"""Compare a saturation-bounded moisture column with the prior RFC 0019 control.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_bounded_column.py SAVE --out DIR
Run only on the pinned A or cold tuning save at density 1, in both wind modes.
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
from pathlib import Path
from time import perf_counter

import numpy as np

import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import (
    build_transport_rates, distribute_rain, solve_bounded_budget, solve_budget,
)
from rfc0019_g3_mixing import RAIN_RATE, SOURCE_DEPTH_M, SOURCE_SECONDS, TUNING_SAVES
from rfc0019_g3_mixing import _coast_distances, _metrics, _save_maps

COLD_COAST_DISTANCE_KM = 500.0


def _downslope_pairs(fields, exposure):
    """Find adjacent windward/leeward samples across strong wind-aligned terrain."""
    windward, leeward, _ = exposure
    row, col = np.where(windward)
    u, v = fields.wind_u[row, col], fields.wind_v[row, col]
    eastward = np.abs(u) >= np.abs(v)
    down_row = row - np.where(eastward, 0, np.sign(v)).astype(int)
    down_col = (col + np.where(eastward, np.sign(u), 0).astype(int)) % windward.shape[1]
    in_bounds = (down_row >= 0) & (down_row < windward.shape[0])
    up, down = (row[in_bounds], col[in_bounds]), (down_row[in_bounds], down_col[in_bounds])
    opposite = leeward[down] & ((fields.wind_u[up] * fields.wind_u[down]
                                + fields.wind_v[up] * fields.wind_v[down]) > 0)
    return tuple(index[opposite] for index in up), tuple(index[opposite] for index in down)


def _pattern_metrics(fields, rain, area, distance, freeze, pairs):
    land = ~fields.is_ocean
    cold_land = land & (fields.air_temperature_c < freeze)
    cold_coast = cold_land & (distance <= COLD_COAST_DISTANCE_KM)
    cold_inland = cold_land & ~cold_coast
    up, down = pairs
    contrast = rain[up] - rain[down]
    cold_regions = {}
    for name, mask in (("cold_coast", cold_coast), ("cold_inland", cold_inland)):
        mean = g3._weighted_mean(rain[mask], area[mask])
        cold_regions[name] = {"area_m2": float(area[mask].sum()), "rain_mm_yr": mean}
    return cold_regions | {
        "adjacent_rain_shadow_pairs": int(contrast.size),
        "pair_median_windward_minus_leeward_mm_yr": float(np.median(contrast)) if contrast.size else None,
        "pair_fraction_windward_wetter": float(np.mean(contrast > 0)) if contrast.size else None,
        "latitude_row_mean_rain_mm_yr": (np.sum(rain * area, axis=1) / np.sum(area, axis=1)).tolist(),
    }


def _arm_ledger(budget, grid, cap, overflow, solves, elapsed):
    column = budget.mass / grid.area
    order = np.argsort(column)
    cumul_area = np.cumsum(grid.area[order])
    cutoff = 0.99 * cumul_area[-1]
    p99_index = np.searchsorted(cumul_area, cutoff)
    tail_share = (cumul_area[p99_index] - cutoff) / grid.area[order[p99_index]]
    tail_rain = budget.rain[order[p99_index + 1:]].sum() + budget.rain[order[p99_index]] * tail_share
    total_rain = float(budget.rain.sum())
    return {
        "source_kg_s": float(sum(values.sum() for values in grid.tagged_source.values())),
        "rain_kg_s": total_rain,
        "column_water_kg_m2": g3._weighted_quantiles(column, grid.area) | {
            "p99": float(column[order[p99_index]]), "max": float(column.max()),
        },
        "rain_fraction_in_highest_column_one_percent": float(tail_rain / total_rain),
        "max_column_to_local_cap": float(np.max(column / cap)),
        "area_fraction_above_cap": float(grid.area[column > cap * (1 + 1e-12)].sum() / grid.area.sum()),
        "overflow_kg_s": float(overflow.sum()),
        "overflow_fraction_of_rain": float(overflow.sum() / total_rain),
        "global_closure": budget.global_relative,
        "summed_cell_closure": budget.summed_cell_relative,
        "free_solves": solves,
        "solve_s_single_call": elapsed,
    }


def _compare(world, wind, out, climate, hydrology, biomes, density):
    world.wind_model = wind
    height, width = climate.grid_dimensions(float(density))
    fields = climate.compute_climate(world, height, width)
    _, _, area = g3._grid_area(height, width)
    _, tags, grid = g3._prepare_source_grid(
        fields, world, climate, hydrology, area, SOURCE_DEPTH_M, SOURCE_SECONDS, density,
    )
    faces = build_transport_rates(grid, fields.wind_u[::-1], fields.wind_v[::-1], diffusivity=0.0)
    source = sum(grid.tagged_source.values())
    rate = np.full_like(source, RAIN_RATE)
    qsat = climate.saturation_specific_humidity_g_kg(fields.air_temperature_c)
    native_cap = qsat[::-1].ravel() * (g3.AIR_DENSITY_KG_M3 * SOURCE_DEPTH_M / 1000.0)
    if not np.all(np.isfinite(native_cap) & (native_cap > 0)):
        raise ValueError("saturation conversion needs positive finite column capacity")
    cap = np.bincount(
        grid.native_ids.ravel(), weights=native_cap * grid.native_area.ravel(), minlength=grid.area.size,
    ) / grid.area

    exposure = g3._windward_leeward_masks(fields, area)
    distance = _coast_distances(fields)
    freeze = hydrology.FREEZE_POINT_C
    pairs = _downslope_pairs(fields, exposure)
    reference = fields.precipitation_mm
    temperature = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = biomes.grid_slope(fields.elevation_m, fields.lat_deg)
    start = perf_counter()
    control = solve_budget(grid.area, source, rate, faces)
    control_time = perf_counter() - start
    start = perf_counter()
    bounded, overflow, solves = solve_bounded_budget(grid.area, source, rate, faces, cap * grid.area)
    bounded_time = perf_counter() - start

    def patterns(rain, ids):
        return (_metrics(fields, rain, area, biomes, ids, exposure, distance, freeze)
                | _pattern_metrics(fields, rain, area, distance, freeze, pairs))

    result = {"wind": wind, "legacy": patterns(reference, fields.biome_ids), "arms": {}}
    panels = [(reference, fields.biome_ids)]
    arrays = {
        "cap_kg_m2": cap, "source_kg_s": source, "lat_deg_north_to_south": fields.lat_deg,
        "wind_u_m_s": fields.wind_u, "wind_v_m_s": fields.wind_v, "is_ocean": fields.is_ocean,
        "air_temperature_c": fields.air_temperature_c, "elevation_m": fields.elevation_m,
        "coast_distance_km": distance, "reference_precipitation_mm_yr": reference,
    }
    for name, budget, extra, steps, elapsed in (
        ("control", control, np.zeros_like(source), 1, control_time),
        ("bounded", bounded, overflow, solves, bounded_time),
    ):
        native_rain_kg_s = distribute_rain(grid, budget.rain, np.ones(grid.native_ids.shape))[::-1]
        rain = native_rain_kg_s / area * g3.SECONDS_PER_YEAR
        ids = biomes.smooth_biome_field(
            temperature, rain, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m,
            lat_deg=fields.lat_deg, axial_tilt_deg=world.axial_tilt_deg,
        )
        panels.append((rain, ids))
        result["arms"][name] = _arm_ledger(budget, grid, cap, extra, steps, elapsed) | {
            "rms_to_current_source": g3._rms_pair(reference, rain, area), "metrics": patterns(rain, ids),
        }
        arrays[f"{name}_mass_kg"] = budget.mass
        arrays[f"{name}_rain_kg_s"] = budget.rain
    np.savez_compressed(out / f"{wind}-columns.npz", area_m2=grid.area, overflow_kg_s=overflow, **arrays)
    _save_maps(out / f"{wind}-rain-and-biomes.png", panels, biomes.BIOME_COLORS,
               labels=("Current source", "Unbounded column", "Saturation-bounded column"),
               difference_label="Bounded − unbounded")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--density", type=int, choices=(1, 4), default=1)
    args = parser.parse_args()
    root, digest, world, biomes, climate, hydrology, resampled, shape = g3._load_tuning_world(
        args.save, parser, TUNING_SAVES, args.density,
    )
    args.out.mkdir(parents=True, exist_ok=False)
    scripts = ("rfc0019_g3_bounded_column.py", "rfc0019_g3_compare.py", "rfc0019_g3_mixing.py",
               "rfc0019_coarse_moisture.py")
    output = {
        "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "scripts": {name: g3._sha256(root / "bin" / "debug" / name) for name in scripts},
        "save": args.save.name,
        "save_sha256": digest,
        "grid_shape": shape,
        "density": args.density,
        "lagged_cache_resampled": resampled,
        "source_depth_m": SOURCE_DEPTH_M,
        "air_density_kg_m3": g3.AIR_DENSITY_KG_M3,
        "cold_coast_distance_km": COLD_COAST_DISTANCE_KM,
        "results": [],
    }
    for wind in ("diagnostic", "cfd"):
        output["results"].append(_compare(world, wind, args.out, climate, hydrology, biomes, args.density))
        output["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        (args.out / "results.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        print(f"{args.save.name}: completed {wind} control and bounded arms", flush=True)


if __name__ == "__main__":
    main()
