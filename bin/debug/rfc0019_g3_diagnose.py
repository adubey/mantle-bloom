"""Run RFC 0019 G3 source, transport, and scoring diagnostics on a pinned tuning save.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_diagnose.py SAVE --out PRIVATE_DIR
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter1d

import rfc0019_g3_compare as g3
from rfc0019_belt_transport import broaden_convergence
from rfc0019_coarse_moisture import build_transport_rates, solve_budget
from rfc0019_g3_mixing import (
    GRID_SHAPE,
    RAIN_RATE,
    SOURCE_DEPTH_M,
    SOURCE_SECONDS,
    TUNING_SAVES,
    _coast_distances,
    _metrics,
)
from rfc0019_rainout_feasibility import bounded_rainout

SIGMA_DEG = (4.0, 8.0)
AMBIENT_RELATIVE_HUMIDITY = 0.8
INTERIOR_BANDS = (0, 500, 1500, np.inf)


def _scoring_groups(fields, reference, area, distance):
    # Include the first row that takes cumulative reference rain past 10%.
    row_rain = np.sum(reference * area, axis=1)
    order = np.argsort(-row_rain, kind="stable")
    count = np.searchsorted(np.cumsum(row_rain[order]), 0.1 * row_rain.sum()) + 1
    wettest = np.zeros(GRID_SHAPE[0], dtype=bool)
    wettest[order[:count]] = True
    lat = fields.lat_deg[:, None]
    groups = {}
    for lo in range(-90, 90, 15):
        band = (lat >= lo) & (lat < lo + 15)
        groups[f"{lo}:{lo + 15}/ocean"] = band & fields.is_ocean
        for index, (start, end) in enumerate(zip(INTERIOR_BANDS[:-1], INTERIOR_BANDS[1:], strict=True)):
            groups[f"{lo}:{lo + 15}/land/{index}"] = band & ~fields.is_ocean & (distance >= start) & (distance < end)
    return groups, wettest


def _aggregate_scores(fields, reference, rain, area):
    native = g3._rms_pair(reference, rain, area)
    scores = {"2": native["raw_rel_rms"]}
    for width in (4, 8):
        # Both bins divide longitude. Latitude keeps symmetric 2° polar caps.
        lat_edges = np.r_[-90, np.arange(-88, 89, width), 90]
        row_ids = np.searchsorted(lat_edges, fields.lat_deg, side="right") - 1
        nlon = 360 // width
        col_ids = np.arange(GRID_SHAPE[1]) // (width // 2)
        ids = (row_ids[:, None] * nlon + col_ids[None, :]).ravel()
        block_area = np.bincount(ids, weights=area.ravel())
        sums = [np.bincount(ids, weights=(values * area).ravel()) for values in (reference, rain)]
        scores[str(width)] = g3._rms_pair(sums[0] / block_area, sums[1] / block_area, block_area)["raw_rel_rms"]
    return native, scores


def _regional_error(reference, rain, area, groups, wettest):
    squared = area * (reference - rain)**2
    total = float(squared.sum())
    by_region = {}
    for label, mask in groups.items():
        if np.any(mask):
            by_region[label] = {
                "area_fraction": float(area[mask].sum() / area.sum()),
                "error_fraction": float(squared[mask].sum() / max(total, 1e-300)),
                "reference_rain_kg_s": float(np.sum(reference[mask] * area[mask]) / g3.SECONDS_PER_YEAR),
                "rain_kg_s": float(np.sum(rain[mask] * area[mask]) / g3.SECONDS_PER_YEAR),
            }
    return float(squared[wettest].sum() / max(total, 1e-300)), by_region


def _latitude_flux(grid, faces, budget, source, ref_kg_s):
    nlat, nlon = grid.coarse_shape
    rows = grid.coarse_cell // nlon
    first_row, second_row = rows[faces.first], rows[faces.second]
    north = second_row > first_row
    modeled = np.bincount(first_row[north], weights=budget.face_flux[north], minlength=nlat - 1)
    normalized_rain = ref_kg_s * (source.sum() / ref_kg_s.sum())
    # Net northward flux through the edge equals north-cap rain minus source.
    row_demand = np.bincount(rows, weights=normalized_rain - source, minlength=nlat)
    required = np.cumsum(row_demand[::-1])[::-1][1:]
    reconstructed = np.cumsum(np.bincount(rows, weights=budget.rain - source, minlength=nlat)[::-1])[::-1][1:]
    np.testing.assert_allclose(modeled, reconstructed, rtol=0, atol=source.sum() * 1e-8)
    edge_deg = np.linspace(-90, 90, nlat + 1)[1:-1]
    return {
        "edge_latitude_deg": edge_deg.tolist(),
        "modeled_northward_kg_s": modeled.tolist(),
        "required_northward_kg_s": required.tolist(),
        "flux_rmse_kg_s": float(np.sqrt(np.mean((modeled - required)**2))),
    }


def _evaporation_only(fields, climate, area, tags):
    sea_sat = climate.saturation_specific_humidity_g_kg(fields.ocean_temperature_c)
    ambient = AMBIENT_RELATIVE_HUMIDITY * climate.saturation_specific_humidity_g_kg(fields.air_temperature_c)
    shape = area * np.hypot(fields.wind_u, fields.wind_v) * np.maximum(sea_sat - ambient, 0.0)
    shape = np.where(fields.is_ocean, shape, 0.0)
    if shape.sum() <= 0.0:
        raise ValueError("bulk evaporation shape has no positive ocean source")
    sea_source = shape * (tags["ocean_target"].sum() / shape.sum())
    np.testing.assert_allclose(sea_source.sum(), tags["ocean_target"].sum(), rtol=1e-14)
    return tags | {"ocean_target": sea_source}


def _arms(v, tags, evaporation):
    mean = v.mean(axis=1, keepdims=True)
    local = v - mean
    yield "control", v, tags
    for sigma in SIGMA_DEG:
        row_sigma = sigma / (180.0 / v.shape[0])
        smooth_mean = gaussian_filter1d(mean, row_sigma, axis=0, mode="nearest", truncate=3)
        smooth_local = gaussian_filter1d(local, row_sigma, axis=0, mode="nearest", truncate=3)
        yield f"zonal_mean_sigma_{sigma:g}", smooth_mean + local, tags
        yield f"local_departures_sigma_{sigma:g}", mean + smooth_local, tags
    yield "evaporation_placement", v, evaporation


def _belt_arms(v, tags):
    correction, geometry = broaden_convergence(v[::-1])
    mean = v.mean(axis=1, keepdims=True)
    smooth_mean = gaussian_filter1d(mean, 4.0, axis=0, mode="nearest", truncate=3)
    yield "control", v, tags, None, None
    yield "zonal_mean_sigma_8", smooth_mean + (v - mean), tags, None, None
    yield "targeted_belts", v, tags, correction, geometry


def _diagnose(world, wind, climate, hydrology, biomes, out, belt_screen=False):
    world.wind_model = wind
    fields = climate.compute_climate(world, *GRID_SHAPE)
    _, _, area = g3._grid_area(*GRID_SHAPE)
    _, tags, _ = g3._prepare_source_grid(
        fields, world, climate, hydrology, area, SOURCE_DEPTH_M, SOURCE_SECONDS,
    )
    exposure = g3._windward_leeward_masks(fields, area)
    distance = _coast_distances(fields)
    freeze = hydrology.FREEZE_POINT_C
    reference = fields.precipitation_mm
    ref_kg_s = (reference * area / g3.SECONDS_PER_YEAR)[::-1].ravel()
    groups, wettest = _scoring_groups(fields, reference, area, distance)
    temperature = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = biomes.grid_slope(fields.elevation_m, fields.lat_deg)
    result = {
        "wind": wind,
        "wettest_latitude_rows": fields.lat_deg[wettest].tolist(),
        "reference_metrics": _metrics(fields, reference, area, biomes, fields.biome_ids, exposure, distance, freeze),
        "arms": {},
    }
    if belt_screen:
        arms = _belt_arms(fields.wind_v, tags)
    else:
        evaporation = _evaporation_only(fields, climate, area, tags)
        arms = ((name, v, arm_tags, None, None) for name, v, arm_tags in _arms(fields.wind_v, tags, evaporation))
    for name, v, arm_tags, delta, geometry in arms:
        grid = g3._aggregate_tagged_sources(fields.is_ocean, arm_tags)
        faces = build_transport_rates(grid, fields.wind_u[::-1], v[::-1], diffusivity=0.0, north_face_delta=delta)
        source = sum(grid.tagged_source.values())
        budget = solve_budget(grid.area, source, np.full_like(source, RAIN_RATE), faces)
        rain = (budget.rain / grid.area).reshape(GRID_SHAPE)[::-1] * g3.SECONDS_PER_YEAR
        ids = biomes.smooth_biome_field(
            temperature, rain, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m,
            lat_deg=fields.lat_deg, axial_tilt_deg=world.axial_tilt_deg,
        )
        native, scales = _aggregate_scores(fields, reference, rain, area)
        wettest_error, partition = _regional_error(reference, rain, area, groups, wettest)
        result["arms"][name] = {
            "native_fit": native,
            "raw_rms_by_grid_deg": scales,
            "wettest_rows_error_fraction": wettest_error,
            "error_partition": partition,
            "latitude_flux": _latitude_flux(grid, faces, budget, source, ref_kg_s),
            "metrics": _metrics(fields, rain, area, biomes, ids, exposure, distance, freeze),
            "global_closure": budget.global_relative,
            "summed_cell_closure": budget.summed_cell_relative,
            "clipped_negative_mass_kg": budget.clipped_negative_mass_kg,
        }
        arrays = {
            "source_kg_s": source, "mass_kg": budget.mass, "rain_kg_s": budget.rain,
            "face_flux_kg_s": budget.face_flux, "face_first": faces.first, "face_second": faces.second,
            "area_m2": grid.area, "is_ocean": grid.is_ocean, "ref_kg_s": ref_kg_s,
        }
        if geometry is not None:
            result["arms"][name]["belt_geometry"] = geometry
            arrays["north_face_delta_m_s"] = delta
        if name == ("targeted_belts" if belt_screen else "control"):
            bounds, witness = bounded_rainout(grid, faces, source, ref_kg_s, budget)
            result["bounded_rainout"] = bounds | {"operator_arm": name}
            if witness is not None:
                arrays["diagnostic_witness_sink_s-1"] = witness
        np.savez_compressed(out / f"{wind}-{name}.npz", **arrays)
        print(f"{wind}: completed {name}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--belt-screen", action="store_true", help="compare targeted convergence against the control and blanket smoothing")
    args = parser.parse_args()
    root, digest, world, biomes, climate, hydrology, resampled, shape = g3._load_tuning_world(
        args.save, parser, TUNING_SAVES,
    )
    args.out.mkdir(parents=True, exist_ok=False)
    scripts = (
        "rfc0019_g3_diagnose.py", "rfc0019_g3_mixing.py", "rfc0019_g3_compare.py",
        "rfc0019_coarse_moisture.py", "rfc0019_rainout_feasibility.py", "rfc0019_belt_transport.py",
    )
    result = {
        "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "scripts": {name: g3._sha256(root / "bin" / "debug" / name) for name in scripts},
        "save": args.save.name,
        "save_sha256": digest,
        "grid_shape": shape,
        "lagged_cache_resampled": resampled,
        "belt_screen": args.belt_screen,
        "results": [],
    }
    for wind in ("diagnostic", "cfd"):
        result["results"].append(_diagnose(world, wind, climate, hydrology, biomes, args.out, args.belt_screen))
        result["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        (args.out / "results.json").write_text(json.dumps(result, allow_nan=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
