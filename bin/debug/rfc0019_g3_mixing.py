"""RFC 0019 G3 transport screen on one pinned tuning world.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_mixing.py SAVE --out PRIVATE_DIR
Compare zero mixing with 30,000 m²/s at density 1 under both saved wind modes.
"""

from __future__ import annotations

import argparse
import json
import resource
import subprocess
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import build_transport_rates, distribute_rain, solve_budget

TUNING_SAVES = {
    "seed910211954-221.6myr.bloom": "1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67",
    "seed910211954-225.1myr.bloom": "c5b1f5d832e83e63f2a90a9f9f7bff03c885e831a0c6e585e22678a718168043",
    "seed910211954-224myr.bloom": "3e0fa48182b86b7d7e293cf73c180c6a4667f1b4c48dac7effaf57c9d2b87ca7",
    "mantle-bloom-seed317739385-200800000y.bloom": "c7e8c4094376b2b57158e36b1066cf529eb475ac0bbf54a6b888b391f294d344",
}
GRID_SHAPE = (90, 180)
DIFFUSIVITIES = (0.0, 30_000.0)
SOURCE_DEPTH_M = 2_000.0
SOURCE_SECONDS = 9.1 * 86_400.0
RAIN_RATE = 1 / (9.1 * 86_400.0)


def _coast_distances(fields):
    land = ~fields.is_ocean
    chord, _ = cKDTree(fields.world_xyz[~land]).query(fields.world_xyz[land])
    distance = np.zeros(land.shape)
    distance[land] = 2 * g3.EARTH_RADIUS_M * np.arcsin(np.minimum(chord / 2, 1)) / 1000
    return distance


def _regional_metrics(fields, rain, area, distance, freeze_point):
    land = ~fields.is_ocean
    total_land_rain = np.sum(rain[land] * area[land])
    cold = land & (fields.air_temperature_c < freeze_point)
    distances = {}
    for name, lower, upper in (
        ("0_to_500", 0, 500), ("500_to_1500", 500, 1500), ("1500_plus", 1500, np.inf),
    ):
        mask = land & (distance >= lower) & (distance < upper)
        distances[name] = {
            "rain_mm_yr": g3._weighted_mean(rain[mask], area[mask]),
            "fraction_of_land": float(area[mask].sum() / area[land].sum()),
        }
    return {
        "land_rain_by_distance_to_ocean_km": distances,
        "below_freezing_land_rain_mm_yr": g3._weighted_mean(rain[cold], area[cold]),
        "fraction_of_land_rain_below_freezing": float(np.sum(rain[cold] * area[cold]) / total_land_rain),
    }


def _face_diagnostics(grid, faces, budget, ocean_origin):
    cross = grid.is_ocean[faces.first] != grid.is_ocean[faces.second]
    signed = np.where(grid.is_ocean[faces.first], budget.face_flux, -budget.face_flux)[cross]
    nlat, nlon = grid.coarse_shape
    row = grid.coarse_cell // nlon
    south_row = row[faces.first]
    northward = row[faces.second] > south_row
    edges = np.linspace(-90, 90, nlat + 1)
    edge_south_row = np.searchsorted(edges, [-60, 60]) - 1
    if not np.allclose(edges[edge_south_row + 1], [-60, 60]):
        raise ValueError("the 60 degree diagnostics require latitude faces at ±60 degrees")
    return {
        "gross_ocean_to_land_kg_s": float(np.maximum(signed, 0).sum()),
        "gross_land_to_ocean_kg_s": float(np.maximum(-signed, 0).sum()),
        "net_ocean_to_land_kg_s": float(signed.sum()),
        "ocean_origin_land_rain_kg_s": float(ocean_origin.rain[~grid.is_ocean].sum()),
        "poleward_across_60s_kg_s": float(-budget.face_flux[northward & (south_row == edge_south_row[0])].sum()),
        "poleward_across_60n_kg_s": float(budget.face_flux[northward & (south_row == edge_south_row[1])].sum()),
    }


def _solve_arm(fields, grid, diffusivity):
    faces = build_transport_rates(grid, fields.wind_u[::-1], fields.wind_v[::-1], diffusivity)
    sink = np.full(grid.area.shape, RAIN_RATE)
    total = sum(grid.tagged_source.values())
    budget = solve_budget(grid.area, total, sink, faces)
    ocean_origin = solve_budget(grid.area, grid.tagged_source["ocean_target"], sink, faces)
    rain = distribute_rain(grid, budget.rain, np.ones(GRID_SHAPE))[::-1]
    rain = rain / grid.native_area[::-1] * g3.SECONDS_PER_YEAR
    return rain, {
        "diffusivity_m2_s": diffusivity,
        "total_source_kg_s": float(total.sum()),
        "total_rain_kg_s": float(budget.rain.sum()),
        "global_closure_relative": budget.global_relative,
        "summed_cell_closure_relative": budget.summed_cell_relative,
        "max_cell_closure_relative": budget.max_cell_relative,
        "face_transport": _face_diagnostics(grid, faces, budget, ocean_origin),
    }


def _metrics(fields, rain, area, biomes, biome_ids, exposure, distance, freeze_point):
    stats = g3._rain_metrics(fields, rain, area, biomes, biome_ids=biome_ids)
    windward, leeward, _ = exposure
    return stats | g3._extended_rain_diagnostics(
        fields, rain, area, windward, leeward, stats["koppen_group_shares"]["B"],
    ) | _regional_metrics(fields, rain, area, distance, freeze_point)


def _save_maps(path, panels, palette, labels=("Current source", "Zero mixing", "30,000 m²/s"),
               difference_label="Mixing − control"):
    height, width = panels[0][0].shape
    maximum = float(np.percentile(np.stack([rain for rain, _ in panels]), 99.5))
    diff = panels[2][0] - panels[1][0]
    diff_max = float(np.percentile(np.abs(diff), 99.5))
    panel_width, panel_height = width * 3, height * 3
    image = Image.new("RGB", (4 * panel_width, 2 * panel_height + 56), (248, 248, 248))
    draw = ImageDraw.Draw(image)
    for index, ((rain, biome), label) in enumerate(zip(panels, labels, strict=True)):
        offset = index * panel_width
        draw.text((offset + 10, 8), label, fill=(20, 20, 20))
        for y, pixels in ((28, g3._rain_map_color(rain, maximum)), (panel_height + 56, palette[biome])):
            raster = Image.fromarray(pixels).resize((panel_width, panel_height), Image.Resampling.NEAREST)
            image.paste(raster, (offset, y))
    draw.text((3 * panel_width + 10, 8), f"{difference_label}, ±{diff_max:.1f} mm/yr", fill=(20, 20, 20))
    change = Image.fromarray(g3._difference_map_color(diff, diff_max))
    image.paste(change.resize((panel_width, panel_height), Image.Resampling.NEAREST), (3 * panel_width, 28))
    draw.text((10, panel_height + 35), f"Shared log rain scale to {maximum:.0f} mm/yr. Biomes below.",
              fill=(20, 20, 20))
    image.save(path)


def _compare(world, wind, out, climate, hydrology, biomes):
    world.wind_model = wind
    fields = climate.compute_climate(world, *GRID_SHAPE)
    _, _, area = g3._grid_area(*GRID_SHAPE)
    _, tags, grid = g3._prepare_source_grid(
        fields, world, climate, hydrology, area, SOURCE_DEPTH_M, SOURCE_SECONDS,
    )
    exposure = g3._windward_leeward_masks(fields, area)
    distance = _coast_distances(fields)
    reference = fields.precipitation_mm
    freeze = hydrology.FREEZE_POINT_C
    current = _metrics(fields, reference, area, biomes, fields.biome_ids, exposure, distance, freeze)
    temperature = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = biomes.grid_slope(fields.elevation_m, fields.lat_deg)
    panels = [(reference, fields.biome_ids)]
    arms = []
    for diffusivity in DIFFUSIVITIES:
        rain, ledger = _solve_arm(fields, grid, diffusivity)
        ids = biomes.smooth_biome_field(
            temperature, rain, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m,
            lat_deg=fields.lat_deg, axial_tilt_deg=world.axial_tilt_deg,
        )
        arms.append(ledger | {
            "rms_to_current_source": g3._rms_pair(reference, rain, area),
            "metrics": _metrics(fields, rain, area, biomes, ids, exposure, distance, freeze),
        })
        panels.append((rain, ids))
    map_path = out / f"{wind}-rain-and-biomes.png"
    _save_maps(map_path, panels, biomes.BIOME_COLORS)
    return {
        "wind": wind,
        "legacy": current,
        "windward_leeward_mask": exposure[2],
        "source_kg_s_by_tag": {name: float(values.sum()) for name, values in grid.tagged_source.items()},
        "arms": arms,
        "rms_mixing_to_control": g3._rms_pair(panels[1][0], panels[2][0], area)["raw_rel_rms"],
        "map": map_path.name,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    root, digest, world, biomes, climate, hydrology, cache_resampled, shape = g3._load_tuning_world(
        args.save, parser, TUNING_SAVES,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    scripts = ("rfc0019_g3_mixing.py", "rfc0019_g3_compare.py", "rfc0019_coarse_moisture.py")
    output = {
        "base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "script_sha256": {name: g3._sha256(root / "bin" / "debug" / name) for name in scripts},
        "save": args.save.name,
        "save_sha256": digest,
        "shape": shape,
        "lagged_climate_cache_resampled": cache_resampled,
        "source_depth_m": SOURCE_DEPTH_M,
        "source_timescale_days": SOURCE_SECONDS / 86_400,
        "rainout_days": 1 / (RAIN_RATE * 86_400),
        "wind_multiplier": 1.0,
        "results": [],
    }
    for wind in ("diagnostic", "cfd"):
        output["results"].append(_compare(world, wind, args.out, climate, hydrology, biomes))
        output["peak_rss_mb_process"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        (args.out / "results.json").write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        print(f"{args.save.name}: completed {wind} control and mixing arms", flush=True)


if __name__ == "__main__":
    main()
