#!/usr/bin/env python3
"""Time the frozen p=0.25 coarse ridge candidate against current climate.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_ridge_runtime.py SAVE --out PRIVATE.json

This is an offline candidate-update proxy: production climate code does not yet call the
coarse moisture solver. It measures one warmup plus three updates per solver and wind mode.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

import rfc0019_g3_compare as g3
import rfc0019_g3_ridge_rainout as ridge
from rfc0019_coarse_moisture import solve_bounded_budget
from rfc0019_g3_mixing import RAIN_RATE, SOURCE_DEPTH_M, SOURCE_SECONDS, TUNING_SAVES

WARMUPS = 1
REPEATS = 3
DENSITY = 4


def _candidate_update(world, climate, biomes, hydrology, solver):
    height, width = climate.grid_dimensions(float(DENSITY))
    update_start = time.perf_counter()
    fields, discarded_zero_rain_biome = g3._compute_candidate_fields(
        world, height, width, climate, biomes,
    )
    _, _, area = g3._grid_area(height, width)
    proxies = g3._surface_source_proxies(fields, world, climate, hydrology)
    tags, _ = g3._flux_sources(proxies, area, SOURCE_DEPTH_M, SOURCE_SECONDS)
    grid = g3._aggregate_tagged_sources(fields.is_ocean, tags, density=DENSITY)

    faces = ridge._candidate_faces(grid, fields.wind_u[::-1], fields.wind_v[::-1], solver)

    surface_temp = fields.air_temperature_c[::-1].ravel()
    qsat = climate.saturation_specific_humidity_g_kg(surface_temp)
    native_cap = qsat * (g3.AIR_DENSITY_KG_M3 * SOURCE_DEPTH_M / 1000.0)
    cap = np.bincount(
        grid.native_ids.ravel(), weights=native_cap * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    source = sum(grid.tagged_source.values())
    sink = np.full_like(source, RAIN_RATE)
    crossings = ridge._ridge_crossings(
        grid, faces,
        np.bincount(
            grid.native_ids.ravel(),
            weights=fields.elevation_m[::-1].ravel() * grid.native_area.ravel(),
            minlength=grid.area.size,
        ) / grid.area,
    )
    interception = ridge._interception(crossings, faces.first.size, 0.25)
    budget, overflow, _ = solve_bounded_budget(
        grid.area, source, sink, faces, cap, interception=interception,
    )
    rain, _ = ridge._place_surface_rain(
        grid, budget, overflow, area, fields, hydrology.FREEZE_POINT_C,
    )
    temperature = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = biomes.grid_slope(fields.elevation_m, fields.lat_deg)
    biomes.smooth_biome_field(
        temperature, rain, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m,
        lat_deg=fields.lat_deg, axial_tilt_deg=world.axial_tilt_deg,
    )
    total = time.perf_counter() - update_start - discarded_zero_rain_biome
    return total, budget.summed_cell_relative


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--solver", choices=("both", "finite_volume", "wind_arrows"), default="both",
        help="run one solver per process to attribute peak RSS",
    )
    args = parser.parse_args()
    digest = g3._sha256(args.save)
    if TUNING_SAVES.get(args.save.name) != digest:
        parser.error("save name or SHA-256 does not match the registered G3 tuning inputs")

    backend = Path(__file__).resolve().parents[2] / "backend"
    sys.path.insert(0, str(backend))
    from app import persistence
    from app.hydroclimate import biomes, climate, hydrology

    world = persistence.load_world_bytes(args.save.read_bytes())
    world.climate_density = float(DENSITY)
    target_shape = climate.grid_dimensions(float(DENSITY))
    if world.climate_cache is not None:
        world.climate_cache, cache_resampled = g3._cache_on_grid(world.climate_cache, target_shape)
    else:
        cache_resampled = False

    results = []
    for wind in ("diagnostic", "cfd"):
        world.wind_model = wind
        current_times = []
        for repeat in range(WARMUPS + REPEATS):
            start = time.perf_counter()
            climate.compute_climate(world, *target_shape)
            elapsed = time.perf_counter() - start
            if repeat >= WARMUPS:
                current_times.append(elapsed)
        solvers = ("finite_volume", "wind_arrows") if args.solver == "both" else (args.solver,)
        for solver in solvers:
            candidate_times = []
            closures = []
            for repeat in range(WARMUPS + REPEATS):
                elapsed, closure = _candidate_update(world, climate, biomes, hydrology, solver)
                if repeat >= WARMUPS:
                    candidate_times.append(elapsed)
                    closures.append(closure)
            baseline = float(np.median(current_times))
            candidate = float(np.median(candidate_times))
            results.append({
                "wind": wind,
                "solver": solver,
                "density": DENSITY,
                "shape": list(target_shape),
                "baseline_current_climate_median_seconds": baseline,
                "candidate_update_proxy_median_seconds": candidate,
                "ratio_to_baseline": candidate / baseline,
                "candidate_samples_seconds": candidate_times,
                "closure_max_relative": max(closures),
            })
        print(json.dumps({"save": args.save.name, "wind": wind, "latest": results[-2:]}), flush=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({
            "save": args.save.name,
            "save_sha256": digest,
            "density": DENSITY,
            "partial": wind != "cfd",
            "results": results,
        }, indent=2) + "\n")

    report = {
        "rfc": "0019",
        "measurement": "offline coarse candidate-update runtime proxy",
        "production_integrated": False,
        "commit": subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=backend, check=True,
            capture_output=True, text=True,
        ).stdout.strip(),
        "save": args.save.name,
        "save_sha256": digest,
        "density": DENSITY,
        "lagged_cache_resampled": cache_resampled,
        "warmup_updates": WARMUPS,
        "timed_updates": REPEATS,
        "solver_process": args.solver,
        "peak_rss_mb_process": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
        "results": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
