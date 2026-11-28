#!/usr/bin/env python3
"""Replay one RFC 0019 G4 arm with the fixed conservative diagnostic atmosphere.

Usage (run separately for each arm and window):
  backend/.venv/bin/python bin/debug/rfc0019_g4_coupling.py CHECKPOINT \
    --anchors ANCHORS.json --myr 0.4 --step-years 100000 \
    --g4-atmosphere conservative --out PRIVATE/RESULT.json [--treatment]

The default atmosphere is the current production climate. The conservative arm replaces only
the moisture calculation during ordinary climate snapshots; CFD forcing requests that skip
moisture still use the original climate function. The RFC 0016 replay retains its normal land,
lake, glacier, river, and target-lake ledgers. No production module is modified.
"""

from __future__ import annotations

from dataclasses import is_dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

import replay_rfc0016_river_water as rfc0016
import rfc0019_g3_compare as g3
import rfc0019_g3_mixing as mixing
import rfc0019_g3_ridge_rainout as ridge
from rfc0019_coarse_moisture import solve_bounded_budget
from rfc0019_g3_mixing import RAIN_RATE, SOURCE_DEPTH_M, SOURCE_SECONDS

ALPHA = 1.0
COAST_LENGTH_KM = 1500.0
RIDGE_INTERCEPTION = 0.25
SOLVER = "finite_volume"
FINGERPRINT_FILES = (
    "bin/debug/rfc0019_g4_coupling.py",
    "bin/debug/replay_rfc0016_river_water.py",
    "bin/debug/replay_lib.py",
    "bin/debug/rfc0016_land_store.py",
    "bin/debug/rfc0019_g3_compare.py",
    "bin/debug/rfc0019_g3_mixing.py",
    "bin/debug/rfc0019_g3_ridge_rainout.py",
    "bin/debug/rfc0019_g3_archive_cap.py",
    "bin/debug/rfc0019_g3_bounded_column.py",
    "bin/debug/rfc0019_coarse_moisture.py",
    "bin/debug/rfc0019_transport_prototypes.py",
    "backend/app/hydroclimate/climate.py",
    "backend/app/hydroclimate/biomes.py",
    "backend/app/hydroclimate/hydrology.py",
    "backend/app/hydroclimate/ice_sheet.py",
    "backend/app/plates.py",
    "backend/app/elevation_lines.py",
    "backend/app/world.py",
)


def _candidate_precipitation(world, fields, climate, hydrology, biomes):
    height, width = fields.is_ocean.shape
    density = height // mixing.GRID_SHAPE[0]
    _, _, native_area = g3._grid_area(height, width)
    proxies = g3._surface_source_proxies(fields, world, climate, hydrology)
    tags, _ = g3._flux_sources(proxies, native_area, SOURCE_DEPTH_M, SOURCE_SECONDS)
    grid = g3._aggregate_tagged_sources(fields.is_ocean, tags, density=density)
    faces = ridge._candidate_faces(grid, fields.wind_u[::-1], fields.wind_v[::-1], SOLVER)

    distance_native = mixing._coast_distances(fields)[::-1].ravel()
    distance_coarse = np.bincount(
        grid.native_ids.ravel(), weights=distance_native * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    sink = RAIN_RATE * np.where(
        grid.is_ocean, 1.0,
        1.0 + ALPHA * np.exp(-distance_coarse / COAST_LENGTH_KM),
    )

    native_cap = (
        climate.saturation_specific_humidity_g_kg(fields.air_temperature_c[::-1]).ravel()
        * (g3.AIR_DENSITY_KG_M3 * SOURCE_DEPTH_M / 1000.0)
    )
    cap = np.bincount(
        grid.native_ids.ravel(), weights=native_cap * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    elevation = fields.elevation_m[::-1].ravel()
    elevation_coarse = np.bincount(
        grid.native_ids.ravel(), weights=elevation * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    crossings = ridge._ridge_crossings(grid, faces, elevation_coarse)
    interception = ridge._interception(
        crossings, faces.first.size, RIDGE_INTERCEPTION,
    )
    source = sum(grid.tagged_source.values())
    started = time.perf_counter()
    budget, overflow, _ = solve_bounded_budget(
        grid.area, source, sink, faces, cap * grid.area, interception=interception,
    )
    rain, placement_residual = ridge._place_surface_rain(
        grid, budget, overflow, native_area[::-1], fields, hydrology.FREEZE_POINT_C,
    )
    elapsed = time.perf_counter() - started

    column_water = budget.mass / grid.area
    native_column_water = column_water[grid.native_ids].reshape(grid.native_ids.shape)[::-1]
    humidity = native_column_water * 1000.0 / (g3.AIR_DENSITY_KG_M3 * SOURCE_DEPTH_M)
    temperature = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = biomes.grid_slope(fields.elevation_m, fields.lat_deg)
    biome_ids = biomes.smooth_biome_field(
        temperature, rain, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m,
        lat_deg=fields.lat_deg, axial_tilt_deg=world.axial_tilt_deg,
    )
    _validate_budget(budget, cap * grid.area, overflow, rain, source, placement_residual)

    diagnostics = {
        "source_kg_s": float(np.sum(source)),
        "budget_rain_kg_s": float(np.sum(budget.rain)),
        "overflow_kg_s": float(np.sum(overflow)),
        "global_relative_residual": float(budget.global_relative),
        "summed_cell_relative_residual": float(budget.summed_cell_relative),
        "placement_relative_residual": float(placement_residual),
        "minimum_column_water_kg_m2": float(np.min(column_water)),
        "maximum_cap_kg_m2": float(np.max(cap)),
        "maximum_cap_fraction": float(np.max(column_water / np.maximum(cap, 1e-300))),
        "finite_nonnegative": True,
        "closure_pass": True,
        "solve_seconds": elapsed,
    }
    return humidity, rain, biome_ids, diagnostics


def _validate_budget(budget, mass_cap, overflow, rain, source, placement_residual):
    mass = budget.mass
    checks = (
        np.all(np.isfinite(mass) & (mass >= 0.0)),
        np.all(mass <= mass_cap + 100.0 * np.finfo(float).eps * mass_cap),
        np.all(np.isfinite(overflow) & (overflow >= 0.0)),
        np.all(np.isfinite(rain) & (rain >= 0.0)),
        budget.global_relative <= 1e-8,
        budget.summed_cell_relative <= 1e-8,
        placement_residual <= 1e-12,
        abs(float(np.sum(budget.rain)) - float(np.sum(source)))
        <= 1e-8 * max(float(np.sum(source)), 1e-300),
    )
    if not all(checks):
        raise RuntimeError(f"G4 atmosphere budget failed validity checks: {checks}")


def _record_climate_fields(fields, climate, biomes):
    height, width = fields.precipitation_mm.shape
    area = g3._grid_area(height, width)[2]
    land = ~fields.is_ocean
    weighted = float(np.sum(area))
    land_area = float(np.sum(area[land]))
    biome_area = np.bincount(
        fields.biome_ids.ravel(), weights=area.ravel(), minlength=len(biomes.BIOME_COLORS),
    )
    return {
        "precipitation_area_mean_mm_yr": float(np.sum(fields.precipitation_mm * area) / weighted),
        "land_precipitation_area_mean_mm_yr": (
            float(np.sum(fields.precipitation_mm[land] * area[land]) / land_area)
            if land_area else None
        ),
        "land_area_fraction": land_area / weighted,
        "humidity_area_p50_g_kg": _weighted_quantile(fields.humidity_g_kg, area, 0.5),
        "humidity_area_p95_g_kg": _weighted_quantile(fields.humidity_g_kg, area, 0.95),
        "biome_area_share": {
            str(index): float(value / weighted)
            for index, value in enumerate(biome_area) if value > 0.0
        },
    }


def _weighted_quantile(values, weights, quantile):
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, quantile * cumulative[-1], side="left")]])


def _erosion_readout(world):
    result = getattr(world, "erosion_cache", None)
    if result is None:
        return None
    budget = getattr(result, "budget", None)
    if budget is not None:
        return {key: float(value) for key, value in budget.items()}
    if is_dataclass(result):
        return {
            "rain_integral": float(np.sum(result.rain)),
            "river_integral": float(np.sum(result.river)),
            "weathering_integral": float(np.sum(result.weathering)),
        }
    return None


def _install_replay_wrapper():
    original_add_arguments = rfc0016.add_arguments
    original_load_replay = rfc0016.load_replay
    original_run_replay = rfc0016.run_replay

    def add_arguments(parser):
        original_add_arguments(parser)
        parser.add_argument("--g4-atmosphere", choices=("legacy", "conservative"), required=True)
        parser.add_argument("--probe-only", action="store_true")

    def load_replay(args, issue):
        world, result = original_load_replay(args, issue)
        result["issue"] = 898
        return world, result

    def run_replay(world, args, result, snapshot, after_step, close_window):
        from app.hydroclimate import biomes, climate, hydrology
        from app import world as world_mod
        from app import plates as plates_mod
        from app.elevation_lines import line_spacing_rad
        from app.hydroclimate import ice_sheet

        original_compute = climate.compute_climate
        original_step_world = world_mod.step_world
        original_compute_hydrology = hydrology.compute_hydrology
        last_fields = [None]
        last_budget = [None]
        last_precipitation = [None]
        climate_call_seconds = [0.0]
        step_climate_seconds = [0.0]

        def capture_frozen_precipitation(world_arg, precipitation_mm, temperature_c, years,
                                         *positional, **keyword):
            hydro = original_compute_hydrology(
                world_arg, precipitation_mm, temperature_c, years, *positional, **keyword,
            )
            amplitude = keyword.get("seasonal_amplitude_at_nodes")
            if amplitude is None and len(positional) >= 2:
                amplitude = positional[1]
            if amplitude is None:
                amplitude = np.zeros_like(temperature_c)
            frozen_fraction = ice_sheet.snow_fraction(temperature_c, amplitude)
            areas = plates_mod.collect_all_accounting_areas_m2(
                hydro.plates_in_order, line_spacing_rad(world_arg.node_density),
            )
            frozen_m3 = float(np.sum(
                np.asarray(precipitation_mm) * frozen_fraction * areas * float(years) / 1000.0
            ))
            liquid_m3 = float(np.sum(
                np.asarray(precipitation_mm) * (1.0 - frozen_fraction) * areas
                * float(years) / 1000.0
            ))
            last_precipitation[0] = {
                "frozen_m3": frozen_m3,
                "liquid_m3": liquid_m3,
                "frozen_share": frozen_m3 / max(frozen_m3 + liquid_m3, 1e-300),
            }
            return hydro

        hydrology.compute_hydrology = capture_frozen_precipitation

        def compute_with_candidate(world_arg, height, width, node_cloud=None, skip_moisture=False):
            started = time.perf_counter()
            if args.g4_atmosphere == "legacy":
                try:
                    fields = original_compute(
                        world_arg, height, width, node_cloud=node_cloud,
                        skip_moisture=skip_moisture,
                    )
                finally:
                    climate_call_seconds[0] += time.perf_counter() - started
                if not skip_moisture:
                    last_fields[0] = fields
                return fields
            try:
                fields = original_compute(
                    world_arg, height, width, node_cloud=node_cloud, skip_moisture=True,
                )
            except Exception:
                climate_call_seconds[0] += time.perf_counter() - started
                raise
            if skip_moisture:
                climate_call_seconds[0] += time.perf_counter() - started
                return fields
            try:
                humidity, rain, biome_ids, diagnostics = _candidate_precipitation(
                    world_arg, fields, climate, hydrology, biomes,
                )
                from dataclasses import replace

                updated = replace(
                    fields, humidity_g_kg=humidity, precipitation_mm=rain, biome_ids=biome_ids,
                )
                last_fields[0], last_budget[0] = updated, diagnostics
                return updated
            finally:
                climate_call_seconds[0] += time.perf_counter() - started

        climate.compute_climate = compute_with_candidate

        def step_world_with_timing(*step_args, **step_kwargs):
            last_fields[0] = None
            last_budget[0] = None
            last_precipitation[0] = None
            climate_call_seconds[0] = 0.0
            try:
                return original_step_world(*step_args, **step_kwargs)
            finally:
                step_climate_seconds[0] = climate_call_seconds[0]

        world_mod.step_world = step_world_with_timing

        result["g4"] = {
            "atmosphere": args.g4_atmosphere,
            "land_store_treatment": bool(args.treatment),
            "step_years": args.step_years,
            "myr": args.myr,
            "checkpoint_myr": args.checkpoint_myr,
            "alpha": ALPHA if args.g4_atmosphere == "conservative" else None,
            "coast_length_km": COAST_LENGTH_KM if args.g4_atmosphere == "conservative" else None,
            "ridge_interception": RIDGE_INTERCEPTION if args.g4_atmosphere == "conservative" else None,
            "solver": SOLVER if args.g4_atmosphere == "conservative" else None,
            "protocol": "analysis/rfc0017/g4-coupling-preregistration.md",
            "anchors_sha256": hashlib.sha256(args.anchors.read_bytes()).hexdigest(),
            "land_store_parameters": {
                "capacity_mm": args.capacity_mm,
                "depletion": args.depletion,
                "demand_factor": args.demand_factor,
                "split_shift": args.split_shift,
                "humidity_source": not args.no_humidity,
            },
            "source_commit": subprocess.run(
                ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2],
                check=True, capture_output=True, text=True,
            ).stdout.strip(),
            "worktree_status": subprocess.run(
                ["git", "status", "--short"], cwd=Path(__file__).resolve().parents[2],
                check=True, capture_output=True, text=True,
            ).stdout.splitlines(),
            "source_hashes": {
                name: hashlib.sha256((Path(__file__).resolve().parents[2] / name).read_bytes()).hexdigest()
                for name in FINGERPRINT_FILES
            },
        }

        if args.probe_only:
            fields = climate.compute_climate(
                world, *climate.grid_dimensions(world.climate_density), skip_moisture=False,
            )
            result["g4_probe"] = {
                "atmosphere": args.g4_atmosphere,
                "fields": _record_climate_fields(fields, climate, biomes),
                "atmosphere_budget": last_budget[0],
            }
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
            climate.compute_climate = original_compute
            hydrology.compute_hydrology = original_compute_hydrology
            world_mod.step_world = original_step_world
            print(f"wrote {args.out}", flush=True)
            return

        def record_step(seconds):
            after_step(seconds)
            step = result["steps"][-1]
            step["climate_compute_seconds"] = step_climate_seconds[0]
            step["precipitation_phase"] = last_precipitation[0]
            if last_fields[0] is not None:
                step["climate_readout"] = _record_climate_fields(last_fields[0], climate, biomes)
            step["erosion_budget_m3"] = _erosion_readout(world)
            if args.g4_atmosphere == "conservative":
                step["atmosphere_budget"] = last_budget[0]
            else:
                step["atmosphere_budget"] = "legacy production climate has no physical source ledger"

        try:
            original_run_replay(world, args, result, snapshot, record_step, close_window)
        finally:
            climate.compute_climate = original_compute
            hydrology.compute_hydrology = original_compute_hydrology
            world_mod.step_world = original_step_world

    rfc0016.add_arguments = add_arguments
    rfc0016.load_replay = load_replay
    rfc0016.run_replay = run_replay


def main():
    _install_replay_wrapper()
    rfc0016.main()


if __name__ == "__main__":
    main()
