#!/usr/bin/env python3
"""Issue #355: replay river and glacier calibration on the same save in two code trees.

Create one starting save per density with the baseline backend, then replay each save on
origin/main and the issue branch. Full outputs belong in scratch space; commit only a compact
comparison to analysis/issue355/.

    backend/.venv/bin/python bin/debug/measure_discharge_calibration.py \
      /tmp/issue355-density1.mbworld --backend /tmp/mantle-bloom-issue355-main/backend \
      --generate-seed 355004 --node-density 1
    backend/.venv/bin/python bin/debug/measure_discharge_calibration.py \
      /tmp/issue355-density1.mbworld --backend /tmp/mantle-bloom-issue355-main/backend \
      --myr 10 --checkpoint-myr 1 --out /tmp/issue355-main-density1.json
    backend/.venv/bin/python bin/debug/measure_discharge_calibration.py \
      /tmp/issue355-density1.mbworld --myr 10 --checkpoint-myr 1 \
      --out /tmp/issue355-branch-density1.json

The per-step erosion arrays are raw stream-power/ice-abrasion capacity before material and
column caps; the erosion budget records total removed/deposited and applied glacial erosion.
Channel percentiles are area-weighted over positive-width/depth land cells, with a separate
land-area coverage fraction. A replay is sequential and reports runtime and peak RSS.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def weighted_percentile(values: np.ndarray, weights: np.ndarray, percentile: float) -> float:
    if not len(values):
        return 0.0
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, percentile / 100.0 * cumulative[-1])]])


def surface_metrics(world, plate_arrays, line_spacing_rad, hydrology) -> dict:
    plates = world.plates
    area = plate_arrays.collect_all_accounting_areas_m2(plates, line_spacing_rad(world.node_density))
    elevation = plate_arrays.collect_all_elevation(plates)
    depth = plate_arrays.collect_all_channel_depth(plates)
    width = plate_arrays.collect_all_channel_width(plates)
    ice = plate_arrays.collect_all_glacier_depth(plates)
    lake = plate_arrays.collect_all_lake_depth(plates)
    hydro = world.hydrology_cache
    land = hydro.is_ocean == 0 if hydro is not None else elevation > world.sea_level_m
    land_area = max(float(area[land].sum()), 1.0)

    def channel(values: np.ndarray) -> dict:
        active = land & (values > 0.0)
        return {
            "coverage_pct_of_land": 100.0 * float(area[active].sum()) / land_area,
            "p50_m": weighted_percentile(values[active], area[active], 50.0),
            "p90_m": weighted_percentile(values[active], area[active], 90.0),
            "p99_m": weighted_percentile(values[active], area[active], 99.0),
            "max_m": float(values[active].max()) if np.any(active) else 0.0,
        }

    result = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "land_pct": 100.0 * float(area[land].sum()) / float(area.sum()),
        "channel_depth": channel(depth),
        "channel_width": channel(width),
        "ice_volume_km3": float(np.dot(area, ice)) / 1e9,
        "ice_cap_pct_of_land": 100.0 * float(area[land & (ice >= hydrology.GLACIER_VISIBLE_DEPTH_M)].sum()) / land_area,
        "lake_volume_km3": float(np.dot(area, lake)) / 1e9,
    }
    erosion = getattr(world, "erosion_cache", None)
    if erosion is not None and len(erosion.is_river_depositing) == len(area):
        result["river_deposition_area_pct_of_land"] = (
            100.0 * float(area[erosion.is_river_depositing & land].sum()) / land_area
        )
    if hydro is not None:
        result["river_area_pct_of_land"] = 100.0 * float(area[hydro.is_river & land].sum()) / land_area
        discharge = getattr(hydro, "discharge_m3_s", None)
        if discharge is not None and len(discharge) == len(area):
            flowing = land & (discharge > 0.0)
            result["flowing_area_pct_of_land"] = 100.0 * float(area[flowing].sum()) / land_area
            result["discharge_p99_m3_s"] = weighted_percentile(discharge[flowing], area[flowing], 99.0)
            result["discharge_max_m3_s"] = float(discharge[flowing].max()) if np.any(flowing) else 0.0
            if erosion is not None and hasattr(hydrology, "discharge_m3_s_to_flow_accum_mm_yr"):
                flow_m = hydrology.discharge_m3_s_to_flow_accum_mm_yr(discharge) / 1000.0
                speed = hydrology.compute_river_speed(erosion.slope, discharge)
                slow = speed < 2.0
                result["deposition_gate_sweep_pct_of_land"] = {
                    str(cutoff): 100.0 * float(area[land & slow & (flow_m > cutoff)].sum()) / land_area
                    for cutoff in (0.005, 0.01, 0.02, 0.03, 0.05)
                }
                result["deposition_speed_sweep_pct_of_land"] = {
                    str(cutoff): 100.0 * float(area[land & (speed < cutoff) & (flow_m > 0.05)].sum()) / land_area
                    for cutoff in (1.0, 2.0, 3.0, 4.0, 8.0)
                }
    return result


def replay(args, world, plate_arrays, line_spacing_rad, erosion_mod, hydrology, world_mod, save_sha: str, commit: str) -> dict:
    result = {
        "issue": 355,
        "backend_commit": commit,
        "save": str(args.save),
        "save_sha256": save_sha,
        "node_density": world.node_density,
        "step_years": args.step_years,
        "river_erosion_coefficient": erosion_mod.RIVER_EROSION_COEFFICIENT,
        "width_growth_coefficient": erosion_mod.WIDTH_GROWTH_COEFFICIENT,
        "checkpoints": [surface_metrics(world, plate_arrays, line_spacing_rad, hydrology)],
    }
    steps = round(args.myr * 1e6 / args.step_years)
    every = max(1, round(args.checkpoint_myr * 1e6 / args.step_years))
    totals = {key: 0.0 for key in (
        "river_erosion_raw_km3", "glacial_erosion_raw_km3", "glacial_erosion_applied_km3",
        "removed_km3", "deposited_km3", "ice_accumulation_km3", "ice_frozen_lake_km3",
        "ice_meltwater_km3", "ice_frozen_river_km3", "ice_calving_km3", "step_seconds",
    )}
    worst_ice_residual_m3 = 0.0
    for step in range(1, steps + 1):
        started = time.perf_counter()
        world_mod.step_world(world, args.step_years)
        totals["step_seconds"] += time.perf_counter() - started
        erosion = world.erosion_cache
        hydro = world.hydrology_cache
        if erosion is not None and hydro is not None:
            area = plate_arrays.collect_all_accounting_areas_m2(
                hydro.plates_in_order, line_spacing_rad(world.node_density)
            )
            land = ~hydro.is_ocean
            totals["river_erosion_raw_km3"] += float(np.dot(erosion.river[land], area[land])) / 1e9
            sliding = hydro.ice.sliding_speed_m_per_yr if hydro.ice is not None else np.zeros(len(area))
            glacier_raw = (
                erosion_mod.GLACIER_ABRASION_COEFFICIENT
                * world.glacier_erosion_multiplier
                * np.power(sliding, erosion_mod.GLACIER_ABRASION_EXPONENT)
                * args.step_years / 1e6
            )
            totals["glacial_erosion_raw_km3"] += float(np.dot(glacier_raw[land], area[land])) / 1e9
            budget = erosion.budget or {}
            for key, budget_key in (
                ("glacial_erosion_applied_km3", "glacial_erosion_m3"),
                ("removed_km3", "removed_m3"),
                ("deposited_km3", "deposited_m3"),
            ):
                totals[key] += budget.get(budget_key, 0.0) / 1e9
            if hydro.ice is not None:
                ledger = hydro.ice.ledger
                totals["ice_accumulation_km3"] += ledger.accumulation_m3 / 1e9
                totals["ice_frozen_lake_km3"] += ledger.frozen_lake_m3 / 1e9
                totals["ice_meltwater_km3"] += ledger.meltwater_m3 / 1e9
                totals["ice_frozen_river_km3"] += ledger.frozen_river_ice_m3 / 1e9
                totals["ice_calving_km3"] += ledger.calved_m3 / 1e9
                worst_ice_residual_m3 = max(worst_ice_residual_m3, abs(ledger.residual_m3))
        if step % every == 0 or step == steps:
            snap = surface_metrics(world, plate_arrays, line_spacing_rad, hydrology)
            snap.update(totals)
            snap["step_seconds_mean"] = totals["step_seconds"] / step
            snap["ice_worst_residual_m3"] = worst_ice_residual_m3
            snap["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
            result["checkpoints"].append(snap)
            print(json.dumps(snap), flush=True)
            args.out.write_text(json.dumps(result, indent=1))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--generate-seed", type=int)
    parser.add_argument("--node-density", type=float)
    parser.add_argument("--myr", type=float, default=10.0)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--checkpoint-myr", type=float, default=1.0)
    parser.add_argument("--river-erosion-coefficient", type=float)
    parser.add_argument("--width-growth-coefficient", type=float)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    backend = args.backend.resolve()
    sys.path.insert(0, str(backend))
    from app import persistence, plates as plate_arrays, world as world_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import erosion, hydrology  # noqa: E402

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=backend, check=True, capture_output=True, text=True).stdout.strip()
    if args.generate_seed is not None:
        if args.node_density is None:
            parser.error("--node-density is required with --generate-seed")
        world = world_mod.generate_world(seed=args.generate_seed, node_density=args.node_density)
        args.save.parent.mkdir(parents=True, exist_ok=True)
        args.save.write_bytes(persistence.save_world_bytes(world))
        print(json.dumps({"save": str(args.save), "seed": world.seed, "node_density": world.node_density, "commit": commit}))
        return
    if args.out is None:
        parser.error("--out is required for a replay")
    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    if args.river_erosion_coefficient is not None:
        erosion.RIVER_EROSION_COEFFICIENT = args.river_erosion_coefficient
    if args.width_growth_coefficient is not None:
        erosion.WIDTH_GROWTH_COEFFICIENT = args.width_growth_coefficient
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result = replay(
        args, world, plate_arrays, line_spacing_rad, erosion, hydrology, world_mod,
        hashlib.sha256(data).hexdigest(), commit,
    )
    args.out.write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
