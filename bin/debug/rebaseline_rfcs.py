#!/usr/bin/env python3
"""Measure RFC 0001–0007 baselines once per seed, serially (issue #403).

Uses compare_cratons.run_job unchanged: default generated quad worlds, cratons on,
5 Myr steps and 100 Myr checkpoints. Read-only taps also capture continental budgets,
each collision merge's craton loss, and attrite's contribution to subducted_m3.
Saves are scratch artifacts; only compact JSON belongs in analysis/issue403/.

Usage:
    backend/.venv/bin/python bin/debug/rebaseline_rfcs.py --out /tmp/tuning.json \
        --save-dir /tmp/issue403
    backend/.venv/bin/python bin/debug/rebaseline_rfcs.py --seeds 456862899 \
        --years 3e9 --out /tmp/long.json --save-dir /tmp/issue403
    backend/.venv/bin/python bin/debug/ledger_history.py \
        /tmp/issue403/seed456862899.mbworld --json /tmp/history.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path

import compare_cratons
import ledger_history
from app import continental_ledger, cratons, merge_split, persistence


def _write(path: Path, result: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=1) + "\n")


def _controls(world) -> dict:
    names = (*compare_cratons.world_mod.TUNING_MULTIPLIER_FIELDS,
             "node_density", "climate_density", "fluid_density", "axial_tilt_deg",
             "fault_deformation_mode", "simulate_climate_biomes", "simulate_plate_movement", "debug_diagnostics")
    return {name: getattr(world, name) for name in names}


def _windows(rows: list[dict]) -> list[dict]:
    """Exact inclusive endpoints, normalized by physical duration and window-start live volume."""
    by_myr = {row["myr"]: row for row in rows}
    windows = []
    for start in range(0, int(rows[-1]["myr"]), 1000):
        first, last = by_myr[start], by_myr[min(start + 1000, rows[-1]["myr"])]
        gyr = (last["myr"] - first["myr"]) / 1000
        volume = first["craton_volume_km3"]
        delta = {key: last["craton_ledger_km3"][key] - first["craton_ledger_km3"][key]
                 for key in (*ledger_history.CRATON_SINKS, *ledger_history.RIFT_PATHS, "formed_m3")}
        initial_material = first["continental_material_ledger"]["surface_continental_derived_m3"]
        final_material = last["continental_material_ledger"]["surface_continental_derived_m3"]
        material_flux = {key: (last["continental_material_ledger"][key] - first["continental_material_ledger"][key]) / 1e18 / gyr
                         for key in continental_ledger.LEDGER_KEYS if key != "initial_continental_m3"}
        windows.append({
            "window_myr": [first["myr"], last["myr"]],
            "start_craton_volume_km3": volume,
            "craton_delta_km3": delta,
            "craton_share_of_start_per_gyr": {key: value / volume / gyr if volume > 0 else None for key, value in delta.items()},
            "events": {key: last["events"].get(key, 0) - first["events"].get(key, 0)
                       for key in compare_cratons.EVENT_KEYS},
            "plate_count_start_end": [first["plate_count"], last["plate_count"]],
            "continental_inventory_start_end_km3": [initial_material / 1e9, final_material / 1e9],
            "continental_net_change_fraction_per_gyr": (final_material / initial_material - 1) / gyr,
            "continental_flux_km3_per_year": material_flux,
            "attrited_craton_share_of_start_per_gyr": (last["attrited_craton_km3"] - first["attrited_craton_km3"]) / volume / gyr
                if volume > 0 else None,
        })
    return windows


def _formation(history: list[dict]) -> dict:
    positive = [last["elapsed_years"] / 1e6 for first, last in zip(history, history[1:])
                if last["craton_ledger_km3"].get("formed_m3", 0.0) > first["craton_ledger_km3"].get("formed_m3", 0.0)]
    return {"last_positive_sample_myr": max(positive, default=None),
            "positive_samples": len(positive), "sample_interval_years": 5e6}


def _history_fields(history: list[dict], rows: list[dict]) -> None:
    """Reuse already recorded stats for RFC 0006's geometry comparator; never step again."""
    by_year = {row["elapsed_years"]: row for row in history}
    for row in rows:
        saved = by_year[row["myr"] * 1e6]
        row["plate_overlap_area_fraction"] = saved["plate_overlap_area_fraction"]
        row["land_fraction_area_weighted"] = saved["land_fraction"]


def _merge_booked_m3(world) -> float:
    """Craton a merge's cap has booked: `delaminated_m3` before #441, `merge_capped_m3` after,
    so the tap measures either backend."""
    ledger = world.craton_ledger
    return ledger.get("delaminated_m3", 0.0) + ledger.get("merge_capped_m3", 0.0)


def _measure(seed: int, years: float, save_dir: Path) -> dict:
    original_snapshot = compare_cratons._snapshot
    original_merge = merge_split._merge_plates_with_craton_audit
    original_attrite = cratons.attrite
    merges = []
    attrited = 0.0
    captured = {}

    def merge(world, keep: int, absorb: int) -> None:
        before = _merge_booked_m3(world)
        live = cratons.live_volume_m3(world)
        original_merge(world, keep, absorb)
        loss = _merge_booked_m3(world) - before
        merges.append({"myr": world.elapsed_years / 1e6, "pair": [keep, absorb],
                       "craton_lost_km3": loss / 1e9,
                       "share_of_live": loss / live if live else 0.0})

    def attrite(*args, **kwargs):
        nonlocal attrited
        lost = original_attrite(*args, **kwargs)
        attrited += lost
        return lost

    def snapshot(world, events):
        captured["world"] = world
        return {**original_snapshot(world, events),
                "continental_material_ledger": continental_ledger.inventories(world),
                "attrited_craton_km3": attrited / 1e9}

    compare_cratons._snapshot = snapshot
    merge_split._merge_plates_with_craton_audit = merge
    cratons.attrite = attrite
    try:
        rows = compare_cratons.run_job(seed, "on", years, 100e6, 5e6, {})
        world = captured["world"]
        _history_fields(world.stats_history, rows)
        save_dir.mkdir(parents=True, exist_ok=True)
        save = save_dir / f"seed{seed}.mbworld"
        data = persistence.save_world_bytes(world)
        save.write_bytes(data)
        return {"seed": seed, "controls": _controls(world), "checkpoints": rows,
                "per_gyr": _windows(rows), "collision_merges": merges,
                "ledger_history": ledger_history.summarize(world.stats_history),
                "formation": _formation(world.stats_history),
                "save_bytes": len(data),
                "save_sha256": hashlib.sha256(data).hexdigest()}
    finally:
        compare_cratons._snapshot = original_snapshot
        merge_split._merge_plates_with_craton_audit = original_merge
        cratons.attrite = original_attrite


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--save-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    model_diff = subprocess.check_output(["git", "diff", "HEAD", "--", "backend/app", "backend/constraints.txt"], cwd=root, text=True)
    if model_diff:
        parser.error("baseline needs an unchanged model and constraints")
    if args.years <= 0 or args.years % 100e6:
        parser.error("--years must be a positive multiple of 100 Myr")
    result = {"baseline_commit": commit, "python": platform.python_version(),
              "constraints_sha256": hashlib.sha256((root / "backend/constraints.txt").read_bytes()).hexdigest(),
              "generation": "world.generate_world(seed), all other arguments default; no control overrides",
              "step_years": 5e6, "every_years": 100e6, "years": args.years,
              "split_min_age_years": merge_split.SPLIT_MIN_AGE_YEARS, "runs": []}
    for seed in map(int, args.seeds.split(",")):
        result["runs"].append(_measure(seed, args.years, args.save_dir))
        _write(args.out, result)


if __name__ == "__main__":
    main()
