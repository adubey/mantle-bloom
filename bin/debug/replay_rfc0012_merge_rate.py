#!/usr/bin/env python3
"""Run one RFC 0012 #547 merge-only replay from a generated world.

This is one arm of the preregistered baseline/candidate matrix in
`analysis/issue547/README.md`. Invoke it once per seed and step size, serially across
baseline and candidate worktrees. It writes checkpoint summaries and cumulative motion,
deformation and ledger measures; it does not save full worlds.

    PYTHONHASHSEED=0 backend/.venv/bin/python bin/debug/replay_rfc0012_merge_rate.py \
        --seed 326494039 --step-myr 0.1 --myr 50 --out <scratch-output>.json
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

backend_root = Path(
    os.environ.get("MANTLE_BLOOM_REPLAY_BACKEND", str(Path(__file__).resolve().parents[2] / "backend"))
).resolve()
sys.path.insert(0, str(backend_root))

from app import continental_ledger, cratons, elevation_lines, hm_ledger, merge_split, world as world_mod
from app.hydroclimate import mobile_cover


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if not len(values):
        return float("nan")
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, q * cumulative[-1])]])


def _snapshot(world) -> dict:
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    hc, hm, area, land = [], [], [], []
    for plate in world.plates:
        if not plate.node_count():
            continue
        hc.append(np.asarray(plate.collect("crustal_thickness_m"), dtype=float))
        hm.append(np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float))
        area.append(np.asarray(plate.accounting_areas_m2(spacing), dtype=float))
        land.append(np.asarray(plate.collect("elevation"), dtype=float) > world.sea_level_m)
    hc, hm, area, land = (np.concatenate(parts) for parts in (hc, hm, area, land))
    total_area = float(area.sum())
    return {
        "myr": world.elapsed_years / 1e6,
        "plate_count": len(world.plates),
        "land_fraction": float(area[land].sum() / total_area) if total_area else 0.0,
        "hc_area_weighted_km": {
            f"p{int(q * 100)}": _weighted_quantile(hc, area, q) / 1e3 for q in (0.5, 0.9, 0.99)
        },
        "hm_area_weighted_km": {
            f"p{int(q * 100)}": _weighted_quantile(hm, area, q) / 1e3 for q in (0.5, 0.9, 0.99)
        },
        "ledger_balance_error_m3": {
            "continental": continental_ledger.balance_error_m3(world),
            "craton": cratons.balance_error_m3(world),
            "mobile_cover": mobile_cover.balance_error_m3(world),
            "hm_by_scope": {},
        },
        "hm_inventory_m3": hm_ledger.inventory_scopes_m3(world),
        "hm_accounts_m3": hm_ledger.cumulative_scopes(world),
    }


def _phase_totals(budget: dict) -> dict:
    totals = {}
    for phase, entry in budget.items():
        scopes = entry.get("scopes")
        if scopes and "all" in scopes:
            row = scopes["all"]
            totals[phase] = {
                "hc_m3": row["hc_volume_after_m3"] - row["hc_volume_before_m3"],
                "hm_m3": row["hm_volume_after_m3"] - row["hm_volume_before_m3"],
            }
    return totals


def _hm_closure(snapshot: dict, initial_inventory: dict) -> dict:
    return {
        scope: snapshot["hm_inventory_m3"][scope]
        - initial_inventory[scope]
        - snapshot["hm_accounts_m3"][scope]["source_m3"]
        + snapshot["hm_accounts_m3"][scope]["sink_m3"]
        - snapshot["hm_accounts_m3"][scope]["reclassification_m3"]
        for scope in initial_inventory
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--step-myr", type=float, choices=(0.1, 1.0, 2.0), required=True)
    parser.add_argument("--myr", type=float, default=50.0)
    parser.add_argument("--checkpoint-myr", type=float, default=10.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=backend_root.parent, check=True, capture_output=True, text=True
    ).stdout.strip()
    world = world_mod.generate_world(seed=args.seed, node_density=4.0)
    world.debug_diagnostics = True
    cratons.ensure_initialized(world)
    continental_ledger.ensure_initialized(world)
    mobile_cover.ensure_ledger(world)
    initial = _snapshot(world)
    initial_hm = initial["hm_inventory_m3"]
    result = {
        "commit": commit,
        "seed": args.seed,
        "node_density": 4.0,
        "step_myr": args.step_myr,
        "horizon_myr": args.myr,
        "checkpoints": [initial],
        "topology_events": [],
        "motion": defaultdict(float),
        "deformation": defaultdict(lambda: defaultdict(float)),
        "deformation_window": defaultdict(lambda: defaultdict(float)),
    }
    seconds = args.step_myr * 1e6
    steps = round(args.myr / args.step_myr)
    checkpoint_steps = round(args.checkpoint_myr / args.step_myr)
    previous_phases = _phase_totals(world.phase_budget)
    original_topology = merge_split.apply_topology_changes

    def capture_topology_events(step_world, years):
        for message in original_topology(step_world, years):
            result["topology_events"].append([step_world.elapsed_years, message])
            yield message

    merge_split.apply_topology_changes = capture_topology_events
    started = time.time()
    for index in range(1, steps + 1):
        motion = []
        world_mod.step_world(world, seconds, motion_diagnostics=motion)
        if motion:
            diag = motion[-1]
            result["motion"]["accumulated_driving_torque_nm_year"] += diag["driving_torque_norm_nm"] * seconds
            result["motion"]["accumulated_residual_torque_nm_year"] += diag["residual_torque_norm_nm"] * seconds
            rates = [row["rate_cm_per_yr"] for row in diag["plate_rates_cm_per_yr"]]
            result["motion"]["rate_sum_cm_per_yr"] += sum(rates)
            result["motion"]["rate_samples"] += len(rates)
            result["motion"]["peak_solved_rate_cm_per_yr"] = max(
                result["motion"]["peak_solved_rate_cm_per_yr"], max(rates, default=0.0)
            )
            contacts = diag["blocked_contact_responses"]
            result["motion"]["contact_samples"] += len(contacts)
            result["motion"]["absolute_closing_sum_m_per_s"] += sum(abs(c["closing_rate_m_per_s"]) for c in contacts)
            result["motion"]["peak_absolute_closing_m_per_s"] = max(
                result["motion"]["peak_absolute_closing_m_per_s"],
                max((abs(c["closing_rate_m_per_s"]) for c in contacts), default=0.0),
            )
            result["motion"]["contact_force_sum_n_per_m"] += sum(c["force_n_per_m"] for c in contacts)
        phases = _phase_totals(world.phase_budget)
        for phase in phases.keys() | previous_phases.keys():
            for measure in ("hc_m3", "hm_m3"):
                delta = phases.get(phase, {}).get(measure, 0.0) - previous_phases.get(phase, {}).get(measure, 0.0)
                result["deformation"][phase][measure] += delta
                result["deformation_window"][phase][measure] += delta
        previous_phases = phases
        if index % checkpoint_steps == 0 or index == steps:
            row = _snapshot(world)
            row["hm_ledger_closure_error_m3"] = _hm_closure(row, initial_hm)
            row["topology_events_so_far"] = len(result["topology_events"])
            row["deformation_since_checkpoint"] = {
                phase: dict(values) for phase, values in result["deformation_window"].items()
            }
            result["deformation_window"] = defaultdict(lambda: defaultdict(float))
            result["checkpoints"].append(row)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            result["runtime_seconds"] = time.time() - started
            result["peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20)
            args.out.write_text(json.dumps(result, indent=1, default=dict) + "\n")
            print(json.dumps({"checkpoint": row, "runtime_seconds": result["runtime_seconds"]}), flush=True)
    motion = result["motion"]
    rate_samples = max(motion.pop("rate_samples"), 1)
    contact_samples = max(motion.pop("contact_samples"), 1)
    motion["mean_solved_rate_cm_per_yr"] = motion.pop("rate_sum_cm_per_yr") / rate_samples
    motion["mean_absolute_closing_m_per_s"] = motion.pop("absolute_closing_sum_m_per_s") / contact_samples
    motion["mean_contact_force_n_per_m"] = motion.pop("contact_force_sum_n_per_m") / contact_samples
    result["runtime_seconds"] = time.time() - started
    result["peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20)
    result["motion"] = dict(motion)
    result["deformation"] = {phase: dict(values) for phase, values in result["deformation"].items()}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({"out": str(args.out), "runtime_seconds": result["runtime_seconds"], "peak_rss_mib": result["peak_rss_mib"]}), flush=True)


if __name__ == "__main__":
    main()
