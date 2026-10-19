#!/usr/bin/env python3
"""Run one G2 replay with a selected craton-formation waiting time.

Run each seed/threshold pair as a separate process. The geometric/Hc eligibility gate is
unchanged; only ``CRATON_FORMATION_MYR`` varies. Example:

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/g2_formation_gate_sensitivity.py \
      --seed 326494039 --threshold-myr 300 --steps 400 --step-myr 1 \
      --out /tmp/g2-gates-326494039-300myr.jsonl

This is a timer sensitivity experiment, not a test of whether the current definition of
"quiet" measures low physical deformation.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from app import cratons
from app.elevation_lines import line_spacing_rad
from app.world import generate_world, step_world


def _weighted_mean(values: np.ndarray, areas: np.ndarray, mask: np.ndarray) -> float | None:
    if not np.any(mask):
        return None
    return float(np.dot(values[mask], areas[mask]) / areas[mask].sum())


def _snapshot(world, step: int, threshold_myr: float, years_myr: float) -> dict:
    spacing = line_spacing_rad(world.node_density)
    hops = cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, spacing)
    areas = []
    hcs = []
    hms = []
    stable = []
    eligible = []
    live_craton = []
    genuine = []
    for plate in world.plates:
        if not plate.node_count():
            continue
        area = cratons._areas(world, plate)
        hc = plate.collect("crustal_thickness_m")
        craton = plate.collect("craton_crust_m")
        host = cratons._genuine_continental(plate)
        interior = cratons._interior_hops(plate, host, hops)
        if interior is None:
            interior = np.zeros(plate.node_count(), dtype=int)
        areas.append(area)
        hcs.append(hc)
        hms.append(plate.collect("mantle_lithosphere_thickness_m"))
        stable.append(plate.collect("stable_continental_myr"))
        genuine.append(host)
        eligible.append(host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & (interior >= hops) & (craton <= 0.0))
        live_craton.append(craton > 0.0)

    if not areas:
        return {}
    area = np.concatenate(areas)
    hc = np.concatenate(hcs)
    hm = np.concatenate(hms)
    stable_clock = np.concatenate(stable)
    eligible_pool = np.concatenate(eligible)
    live_craton_nodes = np.concatenate(live_craton)
    genuine_host = np.concatenate(genuine)
    near_maturity = eligible_pool & (stable_clock >= max(0.0, threshold_myr - years_myr))
    craton_volume_m3 = float(sum(
        np.dot(plate.collect("craton_crust_m"), cratons._areas(world, plate))
        for plate in world.plates if plate.node_count()
    ))
    formed_m3 = float(world.craton_ledger.get("formed_m3", 0.0))
    return {
        "step": step,
        "elapsed_myr": world.elapsed_years / 1e6,
        "threshold_myr": threshold_myr,
        "plate_count": len(world.plates),
        "formed_total_mkm3": formed_m3 / 1e15,
        "live_craton_mkm3": craton_volume_m3 / 1e15,
        "genuine_host_km2": float(area[genuine_host & ~live_craton_nodes].sum() / 1e6),
        "eligible_uncratonized_km2": float(area[eligible_pool].sum() / 1e6),
        "eligible_within_one_step_of_maturity_km2": float(area[near_maturity].sum() / 1e6),
        "eligible_uncratonized_hc_mean_km": _weighted_mean(hc / 1000.0, area, eligible_pool),
        "eligible_uncratonized_hm_mean_km": _weighted_mean(hm / 1000.0, area, eligible_pool),
        "eligible_clock_mean_myr": _weighted_mean(stable_clock, area, eligible_pool),
    }


def run(seed: int, threshold_myr: float, steps: int, step_myr: float,
        snapshot_every: int, out: Path) -> None:
    world = generate_world(seed=seed, node_density=1.0)
    original_threshold = cratons.CRATON_FORMATION_MYR
    original_update = cratons.update
    cratons.CRATON_FORMATION_MYR = threshold_myr
    dt_years = step_myr * 1e6
    rows = []
    started = time.perf_counter()

    def wrapped_update(current_world, years):
        before = float(current_world.craton_ledger.get("formed_m3", 0.0))
        original_update(current_world, years)
        formed_this_step_m3 = float(current_world.craton_ledger.get("formed_m3", 0.0)) - before
        step = int(round(current_world.elapsed_years / years)) if years else 0
        if step % snapshot_every == 0 or formed_this_step_m3 > 0.0 or step == 1 or step == steps:
            row = _snapshot(current_world, step, threshold_myr, step_myr)
            row["formed_this_step_mkm3"] = formed_this_step_m3 / 1e15
            rows.append(row)
            with out.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")

    cratons.update = wrapped_update
    try:
        for index in range(steps):
            step_world(world, dt_years)
            if (index + 1) % 25 == 0:
                print(
                    f"seed={seed} threshold={threshold_myr:g} Myr "
                    f"step={index + 1}/{steps} myr={(index + 1) * step_myr:g} "
                    f"runtime_min={(time.perf_counter() - started) / 60:.1f}",
                    flush=True,
                )
    finally:
        cratons.update = original_update
        cratons.CRATON_FORMATION_MYR = original_threshold

    summary = {
        "done": True,
        "seed": seed,
        "threshold_myr": threshold_myr,
        "steps": steps,
        "step_myr": step_myr,
        "elapsed_myr": steps * step_myr,
        "runtime_min": (time.perf_counter() - started) / 60,
        "checkpoints": len(rows),
        "ledger": world.craton_ledger,
    }
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--threshold-myr", type=float, choices=(100.0, 200.0, 300.0), required=True)
    parser.add_argument("--steps", type=int, default=400)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--snapshot-every", type=int, default=10)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.steps < 1 or args.step_myr <= 0 or args.snapshot_every < 1:
        parser.error("steps, step-myr, and snapshot-every must be positive")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    run(args.seed, args.threshold_myr, args.steps, args.step_myr, args.snapshot_every, args.out)


if __name__ == "__main__":
    main()
