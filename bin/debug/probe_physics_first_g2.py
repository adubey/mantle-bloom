#!/usr/bin/env python3
"""Record one generated world's G2 quasi-static motion diagnostics.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_physics_first_g2.py \
        --seed 326494039 --steps 10 --step-myr 2 --node-density 1 --out /tmp/g2.json

The probe advances complete world steps and records the torque residual, blocked-contact
forces, yielded contacts, plate-rate distribution, and total step time. It is a pilot
diagnostic, not the matched origin/main comparison or the G2 acceptance replay.
"""

from __future__ import annotations

import argparse
import inspect
import json
import statistics
import subprocess
import time
from pathlib import Path

import numpy as np

from app.elevation_lines import line_spacing_rad
from app import mantle
from app.world import generate_world, step_world


def _commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except subprocess.CalledProcessError:
        return "unknown"


def _geology_snapshot(world) -> dict:
    areas = []
    land = []
    hc = []
    hm = []
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        area = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
        areas.append(area)
        land.append(plate.collect("elevation") > 0.0)
        hc.append(plate.collect("crustal_thickness_m"))
        hm.append(plate.collect("mantle_lithosphere_thickness_m"))
    if not areas:
        return {"land_fraction": 0.0, "mean_hc_m": 0.0, "mean_hm_m": 0.0}
    area = np.concatenate(areas)
    land_mask = np.concatenate(land)
    crust = np.concatenate(hc)
    mantle = np.concatenate(hm)
    total_area = float(area.sum())
    return {
        "land_fraction": float(area[land_mask].sum() / total_area),
        "mean_hc_m": float((crust * area).sum() / total_area),
        "mean_hm_m": float((mantle * area).sum() / total_area),
        "plate_count": len(world.plates),
    }


def run(seed: int, steps: int, step_myr: float, node_density: float, snapshot_every: int) -> dict:
    world = generate_world(seed=seed, node_density=node_density)
    motion_rows: list[dict] = []
    runtimes = []
    geology_rows = [{"step": 0, "elapsed_myr": 0.0, **_geology_snapshot(world)}]
    records_motion_diagnostics = "motion_diagnostics" in inspect.signature(step_world).parameters
    for step in range(1, steps + 1):
        started = time.perf_counter()
        if records_motion_diagnostics:
            step_world(world, step_myr * 1_000_000.0, motion_diagnostics=motion_rows)
        else:
            step_world(world, step_myr * 1_000_000.0)
        runtimes.append(time.perf_counter() - started)
        if records_motion_diagnostics:
            motion_rows[-1]["step_seconds"] = runtimes[-1]
        else:
            rates = [mantle.rad_per_yr_to_cm_per_yr(float(np.linalg.norm(p.omega))) for p in world.plates]
            motion_rows.append({
                "step": step,
                "elapsed_years": (step - 1) * step_myr * 1_000_000.0,
                "plates": len(rates),
                "max_plate_rate_cm_per_yr": max(rates, default=0.0),
                "median_plate_rate_cm_per_yr": statistics.median(rates) if rates else 0.0,
                "step_seconds": runtimes[-1],
            })
        if step % snapshot_every == 0 or step == steps:
            geology_rows.append({
                "step": step,
                "elapsed_myr": step * step_myr,
                **_geology_snapshot(world),
            })
    return {
        "commit": _commit(),
        "seed": seed,
        "steps": steps,
        "step_myr": step_myr,
        "node_density": node_density,
        "elapsed_myr": steps * step_myr,
        "median_step_seconds": statistics.median(runtimes) if runtimes else 0.0,
        "motion": motion_rows,
        "geology": geology_rows,
        "ledgers": world.craton_ledger,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--step-myr", type=float, default=2.0)
    parser.add_argument("--node-density", type=float, default=1.0)
    parser.add_argument("--snapshot-every", type=int, default=25)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.steps < 1 or args.step_myr <= 0 or args.node_density <= 0 or args.snapshot_every < 1:
        parser.error("steps, step-myr, node-density, and snapshot-every must be positive")
    result = run(args.seed, args.steps, args.step_myr, args.node_density, args.snapshot_every)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(serialized)
    else:
        print(serialized, end="")


if __name__ == "__main__":
    main()
