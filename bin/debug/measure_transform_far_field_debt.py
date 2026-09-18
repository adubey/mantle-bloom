#!/usr/bin/env python3
"""Investigation for GitHub issue #189's "still open, separate from this fix" note:
`lithosphere_plate.py`'s `deform()` adds `transform_uplift` (transform-boundary pressure-ridge
relief) and `far_field_uplift` (broad far-field collision relief) as bare elevation deltas with
no `crustal_thickness_m` backing -- by design ("local relief without net crustal shortening"),
unlike the fault-relief bug #191 fixed. The issue asks whether this unbacked-by-design drift is
actually a live problem worth bounding over long runs, the same shape of problem as #191's bug
just intentional rather than an oversight.

Same "debt" concept `drill_down_isostasy.py` uses per pinned node (`elevation -
isostatic_elevation(hc, hm, rho_c)`), but sampled across *all* land at several checkpoints of a
run several times longer than #189's original 120 My, not just nodes already pinned at
MAX_ELEVATION_M -- since after #191, transform_uplift/far_field_uplift are the only remaining
source of this debt for any line with real Hc/Hm tracking (every other elevation-moving path
either derives elevation from Hc/Hm directly or -- faults.py, volcanism.py -- bills its own
delta through lithosphere.back_elevation_gain), so debt growth over time is a direct read on
how much these two terms alone are injecting.

Usage: backend/.venv/bin/python bin/debug/measure_transform_far_field_debt.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import SWEEP_SEEDS  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import lithosphere  # noqa: E402
from app.elevation_lines import MAX_ELEVATION_M  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

NODE_DENSITY = 1.0
STEP_YEARS = 10_000_000
CHECKPOINTS_MY = (120, 240, 360, 500)
COLLISION_UPLIFT_MULTIPLIERS = (1.0, 0.25)
PIN_EPS_M = 1e-6


def run_one(seed: int, collision_uplift_multiplier: float) -> list[dict]:
    world = generate_world(seed=seed, node_density=NODE_DENSITY)
    world.collision_uplift_multiplier = collision_uplift_multiplier

    rows = []
    years_done = 0.0
    for checkpoint_my in CHECKPOINTS_MY:
        target = checkpoint_my * 1_000_000
        while years_done < target:
            step = min(STEP_YEARS, target - years_done)
            step_world(world, years=step)
            years_done += step

        for plate in world.plates:
            if not plate.lines:
                continue
            elevation = plate.collect("elevation")
            land = elevation > world.sea_level_m
            if not np.any(land):
                continue
            hc = plate.collect("crustal_thickness_m")[land]
            hm = plate.collect("mantle_lithosphere_thickness_m")[land]
            crust_type_code = plate.collect("crust_type_code")[land]
            rho_c = lithosphere.node_crust_density(crust_type_code, plate.crust_type)
            equilibrium = lithosphere.isostatic_elevation(hc, hm, rho_c)
            debt = elevation[land] - equilibrium
            pinned = elevation[land] >= MAX_ELEVATION_M - PIN_EPS_M

            rows.append({
                "seed": seed,
                "collision_uplift_multiplier": collision_uplift_multiplier,
                "years_my": years_done / 1_000_000,
                "plate_crust_type": plate.crust_type,
                "land_nodes": int(np.sum(land)),
                "pinned_nodes": int(np.sum(pinned)),
                "debt_mean_m": float(debt.mean()),
                "debt_p90_m": float(np.percentile(debt, 90)),
                "debt_max_m": float(debt.max()),
                "positive_debt_nodes": int(np.sum(debt > 100.0)),
            })

    return rows


def main() -> None:
    jobs = [(seed, m) for seed in SWEEP_SEEDS for m in COLLISION_UPLIFT_MULTIPLIERS]
    print(f"{len(jobs)} job(s) to run, checkpoints at {CHECKPOINTS_MY} My.", file=sys.stderr)

    all_rows: list[dict] = []
    t0 = time.perf_counter()
    for n, (seed, m) in enumerate(jobs, start=1):
        try:
            rows = run_one(seed, m)
        except Exception as exc:  # noqa: BLE001 -- keep the sweep going past one bad job
            print(f"FAILED seed={seed} x{m}: {exc!r}", file=sys.stderr)
            continue
        all_rows.extend(rows)
        elapsed = time.perf_counter() - t0
        print(f"[{n}/{len(jobs)}] seed={seed} x{m} done ({elapsed:.0f}s elapsed)", file=sys.stderr)

    out_dir = Path(__file__).resolve().parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "measure_transform_far_field_debt.jsonl"
    with out_path.open("w") as f:
        for row in all_rows:
            f.write(json.dumps(row) + "\n")
    print(f"\nWrote {len(all_rows)} row(s) to {out_path}", file=sys.stderr)

    if not all_rows:
        print("No land nodes found at all -- nothing to analyze.")
        return

    print("\nBy checkpoint (pooled across seeds/multipliers):")
    for cp in CHECKPOINTS_MY:
        cp_rows = [r for r in all_rows if r["years_my"] == cp]
        if not cp_rows:
            continue
        total_land = sum(r["land_nodes"] for r in cp_rows)
        total_pinned = sum(r["pinned_nodes"] for r in cp_rows)
        total_pos_debt = sum(r["positive_debt_nodes"] for r in cp_rows)
        max_debt = max(r["debt_max_m"] for r in cp_rows)
        print(
            f"  {cp:4d} My: land={total_land:6d} pinned={total_pinned:4d} "
            f"({100 * total_pinned / max(total_land, 1):.2f}%) "
            f"nodes w/ >100m debt={total_pos_debt:5d} max_debt={max_debt:8.0f}m"
        )


if __name__ == "__main__":
    main()
