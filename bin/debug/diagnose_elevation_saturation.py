#!/usr/bin/env python3
"""Issue #180 steps 1-2: how widespread/fast is MAX_ELEVATION_M land saturation at the
*baseline* (1.0x, i.e. untouched) collision_uplift_multiplier, across the full sweep_lib.
SWEEP_SEEDS set -- and for nodes pinned at the ceiling, is the underlying crustal thickness
(Hc) also maxed out at lithosphere.MAX_CRUSTAL_THICKNESS_M (a genuine "ran out of room" case
in apply_convergent_deformation's own overflow/delamination path), or does elevation clip
*before* Hc/Hm reach their own ceilings (pointing at the isostasy formula/reference
calibration instead)?

Reuses sweep_lib's exact fresh-generate_world / node_density=1.0 / STEP_YEARS stepping
convention (see its own module docstring) so results are directly comparable to the existing
sweep_results.jsonl, but sidesteps run_one_job/compute_outcome_stats -- those report
aggregate land_fraction/ice_cap_fraction/etc, not per-node Hc/Hm, which is exactly what step 2
needs. `land` here means `elevation > world.sea_level_m`, the same raw definition
stats._total_land_area_and_continental_volume uses (see its own docstring) --
not the climate-grid-resampled land_fraction sweep_lib.compute_outcome_stats reports.

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/diagnose_elevation_saturation.py
    backend/.venv/bin/python bin/debug/diagnose_elevation_saturation.py --seeds 829071382 --workers 1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import (  # noqa: E402
    CHECKPOINT_YEARS,
    NODE_DENSITY,
    STEP_YEARS,
    SWEEP_SEEDS,
    _apply_rotation_rate_overrides,
    BASELINE_AVG_RATE_CM_YR,
    BASELINE_MAX_RATE_CM_YR,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import lithosphere  # noqa: E402
from app.elevation_lines import MAX_ELEVATION_M, effective_is_continental_from_codes  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

# A node's elevation is produced by np.clip(..., MAX_ELEVATION_M) in isostatic_elevation --
# an exact float64 clip, not an accumulated sum, so a pinned node reads back bit-identical to
# the constant rather than merely "close". No tolerance needed for the pin test itself; a
# small one is still used for the "Hc/Hm at their own ceiling" tests below since those values
# arrive via ordinary arithmetic (multiplicative thickening, np.minimum chains), not a bare
# clip against the literal ceiling constant.
NEAR_CAP_FRAC = 0.99


def _seed_report(
    seed: int, node_density: float, checkpoint_years: tuple[int, ...], collision_uplift_multiplier: float = 1.0
) -> list[dict]:
    """One record per checkpoint for this seed. `collision_uplift_multiplier` defaults to 1.0
    (the true baseline, no PARAM_SPECS override at all -- every multiplier left at its
    untouched default, same as sweep_lib.BASELINE_PARAM); pass a different value to reproduce
    the specific point on #171/#180's own collision_uplift_amount sweep."""
    _apply_rotation_rate_overrides(BASELINE_AVG_RATE_CM_YR, BASELINE_MAX_RATE_CM_YR)
    world = generate_world(seed=seed, node_density=node_density)
    world.collision_uplift_multiplier = collision_uplift_multiplier

    records = []
    years_done = 0.0
    for checkpoint in checkpoint_years:
        while years_done < checkpoint:
            step = min(STEP_YEARS, checkpoint - years_done)
            step_world(world, years=step)
            years_done += step

        elevation_chunks, hc_chunks, hm_chunks, continental_chunks = [], [], [], []
        for plate in world.plates:
            elevation = plate.collect("elevation")
            if elevation.size == 0:
                continue
            hc = plate.collect("crustal_thickness_m")
            hm = plate.collect("mantle_lithosphere_thickness_m")
            codes = plate.collect("crust_type_code")
            is_continental = effective_is_continental_from_codes(codes, plate.crust_type == "continental")
            elevation_chunks.append(elevation)
            hc_chunks.append(hc)
            hm_chunks.append(hm)
            continental_chunks.append(is_continental)

        elevation = np.concatenate(elevation_chunks)
        hc = np.concatenate(hc_chunks)
        hm = np.concatenate(hm_chunks)
        is_continental = np.concatenate(continental_chunks)

        land = elevation > world.sea_level_m
        n_land = int(np.count_nonzero(land))
        pinned = land & (elevation >= MAX_ELEVATION_M)
        n_pinned = int(np.count_nonzero(pinned))

        record = {
            "seed": seed,
            "checkpoint_years": checkpoint,
            "node_density": node_density,
            "collision_uplift_multiplier": collision_uplift_multiplier,
            "n_land": n_land,
            "mean_land_elevation_m": float(np.mean(elevation[land])) if n_land else None,
            "pct_land_pinned": (n_pinned / n_land) if n_land else None,
            "pct_pinned_continental": (float(np.mean(is_continental[pinned])) if n_pinned else None),
        }

        if n_pinned:
            hc_pinned = hc[pinned]
            hm_pinned = hm[pinned]
            hc_ratio = hc_pinned / lithosphere.MAX_CRUSTAL_THICKNESS_M
            hm_ratio = hm_pinned / lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
            record.update(
                {
                    "pinned_hc_ratio_mean": float(np.mean(hc_ratio)),
                    "pinned_hc_ratio_p10": float(np.percentile(hc_ratio, 10)),
                    "pinned_hc_ratio_p50": float(np.percentile(hc_ratio, 50)),
                    "pinned_hm_ratio_mean": float(np.mean(hm_ratio)),
                    "pinned_hm_ratio_p50": float(np.percentile(hm_ratio, 50)),
                    # "Hc itself is also maxed out" -- the genuine ran-out-of-room case.
                    "pct_pinned_with_hc_near_cap": float(np.mean(hc_ratio >= NEAR_CAP_FRAC)),
                    # elevation clipped well before Hc got anywhere near its own ceiling --
                    # points at isostasy calibration rather than the thickening rate/cap.
                    "pct_pinned_with_hc_below_half_cap": float(np.mean(hc_ratio < 0.5)),
                }
            )
        else:
            record.update(
                {
                    "pinned_hc_ratio_mean": None,
                    "pinned_hc_ratio_p10": None,
                    "pinned_hc_ratio_p50": None,
                    "pinned_hm_ratio_mean": None,
                    "pinned_hm_ratio_p50": None,
                    "pct_pinned_with_hc_near_cap": None,
                    "pct_pinned_with_hc_below_half_cap": None,
                }
            )

        records.append(record)
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default=",".join(str(s) for s in SWEEP_SEEDS))
    parser.add_argument(
        "--checkpoints", default=",".join(str(y // 1_000_000) for y in CHECKPOINT_YEARS), help="Myr, comma-separated"
    )
    parser.add_argument("--node-density", type=float, default=NODE_DENSITY)
    parser.add_argument("--collision-uplift-multiplier", type=float, default=1.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results" / "elevation_saturation.jsonl")
    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    checkpoint_years = tuple(int(float(m) * 1_000_000) for m in args.checkpoints.split(","))
    node_density = args.node_density
    collision_uplift_multiplier = args.collision_uplift_multiplier

    args.out.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    completed = 0
    with args.out.open("w") as f, ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_seed_report, seed, node_density, checkpoint_years, collision_uplift_multiplier): seed
            for seed in seeds
        }
        for future in as_completed(futures):
            seed = futures[future]
            try:
                records = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the run going past one bad seed
                print(f"FAILED seed={seed}: {exc!r}", file=sys.stderr)
                continue
            for record in records:
                f.write(json.dumps(record) + "\n")
            f.flush()
            completed += 1
            elapsed = time.perf_counter() - t0
            rate = completed / elapsed
            remaining = (len(seeds) - completed) / rate if rate > 0 else float("nan")
            print(f"[{completed}/{len(seeds)}] seed={seed} done ({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)")

    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
