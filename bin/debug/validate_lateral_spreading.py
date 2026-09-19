#!/usr/bin/env python3
"""Before/after sweep-harness validation for the sea-level-biased delamination melt spread
(GitHub issue #176's "lateral spreading" direction), the acceptance test the issue's own
"Suggested next step" calls for: confirm land volume/ice-cap% growth flattens out at high
`collision_uplift_multiplier` values without changing the low-multiplier baseline.

Unlike run_sweep.py's PARAM_SPECS (each a continuous x-baseline multiplier on an existing
World field), the fix isn't its own World field -- it's gated by the module-level
`lithosphere_plate.LATERAL_SPREADING_BIAS_ENABLED` toggle (see that constant's own
comment), so this script flips that directly, same idiom issue #180's now-superseded
validate_crumpling.py used for CRUMPLE_TRANSFER_FRACTION. Reuses sweep_lib.run_one_job as-is
(it already threads `collision_uplift_amount`'s multiplier onto
`world.collision_uplift_multiplier` and returns compute_outcome_stats' four stats), so the
"on"/"off" pair for a given seed+multiplier differs in only the one toggle, from the same
checkout and the same generate_world seed.

Reports, per multiplier, the mean +/- std delta (on - off) of each stat, AND the low-vs-high
multiplier *gap* under each condition -- if the fix works, the 0.25x-vs-4x gap should shrink
under "on" relative to "off" while the 0.25x point itself barely moves.

Usage (run from anywhere; the venv with numpy/scipy/etc is backend/.venv):
    backend/.venv/bin/python bin/debug/validate_lateral_spreading.py
    backend/.venv/bin/python bin/debug/validate_lateral_spreading.py --seeds 829071382,579428537 --workers 4
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

from app import lithosphere_plate  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import NODE_DENSITY, SWEEP_SEEDS, run_one_job  # noqa: E402

MULTIPLIERS = (0.25, 4.0)  # issue #176's own "Run A" low/high ends
CHECKPOINT_YEARS = (120_000_000,)  # issue #176's own headline checkpoint

STAT_KEYS = ("land_fraction", "ice_cap_fraction", "plains_fraction", "land_volume_above_sea_km3")

DEFAULT_OUT = Path(__file__).resolve().parent / "results" / "lateral_spreading_validation.jsonl"


def run_one_toggle_job(seed: int, multiplier: float, bias_on: bool, node_density: float) -> list[dict]:
    """One (seed, multiplier, bias_on) job. Sets the toggle unconditionally at the top (not
    just when it differs from its current value) so a worker process reused across jobs by
    ProcessPoolExecutor never leaks a previous job's setting -- same reasoning as
    sweep_lib._apply_rotation_rate_overrides."""
    lithosphere_plate.LATERAL_SPREADING_BIAS_ENABLED = bias_on
    records = run_one_job(
        "collision_uplift_amount", multiplier, seed, node_density=node_density, checkpoint_years=CHECKPOINT_YEARS
    )
    for r in records:
        r["bias_on"] = bias_on
    return records


def _mean_std(values: list[float]) -> tuple[float | None, float | None]:
    values = [v for v in values if v is not None]
    if not values:
        return None, None
    mean = statistics.fmean(values)
    std = statistics.pstdev(values) if len(values) > 1 else 0.0
    return mean, std


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default=",".join(str(s) for s in SWEEP_SEEDS))
    parser.add_argument("--node-density", type=float, default=NODE_DENSITY)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    node_density = args.node_density

    jobs = [(seed, mult, bias_on) for seed in seeds for mult in MULTIPLIERS for bias_on in (False, True)]
    print(f"{len(jobs)} job(s) ({len(seeds)} seeds x {len(MULTIPLIERS)} multipliers x 2 conditions) at node_density={node_density}.")

    all_records: list[dict] = []
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_one_toggle_job, seed, mult, bias_on, node_density): (seed, mult, bias_on)
            for seed, mult, bias_on in jobs
        }
        completed = 0
        for future in as_completed(futures):
            seed, mult, bias_on = futures[future]
            try:
                records = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the run going past one bad job
                print(f"FAILED seed={seed} multiplier={mult} bias_on={bias_on}: {exc!r}", file=sys.stderr)
                continue
            all_records.extend(records)
            completed += 1
            print(
                f"[{completed}/{len(jobs)}] seed={seed} multiplier={mult} bias_on={bias_on} done "
                f"({time.perf_counter() - t0:.0f}s elapsed)"
            )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for r in all_records:
            f.write(json.dumps(r) + "\n")
    print(f"wrote {len(all_records)} records to {args.out}")

    # Pair each seed's "on"/"off" record at each multiplier and report the delta directly --
    # does the fix move the metric at all, in which direction, per multiplier.
    by_key: dict[tuple[int, float], dict[bool, dict]] = {}
    for r in all_records:
        by_key.setdefault((r["seed"], r["multiplier"]), {})[r["bias_on"]] = r

    print(f"\n{'multiplier':>10}  {'stat':<28} {'mean delta (on - off)':>24} {'std':>10}  n")
    for mult in MULTIPLIERS:
        for stat in STAT_KEYS:
            deltas = []
            for seed in seeds:
                pair = by_key.get((seed, mult))
                if not pair or False not in pair or True not in pair:
                    continue
                off_val, on_val = pair[False][stat], pair[True][stat]
                if off_val is not None and on_val is not None:
                    deltas.append(on_val - off_val)
            mean, std = _mean_std(deltas)
            mean_str = f"{mean:+.6g}" if mean is not None else "n/a"
            std_str = f"{std:.3g}" if std is not None else "n/a"
            print(f"{mult:>9}x  {stat:<28} {mean_str:>24} {std_str:>10}  {len(deltas)}")

    # The actual acceptance test (issue #176's own "Suggested next step"): does the
    # low-multiplier-to-high-multiplier GAP shrink under "on" relative to "off", per stat,
    # without the low end itself moving much.
    print(f"\n{'stat':<28} {'gap off (4x-0.25x)':>20} {'gap on (4x-0.25x)':>20} {'0.25x delta (on-off)':>22}")
    for stat in STAT_KEYS:
        low_off = [by_key[(s, MULTIPLIERS[0])][False][stat] for s in seeds if (s, MULTIPLIERS[0]) in by_key and False in by_key[(s, MULTIPLIERS[0])]]
        high_off = [by_key[(s, MULTIPLIERS[1])][False][stat] for s in seeds if (s, MULTIPLIERS[1]) in by_key and False in by_key[(s, MULTIPLIERS[1])]]
        low_on = [by_key[(s, MULTIPLIERS[0])][True][stat] for s in seeds if (s, MULTIPLIERS[0]) in by_key and True in by_key[(s, MULTIPLIERS[0])]]
        high_on = [by_key[(s, MULTIPLIERS[1])][True][stat] for s in seeds if (s, MULTIPLIERS[1]) in by_key and True in by_key[(s, MULTIPLIERS[1])]]
        gap_off, _ = _mean_std([h - l for h, l in zip(high_off, low_off)])
        gap_on, _ = _mean_std([h - l for h, l in zip(high_on, low_on)])
        low_delta, _ = _mean_std([on - off for on, off in zip(low_on, low_off)])
        gap_off_str = f"{gap_off:+.6g}" if gap_off is not None else "n/a"
        gap_on_str = f"{gap_on:+.6g}" if gap_on is not None else "n/a"
        low_delta_str = f"{low_delta:+.6g}" if low_delta is not None else "n/a"
        print(f"{stat:<28} {gap_off_str:>20} {gap_on_str:>20} {low_delta_str:>22}")


if __name__ == "__main__":
    main()
