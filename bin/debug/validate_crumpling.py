#!/usr/bin/env python3
"""Before/after sweep-harness validation for collision crumpling (GitHub issue #180),
the "not yet done" item called out on PR #183 and in #180's own prototype comment.

Unlike run_sweep.py's PARAM_SPECS (each a continuous x-baseline multiplier on an existing
World field), crumpling isn't wired to a World field at all -- it's gated entirely by the
module-level `lithosphere_plate.CRUMPLE_TRANSFER_FRACTION` constant, so this script toggles
that directly rather than adding it to PARAM_SPECS. Setting the fraction to 0.0 is an exact
behavioral no-op of `_redistribute_crumple_mass` (confirmed by reading it: with the fraction
zeroed, `transfer_hc`/`transfer_hm` are all-zero, so `donors`/`gain_hc`/`receiving` end up
empty and nothing downstream of that point runs) -- i.e. identical to the code as it stood
before #180's commit, on top of the *same* #176 relief-taper baseline used both ways. That's
a tighter isolation of crumpling's own marginal effect than diffing across commits would be,
since both conditions run from the exact same checkout and generate_world seed, differing in
only this one constant.

Runs sweep_lib.compute_outcome_stats's four outcome stats at every sweep_lib.SWEEP_SEEDS x
CHECKPOINT_YEARS combination, once with crumpling on (the shipped 0.3) and once off (0.0), and
reports the paired per-seed delta (on minus off) at each checkpoint -- directly answering the
PR's own open question: does crumpling move land%/ice-cap% in the intended direction, and
is the land% cost the PR worried about (transferring Hc off valley nodes that were visible as
land, onto ridge nodes disproportionately likely to already be capped) small relative to the
land% gain from unsaturating those ridges?

Usage (run from anywhere; the venv with numpy/scipy/etc is backend/.venv):
    backend/.venv/bin/python bin/debug/validate_crumpling.py
    backend/.venv/bin/python bin/debug/validate_crumpling.py --seeds 829071382,579428537 --workers 4
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
from app.world import generate_world, step_world  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import (  # noqa: E402
    CHECKPOINT_YEARS,
    NODE_DENSITY,
    SWEEP_SEEDS,
    compute_outcome_stats,
)

# Captured at import time (before any job mutates it) so "on" always means "whatever this
# checkout actually ships", not a hardcoded 0.3 that would silently go stale if the shipped
# constant is ever retuned.
CRUMPLING_ON_FRACTION = lithosphere_plate.CRUMPLE_TRANSFER_FRACTION
CRUMPLING_OFF_FRACTION = 0.0

STAT_KEYS = ("land_fraction", "ice_cap_fraction", "plains_fraction", "land_volume_above_sea_km3")

DEFAULT_OUT = Path(__file__).resolve().parent / "results" / "crumpling_validation.jsonl"


def run_one_job(seed: int, crumpling_on: bool, node_density: float, checkpoint_years: tuple[int, ...]) -> list[dict]:
    """One (seed, crumpling_on) job, run_sweep.py-style: fresh `generate_world`, stepped out to
    `checkpoint_years`'s last entry, one outcome-stats record per checkpoint. Sets the module
    constant unconditionally at the top (not just when toggling it away from its current value)
    so a worker process reused across jobs by ProcessPoolExecutor never leaks a previous job's
    setting -- same reasoning as sweep_lib._apply_rotation_rate_overrides."""
    lithosphere_plate.CRUMPLE_TRANSFER_FRACTION = CRUMPLING_ON_FRACTION if crumpling_on else CRUMPLING_OFF_FRACTION

    world = generate_world(seed=seed, node_density=node_density)

    records = []
    years_done = 0.0
    for checkpoint in checkpoint_years:
        step_world(world, years=checkpoint - years_done)
        years_done = checkpoint
        record = {
            "seed": seed,
            "crumpling_on": crumpling_on,
            "crumple_transfer_fraction": lithosphere_plate.CRUMPLE_TRANSFER_FRACTION,
            "checkpoint_years": checkpoint,
            "node_density": node_density,
        }
        record.update(compute_outcome_stats(world))
        records.append(record)
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
    parser.add_argument(
        "--checkpoints", default=",".join(str(y // 1_000_000) for y in CHECKPOINT_YEARS), help="Myr, comma-separated"
    )
    parser.add_argument("--node-density", type=float, default=NODE_DENSITY)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    seeds = [int(s) for s in args.seeds.split(",")]
    checkpoint_years = tuple(int(float(m) * 1_000_000) for m in args.checkpoints.split(","))
    node_density = args.node_density

    jobs = [(seed, crumpling_on) for seed in seeds for crumpling_on in (False, True)]
    print(f"{len(jobs)} job(s) ({len(seeds)} seeds x 2 conditions) at node_density={node_density}.")

    all_records: list[dict] = []
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one_job, seed, on, node_density, checkpoint_years): (seed, on) for seed, on in jobs}
        completed = 0
        for future in as_completed(futures):
            seed, on = futures[future]
            try:
                records = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the run going past one bad job
                print(f"FAILED seed={seed} crumpling_on={on}: {exc!r}", file=sys.stderr)
                continue
            all_records.extend(records)
            completed += 1
            print(f"[{completed}/{len(jobs)}] seed={seed} crumpling_on={on} done ({time.perf_counter() - t0:.0f}s elapsed)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as f:
        for r in all_records:
            f.write(json.dumps(r) + "\n")
    print(f"wrote {len(all_records)} records to {args.out}")

    # Pair each seed's "on" and "off" record at each checkpoint and report the delta directly,
    # rather than just each condition's own mean -- the question this script exists to answer
    # is the *effect size* of the toggle, not either condition's absolute level.
    by_key: dict[tuple[int, int], dict[bool, dict]] = {}
    for r in all_records:
        by_key.setdefault((r["seed"], r["checkpoint_years"]), {})[r["crumpling_on"]] = r

    print(f"\n{'checkpoint':>10}  {'stat':<28} {'mean delta (on - off)':>24} {'std':>10}  n")
    for checkpoint in checkpoint_years:
        for stat in STAT_KEYS:
            deltas = []
            for seed in seeds:
                pair = by_key.get((seed, checkpoint))
                if not pair or False not in pair or True not in pair:
                    continue
                off_val, on_val = pair[False][stat], pair[True][stat]
                if off_val is not None and on_val is not None:
                    deltas.append(on_val - off_val)
            mean, std = _mean_std(deltas)
            mean_str = f"{mean:+.6g}" if mean is not None else "n/a"
            std_str = f"{std:.3g}" if std is not None else "n/a"
            print(f"{checkpoint // 1_000_000:>8}My  {stat:<28} {mean_str:>24} {std_str:>10}  {len(deltas)}")


if __name__ == "__main__":
    main()
