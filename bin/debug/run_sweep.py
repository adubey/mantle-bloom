#!/usr/bin/env python3
"""Parameter-sensitivity sweep harness for the "land percent trending low / ice-cap percent
trending high" investigation (plate-collision-driven land loss and runaway plateau glaciation --
see GitHub issue #171). Runs bin/debug/sweep_lib.run_one_job
across every (parameter, multiplier, seed) combination and appends one JSON line per
(job, checkpoint) to --out, so the run can be split across multiple invocations (e.g. one
parameter at a time) and the harness script itself committed independently of any particular
run's results -- see sweep_lib.py's own module docstring for the seeds/checkpoints/step size
this uses and why.

Idempotent/resumable: before building the job list, reads whatever --out already has and skips
any (parameter, multiplier, seed) triple that already has a record for every checkpoint, so
re-running after a partial/interrupted sweep (or adding one more --params later) only computes
what's missing rather than duplicating rows.

Usage (run from anywhere; the venv with numpy/scipy/etc is backend/.venv):
    backend/.venv/bin/python bin/debug/run_sweep.py
    backend/.venv/bin/python bin/debug/run_sweep.py --params rain_erosion,volcanism --workers 4
    backend/.venv/bin/python bin/debug/run_sweep.py --seeds 829071382 --multipliers 1.0
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
    BASELINE_PARAM,
    CHECKPOINT_YEARS,
    MULTIPLIERS,
    PARAM_SPECS,
    SWEEP_SEEDS,
    run_one_job,
)

DEFAULT_OUT = Path(__file__).resolve().parent / "results" / "sweep_results.jsonl"


def _completed_triples(out_path: Path) -> set[tuple[str, float, int]]:
    """(parameter, multiplier, seed) triples that already have a record for every checkpoint in
    `out_path` -- a partially-written triple (e.g. the process was killed mid-job) is treated as
    incomplete and re-run from scratch, since run_one_job only returns a job's records once it
    finishes every checkpoint (there's no partial-job resume finer than "the whole job")."""
    if not out_path.exists():
        return set()
    seen_checkpoints: dict[tuple[str, float, int], set[float]] = {}
    with out_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            key = (row["parameter"], row["multiplier"], row["seed"])
            seen_checkpoints.setdefault(key, set()).add(row["checkpoint_years"])
    needed = set(CHECKPOINT_YEARS)
    return {key for key, checkpoints in seen_checkpoints.items() if needed <= checkpoints}


def build_jobs(params: list[str], multipliers: list[float], seeds: list[int]) -> list[tuple[str, float, int]]:
    """One job per (parameter, multiplier, seed), except multiplier == 1.0 which collapses to a
    single shared BASELINE_PARAM job per seed (see sweep_lib.BASELINE_PARAM) -- every parameter's
    "no change" point is the same run, so it's only computed once regardless of how many `params`
    include 1.0 in their own `multipliers`."""
    jobs: set[tuple[str, float, int]] = set()
    for param in params:
        for m in multipliers:
            for seed in seeds:
                jobs.add((BASELINE_PARAM, 1.0, seed) if m == 1.0 else (param, m, seed))
    return sorted(jobs)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--params", default=",".join(PARAM_SPECS), help=f"comma-separated subset of: {', '.join(PARAM_SPECS)}"
    )
    parser.add_argument("--multipliers", default=",".join(str(m) for m in MULTIPLIERS))
    parser.add_argument("--seeds", default=",".join(str(s) for s in SWEEP_SEEDS))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    params = args.params.split(",")
    for p in params:
        if p not in PARAM_SPECS:
            parser.error(f"unknown parameter {p!r}; choices are {', '.join(PARAM_SPECS)}")
    multipliers = [float(m) for m in args.multipliers.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    done = _completed_triples(args.out)
    jobs = [j for j in build_jobs(params, multipliers, seeds) if j not in done]
    skipped = len(build_jobs(params, multipliers, seeds)) - len(jobs)
    print(f"{len(jobs)} job(s) to run ({skipped} already complete in {args.out}).")
    if not jobs:
        return

    t0 = time.perf_counter()
    completed = 0
    with args.out.open("a") as f, ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one_job, param, multiplier, seed): (param, multiplier, seed) for param, multiplier, seed in jobs}
        for future in as_completed(futures):
            param, multiplier, seed = futures[future]
            try:
                records = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the sweep going past one bad job
                print(f"FAILED {param} x{multiplier} seed={seed}: {exc!r}", file=sys.stderr)
                continue
            for record in records:
                f.write(json.dumps(record) + "\n")
            f.flush()
            completed += 1
            elapsed = time.perf_counter() - t0
            rate = completed / elapsed
            remaining = (len(jobs) - completed) / rate if rate > 0 else float("nan")
            print(
                f"[{completed}/{len(jobs)}] {param} x{multiplier} seed={seed} done "
                f"({elapsed:.0f}s elapsed, ~{remaining:.0f}s remaining)"
            )


if __name__ == "__main__":
    main()
