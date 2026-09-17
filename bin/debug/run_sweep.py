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
    backend/.venv/bin/python bin/debug/run_sweep.py --node-density 2.0 --checkpoints 30,60,90,120,150,180 \
        --params avg_rotation_rate,volcanism,collision_uplift_amount,collision_uplift_distance \
        --out results/sweep_results_density2.jsonl

--node-density/--checkpoints changes are recorded on every output row (see sweep_lib.
run_one_job) and folded into the resume/dedup key below, specifically so a differently-configured
run always lands in (or is recognized as incomplete in) its own --out file rather than silently
averaging together with an incompatible run's rows that happen to share a (parameter, multiplier,
seed) -- point a different configuration at its own --out file rather than reusing one.
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
    NODE_DENSITY,
    PARAM_SPECS,
    SWEEP_SEEDS,
    run_one_job,
)

DEFAULT_OUT = Path(__file__).resolve().parent / "results" / "sweep_results.jsonl"


def _seen_checkpoints(out_path: Path, node_density: float) -> dict[tuple[str, float, int], set[float]]:
    """(parameter, multiplier, seed) -> the set of checkpoint_years already recorded for it at
    this exact `node_density` in `out_path`. Used both to decide which triples are complete
    (see `_completed_triples`) and, row-by-row, to skip re-writing a checkpoint that's already
    on disk -- a triple can be *incomplete* (so it gets re-run from scratch, since run_one_job
    has no partial-job resume finer than "the whole job") while still having *some* of its
    checkpoint rows already flushed, if run_sweep.py itself (not a worker -- worker failures are
    caught below and never write partial records) was killed mid-write between two of a
    finished job's `f.write` calls. Without this, resuming after that exact interruption reruns
    the triple's full checkpoint set and appends duplicate rows for the checkpoints that
    survived. A row from a differently-configured run (different node_density) is ignored
    entirely -- see this module's own docstring for why mixing configurations in one file is
    unsafe."""
    if not out_path.exists():
        return {}
    seen_checkpoints: dict[tuple[str, float, int], set[float]] = {}
    with out_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            row = json.loads(line)
            if row.get("node_density", NODE_DENSITY) != node_density:
                continue
            key = (row["parameter"], row["multiplier"], row["seed"])
            seen_checkpoints.setdefault(key, set()).add(row["checkpoint_years"])
    return seen_checkpoints


def _completed_triples(
    seen_checkpoints: dict[tuple[str, float, int], set[float]], checkpoint_years: list[int]
) -> set[tuple[str, float, int]]:
    """(parameter, multiplier, seed) triples that already have a record for every one of
    `checkpoint_years` -- a partially-written triple (e.g. a worker or the main process was
    killed mid-job) is treated as incomplete and re-run from scratch, since run_one_job only
    returns a job's records once it finishes every checkpoint."""
    needed = set(checkpoint_years)
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
    parser.add_argument(
        "--checkpoints", default=",".join(str(y // 1_000_000) for y in CHECKPOINT_YEARS), help="Myr, comma-separated"
    )
    parser.add_argument("--node-density", type=float, default=NODE_DENSITY)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    params = args.params.split(",")
    for p in params:
        if p not in PARAM_SPECS:
            parser.error(f"unknown parameter {p!r}; choices are {', '.join(PARAM_SPECS)}")
    multipliers = [float(m) for m in args.multipliers.split(",")]
    seeds = [int(s) for s in args.seeds.split(",")]
    checkpoint_years = [int(float(m) * 1_000_000) for m in args.checkpoints.split(",")]
    node_density = args.node_density

    args.out.parent.mkdir(parents=True, exist_ok=True)
    seen_checkpoints = _seen_checkpoints(args.out, node_density)
    done = _completed_triples(seen_checkpoints, checkpoint_years)
    all_jobs = build_jobs(params, multipliers, seeds)
    jobs = [j for j in all_jobs if j not in done]
    print(f"{len(jobs)} job(s) to run ({len(all_jobs) - len(jobs)} already complete in {args.out}).")
    if not jobs:
        return

    t0 = time.perf_counter()
    completed = 0
    with args.out.open("a") as f, ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(run_one_job, param, multiplier, seed, node_density, tuple(checkpoint_years)): (
                param,
                multiplier,
                seed,
            )
            for param, multiplier, seed in jobs
        }
        for future in as_completed(futures):
            param, multiplier, seed = futures[future]
            try:
                records = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the sweep going past one bad job
                print(f"FAILED {param} x{multiplier} seed={seed}: {exc!r}", file=sys.stderr)
                continue
            key = (param, multiplier, seed)
            already = seen_checkpoints.get(key, set())
            for record in records:
                if record["checkpoint_years"] in already:
                    continue  # survived a prior interrupted write of this same triple
                f.write(json.dumps(record) + "\n")
                already.add(record["checkpoint_years"])
            seen_checkpoints[key] = already
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
