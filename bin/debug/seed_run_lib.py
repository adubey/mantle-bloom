"""Shared plumbing for the measurement scripts that generate a fresh world per seed and step it
(measure_breaching.py, measure_cap_overflow.py, measure_collision_denudation.py,
measure_ice_loading.py, measure_meltwater_faults.py, measure_orogenic_relief.py,
measure_seasonal_melt.py): the common command line and the step-and-snapshot loop.

`--backend` points at the backend to import, so one script measures both `origin/main` (a
detached worktree) and a branch. `import_backend` puts it first on `sys.path`, so a script
imports its `app` modules only after calling it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections.abc import Callable
from pathlib import Path


def seed_run_parser(
    doc: str | None, *, steps: int = 20, years: float = 1_000_000.0, every: int = 10, every_help: str = "snapshot every Nth step"
) -> argparse.ArgumentParser:
    """The `--backend --seeds --steps --years --every` options. Scripts add their own on top."""
    parser = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=steps)
    parser.add_argument("--years", type=float, default=years)
    parser.add_argument("--every", type=int, default=every, help=every_help)
    return parser


def import_backend(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(args.backend.resolve()))


def run_seeds(
    args: argparse.Namespace,
    measure: Callable[[object, int, int], dict],
    start_seed: Callable[[object, int], None] | None = None,
) -> list[dict]:
    """Generate a world per `--seeds` and step it `--steps` times. Every `--every` steps and at
    the last, record `{seed, step, seconds, **measure(world, seed, step)}`, print it to 4
    significant figures, and rewrite `--out` (when given) with every row so far.
    `start_seed(world, seed)` runs on each new world before its first step."""
    from app.world import generate_world, step_world  # noqa: E402

    rows: list[dict] = []
    for seed in args.seeds:
        world = generate_world(seed=seed)
        if start_seed is not None:
            start_seed(world, seed)
        start = time.perf_counter()
        for step in range(1, args.steps + 1):
            step_world(world, args.years)
            if step % args.every and step != args.steps:
                continue
            row = {"seed": seed, "step": step, "seconds": round(time.perf_counter() - start, 1), **measure(world, seed, step)}
            rows.append(row)
            print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
            if args.out:
                args.out.parent.mkdir(parents=True, exist_ok=True)
                args.out.write_text(json.dumps(rows, indent=2))
    return rows
