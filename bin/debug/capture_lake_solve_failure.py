#!/usr/bin/env python3
"""Step a saved world and pickle the lake solve's inputs if it raises (#905).

A lake-solver failure deep inside a mature-world step is slow to reach (one
full step) but the solver itself is fast. This driver wraps
``lakes.solve_lakes_prognostic`` so the first call that raises writes its
exact arguments (forest, elevation, prior depth, climate, hypsometry, ...) to a
pickle, which ``--replay`` then re-runs in seconds without stepping the world.

Each step prints its wall time, the last lake solve's lake volume and its ledger residual,
so the same command run in a baseline worktree gives an A/B.

Run from the repository root:

    backend/.venv/bin/python bin/debug/capture_lake_solve_failure.py <save> \\
        --step-years 100000 --steps 3 --out <capture>.pkl
    backend/.venv/bin/python bin/debug/capture_lake_solve_failure.py \\
        --replay <capture>.pkl

Saves and captures stay out of the repo (AGENTS.md).
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app.hydroclimate import lakes  # noqa: E402


def _capture(save: Path, step_years: float, steps: int, out: Path) -> int:
    from app import persistence
    from app.world import step_world

    world = persistence.load_world_bytes(save.read_bytes())
    original = lakes.solve_lakes_prognostic
    captured = {}
    last = {}

    def wrapped(*args, **kwargs):
        snapshot = pickle.dumps((args, kwargs))
        try:
            result = original(*args, **kwargs)
        except Exception:
            if not captured:
                out.write_bytes(snapshot)
                captured["path"] = out
            raise
        last["volume_m3"] = float(result.lake_depth @ args[5].area_m2)
        last["residual_m3"] = float(result.ledger.residual_m3)
        return result

    lakes.solve_lakes_prognostic = wrapped
    for index in range(steps):
        start = time.time()
        try:
            step_world(world, step_years)
        except Exception:
            traceback.print_exc()
            print(f"step {index} failed at {world.elapsed_years / 1e6:.1f} Myr; capture: {captured.get('path')}")
            return 1
        print(
            f"step {index} ok -> {world.elapsed_years / 1e6:.1f} Myr ({time.time() - start:.1f} s) "
            f"lake volume {last.get('volume_m3', float('nan')):.9e} m^3, "
            f"lake ledger residual {last.get('residual_m3', float('nan')):.3e} m^3",
            flush=True,
        )
    return 0


def _replay(capture: Path) -> int:
    args, kwargs = pickle.loads(capture.read_bytes())
    start = time.time()
    result = lakes.solve_lakes_prognostic(*args, **kwargs)
    print(f"solve ok in {time.time() - start:.1f} s; lake volume {float(result.lake_depth @ args[5].area_m2):.6g} m^3")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", nargs="?", type=Path)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--replay", type=Path)
    opts = parser.parse_args()
    if opts.replay:
        return _replay(opts.replay)
    if opts.save is None or opts.out is None:
        parser.error("save and --out are required unless --replay is given")
    return _capture(opts.save, opts.step_years, opts.steps, opts.out)


if __name__ == "__main__":
    raise SystemExit(main())
