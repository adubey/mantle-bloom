"""Shared plumbing for the save-replay measurement scripts (measure_ice_sheet.py,
measure_lake_silt.py): the common command line, loading a save against a chosen backend tree,
and the step-and-checkpoint loop.

`--backend` points at the backend to import, so one script measures both `origin/main` (a
detached worktree) and the branch. `load_replay` puts it first on `sys.path` before anything
from `app` is imported, so a script imports its own `app` modules only after calling it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from collections.abc import Callable
from pathlib import Path


def parse_replay_args(
    doc: str | None, add_arguments: Callable[[argparse.ArgumentParser], None] | None = None
) -> argparse.Namespace:
    """The common command line; `add_arguments` adds a script's own options."""
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    if add_arguments is not None:
        add_arguments(ap)
    ap.add_argument("save", type=Path)
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--myr", type=float, required=True)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--checkpoint-myr", type=float, default=1.0)
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


def load_replay(args: argparse.Namespace, issue: int):
    """Import `app` from `args.backend` and load the save. Returns `(world, result)`, where
    `result` is the output JSON's header (backend, its commit, the save and its SHA-256) with an
    empty `checkpoints` list."""
    backend = args.backend.resolve()
    sys.path.insert(0, str(backend))
    from app import persistence  # noqa: E402

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=backend, check=True, capture_output=True, text=True).stdout.strip()
    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    result = {
        "issue": issue,
        "backend": str(backend),
        "commit": commit,
        "save": str(args.save),
        "save_sha256": hashlib.sha256(data).hexdigest(),
        "step_years": args.step_years,
        "checkpoints": [],
    }
    return world, result


def run_replay(
    world,
    args: argparse.Namespace,
    result: dict,
    snapshot: Callable[[], dict],
    after_step: Callable[[float], None],
    close_window: Callable[[int], dict],
) -> None:
    """Step `world` for `--myr` in `--step-years` steps, writing `result` to `--out` at every
    checkpoint. `snapshot()` measures the world's state. `after_step(seconds)` runs after each
    step with that step's wall time, to add to the script's window. `close_window(steps)` returns
    the window's measures for a checkpoint `steps` steps after the last one, and starts a new
    window."""
    from app import world as world_mod  # noqa: E402

    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result["checkpoints"].append(snapshot())
    print(json.dumps(result["checkpoints"][-1]), flush=True)
    for i in range(1, steps + 1):
        started = time.perf_counter()
        world_mod.step_world(world, args.step_years)
        after_step(time.perf_counter() - started)
        if i % per_checkpoint == 0 or i == steps:
            snap = snapshot()
            snap.update(close_window(per_checkpoint))
            snap["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
            result["checkpoints"].append(snap)
            print(json.dumps(snap), flush=True)
            args.out.write_text(json.dumps(result, indent=1))
    args.out.write_text(json.dumps(result, indent=1))
    print(f"wrote {args.out}")
