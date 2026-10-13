"""Shared plumbing for the save-replay measurement scripts (measure_ice_sheet.py,
measure_lake_silt.py): the common command line, loading a save against a chosen backend tree,
and starting the checkpoint loop.

`--backend` points at the backend to import, so one script measures both `origin/main` (a
detached worktree) and the branch. `load_replay` puts it first on `sys.path` before anything
from `app` is imported, so a script imports its own `app` modules only after calling it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def parse_replay_args(doc: str | None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
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


def start_checkpoints(args: argparse.Namespace, result: dict, first: dict) -> tuple[int, int]:
    """Record and print the starting checkpoint `first`, and create the output directory.
    Returns `(steps, steps_per_checkpoint)`."""
    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    result["checkpoints"].append(first)
    print(json.dumps(first), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    return steps, per_checkpoint
