"""Shared plumbing for the #318 collision-polarity diagnostics (diagnose_polarity_fallbacks.py,
diagnose_polarity_merge_conflicts.py): both replay a save with a `collision_polarity` function
wrapped to record what it decides, and write the records as one JSON file.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Callable
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence, world as world_mod  # noqa: E402

PROGRESS_EVERY_STEPS = 50


def load_polarity_replay(doc: str | None) -> tuple[argparse.Namespace, world_mod.World]:
    """Parse `save --myr --step-years --out` and load the save."""
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--myr", type=float, default=43.1)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    return args, persistence.load_world_bytes(args.save.read_bytes())


def replay(
    world: world_mod.World,
    args: argparse.Namespace,
    progress: Callable[[], str],
    report: Callable[[], dict],
    before_step: Callable[[int], None] | None = None,
) -> None:
    """Step `world` for `--myr`, printing `progress()` every 50 steps, then write
    `{"save": ..., **report()}` to `--out` when given. `before_step(k)` runs before step `k`."""
    steps = int(round(args.myr * 1e6 / args.step_years))
    for k in range(steps):
        if before_step is not None:
            before_step(k)
        world_mod.step_world(world, args.step_years)
        if (k + 1) % PROGRESS_EVERY_STEPS == 0:
            print(f"{world.elapsed_years / 1e6:.1f} Myr: {progress()}", flush=True)
    if args.out is not None:
        args.out.write_text(json.dumps({"save": str(args.save), **report()}, indent=1, default=float))
