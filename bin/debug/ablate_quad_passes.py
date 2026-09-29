"""Issue #228 Phase 4: ablate this branch's quad-pass changes one at a time.

Each variant steps a quad world from the same seed with one change reverted by monkeypatch
and reports crustal volume, land fraction, plate count, exact uncovered/overlapping sphere and
seconds per step, every `--every` steps:

- `all`: every change on;
- `no-overlap`: overlap tracking back on node proximity;
- `no-gapfill`: gap filling back on node distance, with the MIN_GAP_NODES floor on growth;
- `no-vents`: mostly-magmatic new cells don't ignite.

    cd backend && .venv/bin/python ../bin/debug/ablate_quad_passes.py --seed 7 --steps 150 --variant all
"""

from __future__ import annotations

import argparse
import contextlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app import gaps, plates, quad_tectonics  # noqa: E402
from app.sparse_quad_patch import PlateWithSparseQuadPatch  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402
from measure_quad_passes import coverage, fibonacci_sphere, totals  # noqa: E402


@contextlib.contextmanager
def inexact():
    PlateWithSparseQuadPatch.territory_is_exact = False
    try:
        yield
    finally:
        PlateWithSparseQuadPatch.territory_is_exact = True


def without_exact(fn):
    def wrapped(*args, **kwargs):
        with inexact():
            return fn(*args, **kwargs)
    return wrapped


def apply_variant(variant: str) -> None:
    if variant == "no-overlap":
        plates.compute_node_overlap = without_exact(plates.compute_node_overlap)
    elif variant == "no-gapfill":
        gaps.fill_gaps_by_growing_neighbours = without_exact(gaps.fill_gaps_by_growing_neighbours)
        gaps.reconcile_gap_tracks = without_exact(gaps.reconcile_gap_tracks)
        original = quad_tectonics.grow_frontier
        quad_tectonics.grow_frontier = lambda *a, **k: original(*a, **{**k, "standoff": True})
    elif variant == "no-vents":
        quad_tectonics._ignite_early_rift_volcanoes = lambda *a, **k: None
    elif variant != "all":
        raise SystemExit(f"unknown variant {variant}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--steps", type=int, default=150)
    parser.add_argument("--every", type=int, default=50)
    parser.add_argument("--variant", default="all")
    args = parser.parse_args()

    apply_variant(args.variant)
    sample = fibonacci_sphere(100_000)
    world = generate_world(seed=args.seed, node_density=1.0, surface="quad")
    seconds = 0.0
    for step in range(1, args.steps + 1):
        started = time.perf_counter()
        step_world(world, 1_000_000.0)
        seconds += time.perf_counter() - started
        if step % args.every == 0:
            row = {"variant": args.variant, "seed": args.seed, "step": step, **totals(world), **coverage(world, sample), "s_per_step": round(seconds / step, 3)}
            print(json.dumps(row), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
