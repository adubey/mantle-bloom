"""Per-call cost of rebuilding plate outlines on a saved world (#808).

For the largest plates, times a cold `outline_world` (every cache dropped via `_reset_caches`,
as after `_replace_topology`) and a warm `outline_world` (loops cached, only the frame
rotation), and prints a digest of every plate's outline so a change can be shown to leave the
geometry bit-identical. Also prints the `get_neighbours` result for every plate.

Usage (from the repository root, with this worktree's own venv):

    backend/.venv/bin/python bin/debug/bench_plate_outline.py ~/Downloads/<save>.mbworld [--top 6] [--json OUT]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
from app import persistence  # noqa: E402


def _best(fn, repeats: int) -> float:
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn()
        times.append(time.perf_counter() - start)
    return min(times)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("save", type=Path)
    ap.add_argument("--top", type=int, default=6)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    with args.save.open("rb") as handle:
        world = persistence.load_world_file(handle)
    plates = list(world.plates)
    rows = []
    for plate in sorted(plates, key=lambda p: -len(p._keys))[: args.top]:
        def cold(plate=plate):
            plate._reset_caches()
            plate.outline_world()

        cold_s = _best(cold, args.repeats)
        warm_s = _best(plate.outline_world, args.repeats)
        rows.append({"plate": plate.plate_id, "cells": len(plate._keys), "cold_loops_s": cold_s, "warm_outline_s": warm_s})
        print(f"plate {plate.plate_id:4d} cells {len(plate._keys):7d} cold loops {cold_s * 1e3:8.1f} ms  warm outline {warm_s * 1e3:7.2f} ms")
    digest = hashlib.sha256()
    for plate in plates:
        digest.update(np.ascontiguousarray(plate.outline_world()).tobytes())
    neighbours = {p.plate_id: sorted(q.plate_id for q in p.get_neighbours(plates)) for p in plates}
    print("outline sha256", digest.hexdigest())
    if args.json:
        args.json.write_text(json.dumps({"rows": rows, "outline_sha256": digest.hexdigest(), "neighbours": neighbours}, indent=1))


if __name__ == "__main__":
    main()
