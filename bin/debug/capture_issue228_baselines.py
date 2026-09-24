#!/usr/bin/env python3
"""Issue #228 Phase 0a: generate the early/mid-game baseline worlds.

Regenerates the issue #147 profile configuration (seed 0, node_density 4.0, climate_density
4.0, fluid_density 2.0, 100 kyr steps -- see bin/debug/profile_animation_60m.py and
analysis/issue147-profile-20260922/metadata.json) so the representation baselines line up
with the performance baseline the migration is measured against. Saves an .mbworld at every
checkpoint and appends a cheap per-interval time series of line fragmentation and
conservation totals (no climate compute, no rendering) to `seed<N>-timeseries.jsonl`.

The .mbworld checkpoints are large (~100-200 MB each) and are not checked in; the output
directory's own .gitignore excludes them. They are reproducible from this script, and
`characterize_plate_surface.py` records each file's sha256 alongside its results.

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/capture_issue228_baselines.py
    backend/.venv/bin/python bin/debug/capture_issue228_baselines.py --checkpoints-myr 0 10 --out /tmp/x
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import lithosphere, persistence  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

DEFAULT_OUT = REPO_ROOT / "analysis" / "issue228-phase0a" / "worlds"
STEP_YEARS = 100_000
SERIES_EVERY_STEPS = 10


def series_row(world) -> dict:
    """Cheap whole-world totals -- line fragmentation plus the extensive quantities whose
    drift the migration must not worsen. Node area is the nominal constant
    `lithosphere.node_area_m2` (see characterize_plate_surface.py for how far that nominal
    area drifts from the area the nodes actually cover)."""
    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    lines = [line for plate in world.plates for line in plate.lines if len(line) > 0]
    lengths = np.array([len(line) for line in lines], dtype=int)
    hc_cont = hc_all = hm_all = 0.0
    land = 0
    for plate in world.plates:
        hc = plate.collect("crustal_thickness_m")
        hc_all += float(hc.sum())
        hm_all += float(plate.collect("mantle_lithosphere_thickness_m").sum())
        if plate.crust_type == "continental":
            hc_cont += float(hc.sum())
        land += int(np.count_nonzero(plate.collect("elevation") > world.sea_level_m))
    nodes = int(lengths.sum())
    return {
        "elapsed_myr": world.elapsed_years / 1e6,
        "steps": world.steps_taken,
        "plates": len(world.plates),
        "nodes": nodes,
        "lines": int(len(lengths)),
        "one_node_lines": int(np.count_nonzero(lengths == 1)),
        "two_node_lines": int(np.count_nonzero(lengths == 2)),
        "mean_nodes_per_line": float(lengths.mean()) if len(lengths) else 0.0,
        "sea_level_m": float(world.sea_level_m),
        "land_fraction_node": land / nodes if nodes else None,
        "hc_volume_km3": hc_all * area_m2 / 1e9,
        "continental_hc_volume_km3": hc_cont * area_m2 / 1e9,
        "hm_volume_km3": hm_all * area_m2 / 1e9,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--checkpoints-myr", type=float, nargs="+", default=[0, 10, 30, 60])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    world = generate_world(seed=args.seed, node_density=4.0, climate_density=4.0, fluid_density=2.0)
    print(f"generated seed {args.seed} in {time.perf_counter() - started:.1f}s", flush=True)

    checkpoints = sorted(round(m * 1e6 / STEP_YEARS) for m in args.checkpoints_myr)
    series_path = args.out / f"seed{args.seed}-timeseries.jsonl"
    with series_path.open("w") as series:
        series.write(json.dumps(series_row(world)) + "\n")
        for target in checkpoints:
            while world.steps_taken < target:
                step_world(world, STEP_YEARS)
                if world.steps_taken % SERIES_EVERY_STEPS == 0:
                    series.write(json.dumps(series_row(world)) + "\n")
                    series.flush()
            path = args.out / f"seed{args.seed}-{world.steps_taken * STEP_YEARS / 1e6:05.1f}myr.mbworld"
            path.write_bytes(persistence.save_world_bytes(world))
            print(f"{path.name}: step {world.steps_taken}, {time.perf_counter() - started:.0f}s elapsed", flush=True)


if __name__ == "__main__":
    main()
