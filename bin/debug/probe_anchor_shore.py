#!/usr/bin/env python3
"""Record how far a fixed anchor point sits from the nearest lake body at every step (#811).

`trace_lake_anchor.py` reports only the components that contain the anchor, so a dry anchor
gives no information. This probe keeps the nearest wet component whether or not it contains the
anchor, with its volume, area and mean depth, the anchor cell's own ground elevation and water
depth, and the great-circle distance from the anchor cell to the nearest wet cell of that body.
If the anchor is wet only when the body's level clears the anchor cell, the dry steps are shoreline
crossings, not a lake that disappears.

Run from the repository root:

    backend/.venv/bin/python bin/debug/probe_anchor_shore.py <save> \\
        --anchors analysis/rfc0016/target-anchor-window-a.json \\
        --myr 0.4 --step-years 5000 --out /tmp/anchor-shore.json
"""

from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args


def add_arguments(ap) -> None:
    """Add the anchor file to the shared replay CLI."""
    ap.add_argument("--anchors", type=Path, required=True, help="JSON file containing anchor_points")


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=811)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

anchor_points = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
world.ice_age_period_years = 0.0
visible_depth = hydrology.LAKE_MIN_VISIBLE_DEPTH_M


def snapshot(step: int, runtime_s: float = 0.0) -> dict:
    """Describe the wet body nearest the anchor cell, and the anchor cell itself."""
    fields = world.hydrology_cache
    points = np.asarray(fields.points)
    _, nearest = cKDTree(points).query(anchor_points, k=1)
    cell = int(nearest[0])
    depth = np.asarray(fields.lake_depth)
    wet = depth > visible_depth
    groups = hydrology.lake_components(wet, fields.neighbor_idx)
    hyps = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area = hyps.area_m2
    best = None
    for nodes in groups:
        nodes = np.asarray(nodes, dtype=np.int64)
        gap = np.linalg.norm(points[nodes] - points[cell], axis=1).min() * 6371.0
        if best is None or gap < best["gap_km"]:
            best = {
                "gap_km": float(gap),
                "contains_anchor": bool(cell in set(nodes.tolist())),
                "volume_km3": float(np.dot(depth[nodes], area[nodes]) / 1e9),
                "area_km2": float(area[nodes].sum() / 1e6),
                "mean_depth_m": float(np.average(depth[nodes], weights=area[nodes])),
                "max_depth_m": float(depth[nodes].max()),
                "nodes": int(len(nodes)),
            }
    return {
        "step": step,
        "elapsed_years": float(world.elapsed_years),
        "runtime_s": runtime_s,
        "anchor_cell": cell,
        "anchor_cell_elevation_m": float(fields.elevation[cell]),
        "anchor_cell_lake_depth_m": float(depth[cell]),
        "nearest_body": best,
        "world_lake_km3": float(np.dot(depth, area) / 1e9),
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
    }


args.out.parent.mkdir(parents=True, exist_ok=True)
result = {"header": header, "records": [snapshot(0)]}
steps = int(round(args.myr * 1e6 / args.step_years))
for step in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    result["records"].append(snapshot(step, time.perf_counter() - started))
    args.out.write_text(json.dumps(result, indent=1))
    if step % max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years))) == 0 or step == steps:
        print(json.dumps(result["records"][-1]), flush=True)

print(f"wrote {args.out}")
