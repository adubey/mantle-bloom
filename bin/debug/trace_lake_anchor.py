#!/usr/bin/env python3
"""Trace lake components around fixed geographic anchors during a world replay.

Run from the repository root with a saved world and JSON containing an ``anchor_points``
array of unit-sphere XYZ coordinates:

    backend/.venv/bin/python bin/debug/trace_lake_anchor.py <save> \\
        --anchors analysis/issue675/shallow-lake-anchor.json \\
        --myr 0.8 --step-years 10000 --checkpoint-myr 0.01 \\
        --out /tmp/lake-anchor-trace.json

The anchors identify a location, not a persistent lake identity. If they are split across
components, all unique components containing a nearest anchor cell are reported.
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
    """Add anchor-specific arguments to the shared replay CLI."""
    ap.add_argument("--anchors", type=Path, required=True, help="JSON file containing anchor_points")
    ap.add_argument(
        "--anchor-tolerance-km", type=float, default=20.0,
        help="maximum nearest-cell distance accepted for an anchor (default: 20 km)",
    )


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=675)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

anchor_data = json.loads(args.anchors.read_text())
anchor_points = np.asarray(anchor_data["anchor_points"], dtype=float)
if anchor_points.ndim != 2 or anchor_points.shape[1] != 3 or len(anchor_points) == 0:
    raise ValueError("anchor_points must be a nonempty array of XYZ unit-sphere coordinates")
if args.anchor_tolerance_km <= 0.0:
    raise ValueError("anchor tolerance must be positive")

world.ice_age_period_years = 0.0
visible_depth = hydrology.LAKE_MIN_VISIBLE_DEPTH_M


def snapshot(step: int, runtime_s: float = 0.0) -> dict:
    """Record connected water bodies containing the nearest cells to fixed anchors."""
    fields = world.hydrology_cache
    points = np.asarray(fields.points)
    distance, nearest = cKDTree(points).query(anchor_points, k=1)
    distance_km = distance * 6371.0
    close = distance_km <= args.anchor_tolerance_km

    wet = np.asarray(fields.lake_depth) > visible_depth
    components = hydrology.lake_components(wet, fields.neighbor_idx)
    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area_m2 = hypsometry.area_m2
    labels = np.full(len(wet), -1, dtype=np.int32)
    component_data = []
    for component_id, nodes in enumerate(components):
        nodes = np.asarray(nodes, dtype=np.int64)
        labels[nodes] = component_id
        component_data.append(
            {
                "component_id": component_id,
                "area_km2": float(area_m2[nodes].sum() / 1e6),
                "volume_km3": float(np.dot(fields.lake_depth[nodes], area_m2[nodes]) / 1e9),
                "mean_depth_m": float(np.average(fields.lake_depth[nodes], weights=area_m2[nodes])),
            }
        )

    anchor_labels = labels[nearest]
    matched_ids = sorted({int(label) for label, is_close in zip(anchor_labels, close) if is_close and label >= 0})
    matched = [component_data[component_id] for component_id in matched_ids]
    ledger = getattr(fields, "lake_water_ledger", None)
    normalized_residual = None
    if ledger is not None:
        normalized_residual = abs(float(ledger.residual_m3)) / max(float(ledger.external_supply_m3), 1.0)

    return {
        "step": step,
        "elapsed_years": float(world.elapsed_years),
        "runtime_s": runtime_s,
        "anchor_count": int(len(anchor_points)),
        "anchor_points_wet": int(np.count_nonzero((anchor_labels >= 0) & close)),
        "anchor_points_far": int(np.count_nonzero(~close)),
        "nearest_cell_distance_km_median": float(np.median(distance_km)),
        "nearest_cell_distance_km_max": float(np.max(distance_km)),
        "matched_components": matched,
        "anchor_component_area_km2": float(sum(item["area_km2"] for item in matched)),
        "anchor_component_volume_km3": float(sum(item["volume_km3"] for item in matched)),
        "global_area_km2": float(area_m2[wet].sum() / 1e6),
        "global_volume_km3": float(np.dot(fields.lake_depth, area_m2) / 1e9),
        "merge_events": sum(event.kind == "merge" for event in (getattr(fields, "lake_events", []) or [])),
        "split_events": sum(event.kind == "split" for event in (getattr(fields, "lake_events", []) or [])),
        "normalized_ledger_residual": normalized_residual,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
    }


args.out.parent.mkdir(parents=True, exist_ok=True)
result = {"header": header, "anchors": anchor_data, "records": [snapshot(0)]}
steps = int(round(args.myr * 1e6 / args.step_years))
for step in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    record = snapshot(step, time.perf_counter() - started)
    result["records"].append(record)
    args.out.write_text(json.dumps(result, indent=1))
    if step % max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years))) == 0 or step == steps:
        print(json.dumps(record), flush=True)

print(f"wrote {args.out}")
