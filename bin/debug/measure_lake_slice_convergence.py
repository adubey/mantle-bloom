#!/usr/bin/env python3
"""Measure lake-solver slice convergence from a saved checkpoint (#853).

Two arms separate the lake solver's internal slice length from the climate
refresh interval, starting from the same checkpoint:

* ``fixed``: one ``step_world`` of ``--total-years``, with the lake solver
  forced (monkeypatch, this script only) to ``substeps = total / slice``.
  Climate, routing and every other process see the same single step; only
  the lake slice changes.
* ``refreshed``: ``step_world`` repeated in steps of ``--slice-years`` over
  the same total, so climate and routing are recomputed every slice.

Run from the repository root, one arm per invocation:

    backend/.venv/bin/python bin/debug/measure_lake_slice_convergence.py <save> \\
        --anchors analysis/rfc0016/target-anchor-window-d.json \\
        --mode fixed --total-years 25000 --slice-years 200 \\
        --out ~/mantle-bloom-captures/issue853/fixed-25000-slice200.json

The snapshot reports the wet component under the anchor (volume, flooded
area, upstream catchment area), the world lake volume, and the lake
solver's own step ledger (evaporation, sublimation, overflow out, internal
transfers, merge/split events, normalized residual), so a slice effect can
be attributed to a ledger term rather than asserted. Merged inland-sea
state follows #811: anchor catchment above 15M km2.
"""

from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args, upstream_catchment

MERGED_CATCHMENT_KM2 = 15.0e6


def add_arguments(ap) -> None:
    ap.add_argument("--anchors", type=Path, required=True, help="JSON file containing anchor_points")
    ap.add_argument("--mode", choices=["fixed", "refreshed"], required=True)
    ap.add_argument("--total-years", type=float, required=True, help="total interval, years")
    ap.add_argument("--slice-years", type=float, required=True, help="lake slice length, years")
    ap.add_argument("--record-every", type=int, default=1, help="refreshed: record every N steps")


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=853)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

world.ice_age_period_years = 0.0
anchor_points = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
total_years = float(args.total_years)
slice_years = float(args.slice_years)
if slice_years <= 0.0 or total_years <= 0.0 or total_years / slice_years < 1.0:
    raise ValueError("need 0 < slice-years <= total-years")

# The last lake solve's ledger, captured in fixed mode by the wrapper below.
last_solve: dict = {}


def _ledger_record(ledger, events) -> dict:
    return {
        "external_supply_km3": ledger.external_supply_m3 / 1e9,
        "evaporation_km3": ledger.evaporation_m3 / 1e9,
        "sublimation_km3": ledger.sublimation_m3 / 1e9,
        "overflow_out_km3": ledger.overflow_out_m3 / 1e9,
        "overflow_to_ocean_km3": ledger.overflow_to_ocean_m3 / 1e9,
        "internal_transfer_km3": ledger.internal_transfer_m3 / 1e9,
        "storage_change_km3": ledger.storage_change_m3 / 1e9,
        "normalized_residual": abs(ledger.residual_m3) / max(ledger.external_supply_m3, 1.0),
        "merge_events": sum(event.kind == "merge" for event in events),
        "split_events": sum(event.kind == "split" for event in events),
    }


def _wrap_solver(original):
    def forced(forest, elevation, prev_depth, years, climate, hypsometry, **kwargs):
        substeps = max(1, int(round(float(years) / slice_years)))
        kwargs["substeps"] = substeps
        result = original(forest, elevation, prev_depth, years, climate, hypsometry, **kwargs)
        last_solve.clear()
        last_solve.update({"substeps": substeps, **_ledger_record(result.ledger, result.events)})
        return result

    return forced


if args.mode == "fixed":
    lakes.solve_lakes_prognostic = _wrap_solver(lakes.solve_lakes_prognostic)


tracked_members: np.ndarray | None = None


def snapshot(step: int, elapsed: float, runtime_s: float) -> dict:
    global tracked_members
    fields = world.hydrology_cache
    points = np.asarray(fields.points)
    _, nearest = cKDTree(points).query(anchor_points, k=1)
    anchor = int(np.atleast_1d(nearest)[0])
    wet = np.asarray(fields.lake_depth) > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    areas = hypsometry.area_m2
    record: dict = {
        "step": step,
        "elapsed_years": elapsed,
        "runtime_s": runtime_s,
        "anchor_wet": bool(wet[anchor]),
        "world_lake_km3": float(np.dot(fields.lake_depth, areas) / 1e9),
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
    }
    components = [np.asarray(c, dtype=np.int64) for c in hydrology.lake_components(wet, fields.neighbor_idx)]
    members = None
    how = None
    for component in components:
        if anchor in set(component.tolist()):
            members, how = component, "anchor"
            break
    if members is None and tracked_members is not None and len(components):
        # The anchor lake as a body, tracked by largest node overlap with the
        # previous snapshot's body (the anchor is a location, not an identity).
        tracked_set = set(tracked_members.tolist())
        overlaps = [len(tracked_set.intersection(set(c.tolist()))) for c in components]
        best = int(np.argmax(overlaps)) if overlaps else -1
        if best >= 0 and overlaps[best] > 0:
            members, how = components[best], "overlap"
    if members is None and len(components):
        wet_nodes = np.flatnonzero(wet)
        nearest_wet = wet_nodes[cKDTree(points[wet_nodes]).query(points[anchor], k=1)[1]]
        members = next(c for c in components if int(nearest_wet) in set(c.tolist()))
        how = "nearest"
    if members is not None:
        tracked_members = members
        catchment = upstream_catchment(
            np.asarray(fields.flow_target), np.asarray(fields.is_ocean), members
        )
        record.update({
            "tracked_how": how,
            "anchor_volume_km3": float(np.dot(fields.lake_depth[members], areas[members]) / 1e9),
            "anchor_area_km2": float(areas[members].sum() / 1e6),
            "anchor_catchment_km2": float(areas[catchment].sum() / 1e6),
            "merged_state": bool(areas[catchment].sum() / 1e6 > MERGED_CATCHMENT_KM2),
        })
    if last_solve:
        record["lake_solve"] = dict(last_solve)
    else:
        ledger = getattr(fields, "lake_water_ledger", None)
        if ledger is not None:
            record["lake_solve"] = _ledger_record(ledger, getattr(fields, "lake_events", []) or [])
    return record


result = {
    "header": header,
    "mode": args.mode,
    "total_years": total_years,
    "slice_years": slice_years,
    "records": [snapshot(0, 0.0, 0.0)],
}
args.out.parent.mkdir(parents=True, exist_ok=True)

if args.mode == "fixed":
    started = time.perf_counter()
    world_mod.step_world(world, total_years)
    record = snapshot(1, total_years, time.perf_counter() - started)
    result["records"].append(record)
    print(json.dumps(record), flush=True)
else:
    steps = int(round(total_years / slice_years))
    for step in range(1, steps + 1):
        started = time.perf_counter()
        world_mod.step_world(world, slice_years)
        if step % args.record_every == 0 or step == steps:
            record = snapshot(step, step * slice_years, time.perf_counter() - started)
            result["records"].append(record)
            print(json.dumps(record), flush=True)
            args.out.write_text(json.dumps(result, indent=1))

args.out.write_text(json.dumps(result, indent=1))
print(f"wrote {args.out}")
