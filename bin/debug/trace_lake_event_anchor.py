#!/usr/bin/env python3
"""Trace a fixed geographic component and close its persistent-node water budget.

This diagnostic opts into transient per-node source/loss and routed spill-transaction traces.
It initializes a geographic lineage at ``--anchor-step`` (from fixed anchors or the component
nearest a registered volume), follows all overlapping descendants, and separates persistent
node storage changes from storage reclassified when connectivity changes.

Run from the repository root:

    backend/.venv/bin/python bin/debug/trace_lake_event_anchor.py <save> \\
        --anchor-step 4 --reference-volume-km3 136570 --myr 0.5 \\
        --step-years 100000 --out /tmp/seed-event.json
"""

from __future__ import annotations

import json
import os
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args


def add_arguments(ap) -> None:
    ap.add_argument("--anchor-step", type=int, required=True, help="step at which to initialize target lineage")
    ap.add_argument("--anchors", type=Path, help="JSON with fixed XYZ anchor_points")
    ap.add_argument("--reference-volume-km3", type=float, help="target volume used for automatic anchor selection")
    ap.add_argument("--anchor-count", type=int, default=8)


args = parse_replay_args(__doc__, add_arguments)
if bool(args.anchors) == (args.reference_volume_km3 is not None):
    raise ValueError("provide either --anchors or --reference-volume-km3")
os.environ["MANTLE_BLOOM_LAKE_LINEAGE_TRACE"] = "1"
world, header = load_replay(args, issue=675)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

if args.step_years != 100_000:
    raise ValueError("G6 uses the registered 100 kyr schedule")
if args.anchor_step < 1 or args.anchor_count < 1:
    raise ValueError("anchor step and anchor count must be positive")

world.ice_age_period_years = 0.0
visible_depth = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
anchor_points = None
records = []
lineage_ids = None
lineage_lake_ids = None
lineage_wet_ids = None
lineage_component_ids = None
static_anchor_points = (
    np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
    if args.anchors else None
)
if static_anchor_points is not None and (
    static_anchor_points.ndim != 2 or static_anchor_points.shape[1] != 3 or len(static_anchor_points) == 0
):
    raise ValueError("anchor_points must be a nonempty array of XYZ coordinates")


def components(fields):
    wet = np.asarray(fields.lake_depth) > visible_depth
    groups = hydrology.lake_components(wet, fields.neighbor_idx)
    hyps = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area = hyps.area_m2
    return [
        {
            "nodes": np.asarray(nodes, dtype=np.int64),
            "area_km2": float(area[nodes].sum() / 1e6),
            "volume_km3": float(np.dot(fields.lake_depth[nodes], area[nodes]) / 1e9),
        }
        for nodes in groups
    ], area, wet


def snapshot(step: int, runtime_s: float = 0.0) -> dict:
    fields = world.hydrology_cache
    bodies, area, wet = components(fields)
    matched = []
    if anchor_points is not None:
        distance, nearest = cKDTree(fields.points).query(anchor_points, k=1)
        labels = np.full(len(wet), -1, dtype=np.int32)
        for idx, body in enumerate(bodies):
            labels[body["nodes"]] = idx
        anchor_labels = labels[nearest]
        ids = sorted({int(label) for label in anchor_labels if label >= 0})
        matched = []
        for i in ids:
            item = {k: v for k, v in bodies[i].items() if k != "nodes"}
            item["anchor_count"] = int(np.count_nonzero(anchor_labels == i))
            matched.append(item)
        matched_distance_km = (np.asarray(distance) * 6371.0).tolist()
    else:
        matched_distance_km = []
    ledger = getattr(fields, "lake_water_ledger", None)
    ledger_data = None
    if ledger is not None:
        ledger_data = {
            name: float(getattr(ledger, name))
            for name in (
                "external_supply_m3", "wet_precip_m3", "dry_runoff_m3", "meltwater_m3",
                "evaporation_m3", "sublimation_m3", "overflow_out_m3", "storage_change_m3",
                "internal_transfer_m3", "freezing_m3", "melting_m3", "ice_storage_change_m3",
                "residual_m3",
            )
        }
        ledger_data["normalized_residual"] = abs(ledger.residual_m3) / max(ledger.external_supply_m3, 1.0)
    return {
        "step": step,
        "elapsed_years": float(world.elapsed_years),
        "runtime_s": runtime_s,
        "anchor_component": matched,
        "anchor_nearest_distance_km": matched_distance_km,
        "global_liquid_volume_km3": float(np.dot(fields.lake_depth, area) / 1e9),
        "global_lake_area_km2": float(area[wet].sum() / 1e6),
        "component_count": len(bodies),
        "ledger": ledger_data,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
    }


def component_ids_for_points(fields, bodies, points):
    distance, nearest = cKDTree(fields.points).query(points, k=1)
    labels = np.full(len(fields.lake_depth), -1, dtype=np.int32)
    for idx, body in enumerate(bodies):
        labels[body["nodes"]] = idx
    ids = sorted({int(label) for label in labels[nearest] if label >= 0})
    return ids, distance


def id_tuples(node_ids):
    return [tuple(map(int, node_id)) for node_id in np.asarray(node_ids, dtype=np.uint64)]


def storage_snapshot(fields):
    """Copy prior storage arrays before stepping can update their owning world state."""
    depth = np.asarray(fields.lake_depth)
    hypsometry = getattr(fields, "hypsometry", None)
    if hypsometry is None:
        hypsometry = lakes.CellHypsometry.whole_cells(fields.elevation)
    node_ids = getattr(fields, "lake_trace_node_ids", None)
    if node_ids is None:
        node_ids = np.zeros((len(depth), 2), dtype=np.uint64)
    return SimpleNamespace(
        lake_depth=np.array(depth, copy=True),
        lake_ice_m=np.array(getattr(fields, "lake_ice_m", np.zeros_like(depth)), copy=True),
        lake_trace_node_ids=np.array(node_ids, copy=True),
        hypsometry=lakes.CellHypsometry(
            np.array(hypsometry.elevation, copy=True),
            np.array(hypsometry.area_m2, copy=True),
            np.array(hypsometry.relief_m, copy=True),
            float(hypsometry.reference_area_m2),
        ),
    )


def active_basin_members(fields):
    active = lakes.governing_bodies(
        fields.lake_forest, fields.elevation, fields.lake_depth, fields.hypsometry
    )
    return [
        (
            int(body.lake_id),
            set(id_tuples(fields.lake_trace_node_ids[body.members])),
            set(id_tuples(fields.lake_trace_node_ids[body.members][fields.lake_depth[body.members] > visible_depth])),
        )
        for body in active
    ]


def lineage_budget(before, after, old_ids, new_ids, old_component_ids, new_component_ids, old_lake_ids, new_lake_ids):
    """Close persistent surface-node budgets across topology and array reordering."""
    before_ids = id_tuples(before.lake_trace_node_ids)
    after_ids = id_tuples(after.lake_trace_node_ids)
    if len(before_ids) != len(set(before_ids)) or len(after_ids) != len(set(after_ids)):
        raise RuntimeError("lake lineage trace node IDs must be unique across all plates")
    before_index = {node_id: idx for idx, node_id in enumerate(before_ids)}
    after_index = {node_id: idx for idx, node_id in enumerate(after_ids)}
    persistent_ids = old_ids & new_ids & before_index.keys() & after_index.keys()
    old_only_ids = old_ids - persistent_ids
    new_only_ids = new_ids - persistent_ids
    after_lineage = np.asarray([after_index[key] for key in new_ids if key in after_index], dtype=np.int64)
    old_liquid = np.asarray(before.lake_depth) * before.hypsometry.area_m2
    new_liquid = np.asarray(after.lake_depth) * after.hypsometry.area_m2
    old_ice = np.asarray(before.lake_ice_m) * before.hypsometry.area_m2
    new_ice = np.asarray(after.lake_ice_m) * after.hypsometry.area_m2
    old_stock, new_stock = old_liquid + old_ice, new_liquid + new_ice
    old_stock_current_area = np.zeros(len(after.lake_depth))
    old_liquid_current_area = np.zeros(len(after.lake_depth))
    old_ice_current_area = np.zeros(len(after.lake_depth))
    for node_id, before_idx in before_index.items():
        after_idx = after_index.get(node_id)
        if after_idx is not None:
            area = after.hypsometry.area_m2[after_idx]
            old_liquid_current_area[after_idx] = before.lake_depth[before_idx] * area
            old_ice_current_area[after_idx] = before.lake_ice_m[before_idx] * area
            old_stock_current_area[after_idx] = old_liquid_current_area[after_idx] + old_ice_current_area[after_idx]

    def old_storage(node_ids):
        retained = sum(
            old_stock_current_area[after_index[node_id]]
            for node_id in node_ids if node_id in after_index
        )
        removed = sum(
            old_stock[before_index[node_id]]
            for node_id in node_ids if node_id not in after_index and node_id in before_index
        )
        return float(retained + removed)

    after_component = np.asarray([after_index[key] for key in new_component_ids if key in after_index], dtype=np.int64)
    component_old_storage = old_storage(old_component_ids)
    component_new_storage = float(new_stock[after_component].sum())
    basin_old_storage = old_storage(old_ids)
    basin_new_storage = float(new_stock[after_lineage].sum())
    component_delta = component_new_storage - component_old_storage
    basin_membership_delta = basin_new_storage - basin_old_storage
    # Compare storage on the same endpoint basin boundary used for this step's source and spill
    # attribution. This includes cells newly assigned to the basin, but uses their prior stock
    # too; using the previous membership for old storage would mislabel removed-cell stock as
    # physical loss when the hierarchy changes.
    basin_delta = float(
        new_stock[after_lineage].sum() - old_storage(new_ids)
    )
    physical_delta = basin_delta
    connectivity = component_delta - basin_delta
    trace = after.lake_solve_trace
    ice_trace = after.lake_ice_trace
    if trace is None or ice_trace is None:
        raise RuntimeError("lake lineage trace was not recorded")
    node_ids = id_tuples(after.lake_trace_node_ids)
    basin_after = np.zeros(len(after.lake_depth), dtype=bool)
    basin_after[after_lineage] = True
    supply = float(trace.external_supply_m3[basin_after].sum())
    wet_precip = float(trace.wet_precip_m3[basin_after].sum())
    dry_runoff = float(trace.dry_runoff_m3[basin_after].sum())
    meltwater = float(trace.meltwater_m3[basin_after].sum())
    evaporation = float(trace.evaporation_m3[basin_after].sum())
    liquid_sublimation = float(trace.liquid_sublimation_m3[basin_after].sum())
    ice_sublimation = float(ice_trace["sublimation_by_node_m3"][basin_after].sum())
    transfer_net = 0.0
    transfer_edges = []
    transfer_in_by_node = np.zeros(len(after.lake_depth))
    transfer_out_by_node = np.zeros(len(after.lake_depth))

    def boundary_transfer_record(edge, direction):
        source_idx = int(edge.source_node_idx)
        target_idx = int(edge.target_node_idx)
        target_is_ocean = (
            0 <= target_idx < len(after.is_ocean) and bool(after.is_ocean[target_idx])
        )
        if target_is_ocean:
            destination_type = "ocean"
        elif edge.receiver_lake_id is not None:
            destination_type = "lake"
        elif target_idx >= 0:
            destination_type = "non_lake_sink"
        else:
            destination_type = "orphan_export"
        return {
            "source_lake_id": edge.source_lake_id,
            "receiver_lake_id": edge.receiver_lake_id,
            "source_node_id": node_ids[source_idx] if 0 <= source_idx < len(node_ids) else None,
            "target_node_id": node_ids[target_idx] if 0 <= target_idx < len(node_ids) else None,
            "target_is_ocean": target_is_ocean,
            "destination_type": destination_type,
            "volume_m3": edge.volume_m3,
            "direction": direction,
        }

    for edge in trace.transfers:
        if 0 <= edge.source_node_idx < len(after.lake_depth):
            transfer_out_by_node[edge.source_node_idx] += edge.volume_m3
        if edge.receiver_lake_id is not None and 0 <= edge.target_node_idx < len(after.lake_depth):
            transfer_in_by_node[edge.target_node_idx] += edge.volume_m3
        source_inside = (
            0 <= edge.source_node_idx < len(node_ids)
            and node_ids[edge.source_node_idx] in new_ids
        )
        target_inside = (
            edge.receiver_lake_id is not None
            and 0 <= edge.target_node_idx < len(node_ids)
            and node_ids[edge.target_node_idx] in new_ids
        )
        if source_inside and not target_inside:
            transfer_net -= edge.volume_m3
            transfer_edges.append(boundary_transfer_record(edge, "out"))
        elif target_inside and not source_inside:
            transfer_net += edge.volume_m3
            transfer_edges.append(boundary_transfer_record(edge, "in"))
    expected_physical = supply - evaporation - liquid_sublimation - ice_sublimation + transfer_net
    residual = physical_delta - expected_physical
    traced_overflow = sum(edge.volume_m3 for edge in trace.transfers if not edge.internal)
    global_storage_change = float(new_stock.sum() - old_storage(set(before_index)))
    global_trace_residual = (
        float(trace.external_supply_m3.sum())
        - float(trace.evaporation_m3.sum())
        - float(trace.liquid_sublimation_m3.sum())
        - float(ice_trace["sublimation_by_node_m3"].sum())
        - traced_overflow
        - global_storage_change
    )
    local_residual_rows = []
    for node_id in new_ids:
        idx = after_index.get(node_id)
        if idx is None:
            continue
        old_idx = before_index.get(node_id)
        storage_delta = float(new_stock[idx] - old_stock_current_area[idx]) if old_idx is not None else float(new_stock[idx])
        source = float(trace.external_supply_m3[idx])
        losses = float(
            trace.evaporation_m3[idx] + trace.liquid_sublimation_m3[idx]
            + ice_trace["sublimation_by_node_m3"][idx]
        )
        outflow = float(transfer_out_by_node[idx])
        inflow = float(transfer_in_by_node[idx])
        node_residual = storage_delta - source + losses + outflow - inflow
        if abs(node_residual) > 1.0:
            local_residual_rows.append({
                "node_id": list(node_id), "residual_m3": node_residual,
                "storage_change_m3": storage_delta, "supply_m3": source,
                "losses_m3": losses, "spill_out_m3": outflow, "spill_in_m3": inflow,
            })
    top_local_residuals = sorted(local_residual_rows, key=lambda row: -abs(row["residual_m3"]))[:10]
    return {
        "anchor_component_storage_change_m3": component_delta,
        "tracked_basin_storage_change_m3": basin_delta,
        "anchor_component_old_storage_m3": component_old_storage,
        "anchor_component_new_storage_m3": component_new_storage,
        "tracked_basin_old_storage_m3": basin_old_storage,
        "tracked_basin_new_storage_m3": basin_new_storage,
        "tracked_basin_membership_storage_change_m3": basin_membership_delta,
        "liquid_storage_change_m3": float(
            new_liquid[after_lineage].sum() - old_liquid_current_area[after_lineage].sum()
        ),
        "ice_storage_change_m3": float(
            new_ice[after_lineage].sum() - old_ice_current_area[after_lineage].sum()
        ),
        "persistent_storage_change_m3": physical_delta,
        "connectivity_reassignment_m3": connectivity,
        "external_supply_m3": supply,
        "wet_precip_m3": wet_precip,
        "dry_runoff_m3": dry_runoff,
        "meltwater_m3": meltwater,
        "evaporation_m3": evaporation,
        "liquid_sublimation_m3": liquid_sublimation,
        "ice_sublimation_m3": ice_sublimation,
        "internal_transfer_net_m3": transfer_net,
        "transfer_edges": transfer_edges,
        "residual_m3": residual,
        "normalized_residual": abs(residual) / max(supply, 1.0),
        "attribution_fraction": 1.0 - abs(component_delta - (expected_physical + connectivity)) / max(abs(component_delta), 1.0),
        "persistent_node_count": int(len(persistent_ids)),
        "connectivity_added_nodes": int(len(new_only_ids)),
        "connectivity_removed_nodes": int(len(old_only_ids)),
        "tracked_lake_ids_before": sorted(old_lake_ids),
        "tracked_lake_ids_after": sorted(new_lake_ids),
        "trace_totals": {
            "external_supply_m3": float(trace.external_supply_m3.sum()),
            "wet_precip_m3": float(trace.wet_precip_m3.sum()),
            "dry_runoff_m3": float(trace.dry_runoff_m3.sum()),
            "meltwater_m3": float(trace.meltwater_m3.sum()),
            "evaporation_m3": float(trace.evaporation_m3.sum()),
            "liquid_sublimation_m3": float(trace.liquid_sublimation_m3.sum()),
            "ice_sublimation_m3": float(ice_trace["sublimation_by_node_m3"].sum()),
        },
        "trace_vs_whole_world": {
            "supply_error_m3": float(trace.external_supply_m3.sum() - after.lake_water_ledger.external_supply_m3),
            "wet_precip_error_m3": float(trace.wet_precip_m3.sum() - after.lake_water_ledger.wet_precip_m3),
            "dry_runoff_error_m3": float(trace.dry_runoff_m3.sum() - after.lake_water_ledger.dry_runoff_m3),
            "meltwater_error_m3": float(trace.meltwater_m3.sum() - after.lake_water_ledger.meltwater_m3),
            "evaporation_error_m3": float(
                trace.evaporation_m3.sum() + trace.liquid_sublimation_m3.sum()
                - after.lake_water_ledger.evaporation_m3
            ),
            "ice_sublimation_error_m3": float(
                ice_trace["sublimation_by_node_m3"].sum() - after.lake_water_ledger.sublimation_m3
            ),
            "traced_overflow_error_m3": float(traced_overflow - after.lake_water_ledger.overflow_out_m3),
            "global_trace_residual_m3": global_trace_residual,
            "top_basin_node_residuals": top_local_residuals,
        },
    }


args.out.parent.mkdir(parents=True, exist_ok=True)
result = {"header": header, "anchor_selection": None, "records": [snapshot(0)]}
steps = int(round(args.myr * 1e6 / args.step_years))
for step in range(1, steps + 1):
    before = world.hydrology_cache
    before_storage = storage_snapshot(before)
    old_lineage_ids = None if lineage_ids is None else set(lineage_ids)
    old_lineage_lake_ids = None if lineage_lake_ids is None else set(lineage_lake_ids)
    old_lineage_wet_ids = None if lineage_wet_ids is None else set(lineage_wet_ids)
    old_lineage_component_ids = None if lineage_component_ids is None else set(lineage_component_ids)
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    elapsed = time.perf_counter() - started
    fields = world.hydrology_cache
    bodies, _, _ = components(fields)
    if step == args.anchor_step:
        if static_anchor_points is not None:
            anchor_points = static_anchor_points
            selected_ids, distances = component_ids_for_points(fields, bodies, anchor_points)
            if not selected_ids:
                raise RuntimeError("registered anchors do not overlap a wet component at anchor step")
            anchor_selection = {
                "step": step,
                "components": [
                    {key: value for key, value in bodies[idx].items() if key != "nodes"}
                    for idx in selected_ids
                ],
                "anchor_count": int(len(anchor_points)),
                "nearest_distance_km_max": float(np.max(distances) * 6371.0),
                "anchor_points_xyz": anchor_points.tolist(),
            }
        else:
            if not bodies:
                raise RuntimeError("no wet component exists at the requested anchor step")
            selected = min(
                bodies,
                key=lambda body: abs(body["volume_km3"] - args.reference_volume_km3),
            )
            nodes = selected["nodes"]
            node_points = np.asarray(fields.points)[nodes]
            center = node_points.mean(axis=0)
            first = int(np.argmax(np.linalg.norm(node_points - center, axis=1)))
            chosen_local = [first]
            min_distance = np.linalg.norm(node_points - node_points[first], axis=1)
            while len(chosen_local) < min(args.anchor_count, len(nodes)):
                next_idx = int(np.argmax(min_distance))
                chosen_local.append(next_idx)
                next_distance = np.linalg.norm(node_points - node_points[next_idx], axis=1)
                min_distance = np.minimum(min_distance, next_distance)
            chosen = nodes[np.asarray(chosen_local, dtype=int)]
            anchor_points = np.asarray(fields.points)[chosen].copy()
            anchor_selection = {
                "step": step,
                "volume_km3": selected["volume_km3"],
                "area_km2": selected["area_km2"],
                "anchor_node_count": int(len(chosen)),
                "anchor_points_xyz": anchor_points.tolist(),
            }
            selected_ids = [idx for idx, body in enumerate(bodies) if body is selected]
        result["anchor_selection"] = anchor_selection
        component_node_ids = {
            node_id
            for idx in selected_ids
            for node_id in id_tuples(fields.lake_trace_node_ids[bodies[idx]["nodes"]])
        }
        lineage_component_ids = set(component_node_ids)
        selected_basins = [
            item for item in active_basin_members(fields)
            if component_node_ids.intersection(item[2])
        ]
        if not selected_basins:
            raise RuntimeError("anchor component does not overlap an active lake basin")
        result["anchor_selection"]["basins"] = [
            {"lake_id": lake_id, "member_node_count": len(members), "wet_member_node_count": len(wet_ids)}
            for lake_id, members, wet_ids in selected_basins
        ]
        lineage_ids = set().union(*(members for _, members, _ in selected_basins))
        lineage_lake_ids = {lake_id for lake_id, _, _ in selected_basins}
        lineage_wet_ids = set().union(*(wet_ids for _, _, wet_ids in selected_basins))
        result["lineage_budgets"] = []
    elif old_lineage_ids is not None:
        active_basins = active_basin_members(fields)
        selected_basins = [
            item for item in active_basins if old_lineage_wet_ids.intersection(item[2])
        ]
        if anchor_points is not None:
            anchor_ids, _ = component_ids_for_points(fields, bodies, anchor_points)
            anchor_node_ids = {
                node_id for idx in anchor_ids
                for node_id in id_tuples(fields.lake_trace_node_ids[bodies[idx]["nodes"]])
            }
            lineage_component_ids = set(anchor_node_ids)
            selected_by_id = {lake_id: (members, wet_ids) for lake_id, members, wet_ids in selected_basins}
            for lake_id, members, wet_ids in active_basins:
                if anchor_node_ids.intersection(wet_ids):
                    selected_by_id[lake_id] = (members, wet_ids)
            selected_basins = [(lake_id, members, wet_ids) for lake_id, (members, wet_ids) in selected_by_id.items()]
        lineage_ids = set().union(*(members for _, members, _ in selected_basins)) if selected_basins else set()
        lineage_lake_ids = {lake_id for lake_id, _, _ in selected_basins}
        lineage_wet_ids = set().union(*(wet_ids for _, _, wet_ids in selected_basins)) if selected_basins else set()
        event_budget = lineage_budget(
            before_storage, fields, old_lineage_ids, lineage_ids,
            old_lineage_component_ids, lineage_component_ids,
            old_lineage_lake_ids, lineage_lake_ids,
        )
        result.setdefault("lineage_budgets", []).append({"step": step, **event_budget})
        print(json.dumps({"step": step, "lineage_budget": event_budget}), flush=True)
    record = snapshot(step, elapsed)
    records.append(record)
    result["records"] = records
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(record), flush=True)
print(f"wrote {args.out}")
