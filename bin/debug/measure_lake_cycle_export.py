#!/usr/bin/env python3
"""Measure the lake solve's flush passes and where its final export goes (#696).

Replay a save unperturbed and, each step, record the solve's passes against its limit, internal
transfer, the most repeated body pairs, final export, ocean overflow, lake volume, silt, the
wet component at a fixed anchor, and step time.

Final export is routed again on its own through the step's flow targets and transit losses, to
split it into ocean, in-transit loss, and settling on a cycle or a dead end, on wet or dry
ground. Each exporting pair's parent is classified too: what its outlet would be if the export
used that instead of the leaf sink. Changes no model term, so `--backend` compares two trees.

Example from the repository root:

    backend/.venv/bin/python bin/debug/measure_lake_cycle_export.py <checkpoint.mbworld> \\
      --anchors analysis/issue696/anchor-window-a.json --myr 0.5 --out <scratch>/branch-a.json
"""

from __future__ import annotations

import json
import resource
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from replay_lib import load_replay, parse_replay_args

TOP_PAIRS = 6


def add_arguments(parser) -> None:
    parser.add_argument("--anchors", type=Path, required=True, help="anchor JSON with one XYZ unit point")


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=696)
anchor = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)[0]

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

pairs: dict[tuple[int, int], dict] = {}
exports: list[dict] = []
counts = {"active_bodies": 0, "solves": 0}
overflow_route: dict = {}

trace_transfer = lakes._trace_transfer
trace_deferred_export = lakes._trace_deferred_export
active_bodies = lakes._active_bodies
solve_lakes_prognostic = lakes.solve_lakes_prognostic
route_downstream = hydrology.route_downstream


def observing_transfer(accounting, source, receiver, target, volume, *, internal, depth):
    if internal and receiver is not None:
        row = pairs.setdefault((id(source), id(receiver)), {"sends": 0, "total_m3": 0.0})
        row["sends"] += 1
        row["total_m3"] += float(volume)
    return trace_transfer(accounting, source, receiver, target, volume, internal=internal, depth=depth)


def observing_export(accounting, leaf, sink, volume, depth):
    exports.append({"leaf": leaf, "sink": int(sink), "volume_m3": float(volume)})
    return trace_deferred_export(accounting, leaf, sink, volume, depth)


def observing_active_bodies(*body_args, **body_kwargs):
    counts["active_bodies"] += 1
    return active_bodies(*body_args, **body_kwargs)


def observing_solve(*solve_args, **solve_kwargs):
    # A later solve in the same step replaces the earlier one's result, so measure only the last.
    pairs.clear()
    exports.clear()
    counts["active_bodies"] = 0
    counts["solves"] += 1
    return solve_lakes_prognostic(*solve_args, **solve_kwargs)


def observing_route(elevation, is_ocean, flow_target, source_amount, *route_args, **route_kwargs):
    # The overflow route is the step's last one that carries only a transit loss.
    if "loss_fraction" in route_kwargs and "erosion_only_source_amount" not in route_kwargs:
        overflow_route.update(
            elevation=elevation, is_ocean=is_ocean, flow_target=flow_target,
            loss_fraction=route_kwargs["loss_fraction"],
        )
    return route_downstream(elevation, is_ocean, flow_target, source_amount, *route_args, **route_kwargs)


lakes._trace_transfer = observing_transfer
lakes._trace_deferred_export = observing_export
lakes._active_bodies = observing_active_bodies
lakes.solve_lakes_prognostic = observing_solve
hydrology.route_downstream = observing_route


def route_fate(source_m3: np.ndarray, lake_depth: np.ndarray) -> dict:
    """Where `source_m3` ends up if routed alone through this step's overflow route."""
    is_ocean = overflow_route["is_ocean"]
    flow_target = overflow_route["flow_target"]
    _, settled = route_downstream(
        overflow_route["elevation"], is_ocean, flow_target, source_m3,
        loss_fraction=overflow_route["loss_fraction"],
    )
    settled = settled.copy()
    settled[is_ocean] += source_m3[is_ocean]
    land = ~is_ocean
    wet = land & (lake_depth > 0.0)
    on_cycle = land & (flow_target >= 0)
    return {
        "source_m3": float(source_m3.sum()),
        "ocean_m3": float(settled[is_ocean].sum()),
        "transit_loss_m3": float(source_m3.sum() - settled.sum()),
        "settled_wet_m3": float(settled[wet].sum()),
        "settled_dry_m3": float(settled[land & ~wet].sum()),
        "settled_on_flow_cycle_m3": float(settled[on_cycle].sum()),
        "settled_at_dead_end_m3": float(settled[land & ~on_cycle].sum()),
    }


def parent_outlets(fields) -> dict:
    """Classify each exporting body's parent outlet, weighted by exported volume."""
    forest = fields.lake_forest or []
    parent_of = {id(child): lake for lake in lakes.iter_all_lakes(forest) for child in lake.children}
    bodies = lakes.governing_bodies(forest, fields.elevation, fields.lake_depth, fields.hypsometry)
    body_of_leaf = {}
    for body in bodies:
        for leaf in lakes.iter_all_lakes([body]):
            body_of_leaf[id(leaf)] = body
    is_ocean = overflow_route["is_ocean"]
    leaf_nodes = np.zeros(len(is_ocean), dtype=bool)
    for lake in lakes.iter_all_lakes(forest):
        if not lake.children:
            leaf_nodes[lake.members] = True
    classes: dict[str, float] = defaultdict(float)
    source = np.zeros(len(is_ocean))
    for item in exports:
        body = body_of_leaf.get(id(item["leaf"]))
        parent = parent_of.get(id(body)) if body is not None else None
        target = -1 if parent is None else int(parent.outlet_target_idx)
        if parent is None:
            kind = "no_parent"
        elif target < 0:
            kind = "endorheic_parent"
        elif is_ocean[target]:
            kind = "ocean"
        elif leaf_nodes[target]:
            kind = "another_basin"
        else:
            kind = "river"
        classes[kind] += item["volume_m3"]
        if target >= 0:
            source[target] += item["volume_m3"]
    return {"by_class_m3": dict(classes), "routed_from_parent_outlet": route_fate(source, np.asarray(fields.lake_depth))}


def anchor_component(fields) -> dict | None:
    depth = np.asarray(fields.lake_depth)
    area = np.asarray(fields.hypsometry.area_m2)
    best, best_distance = None, np.inf
    for nodes in hydrology.lake_components(depth > 0.0, fields.neighbor_idx):
        distance = float(np.min(np.linalg.norm(fields.points[nodes] - anchor, axis=1)))
        if distance < best_distance:
            best, best_distance = nodes, distance
    if best is None:
        return None
    return {"nodes": len(best), "volume_m3": float(np.dot(depth[best], area[best]))}


def snapshot(step: int, seconds: float) -> dict:
    fields = world.hydrology_cache
    ledger = fields.lake_water_ledger
    area = np.asarray(fields.hypsometry.area_m2)
    depth = np.asarray(fields.lake_depth)
    leaves = sum(1 for lake in lakes.iter_all_lakes(fields.lake_forest or []) if not lake.children)
    # `_active_bodies` runs once before and once after each pass.
    passes = counts["active_bodies"] // 2
    export_by_node = np.zeros(len(depth))
    for item in exports:
        export_by_node[item["sink"]] += item["volume_m3"]
    ranked = sorted(pairs.values(), key=lambda row: -row["total_m3"])
    silt_m3 = float(np.dot(fields.silt_deposited, area))
    return {
        "step": step,
        "age_myr": float(world.elapsed_years) / 1e6,
        "seconds": seconds,
        "solves": counts["solves"],
        "leaves": leaves,
        "passes": passes,
        "pass_limit": 1 + leaves,
        "hit_limit": passes >= 1 + leaves,
        "external_supply_m3": float(ledger.external_supply_m3),
        "internal_transfer_m3": float(ledger.internal_transfer_m3),
        "overflow_out_m3": float(ledger.overflow_out_m3),
        "overflow_to_ocean_m3": float(ledger.overflow_to_ocean_m3),
        "residual_m3": float(ledger.residual_m3),
        "lake_liquid_volume_m3": float(np.dot(depth, area)),
        "silt_m3": silt_m3,
        "silt_mean_mm_yr": silt_m3 / float(area[depth > 0.0].sum()) / args.step_years * 1e3,
        "max_pair_sends": max((row["sends"] for row in pairs.values()), default=0),
        "top_pairs": [{"sends": row["sends"], "total_m3": row["total_m3"]} for row in ranked[:TOP_PAIRS]],
        "final_export_m3": float(export_by_node.sum()),
        "final_export_sinks": len(np.flatnonzero(export_by_node)),
        "final_export_fate": route_fate(export_by_node, depth),
        "parent_outlet": parent_outlets(fields),
        "anchor_component": anchor_component(fields),
    }


world.ice_age_period_years = 0.0
args.out.parent.mkdir(parents=True, exist_ok=True)
result = {**header, "anchor": anchor.tolist(), "records": []}
for step in range(1, int(round(args.myr * 1e6 / args.step_years)) + 1):
    counts["solves"] = 0
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    seconds = time.perf_counter() - started
    result["records"].append(snapshot(step, seconds))
    result["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(f"step {step}: {seconds:.1f}s", flush=True)
print(f"wrote {args.out}")
