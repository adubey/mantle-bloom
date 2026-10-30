#!/usr/bin/env python3
"""Budget every surface pathway at the spill connectors around a lake (#691).

Replay a checkpoint unperturbed and, each step, find the wet component nearest a fixed anchor
and every hierarchy lake containing it. For each of those lakes' rim edges, record the bed and
passage state plus water-derived lake silt, routed sediment and the step's net geomorphic
change, and the terrain change outside erosion measured at the next step. Also record the
component's lake-wide budget and its largest inlets. The #680 preflight read `silt_deposited`
only, so it could not separate that pathway from sediment delivery or from the sill's other
elevation changes. Each record also gives the largest supply any final body's silt implies this
step, for the #693 bound against world external supply, and on code with the #693 count, the
largest supply before that bound.

A rim is a model graph edge between two cells, not a surveyed geomorphic sill. This script
changes no model term and does not test a mechanism.

Example from the repository root:

    backend/.venv/bin/python bin/debug/measure_lake_connector_budget.py \\
      /tmp/seed910211954-221.6myr.mbworld \\
      --anchors analysis/issue691/anchor-window-a.json \\
      --myr 0.4 --step-years 100000 --out /tmp/issue691-window-a.json
"""

from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args

EARTH_RADIUS_KM = 6371.0
MAX_CONTAINING_LAKES = 8
TOP_INLETS = 5


def add_arguments(parser) -> None:
    parser.add_argument("--anchors", type=Path, required=True, help="anchor JSON with one XYZ unit point")


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=691)

from app import plates as plates_mod  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

anchor = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)[0]
step_transfers: list[tuple[int, float]] = []
liquid_depth_after_solve: list[np.ndarray] = []
max_body_implied_supply: list[float] = []
max_body_unbounded_supply: list[float | None] = [None]
trace_transfer = lakes._trace_transfer
solve_lakes_prognostic = lakes.solve_lakes_prognostic


def counting_trace_transfer(accounting, source, receiver, target, volume, *, internal, depth):
    """Observe each internal transfer, then do what the solver would have done."""
    if internal:
        step_transfers.append((id(source), float(volume)))
    return trace_transfer(accounting, source, receiver, target, volume, internal=internal, depth=depth)


def observing_solve(*solve_args, **solve_kwargs):
    """Keep the solver's own depth, before `advance_lake_ice` rewrites it in place."""
    result = solve_lakes_prognostic(*solve_args, **solve_kwargs)
    liquid_depth_after_solve[:] = [result.lake_depth.copy()]
    forest, elevation, hypsometry = solve_args[0], solve_args[1], solve_args[5]
    bodies = lakes.governing_bodies(forest, elevation, result.lake_depth, hypsometry)
    silt_m3 = [
        float(np.dot(result.silt_deposited[body.members], hypsometry.area_m2[body.members]))
        for body in bodies
    ]
    max_body_implied_supply[:] = [max(silt_m3, default=0.0) / lakes.SILT_ACCUMULATION_COEFFICIENT]
    return result


def observing_body_supply(*supply_args):
    """Note the largest body supply before the world bound, then return it unchanged."""
    supply = final_body_silt_supply(*supply_args)
    max_body_unbounded_supply[:] = [max(supply, default=0.0)]
    return supply


lakes._trace_transfer = counting_trace_transfer
lakes.solve_lakes_prognostic = observing_solve
final_body_silt_supply = getattr(lakes, "_final_body_silt_supply", None)
if final_body_silt_supply is not None:
    lakes._final_body_silt_supply = observing_body_supply
world.ice_age_period_years = 0.0
step_years = float(args.step_years)


def surface_km(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.arccos(np.clip(np.dot(a, b), -1.0, 1.0)) * EARTH_RADIUS_KM)


def nearest_component(fields) -> np.ndarray | None:
    wet = np.asarray(fields.lake_depth) > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    best, best_distance = None, np.inf
    for nodes in hydrology.lake_components(wet, fields.neighbor_idx):
        nodes = np.asarray(nodes, dtype=np.int64)
        distance = float(np.min(np.linalg.norm(fields.points[nodes] - anchor, axis=1)))
        if distance < best_distance:
            best, best_distance = nodes, distance
    return best


def node_state(fields, erosion, notch, channel, node: int) -> dict:
    silt = float(fields.silt_deposited[node])
    routed = float(erosion.sediment_deposited[node])
    net = float(erosion.net_elevation_change_m[node])
    return {
        "node": node,
        "xyz": [float(v) for v in fields.points[node]],
        "bed_m": float(fields.elevation[node]),
        "lake_depth_m": float(fields.lake_depth[node]),
        "wet_fraction_positive": bool(fields.lake_depth[node] > 0.0),
        "is_ocean": bool(fields.is_ocean[node]),
        "silt_m": silt,
        # Excludes lake silt (`apply_erosion` books that column apart from total_deposited).
        "sediment_deposited_m": routed,
        "net_geomorphic_m": net,
        # Airy subsidence of the deposit plus any removal: only about 1/6-1/4 of a deposit
        # survives as surface rise (`apply_erosion`), so this is not an erosion measure alone.
        "deposit_minus_net_m": silt + routed - net,
        "breach_notch_m": float(notch[node]),
        "channel_depth_m": float(channel[node]),
    }


def lake_record(fields, erosion, notch, channel, lake, surface: float) -> dict:
    rim, target = int(lake.outlet_node_idx), int(lake.outlet_target_idx)
    spill = None if lake.max_depth is None else float(lake.max_depth)
    record = {
        "lake_id": int(lake.lake_id),
        "member_count": len(lake.members),
        "is_spilling": bool(lake.is_spilling),
        "spill_level_m": spill,
        "clearance_m": None if spill is None else spill - surface,
    }
    if rim >= 0 and target >= 0:
        record["rim"] = node_state(fields, erosion, notch, channel, rim)
        record["target"] = node_state(fields, erosion, notch, channel, target)
        record["edge_km"] = surface_km(fields.points[rim], fields.points[target])
    return record


def anchor_body_chain(fields, erosion, notch, channel, component: np.ndarray) -> dict:
    """The hierarchy around the wet cell nearest the anchor, judged by that cell's own water
    level. A wet component can join bodies standing at different levels, so a whole-component
    surface can read kilometres from the anchor's. Keeps the highest saddle already drowned
    and the next ones above the level; those are the connectors this water could be using."""
    depth = np.asarray(fields.lake_depth)
    node = int(component[np.argmin(np.linalg.norm(fields.points[component] - anchor, axis=1))])
    level = float(fields.hypsometry.level_from_depth(depth[[node]], np.array([node]))[0])
    chain = sorted(
        (lake for lake in lakes.iter_all_lakes(fields.lake_forest or []) if node in lake.members),
        key=lambda lake: len(lake.members),
    )
    resolved = [lake for lake in chain if lake.max_depth is not None]
    below = [lake for lake in resolved if lake.max_depth <= level]
    above = [lake for lake in resolved if lake.max_depth > level]
    selected = below[-2:] + above[:3]
    return {
        "node": node,
        "level_m": level,
        "chain_length": len(chain),
        "lakes": [lake_record(fields, erosion, notch, channel, lake, level) for lake in selected],
    }


def component_levels(fields, component: np.ndarray) -> list[dict]:
    """Water levels in the component, rounded to 1 m, largest volume first."""
    depth = np.asarray(fields.lake_depth)[component]
    area = np.asarray(fields.hypsometry.area_m2)[component]
    levels = np.round(fields.hypsometry.level_from_depth(depth, component))
    rows = [
        {
            "level_m": float(value),
            "nodes": int(np.sum(levels == value)),
            "volume_km3": float(np.dot(depth[levels == value], area[levels == value]) / 1e9),
        }
        for value in np.unique(levels)
    ]
    return sorted(rows, key=lambda row: -row["volume_km3"])[:6]


def transfer_summary() -> dict:
    """Whether internal volume is many one-way hops or the same bodies sending repeatedly."""
    volumes = np.array([volume for _, volume in step_transfers]) if step_transfers else np.zeros(0)
    sends: dict[int, int] = {}
    for source, _ in step_transfers:
        sends[source] = sends.get(source, 0) + 1
    counts = np.array(list(sends.values())) if sends else np.zeros(0)
    return {
        "count": len(volumes),
        "total_m3": float(volumes.sum()),
        "largest_m3": float(volumes.max()) if len(volumes) else 0.0,
        "distinct_sources": len(sends),
        "max_sends_by_one_source": int(counts.max()) if len(counts) else 0,
    }


def zero_depth_silt(fields, area: np.ndarray) -> dict:
    """Silt on cells reading no liquid at step end, split by what the solver itself reported."""
    silted = np.asarray(fields.silt_deposited) > 0.0
    solved = liquid_depth_after_solve[0]
    final = np.asarray(fields.lake_depth)
    zero_final = silted & (final <= 0.0)
    zero_solved = silted & (solved <= 0.0)
    return {
        "nodes_final_zero": int(zero_final.sum()),
        "nodes_solver_zero": int(zero_solved.sum()),
        "silt_m3_final_zero": float(np.dot(fields.silt_deposited[zero_final], area[zero_final])),
        "silt_m3_solver_zero": float(np.dot(fields.silt_deposited[zero_solved], area[zero_solved])),
        "solver_depth_median_m": float(np.median(solved[zero_final])) if zero_final.any() else None,
        "lake_ice_median_m": float(np.median(fields.lake_ice_m[zero_final])) if zero_final.any() else None,
    }


def inlets(fields, erosion, component: np.ndarray, rim_points: list[np.ndarray]) -> list[dict]:
    """The dry cells that drain into the component, largest first, with their distance to the
    nearest rim and the routed sediment on the cell."""
    in_component = np.zeros(len(fields.lake_depth), dtype=bool)
    in_component[component] = True
    target = np.asarray(fields.flow_target)
    feeds = (target >= 0) & ~in_component & ~np.asarray(fields.is_ocean)
    feeds[feeds] = in_component[target[feeds]]
    discharge = np.asarray(fields.discharge_m3_s)
    nodes = np.nonzero(feeds)[0]
    nodes = nodes[np.argsort(-discharge[nodes])][:TOP_INLETS]
    return [
        {
            "node": int(node),
            "discharge_m3_s": float(discharge[node]),
            "sediment_deposited_m": float(erosion.sediment_deposited[node]),
            "distance_to_nearest_rim_km": (
                min(surface_km(fields.points[node], rim) for rim in rim_points) if rim_points else None
            ),
        }
        for node in nodes
    ]


def snapshot(step: int) -> dict:
    fields, erosion = world.hydrology_cache, world.erosion_cache
    area = np.asarray(fields.hypsometry.area_m2)
    record: dict = {"step": step, "age_years": float(world.elapsed_years)}
    component = nearest_component(fields)
    if component is None or erosion is None:
        return record
    notch = plates_mod.collect_all_breach_notch_depth(fields.plates_in_order)
    channel = plates_mod.collect_all_channel_depth(fields.plates_in_order)
    depth = np.asarray(fields.lake_depth)
    surface = float(np.max(fields.hypsometry.level_from_depth(depth[component], component)))
    containing = sorted(
        (lake for lake in lakes.iter_all_lakes(fields.lake_forest or []) if np.isin(component, lake.members).all()),
        key=lambda lake: len(lake.members),
    )[:MAX_CONTAINING_LAKES]
    records = [lake_record(fields, erosion, notch, channel, lake, surface) for lake in containing]
    rim_points = [np.asarray(item["rim"]["xyz"]) for item in records if "rim" in item]
    component_area = float(area[component].sum())
    silt_m3 = float(np.dot(fields.silt_deposited[component], area[component]))
    routed_m3 = float(np.dot(erosion.sediment_deposited[component], area[component]))
    ledger = fields.lake_water_ledger
    record.update({
        "anchor_distance_km": float(np.min(np.linalg.norm(fields.points[component] - anchor, axis=1)) * EARTH_RADIUS_KM),
        "component": {
            "nodes": len(component),
            "area_km2": component_area / 1e6,
            "volume_km3": float(np.dot(depth[component], area[component]) / 1e9),
            "surface_m": surface,
            "silt_m3": silt_m3,
            "sediment_deposited_m3": routed_m3,
            "silt_mean_mm_yr": silt_m3 / component_area / step_years * 1e3,
            "sediment_mean_mm_yr": routed_m3 / component_area / step_years * 1e3,
            "silt_max_m": float(np.max(fields.silt_deposited[component])),
            "sediment_max_m": float(np.max(erosion.sediment_deposited[component])),
            "bed_min_m": float(np.min(fields.elevation[component])),
            "bed_max_m": float(np.max(fields.elevation[component])),
        },
        "world_lake_silt_m3": float(np.dot(fields.silt_deposited, area)),
        "median_edge_km": float(np.median(
            np.linalg.norm(fields.points[:, None, :] - fields.points[fields.neighbor_idx], axis=2)
        ) * EARTH_RADIUS_KM),
        "containing_lakes": records,
        "anchor_body": anchor_body_chain(fields, erosion, notch, channel, component),
        "component_levels": component_levels(fields, component),
        "implied_supply_m3": silt_m3 / lakes.SILT_ACCUMULATION_COEFFICIENT,
        "world_implied_supply_m3": float(np.dot(fields.silt_deposited, area)) / lakes.SILT_ACCUMULATION_COEFFICIENT,
        "max_body_implied_over_world_external": (
            max_body_implied_supply[0] / max(float(ledger.external_supply_m3), 1.0)
        ),
        "max_body_unbounded_supply_m3": max_body_unbounded_supply[0],
        "transfers": transfer_summary(),
        "zero_depth_silt": zero_depth_silt(fields, area),
        # Silt this component would hold if its supply were only the external water settling
        # on its own nodes: a through-flow-free floor, not a proposed rule.
        "component_external_supply_m3": float(
            fields.lake_external_supply_m3_per_yr[component].sum() * step_years
        ),
        "ledger_external_supply_m3": float(ledger.external_supply_m3),
        "ledger_internal_transfer_m3": float(ledger.internal_transfer_m3),
        "ledger_overflow_out_m3": float(ledger.overflow_out_m3),
        "world_silt_on_zero_depth_m3": float(
            np.dot(np.where(depth <= 0.0, fields.silt_deposited, 0.0), area)
        ),
        "inlets": inlets(fields, erosion, component, rim_points),
        "normalized_ledger_residual": abs(float(ledger.residual_m3)) / max(float(ledger.external_supply_m3), 1.0),
    })
    return record


def add_terrain_change_outside_erosion(previous: dict, fields) -> None:
    """Fill each rim edge's change from tectonics, isostasy and everything else outside
    `net_elevation_change_m`: the next step's bed at the nearest node, less this step's
    post-erosion bed. Nodes ride their plates and may be resampled, so the match is geographic
    and its distance is reported."""
    tree = cKDTree(fields.points)
    lakes_seen = previous.get("containing_lakes", []) + previous.get("anchor_body", {}).get("lakes", [])
    for lake in lakes_seen:
        for key in ("rim", "target"):
            node = lake.get(key)
            if node is None:
                continue
            distance, nearest = tree.query(node["xyz"], k=1)
            node["next_bed_m"] = float(fields.elevation[nearest])
            node["next_match_km"] = float(distance * EARTH_RADIUS_KM)
            node["outside_erosion_m"] = float(
                fields.elevation[nearest] - (node["bed_m"] + node["net_geomorphic_m"])
            )


args.out.parent.mkdir(parents=True, exist_ok=True)
result = {**header, "anchor": anchor.tolist(), "records": []}
steps = int(round(args.myr * 1e6 / step_years))
for step in range(1, steps + 1):
    step_transfers.clear()
    world_mod.step_world(world, step_years)
    if result["records"]:
        add_terrain_change_outside_erosion(result["records"][-1], world.hydrology_cache)
    result["records"].append(snapshot(step))
    result["peak_rss_gb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: result["records"][-1].get(k) for k in ("step", "age_years", "component")}), flush=True)
print(f"wrote {args.out}")
