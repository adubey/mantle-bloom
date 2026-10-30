#!/usr/bin/env python3
"""Measure lake sublimation, level spread and cold lakes step by step (#694).

Replay a checkpoint with `ice_age_period_years = 0` and record, each step, the lake ledger, the
wet component at a fixed anchor, how much each solved body's level surface was disturbed after
the solve, cold-lake state and the ice-cap checks. The script only observes. The arm under test
is whichever `--backend` it imports: `origin/main`, or a scratch tree carrying
`analysis/issue694/body-sublimation.patch`.

`--sublimation-scale` multiplies `LAKE_ICE_SUBLIMATION_M_PER_YR`, for the zero-rate control and
the rate-sensitivity arm.

Run one at a time, from the repository root:

    backend/.venv/bin/python bin/debug/measure_lake_sublimation_probe.py \\
      /tmp/seed910211954-221.6myr.mbworld --anchors analysis/issue696/anchor-window-a.json \\
      --backend ../mantle-bloom-issue-694-baseline/backend \\
      --myr 0.5 --step-years 100000 --out /tmp/issue694-a-current-100k.json
"""

from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args

EARTH_RADIUS_KM = 6371.0
ANCHOR_TOLERANCE_KM = 20.0  # As bin/debug/trace_lake_anchor.py.
LARGE_BODY_KM2 = 1.0e4  # RFC 0014's Flicker size.
MIXED_LEVEL_M = 1.0
COLD_COVER = 0.5
SHALLOW_WATER_M = 2.0  # RFC 0014 G5's split of bedfast area.


def add_arguments(parser) -> None:
    parser.add_argument("--anchors", type=Path, required=True, help="anchor JSON with one XYZ unit point")
    parser.add_argument("--sublimation-scale", type=float, default=1.0)


args = parse_replay_args(__doc__, add_arguments)
world, header = load_replay(args, issue=694)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

anchor = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)[0]
hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR *= args.sublimation_scale
world.ice_age_period_years = 0.0
captured: dict = {}
solve_lakes_prognostic = lakes.solve_lakes_prognostic
advance_lake_ice = hydrology.advance_lake_ice


def observing_solve(*solve_args, **solve_kwargs):
    """Keep the solver's bodies and depth, before `advance_lake_ice` replaces the depth in place."""
    result = solve_lakes_prognostic(*solve_args, **solve_kwargs)
    forest, elevation, hypsometry = solve_args[0], solve_args[1], solve_args[5]
    captured.update(
        solved_depth=result.lake_depth.copy(),
        bodies=[body.members for body in lakes.governing_bodies(forest, elevation, result.lake_depth, hypsometry)],
        solver_loss_m3=float(result.ledger.evaporation_m3),
        solver_sublimation_m3=float(result.ledger.sublimation_m3),
    )
    return result


def observing_advance_lake_ice(*ice_args, **ice_kwargs):
    captured["ice_cover"] = np.array(ice_args[5])
    return advance_lake_ice(*ice_args, **ice_kwargs)


lakes.solve_lakes_prognostic = observing_solve
hydrology.advance_lake_ice = observing_advance_lake_ice


def anchor_state(fields, wet: np.ndarray, area: np.ndarray) -> dict:
    """The fixed anchor's wet/dry class (#680) and the wet component nearest it (#691, #693),
    so the event `ΔV` reads the same way as the earlier G8 reports."""
    nearest = int(np.argmin(np.linalg.norm(fields.points - anchor, axis=1)))
    distance_km = float(np.linalg.norm(fields.points[nearest] - anchor) * EARTH_RADIUS_KM)
    state = {
        "nearest_cell_km": distance_km,
        "anchor_wet": bool(wet[nearest]) and distance_km <= ANCHOR_TOLERANCE_KM,
    }
    components = [np.asarray(nodes, dtype=np.int64) for nodes in hydrology.lake_components(wet, fields.neighbor_idx)]
    if not components:
        return state
    distances = [float(np.min(np.linalg.norm(fields.points[nodes] - anchor, axis=1))) for nodes in components]
    nodes = components[int(np.argmin(distances))]
    depth, ice = np.asarray(fields.lake_depth)[nodes], np.asarray(fields.lake_ice_m)[nodes]
    levels = fields.hypsometry.level_from_depth(depth, nodes)
    state.update({
        "component_distance_km": min(distances) * EARTH_RADIUS_KM,
        "component_nodes": len(nodes),
        "component_area_km2": float(area[nodes].sum()) / 1e6,
        "component_volume_km3": float(np.dot(depth, area[nodes])) / 1e9,
        "component_ice_km3": float(np.dot(ice, area[nodes])) / 1e9,
        "component_level_spread_m": float(levels.max() - levels.min()),
    })
    return state


def body_disturbance(fields, area: np.ndarray) -> dict:
    """How far the step's end state stands from the solver's level surfaces (#691). Each body
    leaves the solve at one level, so any spread is what happened to it afterwards."""
    solved, final = captured["solved_depth"], np.asarray(fields.lake_depth)
    spreads, large_areas, emptied_large = [], [], []
    for members in captured["bodies"]:
        wet_members = members[solved[members] > 0.0]
        solved_area = float(area[wet_members].sum()) / 1e6
        if solved_area < LARGE_BODY_KM2:
            continue
        kept = wet_members[final[wet_members] > 0.0]
        levels = fields.hypsometry.level_from_depth(final[kept], kept) if len(kept) else np.zeros(1)
        spreads.append(float(levels.max() - levels.min()))
        large_areas.append(solved_area)
        emptied_large.append(float(area[wet_members[final[wet_members] <= 0.0]].sum()) / 1e6)
    spreads_array, areas_array = np.asarray(spreads), np.asarray(large_areas)
    emptied = (solved > 0.0) & (final <= 0.0)
    return {
        "large_bodies": len(spreads),
        "large_bodies_mixed": int(np.sum(spreads_array > MIXED_LEVEL_M)),
        "large_body_spread_median_m": float(np.median(spreads_array)) if spreads else 0.0,
        "large_body_spread_max_m": float(spreads_array.max()) if spreads else 0.0,
        "large_body_emptied_area_pct": 100.0 * sum(emptied_large) / max(float(areas_array.sum()), 1.0),
        "emptied_cells": int(emptied.sum()),
        "emptied_solved_km3": float(np.dot(solved[emptied], area[emptied])) / 1e9,
    }


def cold_state(fields, wet: np.ndarray, area: np.ndarray, previous_cold: dict | None) -> tuple[dict, dict]:
    """Cold water at step end. Nodes ride plates, so last step's cold wet cells are matched to
    this step's geographically."""
    cover = captured["ice_cover"]
    final = np.asarray(fields.lake_depth)
    cold_wet = wet & (cover >= COLD_COVER)
    record = {
        "cold_wet_area_km2": float(area[cold_wet].sum()) / 1e6,
        "cold_liquid_km3": float(np.dot(final[cold_wet], area[cold_wet])) / 1e9,
        "cold_solved_liquid_km3": float(np.dot(captured["solved_depth"][cover >= COLD_COVER], area[cover >= COLD_COVER])) / 1e9,
        "cold_components": len(hydrology.lake_components(cold_wet, fields.neighbor_idx)),
    }
    if previous_cold is not None and len(previous_cold["points"]):
        _, matched = cKDTree(fields.points).query(previous_cold["points"], k=1)
        record["cold_persistence_pct"] = 100.0 * float(
            np.dot(previous_cold["area"], wet[matched]) / max(float(previous_cold["area"].sum()), 1.0)
        )
    return record, {"points": np.array(fields.points[cold_wet]), "area": area[cold_wet].copy()}


def ice_checks(fields, area: np.ndarray) -> dict:
    """G5's checks. Bedfast means no liquid is left under the ice, split at 2 m of water."""
    liquid = np.clip(np.asarray(fields.lake_depth), 0.0, None)
    ice = np.asarray(fields.lake_ice_m)
    physical = ice / hydrology.LAKE_ICE_DENSITY_RATIO
    water = liquid + ice
    holding = water > 0.0
    bedfast = holding & (liquid <= 0.0)
    shallow = water <= SHALLOW_WATER_M
    return {
        "ice_km3": float(np.dot(ice, area)) / 1e9,
        "ice_max_thickness_m": float(physical.max()) if len(physical) else 0.0,
        "ice_bound_violations": int(np.sum((ice < 0.0) | (physical > hydrology.LAKE_ICE_MAX_THICKNESS_M + 1e-9))),
        "bedfast_deep_pct_of_deep": 100.0 * float(area[bedfast & ~shallow].sum()) / max(float(area[holding & ~shallow].sum()), 1.0),
        "bedfast_shallow_pct_of_water": 100.0 * float(area[bedfast & shallow].sum()) / max(float(area[holding].sum()), 1.0),
    }


def snapshot(step: int, seconds: float, previous_cold: dict | None) -> tuple[dict, dict]:
    fields = world.hydrology_cache
    area = np.asarray(fields.hypsometry.area_m2)
    depth = np.asarray(fields.lake_depth)
    wet = depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    ledger = fields.lake_water_ledger
    supply = max(float(ledger.external_supply_m3), 1.0)
    cold, cold_cells = cold_state(fields, wet, area, previous_cold)
    record = {
        "step": step,
        "age_years": float(world.elapsed_years),
        "seconds": seconds,
        "external_supply_km3": float(ledger.external_supply_m3) / 1e9,
        "evaporation_km3": float(ledger.evaporation_m3) / 1e9,
        "sublimation_km3": float(ledger.sublimation_m3) / 1e9,
        "sublimation_km3_yr": float(ledger.sublimation_m3) / 1e9 / args.step_years,
        "sublimation_pct_of_supply": 100.0 * float(ledger.sublimation_m3) / supply,
        "continuous_sublimation_km3": float(np.dot(
            hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR * np.clip(captured["ice_cover"], 0.0, 1.0) * args.step_years,
            np.where(captured["solved_depth"] > 0.0, area, 0.0),
        )) / 1e9,
        "freezing_km3": float(ledger.freezing_m3) / 1e9,
        "melting_km3": float(ledger.melting_m3) / 1e9,
        "overflow_out_km3": float(ledger.overflow_out_m3) / 1e9,
        "overflow_to_ocean_km3": float(ledger.overflow_to_ocean_m3) / 1e9,
        "normalized_residual": abs(float(ledger.residual_m3)) / supply,
        # Booked once: in the Body arm, ledger evaporation plus sublimation is the solve's whole loss.
        "solver_loss_km3": captured["solver_loss_m3"] / 1e9,
        "solver_sublimation_km3": captured["solver_sublimation_m3"] / 1e9,
        "liquid_km3": float(np.dot(depth, area)) / 1e9,
        "solved_liquid_km3": float(np.dot(captured["solved_depth"], area)) / 1e9,
        "lake_area_km2": float(area[wet].sum()) / 1e6,
        "anchor": anchor_state(fields, wet, area),
        "bodies": body_disturbance(fields, area),
        "cold": cold,
        "ice": ice_checks(fields, area),
        "peak_rss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
    }
    return record, cold_cells


steps = int(round(args.myr * 1e6 / args.step_years))
header.update(sublimation_scale=args.sublimation_scale, sublimation_m_per_yr=hydrology.LAKE_ICE_SUBLIMATION_M_PER_YR)
output = {**header, "records": []}
args.out.parent.mkdir(parents=True, exist_ok=True)
cold_cells = None
for step in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    seconds = time.perf_counter() - started
    record, cold_cells = snapshot(step, seconds, cold_cells)
    output["records"].append(record)
    args.out.write_text(json.dumps(output, indent=1))
    print(json.dumps({key: record[key] for key in ("step", "age_years", "normalized_residual", "peak_rss_gb")}), flush=True)
print(f"wrote {args.out}")
