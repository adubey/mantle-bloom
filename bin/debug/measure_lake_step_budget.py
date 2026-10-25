#!/usr/bin/env python3
"""Issue #588: compare the lake step budget and trajectories on two backends.

Replays a saved world for the requested horizon, recording per-step liquid lake ledger
terms (when available), lake volume/area, large-lake flicker reversals, breaches, floor-drain
toggles and runtime. On a pre-ledger backend, also estimates uncapped balance surplus and
surplus held at the depth cap by tapping `_water_balance`.

Run from `backend/` so the issue worktree's environment is active:

    .venv/bin/python ../bin/debug/measure_lake_step_budget.py <save> --myr 1 \\
        --step-years 100000 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_lake_step_budget.py <save> --myr 1 \\
        --step-years 100000 --backend <baseline>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import inspect
import json
import resource
import sys
import time

import numpy as np
from replay_lib import load_replay, parse_replay_args

args = parse_replay_args(__doc__)
world, result = load_replay(args, issue=588)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

world.ice_age_period_years = 0.0
visible_depth = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
balance_calls: list[tuple] = []
original_balance = lakes._water_balance
balance_signature = inspect.signature(original_balance)


def tap_balance(*positional, **keyword):
    level = original_balance(*positional, **keyword)
    bound = balance_signature.bind(*positional, **keyword)
    bound.apply_defaults()
    p = bound.arguments
    balance_calls.append((p["lake"], level, p["tier_max_depth"], p["climate"], p["hypsometry"]))
    return level


if hasattr(lakes, "_water_balance"):
    lakes._water_balance = tap_balance


def legacy_surplus(fields):
    """Return pre-ledger unbooked surplus and cap-held surplus in m³ for this step."""
    if not balance_calls:
        return 0.0, 0.0
    unbooked = capped = 0.0
    for lake, level, tier_depth, climate, hyps in balance_calls:
        members = lake.members
        if lake.children and level < lake.min_depth:
            continue
        if not np.allclose(fields.lake_depth[members], hyps.mean_depth(level, members), rtol=1e-9, atol=1e-9):
            continue
        sinks = lakes.leaf_sinks(lake)
        if fields.should_spill is not None and (lake.drains_at_floor or (
            lake.max_depth is not None and level >= lake.max_depth - 1e-6
        )) and fields.should_spill[sinks].any():
            continue
        supplied = float(np.maximum(climate.inflow_m3_per_yr[members], 0.0).sum()) * args.step_years
        evaporated = lakes.lake_evaporation_m3_per_yr(level, members, hyps, climate) * args.step_years
        surplus = max(supplied - evaporated, 0.0)
        rim_depth = None if lake.max_depth is None else lake.max_depth - lake.floor_elevation
        at_cap = (rim_depth is None or tier_depth < rim_depth) and level >= lake.floor_elevation + tier_depth - 1e-6
        if at_cap:
            capped += surplus
        else:
            unbooked += surplus
    return unbooked, capped


def components(fields):
    wet = fields.lake_depth > visible_depth
    members = hydrology.lake_components(wet, fields.neighbor_idx)
    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area = hypsometry.area_m2
    records = []
    for nodes in members:
        records.append({
            "nodes": nodes,
            "area_km2": float(area[nodes].sum() / 1e6),
            "volume_m3": float(np.dot(fields.lake_depth[nodes], area[nodes])),
        })
    return records


def match_bodies(previous, current):
    candidates = []
    for old_i, old in enumerate(previous):
        for new_i, new in enumerate(current):
            overlap = np.intersect1d(old["nodes"], new["nodes"], assume_unique=True).size
            if overlap:
                candidates.append((overlap, old_i, new_i))
    candidates.sort(reverse=True)
    old_for_new = {}
    used_old = set()
    for overlap, old_i, new_i in candidates:
        if old_i in used_old or new_i in old_for_new:
            continue
        smaller = min(len(previous[old_i]["nodes"]), len(current[new_i]["nodes"]))
        if overlap >= max(1, smaller // 2):
            old_for_new[new_i] = old_i
            used_old.add(old_i)
    return old_for_new


def lake_volume(fields):
    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    return float(np.dot(fields.lake_depth, hypsometry.area_m2))


def step_record(fields, step, runtime_s, previous_bodies):
    bodies = components(fields)
    old_for_new = match_bodies(previous_bodies, bodies)
    reversals = 0
    eligible = 0
    for new_i, old_i in old_for_new.items():
        old, new = previous_bodies[old_i], bodies[new_i]
        if old["area_km2"] <= 10_000.0 or new["area_km2"] <= 10_000.0:
            continue
        eligible += 1
        prior_delta = old.get("delta_m3")
        delta = new["volume_m3"] - old["volume_m3"]
        if prior_delta is not None and prior_delta * delta < 0.0:
            baseline = max(old["volume_m3"], 1.0)
            if abs(prior_delta) > 0.5 * baseline and abs(delta) > 0.5 * baseline:
                reversals += 1
        new["delta_m3"] = delta

    hypsometry = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area = hypsometry.area_m2
    wet = fields.lake_depth > visible_depth
    ledger = getattr(fields, "lake_water_ledger", None)
    if ledger is not None:
        ledger_data = {
            name: float(getattr(ledger, name))
            for name in (
                "external_supply_m3", "wet_precip_m3", "dry_runoff_m3", "meltwater_m3",
                "evaporation_m3", "overflow_to_ocean_m3", "overflow_out_m3",
                "storage_change_m3", "internal_transfer_m3", "residual_m3",
            )
        }
        ledger_data["normalized_residual"] = abs(ledger.residual_m3) / max(ledger.external_supply_m3, 1.0)
        unbooked, capped = 0.0, 0.0
    else:
        unbooked, capped = legacy_surplus(fields)
        ledger_data = None

    drains = {
        lake.sink_node_idx for lake in lakes.iter_all_lakes(getattr(fields, "lake_forest", []))
        if not lake.children and lake.drains_at_floor
    }
    prior_drains = getattr(step_record, "previous_drains", set())
    step_record.previous_drains = drains
    notch = getattr(fields, "breach_notch_m", np.zeros(0))
    return {
        "step": step,
        "elapsed_years": world.elapsed_years,
        "runtime_s": runtime_s,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (2**20 if sys.platform == "darwin" else 2**10),
        "lake_volume_km3": lake_volume(fields) / 1e9,
        "lake_area_km2": float(area[wet].sum() / 1e6),
        "lake_body_count": len(bodies),
        "merge_events": sum(event.kind == "merge" for event in (getattr(fields, "lake_events", []) or [])),
        "split_events": sum(event.kind == "split" for event in (getattr(fields, "lake_events", []) or [])),
        "large_body_flicker_reversals": reversals,
        "large_body_flicker_comparable": eligible,
        "breached_pits": int(np.count_nonzero(notch > 0.0)) if len(notch) == len(fields.elevation) else 0,
        "floor_drain_leaves": len(drains),
        "floor_drain_toggles": len(drains.symmetric_difference(prior_drains)),
        "unbooked_balance_surplus_m3": unbooked,
        "surplus_held_at_cap_m3": capped,
        "lake_ledger": ledger_data,
    }, bodies


result.update(elapsed_years=world.elapsed_years, nodes=len(world.hydrology_cache.elevation), steps=[])
args.out.parent.mkdir(parents=True, exist_ok=True)
fields = world.hydrology_cache
previous_bodies = components(fields)
print(json.dumps({key: value for key, value in result.items() if key != "steps"}), flush=True)

for step in range(1, int(round(args.myr * 1e6 / args.step_years)) + 1):
    balance_calls.clear()
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    runtime_s = time.perf_counter() - started
    fields = world.hydrology_cache
    record, previous_bodies = step_record(fields, step, runtime_s, previous_bodies)
    result["steps"].append(record)
    print(json.dumps(record), flush=True)
    args.out.write_text(json.dumps(result, indent=1))

