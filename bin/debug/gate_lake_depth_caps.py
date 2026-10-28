#!/usr/bin/env python3
"""RFC 0014 G4: record depth-cap binding by lake body during a saved-world replay.

The probe wraps the production prognostic body integrator. An event is a lake body in a
step whose artificial lake/sea tier depth cap returns positive overflow. Rim-limited
spilling and floor drains are excluded. Event surplus is the volume routed out while
the body is held at its active tier cap; the record also includes the body's external
source supply during the step.

Run from ``backend/`` with the issue worktree's environment:

    .venv/bin/python ../bin/debug/gate_lake_depth_caps.py <save> --myr 1 \\
        --step-years 100000 --out <scratch>/seed.json
"""

from __future__ import annotations

import json
import time

import numpy as np
from replay_lib import load_replay, parse_replay_args

args = parse_replay_args(__doc__)
world, result = load_replay(args, issue=553)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import lakes  # noqa: E402

world.ice_age_period_years = 0.0
records: list[dict] = []
active_step = 0
active_context: dict | None = None
original_advance = lakes._advance_prognostic_body
original_integrate = lakes._integrate_body_volume


def tap_advance(body, dt_step, rates, climate, evaporation, sublimation, depth,
                prev_lake_depth, hypsometry, node_area_km2, node_count, is_ocean,
                accounting, leaves, body_leaf_ids, leaf_of_node, body_for_leaf,
                pending, processed):
    global active_context
    is_sea, tier_depth = lakes._classify_tier(body.members, prev_lake_depth, node_area_km2)
    rim_level = body.max_depth if body.max_depth is not None else np.inf
    cap_level = min(rim_level, body.floor_elevation + tier_depth)
    artificial = cap_level < rim_level - 1e-6
    if body.max_depth is None:
        artificial = True
    sinks = sorted(int(x) for x in lakes.leaf_sinks(body))
    area_km2 = float(hypsometry.area_m2[body.members].sum() / 1e6)
    supply_m3 = float(np.maximum(np.asarray(rates)[body.members], 0.0).sum()) * dt_step
    prior = active_context
    active_context = {
        "step": active_step,
        "sinks": sinks,
        "tier": "sea" if is_sea else "lake",
        "tier_cap_m": float(tier_depth),
        "cap_level_m": float(cap_level),
        "area_km2": area_km2,
        "member_count": int(len(body.members)),
        "supply_m3": supply_m3,
        "dt_years": float(dt_step),
        "artificial_cap": artificial and not body.drains_at_floor,
    }
    try:
        return original_advance(
            body, dt_step, rates, climate, evaporation, sublimation, depth,
            prev_lake_depth, hypsometry, node_area_km2, node_count, is_ocean,
            accounting, leaves, body_leaf_ids, leaf_of_node, body_for_leaf,
            pending, processed,
        )
    finally:
        active_context = prior


def tap_integrate(available, members, cap_level, dry_level, dt_step, evaporation,
                  sublimation, hypsometry):
    output = original_integrate(
        available, members, cap_level, dry_level, dt_step, evaporation,
        sublimation, hypsometry,
    )
    volume, evap_volume, overflow, level = output
    if active_context and active_context["artificial_cap"] and overflow > 0.0:
        records.append({
            **active_context,
            "surplus_m3": float(overflow),
            "surplus_m3_per_year": float(overflow / active_context["dt_years"])
            if active_context["dt_years"] > 0.0 else None,
            "evaporation_m3": float(evap_volume),
            "stored_volume_m3": float(volume),
            "surface_level_m": float(level),
        })
    return output


lakes._advance_prognostic_body = tap_advance
lakes._integrate_body_volume = tap_integrate

steps = int(round(args.myr * 1e6 / args.step_years))
result.update(
    elapsed_years=world.elapsed_years,
    nodes=len(world.hydrology_cache.elevation),
    replay_commit="70b3e904899f22f4cda2d3134ca2e8fb8e898070",
    step_years=args.step_years,
    steps=[],
    cap_events=[],
)
args.out.parent.mkdir(parents=True, exist_ok=True)
print(json.dumps({key: value for key, value in result.items() if key != "cap_events"}), flush=True)

for active_step in range(1, steps + 1):
    first_record = len(records)
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    step_records = records[first_record:]
    by_body: dict[tuple, dict] = {}
    for event in step_records:
        key = (event["tier"], tuple(event["sinks"]))
        row = by_body.setdefault(key, {
            "step": active_step,
            "tier": event["tier"],
            "sinks": event["sinks"],
            "tier_cap_m": event["tier_cap_m"],
            "cap_level_m": event["cap_level_m"],
            "area_km2": event["area_km2"],
            "member_count": event["member_count"],
            "source_supply_m3": 0.0,
            "surplus_m3": 0.0,
            "evaporation_m3": 0.0,
            "integration_calls": 0,
        })
        row["source_supply_m3"] += event["supply_m3"]
        row["surplus_m3"] += event["surplus_m3"]
        row["evaporation_m3"] += event["evaporation_m3"]
        row["integration_calls"] += 1
    for row in by_body.values():
        row["surplus_m3_per_year"] = row["surplus_m3"] / args.step_years
        result["cap_events"].append(row)
    step_record = {
        "step": active_step,
        "elapsed_years": world.elapsed_years,
        "runtime_s": time.perf_counter() - started,
        "binding_bodies": len(by_body),
        "lake_bodies": sum(row["tier"] == "lake" for row in by_body.values()),
        "sea_bodies": sum(row["tier"] == "sea" for row in by_body.values()),
        "surplus_m3": float(sum(row["surplus_m3"] for row in by_body.values())),
    }
    result["steps"].append(step_record)
    print(json.dumps(step_record), flush=True)
    args.out.write_text(json.dumps(result, indent=1))

summary: dict[str, dict] = {}
for row in result["cap_events"]:
    key = row["tier"]
    item = summary.setdefault(key, {"binding_body_steps": 0, "surplus_m3": 0.0, "sinks": set()})
    item["binding_body_steps"] += 1
    item["surplus_m3"] += row["surplus_m3"]
    item["sinks"].update(row["sinks"])
result["summary"] = {}
for key, value in summary.items():
    tier_events = [row for row in result["cap_events"] if row["tier"] == key]
    result["summary"][key] = {
        "binding_body_steps": value["binding_body_steps"],
        "steps_with_binding": len({row["step"] for row in tier_events}),
        "replay_steps": steps,
        "surplus_m3": value["surplus_m3"],
        "surplus_m3_per_year": value["surplus_m3"] / (args.myr * 1e6),
        "distinct_body_signatures": len({tuple(row["sinks"]) for row in tier_events}),
        "distinct_sinks": sorted(value["sinks"]),
    }
result["cap_events"] = sorted(
    result["cap_events"], key=lambda row: (row["step"], row["tier"], row["sinks"])
)
args.out.write_text(json.dumps(result, indent=1))
print(json.dumps({"summary": result["summary"], "cap_event_body_steps": len(result["cap_events"])}))
