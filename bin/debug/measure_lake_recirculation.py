#!/usr/bin/env python3
"""Measure repeated internal transfers in the lake solve's flush passes.

Replay a save unperturbed and, each step, count the solve's passes against its limit of one
plus the leaf count, list the body pairs that send most, and record the surplus still deferred
when the loop ends, which `solve_lakes_prognostic` exports at a leaf sink. Changes no model
term.

Example from the repository root:

    backend/.venv/bin/python bin/debug/measure_lake_recirculation.py \\
      /tmp/seed910211954-221.6myr.mbworld --myr 0.3 --out /tmp/recirculation-a.json
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from replay_lib import load_replay, parse_replay_args

TOP_PAIRS = 6
TOP_EXPORTS = 6

args = parse_replay_args(__doc__)
world, header = load_replay(args, issue=691)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import lakes  # noqa: E402

pairs: dict[tuple[int, int], dict] = {}
exports: list[dict] = []
pass_counts = {"active_bodies": 0}
last_is_ocean: list[np.ndarray] = []

trace_transfer = lakes._trace_transfer
trace_deferred_export = lakes._trace_deferred_export
active_bodies = lakes._active_bodies
solve_lakes_prognostic = lakes.solve_lakes_prognostic


def body_state(body) -> dict:
    tier_floor = float(body.floor_elevation)
    return {
        "lake_id": int(body.lake_id),
        "members": len(body.members),
        "floor_m": tier_floor,
        "rim_m": None if body.max_depth is None else float(body.max_depth),
        "rim_depth_m": None if body.max_depth is None else float(body.max_depth) - tier_floor,
        "level_m": float(body.current_water_elevation),
        "drains_at_floor": bool(body.drains_at_floor),
        "outlet_target": int(body.outlet_target_idx),
    }


def observing_transfer(accounting, source, receiver, target, volume, *, internal, depth):
    if internal and receiver is not None:
        key = (id(source), id(receiver))
        row = pairs.setdefault(key, {"sends": 0, "total_m3": 0.0, "first_m3": float(volume)})
        row.update(sends=row["sends"] + 1, total_m3=row["total_m3"] + float(volume), last_m3=float(volume))
        row["source"], row["receiver"] = source, receiver
    return trace_transfer(accounting, source, receiver, target, volume, internal=internal, depth=depth)


def observing_export(accounting, leaf, sink, volume, depth):
    is_ocean = last_is_ocean[0] if last_is_ocean else None
    exports.append({
        "leaf_id": int(leaf.lake_id),
        "sink": int(sink),
        "volume_m3": float(volume),
        "sink_is_ocean": bool(is_ocean[sink]) if is_ocean is not None and 0 <= sink < len(is_ocean) else None,
    })
    return trace_deferred_export(accounting, leaf, sink, volume, depth)


def observing_active_bodies(*body_args, **body_kwargs):
    pass_counts["active_bodies"] += 1
    return active_bodies(*body_args, **body_kwargs)


def observing_solve(*solve_args, **solve_kwargs):
    last_is_ocean[:] = [np.asarray(solve_kwargs["is_ocean"])] if solve_kwargs.get("is_ocean") is not None else []
    return solve_lakes_prognostic(*solve_args, **solve_kwargs)


lakes._trace_transfer = observing_transfer
lakes._trace_deferred_export = observing_export
lakes._active_bodies = observing_active_bodies
lakes.solve_lakes_prognostic = observing_solve


def snapshot(step: int) -> dict:
    fields = world.hydrology_cache
    ledger = fields.lake_water_ledger
    leaves = sum(1 for lake in lakes.iter_all_lakes(fields.lake_forest or []) if not lake.children)
    ranked = sorted(pairs.values(), key=lambda row: -row["total_m3"])
    reciprocal = {(id(row["source"]), id(row["receiver"])) for row in pairs.values()}
    grouped: dict[str, float] = defaultdict(float)
    for row in pairs.values():
        kind = "reciprocal" if (id(row["receiver"]), id(row["source"])) in reciprocal else "one_way"
        grouped[kind] += row["total_m3"]
    # `_active_bodies` runs once before and once after each pass.
    passes = pass_counts["active_bodies"] // 2
    return {
        "step": step,
        "age_myr": float(world.elapsed_years) / 1e6,
        "leaves": leaves,
        "passes": passes,
        "pass_limit": 1 + leaves,
        "hit_limit": passes >= 1 + leaves,
        "internal_transfer_m3": float(ledger.internal_transfer_m3),
        "external_supply_m3": float(ledger.external_supply_m3),
        "overflow_to_ocean_m3": float(ledger.overflow_to_ocean_m3),
        "transfer_m3_by_kind": dict(grouped),
        "top_pairs": [
            {
                "sends": row["sends"],
                "total_m3": row["total_m3"],
                "first_m3": row["first_m3"],
                "last_m3": row["last_m3"],
                "reciprocal": (id(row["receiver"]), id(row["source"])) in reciprocal,
                "source": body_state(row["source"]),
                "receiver": body_state(row["receiver"]),
            }
            for row in ranked[:TOP_PAIRS]
        ],
        "final_export_m3": float(sum(item["volume_m3"] for item in exports)),
        "final_exports": sorted(exports, key=lambda item: -item["volume_m3"])[:TOP_EXPORTS],
    }


world.ice_age_period_years = 0.0
args.out.parent.mkdir(parents=True, exist_ok=True)
result = {**header, "records": []}
for step in range(1, int(round(args.myr * 1e6 / args.step_years)) + 1):
    pairs.clear()
    exports.clear()
    pass_counts.update(active_bodies=0)
    world_mod.step_world(world, args.step_years)
    result["records"].append(snapshot(step))
    args.out.write_text(json.dumps(result, indent=1) + "\n")
print(f"wrote {args.out}")
