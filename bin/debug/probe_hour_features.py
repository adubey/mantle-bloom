#!/usr/bin/env python3
"""Explain single-step lake drops and discharge spikes in a fixed-step replay (#811).

Repeats `replay_timescale_sweep.py`'s warm-up and stepping (10 day-steps, then the requested
step), and at each `--watch` step records what changed: the lake ledger in km3, the lake merge
and split events, every matched lake body whose volume moved by more than `--min-km3` (with the
share of its former nodes that are ocean afterwards), and the three largest land discharges with
their cell, elevation and lake depth. A discharge times the step length is compared with the
world lake change to test whether a spike is a one-step lake release booked as a rate.

Run from the repository root:

    backend/.venv/bin/python bin/debug/probe_hour_features.py <save> --backend <backend> \\
        --step-years 0.00011408 --myr 1.4e-5 --watch 12 13 14 109 110 111 \\
        --out /tmp/hour-features.json
"""

from __future__ import annotations

import json

import numpy as np

from lake_body_tracking import components, match_bodies
from replay_lib import load_replay, parse_replay_args

DAY_YEARS = 1.0 / 365.25
M3_PER_KM3 = 1e9
SECONDS_PER_YEAR = 365.25 * 86400.0


def add_arguments(ap) -> None:
    """Add the probe's own arguments to the shared replay CLI."""
    ap.add_argument("--warmup-days", type=int, default=10)
    ap.add_argument("--watch", type=int, nargs="+", required=True, help="1-based measured steps to describe")
    ap.add_argument("--min-km3", type=float, default=1.0, help="smallest body volume change reported")
    ap.add_argument("--wind", choices=("diagnostic", "cfd"), default="diagnostic")


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=811)
    world.ice_age_period_years = 0.0
    world.wind_model = args.wind
    from app import world as world_mod  # noqa: E402

    for _ in range(args.warmup_days):
        world_mod.step_world(world, DAY_YEARS)
    previous = components(world.hydrology_cache)
    step_seconds = args.step_years * SECONDS_PER_YEAR
    steps = int(round(args.myr * 1e6 / args.step_years)) if args.myr >= 1e-3 else max(args.watch)
    records = []
    for step in range(1, max(args.watch) + 1):
        world_mod.step_world(world, args.step_years)
        fields = world.hydrology_cache
        current = components(fields)
        if step in args.watch:
            records.append(describe(step, step_seconds, fields, previous, current, args.min_km3))
            print(json.dumps(records[-1])[:400], flush=True)
        previous = current
    result.update(wind_model=args.wind, watch=records)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(f"wrote {args.out} ({steps} requested steps)")


def describe(step, step_seconds, fields, previous, current, min_km3) -> dict:
    """Ledger, topology events, moved bodies and top discharges for one step."""
    ledger = fields.lake_water_ledger
    ledger_km3 = {
        name: float(getattr(ledger, name)) / M3_PER_KM3
        for name in ("external_supply_m3", "evaporation_m3", "overflow_to_ocean_m3", "overflow_out_m3",
                     "storage_change_m3", "internal_transfer_m3", "sublimation_m3")
    }
    events = [
        {"kind": e.kind, "node_count": int(e.node_count), "elevation_m": float(e.elevation_m)}
        for e in (getattr(fields, "lake_events", []) or [])
    ]
    is_ocean = np.asarray(fields.is_ocean)
    moved = []
    matched = match_bodies(previous, current)
    for new_i, old_i in matched.items():
        old, new = previous[old_i], current[new_i]
        delta = (new["volume_m3"] - old["volume_m3"]) / M3_PER_KM3
        if abs(delta) >= min_km3:
            moved.append({
                "delta_km3": delta, "before_km3": old["volume_m3"] / M3_PER_KM3,
                "after_km3": new["volume_m3"] / M3_PER_KM3,
                "area_before_km2": old["area_km2"], "area_after_km2": new["area_km2"],
                "former_nodes_ocean_now": float(is_ocean[old["nodes"]].mean()),
            })
    gone = []
    for old_i, old in enumerate(previous):
        if old_i in matched.values():
            continue
        volume = old["volume_m3"] / M3_PER_KM3
        if volume >= min_km3:
            gone.append({
                "volume_km3": volume, "area_km2": old["area_km2"],
                "former_nodes_ocean_now": float(is_ocean[old["nodes"]].mean()),
            })
    discharge = np.where(is_ocean, 0.0, np.asarray(fields.discharge_m3_s))
    top = []
    for node in np.argsort(discharge)[-3:][::-1]:
        top.append({
            "node": int(node), "discharge_m3_s": float(discharge[node]),
            "flow_km3_over_step": float(discharge[node]) * step_seconds / M3_PER_KM3,
            "elevation_m": float(fields.elevation[node]), "lake_depth_m": float(fields.lake_depth[node]),
        })
    world_change = sum(body["delta_km3"] for body in moved) - sum(body["volume_km3"] for body in gone)
    return {
        "step": step, "ledger_km3": ledger_km3, "lake_events": events,
        "bodies_moved": sorted(moved, key=lambda b: -abs(b["delta_km3"]))[:6],
        "bodies_gone": sorted(gone, key=lambda b: -b["volume_km3"])[:6],
        "net_matched_and_gone_change_km3": world_change, "top_discharges": top,
    }


if __name__ == "__main__":
    main()
