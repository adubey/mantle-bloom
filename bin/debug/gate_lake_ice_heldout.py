#!/usr/bin/env python3
"""Issue #632: held-out RFC 0014 lake-ice acceptance replay.

Advance one preregistered held-out save through ten production 100 kyr world steps.
Record the combined liquid/ice ledger and local lake-ice bounds without tuning. This driver
accepts only the three held-out saves registered on RFC 0014 issue #538.

Follow the compute-budget skill before starting each save:

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_heldout.py \
        <held-out-save> --myr 1 --step-years 100000 --checkpoint-myr 0.1 \
        --out <scratch>/heldout-seed<seed>.json
"""

from __future__ import annotations

import hashlib
import json
import resource
import subprocess
import sys
import time

import numpy as np
from replay_lib import parse_replay_args

args = parse_replay_args(__doc__, lambda ap: ap.set_defaults(myr=1.0, checkpoint_myr=0.1))

from app import persistence, plates, world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from lake_ice_snapshot import lake_ice_snapshot  # noqa: E402

HELDOUT_SHA256 = {
    "9980103ce9bd7733e4671bc3bec0d48359b45a8e457bc76a73dd769f8840f342": 373936961,
    "a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7": 997271774,
    "e6a6c7514ae52601d8d2e3ae541da3c1dcb727cde141a79b0bee39d3c962e865": 456862899,
}
ICE_TO_WATER_DENSITY = 1000.0 / 917.0
RSS_PER_MB = 2**20 if sys.platform == "darwin" else 2**10
LEDGER_TERMS = (
    "external_supply_m3",
    "evaporation_m3",
    "overflow_out_m3",
    "storage_change_m3",
    "sublimation_m3",
    "ice_storage_change_m3",
    "freezing_m3",
    "melting_m3",
    "residual_m3",
)

save_bytes = args.save.read_bytes()
save_sha256 = hashlib.sha256(save_bytes).hexdigest()
if save_sha256 not in HELDOUT_SHA256:
    raise ValueError("Held-out acceptance accepts only the three preregistered RFC 0014 saves")
world = persistence.load_world_bytes(save_bytes)
if world.seed != HELDOUT_SHA256[save_sha256]:
    raise ValueError("save checksum is registered for a different seed")
world.ice_age_period_years = 0.0

steps = int(round(args.myr * 1_000_000.0 / args.step_years))
if steps != 10 or not np.isclose(args.step_years, 100_000.0):
    raise ValueError("RFC 0014 held-out ledger acceptance is exactly 10 steps of 100 kyr")
checkpoint_steps = max(1, int(round(args.checkpoint_myr * 1_000_000.0 / args.step_years)))
if checkpoint_steps != 1:
    raise ValueError("Held-out acceptance records every 100 kyr step")

result = {
    "issue": 632,
    "rfc": "0014",
    "gate": "held-out-lake-ice-acceptance",
    "backend_commit": subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip(),
    "save": args.save.name,
    "save_sha256": save_sha256,
    "seed": int(world.seed),
    "start_age_myr": float(world.elapsed_years / 1_000_000.0),
    "step_years": float(args.step_years),
    "steps_requested": steps,
    "ice_age_period_years": float(world.ice_age_period_years),
    "steps": [],
}
ledger_totals = {name: 0.0 for name in LEDGER_TERMS}
max_residual_fraction = 0.0
zero_supply_steps = 0
step_seconds = []


def snapshot(step: int, runtime_s: float | None) -> dict:
    state = lake_ice_snapshot(world, plates, line_spacing_rad)
    hydro = state["hydro"]
    liquid_m = state["liquid_m"]
    ice_we_m = state["ice_we_m"]
    areas_m2 = state["areas_m2"]
    wet = state["wet"]
    shallow = state["shallow"]
    deep = state["deep"]
    bedfast = state["bedfast"]
    total_water_m = liquid_m + ice_we_m
    ice_m = ice_we_m * ICE_TO_WATER_DENSITY
    wet_area = float(areas_m2[wet].sum())
    shallow_area = float(areas_m2[shallow].sum())
    deep_area = float(areas_m2[deep].sum())
    under_ice_area = float(areas_m2[wet & (ice_we_m > 1e-6) & (liquid_m > 1e-6)].sum())
    ledger = hydro.lake_water_ledger if hydro is not None else None
    return {
        "step": step,
        "elapsed_years": float(world.elapsed_years),
        "runtime_seconds": runtime_s,
        "wet_lake_nodes": int(wet.sum()),
        "wet_lake_area_km2": wet_area / 1e6,
        "liquid_volume_km3": float(np.dot(liquid_m[wet], areas_m2[wet])) / 1e9,
        "ice_water_equivalent_volume_km3": float(np.dot(ice_we_m[wet], areas_m2[wet])) / 1e9,
        "max_physical_ice_m": float(np.max(ice_m[wet], initial=0.0)),
        "max_water_equivalent_excess_m": float(
            np.max(ice_we_m[wet] - total_water_m[wet], initial=0.0)
        ),
        "shallow_bedfast_area_pct": (
            100.0 * float(areas_m2[bedfast & shallow].sum()) / shallow_area if shallow_area else 0.0
        ),
        "deep_bedfast_area_pct": (
            100.0 * float(areas_m2[bedfast & deep].sum()) / deep_area if deep_area else 0.0
        ),
        "liquid_under_ice_area_pct": (
            100.0 * under_ice_area / wet_area if wet_area else 0.0
        ),
        "minimum_liquid_depth_m": float(np.min(liquid_m, initial=0.0)),
        "minimum_ice_water_equivalent_m": float(np.min(ice_we_m, initial=0.0)),
        "ledger": ({
            **{name: float(getattr(ledger, name)) for name in LEDGER_TERMS},
            "residual_pct_of_external_supply": (
                100.0 * abs(float(ledger.residual_m3)) / float(ledger.external_supply_m3)
                if ledger.external_supply_m3 > 0.0
                else None
            ),
        } if ledger is not None else None),
        "ledger_cumulative_m3": dict(ledger_totals),
        "zero_supply_steps": zero_supply_steps,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / RSS_PER_MB,
    }


result["initial"] = snapshot(0, None)
args.out.parent.mkdir(parents=True, exist_ok=True)
print(json.dumps(result["initial"]), flush=True)
for step in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    elapsed = time.perf_counter() - started
    step_seconds.append(elapsed)
    hydro = world.hydrology_cache
    ledger = hydro.lake_water_ledger if hydro is not None else None
    if ledger is not None:
        for name in LEDGER_TERMS:
            ledger_totals[name] += float(getattr(ledger, name))
        supply = float(ledger.external_supply_m3)
        residual = abs(float(ledger.residual_m3))
        if supply > 0.0:
            max_residual_fraction = max(max_residual_fraction, residual / supply)
        else:
            zero_supply_steps += 1
    record = snapshot(step, elapsed)
    record["max_step_residual_pct_of_external_supply"] = 100.0 * max_residual_fraction
    result["steps"].append(record)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(record), flush=True)

result["runtime_seconds"] = {
    "total": float(sum(step_seconds)),
    "median_step": float(np.median(step_seconds)),
    "p90_step": float(np.percentile(step_seconds, 90)),
}
result["max_step_residual_pct_of_external_supply"] = 100.0 * max_residual_fraction
result["acceptance"] = {
    "physical_ice_at_most_2m": all(step["max_physical_ice_m"] <= 2.0 + 1e-9 for step in result["steps"]),
    "nonnegative_reservoirs": all(
        step["minimum_liquid_depth_m"] >= -1e-9
        and step["minimum_ice_water_equivalent_m"] >= -1e-9
        for step in result["steps"]
    ),
    "ice_within_local_water_inventory": all(
        step["max_water_equivalent_excess_m"] <= 1e-9 for step in result["steps"]
    ),
    "deep_bedfast_area_at_most_1pct": all(
        step["deep_bedfast_area_pct"] <= 1.0 for step in result["steps"]
    ),
    "ledger_residual_at_most_0p1pct": max_residual_fraction <= 0.001,
    "zero_supply_steps": zero_supply_steps,
}
args.out.write_text(json.dumps(result, indent=1))
print(f"wrote {args.out}", flush=True)
