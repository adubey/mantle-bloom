#!/usr/bin/env python3
"""Issue #632: production-code RFC 0014 G5 lake-ice replay.

Advance one preregistered diagnosis save through the full production world step and record
the lake ice cap, bedfast-area split, liquid-under-ice storage, lake water ledger, runtime,
and peak memory at each checkpoint. Unlike gate_lake_ice_cap.py, this does not hold climate,
terrain, or lake depth fixed: it calls `world.step_world` and reads the resulting production
state. Never use an RFC 0014 held-out save here.

Run one save at a time, and never alongside a replay or test:

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_production.py \\
        <local-save> --myr 10 --step-years 100000 --checkpoint-myr 0.1 \\
        --out analysis/issue632/production-g5-seed<seed>.json
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

args = parse_replay_args(
    __doc__, lambda ap: ap.set_defaults(checkpoint_myr=0.1)
)

from app import persistence, plates, world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from lake_ice_snapshot import lake_ice_snapshot  # noqa: E402

DIAGNOSIS_SHA256 = {
    "bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2": 910211954,
    "ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045": 635056485,
    "6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd": 317739385,
}
RHO_WATER_TO_ICE = 1000.0 / 917.0
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
if save_sha256 not in DIAGNOSIS_SHA256:
    raise ValueError("Production G5 accepts only the three preregistered diagnosis saves")
world = persistence.load_world_bytes(save_bytes)
if world.seed != DIAGNOSIS_SHA256[save_sha256]:
    raise ValueError("save checksum is registered for a different seed")
world.ice_age_period_years = 0.0

args.out.parent.mkdir(parents=True, exist_ok=True)
steps = int(round(args.myr * 1_000_000.0 / args.step_years))
checkpoint_steps = max(1, int(round(args.checkpoint_myr * 1_000_000.0 / args.step_years)))
result = {
    "issue": 632,
    "rfc": "0014",
    "gate": "G5-production",
    "backend_commit": subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip(),
    "save": args.save.name,
    "save_sha256": save_sha256,
    "seed": int(world.seed),
    "start_age_myr": float(world.elapsed_years / 1_000_000.0),
    "duration_myr": float(args.myr),
    "step_years": float(args.step_years),
    "ice_age_period_years": float(world.ice_age_period_years),
    "checkpoints": [],
}

ledger_totals = {name: 0.0 for name in LEDGER_TERMS}
max_step_residual_fraction = 0.0
zero_supply_steps = 0
step_seconds = []


def snapshot(time_myr: float) -> dict:
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
    physical_ice_m = ice_we_m * RHO_WATER_TO_ICE
    wet_area_m2 = float(areas_m2[wet].sum())
    deep_area_m2 = float(areas_m2[deep].sum())
    if np.any(wet):
        thickness = physical_ice_m[wet]
        thickness_stats = {
            f"p{p}_m": round(float(np.percentile(thickness, p)), 6)
            for p in (50, 90, 99)
        }
        thickness_stats["max_m"] = round(float(np.max(thickness, initial=0.0)), 6)
    else:
        thickness_stats = {"p50_m": 0.0, "p90_m": 0.0, "p99_m": 0.0, "max_m": 0.0}
    ledger = hydro.lake_water_ledger if hydro is not None else None
    return {
        "time_myr": round(time_myr, 4),
        "wet_lake_nodes": int(wet.sum()),
        "wet_lake_area_km2": round(wet_area_m2 / 1e6, 3),
        "liquid_volume_km3": round(float(np.dot(liquid_m[wet], areas_m2[wet])) / 1e9, 6),
        "ice_water_equivalent_volume_km3": round(float(np.dot(ice_we_m[wet], areas_m2[wet])) / 1e9, 6),
        "liquid_under_ice_area_pct": round(
            100.0 * float(areas_m2[wet & (ice_we_m > 1e-6) & (liquid_m > 1e-6)].sum())
            / wet_area_m2, 6
        ) if wet_area_m2 else 0.0,
        "physical_ice_thickness": thickness_stats,
        "max_water_equivalent_excess_m": round(
            float(np.max(ice_we_m[wet] - total_water_m[wet], initial=0.0)), 12
        ) if np.any(wet) else 0.0,
        "shallow_bedfast_area_pct_of_shallow_lake_area": round(
            100.0 * float(areas_m2[bedfast & shallow].sum())
            / float(areas_m2[shallow].sum()), 6
        ) if np.any(shallow) else 0.0,
        "deep_bedfast_area_pct_of_deep_lake_area": round(
            100.0 * float(areas_m2[bedfast & deep].sum()) / deep_area_m2, 6
        ) if deep_area_m2 else 0.0,
        "ledger_last_step": {
            name: float(getattr(ledger, name)) for name in LEDGER_TERMS
        } if ledger is not None else None,
        "ledger_cumulative_m3": dict(ledger_totals),
        "max_step_residual_pct_of_external_supply": round(
            max_step_residual_fraction * 100.0, 9
        ),
        "zero_supply_steps": zero_supply_steps,
        "mean_step_seconds_since_checkpoint": round(
            float(np.mean(step_seconds[-checkpoint_steps:])) if step_seconds else 0.0, 4
        ),
        "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / RSS_PER_MB, 1),
    }


result["checkpoints"].append(snapshot(0.0))
print(json.dumps(result["checkpoints"][-1]), flush=True)
for step in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    elapsed = time.perf_counter() - started
    step_seconds.append(elapsed)
    ledger = world.hydrology_cache.lake_water_ledger if world.hydrology_cache is not None else None
    if ledger is not None:
        for name in LEDGER_TERMS:
            ledger_totals[name] += float(getattr(ledger, name))
        supply = float(ledger.external_supply_m3)
        residual = abs(float(ledger.residual_m3))
        if supply > 0.0:
            max_step_residual_fraction = max(max_step_residual_fraction, residual / supply)
        else:
            zero_supply_steps += 1
    if step % checkpoint_steps == 0 or step == steps:
        point = snapshot((world.elapsed_years - result["start_age_myr"] * 1_000_000.0) / 1_000_000.0)
        result["checkpoints"].append(point)
        args.out.write_text(json.dumps(result, indent=1))
        print(json.dumps(point), flush=True)

result["step_runtime_seconds"] = {
    "total": round(float(sum(step_seconds)), 3),
    "median": round(float(np.median(step_seconds)), 4),
    "p90": round(float(np.percentile(step_seconds, 90)), 4),
}
result["plateau_change_0p1_to_10myr_pct"] = None
if args.myr >= 0.1 and len(result["checkpoints"]) >= 3:
    at_0p1 = result["checkpoints"][1]["ice_water_equivalent_volume_km3"]
    at_end = result["checkpoints"][-1]["ice_water_equivalent_volume_km3"]
    result["plateau_change_0p1_to_10myr_pct"] = round(
        100.0 * (at_end - at_0p1) / at_0p1, 6
    ) if at_0p1 else None
args.out.write_text(json.dumps(result, indent=1))
print(f"wrote {args.out}", flush=True)
