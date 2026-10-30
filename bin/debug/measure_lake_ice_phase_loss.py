#!/usr/bin/env python3
"""Reconcile the liquid lake water `advance_lake_ice` removes in one step (#691).

Step a save once, and compare the solver's liquid depth with the depth after the phase update.
Reports the volume lost, how much sits on silted cells left at exactly zero, and the ice budget
that books it. Changes no model term.

Example from the repository root:

    backend/.venv/bin/python bin/debug/measure_lake_ice_phase_loss.py \\
      /tmp/seed910211954-221.6myr.mbworld --myr 0.1 --out /tmp/issue691-ice-phase.json
"""

from __future__ import annotations

import json

import numpy as np

from replay_lib import load_replay, parse_replay_args

args = parse_replay_args(__doc__)
world, header = load_replay(args, issue=691)

from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology  # noqa: E402

captured: dict = {}
advance_lake_ice = hydrology.advance_lake_ice


def observing_advance_lake_ice(prior_ice, liquid_depth, *rest, **kwargs):
    result = advance_lake_ice(prior_ice, liquid_depth, *rest, **kwargs)
    captured.update(solved=np.array(liquid_depth), ice=np.array(result[0]), budget=result[2])
    return result


hydrology.advance_lake_ice = observing_advance_lake_ice
world.ice_age_period_years = 0.0
world_mod.step_world(world, args.step_years)

fields = world.hydrology_cache
area = np.asarray(fields.hypsometry.area_m2)
solved, final = captured["solved"], np.asarray(fields.lake_depth)
silted_to_zero = (np.asarray(fields.silt_deposited) > 0.0) & (final <= 0.0)
budget = captured["budget"]
result = {
    **header,
    "solved_liquid_km3": float(np.dot(solved, area)) / 1e9,
    "final_liquid_km3": float(np.dot(final, area)) / 1e9,
    "liquid_lost_km3": float(np.dot(np.clip(solved - final, 0.0, None), area)) / 1e9,
    "final_ice_km3": float(np.dot(captured["ice"], area)) / 1e9,
    "silted_cells_at_zero": int(silted_to_zero.sum()),
    "silted_cells_at_zero_solved_km3": float(np.dot(solved[silted_to_zero], area[silted_to_zero])) / 1e9,
    "sublimation_km3": float(budget["sublimation_m3"]) / 1e9,
    "freezing_km3": float(budget["freezing_m3"]) / 1e9,
    "melting_km3": float(budget["melting_m3"]) / 1e9,
}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(result, indent=1) + "\n")
print(json.dumps(result, indent=1))
