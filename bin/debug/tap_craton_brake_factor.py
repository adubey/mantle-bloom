#!/usr/bin/env python3
"""How much RFC 0008 solution 2's craton proximity factor raises the contact brake (#423).

Replays one seed from generation at fixed steps and taps `torque.collision_drag_coefficients`
in this process only. At every call it compares the brake factor g with the proximity factor
c (`torque.contact_drag_pa_s_per_m` with and without `craton`) and counts, over the run:

- band: boundary-band nodes (within `BOUNDARY_FORCE_REACH_MULTIPLIER` spacings of a neighbour);
- braked: nodes with g > 0 (contested or jam-contact nodes, the only nodes c can reach);
- c_pos: braked-eligible nodes with a craton in reach (c > 0);
- c_raises: nodes where c is larger than max(j, relief), so it changes the brake;
- g_sum_without / g_sum_with: the sum of g over all nodes without and with c.

Needs the solution 2 code (branch claude/issue-423-proximity-brake); on main, the brake takes
craton thickness, not c.

    backend/.venv/bin/python bin/debug/tap_craton_brake_factor.py --seed 7 --years 3e8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import lithosphere, torque  # noqa: E402
from app import world as world_mod  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--years", type=float, default=3e8)
    parser.add_argument("--step", type=float, default=5e6)
    args = parser.parse_args()

    counts = dict(band=0, braked=0, c_pos=0, c_raises=0, g_sum_without=0.0, g_sum_with=0.0)
    original = torque.collision_drag_coefficients

    def tapped(plate, inputs, mask, spacing_rad, craton=None, jam=None):
        reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad
        unit = torque.CONTACT_RESISTANCE_N_PER_M / (
            reach * lithosphere.PLANET_RADIUS_M * torque.CONTACT_RESISTANCE_REF_CM_PER_YR * 0.01 / torque.SECONDS_PER_YEAR
        )
        with_c = torque.contact_drag_pa_s_per_m(plate, inputs, mask, spacing_rad, craton, jam) / unit
        without = torque.contact_drag_pa_s_per_m(plate, inputs, mask, spacing_rad, None, jam) / unit
        counts["band"] += int(np.sum(inputs.dist_to_neighbor <= reach))
        counts["braked"] += int(np.sum(with_c > 0.0))
        if craton is not None:
            counts["c_pos"] += int(np.sum(craton > 0.0))
            counts["c_raises"] += int(np.sum(with_c > without + 1e-12))
        counts["g_sum_without"] += float(without.sum())
        counts["g_sum_with"] += float(with_c.sum())
        return original(plate, inputs, mask, spacing_rad, craton, jam)

    torque.collision_drag_coefficients = tapped
    world = world_mod.generate_world(args.seed)
    while world.elapsed_years < args.years - 1.0:
        world_mod.step_world(world, min(args.step, args.years - world.elapsed_years))
    print({k: round(v, 1) if isinstance(v, float) else v for k, v in counts.items()})
    print(
        f"braked share of band {counts['braked'] / max(counts['band'], 1):.1%}; "
        f"mean g at braked nodes {counts['g_sum_without'] / max(counts['braked'], 1):.2f}; "
        f"c raises g at {counts['c_raises'] / max(counts['braked'], 1):.1%} of braked nodes; "
        f"total g +{counts['g_sum_with'] / max(counts['g_sum_without'], 1e-30) - 1:.1%}"
    )


if __name__ == "__main__":
    main()
