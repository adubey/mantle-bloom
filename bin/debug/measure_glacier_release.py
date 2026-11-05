#!/usr/bin/env python3
"""RFC 0016 G0: explain the interim glacier account's release from the model ice's own budget.

The account releases the fraction of opening model ice lost each step, so its size follows the
ice model's budget. This replays a save and records the raw `IceLedger`, splits calving by the
branch that produced it, and locates each step's calving.

Example:

    backend/.venv/bin/python bin/debug/measure_glacier_release.py <checkpoint.mbworld> \\
        --myr 0.4 --step-years 100000 --checkpoint-myr 0.4 --out <scratch>/glacier-release-a.json
"""

from __future__ import annotations

from dataclasses import asdict

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay


def add_arguments(ap) -> None:
    ap.add_argument("--keep-ice-ages", action="store_true", help="keep the save's ice-age cycle")


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=725)
    if not args.keep_ice_ages:
        world.ice_age_period_years = 0.0
    result["ice_age_period_years"] = world.ice_age_period_years
    from app.hydroclimate import ice_sheet  # noqa: E402

    result["steps"] = []
    branch = {"saturated_excess_m3": 0.0, "final_afloat_m3": 0.0}
    original_shed, original_guard = ice_sheet._shed_excess, ice_sheet._guard_body

    def shed(out, body, excess, margin, marine, area):
        before = body.calved_m3
        original_shed(out, body, excess, margin, marine, area)
        branch["saturated_excess_m3"] += body.calved_m3 - before

    def guard(out, body, bed, area, sea_level_m, is_ocean):
        before = body.calved_m3
        original_guard(out, body, bed, area, sea_level_m, is_ocean)
        branch["final_afloat_m3"] += body.calved_m3 - before

    ice_sheet._shed_excess, ice_sheet._guard_body = shed, guard
    original_step = ice_sheet.step_ice
    stray_anatomy = {}

    def step_ice(points, neighbor_idx, area, bed, deflection, sea_level_m, prev_depth, climate, cover, is_ocean, *rest, **kw):
        h0 = np.clip(prev_depth, 0.0, None)
        real_bed = bed - deflection
        ice_free = (h0 > 0.0) & ~ice_sheet.is_grounded(h0, real_bed, sea_level_m, is_ocean)
        volume = h0 * area
        stray_anatomy.update(
            ungrounded_nodes=int(ice_free.sum()),
            ungrounded_m3=float(volume[ice_free].sum()),
            ungrounded_land_nodes=int((ice_free & ~is_ocean).sum()),
            ungrounded_land_m3=float(volume[ice_free & ~is_ocean].sum()),
            ungrounded_ocean_m3=float(volume[ice_free & is_ocean].sum()),
            ungrounded_median_ice_m=float(np.median(h0[ice_free])) if ice_free.any() else 0.0,
            ungrounded_median_bed_m=float(np.median(real_bed[ice_free])) if ice_free.any() else 0.0,
            sea_level_m=float(sea_level_m),
            land_below_sea_level_nodes=int(((real_bed < sea_level_m) & ~is_ocean).sum()),
        )
        return original_step(points, neighbor_idx, area, bed, deflection, sea_level_m, prev_depth, climate, cover, is_ocean, *rest, **kw)

    ice_sheet.step_ice = step_ice

    def snapshot() -> dict:
        return {"elapsed_myr": world.elapsed_years / 1e6}

    def after_step(_seconds: float) -> None:
        hydro = world.hydrology_cache
        ice = hydro.ice
        ledger = asdict(ice.ledger)
        calved = np.asarray(ice.calved_m3)
        total = float(calved.sum())
        order = np.argsort(-calved)
        cumulative = np.cumsum(calved[order])
        hit = calved > 0.0
        result["steps"].append({
            "step": len(result["steps"]) + 1,
            "ledger": ledger,
            "ice_residual_m3": ice.ledger.residual_m3,
            "calving_branches_m3": {
                **branch, "stray_at_step_start_m3": total - sum(branch.values()),
            },
            "calving_nodes": int(hit.sum()),
            "nodes_for_half_of_calving": int(np.searchsorted(cumulative, total / 2.0)) + 1 if total else 0,
            "calving_on_ocean_share": float(calved[hydro.is_ocean].sum()) / max(total, 1.0),
            "largest_calving_point": hydro.points[int(order[0])].tolist(),
            "sea_level_m": world.sea_level_m,
            "bodies": len(ice.bodies),
            "stray_anatomy": dict(stray_anatomy),
            "largest_bodies_calved_m3": sorted((b.calved_m3 for b in ice.bodies), reverse=True)[:5],
        })
        branch.update(saturated_excess_m3=0.0, final_afloat_m3=0.0)

    run_replay(world, args, result, snapshot, after_step, lambda _steps: {})


if __name__ == "__main__":
    main()
