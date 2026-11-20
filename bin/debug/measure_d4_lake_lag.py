"""Measure RFC 0015 D4's opening-ice lag on a diagnosis save.

Run from the repository root with backend/.venv/bin/python
bin/debug/measure_d4_lake_lag.py <save> --opening-ice --myr 10 --step-years 250000 --out <json>.
Stop at the first step that crosses the preregistered 1% threshold.
"""

import json
import time

from replay_lib import load_replay, parse_replay_args


def _add_arguments(parser):
    parser.add_argument("--opening-ice", action="store_true", help="Measure the one-step lag on opening geometry.")


def main():
    args = parse_replay_args(__doc__, _add_arguments)
    world, result = load_replay(args, 990)

    from app import world as world_module
    from app.hydroclimate import hydrology

    if args.opening_ice:
        hydrology.LAKES_USE_CLOSING_GLACIER_ICE = False
    result["lake_geometry"] = "closing" if hydrology.LAKES_USE_CLOSING_GLACIER_ICE else "opening"
    result["save"] = args.save.name
    del result["backend"]
    world.simulate_glacier_march = True
    result["seed"] = world.seed
    result["ice_age_period_years"] = world.ice_age_period_years
    result["start_myr"] = world.elapsed_years / 1e6
    threshold = hydrology.GLACIER_VISIBLE_DEPTH_M
    for step in range(1, round(args.myr * 1e6 / args.step_years) + 1):
        start = time.perf_counter()
        world_module.step_world(world, args.step_years)
        hydro = world.hydrology_cache
        glacier = hydro.glacier_depth >= threshold
        opening = hydro.hypsometry.glacier_mask if not hydrology.LAKES_USE_CLOSING_GLACIER_ICE else None
        covered = glacier if opening is None else glacier & ~opening
        volume = hydro.lake_depth * hydro.hypsometry.area_m2
        lake_m3 = float(volume.sum())
        covered_m3 = float(volume[covered].sum())
        share = covered_m3 / lake_m3 if lake_m3 else 0.0
        record = {
            "step": step,
            "lake_liquid_m3": lake_m3,
            "covered_liquid_m3": covered_m3,
            "share": share,
            "lake_cap_transfer_m3_water": hydro.lake_water_ledger.glacier_transfer_m3,
            "lake_ledger_residual_m3_water": hydro.lake_water_ledger.residual_m3,
            "land_ledger_residual_m3_water": hydro.land_water_ledger.residual_m3,
            "ice_ledger_residual_m3_ice": hydro.ice.ledger.residual_m3,
            "step_seconds": time.perf_counter() - start,
        }
        result["checkpoints"].append(record)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(record), flush=True)
        if share > 0.01:
            break


if __name__ == "__main__":
    main()
