"""Measure RFC 0015 D4's opening-ice lag on a diagnosis save.

Run from the repository root with backend/.venv/bin/python
bin/debug/measure_d4_lake_lag.py <save> --myr 10 --step-years 250000 --out <json>.
Stop at the first step that crosses the preregistered 1% threshold.
"""

import json
import time

from replay_lib import load_replay, parse_replay_args


def main():
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, 990)

    from app import world as world_module
    from app.hydroclimate import hydrology

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
        opening = hydro.hypsometry.glacier_mask
        closing = hydro.glacier_depth >= threshold
        newly_visible = closing if opening is None else closing & ~opening
        volume = hydro.lake_depth * hydro.hypsometry.area_m2
        lake_m3 = float(volume.sum())
        covered_m3 = float(volume[newly_visible].sum())
        share = covered_m3 / lake_m3 if lake_m3 else 0.0
        record = {
            "step": step,
            "lake_liquid_m3": lake_m3,
            "newly_covered_liquid_m3": covered_m3,
            "share": share,
            "step_seconds": time.perf_counter() - start,
        }
        result["checkpoints"].append(record)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(record), flush=True)
        if share > 0.01:
            break


if __name__ == "__main__":
    main()
