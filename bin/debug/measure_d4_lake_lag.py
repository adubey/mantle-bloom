"""Measure RFC 0015 D4's opening-ice lag on a diagnosis save.

Run from the repository root with backend/.venv/bin/python
bin/debug/measure_d4_lake_lag.py <save> --opening-ice --myr 10 --step-years 250000 --out <json>.
Stop at the first step that crosses the preregistered 1% threshold. `--surface-relief`
recreates the pre-review relief calculation as a controlled counterfactual.
"""

import json
import time
from dataclasses import replace

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args


def _add_arguments(parser):
    parser.add_argument("--opening-ice", action="store_true", help="Measure the one-step lag on opening geometry.")
    parser.add_argument(
        "--surface-relief", action="store_true",
        help="Use glacier-inflated surface relief as a counterfactual.",
    )


def main():
    args = parse_replay_args(__doc__, _add_arguments)
    world, result = load_replay(args, 990)

    from app import world as world_module
    from app import plates
    from app.hydroclimate import hydrology, lakes

    if args.surface_relief:
        lake_terrain = hydrology._lake_terrain

        def with_surface_relief(bed, glacier_depth, areas, neighbor_idx):
            hypsometry = lake_terrain(bed, glacier_depth, areas, neighbor_idx)
            return replace(hypsometry, relief_m=lakes.subcell_relief(hypsometry.elevation, neighbor_idx))

        hydrology._lake_terrain = with_surface_relief

    if args.opening_ice:
        hydrology.LAKES_USE_CLOSING_GLACIER_ICE = False
    result["lake_geometry"] = "closing" if hydrology.LAKES_USE_CLOSING_GLACIER_ICE else "opening"
    result["relief_geometry"] = "glacier-surface" if args.surface_relief else "rock-bed"
    result["save"] = args.save.name
    del result["backend"]
    world.simulate_glacier_march = True
    result["seed"] = world.seed
    result["ice_age_period_years"] = world.ice_age_period_years
    result["start_myr"] = world.elapsed_years / 1e6
    threshold = hydrology.GLACIER_VISIBLE_DEPTH_M
    opening_hydro = getattr(world, "hydrology_cache", None)
    opening_depth = getattr(opening_hydro, "glacier_depth", None)
    opening_points, _ = plates.gather_node_positions(world.plates)
    if opening_depth is None:
        opening_elevation = getattr(opening_hydro, "elevation", None)
        if opening_elevation is None:
            opening_elevation = np.concatenate([plate.collect("elevation") for plate in world.plates])
        opening_glacier = np.zeros_like(opening_elevation, dtype=bool)
    else:
        opening_glacier = np.asarray(opening_depth) >= threshold
    for step in range(1, round(args.myr * 1e6 / args.step_years) + 1):
        start = time.perf_counter()
        world_module.step_world(world, args.step_years)
        hydro = world.hydrology_cache
        closing_points, _ = plates.gather_node_positions(world.plates)
        glacier = hydro.glacier_depth >= threshold
        # Tectonic remeshing changes the node count within a step. Compare coverage
        # by mapping each closing node to its nearest opening-state node.
        opening_at_closing = opening_glacier[cKDTree(opening_points).query(closing_points)[1]]
        newly_covered = glacier & ~opening_at_closing
        volume = hydro.lake_depth * hydro.hypsometry.area_m2
        lake_m3 = float(volume.sum())
        covered_m3 = float(volume[glacier].sum())
        newly_covered_m3 = float(volume[newly_covered].sum())
        margin = ~glacier & np.isin(np.arange(len(glacier)), hydro.neighbor_idx[glacier].ravel())
        margin_m3 = float(volume[margin].sum())
        other_m3 = float(volume[~glacier & ~margin].sum())
        record = {
            "step": step,
            "lake_liquid_m3": lake_m3,
            "covered_liquid_m3": covered_m3,
            "newly_covered_liquid_m3": newly_covered_m3,
            "glacier_margin_liquid_m3": margin_m3,
            "other_liquid_m3": other_m3,
            "share": newly_covered_m3 / lake_m3 if lake_m3 else 0.0,
            "closing_ice_share": covered_m3 / lake_m3 if lake_m3 else 0.0,
            "lake_cap_transfer_m3_water": hydro.lake_water_ledger.glacier_transfer_m3,
            "lake_ledger_residual_m3_water": hydro.lake_water_ledger.residual_m3,
            "land_ledger_residual_m3_water": hydro.land_water_ledger.residual_m3,
            "ice_ledger_residual_m3_ice": hydro.ice.ledger.residual_m3,
            "step_seconds": time.perf_counter() - start,
        }
        result["checkpoints"].append(record)
        args.out.write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(record), flush=True)
        if record["share"] > 0.01:
            break


if __name__ == "__main__":
    main()
