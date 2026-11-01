#!/usr/bin/env python3
"""RFC 0016 G0: replay a checkpoint and record the river inventory and rain-only land ledger.

Steps the saved world on the selected backend and records, per step, the river and lake
sub-ledgers, the rain-only land ledger, the registered target lake's supply by source, and how
far river inventory moved between hydrology steps. This does not claim the full G0 gate: RFC 0015
snow and glacier accounts are outside the rain-only boundary.

`--velocity` and `--residence-cap-years` select a labelled sensitivity arm. `--velocity jobson`
swaps Raymond's mean velocity for Jobson's (1996) peak tracer velocity, with mean annual discharge
set to the interval discharge because the model has no such climatology. `--reach-dir` writes each
step's per-reach arrays there as compressed `.npz`; keep those outside the repository.

Example:

    backend/.venv/bin/python bin/debug/replay_rfc0016_river_water.py <checkpoint.mbworld> \\
        --anchors analysis/rfc0016/target-anchor-window-a.json --myr 0.4 --step-years 100000 \\
        --checkpoint-myr 0.1 --out <scratch>/rfc0016-g0-a.json
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from replay_lib import load_replay, parse_replay_args, run_replay

GRAVITY_M_S2 = 9.80665


def add_arguments(ap) -> None:
    ap.add_argument("--anchors", type=Path, required=True, help="JSON file containing anchor_points")
    ap.add_argument("--velocity", choices=("raymond", "jobson"), default="raymond")
    ap.add_argument("--residence-cap-years", type=float, default=1.0)
    ap.add_argument("--reach-dir", type=Path, help="write per-reach arrays for each step here")


def jobson_residence_time_s(length_m, discharge_m3_s, drainage_m2, slope, cap_s):
    """Jobson's peak tracer velocity with `Qa = Q`, as travel time along each reach."""
    q = np.maximum(discharge_m3_s, 1e-12)
    drainage = np.maximum(drainage_m2, 1.0)
    dimensionless_drainage = drainage**1.25 * np.sqrt(GRAVITY_M_S2) / q
    velocity = 0.094 + 0.0143 * dimensionless_drainage**0.919 * np.maximum(slope, 0.0) ** 0.159 * (q / drainage)
    return np.clip(np.maximum(length_m, 1.0) / velocity, 1.0, cap_s)


def use_jobson_velocity(world, hydrology, plates_mod, line_spacing_rad) -> None:
    """Wrap reservoir routing so travel times use each pass's own drainage graph."""
    raymond_route = hydrology.route_river_reservoirs

    def route(points, elevation, is_ocean, flow_target, *rest, **kw):
        _, plates_in_order = plates_mod.gather_node_positions(world.plates)
        areas = plates_mod.collect_all_accounting_areas_m2(
            plates_in_order, line_spacing_rad(world.node_density)
        )
        drainage, _ = hydrology.route_downstream(elevation, is_ocean, flow_target, areas)
        slope = hydrology._slope_to_flow_target(points, elevation, flow_target)
        hydrology.river_residence_time_s = lambda length_m, q: jobson_residence_time_s(
            length_m, np.asarray(q, dtype=float), drainage, slope, hydrology.RIVER_MAX_RESIDENCE_SECONDS
        )
        return raymond_route(points, elevation, is_ocean, flow_target, *rest, **kw)

    hydrology.route_river_reservoirs = route


def target_lake(hydro, areas, anchor_points, visible_depth, lake_components) -> dict:
    """The wet component under the registered anchor, with this step's supply by source. The
    anchor is a location, not a lake identity, so a dry anchor reports no component."""
    _, nearest = cKDTree(hydro.points).query(anchor_points, k=1)
    wet = np.asarray(hydro.lake_depth) > visible_depth
    anchor = int(np.atleast_1d(nearest)[0])
    if not wet[anchor]:
        return {"anchor_wet": False}
    members = next(
        np.asarray(nodes) for nodes in lake_components(wet, hydro.neighbor_idx) if anchor in set(nodes.tolist())
    )
    sources = {
        "wet_precip": hydro.lake_external_wet_precip_m3_per_yr,
        "dry_runoff": hydro.lake_external_dry_runoff_m3_per_yr,
        "meltwater": hydro.lake_external_meltwater_m3_per_yr,
        "river_storage": hydro.lake_external_river_storage_m3_per_yr,
    }
    return {
        "anchor_wet": True,
        "volume_km3": float(np.dot(hydro.lake_depth[members], areas[members])) / 1e9,
        "area_km2": float(areas[members].sum()) / 1e6,
        "supply_m3_per_yr": {name: float(rate[members].sum()) for name, rate in sources.items()},
    }


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=725)
    # The registered #680 windows replay with the ice-age cycle off (measure_lake_step_budget.py).
    world.ice_age_period_years = 0.0
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import hydrology  # noqa: E402

    hydrology.RIVER_MAX_RESIDENCE_SECONDS = args.residence_cap_years * 365.25 * 86400.0
    if args.velocity == "jobson":
        use_jobson_velocity(world, hydrology, plates_mod, line_spacing_rad)
    anchors = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
    result.update(velocity=args.velocity, residence_cap_years=args.residence_cap_years, steps=[])
    if args.reach_dir is not None:
        args.reach_dir.mkdir(parents=True, exist_ok=True)
    previous_final_m3 = [None]

    def snapshot() -> dict:
        return {"elapsed_myr": world.elapsed_years / 1e6}

    def write_reaches(hydro, areas, step: int) -> None:
        np.savez_compressed(
            args.reach_dir / f"reaches-step{step:03d}.npz",
            points=hydro.points, flow_target=hydro.flow_target, is_ocean=hydro.is_ocean,
            lake_depth=hydro.lake_depth, inflow_m3=hydro.river_water_inflow_m3,
            initial_m3=hydro.river_water_initial_m3, final_m3=hydro.river_water_m * areas,
            outflow_m3=hydro.river_water_outflow_m3, routed_m3=hydro.river_water_routed_m3,
            evaporation_m3=hydro.river_water_evaporation_m3, residual_m3=hydro.river_water_residual_m3,
        )

    def after_step(seconds: float) -> None:
        hydro = world.hydrology_cache
        areas = plates_mod.collect_all_accounting_areas_m2(
            hydro.plates_in_order, line_spacing_rad(world.node_density)
        )
        river, land, lake = hydro.river_water_ledger, hydro.land_water_ledger, hydro.lake_water_ledger
        cycles = hydrology.flow_cycles(hydro.flow_target)
        spill_edge_cycles = sum(
            np.any(hydro.flow_target[cycle] == hydro.spill_target[cycle])
            for cycle in cycles
        )
        step = len(result["steps"]) + 1
        drift = None if previous_final_m3[0] is None else river.initial_storage_m3 - previous_final_m3[0]
        previous_final_m3[0] = river.final_storage_m3
        external = land.rain_m3 + land.meltwater_m3
        result["steps"].append({
            "step": step,
            "elapsed_myr": world.elapsed_years / 1e6,
            "runtime_s": seconds,
            "river": asdict(river),
            "land_rain_only": {**asdict(land), "residual_m3": land.residual_m3},
            "lake": {
                key: getattr(lake, key)
                for key in ("external_supply_m3", "evaporation_m3", "overflow_out_m3", "residual_m3")
            },
            "normalized_land_residual": abs(land.residual_m3) / external,
            "normalized_stranded_overflow": land.stranded_overflow_m3 / external,
            "normalized_land_sink_settlement": getattr(land, "land_sink_settlement_m3", 0.0) / external,
            "normalized_lake_residual": abs(lake.residual_m3) / max(lake.external_supply_m3, 1.0),
            "normalized_river_residual": abs(river.residual_m3) / max(river.input_m3, 1.0),
            "ocean_export_km3": land.ocean_export_m3 / 1e9,
            "internal_transfer_km3": lake.internal_transfer_m3 / 1e9,
            "flow_cycles": len(cycles),
            "spill_edge_cycles": int(spill_edge_cycles),
            "river_routing_passes": 2,
            "inter_step_river_drift_m3": drift,
            "river_storage_after_step_km3": float(np.dot(hydro.river_water_m, areas)) / 1e9,
            "reaches_holding_water": int(np.count_nonzero(hydro.river_water_m > 0.0)),
            "lake_volume_km3": float(np.dot(hydro.lake_depth, areas)) / 1e9,
            "target": target_lake(
                hydro, areas, anchors, hydrology.LAKE_MIN_VISIBLE_DEPTH_M, hydrology.lake_components
            ),
        })
        if args.reach_dir is not None:
            write_reaches(hydro, areas, step)

    def close_window(_steps: int) -> dict:
        return {"completed_steps": len(result["steps"])}

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
