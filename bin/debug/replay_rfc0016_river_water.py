#!/usr/bin/env python3
"""RFC 0016 G0: replay a checkpoint and summarize river/lake water inventories.

Runs the saved world on the selected backend and records river and lake storage, transfers,
subsystem residuals, and peak RSS at each checkpoint. This is a diagnostic of the river and
lake sub-ledgers; it does not claim the strict whole-land G0 gate while RFC 0015 snow/glacier
accounting is absent.

Example:

    backend/.venv/bin/python bin/debug/replay_rfc0016_river_water.py <checkpoint.mbworld> --myr 0.4 --step-years 100000 --checkpoint-myr 0.4 --out <scratch>/rfc0016-g0-a.json
"""

from __future__ import annotations

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=725)
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402

    def empty_window() -> dict[str, float]:
        return {
            "step_count": 0.0,
            "river_input_m3": 0.0,
            "river_initial_storage_m3": 0.0,
            "river_final_storage_m3": 0.0,
            "river_evaporation_m3": 0.0,
            "river_lake_transfer_m3": 0.0,
            "river_ocean_export_m3": 0.0,
            "river_sink_settlement_m3": 0.0,
            "river_residual_m3": 0.0,
            "river_max_reach_residual_m3": 0.0,
            "lake_external_supply_m3": 0.0,
            "lake_residual_m3": 0.0,
        }

    window = empty_window()

    def snapshot() -> dict:
        hydro = world.hydrology_cache
        if hydro is None:
            return {"elapsed_myr": world.elapsed_years / 1e6}
        areas = plates_mod.collect_all_accounting_areas_m2(
            hydro.plates_in_order, line_spacing_rad(world.node_density)
        )
        land = ~hydro.is_ocean
        river_volume = float(np.dot(np.clip(hydro.river_water_m, 0.0, None), areas))
        lake_volume = float(np.dot(np.clip(hydro.lake_depth, 0.0, None), areas))
        lake_ice = float(np.dot(np.clip(hydro.lake_ice_m, 0.0, None), areas))
        return {
            "elapsed_myr": world.elapsed_years / 1e6,
            "river_storage_m3": river_volume,
            "river_storage_km3": river_volume / 1e9,
            "lake_liquid_m3": lake_volume,
            "lake_ice_water_equivalent_m3": lake_ice,
            "land_area_m2": float(areas[land].sum()),
            "river_ledger_available": hydro.river_water_ledger is not None,
        }

    def after_step(_seconds: float) -> None:
        hydro = world.hydrology_cache
        if hydro is None:
            return
        river = hydro.river_water_ledger
        if river is not None:
            if window["step_count"] == 0.0:
                window["river_initial_storage_m3"] = river.initial_storage_m3
            window["river_final_storage_m3"] = river.final_storage_m3
            window["river_input_m3"] += river.input_m3
            window["river_evaporation_m3"] += river.evaporation_m3
            window["river_lake_transfer_m3"] += river.lake_transfer_m3
            window["river_ocean_export_m3"] += river.ocean_export_m3
            window["river_sink_settlement_m3"] += river.sink_settlement_m3
            window["river_residual_m3"] += river.residual_m3
            window["river_max_reach_residual_m3"] = max(
                window["river_max_reach_residual_m3"], river.max_reach_residual_m3
            )
        window["step_count"] += 1.0
        lake = hydro.lake_water_ledger
        if lake is not None:
            window["lake_external_supply_m3"] += lake.external_supply_m3
            window["lake_residual_m3"] += lake.residual_m3

    def close_window(steps: int) -> dict:
        nonlocal window
        result_window = {key: float(value) for key, value in window.items()}
        result_window.pop("step_count")
        result_window["steps"] = steps
        window = empty_window()
        return {"river_lake_subledgers": result_window}

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
