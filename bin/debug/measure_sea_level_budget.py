#!/usr/bin/env python3
"""Issue #571: where does sea level go? Replay a save and split each step's sea-level change
into water held on land and ocean-basin capacity.

Eustasy holds the total surface-water budget (`World.ocean_water_volume_m3`) fixed and fills
the connected ocean with what isn't held on land (`eustasy.trapped_water_volume_m3`: ice as
water, and lakes at least `_LAKE_MIN_TRAPPED_DEPTH_M` deep). So the change in sea level between
two checkpoints has two sources, both reported in metres of sea level over the opening ocean
area:
- **trapped water**: the change in ice debit and in trapped lake water;
- **basin capacity**: the change in how much water the connected ocean holds below the opening
  sea level (`eustasy.ocean_water_volume_m3` at that fixed level). Tectonics, sedimentation
  and isostasy reshape the hypsometry.

Each checkpoint also records the budget itself (it should stay constant), all lake water at
any depth, the number and area of trapped lakes, and the ocean area. To split the basin term,
it reports the area below the opening sea level and the mean depth below it. For oceanic and
continental columns (by node crust density), it reports their area, mean elevation, and mean
crustal (Hc) and mantle-lithosphere (Hm) thickness, which set isostatic elevation.

`--backend` picks the tree, so the same script measures `origin/main` and an older commit:

    cd backend
    .venv/bin/python ../bin/debug/measure_sea_level_budget.py <save> --myr 2 --step-years 250000 \\
        --checkpoint-myr 0.25 --out <scratch>/main.json
    .venv/bin/python ../bin/debug/measure_sea_level_budget.py <save> --myr 2 --step-years 250000 \\
        --checkpoint-myr 0.25 --backend <pre-change-tree>/backend --out <scratch>/before.json
"""

from __future__ import annotations

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=571)
    from app import lithosphere  # noqa: E402
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import eustasy  # noqa: E402

    opening_level = float(world.sea_level_m)
    elevations, hyps_areas = eustasy.hypsometry(world)
    opening_capacity = eustasy.ocean_water_volume_m3(world, elevations, hyps_areas, opening_level)
    opening = {}

    def snapshot() -> dict:
        areas = plates_mod.collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
        ice = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
        lake = np.clip(plates_mod.collect_all_lake_depth(world.plates), 0.0, None)
        trapped_lake = lake >= eustasy._LAKE_MIN_TRAPPED_DEPTH_M
        elev, h_areas = eustasy.hypsometry(world)
        capacity = eustasy.ocean_water_volume_m3(world, elev, h_areas, opening_level)
        ocean_area = float(h_areas[elev < world.sea_level_m].sum())
        below = elev < opening_level
        hc = plates_mod.collect_all_crustal_thickness(world.plates)
        hm = plates_mod.collect_all_mantle_lithosphere_thickness(world.plates)
        rho = np.concatenate([lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in world.plates])
        oceanic = rho == lithosphere.RHO_OCEANIC_CRUST
        columns = {}
        for name, mask in (("oceanic", oceanic), ("continental", ~oceanic)):
            w = areas[mask]
            columns[name] = {
                "area_mkm2": round(float(w.sum()) / 1e12, 3),
                "mean_elevation_m": round(float(np.average(elev[mask], weights=w)), 1) if w.sum() > 0 else None,
                "mean_hc_m": round(float(np.average(hc[mask], weights=w)), 1) if w.sum() > 0 else None,
                "mean_hm_m": round(float(np.average(hm[mask], weights=w)), 1) if w.sum() > 0 else None,
            }
        out = {
            "elapsed_myr": round(world.elapsed_years / 1e6, 4),
            "sea_level_m": float(world.sea_level_m),
            "budget_km3": float(world.ocean_water_volume_m3) / 1e9,
            "ice_debit_km3": eustasy._ICE_WATER_EQUIVALENT * float(np.dot(areas, ice)) / 1e9,
            "trapped_lake_km3": float(np.dot(areas, np.where(trapped_lake, lake, 0.0))) / 1e9,
            "all_lake_km3": float(np.dot(areas, lake)) / 1e9,
            "trapped_lake_cells": int(np.count_nonzero(trapped_lake)),
            "trapped_lake_area_km2": float(areas[trapped_lake].sum()) / 1e6,
            "capacity_at_opening_level_km3": capacity / 1e9,
            "ocean_area_km2": ocean_area / 1e6,
            "area_below_opening_level_km2": float(h_areas[below].sum()) / 1e6,
            "mean_depth_below_opening_level_m": float(np.average(opening_level - elev[below], weights=h_areas[below])),
            "columns": columns,
        }
        if not opening:
            opening.update(out, area_m2=ocean_area)
        per_m = opening["area_m2"] / 1e9  # km^3 per metre of sea level over the opening ocean
        out["from_trapped_m"] = -((out["ice_debit_km3"] - opening["ice_debit_km3"]) + (out["trapped_lake_km3"] - opening["trapped_lake_km3"])) / per_m
        out["from_capacity_m"] = -(out["capacity_at_opening_level_km3"] - opening_capacity / 1e9) / per_m
        out["sea_level_change_m"] = out["sea_level_m"] - opening_level
        return out

    run_replay(world, args, result, snapshot, lambda seconds: None, lambda steps: {})


if __name__ == "__main__":
    main()
