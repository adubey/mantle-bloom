#!/usr/bin/env python3
"""Issue #353: replay a saved world and measure its glacier ice, on `main` or on the branch.

Replays `--myr` from the save in `--step-years` steps and records, at every checkpoint:

- land fraction (by area) and the ice cap: the fraction of land under visible ice
  (>= GLACIER_VISIBLE_DEPTH_M), ice volume (km^3), and mean and maximum ice depth over
  visibly glaciated cells;
- glacial erosion, summed over the steps since the last checkpoint (km^3): `raw` is the
  glacier term's detachment capacity over land before any cap (old law: coefficient * slope *
  ice factor; new law: coefficient * sliding speed), comparable across both; `applied` is
  what was taken (branch only, from the erosion budget);
- lake water frozen into the ice (km^3, summed like erosion): the branch's ledger value, or on
  `main` the prior lake depth on nodes main's rule freezes (below freezing, or under 10 m of
  ice);
- the branch's ice ledger (accumulation, melt, calving, overflow, frozen river ice) and its
  worst residual, and the share of step time spent in `ice_sheet.step_ice`;
- step runtime and peak RSS.

`--backend` points at the backend to import, so the same script measures both trees:

    cd backend
    .venv/bin/python ../bin/debug/measure_ice_sheet.py ~/Downloads/<save>.mbworld --myr 10 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_ice_sheet.py ~/Downloads/<save>.mbworld --myr 10 \\
        --backend ../../mantle-bloom-origin-main/backend --out <scratch>/main.json
"""

from __future__ import annotations

import time

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=353)
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import erosion, hydrology  # noqa: E402

    try:
        from app.hydroclimate import ice_sheet  # noqa: E402
    except ImportError:
        ice_sheet = None

    # Taps: the glacier depth erosion sees (the old law's input), the erosion result, and the
    # time spent in the ice solve.
    tap = {"prior_glacier": None, "ice_seconds": 0.0, "lake_frozen_km3": 0.0}
    original_hydrology = hydrology.compute_hydrology

    def tapped_hydrology(world_, precipitation, temperature, *a, **k):
        tap["prior_glacier"] = plates_mod.collect_all_glacier_depth(world_.plates)
        lake = plates_mod.collect_all_lake_depth(world_.plates)
        frozen = (temperature < hydrology.FREEZE_POINT_C) | (tap["prior_glacier"] >= hydrology.GLACIER_VISIBLE_DEPTH_M)
        areas = plates_mod.collect_all_accounting_areas_m2(world_.plates, line_spacing_rad(world_.node_density))
        tap["lake_frozen_km3"] += float(np.dot(np.where(frozen, lake, 0.0), areas)) / 1e9
        return original_hydrology(world_, precipitation, temperature, *a, **k)

    hydrology.compute_hydrology = tapped_hydrology
    if ice_sheet is not None:
        original_step_ice = ice_sheet.step_ice

        def timed_step_ice(*a, **k):
            started = time.perf_counter()
            try:
                return original_step_ice(*a, **k)
            finally:
                tap["ice_seconds"] += time.perf_counter() - started

        ice_sheet.step_ice = timed_step_ice

    def glacial_erosion_km3(res, hydro, years) -> tuple[float, float | None]:
        areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
        land = ~hydro.is_ocean
        dt_myr = years / 1e6
        if hasattr(erosion, "GLACIER_ABRASION_COEFFICIENT"):
            sliding = hydro.ice.sliding_speed_m_per_yr if hydro.ice is not None else np.zeros(len(areas))
            raw = erosion.GLACIER_ABRASION_COEFFICIENT * world.glacier_erosion_multiplier * np.power(sliding, erosion.GLACIER_ABRASION_EXPONENT) * dt_myr
        else:
            ice_factor = np.clip(tap["prior_glacier"] / erosion.GLACIER_EROSION_REFERENCE_DEPTH_M, 0.0, erosion.GLACIER_EROSION_MAX_FACTOR)
            raw = erosion.GLACIER_EROSION_COEFFICIENT * world.glacier_erosion_multiplier * res.slope * ice_factor * dt_myr
        applied = (res.budget or {}).get("glacial_erosion_m3")
        return float(np.sum(raw[land] * areas[land])) / 1e9, None if applied is None else applied / 1e9

    def snapshot() -> dict:
        hydro = world.hydrology_cache
        if hydro is None:
            return {"elapsed_myr": world.elapsed_years / 1e6}
        areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
        land = ~hydro.is_ocean
        ice = hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
        on_land = land & ice
        land_area = max(float(areas[land].sum()), 1.0)
        out = {
            "elapsed_myr": world.elapsed_years / 1e6,
            "land_pct": 100.0 * float(areas[land].sum() / areas.sum()),
            "ice_cap_pct_of_land": 100.0 * float(areas[on_land].sum()) / land_area,
            "ice_volume_km3": float(np.dot(areas, hydro.glacier_depth)) / 1e9,
            "ice_mean_depth_m": float(np.dot(areas[ice], hydro.glacier_depth[ice]) / max(areas[ice].sum(), 1.0)),
            "ice_max_depth_m": float(hydro.glacier_depth.max()) if len(hydro.glacier_depth) else 0.0,
            "sea_level_m": float(world.sea_level_m),
        }
        if getattr(hydro, "ice", None) is not None:
            out["ice_bodies"] = len(hydro.ice.bodies)
            out["ice_thickness_max_m"] = float(hydro.ice.thickness.max())
            out["warm_based_pct_of_ice"] = 100.0 * float(areas[ice & hydro.ice.warm_based].sum()) / max(float(areas[ice].sum()), 1.0)
        return out

    def empty_window() -> dict:
        return {"raw": 0.0, "applied": 0.0, "seconds": 0.0, "ice_seconds": 0.0, "ledger": {}, "worst_residual_m3": 0.0}

    window = empty_window()

    def after_step(seconds: float) -> None:
        window["seconds"] += seconds
        window["ice_seconds"] += tap["ice_seconds"]
        tap["ice_seconds"] = 0.0
        res, hydro = world.erosion_cache, world.hydrology_cache
        if res is not None and hydro is not None:
            raw, applied = glacial_erosion_km3(res, hydro, args.step_years)
            window["raw"] += raw
            window["applied"] = None if applied is None or window["applied"] is None else window["applied"] + applied
            ice = getattr(hydro, "ice", None)
            if ice is not None:
                for key, value in vars(ice.ledger).items():
                    window["ledger"][key] = window["ledger"].get(key, 0.0) + value / 1e9
                window["worst_residual_m3"] = max(window["worst_residual_m3"], abs(ice.ledger.residual_m3))

    def close_window(steps: int) -> dict:
        nonlocal window
        closed = {
            "glacial_erosion_raw_km3": window["raw"],
            "glacial_erosion_applied_km3": window["applied"],
            "lake_frozen_km3": window["ledger"].get("frozen_lake_m3", tap["lake_frozen_km3"]),
            "step_seconds_mean": window["seconds"] / steps,
            "ice_step_share_pct": 100.0 * window["ice_seconds"] / max(window["seconds"], 1e-9),
            "ice_ledger_km3": {k: round(v, 3) for k, v in window["ledger"].items()} or None,
            "ice_worst_residual_m3": window["worst_residual_m3"],
        }
        window = empty_window()
        tap["lake_frozen_km3"] = 0.0
        return closed

    run_replay(world, args, result, snapshot, after_step, close_window)

if __name__ == "__main__":
    main()
