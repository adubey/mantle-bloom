#!/usr/bin/env python3
"""Issue #811: replay a save at one step size and record cross-process metrics per step.

Steps a saved world with the ice-age cycle off (as #754's windows do) at `--step-years`, in
`diagnostic` or `cfd` wind (`--wind`), and records per step:

- whole-world liquid lake volume and lake ice, sea level, and total glacier ice volume
  (`glacier_depth` x accounting area; the ice metrics may change when RFC 0015 D2 lands);
- the largest land discharge (m^3/s) and the mean of the ten largest, as a proxy for the
  biggest river mouths;
- the lake and land water-ledger residuals, normalized as #757's replay does;
- area-weighted mean precipitation and air temperature from the climate cache;
- the CFD atmosphere's simulated seconds (`cfd` only), against the world's elapsed seconds,
  to show the clock mismatch directly;
- step runtime.

The summary gives, per series, the lag-1 autocorrelation of the step change (about -1 is
period-2 alternation), the mean absolute step change from step 3, the largest ratio between
consecutive steps, and the relative drift first-to-last.

It also counts #754's large-body flicker reversals per step (`body_reversals`, over
`body_comparable` matched bodies above 10,000 km^2; see `lake_body_tracking.py`), and the summary
gives their total and the rate per comparable body-step.

Nothing here changes the model. Compare modes with a fixed step count (the issue uses 200).

Example (day mode, 200 steps):

    backend/.venv/bin/python bin/debug/replay_timescale_sweep.py <save>.mbworld \\
        --step-years 0.0027379 --myr 5.4757e-7 --checkpoint-myr 1.4e-7 --wind diagnostic \\
        --out <scratch>/issue811-day-diagnostic.json
"""

from __future__ import annotations

import numpy as np
from lake_body_tracking import components, count_reversals
from replay_lib import load_replay, parse_replay_args, run_replay

M3_PER_KM3 = 1e9
SECONDS_PER_YEAR = 3.15576e7
TOP_RIVERS = 10
DAY_YEARS = 1.0 / 365.25


def add_arguments(ap) -> None:
    ap.add_argument("--wind", choices=("diagnostic", "cfd"), default="diagnostic")
    ap.add_argument(
        "--warmup-days", type=int, default=10,
        help="Day-length steps run before measuring, the same in every mode: a save written before #758 "
        "sheds a lake-sized volume and reclassifies about 2 million km^2 as ocean at its second step (#811)",
    )


def series_summary(values: list[float]) -> dict:
    v = np.asarray(values, dtype=float)
    change = np.diff(v)
    ratios = [max(a, b) / max(min(a, b), 1e-9) for a, b in zip(v[:-1], v[1:]) if max(a, b) > 0.0]
    return {
        "first": float(v[0]) if len(v) else None,
        "last": float(v[-1]) if len(v) else None,
        "relative_drift": float((v[-1] - v[0]) / v[0]) if len(v) and v[0] != 0.0 else None,
        "lag1_autocorrelation_of_change": (
            float(np.corrcoef(change[:-1], change[1:])[0, 1]) if len(change) > 2 and change[:-1].std() > 0 and change[1:].std() > 0 else None
        ),
        "mean_abs_change_from_step3": float(np.abs(change[2:]).mean()) if len(change) > 2 else None,
        "largest_consecutive_ratio": float(max(ratios)) if ratios else None,
    }


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=811)
    world.ice_age_period_years = 0.0
    world.wind_model = args.wind
    result.update(wind_model=args.wind, steps=[])
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402

    def snapshot() -> dict:
        return {"elapsed_years": world.elapsed_years, "sea_level_m": world.sea_level_m}

    def climate_means() -> dict:
        fields = world.climate_cache
        if fields is None:
            return {}
        weight = np.broadcast_to(np.cos(np.radians(fields.lat_deg)).reshape(-1, 1), fields.precipitation_mm.shape)
        return {
            "mean_precipitation_mm": float(np.average(fields.precipitation_mm, weights=weight)),
            "mean_air_temperature_c": float(np.average(fields.air_temperature_c, weights=weight)),
        }

    tracked = {"bodies": None}

    def after_step(seconds: float) -> None:
        hydro = world.hydrology_cache
        bodies = components(hydro)
        reversals, comparable, events = count_reversals(tracked["bodies"], bodies)
        tracked["bodies"] = bodies
        areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
        depth, lake_ice = np.asarray(hydro.lake_depth), np.asarray(hydro.lake_ice_m)
        ice = np.clip(np.asarray(hydro.glacier_depth), 0.0, None)
        discharge = np.asarray(hydro.discharge_m3_s) if hydro.discharge_m3_s is not None else np.zeros(1)
        land_discharge = np.where(np.asarray(hydro.is_ocean), 0.0, discharge)
        top = np.sort(land_discharge)[-TOP_RIVERS:]
        lake, land = hydro.lake_water_ledger, hydro.land_water_ledger
        state = world.atmosphere_cfd_state
        step = {
            "step": len(result["steps"]) + 1,
            "elapsed_years": world.elapsed_years,
            "runtime_s": seconds,
            "sea_level_m": world.sea_level_m,
            "world_lake_km3": float(np.dot(depth, areas)) / M3_PER_KM3,
            "world_lake_ice_km3": float(np.dot(lake_ice, areas)) / M3_PER_KM3,
            "glacier_ice_km3": float(np.dot(ice, areas)) / M3_PER_KM3,
            "max_land_discharge_m3_s": float(top[-1]),
            "top_river_mean_discharge_m3_s": float(top.mean()),
            "normalized_lake_residual": abs(lake.residual_m3) / max(lake.external_supply_m3, 1.0) if lake else None,
            "normalized_land_residual": (
                abs(land.residual_m3) / (land.rain_m3 + land.meltwater_m3) if land and (land.rain_m3 + land.meltwater_m3) > 0 else None
            ),
            "body_reversals": reversals,
            "body_comparable": comparable,
            "body_reversal_events": events,
            "cfd_simulated_seconds": float(state.elapsed_seconds) if state is not None else None,
            "world_elapsed_seconds": world.elapsed_years * SECONDS_PER_YEAR,
        }
        step.update(climate_means())
        result["steps"].append(step)

    def close_window(_steps: int) -> dict:
        keys = (
            "sea_level_m", "world_lake_km3", "world_lake_ice_km3", "glacier_ice_km3",
            "max_land_discharge_m3_s", "top_river_mean_discharge_m3_s", "normalized_lake_residual",
            "normalized_land_residual", "mean_precipitation_mm", "mean_air_temperature_c",
        )
        result["summary"] = {
            key: series_summary([s[key] for s in result["steps"] if s.get(key) is not None]) for key in keys
        }
        reversals = sum(s["body_reversals"] for s in result["steps"])
        comparable = sum(s["body_comparable"] for s in result["steps"])
        result["summary"]["body_flicker"] = {
            "reversals": reversals, "comparable_body_steps": comparable,
            "reversal_rate": reversals / comparable if comparable else None,
        }
        result["summary"]["mean_runtime_s"] = float(np.mean([s["runtime_s"] for s in result["steps"]]))
        return {"completed_steps": len(result["steps"])}

    from app import world as world_mod  # noqa: E402

    for _ in range(args.warmup_days):
        world_mod.step_world(world, DAY_YEARS)
    result["warmup_days"] = args.warmup_days
    tracked["bodies"] = components(world.hydrology_cache)
    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
