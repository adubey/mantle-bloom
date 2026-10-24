#!/usr/bin/env python3
"""RFC 0015 gate G0 (#554): replay a save under one arm of the physical-rate ice prototype.

Gate branch only. `ice_sheet.G0` holds the prototype's switches; every switch off is
origin/main. The arm presets:

- `main`: origin/main, model-unit mass balance (the baseline arm).
- `a`: D1's physical mass balance with today's inputs and root rule.
- `b`: (a) with sea level held at the run's opening value.
- `c`: (a) with the mass-balance surface on the run's opening deflection.
- `e`: (a) with D1's seasonal snow fraction in place of the is_frozen switch.
- `h1`, `h2`: (a) keeping the largest or the smallest stable root.
- `f`: (a) with membership and yield stress re-evaluated at the solved thickness.
- `b2`, `c2`: (a) with sea level (b') or the mass-balance deflection (c') solved within the step.
- `j`: (a) with frictional heat from the current step by outer iteration (combine with --set).
- `g`: every input fix, (b') + (c') + (e) + (f); `g-b2`, `g-c2`, `g-e`, `g-f` leave one out.
  Arm (d), frozen river and lake ice off, is origin/main since #550, so it has no switch.
- `i1`, `i2`: (g) with the largest or the smallest stable root.

`--set key=value` overrides any `G0Config` field, for combined arms. Physical arms use D3's
starting abrasion law (E = 1e-4 u_b^0.65 m/yr; coefficient 100 m/Myr) and drop the ice-age
sea-ice cap (D6 removes it), and meltwater is routed by mass. Every physical arm enumerates
all roots of each body's residual and logs, per step, how many bodies have more than one
stable root and which the root rule kept.

Records a checkpoint every `--checkpoint-myr` (default: every step): ice volume and visible
area, sea level, the ice ledger, the root log, ice-age intensity, and step and ice-solve time.

    cd backend
    .venv/bin/python ../bin/debug/replay_g0.py <save> --arm a --myr 10 --step-years 250000 \\
        --checkpoint-myr 0.25 --out <scratch>/g0-a.json
    .venv/bin/python ../bin/debug/replay_g0.py <save> --arm a --ice-age-period-years 8e6 --myr 18 ...
"""

from __future__ import annotations

import argparse
import time

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay

ARMS = {
    "main": {},
    "a": {"physical": True},
    "b": {"physical": True, "hold_sea_level": True},
    "c": {"physical": True, "hold_deflection": True},
    "e": {"physical": True, "f_snow": True},
    "h1": {"physical": True, "root_rule": "largest"},
    "h2": {"physical": True, "root_rule": "smallest"},
    "f": {"physical": True, "reevaluate_membership": True},
    "b2": {"physical": True, "couple_sea_level": True},
    "c2": {"physical": True, "couple_deflection": True},
    "j": {"physical": True, "current_friction": True},
}
_FIXES = {"b2": "couple_sea_level", "c2": "couple_deflection", "e": "f_snow", "f": "reevaluate_membership"}
ARMS["g"] = {"physical": True, **{field: True for field in _FIXES.values()}}
for _fix, _field in _FIXES.items():
    ARMS[f"g-{_fix}"] = {**ARMS["g"], _field: False}
ARMS["i1"] = {**ARMS["g"], "root_rule": "largest"}
ARMS["i2"] = {**ARMS["g"], "root_rule": "smallest"}


def add_arguments(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--arm", choices=sorted(ARMS), required=True)
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a G0Config field")
    ap.add_argument("--ice-age-period-years", type=float, default=None)
    ap.add_argument("--no-enumerate", action="store_true", help="skip root enumeration on physical arms")


def _parse_value(text: str):
    if text in ("true", "false"):
        return text == "true"
    try:
        return float(text)
    except ValueError:
        return text


def configure(args: argparse.Namespace, world) -> dict:
    """Set ice_sheet.G0 and the arm's patches; returns the configuration used."""
    from app.hydroclimate import erosion, eustasy, ice_sheet  # noqa: E402

    settings = dict(ARMS[args.arm])
    for item in args.set:
        key, value = item.split("=", 1)
        settings[key] = _parse_value(value)
    for key, value in settings.items():
        if not hasattr(ice_sheet.G0, key):
            raise SystemExit(f"unknown G0 field {key}")
        setattr(ice_sheet.G0, key, value)
    g0 = ice_sheet.G0
    if g0.physical:
        g0.enumerate_roots = not args.no_enumerate
        erosion.GLACIER_ABRASION_COEFFICIENT = 100.0
        erosion.GLACIER_ABRASION_EXPONENT = 0.65
    if g0.hold_sea_level:
        eustasy.update_sea_level = lambda world_: None
    if g0.hold_deflection:
        for plate in world.plates:
            plate.set_fields_on_plate(g0_hold_deflection_m=plate.collect("ice_load_deflection_m").copy())
    if args.ice_age_period_years is not None:
        world.ice_age_period_years = args.ice_age_period_years
    return {
        "arm": args.arm,
        **{k: v for k, v in vars(g0).items() if k != "frictional_heat_w_m2"},
        "abrasion_coefficient": erosion.GLACIER_ABRASION_COEFFICIENT,
        "abrasion_exponent": erosion.GLACIER_ABRASION_EXPONENT,
        "ice_age_period_years": float(getattr(world, "ice_age_period_years", 0.0) or 0.0),
    }


def summarise_roots(log: list[dict]) -> dict:
    multi = [r for r in log if r["stable"] > 1]
    return {
        "bodies": len(log),
        "bodies_multi_stable": len(multi),
        "max_stable": max((r["stable"] for r in log), default=0),
        "max_roots": max((r["roots"] for r in log), default=0),
        "walk_not_largest": sum(1 for r in multi if r["walk_rank"] != r["stable"] - 1),
        "walk_on_unstable": sum(1 for r in log if r["walk_rank"] < 0 and not r["saturated"] and r["roots"] > 0),
        "no_stable": sum(1 for r in log if r["stable"] == 0 and not r["saturated"]),
        "saturated": sum(1 for r in log if r["saturated"]),
        "kept_by_rule": sum(1 for r in log if r["kept"] != "walk"),
        "largest_multi": sorted(multi, key=lambda r: -r["opening_km3"])[:3],
    }


def summarise_iterations(log: list[dict]) -> dict | None:
    if not log:
        return None
    return {"solves": len(log), "max_iterations": max(r["iterations"] for r in log), "not_converged": sum(not r["converged"] for r in log)}


def summarise_outer(log: list[dict]) -> dict | None:
    inner = [r for r in log if r["kind"] == "inner"]
    sea = [r for r in log if r["kind"] == "sea_level"]
    if not log:
        return None
    out = {"inner": summarise_iterations(inner)}
    if sea:
        out["sea_level"] = {
            "steps": len(sea),
            "max_evaluations": max(r["evaluations"] for r in sea),
            "unbracketed": sum(not r["bracketed"] for r in sea),
            "max_residual_m": max(r["residual_m"] for r in sea),
            "last": sea[-1],
        }
    return out


def main() -> None:
    args = parse_replay_args(__doc__, add_arguments)
    world, result = load_replay(args, issue=554)
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import climate, hydrology, ice_sheet  # noqa: E402

    result["config"] = configure(args, world)
    tap = {"ice_seconds": 0.0}
    original_step_ice = ice_sheet.step_ice

    def timed_step_ice(*a, **k):
        started = time.perf_counter()
        try:
            return original_step_ice(*a, **k)
        finally:
            tap["ice_seconds"] += time.perf_counter() - started

    ice_sheet.step_ice = timed_step_ice
    window = {"seconds": 0.0, "ice_seconds": 0.0, "ledger": {}, "worst_residual_m3": 0.0, "roots": [], "membership": [], "outer": []}

    def snapshot() -> dict:
        areas = plates_mod.collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
        depth = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
        visible = depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
        return {
            "elapsed_myr": round(world.elapsed_years / 1e6, 4),
            "ice_volume_km3": float(np.dot(areas, depth)) / 1e9,
            "visible_ice_area_km2": float(areas[visible].sum()) / 1e6,
            "ice_max_depth_m": float(depth.max()) if len(depth) else 0.0,
            "sea_level_m": float(world.sea_level_m),
            "ice_age_intensity": round(climate.ice_age_glacial_intensity(world), 4),
        }

    def after_step(seconds: float) -> None:
        window["seconds"] += seconds
        window["ice_seconds"] += tap["ice_seconds"]
        tap["ice_seconds"] = 0.0
        window["roots"].extend(ice_sheet.G0_ROOT_LOG)
        ice_sheet.G0_ROOT_LOG.clear()
        window["membership"].extend(ice_sheet.G0_MEMBERSHIP_LOG)
        ice_sheet.G0_MEMBERSHIP_LOG.clear()
        window["outer"].extend(ice_sheet.G0_OUTER_LOG)
        ice_sheet.G0_OUTER_LOG.clear()
        hydro = world.hydrology_cache
        ice = getattr(hydro, "ice", None) if hydro is not None else None
        if ice is not None:
            for key, value in vars(ice.ledger).items():
                if key != "capped_nodes":
                    window["ledger"][key] = window["ledger"].get(key, 0.0) + value / 1e9
            window["worst_residual_m3"] = max(window["worst_residual_m3"], abs(ice.ledger.residual_m3))
            window["warm_based_nodes"] = int(np.count_nonzero(ice.warm_based))
            window["sliding_max_m_per_yr"] = float(ice.sliding_speed_m_per_yr.max()) if len(ice.sliding_speed_m_per_yr) else 0.0

    def close_window(steps: int) -> dict:
        closed = {
            "step_seconds_mean": window["seconds"] / steps,
            "ice_step_share_pct": 100.0 * window["ice_seconds"] / max(window["seconds"], 1e-9),
            "ice_ledger_km3": {k: round(v, 3) for k, v in window["ledger"].items()},
            "ice_worst_residual_m3": window["worst_residual_m3"],
            "warm_based_nodes": window.get("warm_based_nodes"),
            "sliding_max_m_per_yr": window.get("sliding_max_m_per_yr"),
            "roots": summarise_roots(window["roots"]) if window["roots"] else None,
            "membership": summarise_iterations(window["membership"]),
            "outer": summarise_outer(window["outer"]),
        }
        window.update(seconds=0.0, ice_seconds=0.0, ledger={}, worst_residual_m3=0.0, roots=[], membership=[], outer=[])
        return closed

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
