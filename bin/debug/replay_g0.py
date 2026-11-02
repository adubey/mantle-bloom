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
- `n`: (a) with the climate's per-step random fields (humidity, moisture-flux and current
  noise, drawn from (seed, elapsed_years) in climate.py) frozen at the run's opening draw.
  Outside the RFC's list: it tests whether step-to-step climate noise moves fast ice.
- `m`: (a) marched to equilibrium in short substeps within each step, then held (#554's
  design direction): every fast input is recomputed between substeps.
- `m2`: (m) with the lagged inputs damped (mean of the last two substeps), adaptive substeps
  (250 yr to 10 kyr) and no march cap short of the step, at a tighter settling tolerance.
- `m3`: (m2) at (m)'s settling tolerance (0.1%), halving substeps only for alternations over
  1%, and holding the mean rates of the last six substeps. (m2) never settled: a 0.1-0.3%
  two-substep alternation survives the damping and kept halving its substeps.
- `y`, `q`, `x`: diagnostic arms outside the RFC's list, each (a) less one feedback: yield
  stress fixed at the hard-bed value (y), no frictional heat (q), no glacial erosion (x).

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
import json
import time
from pathlib import Path

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
    "n": {"physical": True},
    "y": {"physical": True, "fixed_yield": True},
    "q": {"physical": True, "no_friction": True},
    "x": {"physical": True},
    "m": {"physical": True, "march": True},
    "m2": {"physical": True, "march": True, "march_damped": True, "march_adaptive": True, "march_tolerance": 2.0e-4},
    "m3": {"physical": True, "march": True, "march_damped": True, "march_adaptive": True,
           "march_flip_fraction": 0.01, "march_hold_window": 6},
}
NO_EROSION_ARMS = {"x"}
FREEZE_NOISE_ARMS = {"n"}
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
    ap.add_argument("--sliding-report", action="store_true", help="report where the sliding is fast, per checkpoint")
    ap.add_argument("--noise-salt-offset", type=int, default=0, help="shift the drifting noise salts (#575 draws)")
    ap.add_argument("--freeze-climate-noise", action="store_true", help="arm (n)'s noise freeze on any arm")
    ap.add_argument("--factorial-myr", type=float, nargs="*", default=[], help="#741: steps (by closing time) to re-march under last step's inputs")
    ap.add_argument("--ensemble-old-climate", action="store_true", help="#741: also run each --factorial-myr ensemble under last step's climate")
    ap.add_argument("--ensemble-max-substep-years", type=float, default=None, help="#741: the march's longest substep in the ensemble cells only")
    ap.add_argument("--ensemble-only", action="store_true", help="#741: replace the factorial cells with the opening-ice ensemble")
    ap.add_argument("--adopt-opening", type=float, nargs=2, default=None, metavar=("MYR", "FACTOR"), help="#741: scale the opening ice of the step closing at MYR, and keep the result")
    ap.add_argument("--hold-noise-from", type=float, default=None, help="#741: hold the drifting noise at this time (Myr) onwards")
    ap.add_argument("--substep-log-myr", type=float, default=None, help="#741: log every substep of this step, and capture any that moves the volume over 15%%")
    ap.add_argument("--capture-dir", type=Path, default=None, help="#741: where --factorial-myr writes its node fields and cells")


def _parse_value(text: str):
    if text in ("true", "false"):
        return text == "true"
    try:
        return float(text)
    except ValueError:
        return text


def freeze_climate_noise(world) -> None:
    """Seed climate.py's per-step generators with the opening time instead of the current one.
    Only calls made from climate.py itself are changed."""
    import sys

    import numpy as np_mod

    original = np_mod.random.default_rng
    opening = round(world.elapsed_years)

    def frozen(seed=None):
        caller = sys._getframe(1).f_globals.get("__name__", "")
        if caller.endswith("hydroclimate.climate") and isinstance(seed, tuple) and len(seed) >= 2 and seed[0] == world.seed:
            seed = (seed[0], opening, *seed[2:])
        return original(seed)

    np_mod.random.default_rng = frozen


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
        # A march's substeps are short, so the nearest root is the one the ice evolves to.
        g0.enumerate_roots = not args.no_enumerate and not g0.march
        erosion.GLACIER_ABRASION_COEFFICIENT = 0.0 if args.arm in NO_EROSION_ARMS else 100.0
        erosion.GLACIER_ABRASION_EXPONENT = 0.65
        if args.arm in NO_EROSION_ARMS:
            erosion.GLACIER_FLATTEN_RATE_PER_MYR = 0.0
    if g0.hold_sea_level:
        eustasy.update_sea_level = lambda world_: None
    if g0.hold_deflection:
        for plate in world.plates:
            plate.set_fields_on_plate(g0_hold_deflection_m=plate.collect("ice_load_deflection_m").copy())
    if args.ice_age_period_years is not None:
        world.ice_age_period_years = args.ice_age_period_years
    if args.arm in FREEZE_NOISE_ARMS or args.freeze_climate_noise:
        freeze_climate_noise(world)
    return {
        "arm": args.arm,
        **{k: v for k, v in vars(g0).items() if k != "frictional_heat_w_m2"},
        "abrasion_coefficient": erosion.GLACIER_ABRASION_COEFFICIENT,
        "abrasion_exponent": erosion.GLACIER_ABRASION_EXPONENT,
        "ice_age_period_years": float(getattr(world, "ice_age_period_years", 0.0) or 0.0),
        "climate_noise_frozen": args.arm in FREEZE_NOISE_ARMS or args.freeze_climate_noise,
    }


def _align(field, old_points, new_points):
    """A previous step's node field on this step's nodes: by index when the node set is the same
    size, else from the nearest old node."""
    if field is None or len(old_points) == len(new_points):
        return field
    from scipy.spatial import cKDTree

    return field[cKDTree(old_points).query(new_points)[1]]


def substep_logger(args: argparse.Namespace, step_ice, areas):
    """#741: wrap step_ice for one march. Writes one JSON line per substep, and every node field
    of the solve for any substep that moves the volume by over 15%, with the substep before it."""
    log = (args.capture_dir / "substeps.jsonl").open("w")
    count = {"n": 0}
    before: dict = {}

    def fields(ice, depth_in, depth_out, melt) -> dict:
        out = {"depth_in": depth_in, "depth_out": depth_out, "melt": melt}
        for key, value in vars(ice).items():
            if isinstance(value, np.ndarray) and value.dtype.kind in "fb" and value.shape[:1] == depth_in.shape:
                out[key] = value
        return {k: np.array(v, dtype=bool if v.dtype.kind == "b" else np.float32) for k, v in out.items()}

    def wrapped(*a, **k):
        depth_in = a[6]
        dt = a[-1]
        new_depth, melt, ice = step_ice(*a, **k)
        volume_in, volume_out = float(np.dot(depth_in, areas)) / 1e9, float(np.dot(new_depth, areas)) / 1e9
        count["n"] += 1
        row = {"substep": count["n"], "dt": dt, "volume_in_km3": volume_in, "volume_out_km3": volume_out,
               "visible_nodes": int(np.count_nonzero(new_depth >= 1.0)), "max_depth_m": float(new_depth.max()),
               "max_sliding": float(ice.sliding_speed_m_per_yr.max()) if len(ice.sliding_speed_m_per_yr) else 0.0,
               "warm_nodes": int(np.count_nonzero(ice.warm_based)),
               "ledger": {key: value / 1e9 for key, value in vars(ice.ledger).items() if isinstance(value, float)}}
        log.write(json.dumps(row) + "\n")
        log.flush()
        now = fields(ice, depth_in, new_depth, melt)
        if abs(volume_out / max(volume_in, 1.0) - 1.0) > 0.15:
            np.savez_compressed(args.capture_dir / f"substep-{count['n']:04d}.npz", **now, **{f"prev_{k}": v for k, v in before.items()})
        before.clear()
        before.update(now)
        return new_depth, melt, ice

    return wrapped


def install_factorial(args: argparse.Namespace, world) -> None:
    """#741: on each step named by --factorial-myr, re-march that step from its real opening ice
    with last step's climate, last step's surface, both, and +-0.1% opening ice. All in-process,
    so a save round trip can't move the state. Cell `new_new` repeats the real march and must
    match it exactly: that is the check that a cell leaves nothing behind."""
    from app.hydroclimate import hydrology, ice_sheet  # noqa: E402

    real_march = hydrology._g0_march_ice
    wanted = [round(t * 1e6) for t in args.factorial_myr]
    if args.capture_dir is not None:
        args.capture_dir.mkdir(parents=True, exist_ok=True)
    previous: dict = {}

    def summary(depth, ice, areas, log) -> dict:
        path = np.array(log["volume_full_km3"], dtype=float)
        visible = depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
        return {
            "volume_km3": float(np.dot(depth, areas)) / 1e9,
            "visible_area_km2": float(areas[visible].sum()) / 1e6,
            "accumulation_km3": ice.ledger.accumulation_m3 / 1e9,
            "surface_melt_km3": ice.ledger.surface_melt_m3 / 1e9,
            "sea_level_m": log["sea_level_m"],
            "substeps": log["substeps"],
            "flips": log["flips"],
            "settled": log["settled"],
            "path_km3": path.tolist(),
            "worst_substep_change": float((path[1:] / path[:-1] - 1.0).min()) if len(path) > 1 else 0.0,
        }

    def wrapped(world_, points, neighbor_idx, elevation, deflection, areas, prev_depth, precipitation,
                temperature, amplitude, heat, mobile_cover, is_ocean, allow_ocean_ice, plates_in_order, years):
        fields = {"points": points, "precipitation": precipitation, "temperature": temperature, "amplitude": amplitude,
                  "elevation": elevation, "deflection": deflection, "mobile_cover": mobile_cover, "is_ocean": is_ocean,
                  "allow_ocean_ice": allow_ocean_ice}
        heat_snapshot = None if heat is None else heat.copy()
        if args.adopt_opening is not None and abs(world_.elapsed_years - args.adopt_opening[0] * 1e6) < 1000:
            prev_depth = prev_depth * args.adopt_opening[1]

        def march(cl, su, opening):
            return real_march(world_, points, neighbor_idx, su["elevation"], su["deflection"], areas, opening.copy(),
                              cl["precipitation"], cl["temperature"], cl["amplitude"],
                              None if heat_snapshot is None else heat_snapshot.copy(), su["mobile_cover"], su["is_ocean"],
                              su["allow_ocean_ice"], plates_in_order, years)

        logging_step = args.substep_log_myr is not None and abs(world_.elapsed_years - args.substep_log_myr * 1e6) < 1000
        plain_step_ice = ice_sheet.step_ice
        if logging_step:
            ice_sheet.step_ice = substep_logger(args, plain_step_ice, areas)
        depth, melt, ice = real_march(world_, points, neighbor_idx, elevation, deflection, areas, prev_depth, precipitation,
                                      temperature, amplitude, heat, mobile_cover, is_ocean, allow_ocean_ice,
                                      plates_in_order, years)
        ice_sheet.step_ice = plain_step_ice
        target = next((t for t in wanted if abs(world_.elapsed_years - t) < 1000), None)
        if target is not None and previous:
            real_log = ice_sheet.G0_MARCH_LOG[-1]
            kept_log = len(ice_sheet.G0_MARCH_LOG)
            kept_heat = ice_sheet.G0.frictional_heat_w_m2
            old = {k: _align(v, previous["points"], points) for k, v in previous.items() if k != "points"}
            new = {k: v for k, v in fields.items() if k != "points"}
            cells = {"real": summary(depth, ice, areas, real_log)}
            saved = {"depth_real": depth.astype(np.float32), "opening": prev_depth.astype(np.float32)}
            plan = {
                "new_new": (new, new, prev_depth), "old_new": (old, new, prev_depth), "new_old": (new, old, prev_depth),
                "old_old": (old, old, prev_depth), "ice_plus": (new, new, prev_depth * 1.001),
                "ice_minus": (new, new, prev_depth * 0.999),
            }
            if args.ensemble_only:
                plan = {"new_new": (new, new, prev_depth)}
                for percent in (0.01, 0.03, 0.1, 0.2, 0.3, 0.5):
                    for sign in (1.0, -1.0):
                        factor = 1.0 + sign * percent / 100.0
                        plan[f"new_{factor:.4f}"] = (new, new, prev_depth * factor)
                        if args.ensemble_old_climate:
                            plan[f"old_{factor:.4f}"] = (old, new, prev_depth * factor)
            longest = ice_sheet.G0.march_max_substep_years
            if args.ensemble_max_substep_years is not None:
                ice_sheet.G0.march_max_substep_years = args.ensemble_max_substep_years
            for name, (cl, su, opening) in plan.items():
                started = time.perf_counter()
                cell_depth, _, cell_ice = march(cl, su, opening)
                cells[name] = summary(cell_depth, cell_ice, areas, ice_sheet.G0_MARCH_LOG[-1])
                cells[name]["seconds"] = round(time.perf_counter() - started, 1)
                saved[f"depth_{name}"] = cell_depth.astype(np.float32)
                del ice_sheet.G0_MARCH_LOG[kept_log:]
                print(json.dumps({"target": target, "cell": name, **{k: v for k, v in cells[name].items() if k != "path_km3"}}), flush=True)
            ice_sheet.G0.march_max_substep_years = longest
            ice_sheet.G0.frictional_heat_w_m2 = kept_heat
            ice_sheet.G0.reference_depth = None
            stem = args.capture_dir / f"step-{round(target / 1e4):05d}"
            np.savez_compressed(
                stem.with_suffix(".npz"), points=points.astype(np.float32), areas=areas.astype(np.float32),
                precipitation=precipitation.astype(np.float32), temperature=temperature.astype(np.float32),
                old_precipitation=old["precipitation"].astype(np.float32), old_temperature=old["temperature"].astype(np.float32),
                elevation=elevation.astype(np.float32), old_elevation=old["elevation"].astype(np.float32),
                deflection=deflection.astype(np.float32), is_ocean=is_ocean, old_is_ocean=old["is_ocean"], **saved,
            )
            stem.with_suffix(".json").write_text(json.dumps({"target_myr": target / 1e6, "opening_sea_level_m": float(world_.sea_level_m), "cells": cells}, indent=1))
        previous.clear()
        previous.update({k: None if v is None else np.array(v, copy=True) for k, v in fields.items()})
        return depth, melt, ice

    hydrology._g0_march_ice = wrapped


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


def sliding_report(ice, areas) -> dict:
    """Where the diagnosed sliding is fast: the iced area above each speed, and the thickness
    and balance flux at the fastest 1% of sliding nodes."""
    iced = ice.thickness > 0.0
    speed = ice.sliding_speed_m_per_yr
    total = float(areas[iced].sum()) or 1.0
    out = {f"area_share_above_{int(v)}_m_per_yr": round(float(areas[iced & (speed > v)].sum()) / total, 5) for v in (100, 1000, 10000)}
    sliding = np.flatnonzero(speed > 0.0)
    if len(sliding):
        top = sliding[np.argsort(-speed[sliding])][: max(1, len(sliding) // 100)]
        out.update(
            top_speed_median_m_per_yr=float(np.median(speed[top])),
            top_thickness_median_m=float(np.median(ice.thickness[top])),
            top_flux_median_m3_per_yr=float(np.median(ice.balance_flux_m3_per_yr[top])),
            top_nodes=int(len(top)),
        )
    return out


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

    for salt in ("MIXING_NOISE_SALT", "HUMIDITY_NOISE_SALT", "MFC_NOISE_SALT"):
        if hasattr(climate, salt):
            setattr(climate, salt, getattr(climate, salt) + args.noise_salt_offset)

    result["config"] = configure(args, world)
    if args.factorial_myr or args.adopt_opening is not None or args.substep_log_myr is not None:
        install_factorial(args, world)
    if args.hold_noise_from is not None:
        drifting = climate._drifting_white_noise
        held = args.hold_noise_from * 1e6
        climate._drifting_white_noise = lambda seed, salt, elapsed, shape: drifting(seed, salt, min(elapsed, held), shape)
    tap = {"ice_seconds": 0.0}
    original_step_ice = ice_sheet.step_ice

    def timed_step_ice(*a, **k):
        started = time.perf_counter()
        try:
            return original_step_ice(*a, **k)
        finally:
            tap["ice_seconds"] += time.perf_counter() - started

    ice_sheet.step_ice = timed_step_ice
    window = {"seconds": 0.0, "ice_seconds": 0.0, "ledger": {}, "worst_residual_m3": 0.0, "roots": [], "membership": [], "outer": [], "march": []}

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
        window["march"].extend(ice_sheet.G0_MARCH_LOG)
        ice_sheet.G0_MARCH_LOG.clear()
        hydro = world.hydrology_cache
        ice = getattr(hydro, "ice", None) if hydro is not None else None
        if ice is not None:
            for key, value in vars(ice.ledger).items():
                if key != "capped_nodes":
                    window["ledger"][key] = window["ledger"].get(key, 0.0) + value / 1e9
            window["worst_residual_m3"] = max(window["worst_residual_m3"], abs(ice.ledger.residual_m3))
            window["warm_based_nodes"] = int(np.count_nonzero(ice.warm_based))
            if args.sliding_report:
                areas_now = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
                window["sliding"] = sliding_report(ice, areas_now)
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
            "march": window["march"][-1] if window["march"] else None,
            "march_unsettled_steps": sum(not r["settled"] for r in window["march"]),
            "sliding": window.pop("sliding", None),
        }
        window.update(seconds=0.0, ice_seconds=0.0, ledger={}, worst_residual_m3=0.0, roots=[], membership=[], outer=[], march=[])
        return closed

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
