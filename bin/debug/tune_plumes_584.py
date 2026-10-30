#!/usr/bin/env python3
"""RFC 0007 tuning + sensitivity harness (issue #584).

One world at a time: baselines on seeds 7/11/2024 to 1 Gyr, regenerated seed
456862899 to 3 Gyr, and single-knob sensitivity cases (seed 7, 500 Myr unless
noted). Every tap wraps, records, delegates, and is restored; taps never
mutate world state. Split counting matches #407/#582 exactly
(`maybe_split_plate` successes only).

Usage:
    backend/.venv/bin/python bin/debug/tune_plumes_584.py --list
    backend/.venv/bin/python bin/debug/tune_plumes_584.py --case baseline_7 \\
        --out analysis/issue584/tuning.json
    backend/.venv/bin/python bin/debug/tune_plumes_584.py --all \\
        --out analysis/issue584/tuning.json   # serial, hours-long

One world runs at a time (replays are RAM-heavy). Only the compact JSON is
reviewed; no saves are written.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
from app import continental_ledger, merge_split, plumes, stats  # noqa: E402
from app import terranes as terranes_mod  # noqa: E402
from app import world as world_mod  # noqa: E402

CHECKPOINT_EVERY_YEARS = 1e8

SCRIPT_NAME = Path(__file__).name
STEP_YEARS = 5e6

CASES: list[dict] = [
    {"name": "baseline_7", "seed": 7, "years": 1e9, "knobs": {}},
    {"name": "baseline_11", "seed": 11, "years": 1e9, "knobs": {}},
    {"name": "baseline_2024", "seed": 2024, "years": 1e9, "knobs": {}},
    {"name": "regen_456862899_3gyr", "seed": 456862899, "years": 3e9, "knobs": {}},
    # Sensitivity: one knob at a time off the seed-7 baseline, 500 Myr.
    {"name": "sens_rate_0.05", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_RATE_PER_MYR": 0.05}},
    {"name": "sens_rate_0.2", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_RATE_PER_MYR": 0.2}},
    {"name": "sens_cox_on", "seed": 7, "years": 5e8,
     "knobs": {"world.plume_cox_modulation": True}},
    {"name": "sens_vol_lognormal", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_VOLUME_DISTRIBUTION": "lognormal"}},
    {"name": "sens_vol_powerlaw", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_VOLUME_DISTRIBUTION": "power_law"}},
    {"name": "sens_dock_easy", "seed": 7, "years": 5e8,
     "knobs": {"terranes.DOCK_E_YOUNG_M": 2500.0, "terranes.DOCK_E_OLD_M": 10000.0}},
    {"name": "sens_dock_hard", "seed": 7, "years": 5e8,
     "knobs": {"terranes.DOCK_E_YOUNG_M": 8000.0, "terranes.DOCK_E_OLD_M": 20000.0}},
    {"name": "sens_dock_loss_0.4", "seed": 7, "years": 5e8,
     "knobs": {"terranes.DOCK_LOWER_CRUST_LOSS": 0.4}},
    {"name": "sens_plateau_0.5", "seed": 7, "years": 5e8,
     "knobs": {"terranes.PLATEAU_PRESERVED_FRACTION": 0.5}},
    {"name": "sens_heat_5myr", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAT_MYR": 5.0}},
    {"name": "sens_heat_50myr", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAT_MYR": 50.0}},
    {"name": "sens_weaken_0.5", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_CRATON_WEAKENING": 0.5}},
    {"name": "sens_weaken_1.0", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_CRATON_WEAKENING": 1.0}},
    {"name": "sens_stress_0", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_STRESS": 0.0}},
    {"name": "sens_stress_0.75", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_STRESS": 0.75}},
    {"name": "sens_stress_2", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_STRESS": 2.0}},
    {"name": "sens_stress_3", "seed": 7, "years": 5e8,
     "knobs": {"plumes.PLUME_HEAD_STRESS": 3.0}},
]

_MODULES = {"plumes": plumes, "terranes": terranes_mod}


def _arrivals_since(world, t0_years: float, t1_years: float) -> list[dict]:
    """Head arrivals born in (t0, t1], each with host judgment and volume."""
    live = [plate for plate in world.plates if plate.node_count() > 0]
    arrivals = []
    for plume in world.mantle_plumes:
        if plume.is_tail_only:
            continue
        if not t0_years < plume.born_years <= t1_years:
            continue
        continental, host = plumes._center_host(live, plume.position)
        arrivals.append(
            {
                "born": float(plume.born_years),
                "host": host,
                "continental": bool(continental),
                "volume_km3": float(plume.head_volume_m3) / 1e9,
            }
        )
    return arrivals


class TuningTaps:
    """Record head arrivals (volume, host type) and split parent heat cover."""

    def __init__(self) -> None:
        self.heads = []
        self.splits = []  # each: time, parent, parent_heated_fraction
        self._saved = []

    def install(self) -> None:
        heads, splits, saved = self.heads, self.splits, self._saved

        def replace(obj, name, make_wrapper):
            target = getattr(obj, name)
            saved.append((obj, name, target))
            setattr(obj, name, make_wrapper(target))

        def apply_wrapper(original):
            def apply(world, years, *args, **kwargs):
                out = original(world, years, *args, **kwargs)
                # Arrivals are derived after the step from the world's own
                # clock, so no pre-step timestamp needs capturing.
                heads.extend(
                    _arrivals_since(world, world.elapsed_years - years, world.elapsed_years)
                )
                return out

            return apply

        def split_wrapper(original):
            def maybe_split_plate(world, plate, *args, **kwargs):
                heated_fraction = _heated_fraction(world, plate)
                out = original(world, plate, *args, **kwargs)
                if out is not None:
                    splits.append(
                        {
                            "time": float(world.elapsed_years),
                            "parent": plate.plate_id,
                            "parent_heated_fraction": heated_fraction,
                        }
                    )
                return out

            return maybe_split_plate

        replace(plumes, "apply", apply_wrapper)
        replace(merge_split, "maybe_split_plate", split_wrapper)

    def uninstall(self) -> None:
        while self._saved:
            obj, name, original = self._saved.pop()
            setattr(obj, name, original)


def _heated_fraction(world, plate) -> float | None:
    """Share of the parent plate's area under live plume heat at split time."""
    try:
        spacing = stats.line_spacing_rad(world.node_density)
        areas = plate.accounting_areas_m2(spacing)
        heat = plate.collect("plume_heat_myr")
        total = float(sum(areas))
        if total <= 0.0:
            return None
        return round(float(sum(a for a, h in zip(areas, heat) if h > 0.0)) / total, 4)
    except Exception:
        return None


def _apply_knobs(knobs: dict) -> list[tuple[object, str, object]]:
    """setattr overrides; returns restore records (world.* handled by caller)."""
    saved = []
    for dotted, value in knobs.items():
        mod_name, attr = dotted.split(".", 1)
        if mod_name == "world":
            continue
        mod = _MODULES[mod_name]
        saved.append((mod, attr, getattr(mod, attr)))
        setattr(mod, attr, value)
    return saved


def _restore(saved) -> None:
    for obj, attr, value in reversed(saved):
        setattr(obj, attr, value)


def _checkpoint(world, taps) -> dict:
    """100 Myr snapshot: active hotspots, plate count, heat/plateau area."""
    census = plumes.diagnostics(world)
    spacing = stats.line_spacing_rad(world.node_density)
    heat_area_m2 = plateau_area_m2 = 0.0
    plates = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        plates += 1
        areas = plate.accounting_areas_m2(spacing)
        heat = plate.collect("plume_heat_myr")
        if np.any(heat > 0.0):
            heat_area_m2 += float(np.sum(areas[heat > 0.0]))
        plateau, _, _ = terranes_mod.plate_tracers(plate)
        if np.any(plateau > 0.0):
            plateau_area_m2 += float(np.sum(areas[plateau > 0.0]))
    return {
        "elapsed_years": float(world.elapsed_years),
        "active_heads": census["heads"],
        "active_tails": census["tails"],
        "plates": plates,
        "heads_so_far": len(taps.heads),
        "splits_so_far": len(taps.splits),
        "heat_area_km2": round(heat_area_m2 / 1.0e6, 1),
        "plateau_area_km2": round(plateau_area_m2 / 1.0e6, 1),
    }


def _run_case(case: dict) -> dict:
    seed, years, knobs = case["seed"], case["years"], case["knobs"]
    saved = _apply_knobs(knobs)
    taps = TuningTaps()
    taps.install()
    try:
        world = world_mod.generate_world(seed)
        for dotted, value in knobs.items():
            mod_name, attr = dotted.split(".", 1)
            if mod_name == "world":
                setattr(world, attr, value)
        t0 = time.time()
        step_times = []
        checkpoints = [_checkpoint(world, taps)]
        next_checkpoint = CHECKPOINT_EVERY_YEARS
        while world.elapsed_years < years - 1.0:
            s0 = time.time()
            world_mod.step_world(world, min(STEP_YEARS, years - world.elapsed_years))
            step_times.append(time.time() - s0)
            if world.elapsed_years + 1.0 >= next_checkpoint and next_checkpoint < years:
                checkpoints.append(_checkpoint(world, taps))
                next_checkpoint += CHECKPOINT_EVERY_YEARS
        wall = time.time() - t0
        full = stats.compute_stats(world)
        census = plumes.diagnostics(world)
        checkpoints.append(_checkpoint(world, taps))
    finally:
        taps.uninstall()
        _restore(saved)
    gyr = years / 1e9
    cont = [h for h in taps.heads if h["continental"]]
    ocn = [h for h in taps.heads if not h["continental"]]
    heated = [s["parent_heated_fraction"] for s in taps.splits
              if s["parent_heated_fraction"] is not None]
    volumes = sorted(h["volume_km3"] for h in taps.heads)
    deciles = (
        [round(float(np.percentile(volumes, q)), 1) for q in (0, 10, 25, 50, 75, 90, 100)]
        if volumes else None
    )
    juveniles = {
        key: round(float(world.continental_material_ledger.get(key, 0.0)) / 1.0e9, 3)
        for key in continental_ledger.JUVENILE_SUB_ACCOUNTS
    }
    return {
        "case": case["name"],
        "seed": seed,
        "years": years,
        "knobs": knobs,
        "wall_s": round(wall, 1),
        "mean_step_s": round(sum(step_times) / max(1, len(step_times)), 2),
        "lip_events_per_gyr": round(len(taps.heads) / gyr, 2),
        "lip_continental": len(cont),
        "lip_oceanic": len(ocn),
        "lip_volume_km3_per_gyr": round(sum(h["volume_km3"] for h in taps.heads) / gyr, 1),
        "lip_volume_continental_km3_per_gyr": round(sum(h["volume_km3"] for h in cont) / gyr, 1),
        "lip_volume_oceanic_km3_per_gyr": round(sum(h["volume_km3"] for h in ocn) / gyr, 1),
        "active_heads_end": census["heads"],
        "active_tails_end": census["tails"],
        "splits": len(taps.splits),
        "splits_parent_heated_fraction_mean": (
            round(sum(heated) / len(heated), 4) if heated else None
        ),
        "splits_parent_fully_heated": sum(1 for f in heated if f > 0.5),
        "splits_parent_unheated": sum(1 for f in heated if f == 0.0),
        "plume_lip_events_stat": full.get("plume_lip_events"),
        "plume_plateau_area_km2": full.get("plume_plateau_area_km2"),
        "plume_heat_area_km2": full.get("plume_heat_area_km2"),
        "buoyant_crust_ledger_km3": full.get("buoyant_crust_ledger_km3"),
        "juvenile_continental_km3": juveniles,
        "head_volume_deciles_km3": deciles,
        "checkpoints": checkpoints,
        "craton_keys": {k: v for k, v in full.items() if k.startswith("craton")},
    }


def _provenance() -> dict:
    root = Path(__file__).resolve().parents[2]

    def git(*args):
        try:
            return subprocess.run(
                ["git", "-C", str(root), *args],
                capture_output=True, text=True, timeout=30,
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return "unavailable"

    return {
        "script": SCRIPT_NAME,
        "git_head": git("rev-parse", "HEAD"),
        "git_status_short": git("status", "--short"),
        "python": platform.python_version(),
        "controls": {"steps": "default generate_world, cratons on, 5 Myr steps"},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--list", action="store_true")
    group.add_argument("--case", default=None)
    group.add_argument("--all", action="store_true")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    if args.list:
        for case in CASES:
            print(f"{case['name']}: seed={case['seed']} years={case['years']:.0e} knobs={case['knobs']}")
        return

    wanted = CASES if args.all else [c for c in CASES if c["name"] == args.case]
    if not wanted:
        raise SystemExit(f"unknown case: {args.case}")
    results = []
    for case in wanted:
        print(f"[{time.strftime('%H:%M:%S')}] start {case['name']}", flush=True)
        result = _run_case(case)
        results.append(result)
        print(
            f"[{time.strftime('%H:%M:%S')}] done {case['name']}: "
            f"LIPs/Gyr={result['lip_events_per_gyr']} splits={result['splits']} "
            f"wall={result['wall_s']}s",
            flush=True,
        )
    payload = {"provenance": _provenance(), "results": results}
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        if out.exists():
            prior = json.loads(out.read_text()).get("results", [])
            seen = {r["case"] for r in results}
            results = [r for r in prior if r["case"] not in seen] + results
            payload["results"] = results
        out.write_text(json.dumps(payload, indent=2))
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
