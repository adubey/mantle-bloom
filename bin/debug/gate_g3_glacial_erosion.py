#!/usr/bin/env python3
"""RFC 0015 gate G3 follow-up (#973): integrated glacial abrasion and its limits.

Replays a save with `step_world` under one of three glacial abrasion laws:

- `legacy-law`: D2's legacy coefficient expressed in the new m/yr law, exponent = 1.0
- `cook2020`: D3 start (`K_G = 1e-4`, `l = 0.65`, Cook et al., 2020)
- `humphrey1994`: D3 alternative (`K_G = 1e-4`, `l = 1.0`, Humphrey & Raymond, 1994)

On each step, records:
- raw glacial abrasion volume (`m3`), applied glacial abrasion volume (`m3`), and their ratio;
- glacial flattening volume (`m3`);
- total subaerial erosion volume (`m3`, from `removed_m3` minus marine/coastal/leveling/flatten,
  or `erosion_amount`);
- end-of-step ice volume (`km3`), visible ice area (`km2`), land ice area (`km2`),
  `land_volume_above_sea_km3`, and `sea_level_m`.

Supports `--multiplier` (`World.glacier_erosion_multiplier`) and an optional
`--flatten-multiplier` override to separate flattening's share of the knob's effect.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_g3_glacial_erosion.py \
        <save> --arm cook2020 --steps 40 --myr 0.25 --out <captures>/run.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

BACKEND = Path(os.environ.get("MANTLE_BLOOM_BACKEND", str(Path(__file__).resolve().parents[2] / "backend")))
DEBUG_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BACKEND))
sys.path.insert(0, str(DEBUG_DIR))

from app import cratons, persistence, stats  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.hydroclimate import erosion, glacier_erosion, hydrology, ice_sheet  # noqa: E402
from sweep_lib import _land_volume_above_sea_km3  # noqa: E402

ARMS: dict[str, tuple[float, float]] = {
    "legacy-law": (430.0 * hydrology.GLACIER_ACCUMULATION_RATE / (1000.0 / 917.0 * 1.0e-3 * 1.0e6) / 1.0e6, 1.0),
    "cook2020": (1.0e-4, 0.65),
    "humphrey1994": (1.0e-4, 1.0),
}


def _snapshot_state(world) -> dict:
    """Measure ice-cap area, ice volume, land volume above sea level, and sea level."""
    spacing_rad = line_spacing_rad(world.node_density)
    weights = stats.overlap_area_weights(world)
    ice_vol_m3 = 0.0
    vis_ice_area_m2 = 0.0
    land_ice_area_m2 = 0.0
    land_area_m2 = 0.0
    for plate in world.plates:
        if plate.plate_id not in weights:
            continue
        _, elevation = plate.all_points_and_elevation()
        glacier = np.clip(plate.collect("glacier_depth"), 0.0, None)
        area_w = plate.accounting_areas_m2(spacing_rad) * weights[plate.plate_id]
        ice_vol_m3 += float(np.dot(glacier, area_w))
        vis = glacier >= hydrology.GLACIER_VISIBLE_DEPTH_M
        land = elevation > world.sea_level_m
        vis_ice_area_m2 += float(area_w[vis].sum())
        land_ice_area_m2 += float(area_w[vis & land].sum())
        land_area_m2 += float(area_w[land].sum())
    return {
        "age_myr": float(world.elapsed_years) / 1.0e6,
        "sea_level_m": float(world.sea_level_m),
        "land_volume_above_sea_km3": _land_volume_above_sea_km3(world),
        "land_area_km2": land_area_m2 / 1.0e6,
        "ice_volume_km3": ice_vol_m3 / 1.0e9,
        "visible_ice_area_km2": vis_ice_area_m2 / 1.0e6,
        "land_ice_area_km2": land_ice_area_m2 / 1.0e6,
    }


def _run_replay(ns: argparse.Namespace) -> dict:
    raw = ns.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    if ns.ice_age_period_years is not None:
        world.ice_age_period_years = ns.ice_age_period_years
    world.glacier_erosion_multiplier = ns.multiplier

    arm_cfg = ARMS[ns.arm]
    orig_coeff = glacier_erosion.GLACIAL_ABRASION_COEFFICIENT
    orig_exp = glacier_erosion.GLACIAL_ABRASION_EXPONENT
    glacier_erosion.GLACIAL_ABRASION_COEFFICIENT, glacier_erosion.GLACIAL_ABRASION_EXPONENT = arm_cfg

    preset = next((k for k, v in ice_sheet.MARCH_COST_PRESETS.items() if v == ice_sheet.MARCH_COST), "custom")
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=BACKEND, check=True, capture_output=True, text=True
    ).stdout.strip()

    step_capture: dict[str, float] = {}
    orig_flatten = erosion._flatten
    orig_apply_erosion = erosion.apply_erosion

    def _wrapped_flatten(hydro, ice_factor, years, removable_m=None, multiplier=1.0):
        eff_mult = ns.flatten_multiplier if ns.flatten_multiplier is not None else multiplier
        send = orig_flatten(hydro, ice_factor, years, removable_m=removable_m, multiplier=eff_mult)
        plates_in_order = [p for p in world.plates if p.node_count() > 0]
        area = erosion._gather_areas(world, plates_in_order)
        craton_keep = 1.0 - cratons.CRATON_EROSION_RESISTANCE * cratons.strength(
            np.concatenate([p.collect("craton_crust_m") for p in plates_in_order])
        )
        flatten_removed = (send * craton_keep[:, None]).sum(axis=1)
        raw_glacier = (
            np.asarray(hydro.ice.raw_abrasion_depth_m, dtype=float)
            if hydro.ice is not None and hydro.ice.raw_abrasion_depth_m is not None
            else np.zeros(len(hydro.elevation))
        )
        step_capture["raw_glacial_erosion_m3"] = float(np.sum(raw_glacier * area))
        step_capture["glacial_flatten_m3"] = float(np.sum(flatten_removed * area))
        return send

    def _wrapped_apply_erosion(w, years, **kwargs):
        res = orig_apply_erosion(w, years, **kwargs)
        step_capture["applied_glacial_erosion_m3"] = float(res.budget["glacial_erosion_m3"])
        step_capture["unapplied_glacial_erosion_m3"] = float(res.budget["glacial_erosion_unapplied_m3"])
        step_capture["total_removed_m3"] = float(res.budget["removed_m3"])
        return res

    erosion._flatten = _wrapped_flatten
    erosion.apply_erosion = _wrapped_apply_erosion

    try:
        opening = _snapshot_state(world)
        steps = []
        t0 = time.perf_counter()
        for i in range(ns.steps):
            step_capture.clear()
            s0 = time.perf_counter()
            world_mod.step_world(world, ns.myr * 1.0e6)
            wall = time.perf_counter() - s0
            state = _snapshot_state(world)
            raw_m3 = step_capture["raw_glacial_erosion_m3"]
            app_m3 = step_capture["applied_glacial_erosion_m3"]
            flat_m3 = step_capture["glacial_flatten_m3"]
            ratio = app_m3 / raw_m3 if raw_m3 > 0.0 else 1.0
            rec = {
                "step": i,
                "wall_s": round(wall, 2),
                **state,
                "raw_glacial_erosion_m3": raw_m3,
                "applied_glacial_erosion_m3": app_m3,
                "applied_raw_ratio": ratio,
                "glacial_flatten_m3": flat_m3,
                "total_removed_m3": step_capture["total_removed_m3"],
            }
            steps.append(rec)
            print(
                f"step {i:2d}/{ns.steps}: {wall:4.1f}s age={state['age_myr']:6.2f} "
                f"raw={raw_m3/1e9:9.2f} app={app_m3/1e9:8.2f} km3 "
                f"ratio={ratio:.4f} flat={flat_m3/1e9:8.2f} km3 "
                f"ice_area={state['land_ice_area_km2']/1e6:.2f}M land_vol={state['land_volume_above_sea_km3']/1e6:.2f}M",
                flush=True,
            )
    finally:
        erosion._flatten = orig_flatten
        erosion.apply_erosion = orig_apply_erosion
        glacier_erosion.GLACIAL_ABRASION_COEFFICIENT = orig_coeff
        glacier_erosion.GLACIAL_ABRASION_EXPONENT = orig_exp

    total_raw = sum(s["raw_glacial_erosion_m3"] for s in steps)
    total_app = sum(s["applied_glacial_erosion_m3"] for s in steps)
    total_flat = sum(s["glacial_flatten_m3"] for s in steps)
    closing = steps[-1] if steps else opening
    rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        rss_bytes *= 1024
    return {
        "issue": 973,
        "commit": commit,
        "save": ns.save.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "march_preset": preset,
        "arm": ns.arm,
        "coefficient_m_per_yr": arm_cfg[0],
        "exponent": arm_cfg[1],
        "glacier_erosion_multiplier": ns.multiplier,
        "flatten_multiplier": ns.flatten_multiplier if ns.flatten_multiplier is not None else ns.multiplier,
        "ice_age_period_years": float(world.ice_age_period_years),
        "steps": ns.steps,
        "myr_per_step": ns.myr,
        "total_myr": ns.steps * ns.myr,
        "elapsed_s": round(time.perf_counter() - t0, 2),
        "peak_rss_gb": round(rss_bytes / 1e9, 2),
        "summary": {
            "cumulative_raw_glacial_erosion_km3": total_raw / 1.0e9,
            "cumulative_applied_glacial_erosion_km3": total_app / 1.0e9,
            "cumulative_applied_raw_ratio": total_app / total_raw if total_raw > 0.0 else 1.0,
            "cumulative_glacial_flatten_km3": total_flat / 1.0e9,
            "flatten_to_applied_abrasion_ratio": total_flat / total_app if total_app > 0.0 else None,
            "opening_land_ice_area_km2": opening["land_ice_area_km2"],
            "closing_land_ice_area_km2": closing["land_ice_area_km2"],
            "opening_visible_ice_area_km2": opening["visible_ice_area_km2"],
            "closing_visible_ice_area_km2": closing["visible_ice_area_km2"],
            "opening_land_volume_above_sea_km3": opening["land_volume_above_sea_km3"],
            "closing_land_volume_above_sea_km3": closing["land_volume_above_sea_km3"],
            "opening_ice_volume_km3": opening["ice_volume_km3"],
            "closing_ice_volume_km3": closing["ice_volume_km3"],
        },
        "opening": opening,
        "per_step": steps,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--arm", choices=tuple(ARMS), required=True)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--myr", type=float, default=1.0)
    ap.add_argument("--multiplier", type=float, default=1.0)
    ap.add_argument("--flatten-multiplier", type=float, default=None)
    ap.add_argument("--ice-age-period-years", type=float, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ns = ap.parse_args()

    ns.out.parent.mkdir(parents=True, exist_ok=True)
    result = _run_replay(ns)
    ns.out.write_text(json.dumps(result, indent=1) + "\n")
    s = result["summary"]
    print(
        f"wrote {ns.out} | app={s['cumulative_applied_glacial_erosion_km3']:.1f} km3 "
        f"raw={s['cumulative_raw_glacial_erosion_km3']:.1f} km3 "
        f"app/raw={s['cumulative_applied_raw_ratio']:.4f} "
        f"flat={s['cumulative_glacial_flatten_km3']:.1f} km3 "
        f"rss={result['peak_rss_gb']} GB"
    )


if __name__ == "__main__":
    main()
