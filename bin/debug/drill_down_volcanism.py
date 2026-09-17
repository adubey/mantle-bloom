#!/usr/bin/env python3
"""Single-seed, step-by-step drill-down for issue #173 ("volcanism multiplier trends inverse:
more volcanism -> less land"). Runs one seed at volcanism_multiplier=1.0 and 4.0 side by side,
instrumenting erosion.apply_erosion and volcanism.apply_volcanic_activity to record, every
step: total_continental_crust_volume_km3, land_volume_above_sea_km3, how much elevation
volcanism added this step, and how much net erosion landed on is_volcano-tagged nodes vs
everywhere else (is_volcano is permanent provenance -- see elevation_lines.py -- so this
isolates "does land erosion disproportionately target ground a volcano has ever touched").

Usage: backend/.venv/bin/python bin/debug/drill_down_volcanism.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import _land_volume_above_sea_km3  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import erosion, lithosphere, plates as plates_mod, stats, volcanism  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

SEED = 829071382
NODE_DENSITY = 1.0
STEP_YEARS = 10_000_000
N_STEPS = 12  # 120 My


def _total_isostatic_debt_km3(world) -> float:
    """Total volume (km^3, positive part only) by which every node's actual `elevation`
    exceeds what its own crustal_thickness_m/mantle_lithosphere_thickness_m would isostatically
    support -- i.e. "surface elevation not backed by any crustal mass." volcanism.py only ever
    writes `elevation` (see apply_volcanic_activity/_spread_volcanic_plains), never
    crustal_thickness_m, so every eruption/apron bump is untethered debt of exactly this kind
    until something either raises Hc under it or erosion planes the surface back down to match.
    Reuses lithosphere.isostatic_elevation, the same Airy formula erosion.py's own Hc-driven
    correction is derived from."""
    from app.elevation_lines import line_spacing_rad

    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    total_m3 = 0.0
    for plate in world.plates:
        if not plate.lines:
            continue
        hc = plate.collect("crustal_thickness_m")
        hm = plate.collect("mantle_lithosphere_thickness_m")
        elevation = plate.collect("elevation")
        rho_c = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
        equilibrium = lithosphere.isostatic_elevation(hc, hm, rho_c)
        excess = np.clip(elevation - equilibrium, 0.0, None)
        total_m3 += float(np.sum(excess)) * area_m2
    return total_m3 / 1.0e9


def run_one(multiplier: float) -> list[dict]:
    world = generate_world(seed=SEED, node_density=NODE_DENSITY)
    world.volcanism_multiplier = multiplier

    telemetry = {"volcanic_added_m": 0.0, "erosion_volcano_m": 0.0, "erosion_other_m": 0.0,
                 "n_volcano": 0, "n_other": 0}

    orig_apply_erosion = erosion.apply_erosion
    orig_apply_volcanic_activity = volcanism.apply_volcanic_activity

    def wrapped_apply_erosion(world_, years, node_cloud=None):
        if node_cloud is not None:
            points, plates_in_order = node_cloud
            is_volcano = plates_mod.collect_all_is_volcano(plates_in_order)
        else:
            is_volcano = None
        result = orig_apply_erosion(world_, years, node_cloud=node_cloud)
        if result is not None and is_volcano is not None and len(is_volcano) == len(result.net_elevation_change_m):
            net = result.net_elevation_change_m
            lowering = np.clip(-net, 0.0, None)  # positive where net erosion, zero where net deposition
            telemetry["erosion_volcano_m"] += float(lowering[is_volcano].sum())
            telemetry["erosion_other_m"] += float(lowering[~is_volcano].sum())
            telemetry["n_volcano"] = int(is_volcano.sum())
            telemetry["n_other"] = int((~is_volcano).sum())
        return result

    def wrapped_apply_volcanic_activity(world_, years):
        before = plates_mod.collect_all_elevation(world_.plates)
        orig_apply_volcanic_activity(world_, years)
        after = plates_mod.collect_all_elevation(world_.plates)
        telemetry["volcanic_added_m"] += float(np.clip(after - before, 0.0, None).sum())

    erosion.apply_erosion = wrapped_apply_erosion
    volcanism.apply_volcanic_activity = wrapped_apply_volcanic_activity

    records = []
    try:
        for step in range(1, N_STEPS + 1):
            telemetry["volcanic_added_m"] = 0.0
            telemetry["erosion_volcano_m"] = 0.0
            telemetry["erosion_other_m"] = 0.0
            step_world(world, years=STEP_YEARS)
            snapshot = stats.compute_stats(world)
            records.append({
                "step": step,
                "years": step * STEP_YEARS,
                "total_continental_crust_volume_km3": snapshot["total_continental_crust_volume_km3"],
                "land_volume_above_sea_km3": _land_volume_above_sea_km3(world),
                "isostatic_debt_km3": _total_isostatic_debt_km3(world),
                "land_fraction": snapshot["land_fraction"],
                "volcanic_added_m": telemetry["volcanic_added_m"],
                "erosion_volcano_m": telemetry["erosion_volcano_m"],
                "erosion_other_m": telemetry["erosion_other_m"],
                "n_volcano": telemetry["n_volcano"],
                "n_other": telemetry["n_other"],
            })
    finally:
        erosion.apply_erosion = orig_apply_erosion
        volcanism.apply_volcanic_activity = orig_apply_volcanic_activity

    return records


def main() -> None:
    results = {}
    for m in (1.0, 4.0):
        print(f"=== running seed={SEED} volcanism_multiplier={m} ===", file=sys.stderr)
        results[m] = run_one(m)

    header = (
        f"{'step':>4} {'crust_vol_km3':>14} {'land_vol_km3':>14} {'isostatic_debt_km3':>18} {'land%':>6} "
        f"{'volc_added_m':>13} {'eros@volc_m':>12} {'eros@other_m':>13} "
        f"{'mean_eros@volc':>15} {'mean_eros@other':>16} {'n_volc':>7} {'n_other':>8}"
    )
    for m in (1.0, 4.0):
        print(f"\n--- multiplier {m} ---")
        print(header)
        for r in results[m]:
            mean_v = r["erosion_volcano_m"] / r["n_volcano"] if r["n_volcano"] else 0.0
            mean_o = r["erosion_other_m"] / r["n_other"] if r["n_other"] else 0.0
            print(
                f"{r['step']:>4} {r['total_continental_crust_volume_km3']:>14.4e} "
                f"{r['land_volume_above_sea_km3']:>14.4e} {r['isostatic_debt_km3']:>18.4e} {r['land_fraction']*100:>5.1f}% "
                f"{r['volcanic_added_m']:>13.3e} {r['erosion_volcano_m']:>12.3e} {r['erosion_other_m']:>13.3e} "
                f"{mean_v:>15.4f} {mean_o:>16.4f} {r['n_volcano']:>7} {r['n_other']:>8}"
            )

    import json
    out_path = Path(__file__).resolve().parent / "results" / "drill_down_volcanism.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump({str(k): v for k, v in results.items()}, f, indent=2)
    print(f"\nWrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
