#!/usr/bin/env python3
"""Issue #198 follow-up to #178: is a volcano's own tagged terrain actually landing where
ERUPTION_ELEVATION_M/VOLCANIC_PLAIN_ELEVATION_M were tuned for it to (elevation_lines.py's own
2026-09-04 comment: "a volcano that does erupt a few times actually builds a lasting peak", the
fix for the pre-tuning state where "across a 188 My save, volcano nodes averaged *below* sea
level")? #178 backed every eruption/apron elevation gain with real crustal_thickness_m so it's
no longer eaten away by erosion's isostatic correction the way unbacked relief was -- this
checks, on current main, whether the terrain that actually erupted now reads as the intended
"lasting peak" (mean elevation comfortably above sea level, matching the up-to-1500m-gross-
relief-per-volcano the constants describe) or as something that's overshot that intent.

`is_volcano` (plates.collect_all_is_volcano) is *permanent field provenance* -- every node
that was ever part of a rift-spawned volcanic field, the large majority of which is ordinary
subsided oceanic crust that never itself erupted (most nodes near a rift's volcanic field
never roll the per-node eruption chance in volcanism.py before the field's plate moves on) --
not "this specific node erupted." A first pass of this script measured is_volcano's mean
elevation directly and got ~-5200m, deeply submerged -- but that's dominated by that ordinary
subsided-crust population, not eruption sites, and doesn't speak to whether
ERUPTION_ELEVATION_M/VOLCANIC_PLAIN_ELEVATION_M themselves are well-tuned. `mineral_deposit_m`
only grows on the exact per-node eruption roll that also adds ERUPTION_ELEVATION_M (see
volcanism.py's MINERAL_DEPOSIT_PER_ERUPTION_M), so `mineral_deposit_m > 0` isolates "has this
node itself actually erupted at least once" -- the right population to check the constants
against. `volcano_active_years_remaining > 0` (currently active, may not have erupted yet) is
tracked alongside for context.

Usage: backend/.venv/bin/python bin/debug/measure_volcano_node_elevation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import SWEEP_SEEDS  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import lithosphere, plates as plates_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

NODE_DENSITY = 1.0
STEP_YEARS = 10_000_000
CHECKPOINTS = (30_000_000, 60_000_000, 90_000_000, 120_000_000)


def measure_one(seed: int) -> list[dict]:
    world = generate_world(seed=seed, node_density=NODE_DENSITY)
    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    records = []
    years_done = 0.0
    for checkpoint in CHECKPOINTS:
        while years_done < checkpoint:
            step = min(STEP_YEARS, checkpoint - years_done)
            step_world(world, years=step)
            years_done += step
        is_volcano = plates_mod.collect_all_is_volcano(world.plates)
        elevation = plates_mod.collect_all_elevation(world.plates)
        mineral_deposit = plates_mod.collect_all_mineral_deposit(world.plates)
        active_remaining = plates_mod.collect_all_volcano_active_years_remaining(world.plates)
        erupted = mineral_deposit > 0.0
        active = is_volcano & (active_remaining > 0.0)

        def _summary(mask: np.ndarray) -> dict:
            sub = elevation[mask]
            above = sub[sub > world.sea_level_m]
            return {
                "n": int(mask.sum()),
                "mean_elevation_m": float(np.mean(sub)) if len(sub) else float("nan"),
                "median_elevation_m": float(np.median(sub)) if len(sub) else float("nan"),
                "frac_above_sea": float(len(above) / len(sub)) if len(sub) else float("nan"),
            }

        # Land-volume decomposition: how much of total land_volume_above_sea_km3 traces to
        # is_volcano-tagged terrain (any node a volcanic field ever touched) vs everything
        # else, so a change in total volume can be attributed to volcanism specifically or to
        # something broader (e.g. #178's erosion-budget fix affecting non-volcanic terrain too).
        above_sea = elevation > world.sea_level_m
        volcano_land_m3 = float(np.sum((elevation[above_sea & is_volcano] - world.sea_level_m))) * area_m2
        other_land_m3 = float(np.sum((elevation[above_sea & ~is_volcano] - world.sea_level_m))) * area_m2

        records.append({
            "seed": seed,
            "checkpoint_years": checkpoint,
            "n_total_nodes": int(len(elevation)),
            "mean_all_elevation_m": float(np.mean(elevation)),
            "sea_level_m": float(world.sea_level_m),
            "is_volcano": _summary(is_volcano),
            "erupted_at_least_once": _summary(erupted),
            "currently_active": _summary(active),
            "volcano_land_volume_km3": volcano_land_m3 / 1.0e9,
            "other_land_volume_km3": other_land_m3 / 1.0e9,
            "total_land_volume_km3": (volcano_land_m3 + other_land_m3) / 1.0e9,
        })
    return records


def main() -> None:
    all_records = []
    for seed in SWEEP_SEEDS:
        print(f"=== seed={seed} ===", file=sys.stderr)
        all_records.extend(measure_one(seed))

    header = (
        f"{'seed':>10} {'ckpt_My':>8} "
        f"{'n_erupted':>9} {'mean_erupt_elev':>15} {'frac_above_sea':>15} "
        f"{'n_active':>8} {'mean_active_elev':>16} "
        f"{'n_is_volc':>9} {'mean_is_volc_elev':>17} {'mean_all_elev_m':>16}"
    )
    print(header)
    for r in all_records:
        e, a, v = r["erupted_at_least_once"], r["currently_active"], r["is_volcano"]
        print(
            f"{r['seed']:>10} {r['checkpoint_years'] // 1_000_000:>8} "
            f"{e['n']:>9} {e['mean_elevation_m']:>15.1f} {e['frac_above_sea']:>15.3f} "
            f"{a['n']:>8} {a['mean_elevation_m']:>16.1f} "
            f"{v['n']:>9} {v['mean_elevation_m']:>17.1f} {r['mean_all_elevation_m']:>16.1f}"
        )

    print("\n=== aggregate by checkpoint (mean across seeds) ===")
    for checkpoint in CHECKPOINTS:
        rows = [r for r in all_records if r["checkpoint_years"] == checkpoint]
        for key, label in (
            ("erupted_at_least_once", "erupted>=1x"),
            ("currently_active", "currently_active"),
            ("is_volcano", "is_volcano (all-time field provenance)"),
        ):
            mean_elev = np.nanmean([r[key]["mean_elevation_m"] for r in rows])
            mean_frac_above = np.nanmean([r[key]["frac_above_sea"] for r in rows])
            mean_n = np.mean([r[key]["n"] for r in rows])
            print(
                f"{checkpoint // 1_000_000:>4} My [{label:>38}]: mean_elevation_m={mean_elev:>9.1f} "
                f"frac_above_sea={mean_frac_above:.3f} mean_n={mean_n:.1f}"
            )

    print("\n=== land-volume decomposition (mean across seeds, km^3) ===")
    for checkpoint in CHECKPOINTS:
        rows = [r for r in all_records if r["checkpoint_years"] == checkpoint]
        mean_total = np.mean([r["total_land_volume_km3"] for r in rows])
        mean_volcano = np.mean([r["volcano_land_volume_km3"] for r in rows])
        mean_other = np.mean([r["other_land_volume_km3"] for r in rows])
        pct_volcano = 100.0 * mean_volcano / mean_total if mean_total else float("nan")
        print(
            f"{checkpoint // 1_000_000:>4} My: total={mean_total:.3e} "
            f"volcano_tagged={mean_volcano:.3e} ({pct_volcano:.1f}%) other={mean_other:.3e}"
        )

    import json
    out_path = Path(__file__).resolve().parent / "results" / "volcano_node_elevation.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        json.dump(all_records, f, indent=2)
    print(f"\nWrote {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
