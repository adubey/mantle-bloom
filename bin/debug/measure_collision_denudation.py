#!/usr/bin/env python3
"""Issue #275 phase 4: measure Hc-cap saturation and collision-belt denudation on real runs.

Generates a world per seed, steps it, and every few steps records:

- `hc_at_max_fraction`, the same node-count metric `/world/stats` reports, plus an
  area-weighted version;
- land fraction and every continental-material ledger reservoir;
- for columns at the Hc cap after the step: the share that is land, their median slope, and
  the share erosion took nothing from this step (a capped column erosion never touches stays
  capped);
- per-step erosion totals summed since the last snapshot: removed, collision-belt removal
  (columns at or above `--belt-fraction` of the cap), and `hc_cap_overflow_m3`.

Runs against any backend, so the same script gives a `main` baseline:

    cd backend
    .venv/bin/python ../bin/debug/measure_collision_denudation.py --seeds 1 2 3 --steps 60
    .venv/bin/python ../bin/debug/measure_collision_denudation.py --backend ../../mantle-bloom-origin-main/backend --seeds 1 2 3 --steps 60
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

LEDGER_KEYS = (
    "surface_continental_derived_m3",
    "continental_sediment_on_oceanic_hosts_m3",
    "delaminated_lower_crust_m3",
    "deeply_subducted_m3",
    "numerical_unplaced_m3",
    "discarded_marine_sediment_m3",
    "juvenile_additions_m3",
    "balance_error_m3",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10, help="snapshot every Nth step")
    parser.add_argument("--belt-fraction", type=float, default=0.9, help="Hc / cap that counts as a collision belt")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import continental_ledger, lithosphere
    from app.hydroclimate import erosion  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.plates import collect_all_accounting_areas_m2, collect_all_crustal_thickness  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    cap = lithosphere.MAX_CRUSTAL_THICKNESS_M
    belt_hc = args.belt_fraction * cap
    window: dict[str, float] = {}
    last_removed_m: dict[str, np.ndarray] = {}
    original_apply = erosion.apply_erosion

    def wrapped(world, years, node_cloud=None):
        # Erosion never changes node topology, so before/after Hc line up node for node.
        order = node_cloud[1] if node_cloud is not None else [p for p in world.plates if p.node_count() > 0]
        hc_before = collect_all_crustal_thickness(order)
        result = original_apply(world, years, node_cloud=node_cloud)
        if result is None:
            return result
        areas = collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
        removed_m = np.clip(hc_before - collect_all_crustal_thickness(order), 0.0, None)
        last_removed_m["m"] = removed_m
        budget = result.budget or {}
        window["removed_m3"] = window.get("removed_m3", 0.0) + budget.get("removed_m3", 0.0)
        window["belt_net_hc_loss_m3"] = window.get("belt_net_hc_loss_m3", 0.0) + float(np.sum((removed_m * areas)[hc_before >= belt_hc]))
        window["hc_cap_overflow_m3"] = window.get("hc_cap_overflow_m3", 0.0) + budget.get("hc_cap_overflow_m3", 0.0)
        for key, value in budget.items():
            if key.startswith("collision_") or key.startswith("mass_wasting_"):
                window[key] = window.get(key, 0.0) + value
        return result

    erosion.apply_erosion = wrapped

    rows = []
    for seed in args.seeds:
        world = generate_world(seed=seed)
        start = time.perf_counter()
        window.clear()
        for step in range(1, args.steps + 1):
            step_world(world, args.years)
            if step % args.every and step != args.steps:
                continue
            hydro = world.hydrology_cache
            order = hydro.plates_in_order
            hc = collect_all_crustal_thickness(order)
            areas = collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
            at_cap = hc >= cap - 1e-6
            land = ~hydro.is_ocean
            slope = world.erosion_cache.slope
            removed = last_removed_m.get("m", np.zeros_like(hc))
            inv = continental_ledger.inventories(world)
            row = {
                "seed": seed,
                "step": step,
                "seconds": round(time.perf_counter() - start, 1),
                "land_fraction": float(np.sum(areas[land]) / np.sum(areas)),
                "hc_at_max_fraction": float(np.mean(at_cap)),
                "hc_at_max_area_fraction": float(np.sum(areas[at_cap]) / np.sum(areas)),
                "belt_area_fraction": float(np.sum(areas[hc >= belt_hc]) / np.sum(areas)),
                "capped_land_share": float(np.mean(land[at_cap])) if at_cap.any() else 0.0,
                "capped_median_slope": float(np.median(slope[at_cap])) if at_cap.any() else 0.0,
                "capped_untouched_share": float(np.mean(removed[at_cap] <= 1e-6)) if at_cap.any() else 0.0,
                **{k: v for k, v in window.items()},
                **{f"ledger_{k}": float(inv[k]) for k in LEDGER_KEYS if k in inv},
            }
            window.clear()
            rows.append(row)
            print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
