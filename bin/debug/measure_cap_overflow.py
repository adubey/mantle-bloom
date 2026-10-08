#!/usr/bin/env python3
"""Issue #288: measure depositional Hc-cap overflow and where it ends up on real runs.

Generates a world per seed, steps it, and every `--every` steps records:

- Hc-cap saturation (`hc_at_max_area_fraction`) and crustal-thickening indicators: the
  area-weighted mean Hc of continental hosts, and the area share of columns at or above
  `--belt-fraction` of the cap -- a runaway would show as both climbing without bound;
- erosion's overflow diagnostics summed since the last snapshot: cap-hit node count,
  attempted, placed-onward and terminal volume, total and per depositional pathway, plus the
  rock the overflow-laden ice scoured;
- every continental-material loss account, kept apart so depositional overflow
  (`numerical_unplaced_m3` / `overloaded_root_delaminated_m3`) can be told from
  suture-accretion delamination (`delaminated_lower_crust_m3`) and divergent thinning
  (`rift_thinned_m3`), and `numerical_unplaced_m3` as a share of the initial inventory.

Runs against any backend (missing keys are reported as 0), so the same script gives a `main`
baseline:

    cd backend
    .venv/bin/python ../bin/debug/measure_cap_overflow.py --seeds 1 2 3 --steps 100
    .venv/bin/python ../bin/debug/measure_cap_overflow.py --backend ../../mantle-bloom-origin-main/backend --seeds 1 2 3 --steps 100
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

LEDGER_KEYS = (
    "initial_continental_m3",
    "surface_continental_derived_m3",
    "numerical_unplaced_m3",
    "overloaded_root_delaminated_m3",
    "delaminated_lower_crust_m3",
    "rift_thinned_m3",
    "deeply_subducted_m3",
    "discarded_marine_sediment_m3",
    "balance_error_m3",
)
PATHWAYS = ("river", "lake", "lake_silt", "glacial", "aeolian", "mass_wasting", "coastal_leveling", "marine")
TOTALS = (
    "hc_cap_overflow_nodes",
    "hc_cap_overflow_m3",
    "hc_cap_overflow_placed_m3",
    "hc_cap_overflow_terminal_m3",
    "hc_cap_overflow_ice_scour_m3",
    "continental_overflow_m3",
    "continental_overflow_terminal_m3",
    "continental_unplaced_m3",
    "hc_clip_residual_m3",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10, help="snapshot every Nth step")
    parser.add_argument("--belt-fraction", type=float, default=0.9, help="Hc / cap that counts as thickened")
    parser.add_argument("--per-pathway", action="store_true", help="also print per-pathway columns")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import continental_ledger, erosion, lithosphere  # noqa: E402
    from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
    from app.plates import collect_all_accounting_areas_m2, collect_all_crustal_thickness  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    cap = lithosphere.MAX_CRUSTAL_THICKNESS_M
    window: dict[str, float] = {}
    original_apply = erosion.apply_erosion

    def wrapped(world, years, node_cloud=None):
        result = original_apply(world, years, node_cloud=node_cloud)
        budget = (result.budget if result is not None else None) or {}
        for key in TOTALS:
            window[key] = window.get(key, 0.0) + budget.get(key, 0.0)
        for pathway in PATHWAYS:
            for suffix in ("nodes", "m3", "placed_m3", "terminal_m3"):
                key = f"hc_cap_overflow_{pathway}_{suffix}"
                window[key] = window.get(key, 0.0) + budget.get(key, 0.0)
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
            order = world.hydrology_cache.plates_in_order
            hc = collect_all_crustal_thickness(order)
            areas = collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
            continental = np.concatenate(
                [effective_is_continental_from_codes(p.collect("crust_type_code"), p.crust_type == "continental") for p in order]
            )
            inv = continental_ledger.inventories(world)
            initial = max(float(inv.get("initial_continental_m3", 0.0)), 1.0)
            row = {
                "seed": seed,
                "step": step,
                "myr": step * args.years / 1e6,
                "seconds": round(time.perf_counter() - start, 1),
                "hc_at_max_area_fraction": float(np.sum(areas[hc >= cap - 1e-6]) / np.sum(areas)),
                "thickened_area_fraction": float(np.sum(areas[hc >= args.belt_fraction * cap]) / np.sum(areas)),
                "continental_mean_hc_km": float(np.sum(hc[continental] * areas[continental]) / max(np.sum(areas[continental]), 1.0)) / 1e3,
                "unplaced_share_of_initial": float(inv.get("numerical_unplaced_m3", 0.0)) / initial,
                **{k: v for k, v in window.items() if args.per_pathway or not any(f"_{p}_" in k for p in PATHWAYS)},
                **{f"ledger_{k}": float(inv.get(k, 0.0)) for k in LEDGER_KEYS},
            }
            window.clear()
            rows.append(row)
            print(json.dumps({k: (float(f"{v:.4g}") if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
            if args.out:
                args.out.parent.mkdir(parents=True, exist_ok=True)
                args.out.write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
