#!/usr/bin/env python3
"""Issue #275 phase 3: measure the isostatic ice-load deflection on real runs.

Generates a world per seed, steps it, and every few steps records:

- land fraction and total ice volume (km^3);
- the stored `ice_load_deflection_m`: deepest node, area-weighted mean over glaciated land,
  and the fraction of land carrying more than 10 m of depression;
- `drift_m`: the largest gap between the stored deflection and the one the current ice and
  column imply. Tectonics runs after erosion each step, so this isn't zero, but it should
  stay small next to the deflection itself.

Deflection columns are absent on a backend without the field, so the same script also gives
a `main` baseline for land fraction and ice volume:

    cd backend
    .venv/bin/python ../bin/debug/measure_ice_loading.py --seeds 1 2 3 --steps 60
    .venv/bin/python ../bin/debug/measure_ice_loading.py --backend ../../mantle-bloom-origin-main/backend --seeds 1 2 3 --steps 60
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def snapshot(world, modules) -> dict[str, float]:
    hydrology, lithosphere, plates_mod, line_spacing_rad = modules
    hydro = world.hydrology_cache
    if hydro is None:
        return {}
    order = hydro.plates_in_order
    land = ~hydro.is_ocean
    ice = hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
    areas = plates_mod.collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
    land_area = max(float(np.sum(areas[land])), 1.0)
    out = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "land_fraction": float(np.sum(areas[land]) / np.sum(areas)),
        "ice_fraction_of_land": float(np.sum(areas[land & ice]) / land_area),
        "ice_volume_km3": float(np.dot(areas, hydro.glacier_depth)) / 1e9,
    }
    if "ice_load_deflection_m" not in plates_mod.SURFACE_FIELDS:
        return out
    # The world's current nodes (tectonics may have changed them since hydrology ran).
    current = list(world.plates)
    w = np.concatenate([p.collect("ice_load_deflection_m") for p in current])
    glacier = np.concatenate([p.collect("glacier_depth") for p in current])
    elevation = np.concatenate([p.collect("elevation") for p in current])
    hc = np.concatenate([p.collect("crustal_thickness_m") for p in current])
    hm = np.concatenate([p.collect("mantle_lithosphere_thickness_m") for p in current])
    rho = np.concatenate([lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in current])
    cur_areas = plates_mod.collect_all_accounting_areas_m2(current, line_spacing_rad(world.node_density))
    cur_land = elevation > world.sea_level_m
    load = lithosphere.grounded_ice_load_kg_m2(glacier, elevation - w, world.sea_level_m)
    implied = lithosphere.ice_load_deflection(hc, hm, rho, load)
    glaciated_land = cur_land & (glacier >= hydrology.GLACIER_VISIBLE_DEPTH_M)
    out.update(
        {
            "deflection_min_m": float(w.min()),
            "deflection_mean_glaciated_land_m": float(np.dot(cur_areas[glaciated_land], w[glaciated_land]) / max(np.sum(cur_areas[glaciated_land]), 1.0)),
            "land_depressed_over_10m_fraction": float(np.sum(cur_areas[cur_land & (w < -10.0)]) / max(np.sum(cur_areas[cur_land]), 1.0)),
            "drift_m": float(np.abs(w - implied).max()),
        }
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10, help="print every Nth step")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import lithosphere, plates
    from app.hydroclimate import hydrology  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    modules = (hydrology, lithosphere, plates, line_spacing_rad)
    rows = []
    for seed in args.seeds:
        world = generate_world(seed=seed)
        start = time.perf_counter()
        for step in range(1, args.steps + 1):
            step_world(world, args.years)
            if step % args.every == 0 or step == args.steps:
                row = {"seed": seed, "step": step, "seconds": round(time.perf_counter() - start, 1), **snapshot(world, modules)}
                rows.append(row)
                print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()}), flush=True)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
