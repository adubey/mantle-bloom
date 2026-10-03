#!/usr/bin/env python3
"""Issue #275 phase 2: measure glacier extent, to compare seasonal melt against annual-mean melt.

Generates a world per seed, steps it, and after every step records, from the world's own
hydrology/climate caches:

- the fraction of land nodes carrying visible ice (>= GLACIER_VISIBLE_DEPTH_M);
- total ice volume (km^3, ice depth * accounting area);
- how that glaciated land splits by annual-mean temperature band, so ice held only because
  the annual mean is below freezing (but summer is warm) shows up in the -15..0 C bands.

`--backend` points at the backend to import, so the same script measures `main` and a branch:

    cd backend
    .venv/bin/python ../bin/debug/measure_seasonal_melt.py --seeds 1 2 --steps 30
    .venv/bin/python ../bin/debug/measure_seasonal_melt.py --backend ~/projects/mantle-bloom/backend --seeds 1 2 --steps 30
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

BANDS_C = (-np.inf, -25.0, -15.0, -10.0, -5.0, 0.0, np.inf)


def snapshot(world, modules) -> dict[str, float]:
    erosion, hydrology, plates_mod, line_spacing_rad = modules
    hydro = world.hydrology_cache
    fields = world.climate_cache
    if hydro is None or fields is None:
        return {}
    land = ~hydro.is_ocean
    ice = hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
    areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
    row, col = erosion.climate_grid_indices(hydro.points, *fields.air_temperature_c.shape)
    temp = fields.air_temperature_c[row, col]
    out = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "land_fraction": float(np.sum(areas[land]) / np.sum(areas)),
        "ice_fraction_of_land": float(np.sum(areas[land & ice]) / max(np.sum(areas[land]), 1.0)),
        "ice_volume_km3": float(np.dot(areas, hydro.glacier_depth)) / 1e9,
    }
    for lo, hi in zip(BANDS_C[:-1], BANDS_C[1:]):
        band = land & ice & (temp >= lo) & (temp < hi)
        out[f"ice_land_frac_T[{lo:g},{hi:g})"] = float(np.sum(areas[band]) / max(np.sum(areas[land]), 1.0))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", default=[1])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--surface", default="quad")
    parser.add_argument("--every", type=int, default=5, help="print every Nth step")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    from app import erosion, hydrology, plates  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    modules = (erosion, hydrology, plates, line_spacing_rad)
    rows = []
    for seed in args.seeds:
        world = generate_world(seed=seed, surface=args.surface)
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
