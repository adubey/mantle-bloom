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

from pathlib import Path

import numpy as np
from seed_run_lib import import_backend, run_seeds, seed_run_parser

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
    parser = seed_run_parser(__doc__, every=5, every_help="print every Nth step")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    import_backend(args)
    from app import plates
    from app.hydroclimate import erosion, hydrology  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402

    modules = (erosion, hydrology, plates, line_spacing_rad)
    run_seeds(args, lambda world, seed, step: snapshot(world, modules))


if __name__ == "__main__":
    main()
