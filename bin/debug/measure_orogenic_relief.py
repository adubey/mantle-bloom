#!/usr/bin/env python3
"""Issue #290: where collision crust goes, and whether continents hold their area and volume.

Generates a world per seed, steps it, and every few steps records:

- land fraction, continental-host area, and live continental-derived material;
- cumulative `delaminated_lower_crust_m3`, as a share of the starting inventory;
- Hc-cap saturation (area fraction at the cap), and the continental Hc mean / p99;
- continental Hc volume against the first snapshot's (runaway thickening shows here even
  where the material tracer doesn't, e.g. conserved ceiling overflow);
- with phase 2's fields: restite volume, and the area-weighted mean Moho lag on continent;
- with this branch's `orogeny.py`: the area in each crust state and every
  `world.orogenic_relief_budget` account (suture belts, escape, delamination, far field,
  no-outlet remainder, collapse, ductile flow).

Runs against any backend, so the same script gives a `main` baseline:

    cd backend
    .venv/bin/python ../bin/debug/measure_orogenic_relief.py --seeds 1 2 3 --steps 100
    .venv/bin/python ../bin/debug/measure_orogenic_relief.py --backend ../../mantle-bloom-origin-main/backend --seeds 1 2 3 --steps 100
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from seed_run_lib import import_backend, run_seeds, seed_run_parser


def main() -> None:
    parser = seed_run_parser(__doc__)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    import_backend(args)
    from app import continental_ledger, lithosphere  # noqa: E402
    from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
    from app.surface_fields import SURFACE_FIELDS  # noqa: E402

    anatexis_fields = "restite_m" in SURFACE_FIELDS

    try:
        from app import orogeny  # noqa: E402
    except ImportError:
        orogeny = None

    cap = lithosphere.MAX_CRUSTAL_THICKNESS_M
    initial_by_seed: dict[int, float] = {}
    first_hc_volume: dict[int, float] = {}

    def start_seed(world, seed) -> None:
        continental_ledger.ensure_initialized(world)
        initial_by_seed[seed] = world.continental_material_ledger["initial_continental_m3"]

    def measure(world, seed, step) -> dict:
        initial = initial_by_seed[seed]
        spacing = line_spacing_rad(world.node_density)
        hc_parts, area_parts, cont_parts, restite_parts, lag_parts = [], [], [], [], []
        for plate in world.plates:
            if plate.node_count() == 0:
                continue
            hc_parts.append(plate.collect("crustal_thickness_m"))
            area_parts.append(plate.accounting_areas_m2(spacing))
            cont_parts.append(
                effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
            )
            zeros = np.zeros(plate.node_count())
            restite_parts.append(plate.collect("restite_m") if anatexis_fields else zeros)
            lag_parts.append(plate.collect("moho_thermal_lag_c") if anatexis_fields else zeros)
        hc, areas, continental, restite, lag = (
            np.concatenate(p) for p in (hc_parts, area_parts, cont_parts, restite_parts, lag_parts)
        )
        total_area = float(areas.sum())
        hydro = world.hydrology_cache
        land_area = float(np.sum(areas_in_order(hydro, spacing)[~hydro.is_ocean])) if hydro is not None else float("nan")
        inv = continental_ledger.inventories(world)
        ledger = world.continental_material_ledger
        cont_hc = hc[continental]
        cont_hc_volume = float(cont_hc @ areas[continental])
        first_hc_volume.setdefault(seed, cont_hc_volume)
        row = {
            "land_fraction": land_area / total_area,
            "continental_area_fraction": float(areas[continental].sum()) / total_area,
            "continental_material_m3": inv["surface_continental_derived_m3"],
            "continental_material_vs_initial": inv["surface_continental_derived_m3"] / initial,
            "delaminated_m3": ledger["delaminated_lower_crust_m3"],
            "delaminated_vs_initial": ledger["delaminated_lower_crust_m3"] / initial,
            "deeply_subducted_vs_initial": ledger["deeply_subducted_m3"] / initial,
            "numerical_unplaced_vs_initial": ledger["numerical_unplaced_m3"] / initial,
            "hc_at_max_area_fraction": float(areas[hc >= cap - 1e-6].sum()) / total_area,
            "continental_hc_mean_m": float(np.average(cont_hc, weights=areas[continental])) if cont_hc.size else 0.0,
            "continental_hc_p99_m": float(np.percentile(cont_hc, 99)) if cont_hc.size else 0.0,
            "continental_hc_volume_vs_first": cont_hc_volume / first_hc_volume[seed],
            "restite_vs_initial": float(restite @ areas) / initial,
            "continental_mean_moho_lag_c": float(np.average(lag[continental], weights=areas[continental])) if cont_hc.size else 0.0,
        }
        if orogeny is not None:
            states = orogeny.crust_state_areas_m2(world)
            row.update({f"state_{name}_fraction": value / total_area for name, value in states.items()})
            row.update({f"budget_{key}": value for key, value in orogeny.ensure_budget(world).items()})
        return row

    run_seeds(args, measure, start_seed)

def areas_in_order(hydro, spacing: float) -> np.ndarray:
    """Accounting areas in the hydrology cache's node order, which `is_ocean` follows."""
    return np.concatenate([p.accounting_areas_m2(spacing) for p in hydro.plates_in_order])


if __name__ == "__main__":
    main()
