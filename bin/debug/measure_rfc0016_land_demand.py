#!/usr/bin/env python3
"""RFC 0016 G1a: scale-check the preregistered land demand on a checkpoint's saved climate.

Reads `world.climate_cache` and never steps the world, so nothing here is a paired run. On the
climate grid's unfrozen land it evaluates Oudin demand from monthly temperature and insolation,
runs the seasonal bucket to its repeating year, and compares that year's water partition and
transpiration with the control's Budyko partition and rainfall-proxy humidity source. Results
are land-area weighted. They check scale; they don't select a parameter.

Example:

    backend/.venv/bin/python bin/debug/measure_rfc0016_land_demand.py <checkpoint.mbworld> \\
        --out <scratch>/rfc0016-g1a-window-d.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

BACKEND = Path(__file__).resolve().parents[2] / "backend"
sys.path.insert(0, str(BACKEND))

from app import persistence  # noqa: E402
from app.hydroclimate import biomes, breaching, climate  # noqa: E402
from rfc0016_land_store import MONTHS, monthly_forcing, transpiration_split_table  # noqa: E402

SUBSTEPS = 30
SPIN_YEARS = 40


def forcing(fields, world, mask):
    """Monthly rain and Oudin demand on the masked cells, each (months, cells)."""
    lat = np.broadcast_to(fields.lat_deg[:, None], mask.shape)[mask]
    amplitude = biomes.seasonal_temp_amplitude(
        lat, biomes.grid_continentality(fields.is_ocean)[mask], world.axial_tilt_deg
    )
    return monthly_forcing(
        lat, fields.air_temperature_c[mask], amplitude, fields.precipitation_mm[mask],
        world.axial_tilt_deg, climate.SOLAR_CONSTANT_W_M2, biomes,
    )


def repeating_year(rain, demand, capacity_mm, depletion_p):
    """Actual ET and saturation excess, in mm/yr, for the bucket's repeating year. Sub-steps stand
    in for the closed-form piecewise update here; this measures scale, not closure."""
    store = np.full(rain.shape[1], 0.5 * capacity_mm)
    critical = max((1.0 - depletion_p) * capacity_mm, 1e-9)
    et_year = excess_year = None
    for _ in range(SPIN_YEARS):
        et_year = np.zeros_like(store)
        excess_year = np.zeros_like(store)
        for month in range(MONTHS):
            for _ in range(SUBSTEPS):
                available = store + rain[month] / SUBSTEPS
                stress = np.clip(available / critical, 0.0, 1.0)
                et = np.minimum(demand[month] / SUBSTEPS * stress, available)
                store = available - et
                excess = np.clip(store - capacity_mm, 0.0, None)
                store -= excess
                et_year += et
                excess_year += excess
    return et_year, excess_year


def measure(fields, world, mask, weights, demand_factor, capacity_mm, depletion_p, split):
    rain, demand = forcing(fields, world, mask)
    et, excess = repeating_year(rain, demand_factor * demand, capacity_mm, depletion_p)
    precipitation = fields.precipitation_mm[mask]
    control_runoff = breaching.runoff_mm(
        precipitation, breaching.open_water_evaporation_mm(fields.air_temperature_c[mask])
    )
    weighted = lambda values: float(np.average(values, weights=weights))  # noqa: E731
    transpiration = et * split
    control_source_mm = (
        climate.VEGETATION_TRANSPIRATION_BY_BIOME[fields.biome_ids[mask]]
        * climate.VEGETATION_RECYCLING_FRACTION
        * precipitation
    )
    return {
        "demand_mm": weighted(demand_factor * demand.sum(axis=0)),
        "treatment_et_mm": weighted(et),
        "treatment_excess_mm": weighted(excess),
        "control_runoff_mm": weighted(control_runoff),
        "excess_change_percent": 100.0 * (weighted(excess) / weighted(control_runoff) - 1.0),
        "treatment_transpiration_mm": weighted(transpiration),
        "control_humidity_source_mm": weighted(control_source_mm),
        "source_ratio": weighted(transpiration) / max(weighted(control_source_mm), 1e-9),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    world = persistence.load_world_bytes(args.save.read_bytes())
    fields = world.climate_cache
    mask = ~fields.is_ocean & (fields.air_temperature_c >= 0.0)
    weights = np.broadcast_to(np.cos(np.radians(fields.lat_deg))[:, None], mask.shape)[mask]
    split = transpiration_split_table(biomes)[fields.biome_ids[mask]]
    precipitation = fields.precipitation_mm[mask]
    result = {
        "save": args.save.name,
        "elapsed_years": world.elapsed_years,
        "land_precipitation_mm": float(np.average(precipitation, weights=weights)),
        "control_open_water_demand_mm": float(
            np.average(breaching.open_water_evaporation_mm(fields.air_temperature_c[mask]), weights=weights)
        ),
        "arms": {},
    }
    arms = {
        "central": (1.0, 300.0, 0.5),
        "capacity-150": (1.0, 150.0, 0.5),
        "capacity-1000": (1.0, 1000.0, 0.5),
        "demand-0.85": (0.85, 300.0, 0.5),
        "demand-1.15": (1.15, 300.0, 0.5),
        "depletion-0.3": (1.0, 300.0, 0.3),
        "depletion-0.7": (1.0, 300.0, 0.7),
    }
    for name, (factor, capacity, depletion) in arms.items():
        result["arms"][name] = measure(fields, world, mask, weights, factor, capacity, depletion, split)
        print(name, json.dumps(result["arms"][name]), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
