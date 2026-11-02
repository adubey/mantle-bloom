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

MONTHS = 12
SUBSTEPS = 30
SPIN_YEARS = 40
LATENT_HEAT_MJ_KG = 2.45
W_M2_TO_MJ_M2_DAY = 0.0864
DAYS_PER_MONTH = 365.25 / MONTHS

# Wei et al. (2017) table 1, T/(E+T+I), by Köppen group. Arid and polar classes have no Wei class.
SPLIT_BY_GROUP = {"A": 0.61, "C": 0.55, "D": 0.55, "BS": 0.65, "ET": 0.65, "BW": 0.25, "EF": 0.0}


def transpiration_split(biome_ids: np.ndarray) -> np.ndarray:
    """The potential transpiration share of demand for each cell's Köppen class."""
    by_code = []
    for code in biomes.KOPPEN_CODES:
        group = code if code in SPLIT_BY_GROUP else code[:2] if code[:2] in SPLIT_BY_GROUP else code[0]
        by_code.append(SPLIT_BY_GROUP[group])
    table = np.array(by_code + [0.0] * len(biomes.PELAGIC_NAMES))
    return table[biome_ids]


def monthly_insolation_w_m2(lat_deg: np.ndarray, tilt_deg: float) -> np.ndarray:
    """Daily-mean top-of-atmosphere flux, (months, cells), by the formula of
    `climate.compute_insolation_w_m2`. Month 0 starts at the northern spring equinox."""
    lat = np.radians(lat_deg)[None, :]
    longitude = (np.arange(MONTHS) + 0.5) * 2.0 * np.pi / MONTHS
    declination = np.arcsin(np.sin(np.radians(tilt_deg)) * np.sin(longitude))[:, None]
    hour = np.arccos(np.clip(-np.tan(lat) * np.tan(declination), -1.0, 1.0))
    daily = hour * np.sin(lat) * np.sin(declination) + np.cos(lat) * np.cos(declination) * np.sin(hour)
    return np.clip(daily, 0.0, None) * climate.SOLAR_CONSTANT_W_M2 / np.pi


def monthly_forcing(fields, world, mask):
    """Monthly (temperature °C, rain mm, Oudin demand mm), each (months, cells), for the masked
    cells. Temperature peaks a month after the solstice; rain falls by half-year shares."""
    lat = np.broadcast_to(fields.lat_deg[:, None], mask.shape)[mask]
    temperature = fields.air_temperature_c[mask]
    rain = fields.precipitation_mm[mask]
    amplitude = biomes.seasonal_temp_amplitude(lat, biomes.grid_continentality(fields.is_ocean)[mask], world.axial_tilt_deg)
    north = lat >= 0.0
    warm_peak = np.where(north, 3.5, 9.5)
    phase = 2.0 * np.pi * ((np.arange(MONTHS) + 0.5)[:, None] - warm_peak) / MONTHS
    monthly_temperature = temperature + amplitude * np.cos(phase)
    strength = np.clip(world.axial_tilt_deg / biomes._SEASON_TILT_REFERENCE_DEG, 0.0, 1.3)
    summer_share, _ = biomes._precip_season(lat, strength)
    warm_half = np.cos(phase) >= 0.0
    monthly_rain = rain * np.where(warm_half, summer_share, 1.0 - summer_share) / 6.0
    insolation = monthly_insolation_w_m2(lat, world.axial_tilt_deg) * W_M2_TO_MJ_M2_DAY
    daily_demand = insolation / LATENT_HEAT_MJ_KG * np.clip(monthly_temperature + 5.0, 0.0, None) / 100.0
    return monthly_temperature, monthly_rain, daily_demand * DAYS_PER_MONTH


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
    _, rain, demand = monthly_forcing(fields, world, mask)
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
    split = transpiration_split(fields.biome_ids[mask])
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
