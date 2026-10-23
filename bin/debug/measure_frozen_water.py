#!/usr/bin/env python3
"""RFC 0015: measure how a saved world represents frozen water, without stepping it.

Reads one `.mbworld` save and reports, as JSON:

- **Land ice on the map** (the Biome/Combined views' fine grid, area-weighted by cos(lat)):
  land under visible ice (>= GLACIER_VISIBLE_DEPTH_M), the share of it the Biome view
  classifies as Köppen EF (Ice Cap) or ET (Tundra), the EF land that carries no ice, the
  mean ice thickness on visibly iced land, and the share of iced land whose warmest month is
  above 0 C at the bed but would fall below it at the ice surface (LAPSE_RATE_C_PER_KM per km
  of ice).
- **Snowfall on frozen land** (same grid): hydrology treats all precipitation on land whose
  annual mean is below FREEZE_POINT_C as snow. Under the model's own sinusoidal annual
  cycle (`biomes.seasonal_temp_amplitude`), the precipitation-weighted share of the year
  above 0 C on that land is the share it counts as snow but that falls in thawing months.
- **Polar ocean** (same grid, poleward of 60°): the ocean temperature's floor
  (`climate.WATER_TEMP_MIN_C`) and how much polar ocean falls below seawater's freezing point
  in the annual mean and in the warmest month, against the air temperature over the same
  ocean, to show whether either can mark perennial sea ice.
- **Sea-ice balance by latitude** (RFC 0015 D6, this world's tilt and sunlight): the
  perennial, seasonal or ice-free state from freezing and thawing degree-days of the
  land-temperature construct at sea level, for two ice degree-day factors.
- **Ice and sea level** (plate nodes, exact areas): the ice volume `eustasy` debits from the
  ocean (all `glacier_depth`, x 0.917), split into floating ice, grounded marine ice and land
  ice, against the water equivalent of the ice's volume above flotation.

Usage:

    cd backend
    .venv/bin/python ../bin/debug/measure_frozen_water.py <save>.mbworld --out <scratch>/frozen.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from replay_lib import load_replay

POLAR_LAT_DEG = 60.0
SEAWATER_FREEZING_C = -1.8
SEA_ICE_LATITUDES_DEG = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 89.0)
SEA_ICE_OCEAN_FLUX_W_M2 = 2.0
SEA_ICE_DDF_M_PER_DAY = (0.004, 0.007)


def area_weights(lat_deg: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    return np.broadcast_to(np.cos(np.radians(lat_deg)).reshape(-1, 1), shape)


def share(mask: np.ndarray, within: np.ndarray, weights: np.ndarray) -> float:
    total = float(weights[within].sum())
    return float(weights[within & mask].sum()) / total if total > 0 else 0.0


def polar_ocean_measures(world, lat_grid, is_ocean, air_temp, ocean_temp, w) -> dict:
    """Whether the ocean temperature can show perennial sea ice, against the air above it,
    on ocean poleward of POLAR_LAT_DEG. The warmest and coldest months use the maritime
    amplitude (`seasonal_temp_amplitude` at continentality 0)."""
    from app.hydroclimate import biomes
    from app.hydroclimate.climate import WATER_TEMP_MIN_C

    polar = is_ocean & (np.abs(lat_grid) >= POLAR_LAT_DEG)
    amp = biomes.seasonal_temp_amplitude(lat_grid, np.zeros_like(lat_grid), world.axial_tilt_deg)
    t_f = SEAWATER_FREEZING_C
    out = {
        "polar_ocean_share_of_ocean": share(polar, is_ocean, w),
        "sst_min_c": float(ocean_temp[polar].min()) if polar.any() else None,
        "sst_at_floor": share(ocean_temp <= WATER_TEMP_MIN_C + 0.01, polar, w),
        "sst_annual_below_freezing": share(ocean_temp <= t_f, polar, w),
        "sst_warmest_month_below_freezing": share(ocean_temp + amp <= t_f, polar, w),
        "air_finite": share(np.isfinite(air_temp), polar, w),
    }
    if polar.any() and np.isfinite(air_temp[polar]).all():
        out.update(
            air_mean_c=float(np.average(air_temp[polar], weights=w[polar])),
            air_annual_below_freezing=share(air_temp <= t_f, polar, w),
            air_warmest_month_below_freezing=share(air_temp + amp <= t_f, polar, w),
            air_coldest_month_below_freezing=share(air_temp - amp <= t_f, polar, w),
        )
    return out


def sea_ice_balance_by_latitude(world) -> list[dict]:
    """RFC 0015 D6's sea-ice state by latitude at sea level, from the climate's own
    constructs for this world's tilt and sunlight: the land-temperature construct at sea
    level as the surface over ice, the maritime amplitude, Stefan growth from freezing
    degree-days, and loss from DDF_ice x thawing degree-days plus an ocean heat flux of
    SEA_ICE_OCEAN_FLUX_W_M2. Reported for two ice degree-day factors."""
    from app.hydroclimate import biomes, climate, hydrology
    from app.hydroclimate.ice_sheet import ICE_LATENT_HEAT_J_KG, ICE_THERMAL_CONDUCTIVITY_W_M_K
    from app.lithosphere import RHO_ICE, RHO_WATER

    lat = np.array(SEA_ICE_LATITUDES_DEG)
    insolation = climate.compute_heating_index(lat, world.axial_tilt_deg, world.solar_multiplier)
    t_ice = climate.LAND_TEMP_MIN_C + climate.LAND_TEMP_RANGE_C * insolation
    sst = climate.WATER_TEMP_MIN_C + climate.WATER_TEMP_RANGE_C * insolation
    amp = biomes.seasonal_temp_amplitude(lat, np.zeros_like(lat), world.axial_tilt_deg)
    fdd = 365.25 * hydrology.seasonal_positive_degrees(-t_ice, amp, -SEAWATER_FREEZING_C)
    tdd = 365.25 * hydrology.seasonal_positive_degrees(t_ice, amp, 0.0)
    growth = np.sqrt(2.0 * ICE_THERMAL_CONDUCTIVITY_W_M_K * fdd * 86_400.0 / (RHO_ICE * ICE_LATENT_HEAT_J_KG))
    bottom = SEA_ICE_OCEAN_FLUX_W_M2 * 3.15576e7 / (RHO_ICE * ICE_LATENT_HEAT_J_KG)
    rows = []
    for i, latitude in enumerate(lat):
        row = {"lat_deg": float(latitude), "sst_baseline_c": round(float(sst[i]), 2), "t_ice_c": round(float(t_ice[i]), 2),
               "amplitude_c": round(float(amp[i]), 2), "fdd": round(float(fdd[i])), "tdd": round(float(tdd[i])),
               "growth_m": round(float(growth[i]), 2)}
        for ddf in SEA_ICE_DDF_M_PER_DAY:
            loss = ddf * tdd[i] * RHO_WATER / RHO_ICE + bottom
            state = "ice-free" if fdd[i] <= 0 else ("perennial" if growth[i] > loss else "seasonal")
            row[f"state_ddf_{ddf * 1000:g}mm"] = state
        rows.append(row)
    return rows


def map_measures(world) -> dict:
    """Land ice and snowfall on the Biome view's own grid and classification."""
    from app import render_image
    from app.hydroclimate import biomes, hydrology
    from app.hydroclimate.climate import LAPSE_RATE_C_PER_KM

    fields, winner_ids, _share, _runner = render_image._classify_biome_grid(world)
    lat_deg, _lon, _xyz, elevation, is_ocean, air_temp, ocean_temp, precip, lake_depth, glacier_depth, _is_sea, _hs = fields
    shape = elevation.shape
    w = area_weights(lat_deg, shape)
    biome = winner_ids.reshape(shape)
    land = ~is_ocean & (lake_depth <= hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
    iced = land & (glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M)
    lat_grid = np.broadcast_to(lat_deg.reshape(-1, 1), shape)
    amplitude = biomes.seasonal_temp_amplitude(lat_grid, biomes.grid_continentality(is_ocean), world.axial_tilt_deg)
    warm_bed = air_temp + amplitude
    warm_ice = warm_bed - LAPSE_RATE_C_PER_KM * glacier_depth / 1000.0
    frozen = land & (air_temp < hydrology.FREEZE_POINT_C)
    d = np.clip(air_temp / np.maximum(amplitude, 1e-9), -1.0, 1.0)
    thaw_share = np.arccos(-d) / np.pi  # share of the year above 0 C under the sinusoid
    pw = w * precip
    return {
        "polar_ocean": polar_ocean_measures(world, lat_grid, is_ocean, air_temp, ocean_temp, w),
        "land_area_share_of_globe": float(w[land].sum() / w.sum()),
        "iced_share_of_land": share(iced, land, w),
        "iced_land_classified_EF": share(biome == biomes.ICE_CAP, iced, w),
        "iced_land_classified_ET": share(biome == biomes.koppen_index("ET"), iced, w),
        "EF_land_without_ice": share(~iced, land & (biome == biomes.ICE_CAP), w),
        "EF_share_of_land": share(biome == biomes.ICE_CAP, land, w),
        "iced_land_mean_thickness_m": float(np.average(glacier_depth[iced], weights=w[iced])) if iced.any() else 0.0,
        "iced_land_warm_month_above_0_at_bed_below_0_at_surface": share((warm_bed >= 0.0) & (warm_ice < 0.0), iced, w),
        "frozen_share_of_land": share(frozen, land, w),
        "frozen_land_precip_in_thawing_months": float((pw[frozen] * thaw_share[frozen]).sum() / max(pw[frozen].sum(), 1e-9)),
    }


def sea_level_measures(world) -> dict:
    """The ice eustasy debits from the ocean, by where it sits, against its volume above flotation."""
    from app import plates as plates_mod
    from app.elevation_lines import line_spacing_rad
    from app.hydroclimate import ice_sheet
    from app.lithosphere import RHO_ICE, RHO_WATER

    area = plates_mod.collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
    depth = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
    bed = plates_mod.collect_all_elevation(world.plates) - plates_mod.collect_all_ice_load_deflection(world.plates)
    sea = float(world.sea_level_m)
    water_depth = np.clip(sea - bed, 0.0, None)
    grounded = ice_sheet.is_grounded(depth, bed, sea)
    marine = water_depth > 0.0
    above_flotation = np.where(grounded, depth - ice_sheet.flotation_thickness_m(bed, sea), 0.0)
    ratio = RHO_ICE / RHO_WATER

    def km3_we(mask: np.ndarray, column: np.ndarray = depth) -> float:
        return ratio * float(np.dot(np.where(mask, column, 0.0), area)) / 1e9

    return {
        "sea_level_m": sea,
        "ice_debited_km3_we": km3_we(np.ones_like(grounded)),
        "floating_ice_km3_we": km3_we(~grounded & (depth > 0.0)),
        "grounded_marine_ice_km3_we": km3_we(grounded & marine),
        "land_ice_km3_we": km3_we(grounded & ~marine),
        "above_flotation_km3_we": km3_we(grounded, np.clip(above_flotation, 0.0, None)),
        "ocean_area_km2": float(area[marine & ~grounded].sum()) / 1e6,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.step_years = None  # load_replay's header field; this probe doesn't step
    world, result = load_replay(args, issue=368)
    result.update(
        save=args.save.name, elapsed_myr=world.elapsed_years / 1e6,
        map=map_measures(world), sea_level=sea_level_measures(world),
        sea_ice_balance=sea_ice_balance_by_latitude(world),
    )
    del result["checkpoints"], result["step_years"], result["backend"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
