#!/usr/bin/env python3
"""RFC 0015 gate G5: sea ice from the air above it, on saved worlds and at three tilts (#548).

D6 gives each ocean cell a derived sea-ice state from a degree-day balance:
- the surface over ice, `T_si`, is the climate's land-temperature construct at sea level
  (ice-age cooling included) plus the cell's ocean-current anomaly
  (`ocean_temperature_c` - the zonal ocean baseline), with the maritime half-amplitude
  (`biomes.seasonal_temp_amplitude` at continentality 0);
- freezing degree-days below SEAWATER_FREEZING_C and thawing degree-days above 0 C, in closed
  form from that sinusoid;
- annual growth by Stefan's law, sqrt(2 k FDD 86400 / (rho_i L_f)); annual loss
  DDF_ice x TDD x rho_w / rho_i plus bottom melt F_w / (rho_i L_f) per year;
- *perennial* where FDD > 0 and growth exceeds loss, *seasonal* where FDD > 0 otherwise, and
  the ice-covered share of the year `f_ice` = 1 if perennial, else the share below T_f.

For each save, without stepping it, the climate is computed twice: as saved with ice ages
off, and at the glacial maximum (`ice_age_glacial_intensity` = 1, by setting the period to
twice the elapsed time). Ice-age cooling is uniform in both baselines (`ice_age_cooling_c`),
so the saved state at full intensity is the ice-age arm; no replay is needed. Each arm reports,
for every (F_w, DDF_ice) setting, `sea_ice_extent_stats.hemisphere_summary` with and without
the current anomaly, how far the anomaly moves the per-longitude perennial edge, and the Polar
Sea Ice class area today (`biomes.classify_pelagic`'s SST cutoff) against f_ice >= 0.5.

`--zonal` adds the pure-function bands at tilts of 15, 23.5 and 35 degrees on a 0.25 degree
grid: the latitudes where ice becomes perennial and where freezing starts.

Usage (follow the compute-budget skill before each save; peak RSS about 2.5 GB):

    cd backend
    .venv/bin/python ../bin/debug/gate_sea_ice.py <save>.mbworld --out <scratch>/g5-<seed>.json
    .venv/bin/python ../bin/debug/gate_sea_ice.py --zonal-only --out <scratch>/g5-zonal.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from sea_ice_extent_stats import cell_areas_km2, hemisphere_summary, per_longitude_edge

SEAWATER_FREEZING_C = -1.8
SECONDS_PER_YEAR = 3.15576e7
OCEAN_HEAT_FLUX_W_M2 = (1.3, 2.0, 2.7)  # Maykut & Untersteiner (1971)
# Hock (2003) Table 1: multi-season fits on bare ice span 5.5-8.2 mm w.e. d-1 C-1, and
# Braithwaite (1995) finds about 8 in West Greenland; 20.0 is Table 1's largest single-season
# value (GIMEX, 1028 m), kept as an extreme.
DDF_ICE_M_PER_DAY = (0.0055, 0.008, 0.020)
CENTRAL = (2.0, 0.008)
ANOMALY_MOVE_DEG = 2.0
ZONAL_TILTS_DEG = (15.0, 23.5, 35.0)
ZONAL_STEP_DEG = 0.25


def sea_ice_state(t_si: np.ndarray, amplitude: np.ndarray, ddf_m_per_day: float, ocean_flux_w_m2: float):
    """D6's (perennial, seasonal, f_ice) for surface temperature `t_si` and half-amplitude
    `amplitude`, both arrays of the same shape."""
    from app.hydroclimate import hydrology
    from app.hydroclimate.ice_sheet import ICE_LATENT_HEAT_J_KG, ICE_THERMAL_CONDUCTIVITY_W_M_K
    from app.lithosphere import RHO_ICE, RHO_WATER

    fdd = 365.25 * hydrology.seasonal_positive_degrees(-t_si, amplitude, -SEAWATER_FREEZING_C)
    tdd = 365.25 * hydrology.seasonal_positive_degrees(t_si, amplitude, 0.0)
    growth = np.sqrt(2.0 * ICE_THERMAL_CONDUCTIVITY_W_M_K * fdd * 86_400.0 / (RHO_ICE * ICE_LATENT_HEAT_J_KG))
    loss = ddf_m_per_day * tdd * RHO_WATER / RHO_ICE + ocean_flux_w_m2 * SECONDS_PER_YEAR / (RHO_ICE * ICE_LATENT_HEAT_J_KG)
    freezes = fdd > 0.0
    perennial = freezes & (growth > loss)
    seasonal = freezes & ~perennial
    below = 0.5 + np.arcsin(np.clip((SEAWATER_FREEZING_C - t_si) / np.maximum(amplitude, 1e-9), -1.0, 1.0)) / np.pi
    f_ice = np.where(perennial, 1.0, np.where(freezes, below, 0.0))
    return perennial, seasonal, f_ice


def _setting_name(flux: float, ddf: float) -> str:
    return f"Fw{flux:g}_DDF{ddf * 1000:g}"


def _anomaly_move(lat: np.ndarray, ocean: np.ndarray, with_anomaly: np.ndarray, zonal: np.ndarray) -> dict:
    out = {}
    for name, north in (("north", True), ("south", False)):
        a = per_longitude_edge(lat, ocean, with_anomaly, north)
        z = per_longitude_edge(lat, ocean, zonal, north)
        # Only columns with perennial ice both ways: a column that loses or gains all its ice
        # reads as an edge of 90 degrees on one side, which isn't a move of the edge.
        both = (a < 90.0) & (z < 90.0)
        moved = np.abs(a - z)[both]
        out[name] = None if moved.size == 0 else {
            "max_deg": round(float(moved.max()), 2),
            "columns_moved_ge_2deg": int((moved >= ANOMALY_MOVE_DEG).sum()),
            "columns": int(moved.size),
            "columns_with_ice_one_way_only": int(((a < 90.0) ^ (z < 90.0)).sum()),
        }
    return out


def measure_arm(world) -> dict:
    """Every G5 measure for the world's current climate."""
    from app.hydroclimate import biomes, climate

    fields = climate.compute_climate(world)
    lat = fields.lat_deg
    ocean = fields.is_ocean
    height, width = ocean.shape
    cooling = climate.ice_age_cooling_c(world)
    insolation = climate.compute_heating_index(lat, world.axial_tilt_deg, world.solar_multiplier)
    baseline = climate.compute_ocean_temperature_baseline(insolation, height, width) - cooling
    anomaly = np.where(ocean, fields.ocean_temperature_c - baseline, 0.0)
    # On ocean cells the land construct has no lapse (elevation <= 0), so it is the
    # construct at sea level, ice-age cooling included.
    t_zonal = fields.land_temperature_c
    t_si = t_zonal + anomaly
    lat_grid = np.broadcast_to(lat[:, None], ocean.shape)
    amplitude = biomes.seasonal_temp_amplitude(lat_grid, np.zeros_like(lat_grid), world.axial_tilt_deg)
    area = cell_areas_km2(lat, width)

    settings = {}
    for flux in OCEAN_HEAT_FLUX_W_M2:
        for ddf in DDF_ICE_M_PER_DAY:
            perennial, seasonal, f_ice = sea_ice_state(t_si, amplitude, ddf, flux)
            z_perennial, z_seasonal, _ = sea_ice_state(t_zonal, amplitude, ddf, flux)
            perennial, seasonal, z_perennial, z_seasonal = (m & ocean for m in (perennial, seasonal, z_perennial, z_seasonal))
            with_anomaly = hemisphere_summary(lat, area, ocean, perennial, seasonal)
            without_anomaly = hemisphere_summary(lat, area, ocean, z_perennial, z_seasonal)
            if (flux, ddf) != CENTRAL:  # keep the JSON compact: bands only at the central setting
                for summary in (*with_anomaly.values(), *without_anomaly.values()):
                    del summary["bands"]
            settings[_setting_name(flux, ddf)] = {
                "with_anomaly": with_anomaly,
                "without_anomaly": without_anomaly,
                "anomaly_moves_perennial_edge": _anomaly_move(lat, ocean, perennial, z_perennial),
                "polar_sea_ice_class_after_mkm2": {
                    name: round(float(area[ocean & (f_ice >= 0.5) & side[:, None]].sum()) / 1e6, 3)
                    for name, side in (("north", lat > 0), ("south", lat < 0))
                },
            }
    anomaly_polar = ocean & (np.abs(lat_grid) >= 50.0)
    return {
        "ice_age_cooling_c": round(cooling, 3),
        "polar_sea_ice_class_before_mkm2": round(float(area[ocean & (fields.biome_ids == biomes.POLAR_SEA_ICE)].sum()) / 1e6, 3),
        "anomaly_poleward_50_c": {
            "min": round(float(anomaly[anomaly_polar].min()), 2) if anomaly_polar.any() else None,
            "max": round(float(anomaly[anomaly_polar].max()), 2) if anomaly_polar.any() else None,
        },
        "settings": settings,
    }


def measure_save(world) -> dict:
    saved_period = getattr(world, "ice_age_period_years", 0.0)
    arms = {}
    try:
        world.ice_age_period_years = 0.0
        arms["no_ice_age"] = measure_arm(world)
        world.ice_age_period_years = 2.0 * world.elapsed_years  # half a period: full intensity
        arms["ice_age_max"] = measure_arm(world)
    finally:
        world.ice_age_period_years = saved_period
    return {
        "seed": int(world.seed),
        "elapsed_myr": round(world.elapsed_years / 1e6, 3),
        "node_density": world.node_density,
        "axial_tilt_deg": float(world.axial_tilt_deg),
        "solar_multiplier": float(world.solar_multiplier),
        "wind_model": getattr(world, "wind_model", None),
        "arms": arms,
    }


def zonal_bands() -> dict:
    """Perennial and freezing limits by tilt, without the current anomaly or ice ages."""
    from app.hydroclimate import biomes, climate

    lat = np.arange(ZONAL_STEP_DEG / 2.0, 90.0, ZONAL_STEP_DEG)
    out = {}
    for tilt in ZONAL_TILTS_DEG:
        t_si = climate.LAND_TEMP_MIN_C + climate.LAND_TEMP_RANGE_C * climate.compute_heating_index(lat, tilt)
        amplitude = biomes.seasonal_temp_amplitude(lat, np.zeros_like(lat), tilt)
        rows = {}
        for flux in OCEAN_HEAT_FLUX_W_M2:
            for ddf in DDF_ICE_M_PER_DAY:
                perennial, seasonal, _ = sea_ice_state(t_si, amplitude, ddf, flux)
                rows[_setting_name(flux, ddf)] = {
                    "perennial_from_deg": float(lat[perennial].min()) if perennial.any() else None,
                    "freezing_from_deg": float(lat[perennial | seasonal].min()) if (perennial | seasonal).any() else None,
                }
        out[f"tilt_{tilt:g}"] = rows
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path, nargs="?")
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--zonal-only", action="store_true", help="only the pure-function tilt bands; no save")
    args = ap.parse_args()
    if args.zonal_only:
        sys.path.insert(0, str(args.backend.resolve()))
        result = {"zonal_bands": zonal_bands()}
    else:
        from replay_lib import load_replay

        args.step_years = None  # load_replay's header field; this gate doesn't step
        world, result = load_replay(args, issue=548)
        result.update(save=args.save.name, **measure_save(world))
        del result["checkpoints"], result["step_years"], result["backend"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result)[:2000])


if __name__ == "__main__":
    main()
