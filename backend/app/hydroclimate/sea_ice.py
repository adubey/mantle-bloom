"""RFC 0015 D6's floating sea ice, derived for display from the annual climate."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..lithosphere import RHO_ICE, RHO_WATER
from . import biomes, climate, hydrology, ice_sheet

if TYPE_CHECKING:
    from ..world import World

SEAWATER_FREEZING_C = -1.8
OCEAN_HEAT_FLUX_W_M2 = 2.0
POLAR_SEA_ICE_FRACTION = 0.5


@dataclass(frozen=True)
class SeaIce:
    temperature_c: np.ndarray
    perennial: np.ndarray
    seasonal: np.ndarray
    fraction: np.ndarray


def annual_balance(temperature_c: np.ndarray, amplitude_c: np.ndarray, is_ocean: np.ndarray) -> SeaIce:
    """Stefan growth from winter degree-days versus surface and basal loss (G5, #548)."""
    freezing_days = 365.25 * hydrology.seasonal_freezing_degrees(temperature_c, amplitude_c, SEAWATER_FREEZING_C)
    thawing_days = 365.25 * hydrology.seasonal_positive_degrees(temperature_c, amplitude_c, 0.0)
    growth_m = np.sqrt(
        2.0 * ice_sheet.ICE_THERMAL_CONDUCTIVITY_W_M_K * freezing_days * 86_400.0
        / (RHO_ICE * ice_sheet.ICE_LATENT_HEAT_J_KG)
    )
    surface_loss_m = ice_sheet.DDF_ICE_MM_WATER_PER_DAY_K * 1e-3 * thawing_days * RHO_WATER / RHO_ICE
    bottom_loss_m = OCEAN_HEAT_FLUX_W_M2 * ice_sheet.W_M2_TO_M_ICE_PER_YR
    freezes = is_ocean & (freezing_days > 0.0)
    perennial = freezes & (growth_m > surface_loss_m + bottom_loss_m)
    seasonal = freezes & ~perennial
    winter_fraction = hydrology.seasonal_fraction_above(-temperature_c, amplitude_c, -SEAWATER_FREEZING_C)
    fraction = np.where(perennial, 1.0, np.where(seasonal, winter_fraction, 0.0))
    return SeaIce(temperature_c, perennial, seasonal, fraction)


def climate_grid_balance(
    world: World, latitude: np.ndarray, ocean_temperature_c: np.ndarray, is_ocean: np.ndarray,
) -> SeaIce:
    """Diagnose D6 sea ice from the finished climate fields, for display and reporting."""
    heating = climate.compute_heating_index(latitude, world.axial_tilt_deg, world.solar_multiplier)
    ocean_baseline_c = climate.compute_ocean_temperature_baseline(heating, len(latitude), 1)
    sea_level_air_c = climate.compute_land_temperature(heating, np.zeros((len(latitude), 1)))
    # Climate cools both baselines equally. Their difference stays fixed, so the sampled
    # SST supplies both the current anomaly and the cooling of this climate snapshot.
    temperature_c = ocean_temperature_c + (sea_level_air_c - ocean_baseline_c)
    amplitude_c = biomes.seasonal_temp_amplitude(latitude, np.zeros_like(latitude), world.axial_tilt_deg)
    return annual_balance(temperature_c, amplitude_c[:, None], is_ocean)


def sample(
    world: World, fields: climate.ClimateFields, row: np.ndarray, col: np.ndarray, is_ocean: np.ndarray,
) -> SeaIce:
    """Read the proxy at the corresponding climate cells and gate it by the display ocean mask.

    Where a finer map resolves ocean within a coarse land cell, the climate grid has zero
    current anomaly there. The sea-level temperature and maritime seasonality still apply.
    """
    balance = climate_grid_balance(world, fields.lat_deg, fields.ocean_temperature_c, np.ones_like(fields.is_ocean))
    return SeaIce(
        balance.temperature_c[row, col],
        is_ocean & balance.perennial[row, col],
        is_ocean & balance.seasonal[row, col],
        np.where(is_ocean, balance.fraction[row, col], 0.0),
    )
