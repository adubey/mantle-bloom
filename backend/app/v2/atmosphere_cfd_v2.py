"""HEALPix port of `atmosphere_cfd.py`'s shallow-water wind solver -- same equations, same
constants, same `init_*`/`refresh_forcing`/`step_*` structure, `(npix,)` flat arrays instead
of `(H, W)`, `fluid_dynamics_healpix` primitives instead of `fluid_dynamics`, and no polar
zonal filter/sponge drag (nothing to compensate for on an equal-area grid). See that module's
own docstring for the physical reasoning behind every term -- unchanged here, only the grid
substrate moved.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .. import climate, fluid_dynamics
from . import fluid_dynamics_healpix as fdh
from . import healpix_grid

if TYPE_CHECKING:
    from .world_v2 import WorldV2

REDUCED_GRAVITY_M_S2 = 0.5
EFFECTIVE_TROPOSPHERE_DEPTH_M = 8000.0
ETA_TEMPERATURE_COUPLING_M_PER_C = 8.0
ETA_THERMAL_RELAXATION_PER_S = 3.0e-5

SURFACE_DRAG_BASE_PER_S = 1.0e-5
SURFACE_DRAG_OROGRAPHIC_PER_S = 8.0e-6
VISCOSITY_M2_S = 6.0e4

MOUNTAIN_OBSTACLE_ELEVATION_M = 2000.0
MOUNTAIN_OBSTACLE_RAMP_M = 800.0
MOUNTAIN_DEFLECTION_RATE_PER_S = 3.0e-4

TEMPERATURE_DIFFUSIVITY_M2_S = 3.0e4
RADIATIVE_RELAXATION_PER_S = 3.0e-6

HUMIDITY_DIFFUSIVITY_M2_S = 2.0e4
OCEAN_EVAPORATION_SOURCE_PER_S = 4.0e-7
LAND_EVAPORATION_SOURCE_PER_S = 0.6e-7
CONDENSATION_RATE_PER_S = 2.0e-6
PRECIP_CONDENSATION_TO_MM = 4.0e6

MAX_SUBSTEPS_PER_STEP = 2000
SECONDS_PER_TECTONIC_STEP = 86400.0  # one simulated day, same as v1


@dataclass
class AtmosphereCFDStateV2:
    grid: healpix_grid.HealpixGrid
    is_ocean: np.ndarray  # (npix,) bool -- refreshed once per tectonics step
    elevation_m: np.ndarray  # (npix,) -- refreshed once per tectonics step
    u: np.ndarray  # (npix,) eastward wind, m/s
    v: np.ndarray  # (npix,) northward wind, m/s
    eta: np.ndarray  # (npix,) geopotential-height anomaly, m
    temperature_c: np.ndarray
    equilibrium_temperature_c: np.ndarray  # refreshed once per tectonics step
    humidity: np.ndarray
    precipitation_mm: np.ndarray
    elapsed_seconds: float = 0.0

    def resample_uv_to_equirect(self, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
        """The seam `climate.py`'s `compute_climate` calls polymorphically against either a
        v1 (equirectangular) or v2 (HEALPix) CFD state -- see climate.py's own small edit."""
        u = healpix_grid.resample_to_equirect(self.grid, self.u, height, width)
        v = healpix_grid.resample_to_equirect(self.grid, self.v, height, width)
        return u, v


def _equilibrium_temperature(world: "WorldV2", grid: healpix_grid.HealpixGrid, is_ocean: np.ndarray, elevation_m: np.ndarray) -> np.ndarray:
    """`climate.compute_land_temperature`/`compute_ocean_temperature_baseline` both assume
    the classic "(H,) insolation broadcasting against (H, W) cells" shape -- called here with
    `height=npix, width=1` (one pixel per "row", each its own independent insolation value)
    rather than trying to force HEALPix's genuinely 1D `(npix,)` layout into a `(1, npix)`
    single shared-insolation row, which would silently broadcast wrong (every pixel getting
    every *other* pixel's insolation too, an (npix, npix) result). Squeezed back to `(npix,)`
    after."""
    lat_deg = np.degrees(grid.lat_rad)
    insolation = climate.compute_insolation(lat_deg, world.axial_tilt_deg, world.solar_multiplier)
    land = climate.compute_land_temperature(insolation, elevation_m[:, None])[:, 0]
    ocean = climate.compute_ocean_temperature_baseline(insolation, grid.npix, 1)[:, 0]
    return np.where(is_ocean, ocean, land)


def init_atmosphere_cfd(world: "WorldV2", terrain: climate.ClimateFields) -> AtmosphereCFDStateV2:
    """Bootstraps this world's permanent `WorldV2.atmosphere_cfd_state` from `terrain` (v1's
    own `climate.compute_climate` diagnostic snapshot, on the equirectangular render/erosion
    grid) resampled onto a fresh HEALPix grid -- the HEALPix analogue of v1's own
    `init_atmosphere_cfd`, which bootstraps from the same diagnostic pipeline."""
    grid = healpix_grid.build(healpix_grid.nside_for_density(world.fluid_density))
    lat_rows = terrain.lat_deg
    is_ocean = healpix_grid.resample_from_equirect(grid, terrain.is_ocean.astype(np.float64), lat_rows) > 0.5
    elevation_m = healpix_grid.resample_from_equirect(grid, terrain.elevation_m, lat_rows)
    u0 = healpix_grid.resample_from_equirect(grid, terrain.wind_u, lat_rows)
    v0 = healpix_grid.resample_from_equirect(grid, terrain.wind_v, lat_rows)
    temperature0 = healpix_grid.resample_from_equirect(grid, np.where(terrain.is_ocean, terrain.ocean_temperature_c, terrain.air_temperature_c), lat_rows)
    humidity0 = healpix_grid.resample_from_equirect(grid, terrain.humidity, lat_rows)
    precip0 = healpix_grid.resample_from_equirect(grid, terrain.precipitation_mm, lat_rows)
    equilibrium_temperature_c = _equilibrium_temperature(world, grid, is_ocean, elevation_m)

    return AtmosphereCFDStateV2(
        grid=grid,
        is_ocean=is_ocean,
        elevation_m=elevation_m.astype(np.float32),
        u=u0.astype(np.float32),
        v=v0.astype(np.float32),
        eta=np.zeros(grid.npix, dtype=np.float32),
        temperature_c=temperature0.astype(np.float32),
        equilibrium_temperature_c=equilibrium_temperature_c.astype(np.float32),
        humidity=humidity0.astype(np.float32),
        precipitation_mm=precip0.astype(np.float32),
    )


def refresh_forcing(world: "WorldV2", state: AtmosphereCFDStateV2, terrain: climate.ClimateFields) -> None:
    lat_rows = terrain.lat_deg
    state.is_ocean = healpix_grid.resample_from_equirect(state.grid, terrain.is_ocean.astype(np.float64), lat_rows) > 0.5
    state.elevation_m = healpix_grid.resample_from_equirect(state.grid, terrain.elevation_m, lat_rows).astype(np.float32)
    state.equilibrium_temperature_c = _equilibrium_temperature(world, state.grid, state.is_ocean, state.elevation_m).astype(np.float32)


def _mountain_deflection_geometry(elevation_m: np.ndarray, grid: healpix_grid.HealpixGrid):
    gx, gy = fdh.gradient(elevation_m, grid)
    magnitude = np.hypot(gx, gy)
    uphill_x = np.divide(gx, magnitude, out=np.zeros_like(gx), where=magnitude > 1e-12)
    uphill_y = np.divide(gy, magnitude, out=np.zeros_like(gy), where=magnitude > 1e-12)
    normal_x, normal_y = -uphill_x, -uphill_y
    tangent_x, tangent_y = -normal_y, normal_x
    ramp = np.clip((elevation_m - (MOUNTAIN_OBSTACLE_ELEVATION_M - MOUNTAIN_OBSTACLE_RAMP_M)) / (2.0 * MOUNTAIN_OBSTACLE_RAMP_M), 0.0, 1.0)
    return normal_x, normal_y, tangent_x, tangent_y, ramp


def _mountain_deflection_tendency(u: np.ndarray, v: np.ndarray, geom) -> tuple[np.ndarray, np.ndarray]:
    normal_x, normal_y, tangent_x, tangent_y, ramp = geom
    into_slope = np.clip(-(u * normal_x + v * normal_y), 0.0, None)
    rate = MOUNTAIN_DEFLECTION_RATE_PER_S * ramp * into_slope
    ax = rate * (normal_x + tangent_x)
    ay = rate * (normal_y + tangent_y)
    return ax, ay


def step_atmosphere_cfd(world: "WorldV2", state: AtmosphereCFDStateV2, seconds: float) -> None:
    del world
    grid = state.grid
    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * EFFECTIVE_TROPOSPHERE_DEPTH_M))
    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, fdh.min_spacing_m(grid), wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)

    orographic_drag = SURFACE_DRAG_BASE_PER_S + SURFACE_DRAG_OROGRAPHIC_PER_S * np.clip(state.elevation_m / MOUNTAIN_OBSTACLE_ELEVATION_M, 0.0, 1.0)
    mountain_geom = _mountain_deflection_geometry(state.elevation_m, grid)
    f = fluid_dynamics.coriolis_parameter(np.degrees(grid.lat_rad))

    u, v, eta = state.u, state.v, state.eta
    temperature_c = state.temperature_c
    humidity = state.humidity
    precipitation_mm = state.precipitation_mm

    for _ in range(n_substeps):
        deta_dx, deta_dy = fdh.gradient(eta, grid)

        deflect_ax, deflect_ay = _mountain_deflection_tendency(u, v, mountain_geom)
        du_dt = f * v - REDUCED_GRAVITY_M_S2 * deta_dx + VISCOSITY_M2_S * fdh.laplacian(u, grid) + deflect_ax
        dv_dt = -f * u - REDUCED_GRAVITY_M_S2 * deta_dy + VISCOSITY_M2_S * fdh.laplacian(v, grid) + deflect_ay
        u = (u + dt_s * du_dt) / (1.0 + dt_s * orographic_drag)
        v = (v + dt_s * dv_dt) / (1.0 + dt_s * orographic_drag)
        u = fdh.grid_noise_filter(u, grid)
        v = fdh.grid_noise_filter(v, grid)

        eta_target = ETA_TEMPERATURE_COUPLING_M_PER_C * (temperature_c - float(temperature_c.mean()))
        flux_divergence = fdh.divergence(EFFECTIVE_TROPOSPHERE_DEPTH_M * u, EFFECTIVE_TROPOSPHERE_DEPTH_M * v, grid)
        eta = eta - dt_s * flux_divergence + dt_s * ETA_THERMAL_RELAXATION_PER_S * (eta_target - eta)
        eta = fdh.grid_noise_filter(eta, grid)

        temperature_c = fdh.semi_lagrangian_advect(temperature_c, u, v, dt_s, grid)
        temperature_c = temperature_c + dt_s * TEMPERATURE_DIFFUSIVITY_M2_S * fdh.laplacian(temperature_c, grid)
        temperature_c = temperature_c + dt_s * RADIATIVE_RELAXATION_PER_S * (state.equilibrium_temperature_c - temperature_c)

        evap_source = np.where(state.is_ocean, np.float32(OCEAN_EVAPORATION_SOURCE_PER_S), np.float32(LAND_EVAPORATION_SOURCE_PER_S))
        saturation_ceiling = np.clip(temperature_c / climate.EVAPORATION_REFERENCE_TEMP_C, climate.MIN_EVAPORATION_CEILING, climate.MAX_EVAPORATION_CEILING)
        excess = np.clip(humidity - saturation_ceiling, 0.0, None)
        condensed = CONDENSATION_RATE_PER_S * excess

        humidity = fdh.semi_lagrangian_advect(humidity, u, v, dt_s, grid)
        humidity = humidity + dt_s * HUMIDITY_DIFFUSIVITY_M2_S * fdh.laplacian(humidity, grid)
        humidity = np.clip(humidity + dt_s * (evap_source - condensed), 0.0, None)
        precipitation_mm = condensed * PRECIP_CONDENSATION_TO_MM

        state.elapsed_seconds += dt_s

    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.humidity = humidity
    state.precipitation_mm = precipitation_mm
