"""HEALPix port of `ocean_cfd.py`'s shallow-water ocean-current solver -- same equations,
same constants, same structure as `atmosphere_cfd_v2.py`'s own port of `atmosphere_cfd.py`.
See that module and `ocean_cfd.py` for the physical reasoning; unchanged here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .. import climate, fluid_dynamics
from . import atmosphere_cfd_v2
from . import fluid_dynamics_healpix as fdh
from . import healpix_grid

if TYPE_CHECKING:
    from .world_v2 import WorldV2

REDUCED_GRAVITY_M_S2 = 0.03
MIN_DEPTH_M = 200.0
MAX_DEPTH_M = 6000.0

WIND_STRESS_COEFFICIENT = 1.5e-6
MIXED_LAYER_DEPTH_M = 50.0
BOTTOM_DRAG_PER_S = 2.0e-5
VISCOSITY_M2_S = 4.0e4
TEMPERATURE_DIFFUSIVITY_M2_S = 2.0e4
TEMPERATURE_RELAXATION_PER_S = 1.0e-6

SEDIMENT_DIFFUSIVITY_M2_S = 1.5e4
SEDIMENT_PICKUP_COEFFICIENT = 2.0e-4
SEDIMENT_SETTLING_SPEED_THRESHOLD_M_S = 0.08
SEDIMENT_SETTLING_RATE_PER_S = 3.0e-6
SEDIMENT_DEPOSIT_DEPTH_COEFFICIENT = 0.5

MAX_SUBSTEPS_PER_STEP = 2000
SECONDS_PER_TECTONIC_STEP = 7 * 86400.0  # one simulated week, same as v1


@dataclass
class OceanCFDStateV2:
    grid: healpix_grid.HealpixGrid
    is_ocean: np.ndarray
    elevation_m: np.ndarray
    depth_m: np.ndarray
    u: np.ndarray
    v: np.ndarray
    eta: np.ndarray
    temperature_c: np.ndarray
    baseline_temperature_c: np.ndarray
    sediment_concentration: np.ndarray
    sediment_deposited_m: np.ndarray
    wind_u: np.ndarray
    wind_v: np.ndarray
    elapsed_seconds: float = 0.0

    def resample_uv_to_equirect(self, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
        u = healpix_grid.resample_to_equirect(self.grid, self.u, height, width)
        v = healpix_grid.resample_to_equirect(self.grid, self.v, height, width)
        return u, v


def init_ocean_cfd(world: "WorldV2", terrain: climate.ClimateFields, atmosphere_state: atmosphere_cfd_v2.AtmosphereCFDStateV2) -> OceanCFDStateV2:
    """Same HEALPix grid `atmosphere_cfd_v2.init_atmosphere_cfd` already built for this
    world (`atmosphere_state.grid`) -- reused directly rather than building a second,
    independent one, so `wind_u`/`wind_v` need no resample (both states share one grid,
    matching v1's own same-`fluid_density`-grid reasoning)."""
    grid = atmosphere_state.grid
    lat_rows = terrain.lat_deg
    is_ocean = healpix_grid.resample_from_equirect(grid, terrain.is_ocean.astype(np.float64), lat_rows) > 0.5
    elevation_m = healpix_grid.resample_from_equirect(grid, terrain.elevation_m, lat_rows)
    depth_m = np.where(is_ocean, np.clip(-elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0).astype(np.float32)
    ocean_temp0 = healpix_grid.resample_from_equirect(grid, terrain.ocean_temperature_c, lat_rows).astype(np.float32)
    zeros = np.zeros(grid.npix, dtype=np.float32)

    return OceanCFDStateV2(
        grid=grid,
        is_ocean=is_ocean,
        elevation_m=elevation_m,
        depth_m=depth_m,
        u=zeros.copy(),
        v=zeros.copy(),
        eta=zeros.copy(),
        temperature_c=ocean_temp0.copy(),
        baseline_temperature_c=ocean_temp0.copy(),
        sediment_concentration=zeros.copy(),
        sediment_deposited_m=zeros.copy(),
        wind_u=atmosphere_state.u.copy(),
        wind_v=atmosphere_state.v.copy(),
    )


def refresh_forcing(world: "WorldV2", state: OceanCFDStateV2, terrain: climate.ClimateFields) -> None:
    lat_rows = terrain.lat_deg
    state.is_ocean = healpix_grid.resample_from_equirect(state.grid, terrain.is_ocean.astype(np.float64), lat_rows) > 0.5
    state.elevation_m = healpix_grid.resample_from_equirect(state.grid, terrain.elevation_m, lat_rows)
    state.depth_m = np.where(state.is_ocean, np.clip(-state.elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0).astype(np.float32)
    state.wind_u = world.atmosphere_cfd_state.u.copy()
    state.wind_v = world.atmosphere_cfd_state.v.copy()


def step_ocean_cfd(world: "WorldV2", state: OceanCFDStateV2, seconds: float) -> None:
    del world
    grid = state.grid
    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * max(float(state.depth_m.max()), MIN_DEPTH_M)))
    coastal = fdh.coastal_ocean_mask(state.is_ocean, grid)

    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, fdh.min_spacing_m(grid), wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)
    drag = BOTTOM_DRAG_PER_S
    f = fluid_dynamics.coriolis_parameter(np.degrees(grid.lat_rad))

    u, v, eta = state.u, state.v, state.eta
    temperature_c = state.temperature_c
    sediment_concentration = state.sediment_concentration
    sediment_deposited_m = state.sediment_deposited_m

    for _ in range(n_substeps):
        deta_dx, deta_dy = fdh.gradient(eta, grid)

        wind_speed = np.hypot(state.wind_u, state.wind_v)
        tau_u = (WIND_STRESS_COEFFICIENT * wind_speed * state.wind_u) / MIXED_LAYER_DEPTH_M
        tau_v = (WIND_STRESS_COEFFICIENT * wind_speed * state.wind_v) / MIXED_LAYER_DEPTH_M

        du_dt = f * v - REDUCED_GRAVITY_M_S2 * deta_dx + tau_u + VISCOSITY_M2_S * fdh.laplacian(u, grid)
        dv_dt = -f * u - REDUCED_GRAVITY_M_S2 * deta_dy + tau_v + VISCOSITY_M2_S * fdh.laplacian(v, grid)

        u = (u + dt_s * du_dt) / (1.0 + dt_s * drag)
        v = (v + dt_s * dv_dt) / (1.0 + dt_s * drag)
        u = np.where(state.is_ocean, fdh.grid_noise_filter(u, grid), 0.0)
        v = np.where(state.is_ocean, fdh.grid_noise_filter(v, grid), 0.0)

        flux_divergence = fdh.divergence(state.depth_m * u, state.depth_m * v, grid)
        eta = eta - dt_s * flux_divergence
        eta = np.where(state.is_ocean, fdh.grid_noise_filter(eta, grid), 0.0)

        temperature_c = fdh.semi_lagrangian_advect(temperature_c, u, v, dt_s, grid)
        temperature_c = temperature_c + dt_s * TEMPERATURE_DIFFUSIVITY_M2_S * fdh.laplacian(temperature_c, grid)
        temperature_c = temperature_c + dt_s * TEMPERATURE_RELAXATION_PER_S * (state.baseline_temperature_c - temperature_c)

        speed = np.hypot(u, v)
        pickup = np.where(coastal, SEDIMENT_PICKUP_COEFFICIENT * speed, 0.0)
        settling_factor = np.clip(1.0 - speed / SEDIMENT_SETTLING_SPEED_THRESHOLD_M_S, 0.0, 1.0)
        settle = SEDIMENT_SETTLING_RATE_PER_S * sediment_concentration * settling_factor

        sediment_concentration = fdh.semi_lagrangian_advect(sediment_concentration, u, v, dt_s, grid)
        sediment_concentration = sediment_concentration + dt_s * SEDIMENT_DIFFUSIVITY_M2_S * fdh.laplacian(sediment_concentration, grid)
        sediment_concentration = np.clip(sediment_concentration + dt_s * (pickup - settle), 0.0, None)
        sediment_concentration = np.where(state.is_ocean, sediment_concentration, 0.0)

        sediment_deposited_m = sediment_deposited_m + dt_s * settle * SEDIMENT_DEPOSIT_DEPTH_COEFFICIENT

        state.elapsed_seconds += dt_s

    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.sediment_concentration = sediment_concentration
    state.sediment_deposited_m = sediment_deposited_m
