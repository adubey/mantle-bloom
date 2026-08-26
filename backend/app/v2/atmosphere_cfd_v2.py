"""HEALPix port of `atmosphere_cfd.py`'s shallow-water wind solver -- same equations, same
constants, same `init_*`/`refresh_forcing`/`step_*` structure, `(npix,)` flat arrays instead
of `(H, W)`, `fluid_dynamics_healpix` primitives instead of `fluid_dynamics`, and no polar
zonal filter/sponge drag (nothing to compensate for on an equal-area grid). See that module's
own docstring for the physical reasoning behind every term -- unchanged here, only the grid
substrate moved.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from .. import climate, fluid_dynamics
from . import fluid_dynamics_healpix as fdh
from . import healpix_grid
from .healpix_grid import ang2pix_nest_scalar

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


@njit(cache=True, fastmath=True)
def _atmosphere_substep_loop_kernel(
    u,
    v,
    eta,
    temperature_c,
    humidity,
    precipitation_mm,
    is_ocean,
    equilibrium_temperature_c,
    orographic_drag,
    normal_x,
    normal_y,
    tangent_x,
    tangent_y,
    ramp,
    f,
    neighbours,
    neighbour_valid,
    neighbour_dx_m,
    neighbour_dy_m,
    gradient_inv,
    world_xyz,
    east,
    north,
    radius_m,
    nside,
    order,
    dt_s,
    n_substeps,
    reduced_gravity,
    trop_depth,
    eta_temp_coupling,
    eta_thermal_relaxation,
    viscosity,
    mountain_deflection_rate,
    temp_diffusivity,
    radiative_relaxation,
    ocean_evap_source,
    land_evap_source,
    evap_reference_temp,
    min_evap_ceiling,
    max_evap_ceiling,
    condensation_rate,
    humidity_diffusivity,
    precip_conversion,
):
    """`step_atmosphere_cfd`'s entire per-substep body (and the `n_substeps` loop itself)
    fused into one compiled function -- the atmosphere-solver twin of
    `ocean_cfd_v2._ocean_substep_loop_kernel`; see that function's own docstring for why this
    exists and what it fuses. Structural differences from the ocean kernel: no `is_ocean`
    masking on `u`/`v`/`eta` (wind exists over land too, unlike ocean currents), a per-pixel
    `orographic_drag` array instead of ocean's scalar bottom drag, the mountain-deflection
    tendency (purely local given `u`/`v` and the once-per-step `mountain_geom` arrays) riding
    along in the same pass as the viscosity/Coriolis momentum update, one extra whole-grid
    reduction (`temperature_c`'s own mean, needed for `eta_target`) folded into the
    first-derivative gather pass since that's already a full traversal, and humidity/
    evaporation/condensation/precipitation in place of ocean's sediment pickup/settling."""
    npix = u.shape[0]
    n_neighbours = neighbours.shape[1]

    gx_eta = np.empty(npix, dtype=np.float32)
    gy_eta = np.empty(npix, dtype=np.float32)
    gx_u = np.empty(npix, dtype=np.float32)
    gy_u = np.empty(npix, dtype=np.float32)
    gx_v = np.empty(npix, dtype=np.float32)
    gy_v = np.empty(npix, dtype=np.float32)
    u_raw = np.empty(npix, dtype=np.float32)
    v_raw = np.empty(npix, dtype=np.float32)
    filt_u = np.empty(npix, dtype=np.float32)
    filt_v = np.empty(npix, dtype=np.float32)
    flux_u = np.empty(npix, dtype=np.float32)
    flux_v = np.empty(npix, dtype=np.float32)
    eta_raw = np.empty(npix, dtype=np.float32)
    advected_temp = np.empty(npix, dtype=np.float32)
    advected_hum = np.empty(npix, dtype=np.float32)
    gx_temp = np.empty(npix, dtype=np.float32)
    gy_temp = np.empty(npix, dtype=np.float32)
    gx_hum = np.empty(npix, dtype=np.float32)
    gy_hum = np.empty(npix, dtype=np.float32)

    for _substep in range(n_substeps):
        # Pass A: first-derivative gather of eta, u, v together, plus temperature_c's own
        # mean (needed for eta_target below) -- folded into this pass since it's already a
        # full traversal over npix.
        temp_sum = 0.0
        for i in range(npix):
            ei = eta[i]
            ui = u[i]
            vi = v[i]
            atb_x_eta = 0.0
            atb_y_eta = 0.0
            atb_x_u = 0.0
            atb_y_u = 0.0
            atb_x_v = 0.0
            atb_y_v = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    de = eta[nb] - ei
                    du = u[nb] - ui
                    dv = v[nb] - vi
                    atb_x_eta += dx * de
                    atb_y_eta += dy * de
                    atb_x_u += dx * du
                    atb_y_u += dy * du
                    atb_x_v += dx * dv
                    atb_y_v += dy * dv
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            gx_eta[i] = g00 * atb_x_eta + g01 * atb_y_eta
            gy_eta[i] = g10 * atb_x_eta + g11 * atb_y_eta
            gx_u[i] = g00 * atb_x_u + g01 * atb_y_u
            gy_u[i] = g10 * atb_x_u + g11 * atb_y_u
            gx_v[i] = g00 * atb_x_v + g01 * atb_y_v
            gy_v[i] = g10 * atb_x_v + g11 * atb_y_v
            temp_sum += temperature_c[i]
        temp_mean = temp_sum / npix

        # Pass B: second derivative (laplacian) of u and v, fused with the momentum update
        # (Coriolis/pressure-gradient/viscosity/mountain-deflection, all purely local).
        for i in range(npix):
            gxu_i = gx_u[i]
            gyu_i = gy_u[i]
            gxv_i = gx_v[i]
            gyv_i = gy_v[i]
            a_u = 0.0
            b_u = 0.0
            c_u = 0.0
            d_u = 0.0
            a_v = 0.0
            b_v = 0.0
            c_v = 0.0
            d_v = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dgxu = gx_u[nb] - gxu_i
                    dgyu = gy_u[nb] - gyu_i
                    dgxv = gx_v[nb] - gxv_i
                    dgyv = gy_v[nb] - gyv_i
                    a_u += dx * dgxu
                    b_u += dy * dgxu
                    c_u += dx * dgyu
                    d_u += dy * dgyu
                    a_v += dx * dgxv
                    b_v += dy * dgxv
                    c_v += dx * dgyv
                    d_v += dy * dgyv
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            lap_u_i = (g00 * a_u + g01 * b_u) + (g10 * c_u + g11 * d_u)
            lap_v_i = (g00 * a_v + g01 * b_v) + (g10 * c_v + g11 * d_v)

            ui = u[i]
            vi = v[i]
            into_slope = -(ui * normal_x[i] + vi * normal_y[i])
            if into_slope < 0.0:
                into_slope = 0.0
            rate = mountain_deflection_rate * ramp[i] * into_slope
            deflect_ax = rate * (normal_x[i] + tangent_x[i])
            deflect_ay = rate * (normal_y[i] + tangent_y[i])

            du_dt = f[i] * vi - reduced_gravity * gx_eta[i] + viscosity * lap_u_i + deflect_ax
            dv_dt = -f[i] * ui - reduced_gravity * gy_eta[i] + viscosity * lap_v_i + deflect_ay

            u_raw[i] = (ui + dt_s * du_dt) / (1.0 + dt_s * orographic_drag[i])
            v_raw[i] = (vi + dt_s * dv_dt) / (1.0 + dt_s * orographic_drag[i])

        # Pass C: grid_noise_filter of u_raw/v_raw (no is_ocean masking -- wind exists over
        # land too), fused with the tropospheric mass-flux prep.
        for i in range(npix):
            total_u = 0.0
            total_v = 0.0
            count = 0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    total_u += u_raw[nb]
                    total_v += v_raw[nb]
                    count += 1
            if count == 0:
                count = 1
            avg_u = total_u / count
            avg_v = total_v / count
            fu = u_raw[i] + 0.05 * (avg_u - u_raw[i])
            fv = v_raw[i] + 0.05 * (avg_v - v_raw[i])
            filt_u[i] = fu
            filt_v[i] = fv
            flux_u[i] = trop_depth * fu
            flux_v[i] = trop_depth * fv

        # Pass D: divergence of the tropospheric mass flux, fused with the continuity/eta
        # update (eta_target uses temperature_c's OLD value and Pass A's temp_mean).
        for i in range(npix):
            fui = flux_u[i]
            fvi = flux_v[i]
            atbx_fu = 0.0
            atby_fu = 0.0
            atbx_fv = 0.0
            atby_fv = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dfu = flux_u[nb] - fui
                    dfv = flux_v[nb] - fvi
                    atbx_fu += dx * dfu
                    atby_fu += dy * dfu
                    atbx_fv += dx * dfv
                    atby_fv += dy * dfv
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            du_dx = g00 * atbx_fu + g01 * atby_fu
            dv_dy = g10 * atbx_fv + g11 * atby_fv
            flux_divergence = du_dx + dv_dy

            eta_target = eta_temp_coupling * (temperature_c[i] - temp_mean)
            eta_raw[i] = eta[i] - dt_s * flux_divergence + dt_s * eta_thermal_relaxation * (eta_target - eta[i])

        # Pass E: grid_noise_filter(eta_raw) -> this substep's final eta.
        for i in range(npix):
            total = 0.0
            count = 0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    total += eta_raw[neighbours[i, k]]
                    count += 1
            if count == 0:
                count = 1
            avg = total / count
            eta[i] = eta_raw[i] + 0.05 * (avg - eta_raw[i])

        # Pass F: semi_lagrangian_advect of temperature_c and humidity together.
        for i in range(npix):
            ox = -filt_u[i] * dt_s
            oy = -filt_v[i] * dt_s
            sx = world_xyz[i, 0] + (ox * east[i, 0] + oy * north[i, 0]) / radius_m
            sy = world_xyz[i, 1] + (ox * east[i, 1] + oy * north[i, 1]) / radius_m
            sz = world_xyz[i, 2] + (ox * east[i, 2] + oy * north[i, 2]) / radius_m
            norm = math.sqrt(sx * sx + sy * sy + sz * sz)
            sx /= norm
            sy /= norm
            sz /= norm
            if sz > 1.0:
                sz = 1.0
            elif sz < -1.0:
                sz = -1.0
            src_lon = math.atan2(sy, sx)
            src_lat = math.asin(sz)
            src_pix = ang2pix_nest_scalar(nside, order, src_lon, src_lat)
            advected_temp[i] = temperature_c[src_pix]
            advected_hum[i] = humidity[src_pix]

        # Pass G: first-derivative gather of advected_temp and advected_hum together.
        for i in range(npix):
            ti = advected_temp[i]
            hi = advected_hum[i]
            atbx_t = 0.0
            atby_t = 0.0
            atbx_h = 0.0
            atby_h = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dt_ = advected_temp[nb] - ti
                    dh_ = advected_hum[nb] - hi
                    atbx_t += dx * dt_
                    atby_t += dy * dt_
                    atbx_h += dx * dh_
                    atby_h += dy * dh_
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            gx_temp[i] = g00 * atbx_t + g01 * atby_t
            gy_temp[i] = g10 * atbx_t + g11 * atby_t
            gx_hum[i] = g00 * atbx_h + g01 * atby_h
            gy_hum[i] = g10 * atbx_h + g11 * atby_h

        # Pass H: second derivative (laplacian) of advected_temp/advected_hum, fused with
        # radiative relaxation, evaporation/condensation, and precipitation -- this
        # substep's final write. `humidity[i]` here is still the pre-advection value
        # (matches the original loop, which computes excess/condensed before overwriting
        # `humidity` via advection).
        for i in range(npix):
            gxt_i = gx_temp[i]
            gyt_i = gy_temp[i]
            gxh_i = gx_hum[i]
            gyh_i = gy_hum[i]
            a_t = 0.0
            b_t = 0.0
            c_t = 0.0
            d_t = 0.0
            a_h = 0.0
            b_h = 0.0
            c_h = 0.0
            d_h = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dgxt = gx_temp[nb] - gxt_i
                    dgyt = gy_temp[nb] - gyt_i
                    dgxh = gx_hum[nb] - gxh_i
                    dgyh = gy_hum[nb] - gyh_i
                    a_t += dx * dgxt
                    b_t += dy * dgxt
                    c_t += dx * dgyt
                    d_t += dy * dgyt
                    a_h += dx * dgxh
                    b_h += dy * dgxh
                    c_h += dx * dgyh
                    d_h += dy * dgyh
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            lap_t = (g00 * a_t + g01 * b_t) + (g10 * c_t + g11 * d_t)
            lap_h = (g00 * a_h + g01 * b_h) + (g10 * c_h + g11 * d_h)

            diffused_temp = advected_temp[i] + dt_s * temp_diffusivity * lap_t
            final_temp = diffused_temp + dt_s * radiative_relaxation * (equilibrium_temperature_c[i] - diffused_temp)

            evap_source = ocean_evap_source if is_ocean[i] else land_evap_source
            saturation_ceiling = final_temp / evap_reference_temp
            if saturation_ceiling < min_evap_ceiling:
                saturation_ceiling = min_evap_ceiling
            elif saturation_ceiling > max_evap_ceiling:
                saturation_ceiling = max_evap_ceiling
            excess = humidity[i] - saturation_ceiling
            if excess < 0.0:
                excess = 0.0
            condensed = condensation_rate * excess

            diffused_hum = advected_hum[i] + dt_s * humidity_diffusivity * lap_h
            final_hum = diffused_hum + dt_s * (evap_source - condensed)
            if final_hum < 0.0:
                final_hum = 0.0

            temperature_c[i] = final_temp
            humidity[i] = final_hum
            precipitation_mm[i] = condensed * precip_conversion
            u[i] = filt_u[i]
            v[i] = filt_v[i]

    return u, v, eta, temperature_c, humidity, precipitation_mm


def step_atmosphere_cfd(world: "WorldV2", state: AtmosphereCFDStateV2, seconds: float) -> None:
    del world
    grid = state.grid
    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * EFFECTIVE_TROPOSPHERE_DEPTH_M))
    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, fdh.min_spacing_m(grid), wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)

    orographic_drag = SURFACE_DRAG_BASE_PER_S + SURFACE_DRAG_OROGRAPHIC_PER_S * np.clip(state.elevation_m / MOUNTAIN_OBSTACLE_ELEVATION_M, 0.0, 1.0)
    normal_x, normal_y, tangent_x, tangent_y, ramp = _mountain_deflection_geometry(state.elevation_m, grid)
    f = fluid_dynamics.coriolis_parameter(np.degrees(grid.lat_rad))

    u, v, eta, temperature_c, humidity, precipitation_mm = _atmosphere_substep_loop_kernel(
        state.u.copy(),
        state.v.copy(),
        state.eta.copy(),
        state.temperature_c.copy(),
        state.humidity.copy(),
        state.precipitation_mm.copy(),
        state.is_ocean,
        state.equilibrium_temperature_c,
        orographic_drag.astype(np.float32),
        normal_x,
        normal_y,
        tangent_x,
        tangent_y,
        ramp,
        f,
        grid.neighbours,
        grid.neighbour_valid,
        grid.neighbour_dx_m,
        grid.neighbour_dy_m,
        grid.gradient_inv,
        grid.world_xyz,
        grid.east,
        grid.north,
        fdh.PLANET_RADIUS_M,
        grid.nside,
        grid.order,
        dt_s,
        n_substeps,
        REDUCED_GRAVITY_M_S2,
        EFFECTIVE_TROPOSPHERE_DEPTH_M,
        ETA_TEMPERATURE_COUPLING_M_PER_C,
        ETA_THERMAL_RELAXATION_PER_S,
        VISCOSITY_M2_S,
        MOUNTAIN_DEFLECTION_RATE_PER_S,
        TEMPERATURE_DIFFUSIVITY_M2_S,
        RADIATIVE_RELAXATION_PER_S,
        OCEAN_EVAPORATION_SOURCE_PER_S,
        LAND_EVAPORATION_SOURCE_PER_S,
        climate.EVAPORATION_REFERENCE_TEMP_C,
        climate.MIN_EVAPORATION_CEILING,
        climate.MAX_EVAPORATION_CEILING,
        CONDENSATION_RATE_PER_S,
        HUMIDITY_DIFFUSIVITY_M2_S,
        PRECIP_CONDENSATION_TO_MM,
    )
    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.humidity = humidity
    state.precipitation_mm = precipitation_mm
    state.elapsed_seconds += dt_s * n_substeps
