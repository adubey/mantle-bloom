"""Shallow-water ocean-current solver over the world's HEALPix grid -- same structure as
`atmosphere_cfd.py`'s own wind solver (`(npix,)` flat arrays via `fluid_dynamics_healpix`
primitives), plus its own sediment pickup/settling tracking."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numba import njit

from . import climate, fluid_dynamics
from . import atmosphere_cfd
from . import fluid_dynamics_healpix as fdh
from . import healpix_grid
from .healpix_grid import ang2pix_nest_scalar

if TYPE_CHECKING:
    from .world import World

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
class OceanCFDState:
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

    def resample_scalar_to_equirect(self, field: np.ndarray, height: int, width: int) -> np.ndarray:
        """See atmosphere_cfd.AtmosphereCFDState's own matching method."""
        return healpix_grid.resample_to_equirect(self.grid, field, height, width)

    def resample_uv_to_equirect(self, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
        u = self.resample_scalar_to_equirect(self.u, height, width)
        v = self.resample_scalar_to_equirect(self.v, height, width)
        return u, v


def init_ocean_cfd(world: "World", terrain: climate.ClimateFields, atmosphere_state: atmosphere_cfd.AtmosphereCFDState) -> OceanCFDState:
    """Same HEALPix grid `atmosphere_cfd.init_atmosphere_cfd` already built for this
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

    return OceanCFDState(
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


def refresh_forcing(world: "World", state: OceanCFDState, terrain: climate.ClimateFields) -> None:
    lat_rows = terrain.lat_deg
    state.is_ocean = healpix_grid.resample_from_equirect(state.grid, terrain.is_ocean.astype(np.float64), lat_rows) > 0.5
    state.elevation_m = healpix_grid.resample_from_equirect(state.grid, terrain.elevation_m, lat_rows)
    state.depth_m = np.where(state.is_ocean, np.clip(-state.elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0).astype(np.float32)
    state.wind_u = world.atmosphere_cfd_state.u.copy()
    state.wind_v = world.atmosphere_cfd_state.v.copy()


@njit(cache=True, fastmath=True)
def _ocean_substep_loop_kernel(
    u,
    v,
    eta,
    temperature_c,
    sediment_concentration,
    sediment_deposited_m,
    wind_u,
    wind_v,
    depth_m,
    is_ocean,
    coastal,
    baseline_temperature_c,
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
    fast_substeps_per_block,
    drag,
    reduced_gravity,
    wind_stress_coeff,
    mixed_layer_depth,
    viscosity,
    temp_diffusivity,
    temp_relaxation,
    sediment_diffusivity,
    sediment_pickup_coeff,
    sediment_settle_speed_thresh,
    sediment_settle_rate,
    sediment_deposit_coeff,
):
    """Subcycled fast/slow split of `step_ocean_cfd`'s substep loop, fused into one compiled
    function that also owns the outer-block/fast-substep looping itself.

    The gravity-wave system (pressure-gradient + Coriolis + continuity) is the only thing
    that actually needs `dt_s` (the wave-speed CFL dt `cfl_substeps` already computes) --
    viscosity and the temperature/sediment advection-diffusion are stable at a dt up to
    `fast_substeps_per_block`x larger (`fluid_dynamics.diffusion_stable_dt`, sized from this
    solver's own diffusivity constants; see `step_ocean_cfd`). So each of the
    `ceil(n_substeps / fast_substeps_per_block)` outer blocks below runs:
      - a FAST sub-cycle, `dt_s` at a time, of exactly the terms the wave CFL bounds
        (eta/u/v's own numbers are therefore identical to the pre-split kernel's -- the
        part of the scheme that needs to be numerically careful is untouched);
      - one SLOW update -- viscosity, then advection, then diffusion/relaxation/pickup-settle
        of temperature_c/sediment_concentration -- at `dt_outer` (that block's actual fast-
        substep count times `dt_s`, so the last, possibly-short, block still lands on exactly
        `n_substeps * dt_s` total), using the fast sub-cycle's freshly-updated u/v (slow-
        after-fast operator splitting, so tracer advection reads the current velocity, not a
        stale one).
    Every per-pixel expression below is unchanged from the original single-cadence kernel --
    only which cadence (`dt_s` vs `dt_outer`) each term is evaluated at, and how the momentum
    equation's viscosity term separates from its (still-fast) drag term, changed."""
    npix = u.shape[0]
    n_neighbours = neighbours.shape[1]

    gx_eta = np.empty(npix, dtype=np.float32)
    gy_eta = np.empty(npix, dtype=np.float32)
    u_raw = np.empty(npix, dtype=np.float32)
    v_raw = np.empty(npix, dtype=np.float32)
    filt_u = np.empty(npix, dtype=np.float32)
    filt_v = np.empty(npix, dtype=np.float32)
    flux_u = np.empty(npix, dtype=np.float32)
    flux_v = np.empty(npix, dtype=np.float32)
    eta_raw = np.empty(npix, dtype=np.float32)
    gx_u = np.empty(npix, dtype=np.float32)
    gy_u = np.empty(npix, dtype=np.float32)
    gx_v = np.empty(npix, dtype=np.float32)
    gy_v = np.empty(npix, dtype=np.float32)
    pickup = np.empty(npix, dtype=np.float32)
    settle = np.empty(npix, dtype=np.float32)
    advected_temp = np.empty(npix, dtype=np.float32)
    advected_sed = np.empty(npix, dtype=np.float32)
    gx_temp = np.empty(npix, dtype=np.float32)
    gy_temp = np.empty(npix, dtype=np.float32)
    gx_sed = np.empty(npix, dtype=np.float32)
    gy_sed = np.empty(npix, dtype=np.float32)

    substeps_done = 0
    while substeps_done < n_substeps:
        block_count = n_substeps - substeps_done
        if block_count > fast_substeps_per_block:
            block_count = fast_substeps_per_block
        dt_outer = dt_s * block_count

        for _fast in range(block_count):
            # Fast Pass 1: gradient of eta only, fused with the local momentum update
            # (Coriolis + pressure-gradient + wind stress + drag -- all purely local given
            # gx_eta/gy_eta, no extra neighbour access) -- viscosity is *not* here any more,
            # it moved to the slow block below.
            for i in range(npix):
                ei = eta[i]
                atb_x = 0.0
                atb_y = 0.0
                for k in range(n_neighbours):
                    if neighbour_valid[i, k]:
                        nb = neighbours[i, k]
                        dx = neighbour_dx_m[i, k]
                        dy = neighbour_dy_m[i, k]
                        de = eta[nb] - ei
                        atb_x += dx * de
                        atb_y += dy * de
                g00 = gradient_inv[i, 0, 0]
                g01 = gradient_inv[i, 0, 1]
                g10 = gradient_inv[i, 1, 0]
                g11 = gradient_inv[i, 1, 1]
                gxe_i = g00 * atb_x + g01 * atb_y
                gye_i = g10 * atb_x + g11 * atb_y

                wind_speed_i = math.sqrt(wind_u[i] * wind_u[i] + wind_v[i] * wind_v[i])
                tau_u_i = (wind_stress_coeff * wind_speed_i * wind_u[i]) / mixed_layer_depth
                tau_v_i = (wind_stress_coeff * wind_speed_i * wind_v[i]) / mixed_layer_depth

                du_dt = f[i] * v[i] - reduced_gravity * gxe_i + tau_u_i
                dv_dt = -f[i] * u[i] - reduced_gravity * gye_i + tau_v_i

                u_raw[i] = (u[i] + dt_s * du_dt) / (1.0 + dt_s * drag)
                v_raw[i] = (v[i] + dt_s * dv_dt) / (1.0 + dt_s * drag)

            # Fast Pass 2: grid_noise_filter of u_raw/v_raw (masked by is_ocean) -> this
            # substep's final u/v, fused with flux_u/flux_v prep for the continuity update.
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
                if not is_ocean[i]:
                    fu = 0.0
                    fv = 0.0
                filt_u[i] = fu
                filt_v[i] = fv
                flux_u[i] = depth_m[i] * fu
                flux_v[i] = depth_m[i] * fv

            # Fast Pass 3: divergence of (depth_m * u, depth_m * v), fused with the
            # continuity update and writing this substep's final u/v.
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
                eta_raw[i] = eta[i] - dt_s * (du_dx + dv_dy)
                u[i] = filt_u[i]
                v[i] = filt_v[i]

            # Fast Pass 4: grid_noise_filter(eta_raw), masked -> this substep's final eta.
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
                e = eta_raw[i] + 0.05 * (avg - eta_raw[i])
                eta[i] = e if is_ocean[i] else 0.0

        # Slow Pass 1: first-derivative gather of u, v (post fast-sub-cycle).
        for i in range(npix):
            ui = u[i]
            vi = v[i]
            atb_x_u = 0.0
            atb_y_u = 0.0
            atb_x_v = 0.0
            atb_y_v = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    du = u[nb] - ui
                    dv = v[nb] - vi
                    atb_x_u += dx * du
                    atb_y_u += dy * du
                    atb_x_v += dx * dv
                    atb_y_v += dy * dv
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            gx_u[i] = g00 * atb_x_u + g01 * atb_y_u
            gy_u[i] = g10 * atb_x_u + g11 * atb_y_u
            gx_v[i] = g00 * atb_x_v + g01 * atb_y_v
            gy_v[i] = g10 * atb_x_v + g11 * atb_y_v

        # Slow Pass 2: second derivative (laplacian) of u/v -> explicit viscosity increment
        # over dt_outer, fused with the sediment pickup/settle terms (purely local given the
        # now-final u/v for this block and the pre-advection sediment_concentration).
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

            new_u = u[i] + dt_outer * viscosity * lap_u_i
            new_v = v[i] + dt_outer * viscosity * lap_v_i
            if not is_ocean[i]:
                new_u = 0.0
                new_v = 0.0
            u[i] = new_u
            v[i] = new_v

            speed_i = math.sqrt(new_u * new_u + new_v * new_v)
            pickup[i] = sediment_pickup_coeff * speed_i if coastal[i] else 0.0
            settling_factor = 1.0 - speed_i / sediment_settle_speed_thresh
            if settling_factor < 0.0:
                settling_factor = 0.0
            elif settling_factor > 1.0:
                settling_factor = 1.0
            settle[i] = sediment_settle_rate * sediment_concentration[i] * settling_factor

        # Slow Pass 3: semi_lagrangian_advect of temperature_c and sediment_concentration
        # together, over dt_outer, using this block's final u/v (unconditionally stable
        # regardless of dt -- it's a backward trace via ang2pix, not a stencil, so coarsening
        # its cadence to dt_outer costs nothing it needed).
        for i in range(npix):
            ox = -u[i] * dt_outer
            oy = -v[i] * dt_outer
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
            advected_sed[i] = sediment_concentration[src_pix]

        # Slow Pass 4: first-derivative gather of advected_temp and advected_sed together.
        for i in range(npix):
            ti = advected_temp[i]
            si = advected_sed[i]
            atbx_t = 0.0
            atby_t = 0.0
            atbx_s = 0.0
            atby_s = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dt_ = advected_temp[nb] - ti
                    ds_ = advected_sed[nb] - si
                    atbx_t += dx * dt_
                    atby_t += dy * dt_
                    atbx_s += dx * ds_
                    atby_s += dy * ds_
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            gx_temp[i] = g00 * atbx_t + g01 * atby_t
            gy_temp[i] = g10 * atbx_t + g11 * atby_t
            gx_sed[i] = g00 * atbx_s + g01 * atby_s
            gy_sed[i] = g10 * atbx_s + g11 * atby_s

        # Slow Pass 5: second derivative (laplacian) of advected_temp/advected_sed, fused with
        # every remaining local term (diffusion, radiative relaxation, pickup/settle, the
        # is_ocean mask, and sediment deposit accumulation) -- this block's final write, all
        # over dt_outer.
        for i in range(npix):
            gxt_i = gx_temp[i]
            gyt_i = gy_temp[i]
            gxs_i = gx_sed[i]
            gys_i = gy_sed[i]
            a_t = 0.0
            b_t = 0.0
            c_t = 0.0
            d_t = 0.0
            a_s = 0.0
            b_s = 0.0
            c_s = 0.0
            d_s = 0.0
            for k in range(n_neighbours):
                if neighbour_valid[i, k]:
                    nb = neighbours[i, k]
                    dx = neighbour_dx_m[i, k]
                    dy = neighbour_dy_m[i, k]
                    dgxt = gx_temp[nb] - gxt_i
                    dgyt = gy_temp[nb] - gyt_i
                    dgxs = gx_sed[nb] - gxs_i
                    dgys = gy_sed[nb] - gys_i
                    a_t += dx * dgxt
                    b_t += dy * dgxt
                    c_t += dx * dgyt
                    d_t += dy * dgyt
                    a_s += dx * dgxs
                    b_s += dy * dgxs
                    c_s += dx * dgys
                    d_s += dy * dgys
            g00 = gradient_inv[i, 0, 0]
            g01 = gradient_inv[i, 0, 1]
            g10 = gradient_inv[i, 1, 0]
            g11 = gradient_inv[i, 1, 1]
            lap_t = (g00 * a_t + g01 * b_t) + (g10 * c_t + g11 * d_t)
            lap_s = (g00 * a_s + g01 * b_s) + (g10 * c_s + g11 * d_s)

            diffused_temp = advected_temp[i] + dt_outer * temp_diffusivity * lap_t
            temperature_c[i] = diffused_temp + dt_outer * temp_relaxation * (baseline_temperature_c[i] - diffused_temp)

            diffused_sed = advected_sed[i] + dt_outer * sediment_diffusivity * lap_s
            final_sed = diffused_sed + dt_outer * (pickup[i] - settle[i])
            if final_sed < 0.0:
                final_sed = 0.0
            sediment_concentration[i] = final_sed if is_ocean[i] else 0.0

            sediment_deposited_m[i] += dt_outer * settle[i] * sediment_deposit_coeff

        substeps_done += block_count

    return u, v, eta, temperature_c, sediment_concentration, sediment_deposited_m


def step_ocean_cfd(world: "World", state: OceanCFDState, seconds: float) -> None:
    del world
    grid = state.grid
    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * max(float(state.depth_m.max()), MIN_DEPTH_M)))
    coastal = fdh.coastal_ocean_mask(state.is_ocean, grid)

    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, fdh.min_spacing_m(grid), wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)
    f = fluid_dynamics.coriolis_parameter(np.degrees(grid.lat_rad))

    max_diffusivity = max(VISCOSITY_M2_S, TEMPERATURE_DIFFUSIVITY_M2_S, SEDIMENT_DIFFUSIVITY_M2_S)
    dt_outer_limit = fluid_dynamics.diffusion_stable_dt(fdh.min_spacing_m(grid), max_diffusivity)
    fast_substeps_per_block = max(1, int(dt_outer_limit // dt_s))

    u, v, eta, temperature_c, sediment_concentration, sediment_deposited_m = _ocean_substep_loop_kernel(
        state.u.copy(),
        state.v.copy(),
        state.eta.copy(),
        state.temperature_c.copy(),
        state.sediment_concentration.copy(),
        state.sediment_deposited_m.copy(),
        state.wind_u,
        state.wind_v,
        state.depth_m,
        state.is_ocean,
        coastal,
        state.baseline_temperature_c,
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
        fast_substeps_per_block,
        BOTTOM_DRAG_PER_S,
        REDUCED_GRAVITY_M_S2,
        WIND_STRESS_COEFFICIENT,
        MIXED_LAYER_DEPTH_M,
        VISCOSITY_M2_S,
        TEMPERATURE_DIFFUSIVITY_M2_S,
        TEMPERATURE_RELAXATION_PER_S,
        SEDIMENT_DIFFUSIVITY_M2_S,
        SEDIMENT_PICKUP_COEFFICIENT,
        SEDIMENT_SETTLING_SPEED_THRESHOLD_M_S,
        SEDIMENT_SETTLING_RATE_PER_S,
        SEDIMENT_DEPOSIT_DEPTH_COEFFICIENT,
    )
    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.sediment_concentration = sediment_concentration
    state.sediment_deposited_m = sediment_deposited_m
    state.elapsed_seconds += dt_s * n_substeps

    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.sediment_concentration = sediment_concentration
    state.sediment_deposited_m = sediment_deposited_m
