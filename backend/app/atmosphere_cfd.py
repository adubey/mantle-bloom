"""Atmospheric Fluid Dynamics mode: a genuine time-integrated shallow-water simulation of
wind, replacing climate.py's own diagnostic latitude-banded wind heuristic with real
prognostic state that persists and evolves step to step -- see
docs/simulation-model.md#ocean-atmospheric-fluid-dynamics for the shared design rationale
(reduced gravity, freeze-on-entry, substepping) with ocean_cfd.py.

**Inputs, matching the user's own spec.** Coriolis force and elevation (both baked directly
into the momentum equation -- Coriolis in the usual way, elevation via orographic deflection
around high terrain and a lapse-rate-cooled radiative-equilibrium target temperature, see
`_deflect_around_mountains`/`_equilibrium_temperature`), plus humidity as a genuinely
prognostic, wind-advected field (evaporation source over moist surfaces, a precipitation sink
where locally saturated) rather than climate.py's own one-shot diagnostic formula."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from . import climate, fluid_dynamics

if TYPE_CHECKING:
    from .world import World

# Reduced gravity for the atmosphere's own equivalent-barotropic layer -- larger than the
# ocean's own REDUCED_GRAVITY_M_S2 (see ocean_cfd.py) since real atmospheric circulation
# genuinely moves faster than ocean currents; tuned so gravity-wave speed lands in the
# few-tens-of-m/s range (a realistic large-scale atmospheric wave speed) rather than either
# the ocean's much gentler ~1 m/s or real dry-atmosphere gravity waves' ~300 m/s (which would
# force substep counts back up to the same range reduced gravity exists to avoid).
REDUCED_GRAVITY_M_S2 = 0.5
EFFECTIVE_TROPOSPHERE_DEPTH_M = 8000.0
# Relates local temperature anomaly to geopotential-height anomaly -- the real atmospheric
# "hypsometric" relationship (a warmer air column is thicker, so upper-level geopotential
# height reads higher over warm regions), the actual thermal forcing driving this mode's
# pressure-gradient-like `eta` field, standing in for a full 3D pressure/density solve the
# same way climate.py's own wind model stands in for one (see that module's own docstring).
ETA_TEMPERATURE_COUPLING_M_PER_C = 8.0
ETA_THERMAL_RELAXATION_PER_S = 3.0e-5

# Tuned the same way ocean_cfd.BOTTOM_DRAG_PER_S was (see its own comment): strong enough to
# damp the lightly-damped, Coriolis-dominated near-inertial transient at high latitude within
# a real run's own timescale (confirmed directly: max wind speed over a 10-step, 60-hour run
# settles around ~1 m/s at this value, versus a still-ringing ~2.3 m/s at the original 1.5e-6),
# applied semi-implicitly (see step_atmosphere_cfd) so, like ocean_cfd's own drag, raising it
# further is never a stability risk, only a "how energetic does the wind look" tuning choice.
SURFACE_DRAG_BASE_PER_S = 1.0e-5
SURFACE_DRAG_OROGRAPHIC_PER_S = 8.0e-6
# Extra damping poleward of fluid_dynamics.POLAR_FILTER_START_LAT_DEG -- see
# fluid_dynamics.polar_sponge_drag_per_s's own docstring and ocean_cfd.py's matching constant
# for why the polar cap needs this on top of ordinary surface drag.
POLAR_SPONGE_MAX_DRAG_PER_S = 5.0e-4
VISCOSITY_M2_S = 6.0e4

MOUNTAIN_OBSTACLE_ELEVATION_M = 2000.0
MOUNTAIN_OBSTACLE_RAMP_M = 800.0
# A genuine per-second damping *rate* (see _mountain_deflection_tendency), not a one-shot
# multiplicative state overwrite -- an earlier version applied "cancel the into-slope
# component and redirect it tangentially, sped up by 1.4x" directly to (u, v) every substep,
# which is stable as a single diagnostic transform (climate.py's own _mountain_deflection
# only ever runs once per climate computation) but compounds geometrically when the same
# transform re-applies every substep to an already-redirected velocity -- confirmed directly:
# an early build blew up to absurd speeds within a single UI "Step" from exactly this.
# Framed as a *tendency* instead (added into du_dt/dv_dt alongside Coriolis/pressure-
# gradient/drag), it only ever damps a persistent into-slope flow, the same stability
# property ordinary bottom drag already has.
MOUNTAIN_DEFLECTION_RATE_PER_S = 3.0e-4

TEMPERATURE_DIFFUSIVITY_M2_S = 3.0e4
RADIATIVE_RELAXATION_PER_S = 3.0e-6

HUMIDITY_DIFFUSIVITY_M2_S = 2.0e4
OCEAN_EVAPORATION_SOURCE_PER_S = 4.0e-7
LAND_EVAPORATION_SOURCE_PER_S = 0.6e-7
# Above this, excess humidity condenses out as precipitation each substep -- same evaporation-
# ceiling shape climate.py's own compute_humidity uses (see EVAPORATION_REFERENCE_TEMP_C/
# MIN/MAX_EVAPORATION_CEILING there), reused here as a real prognostic saturation limit
# instead of a one-shot diagnostic formula.
CONDENSATION_RATE_PER_S = 2.0e-6
PRECIP_CONDENSATION_TO_MM = 4.0e6

MAX_SUBSTEPS_PER_STEP = 2000


@dataclass
class AtmosphereCFDState:
    lat_deg: np.ndarray  # (H,)
    lon_deg: np.ndarray  # (W,)
    world_xyz: np.ndarray  # (H, W, 3)
    is_ocean: np.ndarray  # (H, W) bool -- frozen at mode entry, drives evaporation source only
    elevation_m: np.ndarray  # (H, W) -- frozen at mode entry
    u: np.ndarray  # (H, W) eastward wind, m/s
    v: np.ndarray  # (H, W) northward wind, m/s
    eta: np.ndarray  # (H, W) geopotential-height anomaly, m
    temperature_c: np.ndarray  # (H, W)
    equilibrium_temperature_c: np.ndarray  # (H, W) -- fixed radiative-equilibrium target, see module docstring
    humidity: np.ndarray  # (H, W), roughly [0, climate.MAX_EVAPORATION_CEILING]
    precipitation_mm: np.ndarray  # (H, W) -- latest substep's condensation rate, a display field
    elapsed_seconds: float = 0.0


def _equilibrium_temperature(world: "World", fields: climate.ClimateFields) -> np.ndarray:
    """The radiative-equilibrium temperature this mode's own temperature field relaxes
    toward every substep (see RADIATIVE_RELAXATION_PER_S) -- land's solar-heating-plus-lapse-
    rate baseline over land, the zonal water baseline over ocean, reusing climate.py's own
    public formulas directly (same insolation/elevation this world already has, frozen at
    mode entry) rather than re-deriving them."""
    insolation_row = climate.compute_insolation(fields.lat_deg, world.axial_tilt_deg, world.solar_multiplier)
    land = climate.compute_land_temperature(insolation_row, fields.elevation_m)
    ocean = climate.compute_ocean_temperature_baseline(insolation_row, *fields.elevation_m.shape)
    return np.where(fields.is_ocean, ocean, land)


def init_atmosphere_cfd(world: "World") -> AtmosphereCFDState:
    """Snapshots the world's current elevation/temperature/humidity (via climate.py's own
    public pipeline) and starts the atmosphere's wind from World.remembered_wind_u/v (see
    their own docstring) when a prior "atmosphere_cfd" or "ocean_cfd" session left one behind,
    else from that diagnostic snapshot's own wind -- unlike ocean_cfd.py's ocean-at-rest
    start, beginning from a plausible wind field (remembered or diagnostic) avoids a jarring
    dead-calm-to-storm transient the first few substeps would otherwise show."""
    height, width = climate.grid_dimensions(world.climate_density)
    fields = climate.compute_climate(world, height, width)
    equilibrium_temperature_c = _equilibrium_temperature(world, fields)

    return AtmosphereCFDState(
        lat_deg=fields.lat_deg,
        lon_deg=fields.lon_deg,
        world_xyz=fields.world_xyz,
        is_ocean=fields.is_ocean,
        elevation_m=fields.elevation_m,
        u=fields.wind_u.copy() if world.remembered_wind_u is None else world.remembered_wind_u.copy(),
        v=fields.wind_v.copy() if world.remembered_wind_v is None else world.remembered_wind_v.copy(),
        eta=np.zeros((height, width)),
        temperature_c=np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c),
        equilibrium_temperature_c=equilibrium_temperature_c,
        humidity=fields.humidity.copy(),
        precipitation_mm=fields.precipitation_mm.copy(),
    )


def remember_atmosphere_state(world: "World", state: AtmosphereCFDState) -> None:
    """Snapshots this session's final wind onto `world.remembered_wind_u/v` (see their own
    docstring) so a later switch into "ocean_cfd" or back into "atmosphere_cfd" can resume
    from it instead of a fresh climate.py diagnostic."""
    world.remembered_wind_u = state.u.copy()
    world.remembered_wind_v = state.v.copy()


def _mountain_deflection_geometry(elevation_m: np.ndarray, dx_m: np.ndarray, dy_m: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The terrain-derived part of _mountain_deflection_tendency's own math -- the slope
    normal/tangent directions and the obstacle ramp -- which depends only on elevation_m, not
    on wind (u, v). elevation_m is frozen at mode entry (see module docstring) and never
    changes across a step_atmosphere_cfd call's substeps, so callers compute this once per
    step call rather than _mountain_deflection_tendency redundantly rebuilding it (a
    gradient_m call plus several full-grid ops) every single substep."""
    gx, gy = fluid_dynamics.gradient_m(elevation_m, dx_m, dy_m)
    magnitude = np.hypot(gx, gy)
    uphill_x = np.divide(gx, magnitude, out=np.zeros_like(gx), where=magnitude > 1e-12)
    uphill_y = np.divide(gy, magnitude, out=np.zeros_like(gy), where=magnitude > 1e-12)
    # Outward-from-slope normal, and one of its two tangents (a fixed rotation sense -- unlike
    # boundary-current deflection, mesoscale mountain flow-splitting has no Coriolis-preferred
    # side, so either tangent is an equally valid simplification here).
    normal_x, normal_y = -uphill_x, -uphill_y
    tangent_x, tangent_y = -normal_y, normal_x
    ramp = np.clip((elevation_m - (MOUNTAIN_OBSTACLE_ELEVATION_M - MOUNTAIN_OBSTACLE_RAMP_M)) / (2.0 * MOUNTAIN_OBSTACLE_RAMP_M), 0.0, 1.0)
    return normal_x, normal_y, tangent_x, tangent_y, ramp


def _mountain_deflection_tendency(
    u: np.ndarray, v: np.ndarray, geometry: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    """An acceleration (m/s^2) damping wind's into-slope component near high terrain and
    redirecting a matching amount tangentially -- same "cancel and redirect" *shape* as
    climate.py's own `_mountain_deflection`, but expressed here as a per-second rate applied
    to du_dt/dv_dt (see MOUNTAIN_DEFLECTION_RATE_PER_S's own docstring for why a tendency,
    not a direct state overwrite, is required for stability under repeated substepping).
    `geometry` is this step's own _mountain_deflection_geometry(elevation_m, dx_m, dy_m) --
    see its docstring for why callers precompute it once rather than every substep."""
    normal_x, normal_y, tangent_x, tangent_y, ramp = geometry
    into_slope = np.clip(-(u * normal_x + v * normal_y), 0.0, None)
    rate = MOUNTAIN_DEFLECTION_RATE_PER_S * ramp * into_slope
    # Redirects exactly the damped magnitude into the tangent direction (no speedup factor >
    # 1) -- energy-neutral by construction, not a net source, unlike an earlier version (see
    # MOUNTAIN_DEFLECTION_RATE_PER_S's own docstring).
    ax = rate * (normal_x + tangent_x)
    ay = rate * (normal_y + tangent_y)
    return ax, ay


def step_atmosphere_cfd(world: "World", state: AtmosphereCFDState, seconds: float) -> None:
    """Advances `state` by `seconds` of real time, in as many CFL-stable substeps as the
    current grid/reduced-gravity wave speed demand (see fluid_dynamics.cfl_substeps). `world`
    is accepted (unused directly) to match ocean_cfd.step_ocean_cfd's own `(world, ...)`
    calling convention -- see that function's docstring for why."""
    del world
    height, width = state.elevation_m.shape
    dx_m, dy_m = fluid_dynamics.grid_spacing_m(state.lat_deg, height, width)

    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * EFFECTIVE_TROPOSPHERE_DEPTH_M))
    min_spacing_m = fluid_dynamics.stable_min_spacing_m(dx_m, dy_m)
    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, min_spacing_m, wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)

    orographic_drag = (
        SURFACE_DRAG_BASE_PER_S
        + SURFACE_DRAG_OROGRAPHIC_PER_S * np.clip(state.elevation_m / MOUNTAIN_OBSTACLE_ELEVATION_M, 0.0, 1.0)
        + fluid_dynamics.polar_sponge_drag_per_s(state.lat_deg, POLAR_SPONGE_MAX_DRAG_PER_S)
    )
    # elevation_m/lat_deg/width are all fixed for the whole step call, so both of these
    # (unlike u/v/dt_s below) are the same every substep -- see their own docstrings.
    mountain_geom = _mountain_deflection_geometry(state.elevation_m, dx_m, dy_m)
    advect_geom = fluid_dynamics.advection_geometry(state.lat_deg, width)
    f = fluid_dynamics.coriolis_parameter(state.lat_deg)[:, None]

    u, v, eta = state.u, state.v, state.eta
    temperature_c = state.temperature_c
    humidity = state.humidity
    precipitation_mm = state.precipitation_mm

    for _ in range(n_substeps):
        deta_dx, deta_dy = fluid_dynamics.gradient_m(eta, dx_m, dy_m)

        # Drag (surface/orographic + the polar sponge) is applied semi-implicitly, same
        # unconditionally-stable-regardless-of-dt*drag reasoning as ocean_cfd.step_ocean_cfd's
        # own matching comment.
        deflect_ax, deflect_ay = _mountain_deflection_tendency(u, v, mountain_geom)
        du_dt = f * v - REDUCED_GRAVITY_M_S2 * deta_dx + VISCOSITY_M2_S * fluid_dynamics.laplacian_m(u, dx_m, dy_m) + deflect_ax
        dv_dt = -f * u - REDUCED_GRAVITY_M_S2 * deta_dy + VISCOSITY_M2_S * fluid_dynamics.laplacian_m(v, dx_m, dy_m) + deflect_ay
        u = (u + dt_s * du_dt) / (1.0 + dt_s * orographic_drag)
        v = (v + dt_s * dv_dt) / (1.0 + dt_s * orographic_drag)
        u = fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(u), state.lat_deg)
        v = fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(v), state.lat_deg)

        eta_target = ETA_TEMPERATURE_COUPLING_M_PER_C * (temperature_c - float(temperature_c.mean()))
        flux_divergence = fluid_dynamics.divergence_m(EFFECTIVE_TROPOSPHERE_DEPTH_M * u, EFFECTIVE_TROPOSPHERE_DEPTH_M * v, dx_m, dy_m)
        eta = eta - dt_s * flux_divergence + dt_s * ETA_THERMAL_RELAXATION_PER_S * (eta_target - eta)
        eta = fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(eta), state.lat_deg)

        temperature_c = fluid_dynamics.semi_lagrangian_advect(temperature_c, u, v, dt_s, advect_geom)
        temperature_c = temperature_c + dt_s * TEMPERATURE_DIFFUSIVITY_M2_S * fluid_dynamics.laplacian_m(temperature_c, dx_m, dy_m)
        temperature_c = temperature_c + dt_s * RADIATIVE_RELAXATION_PER_S * (state.equilibrium_temperature_c - temperature_c)

        evap_source = np.where(state.is_ocean, OCEAN_EVAPORATION_SOURCE_PER_S, LAND_EVAPORATION_SOURCE_PER_S)
        saturation_ceiling = np.clip(temperature_c / climate.EVAPORATION_REFERENCE_TEMP_C, climate.MIN_EVAPORATION_CEILING, climate.MAX_EVAPORATION_CEILING)
        excess = np.clip(humidity - saturation_ceiling, 0.0, None)
        condensed = CONDENSATION_RATE_PER_S * excess

        humidity = fluid_dynamics.semi_lagrangian_advect(humidity, u, v, dt_s, advect_geom)
        humidity = humidity + dt_s * HUMIDITY_DIFFUSIVITY_M2_S * fluid_dynamics.laplacian_m(humidity, dx_m, dy_m)
        humidity = np.clip(humidity + dt_s * (evap_source - condensed), 0.0, None)
        precipitation_mm = condensed * PRECIP_CONDENSATION_TO_MM

        state.elapsed_seconds += dt_s

    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.humidity = humidity
    state.precipitation_mm = precipitation_mm
