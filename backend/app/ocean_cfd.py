"""Ocean Fluid Dynamics mode: a genuine time-integrated shallow-water simulation of ocean
currents, replacing climate.py's own diagnostic (recomputed-from-scratch-every-call) Ekman-
current heuristic with real prognostic state that persists and evolves step to step -- see
docs/simulation-model.md#ocean-atmospheric-fluid-dynamics for the full design rationale
(reduced gravity, freeze-on-entry, substepping) shared with atmosphere_cfd.py.

**Inputs, matching the user's own spec.** Coriolis force and land (both baked into the
momentum equation itself -- see step_ocean_cfd), wind (sampled once from the world's current
climate at mode-entry and held fixed for the session -- the atmosphere isn't independently
simulated in this mode, so there's no other source for it), and temperature (both an input,
via the wind-driven baseline this mode starts from, and an output the currents themselves
advect).

**Currents move water along with temperature and sediment.** Every substep advects
`temperature_c` and `sediment_concentration` along the same velocity field that's just been
updated -- literally "carried along," not two separately-evolving fields that merely happen
to share a grid. Sediment is a visual tracer only (advected/diffused, with a coastal pickup
source and a slow-flow settling sink) -- it never mutates `world.plates`' persistent
elevation, but `sediment_deposited_m` tracks cumulative settling per cell for its own display
layer, per the user's explicit "keep track of where it would deposit" ask."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from . import climate, fluid_dynamics

if TYPE_CHECKING:
    from .world import World

# Reduced-gravity shallow water: using real ocean depth (~4000m) and g=9.8 gives a gravity-
# wave speed near 200 m/s, which at planet-grid resolution forces the CFL-stable substep down
# near 100 seconds -- correct, but far more substeps per UI "Step" than the barotropic
# circulation patterns this mode actually displays (gyres, boundary currents, coastal
# upwelling) need to look right. Reduced gravity is the standard technique basin-circulation
# models use to keep wave speed -- and thus substep count -- tractable while the solver still
# genuinely integrates Coriolis/wind-stress/topography forcing every substep, not a physics
# shortcut being hidden (see module docstring's own doc reference).
REDUCED_GRAVITY_M_S2 = 0.03
MIN_DEPTH_M = 200.0
MAX_DEPTH_M = 6000.0

# Bulk kinematic wind-stress coefficient (Cd * rho_air/rho_water lumped into one constant --
# this codebase doesn't track real fluid densities anywhere else either, e.g. climate.py's
# own SUNLIGHT/MERIDIONAL_BASE_SPEED constants are similarly tuned rather than derived from
# first principles) -- `WIND_STRESS_COEFFICIENT * wind_speed * wind_component` is a kinematic
# stress (m^2/s^2), which still needs dividing by a real mixed-layer depth to become a
# surface-water *acceleration*: real wind stress drives roughly the top MIXED_LAYER_DEPTH_M
# of the ocean directly, not the full water column, so spreading the same force over a
# realistic depth (rather than none at all) is what keeps equilibrium current speeds in the
# few-tens-of-cm/s range real wind-driven surface currents show, instead of blowing up to
# absurd speeds within the first simulated day.
WIND_STRESS_COEFFICIENT = 1.5e-6
MIXED_LAYER_DEPTH_M = 50.0
# At high latitude (where the real Coriolis parameter f is largest), a weak drag leaves the
# system in a lightly-damped, Coriolis-dominated regime -- confirmed directly during
# development: at BOTTOM_DRAG_PER_S = 5e-6 (a ~2.3-day spin-down), a persistent near-inertial
# oscillation at ~70 deg latitude (well equatorward of the polar filter zone, so unrelated to
# POLAR_SPONGE_MAX_DRAG_PER_S below) still hadn't decayed after 10 simulated days, reaching
# tens of m/s -- an order of magnitude past any real wind-driven surface current. This value
# (a ~14-hour spin-down) damps that transient fast enough to stay physically plausible
# (confirmed: max speed drops to ~0.1 m/s over the same 10-day run) without needing per-
# latitude tuning, since the drag term is applied semi-implicitly (see step_ocean_cfd) and so
# is unconditionally stable regardless of its magnitude.
BOTTOM_DRAG_PER_S = 2.0e-5
# Extra damping poleward of fluid_dynamics.POLAR_FILTER_START_LAT_DEG, on top of
# BOTTOM_DRAG_PER_S -- see fluid_dynamics.polar_sponge_drag_per_s's own docstring for why the
# polar cap specifically (as opposed to the merely-high-latitude band BOTTOM_DRAG_PER_S above
# already handles) needs its own extra damping.
POLAR_SPONGE_MAX_DRAG_PER_S = 2.0e-4
VISCOSITY_M2_S = 4.0e4
TEMPERATURE_DIFFUSIVITY_M2_S = 2.0e4
# Newtonian relaxation toward the wind-driven baseline temperature this mode started from --
# without a real heat-flux/radiation model, temperature would otherwise drift unboundedly
# over many steps; this is the same "accepted simplification, not a bug" this codebase
# already applies elsewhere (see climate.py's own module docstring).
TEMPERATURE_RELAXATION_PER_S = 1.0e-6

SEDIMENT_DIFFUSIVITY_M2_S = 1.5e4
# Classic transport-capacity model: pickup scales with how fast the water is moving past an
# erodible coast, settling kicks in once the local current drops below a threshold speed.
SEDIMENT_PICKUP_COEFFICIENT = 2.0e-4
SEDIMENT_SETTLING_SPEED_THRESHOLD_M_S = 0.08
SEDIMENT_SETTLING_RATE_PER_S = 3.0e-6
# Converts settled concentration (arbitrary units) into an equivalent deposited depth --
# tuned only so sediment_deposited_m reads as a plausible, slowly-growing meters-scale field
# over a real session, not calibrated against a real sediment density.
SEDIMENT_DEPOSIT_DEPTH_COEFFICIENT = 0.5

MAX_SUBSTEPS_PER_STEP = 2000


@dataclass
class OceanCFDState:
    lat_deg: np.ndarray  # (H,)
    lon_deg: np.ndarray  # (W,)
    world_xyz: np.ndarray  # (H, W, 3)
    is_ocean: np.ndarray  # (H, W) bool
    elevation_m: np.ndarray  # (H, W) -- frozen at mode entry
    depth_m: np.ndarray  # (H, W) -- effective shallow-water layer depth, zero on land
    u: np.ndarray  # (H, W) eastward current, m/s
    v: np.ndarray  # (H, W) northward current, m/s
    eta: np.ndarray  # (H, W) sea-surface-height anomaly, m
    temperature_c: np.ndarray  # (H, W)
    baseline_temperature_c: np.ndarray  # (H, W) -- fixed relaxation target, see module docstring
    sediment_concentration: np.ndarray  # (H, W), arbitrary units >= 0
    sediment_deposited_m: np.ndarray  # (H, W), cumulative, tracking-only -- never touches world.plates
    wind_u: np.ndarray  # (H, W) eastward, fixed for the session
    wind_v: np.ndarray  # (H, W) northward, fixed for the session
    elapsed_seconds: float = 0.0


def init_ocean_cfd(world: "World") -> OceanCFDState:
    """Snapshots the world's current elevation/wind/temperature (via climate.py's own public
    pipeline -- reusing its grid construction and elevation resampling rather than
    duplicating them). Wind and the starting current/sea-surface state each independently
    resume from World.remembered_wind_u/World.remembered_ocean_u (etc., see their own
    docstrings) when set -- wind from a prior "atmosphere_cfd" session (or this same mode's
    own prior session, which never changes it -- see module docstring), current/eta/
    temperature/sediment from this mode's own prior session -- falling back to climate.py's
    fresh diagnostic wind and an ocean at rest (u = v = eta = 0, no sediment in suspension)
    for whichever of those has nothing to resume from."""
    height, width = climate.grid_dimensions(world.climate_density)
    fields = climate.compute_climate(world, height, width)

    depth_m = np.where(fields.is_ocean, np.clip(-fields.elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0)
    zeros = np.zeros((height, width))
    has_remembered_current = world.remembered_ocean_u is not None

    return OceanCFDState(
        lat_deg=fields.lat_deg,
        lon_deg=fields.lon_deg,
        world_xyz=fields.world_xyz,
        is_ocean=fields.is_ocean,
        elevation_m=fields.elevation_m,
        depth_m=depth_m,
        u=world.remembered_ocean_u.copy() if has_remembered_current else zeros.copy(),
        v=world.remembered_ocean_v.copy() if has_remembered_current else zeros.copy(),
        eta=world.remembered_ocean_eta.copy() if has_remembered_current else zeros.copy(),
        temperature_c=world.remembered_ocean_temperature_c.copy() if has_remembered_current else fields.ocean_temperature_c.copy(),
        baseline_temperature_c=fields.ocean_temperature_c.copy(),
        sediment_concentration=world.remembered_ocean_sediment_concentration.copy() if has_remembered_current else zeros.copy(),
        sediment_deposited_m=world.remembered_ocean_sediment_deposited_m.copy() if has_remembered_current else zeros.copy(),
        wind_u=fields.wind_u.copy() if world.remembered_wind_u is None else world.remembered_wind_u.copy(),
        wind_v=fields.wind_v.copy() if world.remembered_wind_v is None else world.remembered_wind_v.copy(),
    )


def remember_ocean_state(world: "World", state: OceanCFDState) -> None:
    """Snapshots this session's final current/sea-surface state onto `world`'s
    remembered_ocean_* fields (see World.remembered_ocean_u's own docstring) so a later
    switch back into "ocean_cfd" can resume from it instead of starting the ocean at rest
    again. Leaves world.remembered_wind_u/v untouched -- this mode never changes wind (see
    module docstring), so there's nothing new to remember there."""
    world.remembered_ocean_u = state.u.copy()
    world.remembered_ocean_v = state.v.copy()
    world.remembered_ocean_eta = state.eta.copy()
    world.remembered_ocean_temperature_c = state.temperature_c.copy()
    world.remembered_ocean_sediment_concentration = state.sediment_concentration.copy()
    world.remembered_ocean_sediment_deposited_m = state.sediment_deposited_m.copy()


def step_ocean_cfd(world: "World", state: OceanCFDState, seconds: float) -> None:
    """Advances `state` by `seconds` of real time, in as many CFL-stable substeps as the
    current grid/reduced-gravity wave speed demand (see fluid_dynamics.cfl_substeps). `world`
    is accepted (unused directly) to match step_world's own `(world, ...)` calling
    convention -- main.py dispatches on `world.fluid_mode` without needing to know which
    signature each mode's step function has."""
    del world  # not otherwise needed -- state is fully self-contained once initialized
    height, width = state.is_ocean.shape
    dx_m, dy_m = fluid_dynamics.grid_spacing_m(state.lat_deg, height, width)

    wave_speed = float(np.sqrt(REDUCED_GRAVITY_M_S2 * max(state.depth_m.max(), MIN_DEPTH_M)))
    coastal = fluid_dynamics.coastal_ocean_mask(state.is_ocean)

    min_spacing_m = fluid_dynamics.stable_min_spacing_m(dx_m, dy_m)
    current_speed = float(np.hypot(state.u, state.v).max(initial=0.0))
    n_substeps, dt_s = fluid_dynamics.cfl_substeps(seconds, min_spacing_m, wave_speed, current_speed, MAX_SUBSTEPS_PER_STEP)
    drag = BOTTOM_DRAG_PER_S + fluid_dynamics.polar_sponge_drag_per_s(state.lat_deg, POLAR_SPONGE_MAX_DRAG_PER_S)
    # state.lat_deg/width are fixed for the whole session, so this (unlike u/v/dt_s below) is
    # the same every substep -- see fluid_dynamics.advection_geometry's own docstring.
    advect_geom = fluid_dynamics.advection_geometry(state.lat_deg, width)
    f = fluid_dynamics.coriolis_parameter(state.lat_deg)[:, None]

    u, v, eta = state.u, state.v, state.eta
    temperature_c = state.temperature_c
    sediment_concentration = state.sediment_concentration
    sediment_deposited_m = state.sediment_deposited_m

    for _ in range(n_substeps):
        deta_dx, deta_dy = fluid_dynamics.gradient_m(eta, dx_m, dy_m)

        wind_speed = np.hypot(state.wind_u, state.wind_v)
        tau_u = (WIND_STRESS_COEFFICIENT * wind_speed * state.wind_u) / MIXED_LAYER_DEPTH_M
        tau_v = (WIND_STRESS_COEFFICIENT * wind_speed * state.wind_v) / MIXED_LAYER_DEPTH_M

        # Drag (bottom drag + the polar sponge, see `drag` above) is applied semi-implicitly
        # (backward Euler for just this term: divide by (1 + dt*drag) rather than subtract
        # dt*drag*u outright) -- unconditionally stable regardless of how large dt*drag gets,
        # unlike every other term here. This matters specifically for the polar sponge, whose
        # whole point is a strong damping rate right where dt is least able to resolve it
        # (see fluid_dynamics.polar_sponge_drag_per_s's own docstring); an explicit update
        # would need dt*drag comfortably under 1 to stay stable, confirmed directly during
        # development to blow up (not merely stay too energetic) once the sponge was strong
        # enough to actually fix the polar over-acceleration this exists to prevent.
        du_dt = f * v - REDUCED_GRAVITY_M_S2 * deta_dx + tau_u + VISCOSITY_M2_S * fluid_dynamics.laplacian_m(u, dx_m, dy_m)
        dv_dt = -f * u - REDUCED_GRAVITY_M_S2 * deta_dy + tau_v + VISCOSITY_M2_S * fluid_dynamics.laplacian_m(v, dx_m, dy_m)

        u = (u + dt_s * du_dt) / (1.0 + dt_s * drag)
        v = (v + dt_s * dv_dt) / (1.0 + dt_s * drag)
        u = np.where(state.is_ocean, fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(u), state.lat_deg), 0.0)
        v = np.where(state.is_ocean, fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(v), state.lat_deg), 0.0)

        # Forward-backward scheme: continuity uses the *just-updated* velocity, a standard
        # semi-implicit trick for shallow water that's meaningfully more stable than plain
        # forward Euler (updating eta from the *old* velocity) at the same substep size.
        flux_divergence = fluid_dynamics.divergence_m(state.depth_m * u, state.depth_m * v, dx_m, dy_m)
        eta = eta - dt_s * flux_divergence
        eta = np.where(state.is_ocean, fluid_dynamics.polar_zonal_filter(fluid_dynamics.grid_noise_filter(eta), state.lat_deg), 0.0)

        temperature_c = fluid_dynamics.semi_lagrangian_advect(temperature_c, u, v, dt_s, advect_geom)
        temperature_c = temperature_c + dt_s * TEMPERATURE_DIFFUSIVITY_M2_S * fluid_dynamics.laplacian_m(temperature_c, dx_m, dy_m)
        temperature_c = temperature_c + dt_s * TEMPERATURE_RELAXATION_PER_S * (state.baseline_temperature_c - temperature_c)

        speed = np.hypot(u, v)
        pickup = np.where(coastal, SEDIMENT_PICKUP_COEFFICIENT * speed, 0.0)
        settling_factor = np.clip(1.0 - speed / SEDIMENT_SETTLING_SPEED_THRESHOLD_M_S, 0.0, 1.0)
        settle = SEDIMENT_SETTLING_RATE_PER_S * sediment_concentration * settling_factor

        sediment_concentration = fluid_dynamics.semi_lagrangian_advect(sediment_concentration, u, v, dt_s, advect_geom)
        sediment_concentration = sediment_concentration + dt_s * SEDIMENT_DIFFUSIVITY_M2_S * fluid_dynamics.laplacian_m(sediment_concentration, dx_m, dy_m)
        sediment_concentration = np.clip(sediment_concentration + dt_s * (pickup - settle), 0.0, None)
        sediment_concentration = np.where(state.is_ocean, sediment_concentration, 0.0)

        # settle is already zero on land (sediment_concentration is kept at zero there, see
        # above), so sediment_deposited_m only ever grows over ocean cells with no separate
        # land mask needed here.
        sediment_deposited_m = sediment_deposited_m + dt_s * settle * SEDIMENT_DEPOSIT_DEPTH_COEFFICIENT

        state.elapsed_seconds += dt_s

    state.u, state.v, state.eta = u, v, eta
    state.temperature_c = temperature_c
    state.sediment_concentration = sediment_concentration
    state.sediment_deposited_m = sediment_deposited_m
