"""Ocean Fluid Dynamics: a genuine time-integrated shallow-water simulation of ocean currents,
replacing climate.py's own diagnostic (recomputed-from-scratch-every-call) Ekman-current
heuristic with real prognostic state that persists and evolves continuously -- see
docs/simulation-model.md#ocean-atmospheric-fluid-dynamics for the full design rationale
(reduced gravity, substepping) shared with atmosphere_cfd.py.

**Always-on, not a mode.** `World.ocean_cfd_state` is created once, by `init_ocean_cfd`,
during `generate_world` (after `World.atmosphere_cfd_state`, which it seeds its own wind
forcing from), and never re-initialized after that -- every `/world/step` call advances it by
a fixed `SECONDS_PER_TECTONIC_STEP` (1 simulated week) via `step_ocean_cfd`, regardless of how
many tectonic years that step covers (see world.py's `step_world`), gated on
`World.simulate_climate_biomes` the same way erosion/hydrology already are. `refresh_forcing`
re-samples the terrain/wind this state reacts to (`elevation_m`/`is_ocean`/`depth_m`/
`wind_u`/`wind_v`) once per tectonics step, right before that step's `step_ocean_cfd` call --
terrain changes slowly relative to one tectonics step, but not never, and wind now comes from
the world's own continuously-evolving atmosphere_cfd_state rather than a value frozen forever
-- while leaving every genuinely prognostic field (u/v/eta/temperature_c/sediment_*) untouched,
so currents keep evolving continuously rather than resetting.

**Inputs, matching the user's own spec.** Coriolis force and land (both baked into the
momentum equation itself -- see step_ocean_cfd), wind (see `refresh_forcing` above -- the
world's own continuously-evolving atmosphere_cfd_state, not an independent simulation of its
own), and temperature (both an input, via the wind-driven baseline this state starts from, and
an output the currents themselves advect).

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

# The fixed real-time increment step_world's own _advance_fluid_dynamics advances this state
# by every tectonics step, regardless of how many tectonic years that step covers (see module
# docstring's "Always-on, not a mode") -- one simulated week (longer than atmosphere_cfd's own
# SECONDS_PER_TECTONIC_STEP -- ocean currents evolve on a slower timescale than wind).
SECONDS_PER_TECTONIC_STEP = 7 * 86400.0


@dataclass
class OceanCFDState:
    lat_deg: np.ndarray  # (H,)
    lon_deg: np.ndarray  # (W,)
    world_xyz: np.ndarray  # (H, W, 3)
    is_ocean: np.ndarray  # (H, W) bool -- refreshed once per tectonics step, see refresh_forcing
    elevation_m: np.ndarray  # (H, W) -- refreshed once per tectonics step, see refresh_forcing
    depth_m: np.ndarray  # (H, W) -- effective shallow-water layer depth, zero on land, refreshed once per tectonics step
    u: np.ndarray  # (H, W) eastward current, m/s
    v: np.ndarray  # (H, W) northward current, m/s
    eta: np.ndarray  # (H, W) sea-surface-height anomaly, m
    temperature_c: np.ndarray  # (H, W)
    baseline_temperature_c: np.ndarray  # (H, W) -- fixed relaxation target, see module docstring
    sediment_concentration: np.ndarray  # (H, W), arbitrary units >= 0
    sediment_deposited_m: np.ndarray  # (H, W), cumulative, tracking-only -- never touches world.plates
    wind_u: np.ndarray  # (H, W) eastward, refreshed once per tectonics step from World.atmosphere_cfd_state, see refresh_forcing
    wind_v: np.ndarray  # (H, W) northward, refreshed once per tectonics step from World.atmosphere_cfd_state, see refresh_forcing
    elapsed_seconds: float = 0.0

    def resample_scalar_to_equirect(self, field: np.ndarray, height: int, width: int) -> np.ndarray:
        """See atmosphere_cfd.AtmosphereCFDState's own matching method -- the same
        grid-agnostic seam `climate.compute_climate` reads temperature through."""
        return fluid_dynamics.resample_to_grid(field, height, width)

    def resample_uv_to_equirect(self, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
        """See atmosphere_cfd.AtmosphereCFDState's own matching method -- the same
        grid-agnostic seam `climate.compute_climate` reads currents through."""
        u = self.resample_scalar_to_equirect(self.u, height, width)
        v = self.resample_scalar_to_equirect(self.v, height, width)
        return u, v


def init_ocean_cfd(world: "World") -> OceanCFDState:
    """Called exactly once, by generate_world (after atmosphere_cfd.init_atmosphere_cfd, whose
    result this seeds its own wind forcing from), to seed this world's permanent
    World.ocean_cfd_state. Snapshots the world's current elevation/temperature (via climate.py's
    own public pipeline -- reusing its grid construction and elevation resampling rather than
    duplicating them); starts the ocean at rest (u = v = eta = 0, no sediment in suspension) --
    there's no real ocean-current analog to atmosphere_cfd's own diagnostic-wind bootstrap to
    start from instead. Sized by World.fluid_density, not World.climate_density -- see the
    former's own docstring for why this gets its own, independently choosable grid
    resolution."""
    height, width = climate.grid_dimensions(world.fluid_density)
    fields = climate.compute_climate(world, height, width)

    depth_m = np.where(fields.is_ocean, np.clip(-fields.elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0).astype(np.float32)
    zeros = np.zeros((height, width), dtype=np.float32)

    # Every substep-loop field below is float32 (not the rest of this codebase's usual
    # float64) -- these are the same memory-bandwidth-bound elementwise array ops
    # (np.roll-based gradients/Laplacians, the semi-Lagrangian gather, ...) every substep,
    # up to MAX_SUBSTEPS_PER_STEP times per call, so halving each array's footprint roughly
    # halves that work; profiling confirmed close to a 2x wall-clock win. Precision loss
    # (~7 vs ~15 significant digits) is a non-issue against this solver's own already-tuned-
    # not-derived constants (see module docstring). state.lat_deg is float32 too, since it
    # feeds every one of those per-substep calls (grid spacing, Coriolis, polar filtering,
    # advection geometry) -- NumPy would silently upcast the whole chain back to float64 the
    # moment a float32 field met a float64 lat_deg-derived array. lon_deg/world_xyz stay
    # float64 -- render_image.py's own projection math, not this module's substep loop.
    return OceanCFDState(
        lat_deg=fields.lat_deg.astype(np.float32),
        lon_deg=fields.lon_deg,
        world_xyz=fields.world_xyz,
        is_ocean=fields.is_ocean,
        elevation_m=fields.elevation_m,
        depth_m=depth_m,
        u=zeros.copy(),
        v=zeros.copy(),
        eta=zeros.copy(),
        temperature_c=fields.ocean_temperature_c.astype(np.float32),
        baseline_temperature_c=fields.ocean_temperature_c.astype(np.float32),
        sediment_concentration=zeros.copy(),
        sediment_deposited_m=zeros.copy(),
        wind_u=world.atmosphere_cfd_state.u.copy(),
        wind_v=world.atmosphere_cfd_state.v.copy(),
    )


def refresh_forcing(world: "World", state: OceanCFDState, terrain: climate.ClimateFields) -> None:
    """Re-samples the terrain (`elevation_m`/`is_ocean`/`depth_m`) and wind forcing
    (`wind_u`/`wind_v`, from the world's own continuously-evolving
    `World.atmosphere_cfd_state` -- both already at the same `World.fluid_density` grid, no
    resample needed) this state's substep loop reads every substep, from `terrain` -- this
    tectonics step's own current climate.compute_climate snapshot, see world.py's
    `_advance_fluid_dynamics` -- while leaving every genuinely prognostic field
    (u/v/eta/temperature_c/sediment_*) untouched, so currents/temperature/sediment keep
    evolving continuously across tectonics steps rather than resetting. Must be called after
    `atmosphere_cfd.step_atmosphere_cfd` has already advanced this step's wind (see
    `_advance_fluid_dynamics`'s own call order), so `wind_u`/`wind_v` here reflect *this*
    step's wind, not last step's."""
    state.elevation_m = terrain.elevation_m
    state.is_ocean = terrain.is_ocean
    state.depth_m = np.where(terrain.is_ocean, np.clip(-terrain.elevation_m, MIN_DEPTH_M, MAX_DEPTH_M), 0.0).astype(np.float32)
    state.wind_u = world.atmosphere_cfd_state.u.copy()
    state.wind_v = world.atmosphere_cfd_state.v.copy()


def step_ocean_cfd(world: "World", state: OceanCFDState, seconds: float) -> None:
    """Advances `state` by `seconds` of real time, in as many CFL-stable substeps as the
    current grid/reduced-gravity wave speed demand (see fluid_dynamics.cfl_substeps). Called
    by world.py's `_advance_fluid_dynamics` with `seconds=SECONDS_PER_TECTONIC_STEP`, right
    after `refresh_forcing`. `world` is accepted (unused directly) to match
    atmosphere_cfd.step_atmosphere_cfd's own `(world, ...)` calling convention."""
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
