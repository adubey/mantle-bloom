import numpy as np

from app import atmosphere_cfd
from app.world import generate_world

# Coarse settings throughout -- same rationale as test_world_smoke.py's own
# _COARSE_KWARGS: these are regression/boundedness checks, not physics-precision tests, so
# a small fast-to-run grid is preferable.
_COARSE_KWARGS = dict(node_density=0.5, climate_density=0.5, fluid_density=0.5, num_plates=6)


def _world(seed=1):
    return generate_world(seed=seed, **_COARSE_KWARGS)


def test_init_atmosphere_cfd_produces_correctly_shaped_finite_state():
    world = _world()
    state = world.atmosphere_cfd_state
    npix = state.grid.npix
    assert state.u.shape == (npix,)
    assert state.v.shape == (npix,)
    for field in (state.u, state.v, state.eta, state.temperature_c, state.equilibrium_temperature_c, state.humidity):
        assert np.all(np.isfinite(field))
    assert np.all(state.eta == 0.0)  # geopotential anomaly starts flat


def test_step_atmosphere_cfd_produces_no_nan_or_inf():
    world = _world()
    state = world.atmosphere_cfd_state
    for _ in range(3):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
        for field in (state.u, state.v, state.eta, state.temperature_c, state.humidity, state.precipitation_mm):
            assert np.all(np.isfinite(field))


def test_resample_scalar_to_equirect_matches_resample_uv_to_equirect_for_a_component():
    # resample_uv_to_equirect is defined in terms of resample_scalar_to_equirect (see its own
    # docstring) -- calling the scalar seam directly on `state.u` should give exactly the same
    # result as the first element of the (u, v) pair. Stepped first so u isn't trivially all
    # zero.
    world = _world()
    state = world.atmosphere_cfd_state
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
    height, width = 30, 60
    u_only = state.resample_scalar_to_equirect(state.u, height, width)
    u_from_pair, _ = state.resample_uv_to_equirect(height, width)
    assert np.array_equal(u_only, u_from_pair)


def test_step_atmosphere_cfd_advances_elapsed_seconds_by_exactly_the_requested_amount():
    world = _world()
    state = world.atmosphere_cfd_state
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0)
    assert state.elapsed_seconds == 3600.0
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=1800.0)
    assert state.elapsed_seconds == 5400.0


def test_humidity_never_goes_negative():
    world = _world()
    state = world.atmosphere_cfd_state
    for _ in range(5):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
        assert np.all(state.humidity >= 0.0)
        assert np.all(state.precipitation_mm >= 0.0)


def test_atmosphere_cfd_never_mutates_world_plates():
    world = _world()
    elevations_before = [line.elevation.copy() for plate in world.plates for line in plate.lines]
    state = world.atmosphere_cfd_state
    for _ in range(3):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
    elevations_after = [line.elevation for plate in world.plates for line in plate.lines]
    assert len(elevations_before) == len(elevations_after)
    for before, after in zip(elevations_before, elevations_after):
        assert np.array_equal(before, after)


def test_step_atmosphere_cfd_stays_bounded_over_many_steps():
    # Same "doesn't blow up" regression guard as ocean_cfd's own version -- guards against
    # the fast/slow subcycled split (atmosphere_cfd's own docstring) diverging, not an
    # exact physical value.
    world = generate_world(seed=1, node_density=1.0, climate_density=1.0, fluid_density=1.0, num_plates=6)
    state = world.atmosphere_cfd_state
    for _ in range(10):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
    speed = np.hypot(state.u, state.v)
    assert np.all(np.isfinite(speed))
    assert speed.max() < 200.0  # a real, if strong, large-scale wind speed ceiling
