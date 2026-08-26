import numpy as np
from app import atmosphere_cfd
from app.world import generate_world, step_world


def _world(seed=1, num_plates=8, climate_density=0.5, steps=1, years=5_000_000):
    world = generate_world(seed, num_plates=num_plates, continental_fraction=0.6, land_fraction=0.35, climate_density=climate_density)
    for _ in range(steps):
        step_world(world, years)
    return world


def test_init_atmosphere_cfd_produces_correctly_shaped_finite_state():
    world = _world()
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    height, width = state.elevation_m.shape
    assert state.u.shape == (height, width)
    assert state.v.shape == (height, width)
    for field in (state.u, state.v, state.eta, state.temperature_c, state.equilibrium_temperature_c, state.humidity):
        assert np.all(np.isfinite(field))
    assert np.all(state.eta == 0.0)  # geopotential anomaly starts flat


def test_step_atmosphere_cfd_produces_no_nan_or_inf():
    world = _world()
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    for _ in range(3):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
        for field in (state.u, state.v, state.eta, state.temperature_c, state.humidity, state.precipitation_mm):
            assert np.all(np.isfinite(field))


def test_step_atmosphere_cfd_advances_elapsed_seconds_by_exactly_the_requested_amount():
    world = _world()
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0)
    assert state.elapsed_seconds == 3600.0
    atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=1800.0)
    assert state.elapsed_seconds == 5400.0


def test_humidity_never_goes_negative():
    world = _world()
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    for _ in range(5):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
        assert np.all(state.humidity >= 0.0)
        assert np.all(state.precipitation_mm >= 0.0)


def test_atmosphere_cfd_never_mutates_world_plates():
    world = _world()
    elevations_before = [line.elevation.copy() for plate in world.plates for line in plate.lines]
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    for _ in range(3):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
    elevations_after = [line.elevation for plate in world.plates for line in plate.lines]
    for before, after in zip(elevations_before, elevations_after):
        assert np.array_equal(before, after)


def test_step_atmosphere_cfd_stays_bounded_over_many_steps():
    # Same "doesn't blow up" regression guard as ocean_cfd's own version -- an earlier build's
    # mountain-deflection term (applied as a direct state overwrite every substep instead of a
    # bounded tendency) compounded geometrically to absurd wind speeds within a few simulated
    # days; this catches that class of bug, not an exact physical value.
    world = _world(climate_density=1.0)
    state = atmosphere_cfd.init_atmosphere_cfd(world)
    for _ in range(10):
        atmosphere_cfd.step_atmosphere_cfd(world, state, seconds=3600.0 * 6)
    speed = np.hypot(state.u, state.v)
    assert np.all(np.isfinite(speed))
    assert speed.max() < 200.0  # a real, if strong, large-scale wind speed ceiling
