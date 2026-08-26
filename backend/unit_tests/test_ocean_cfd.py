import numpy as np
from app import ocean_cfd
from app.world import generate_world


def _world(seed=1, num_plates=8, climate_density=0.5, steps=1, years=5_000_000):
    world = generate_world(seed, num_plates=num_plates, continental_fraction=0.6, land_fraction=0.35, climate_density=climate_density)
    from app.world import step_world

    for _ in range(steps):
        step_world(world, years)
    return world


def test_init_ocean_cfd_produces_correctly_shaped_state_at_rest():
    world = _world()
    state = ocean_cfd.init_ocean_cfd(world)
    height, width = state.is_ocean.shape
    assert state.u.shape == (height, width)
    assert state.v.shape == (height, width)
    assert state.eta.shape == (height, width)
    # At rest -- nothing has been forced yet.
    assert np.all(state.u == 0.0)
    assert np.all(state.v == 0.0)
    assert np.all(state.eta == 0.0)
    assert np.all(state.sediment_concentration == 0.0)
    assert np.all(state.sediment_deposited_m == 0.0)
    assert np.all(np.isfinite(state.temperature_c))
    assert np.all(np.isfinite(state.depth_m))
    # depth_m is zero on land -- the shallow-water layer only exists over ocean.
    assert np.all(state.depth_m[~state.is_ocean] == 0.0)


def test_step_ocean_cfd_keeps_land_cells_at_zero_velocity():
    world = _world()
    state = ocean_cfd.init_ocean_cfd(world)
    ocean_cfd.step_ocean_cfd(world, state, seconds=3600.0)
    assert np.all(state.u[~state.is_ocean] == 0.0)
    assert np.all(state.v[~state.is_ocean] == 0.0)
    assert np.all(state.eta[~state.is_ocean] == 0.0)


def test_step_ocean_cfd_produces_no_nan_or_inf():
    world = _world()
    state = ocean_cfd.init_ocean_cfd(world)
    for _ in range(3):
        ocean_cfd.step_ocean_cfd(world, state, seconds=3600.0)
        for field in (state.u, state.v, state.eta, state.temperature_c, state.sediment_concentration, state.sediment_deposited_m):
            assert np.all(np.isfinite(field))


def test_step_ocean_cfd_advances_elapsed_seconds_by_exactly_the_requested_amount():
    world = _world()
    state = ocean_cfd.init_ocean_cfd(world)
    ocean_cfd.step_ocean_cfd(world, state, seconds=7200.0)
    assert state.elapsed_seconds == 7200.0
    ocean_cfd.step_ocean_cfd(world, state, seconds=3600.0)
    assert state.elapsed_seconds == 10800.0


def test_sediment_deposited_m_only_grows_and_stays_nonnegative():
    world = _world()
    state = ocean_cfd.init_ocean_cfd(world)
    previous = state.sediment_deposited_m.copy()
    for _ in range(5):
        ocean_cfd.step_ocean_cfd(world, state, seconds=86400.0)
        assert np.all(state.sediment_deposited_m >= previous - 1e-12)  # monotonically non-decreasing
        assert np.all(state.sediment_deposited_m >= 0.0)
        previous = state.sediment_deposited_m.copy()
    # Some real settling should have accumulated somewhere over 5 simulated days.
    assert state.sediment_deposited_m.max() > 0.0


def test_ocean_cfd_never_mutates_world_plates():
    world = _world()
    elevations_before = [line.elevation.copy() for plate in world.plates for line in plate.lines]
    state = ocean_cfd.init_ocean_cfd(world)
    for _ in range(3):
        ocean_cfd.step_ocean_cfd(world, state, seconds=86400.0)
    elevations_after = [line.elevation for plate in world.plates for line in plate.lines]
    assert len(elevations_before) == len(elevations_after)
    for before, after in zip(elevations_before, elevations_after):
        assert np.array_equal(before, after)


def test_step_ocean_cfd_stays_bounded_over_many_steps():
    # A coarse "doesn't blow up" check -- confirmed directly during development that an
    # earlier build's numerics diverged to absurd speeds within the first simulated day at
    # this same (finer) grid; this guards against that class of regression, not an exact
    # physical value.
    world = _world(climate_density=1.0)
    state = ocean_cfd.init_ocean_cfd(world)
    for _ in range(10):
        ocean_cfd.step_ocean_cfd(world, state, seconds=86400.0)
    speed = np.hypot(state.u, state.v)
    assert np.all(np.isfinite(speed))
    assert speed.max() < 10.0  # a real wind-driven surface current is well under this
