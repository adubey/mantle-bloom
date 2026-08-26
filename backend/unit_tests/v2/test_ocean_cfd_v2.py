import numpy as np

from app.v2 import ocean_cfd_v2
from app.v2.world_v2 import generate_world_v2

# Coarse settings throughout -- same rationale as test_world_v2_smoke.py's own
# _COARSE_KWARGS: these are regression/boundedness checks, not physics-precision tests, so
# a small fast-to-run grid is preferable.
_COARSE_KWARGS = dict(node_density=0.5, climate_density=0.5, fluid_density=0.5, num_plates=6)


def _world(seed=1):
    return generate_world_v2(seed=seed, **_COARSE_KWARGS)


def test_init_ocean_cfd_v2_produces_correctly_shaped_state_at_rest():
    world = _world()
    state = world.ocean_cfd_state
    npix = state.grid.npix
    assert state.u.shape == (npix,)
    assert state.v.shape == (npix,)
    assert state.eta.shape == (npix,)
    # At rest -- nothing has stepped yet.
    assert np.all(state.u == 0.0)
    assert np.all(state.v == 0.0)
    assert np.all(state.eta == 0.0)
    assert np.all(state.sediment_concentration == 0.0)
    assert np.all(state.sediment_deposited_m == 0.0)
    assert np.all(np.isfinite(state.temperature_c))
    assert np.all(np.isfinite(state.depth_m))
    # depth_m is zero on land -- the shallow-water layer only exists over ocean.
    assert np.all(state.depth_m[~state.is_ocean] == 0.0)


def test_step_ocean_cfd_v2_keeps_land_cells_at_zero_velocity():
    world = _world()
    state = world.ocean_cfd_state
    ocean_cfd_v2.step_ocean_cfd(world, state, seconds=3600.0)
    assert np.all(state.u[~state.is_ocean] == 0.0)
    assert np.all(state.v[~state.is_ocean] == 0.0)
    assert np.all(state.eta[~state.is_ocean] == 0.0)


def test_resample_scalar_to_equirect_matches_resample_uv_to_equirect_for_a_component():
    # resample_uv_to_equirect is defined in terms of resample_scalar_to_equirect (see its own
    # docstring) -- calling the scalar seam directly on `state.u` should give exactly the same
    # result as the first element of the (u, v) pair. Stepped first so u isn't trivially all
    # zero (see test_init_ocean_cfd_v2_produces_correctly_shaped_state_at_rest).
    world = _world()
    state = world.ocean_cfd_state
    ocean_cfd_v2.step_ocean_cfd(world, state, seconds=3600.0)
    height, width = 30, 60
    u_only = state.resample_scalar_to_equirect(state.u, height, width)
    u_from_pair, _ = state.resample_uv_to_equirect(height, width)
    assert np.array_equal(u_only, u_from_pair)


def test_step_ocean_cfd_v2_produces_no_nan_or_inf():
    world = _world()
    state = world.ocean_cfd_state
    for _ in range(3):
        ocean_cfd_v2.step_ocean_cfd(world, state, seconds=3600.0)
        for field in (state.u, state.v, state.eta, state.temperature_c, state.sediment_concentration, state.sediment_deposited_m):
            assert np.all(np.isfinite(field))


def test_step_ocean_cfd_v2_advances_elapsed_seconds_by_exactly_the_requested_amount():
    world = _world()
    state = world.ocean_cfd_state
    ocean_cfd_v2.step_ocean_cfd(world, state, seconds=7200.0)
    assert state.elapsed_seconds == 7200.0
    ocean_cfd_v2.step_ocean_cfd(world, state, seconds=3600.0)
    assert state.elapsed_seconds == 10800.0


def test_sediment_deposited_m_only_grows_and_stays_nonnegative():
    world = _world()
    state = world.ocean_cfd_state
    previous = state.sediment_deposited_m.copy()
    for _ in range(5):
        ocean_cfd_v2.step_ocean_cfd(world, state, seconds=86400.0)
        assert np.all(state.sediment_deposited_m >= previous - 1e-12)  # monotonically non-decreasing
        assert np.all(state.sediment_deposited_m >= 0.0)
        previous = state.sediment_deposited_m.copy()
    # Some real settling should have accumulated somewhere over 5 simulated days.
    assert state.sediment_deposited_m.max() > 0.0


def test_ocean_cfd_v2_never_mutates_world_plates():
    world = _world()
    elevations_before = [line.elevation.copy() for plate in world.plates for line in plate.lines]
    state = world.ocean_cfd_state
    for _ in range(3):
        ocean_cfd_v2.step_ocean_cfd(world, state, seconds=86400.0)
    elevations_after = [line.elevation for plate in world.plates for line in plate.lines]
    assert len(elevations_before) == len(elevations_after)
    for before, after in zip(elevations_before, elevations_after):
        assert np.array_equal(before, after)


def test_step_ocean_cfd_v2_stays_bounded_over_many_steps():
    # A coarse "doesn't blow up" check -- guards against the fast/slow subcycled split
    # (ocean_cfd_v2's own docstring) diverging, not an exact physical value.
    world = generate_world_v2(seed=1, node_density=1.0, climate_density=1.0, fluid_density=1.0, num_plates=6)
    state = world.ocean_cfd_state
    for _ in range(10):
        ocean_cfd_v2.step_ocean_cfd(world, state, seconds=86400.0)
    speed = np.hypot(state.u, state.v)
    assert np.all(np.isfinite(speed))
    assert speed.max() < 10.0  # a real wind-driven surface current is well under this
