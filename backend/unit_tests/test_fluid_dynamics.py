import numpy as np
from app import fluid_dynamics


def test_coriolis_parameter_flips_sign_by_hemisphere():
    lat_deg = np.array([-45.0, 0.0, 45.0])
    f = fluid_dynamics.coriolis_parameter(lat_deg)
    assert f[0] < 0 < f[2]
    assert f[1] == 0.0


def test_grid_spacing_m_shrinks_toward_equator_to_pole_but_clamps_near_pole():
    lat_deg = 90.0 - (np.arange(36) + 0.5) * (180.0 / 36)
    dx_m, dy_m = fluid_dynamics.grid_spacing_m(lat_deg, 36, 72)
    # Real east-west spacing shrinks approaching the pole (cos(lat) -> 0)...
    assert dx_m[0] < dx_m[len(dx_m) // 2]
    # ...but grid_spacing_m clamps it at its own value at POLAR_FILTER_START_LAT_DEG, so it
    # never keeps shrinking arbitrarily close to zero the way a raw cos(lat) formula would --
    # see that function's own docstring for why (gradient/Laplacian 1/dx^2 blow-up otherwise).
    unclamped_pole_dx = (2.0 * np.pi * fluid_dynamics.plates.PLANET_RADIUS_KM * 1000.0 / 72) * np.cos(np.radians(lat_deg[0]))
    assert dx_m[0] > unclamped_pole_dx
    assert dy_m > 0


def test_stable_min_spacing_m_is_the_grids_own_minimum():
    lat_deg = 90.0 - (np.arange(10) + 0.5) * (180.0 / 10)
    dx_m, dy_m = fluid_dynamics.grid_spacing_m(lat_deg, 10, 20)
    assert fluid_dynamics.stable_min_spacing_m(dx_m, dy_m) == min(float(dx_m.min()), dy_m)


def test_polar_zonal_filter_blends_toward_zonal_mean_by_the_documented_ramp():
    lat_deg = 90.0 - (np.arange(18) + 0.5) * (180.0 / 18)
    width = 36
    rng = np.random.default_rng(0)
    field = rng.random((18, width))
    filtered = fluid_dynamics.polar_zonal_filter(field, lat_deg)
    # Exactly the documented linear ramp -- 0 strength (no blending at all) equatorward of
    # POLAR_FILTER_START_LAT_DEG, blending toward the zonal mean poleward of it.
    abs_lat = np.abs(lat_deg)
    span = 90.0 - fluid_dynamics.POLAR_FILTER_START_LAT_DEG
    strength = np.clip((abs_lat - fluid_dynamics.POLAR_FILTER_START_LAT_DEG) / span, 0.0, 1.0)[:, None]
    expected = field * (1.0 - strength) + field.mean(axis=1, keepdims=True) * strength
    assert np.allclose(filtered, expected)
    # The row nearest the equator (|lat| well under POLAR_FILTER_START_LAT_DEG) should be
    # untouched entirely (strength == 0 there).
    equator_row = len(lat_deg) // 2
    assert np.allclose(filtered[equator_row], field[equator_row])
    # The row nearest the pole should be blended meaningfully toward its own zonal mean.
    assert not np.allclose(filtered[0], field[0])


def test_cfl_substeps_shrinks_dt_for_a_faster_wave_speed():
    _, dt_slow = fluid_dynamics.cfl_substeps(86400.0, 10_000.0, wave_speed_m_s=1.0, max_advect_speed_m_s=0.0, max_substeps=10_000)
    _, dt_fast = fluid_dynamics.cfl_substeps(86400.0, 10_000.0, wave_speed_m_s=100.0, max_advect_speed_m_s=0.0, max_substeps=10_000)
    assert dt_fast < dt_slow


def test_cfl_substeps_respects_the_max_substeps_ceiling():
    n, dt = fluid_dynamics.cfl_substeps(86400.0, 1.0, wave_speed_m_s=1000.0, max_advect_speed_m_s=0.0, max_substeps=5)
    assert n == 5
    assert dt == 86400.0 / 5


def test_cfl_substeps_dt_always_evenly_divides_seconds():
    n, dt = fluid_dynamics.cfl_substeps(3600.0, 5000.0, wave_speed_m_s=13.0, max_advect_speed_m_s=2.0, max_substeps=2000)
    assert n * dt == 3600.0


def test_semi_lagrangian_advect_is_identity_at_zero_velocity():
    lat_deg = 90.0 - (np.arange(10) + 0.5) * (180.0 / 10)
    field = np.arange(10 * 20, dtype=float).reshape(10, 20)
    zero = np.zeros_like(field)
    geometry = fluid_dynamics.advection_geometry(lat_deg, width=20)
    advected = fluid_dynamics.semi_lagrangian_advect(field, zero, zero, dt_s=3600.0, geometry=geometry)
    assert np.array_equal(advected, field)


def _reference_gradient_m(field, dx_m, dy_m):
    """The pre-Numba formula, reimplemented independently here (not calling production code)
    so these tests catch a regression in the jitted rewrite rather than merely re-asserting
    whatever it currently computes."""
    gx = (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * dx_m[:, None])
    gy = (np.roll(field, 1, axis=0) - np.roll(field, -1, axis=0)) / (2.0 * dy_m)
    return gx, gy


def _reference_laplacian_m(field, dx_m, dy_m):
    d2x = (np.roll(field, -1, axis=1) + np.roll(field, 1, axis=1) - 2.0 * field) / (dx_m[:, None] ** 2)
    d2y = (np.roll(field, 1, axis=0) + np.roll(field, -1, axis=0) - 2.0 * field) / (dy_m**2)
    return d2x + d2y


def _reference_divergence_m(u, v, dx_m, dy_m):
    gx, _ = _reference_gradient_m(u, dx_m, dy_m)
    _, gy = _reference_gradient_m(v, dx_m, dy_m)
    return gx + gy


def _reference_grid_noise_filter(field, weight):
    neighbor_avg = 0.25 * (
        np.roll(field, -1, axis=1) + np.roll(field, 1, axis=1) + np.roll(field, -1, axis=0) + np.roll(field, 1, axis=0)
    )
    return field + weight * (neighbor_avg - field)


def _random_field_and_spacing(rng, height, width):
    field = rng.normal(size=(height, width))
    dx_m = rng.uniform(10_000.0, 100_000.0, size=height)
    dy_m = 12_000.0
    return field, dx_m, dy_m


def test_gradient_m_matches_reference_formula():
    rng = np.random.default_rng(1)
    field, dx_m, dy_m = _random_field_and_spacing(rng, 12, 24)
    gx, gy = fluid_dynamics.gradient_m(field, dx_m, dy_m)
    expected_gx, expected_gy = _reference_gradient_m(field, dx_m, dy_m)
    assert np.allclose(gx, expected_gx)
    assert np.allclose(gy, expected_gy)


def test_laplacian_m_matches_reference_formula():
    rng = np.random.default_rng(2)
    field, dx_m, dy_m = _random_field_and_spacing(rng, 12, 24)
    result = fluid_dynamics.laplacian_m(field, dx_m, dy_m)
    assert np.allclose(result, _reference_laplacian_m(field, dx_m, dy_m))


def test_divergence_m_matches_reference_formula():
    rng = np.random.default_rng(3)
    u, dx_m, dy_m = _random_field_and_spacing(rng, 12, 24)
    v, _, _ = _random_field_and_spacing(rng, 12, 24)
    result = fluid_dynamics.divergence_m(u, v, dx_m, dy_m)
    assert np.allclose(result, _reference_divergence_m(u, v, dx_m, dy_m))


def test_grid_noise_filter_matches_reference_formula():
    rng = np.random.default_rng(4)
    field, _, _ = _random_field_and_spacing(rng, 12, 24)
    result = fluid_dynamics.grid_noise_filter(field, weight=0.05)
    assert np.allclose(result, _reference_grid_noise_filter(field, 0.05))


def test_semi_lagrangian_advect_matches_reference_formula_at_nonzero_velocity():
    """test_semi_lagrangian_advect_is_identity_at_zero_velocity above only exercises the
    zero-velocity path (where every backward-trace lands back on its own cell); this exercises
    the actual round/clip/mod/gather chain the Numba kernel now performs instead of NumPy."""
    height, width = 12, 24
    lat_deg = (90.0 - (np.arange(height) + 0.5) * (180.0 / height)).astype(np.float32)
    rng = np.random.default_rng(5)
    field = rng.normal(size=(height, width)).astype(np.float32)
    u = rng.normal(0.0, 5.0, size=(height, width)).astype(np.float32)
    v = rng.normal(0.0, 5.0, size=(height, width)).astype(np.float32)
    dt_s = 3600.0

    radius_m = fluid_dynamics.plates.PLANET_RADIUS_KM * 1000.0
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_deg_row = -180.0 + (np.arange(width, dtype=np.float32) + 0.5) * (360.0 / width)
    cos_lat = np.clip(np.cos(np.radians(lat_grid)), 0.15, 1.0)
    meters_per_deg_lat = (np.pi * radius_m) / 180.0
    src_lat = lat_grid - (v * dt_s) / meters_per_deg_lat
    src_lon = lon_deg_row[None, :] - (u * dt_s) / (meters_per_deg_lat * cos_lat)
    src_row = np.clip(np.round((90.0 - src_lat) / (180.0 / height) - 0.5).astype(np.int64), 0, height - 1)
    src_col = np.round((src_lon + 180.0) / (360.0 / width) - 0.5).astype(np.int64) % width
    expected = field[src_row, src_col]

    geometry = fluid_dynamics.advection_geometry(lat_deg, width)
    result = fluid_dynamics.semi_lagrangian_advect(field, u, v, dt_s, geometry)
    assert np.array_equal(result, expected)


def test_coastal_ocean_mask_excludes_open_ocean_and_all_land():
    is_ocean = np.array(
        [
            [True, True, True, False],
            [True, True, True, False],
        ]
    )
    coastal = fluid_dynamics.coastal_ocean_mask(is_ocean)
    # Column 2 (ocean, adjacent to land at column 3) is coastal; column 0 (wraps around to
    # column 3, also land) is coastal too; column 1 (ocean on every side) is not.
    assert coastal[0, 2] and coastal[0, 0]
    assert not coastal[0, 1]
    assert not np.any(coastal[:, 3])  # land itself is never "coastal ocean"
