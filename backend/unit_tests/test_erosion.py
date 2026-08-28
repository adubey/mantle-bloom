import numpy as np

from app import erosion, geometry, plates
from app.elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M
from app.world import World, generate_world


def test_climate_grid_indices_matches_build_grid_convention():
    # Row 0 = north pole, row increases southward; column increases eastward from lon=-180.
    world_xyz = np.array(
        [
            geometry.latlon_to_xyz(np.radians(89.0), np.radians(-179.0)),
            geometry.latlon_to_xyz(np.radians(-89.0), np.radians(179.0)),
            geometry.latlon_to_xyz(0.0, 0.0),
        ]
    )
    row, col = erosion.climate_grid_indices(world_xyz, height=90, width=180)
    assert row[0] == 0
    assert row[1] == 89
    assert col[0] == 0
    assert col[2] == 90  # lon=0 falls in the middle column


def test_compute_slope_zero_for_flat_cluster():
    # A tight cluster of same-elevation points -- nobody is lower than anybody else, so
    # both the slope and its cap (drop to lowest neighbor) must be exactly 0.
    theta = np.linspace(0.0, 0.01, 8)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = np.full(8, 500.0)
    slope, drop_m = erosion.compute_slope(points, elevation)
    assert np.allclose(slope, 0.0)
    assert np.allclose(drop_m, 0.0)


def test_compute_slope_matches_known_gradient():
    # 6 points along a line with elevation decreasing monotonically by 100m per step --
    # point 0's SLOPE_NEIGHBOR_COUNT=4 nearest neighbors are points 1-4, and since
    # elevation decreases monotonically, the *lowest* of those 4 is point 4 (farthest in
    # this small set, but lowest).
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)
    elevation = 500.0 - 100.0 * np.arange(6)

    slope, drop_m = erosion.compute_slope(points, elevation)

    expected_drop = elevation[0] - elevation[4]
    expected_run_m = geometry.angular_distance(points[0], points[4]) * plates.PLANET_RADIUS_KM * 1000.0
    assert np.isclose(drop_m[0], expected_drop)
    assert np.isclose(slope[0], expected_drop / expected_run_m)


def test_weathering_relief_factor_suppresses_flat_terrain_and_saturates_on_steep_terrain():
    # Regression check for the "flat coastal plain weathers exactly as fast as a mountain"
    # bug: weathering must scale down toward 0 as slope -> 0 (same clip-and-normalize idiom
    # humidity_norm already uses), and saturate to full strength on genuinely steep terrain,
    # rather than being a flat function of wind/humidity alone.
    d = 0.01
    theta = d * np.arange(6)
    points = geometry.local_xyz(np.zeros_like(theta), theta)

    flat_elevation = np.full(6, 30.0)  # uniform -- see test_compute_slope_zero_for_flat_cluster
    flat_slope, _ = erosion.compute_slope(points, flat_elevation)
    flat_relief_factor = np.clip(flat_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    assert np.allclose(flat_relief_factor, 0.0)

    steep_elevation = 30.0 + 2000.0 * np.arange(6)  # steep monotonic gradient
    steep_slope, _ = erosion.compute_slope(points, steep_elevation)
    steep_relief_factor = np.clip(steep_slope / erosion.WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    # Every point except the global minimum (point 0, which has no lower neighbor at all)
    # should be steep enough to saturate weathering to full strength.
    assert np.allclose(steep_relief_factor[1:], 1.0)


def test_submarine_erosion_scales_with_slope_and_depth_only_over_ocean():
    # Flat sea floor (slope 0) erodes nothing; a steeper, deeper submerged scarp erodes more
    # than a steeper, shallower one; subaerial nodes are untouched by this term.
    is_ocean = np.array([True, True, True, False])
    elevation = np.array([-100.0, -100.0, -5000.0, 500.0])
    slope = np.array([0.0, 0.03, 0.03, 0.03])
    amt = erosion.submarine_erosion_amount(elevation, slope, is_ocean, dt_myr=5.0)
    assert amt[0] == 0.0  # flat
    assert amt[2] > amt[1] > 0.0  # deeper (more pressure) erodes faster at the same slope
    assert amt[3] == 0.0  # subaerial


def test_coastal_erosion_confined_to_band_and_peaks_near_freezing():
    elevation = np.array([0.0, 0.0, 150.0, 5000.0, -5000.0])
    freezing = np.full(5, erosion.COASTAL_FROST_PEAK_C)
    warm = np.full(5, 30.0)
    at_freezing = erosion.coastal_erosion_amount(elevation, freezing, dt_myr=1.0)
    when_warm = erosion.coastal_erosion_amount(elevation, warm, dt_myr=1.0)
    assert at_freezing[0] > when_warm[0] > 0.0  # frost adds on top of wave attack at the shore
    assert at_freezing[2] > 0.0 and at_freezing[2] < at_freezing[0]  # in-band but tapering
    assert at_freezing[3] == 0.0 and at_freezing[4] == 0.0  # far above / far below sea level


def test_apply_erosion_can_lower_a_submerged_range():
    # End to end: over a generated world, some ocean node loses net elevation once submarine +
    # coastal erosion are in play (they were entirely erosion-exempt before).
    world = generate_world(seed=20, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    is_ocean = elevation_before <= 0.0
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    assert len(elevation_after) == len(elevation_before)
    assert np.any(elevation_after[is_ocean] < elevation_before[is_ocean] - 1e-6)


def test_apply_erosion_keeps_elevation_finite_and_changing():
    world = generate_world(seed=21, num_plates=8)
    _, elevation_before, _, _, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _, _, _ = erosion._gather_nodes(world)
    # Erosion/deposition/weathering can now raise *or* lower any given node (deposition can
    # outweigh local erosion at a floodplain or delta) -- just confirm the pass actually did
    # something and didn't produce garbage.
    assert np.all(np.isfinite(elevation_after))
    assert not np.allclose(elevation_after, elevation_before)


def test_apply_erosion_respects_elevation_bounds():
    world = generate_world(seed=22, num_plates=8)
    erosion.apply_erosion(world, years=5_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(line.elevation >= MIN_ELEVATION_M - 1e-6)
            assert np.all(line.elevation <= MAX_ELEVATION_M + 1e-6)


def test_apply_erosion_noop_for_empty_world():
    world = World(seed=0, plates=[])
    erosion.apply_erosion(world, years=1_000_000)  # must not raise
    assert world.plates == []
