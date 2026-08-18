import numpy as np

from app import erosion, geometry, plates
from app.boundary import MAX_ELEVATION_M, MIN_ELEVATION_M
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
    row, col = erosion._climate_grid_indices(world_xyz, height=90, width=180)
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
    slope, drop_m = erosion._compute_slope(points, elevation)
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

    slope, drop_m = erosion._compute_slope(points, elevation)

    expected_drop = elevation[0] - elevation[4]
    expected_run_m = geometry.angular_distance(points[0], points[4]) * plates.PLANET_RADIUS_KM * 1000.0
    assert np.isclose(drop_m[0], expected_drop)
    assert np.isclose(slope[0], expected_drop / expected_run_m)


def test_apply_erosion_never_erodes_ocean_nodes_directly():
    # Erosion itself (rain/river/weathering) is still zeroed over ocean nodes -- but unlike
    # before river deposition existed, an ocean node's elevation *can* now rise, when a
    # river empties sediment at its mouth (a delta) -- see test_deposition_can_raise_a
    # downstream node below for a synthetic case proving that actually happens.
    world = generate_world(seed=20, num_plates=8)
    _, elevation_before, _, _ = erosion._gather_nodes(world)
    is_ocean = elevation_before <= 0.0

    erosion.apply_erosion(world, years=5_000_000)

    _, elevation_after, _, _ = erosion._gather_nodes(world)
    assert len(elevation_after) == len(elevation_before)
    assert np.all(elevation_after[is_ocean] >= elevation_before[is_ocean] - 1e-6)


def test_apply_erosion_keeps_elevation_finite_and_changing():
    world = generate_world(seed=21, num_plates=8)
    _, elevation_before, _, _ = erosion._gather_nodes(world)
    erosion.apply_erosion(world, years=5_000_000)
    _, elevation_after, _, _ = erosion._gather_nodes(world)
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
