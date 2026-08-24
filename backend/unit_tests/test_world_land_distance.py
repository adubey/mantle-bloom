import numpy as np
from app import geometry
from app.plates import ElevationLine, PlateWithLines
from app.world import World, generate_world, step_world


def _plate(plate_id, crust_type, theta, base_elevation):
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), base_elevation))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line])


def test_distance_from_land_approx_zero_at_a_land_node_itself():
    land = _plate(0, "continental", [0.0, 0.01], base_elevation=200.0)
    world = World(seed=0, plates=[land])

    land_point = land.lines[0].world_xyz(land.frame)[0]
    dist = world.distance_from_land_approx(np.array([land_point]))
    assert np.allclose(dist, 0.0, atol=1e-9)


def test_distance_from_land_approx_matches_a_manual_kdtree():
    from scipy.spatial import cKDTree

    land = _plate(0, "continental", [0.0, 0.01, 0.02], base_elevation=200.0)
    submerged = _plate(1, "continental", [0.05, 0.1], base_elevation=-500.0)
    world = World(seed=0, plates=[land, submerged])

    query_points = submerged.lines[0].world_xyz(submerged.frame)
    dist = world.distance_from_land_approx(query_points)

    land_points = land.lines[0].world_xyz(land.frame)
    expected, _ = cKDTree(land_points).query(query_points)
    assert np.allclose(dist, expected)


def test_distance_from_land_approx_is_inf_with_no_land():
    submerged = _plate(0, "continental", [0.0, 0.1], base_elevation=-500.0)
    world = World(seed=0, plates=[submerged])

    dist = world.distance_from_land_approx(submerged.lines[0].world_xyz(submerged.frame))
    assert np.all(np.isinf(dist))


def test_distance_from_land_approx_empty_points_returns_empty_array():
    land = _plate(0, "continental", [0.0], base_elevation=200.0)
    world = World(seed=0, plates=[land])
    dist = world.distance_from_land_approx(np.zeros((0, 3)))
    assert dist.shape == (0,)


def test_distance_from_land_approx_caches_the_kdtree_until_invalidated():
    land = _plate(0, "continental", [0.0, 0.01], base_elevation=200.0)
    world = World(seed=0, plates=[land])

    query = np.array([land.lines[0].world_xyz(land.frame)[0]])
    world.distance_from_land_approx(query)
    cached = world.land_kdtree_cache
    assert cached is not None

    # Elevation changing afterward shouldn't rebuild the tree until it's explicitly reset --
    # the same "up to one step stale" tolerance World.climate_cache/hydrology_cache document.
    land.replace_line(0, land.lines[0].replace(elevation=np.full(2, -500.0)))
    world.distance_from_land_approx(query)
    assert world.land_kdtree_cache is cached


def test_step_world_rebuilds_land_kdtree_cache_each_step():
    world = generate_world(seed=20, num_plates=8, continental_fraction=0.5, land_fraction=0.35)
    world.distance_from_land_approx(np.zeros((0, 3)))  # force a build this "step"
    stale_cache = world.land_kdtree_cache
    assert stale_cache is not None

    step_world(world, years=1_000_000)
    # bathymetry.py/geology.py both read distance_from_land_approx during the step, so the
    # cache is rebuilt from scratch rather than carrying over the previous step's now-stale
    # tree (see World.land_kdtree_cache's own docstring on the reset-then-lazily-rebuild
    # contract).
    assert world.land_kdtree_cache is not None
    assert world.land_kdtree_cache is not stale_cache
