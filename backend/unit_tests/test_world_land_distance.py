import numpy as np
from app.world import World, generate_world, step_world

from .quad_fixtures import node_points, quad_plate


def _plate(plate_id, crust_type, columns, base_elevation):
    return quad_plate(plate_id, crust_type, columns=columns, elevation=base_elevation)


def test_distance_from_land_approx_zero_at_a_land_node_itself():
    land = _plate(0, "continental", [0, 1], base_elevation=200.0)
    world = World(seed=0, plates=[land])

    land_point = node_points(land)[0]
    dist = world.distance_from_land_approx(np.array([land_point]))
    assert np.allclose(dist, 0.0, atol=1e-9)


def test_distance_from_land_approx_matches_a_manual_kdtree():
    from scipy.spatial import cKDTree

    land = _plate(0, "continental", [0, 1, 2], base_elevation=200.0)
    submerged = _plate(1, "continental", [4, 6], base_elevation=-500.0)
    world = World(seed=0, plates=[land, submerged])

    query_points = node_points(submerged)
    dist = world.distance_from_land_approx(query_points)

    land_points = node_points(land)
    expected, _ = cKDTree(land_points).query(query_points)
    assert np.allclose(dist, expected)


def test_distance_from_land_approx_is_inf_with_no_land():
    submerged = _plate(0, "continental", [0, 5], base_elevation=-500.0)
    world = World(seed=0, plates=[submerged])

    dist = world.distance_from_land_approx(node_points(submerged))
    assert np.all(np.isinf(dist))


def test_distance_from_land_approx_empty_points_returns_empty_array():
    land = _plate(0, "continental", [0], base_elevation=200.0)
    world = World(seed=0, plates=[land])
    dist = world.distance_from_land_approx(np.zeros((0, 3)))
    assert dist.shape == (0,)


def test_distance_from_land_approx_caches_the_kdtree_until_invalidated():
    land = _plate(0, "continental", [0, 1], base_elevation=200.0)
    world = World(seed=0, plates=[land])

    query = np.array([node_points(land)[0]])
    world.distance_from_land_approx(query)
    cached = world.land_kdtree_cache
    assert cached is not None

    # Elevation changing afterward shouldn't rebuild the tree until it's explicitly reset --
    # the same "up to one step stale" tolerance World.climate_cache/hydrology_cache document.
    land.set_fields_on_plate(elevation=np.full(2, -500.0))
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
