import numpy as np

from app.rtree_index import RTree


def _brute_force_box(points, min_xy, max_xy):
    return set(np.nonzero(np.all((points >= min_xy) & (points <= max_xy), axis=1))[0].tolist())


def _brute_force_nearest(points, query_xy, exclude_index=None):
    dist = np.hypot(*(points - query_xy).T)
    if exclude_index is not None:
        dist = dist.copy()
        dist[exclude_index] = np.inf
    idx = int(np.argmin(dist))
    return idx, float(dist[idx])


def test_query_box_matches_brute_force_on_random_clouds():
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = int(rng.integers(0, 400))
        points = rng.uniform(-10.0, 10.0, size=(n, 2))
        tree = RTree.build(points)
        for _ in range(5):
            lo = rng.uniform(-10.0, 10.0, size=2)
            hi = lo + rng.uniform(0.0, 6.0, size=2)
            assert set(tree.query_box(lo, hi).tolist()) == _brute_force_box(points, lo, hi)


def test_nearest_one_matches_brute_force_on_random_clouds():
    rng = np.random.default_rng(1)
    for _ in range(20):
        n = int(rng.integers(1, 400))
        points = rng.uniform(-10.0, 10.0, size=(n, 2))
        tree = RTree.build(points)
        for _ in range(5):
            q = rng.uniform(-10.0, 10.0, size=2)
            got_idx, got_dist = tree.nearest_one(q)
            brute_idx, brute_dist = _brute_force_nearest(points, q)
            assert np.isclose(got_dist, brute_dist, atol=1e-9)
            # Distinct points can legitimately tie for nearest -- only the distance (and that
            # the reported index actually achieves it) is guaranteed, not which tied index
            # wins.
            assert np.isclose(np.hypot(*(points[got_idx] - q)), brute_dist, atol=1e-9)


def test_nearest_one_excludes_given_index():
    rng = np.random.default_rng(2)
    points = rng.uniform(-10.0, 10.0, size=(50, 2))
    tree = RTree.build(points)
    for i in range(len(points)):
        got_idx, got_dist = tree.nearest_one(points[i], exclude_index=i)
        assert got_idx != i
        brute_idx, brute_dist = _brute_force_nearest(points, points[i], exclude_index=i)
        assert np.isclose(got_dist, brute_dist, atol=1e-9)


def test_empty_tree_returns_nothing():
    tree = RTree.build(np.zeros((0, 2)))
    assert tree.nearest_one(np.array([0.0, 0.0])) is None
    assert len(tree.query_box(np.array([-1.0, -1.0]), np.array([1.0, 1.0]))) == 0
    assert tree.count_in_box(np.array([-1.0, -1.0]), np.array([1.0, 1.0])) == 0


def test_single_point_tree():
    tree = RTree.build(np.array([[3.0, -2.0]]))
    assert tree.nearest_one(np.array([0.0, 0.0])) == (0, np.hypot(3.0, 2.0))
    assert tree.query_box(np.array([2.0, -3.0]), np.array([4.0, -1.0])).tolist() == [0]
    assert tree.query_box(np.array([10.0, 10.0]), np.array([11.0, 11.0])).tolist() == []


def test_count_in_box_matches_query_box_length():
    rng = np.random.default_rng(3)
    points = rng.uniform(-5.0, 5.0, size=(200, 2))
    tree = RTree.build(points)
    lo, hi = np.array([-2.0, -2.0]), np.array([2.0, 2.0])
    assert tree.count_in_box(lo, hi) == len(tree.query_box(lo, hi))


def test_points_beyond_max_entries_per_node_build_a_multi_level_tree():
    # Forces at least a couple of internal levels (MAX_ENTRIES_PER_NODE == 16) -- exercises
    # the parent-level STR packing, not just a single flat leaf.
    rng = np.random.default_rng(4)
    points = rng.uniform(-100.0, 100.0, size=(2000, 2))
    tree = RTree.build(points)
    for _ in range(10):
        q = rng.uniform(-100.0, 100.0, size=2)
        got_idx, got_dist = tree.nearest_one(q)
        brute_idx, brute_dist = _brute_force_nearest(points, q)
        assert np.isclose(got_dist, brute_dist, atol=1e-9)
