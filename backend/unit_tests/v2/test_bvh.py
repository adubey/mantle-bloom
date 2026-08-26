import numpy as np

from app import geometry
from app.v2 import bvh


def _random_unit_points(rng, n):
    pts = rng.normal(size=(n, 3))
    return pts / np.linalg.norm(pts, axis=-1, keepdims=True)


def test_query_nearest_matches_brute_force():
    rng = np.random.default_rng(0)
    points = _random_unit_points(rng, 400)
    tree = bvh.build(points, leaf_size=8)
    queries = _random_unit_points(rng, 60)

    idx, dist = bvh.query_nearest(tree, points, queries)

    brute_idx = np.array([int(np.argmin(geometry.angular_distance(points, q[None, :]))) for q in queries])
    brute_dist = np.array([float(geometry.angular_distance(points, q[None, :]).min()) for q in queries])
    assert np.array_equal(idx, brute_idx)
    assert np.allclose(dist, brute_dist, atol=1e-12)


def test_query_nearest_empty_tree():
    idx, dist = bvh.query_nearest(None, np.zeros((0, 3)), np.zeros((3, 3)))
    assert np.all(idx == -1)
    assert np.all(np.isinf(dist))


def test_query_nearest_cross_matches_brute_force():
    rng = np.random.default_rng(1)
    points_a = _random_unit_points(rng, 250)
    points_b = _random_unit_points(rng, 180)
    tree_a = bvh.build(points_a, leaf_size=8)
    tree_b = bvh.build(points_b, leaf_size=8)
    reach = 0.05

    pairs = bvh.query_nearest_cross(tree_a, points_a, tree_b, points_b, reach)
    found = {(i, j) for i, j, _ in pairs}

    d = geometry.angular_distance(points_a[:, None, :], points_b[None, :, :])
    expected = {(int(i), int(j)) for i, j in zip(*np.where(d <= reach))}
    assert found == expected


def test_bounding_cap_actually_bounds_its_members():
    rng = np.random.default_rng(2)
    points = _random_unit_points(rng, 100)
    tree = bvh.build(points, leaf_size=8)

    def check(node):
        if node is None:
            return
        if node.is_leaf:
            member_pts = points[node.indices]
            dist = geometry.angular_distance(member_pts, node.center[None, :])
            assert np.all(dist <= node.radius + 1e-9)
        else:
            check(node.left)
            check(node.right)

    check(tree)
