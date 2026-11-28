"""`SpatialIndex` matches `cKDTree` bit for bit on both backends (#1042)."""

import copy
import os
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

from app import spatial_index
from app.spatial_index import SpatialIndex

needs_rust = pytest.mark.skipif(not spatial_index.RUST_AVAILABLE, reason="Rust extension unavailable")


def _sphere(rng, n):
    points = rng.normal(size=(n, 3))
    return points / np.linalg.norm(points, axis=1)[:, None]


@pytest.fixture(params=["rust", "scipy"])
def backend(request, monkeypatch):
    if request.param == "rust":
        if not spatial_index.RUST_AVAILABLE:
            pytest.skip("Rust extension unavailable")
    else:
        monkeypatch.setattr(spatial_index, "_KdIndex", None)
    return request.param


@pytest.fixture
def cloud():
    rng = np.random.default_rng(1042)
    # Enough points for several tree levels, and enough queries (>= 1024) to take the parallel path.
    return _sphere(rng, 3000), _sphere(rng, 2500)


def test_backend_fixture_selects_backend(backend, cloud):
    assert SpatialIndex(cloud[0]).backend == backend


@pytest.mark.parametrize("k", [1, 2, 9])
@pytest.mark.parametrize("bound", [np.inf, 0.03])
def test_query_matches_ckdtree(backend, cloud, k, bound):
    points, queries = cloud
    expected = cKDTree(points).query(queries, k=k, distance_upper_bound=bound)
    actual = SpatialIndex(points).query(queries, k=k, distance_upper_bound=bound)
    for got, want in zip(actual, expected, strict=True):
        assert got.dtype == want.dtype and got.shape == want.shape
        np.testing.assert_array_equal(got, want)


def test_query_shapes_and_scalars(backend, cloud):
    points, queries = cloud
    tree, index = cKDTree(points), SpatialIndex(points)
    d, i = index.query(queries[0])
    assert (type(d), type(i)) == (float, int)
    assert (d, i) == tree.query(queries[0])
    for x, k in ((queries[0], 3), (queries[:12].reshape(3, 4, 3), 1), (queries[:12].reshape(3, 4, 3), 2), (queries[:5], len(points) + 4)):
        for got, want in zip(index.query(x, k=k), tree.query(x, k=k), strict=True):
            np.testing.assert_array_equal(got, want)


def test_query_misses_are_inf_and_n(backend, cloud):
    points, queries = cloud
    d, i = SpatialIndex(points).query(queries, distance_upper_bound=1e-9)
    assert np.all(np.isinf(d)) and np.all(i == len(points))


def test_ball_point_matches_ckdtree(backend, cloud):
    points, queries = cloud
    tree, index = cKDTree(points), SpatialIndex(points)
    radii = np.random.default_rng(7).uniform(0.0, 0.08, len(queries))
    for r in (0.05, radii):
        expected = tree.query_ball_point(queries, r)
        actual = index.query_ball_point(queries, r)
        assert actual.dtype == object and actual.shape == expected.shape
        assert [list(row) for row in actual] == [sorted(row) for row in expected]
        indptr, indices = index.query_ball_point_csr(queries, r)
        assert indptr.dtype == indices.dtype == np.int64
        np.testing.assert_array_equal(np.diff(indptr), [len(row) for row in expected])
        np.testing.assert_array_equal(indices, np.concatenate([sorted(row) for row in expected]))
        np.testing.assert_array_equal(index.query_ball_point(queries, r, return_length=True), tree.query_ball_point(queries, r, return_length=True))


def test_ball_point_single_point_is_sorted_list(backend, cloud):
    points, queries = cloud
    expected = cKDTree(points).query_ball_point(queries[0], 0.1)
    actual = SpatialIndex(points).query_ball_point(queries[0], 0.1)
    assert isinstance(actual, list) and actual == sorted(expected)
    assert SpatialIndex(points).query_ball_point(queries[0], 0.1, return_length=True) == len(expected)


def test_ball_point_boundary_is_inclusive(backend):
    points = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    r = float(np.linalg.norm(points[0] - points[1]))
    assert SpatialIndex(points).query_ball_point(points[0], r) == cKDTree(points).query_ball_point(points[0], r) == [0, 1, 2]


def test_query_pairs_matches_ckdtree(backend, cloud):
    points, _ = cloud
    tree, index = cKDTree(points), SpatialIndex(points)
    assert index.query_pairs(0.04) == tree.query_pairs(0.04)
    pairs = tree.query_pairs(0.04, output_type="ndarray")
    np.testing.assert_array_equal(index.query_pairs(0.04, output_type="ndarray"), pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))])
    assert index.query_pairs(1e-9, output_type="ndarray").shape == (0, 2)


def test_empty_index(backend):
    index = SpatialIndex(np.zeros((0, 3)))
    d, i = index.query(np.eye(3))
    np.testing.assert_array_equal(d, np.full(3, np.inf))
    np.testing.assert_array_equal(i, np.zeros(3, dtype=np.int64))
    assert index.query_pairs(1.0) == set()


@needs_rust
def test_exact_ties_take_the_lowest_index():
    rng = np.random.default_rng(3)
    base = _sphere(rng, 400)
    order = rng.permutation(800)
    points = np.empty((800, 3))
    points[order[:400]] = base
    points[order[400:]] = base
    lowest = np.minimum(order[:400], order[400:])
    d, i = SpatialIndex(points).query(base)
    np.testing.assert_array_equal(i, lowest)
    np.testing.assert_array_equal(d, np.zeros(400))
    _, pair = SpatialIndex(points).query(base, k=2)
    np.testing.assert_array_equal(pair, np.sort(np.stack([order[:400], order[400:]], axis=1), axis=1))


@needs_rust
def test_probe_reproduces_scipy_rounding():
    """The probed rounding reproduces scipy exactly on queries where the two roundings differ."""
    rng = np.random.default_rng(11)
    points, queries = _sphere(rng, 2000), _sphere(rng, 20000)
    fused, plain = (spatial_index._KdIndex(points, flag).query(queries, 1, np.inf)[0] for flag in (True, False))
    assert np.count_nonzero(fused != plain) > 0
    expected, _ = cKDTree(points).query(queries)
    np.testing.assert_array_equal(fused if spatial_index.FUSED else plain, expected)


def test_unsupported_queries_raise(backend, cloud):
    index = SpatialIndex(cloud[0])
    with pytest.raises(ValueError):
        index.query(cloud[1], p=1)
    with pytest.raises(ValueError):
        index.query(cloud[1], eps=0.1)
    with pytest.raises(ValueError):
        index.query(cloud[1], k=0)
    with pytest.raises(ValueError):
        SpatialIndex(np.zeros((4, 2)))


def test_data_is_read_only_and_copies_share(backend, cloud):
    index = SpatialIndex(cloud[0])
    assert not index.data.flags.writeable and index.n == len(cloud[0])
    assert copy.copy(index) is index and copy.deepcopy(index) is index
    restored = pickle.loads(pickle.dumps(index))
    np.testing.assert_array_equal(restored.data, index.data)
    for got, want in zip(restored.query(cloud[1]), index.query(cloud[1]), strict=True):
        np.testing.assert_array_equal(got, want)


def test_disable_rust_environment():
    env = {**os.environ, "MANTLE_BLOOM_DISABLE_RUST": "1"}
    subprocess.run(
        [sys.executable, "-c", "from app import spatial_index as s; assert not s.RUST_AVAILABLE and s.SpatialIndex([[1.0, 0, 0]]).backend == 'scipy'"],
        cwd=Path(__file__).resolve().parents[1], env=env, check=True,
    )
