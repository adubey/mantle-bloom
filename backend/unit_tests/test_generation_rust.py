"""Parity checks for native world-generation kernels (#1042)."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import coo_matrix

from app import _generation_rust, plate_partition, synthetic_relief


def _unit(rng, n):
    xyz = rng.normal(size=(n, 3))
    return xyz / np.linalg.norm(xyz, axis=1, keepdims=True)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_priority_flood_matches_reference(monkeypatch):
    rng = np.random.default_rng(1042)
    for n in (1, 7, 200):
        src = rng.integers(0, n, n * 4)
        dst = rng.integers(0, n, n * 4)
        edges = coo_matrix((rng.random(n * 4), (dst, src)), shape=(n, n)).tocsr()
        edges.sort_indices()
        elevation = rng.uniform(-100.0, 1000.0, n)
        seeds = np.where(rng.random(n) < 0.3, elevation, np.inf)
        for gradient in (None, edges.data * 0.01):
            actual = synthetic_relief._flood(edges, seeds.copy(), elevation, gradient)
            monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
            expected = synthetic_relief._flood(edges, seeds.copy(), elevation, gradient)
            monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
            np.testing.assert_array_equal(actual, expected)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_nearest_distance_matches_reference(monkeypatch):
    rng = np.random.default_rng(1042)
    sources, points = _unit(rng, 37), _unit(rng, 500)
    actual_dist, actual_idx = synthetic_relief._distance_km(sources, points)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
    expected_dist, expected_idx = synthetic_relief._distance_km(sources, points)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
    np.testing.assert_array_equal(actual_idx, expected_idx)
    # Distances agree to the last ulp: scipy accumulates the chord in a
    # different order, so a few entries differ by ~2e-12 km (measured max).
    # Indices are exact; a comparison can only flip within that distance.
    np.testing.assert_allclose(actual_dist, expected_dist, rtol=0.0, atol=1e-11)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_nearest_seed_labels_matches_kdtree():
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(1042)
    for n_seeds, n_points in ((1, 50), (37, 800), (900, 5000)):
        xyz, seeds = _unit(rng, n_points), _unit(rng, n_seeds)
        owners = rng.integers(0, 9, n_seeds)
        _, nearest = cKDTree(seeds).query(xyz, workers=-1)
        np.testing.assert_array_equal(_generation_rust.nearest_seed_labels(xyz, seeds, owners), owners[nearest])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_nearest_seed_labels_ties_take_lowest_seed():
    # Duplicated seeds are exact ties: the native tree resolves them to the
    # lowest seed index (`cKDTree` promises no tie order).
    seeds = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    owners = np.array([5, 7, 9])
    labels = _generation_rust.nearest_seed_labels(np.eye(3), seeds, owners)
    np.testing.assert_array_equal(labels, [5, 9, 5])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_seed_labels_tie_contract_matches_rust():
    # Exact ties (duplicated seeds, symmetric points) must resolve to the
    # lowest seed index in the Python fallback just as in the native tree, so
    # a seeded partition does not depend on the extension being installed.
    seeds = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    owners = np.array([5, 7, 9, 11])
    points = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    expected = _generation_rust.nearest_seed_labels(points, seeds, owners)
    np.testing.assert_array_equal(plate_partition._seed_labels_python(points, seeds, owners), expected)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_polygon_seed_partition_matches_reference(monkeypatch):
    from app.elevation_lines import line_spacing_rad
    from app.layout_lattice import lattice

    lat = lattice(line_spacing_rad(0.25))
    core = np.full(lat.size, -1, dtype=np.int64)
    core[::37] = 0  # a thin spread of craton cells, so polygon seeds join in

    def run():
        return plate_partition._polygon_seed_partition(np.random.default_rng(1042), lat, 9, core)

    actual_labels, actual_sites = run()
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
    expected_labels, expected_sites = run()
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
    np.testing.assert_array_equal(actual_labels, expected_labels)
    np.testing.assert_array_equal(actual_sites, expected_sites)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_label_components_matches_reference(monkeypatch):
    from app.elevation_lines import line_spacing_rad
    from app.layout_lattice import lattice

    lat = lattice(line_spacing_rad(0.25))
    rng = np.random.default_rng(1042)
    labels = rng.integers(0, 5, lat.size)
    actual = plate_partition.label_components(lat, labels)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
    expected = plate_partition.label_components(lat, labels)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
    assert actual == expected


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_drain_land_matches_reference(monkeypatch):
    rng = np.random.default_rng(1042)
    points = _unit(rng, 300)
    elevation = rng.uniform(-500.0, 1500.0, len(points))
    land = elevation > 0.0
    area = np.full(len(points), 100.0)
    actual = synthetic_relief.drain_land(points, elevation.copy(), land, area, seed=7)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
    expected = synthetic_relief.drain_land(points, elevation.copy(), land, area, seed=7)
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_nearest_distance_empty_single_and_ties():
    points = np.eye(3)
    dist, index = _generation_rust.nearest_distance(np.zeros((0, 3)), points)
    assert np.all(np.isinf(dist)) and not index.any()
    dist, index = _generation_rust.nearest_distance(np.array([[0.0, 1.0, 0.0]]), points)
    np.testing.assert_array_equal(index, [0, 0, 0])
    assert dist[1] == 0.0
    # Exact ties (duplicated sources) resolve to the lowest source index;
    # `cKDTree` promises no tie order, so this pins the native contract.
    sources = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    _, index = _generation_rust.nearest_distance(sources, points)
    np.testing.assert_array_equal(index, [0, 2, 0])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_nearest_gate_handover(monkeypatch):
    sources = _unit(np.random.default_rng(1042), 40)
    assert _generation_rust.use_nearest(sources)
    monkeypatch.setattr(_generation_rust, "NEAREST_RUST_MAX_SOURCES", 8)
    assert not _generation_rust.use_nearest(sources)
    points = _unit(np.random.default_rng(7), 300)
    gated = synthetic_relief._distance_km(sources, points)  # gate open: native path
    monkeypatch.setattr(_generation_rust, "NEAREST_RUST_MAX_SOURCES", 8)
    tree = synthetic_relief._distance_km(sources, points)  # gate closed: cKDTree
    monkeypatch.setattr(_generation_rust, "NEAREST_RUST_MAX_SOURCES", 512)
    np.testing.assert_allclose(gated[0], tree[0], rtol=0.0, atol=1e-11)
    np.testing.assert_array_equal(gated[1], tree[1])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_priority_flood_rejects_nan():
    with pytest.raises(ValueError, match="NaN"):
        _generation_rust.priority_flood([0, 0], [], None, [np.nan], [0.0])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_kernels_reject_malformed_csr_and_nan():
    with pytest.raises(ValueError, match="CSR"):
        _generation_rust.priority_flood([0, 2, 1], [0], None, [0.0, 0.0], [0.0, 0.0])
    with pytest.raises(ValueError, match="CSR"):
        _generation_rust.label_component_count([1, 1], [0], [0])
    with pytest.raises(ValueError, match="NaN"):
        _generation_rust.nearest_distance([[np.nan, 0.0, 0.0]], [[1.0, 0.0, 0.0]])
    with pytest.raises(ValueError, match="NaN"):
        _generation_rust.nearest_seed_labels([[1.0, 0.0, 0.0]], [[np.nan, 0.0, 0.0]], [0])


@pytest.mark.skipif(not _generation_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_build_geography_matches_reference(monkeypatch):
    from app import continent_layout
    from app.elevation_lines import line_spacing_rad

    def run():
        geo = continent_layout.build_geography(np.random.default_rng(1042), line_spacing_rad(0.25), 0.3)
        return geo.continent, geo.core, geo.belt, geo.margin, geo.nuclei, geo.nucleus_craton

    actual = run()
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", False)
    expected = run()
    monkeypatch.setattr(_generation_rust, "RUST_AVAILABLE", True)
    for left, right in zip(actual, expected, strict=True):
        np.testing.assert_array_equal(left, right)


def test_disable_rust_environment():
    env = {**os.environ, "MANTLE_BLOOM_DISABLE_RUST": "1"}
    subprocess.run(
        [sys.executable, "-c", "from app import _generation_rust; assert not _generation_rust.RUST_AVAILABLE"],
        cwd=Path(__file__).resolve().parents[1], env=env, check=True,
    )


def test_stale_generation_extension_warns(caplog):
    _generation_rust._warn_missing_exports(object())
    assert "stale build" in caplog.text
    assert "priority_flood_rs" in caplog.text
