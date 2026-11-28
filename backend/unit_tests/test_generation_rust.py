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
    np.testing.assert_allclose(actual_dist, expected_dist, rtol=0.0, atol=1e-9)


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
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1e-9)


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
