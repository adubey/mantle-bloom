"""Parity checks for native hydrology routing (#1041)."""

import numpy as np
import pytest
import os
import subprocess
import sys
from pathlib import Path

from app.hydroclimate import _routing_rust, hydrology


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_flow_direction_matches_reference(monkeypatch):
    rng = np.random.default_rng(1041)
    for n in (1, 20, 300):
        elevation = rng.integers(-5, 6, n).astype(float)
        ocean = rng.random(n) < 0.2
        neighbors = rng.integers(0, n, (n, 8))
        offset = rng.uniform(-2, 2, n)
        filled = elevation + rng.integers(-1, 2, n)
        spill = rng.integers(-1, n, n)
        channel = rng.choice([0.0, 5.0, 6.0, 10.0], n)
        args = elevation, ocean, neighbors, offset, filled, spill, channel
        actual = hydrology._compute_flow_direction(*args)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
        expected = hydrology._compute_flow_direction(*args)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", True)
        for left, right in zip(actual, expected, strict=True):
            np.testing.assert_array_equal(left, right)


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_flow_direction_nonfinite_and_strided_inputs(monkeypatch):
    elevation = np.array([10.0, 0, -np.inf, 0, 3.0, 0, np.nan, 0])[::2]
    ocean = np.array([False, True, False, False])
    neighbors = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])[:, ::-1]
    args = elevation, ocean, neighbors, np.zeros(4), elevation.copy(), np.arange(4), np.array([0., 7., np.inf, 8.])
    actual = hydrology._compute_flow_direction(*args)
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    expected = hydrology._compute_flow_direction(*args)
    for left, right in zip(actual, expected, strict=True):
        np.testing.assert_array_equal(left, right)


def test_disable_rust_environment():
    env = {**os.environ, "MANTLE_BLOOM_DISABLE_RUST": "1"}
    subprocess.run(
        [sys.executable, "-c", "from app.hydroclimate import _routing_rust; assert not _routing_rust.RUST_AVAILABLE"],
        cwd=Path(__file__).resolve().parents[1], env=env, check=True,
    )
