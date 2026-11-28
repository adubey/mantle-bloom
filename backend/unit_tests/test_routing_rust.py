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


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
@pytest.mark.parametrize("with_erosion", [False, True])
def test_downstream_matches_reference_with_cycles_and_losses(monkeypatch, with_erosion):
    rng = np.random.default_rng(1046)
    for _ in range(20):
        n = 80
        elevation = rng.integers(0, 10, n).astype(float)
        ocean = rng.random(n) < 0.15
        target = rng.integers(-1, n, n)
        source, retain, loss, erosion = (rng.random(n) for _ in range(4))
        args = elevation, ocean, target, source, retain, loss, erosion if with_erosion else None
        actual = hydrology.route_downstream(*args)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
        expected = hydrology.route_downstream(*args)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", True)
        for name in ("through_flux", "deposited", "retained", "sink_settlement",
                     "ocean_export", "cycle_settlement", "evaporated"):
            np.testing.assert_array_equal(getattr(actual, name), getattr(expected, name))
        if with_erosion:
            np.testing.assert_array_equal(actual.erosion_only_flux, expected.erosion_only_flux)
        assert actual.residual == expected.residual


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_routing_order_matches_reference_for_uphill_edges_and_ties(monkeypatch):
    rng = np.random.default_rng(1046)
    for _ in range(50):
        n = 100
        elevation = rng.integers(-5, 5, n).astype(float)
        ocean = rng.random(n) < 0.1
        target = rng.integers(-1, n, n)
        actual = hydrology._routing_order(elevation, ocean, target)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
        expected = hydrology._routing_order(elevation, ocean, target)
        monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", True)
        assert actual == expected


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
@pytest.mark.parametrize("dt", [0.0, 100.0])
def test_river_reach_sweep_matches_reference(monkeypatch, dt):
    rng = np.random.default_rng(1046)
    n = 100
    order = np.arange(n, dtype=np.int64)
    target = np.minimum(order + 1, n - 1)
    target[-1] = -1
    ocean = np.zeros(n, dtype=bool)
    ocean[-2] = True
    stores = rng.random(n) < 0.6
    tau = rng.uniform(10.0, 1000.0, n)
    storage, runoff, melt, prior, runoff_loss, melt_loss = (rng.random(n) for _ in range(6))

    def run():
        inputs = [runoff.copy(), melt.copy(), prior.copy()]
        result = hydrology._sweep_river_reaches(
            order, target, ocean, stores, tau, dt, storage, *inputs, runoff_loss, melt_loss,
        )
        return result, inputs

    actual, actual_inputs = run()
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    expected, expected_inputs = run()
    for a, b in zip(list(actual) + actual_inputs, list(expected) + expected_inputs, strict=True):
        np.testing.assert_array_equal(a, b)


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_river_reach_sweep_falls_back_for_strided_mutable_inputs(monkeypatch):
    n = 8
    order = np.arange(n, dtype=np.int64)
    target = np.full(n, -1, dtype=np.int64)
    ocean = np.zeros(n, dtype=bool)
    stores = np.zeros(n, dtype=bool)
    zeros = np.zeros(n)
    runoff = np.ones(n * 2)[::2]
    args = order, target, ocean, stores, np.ones(n), 1.0, zeros, runoff, zeros.copy(), zeros.copy(), zeros, zeros
    monkeypatch.setattr(_routing_rust._rust, "sweep_river_reaches_rs", lambda *args: pytest.fail("native path used"))
    result = hydrology._sweep_river_reaches(*args)
    np.testing.assert_array_equal(result[3], runoff)


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
@pytest.mark.parametrize("elevation", [
    [1.0, 2.0, np.inf, -np.inf, np.nan],
    [-np.inf, np.inf, 0.0, 1.0, np.inf],
])
def test_order_nonfinite_ocean_and_land_infinities_use_native(monkeypatch, elevation):
    elevation = np.array(elevation)
    ocean = np.array([False, False, False, False, True])
    target = np.array([1, -1, 3, -1, -1])
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    expected = hydrology._routing_order(elevation, ocean, target)
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", True)
    calls = []
    original = _routing_rust.routing_order
    def native(*args):
        calls.append(True)
        return original(*args)
    monkeypatch.setattr(_routing_rust, "routing_order", native)
    assert hydrology._routing_order(elevation, ocean, target) == expected
    assert calls


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_order_nan_land_uses_python_and_native_rejects_nan(monkeypatch):
    elevation = np.array([1.0, 2.0, np.nan])
    ocean = np.zeros(3, dtype=bool)
    target = np.array([1, -1, -1])
    with pytest.raises(ValueError, match="NaN"):
        _routing_rust.routing_order(elevation, ocean, target)
    monkeypatch.setattr(_routing_rust, "routing_order", lambda *args: pytest.fail("native path used"))
    result = hydrology._routing_order(elevation, ocean, target)
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    assert result == hydrology._routing_order(elevation, ocean, target)


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_river_reach_aliased_mutable_input_uses_fallback(monkeypatch):
    n = 4
    order = np.arange(n, dtype=np.int64)
    target = np.full(n, -1, dtype=np.int64)
    mask = np.zeros(n, dtype=bool)
    shared = np.ones(n)
    zeros = np.zeros(n)
    monkeypatch.setattr(_routing_rust, "sweep_river_reaches", lambda *args: pytest.fail("native path used"))
    actual = hydrology._sweep_river_reaches(order, target, mask, mask.copy(), np.ones(n), 1.0,
                                           zeros, shared, shared, zeros.copy(), zeros, zeros)
    np.testing.assert_array_equal(actual[3], np.ones(n))
    np.testing.assert_array_equal(actual[4], np.ones(n))


def test_stale_routing_extension_warns(caplog):
    _routing_rust._warn_missing_exports(object())
    assert "stale build" in caplog.text
    assert "routing_order_rs" in caplog.text


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_native_avoids_erosion_array_when_not_requested():
    zeros = np.zeros(2)
    arrays = _routing_rust.route_downstream([0, 1], [], [-1, -1], [False, False],
                                           [1.0, 2.0], zeros, zeros, None)
    assert arrays[-1].size == 0
