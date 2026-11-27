"""Large-graph parity for native hydrology routing (#1046)."""

import numpy as np
import pytest

from app.hydroclimate import _routing_rust, hydrology


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_large_routing_graph_matches_python(monkeypatch):
    rng = np.random.default_rng(1046)
    n = 300_000
    elevation = rng.integers(0, 100, n).astype(float)
    ocean = rng.random(n) < 0.1
    target = rng.integers(-1, n, n)
    source, retain, loss, erosion = (rng.random(n) for _ in range(4))
    args = elevation, ocean, target, source, retain, loss, erosion
    actual_order = hydrology._routing_order(elevation, ocean, target)
    actual = hydrology.route_downstream(*args)
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    expected_order = hydrology._routing_order(elevation, ocean, target)
    expected = hydrology.route_downstream(*args)
    assert actual_order == expected_order
    for name in ("through_flux", "deposited", "retained", "sink_settlement", "ocean_export",
                 "cycle_settlement", "evaporated", "erosion_only_flux"):
        np.testing.assert_array_equal(getattr(actual, name), getattr(expected, name))
    assert actual.residual == expected.residual


@pytest.mark.skipif(not _routing_rust.RUST_AVAILABLE, reason="Rust extension unavailable")
def test_large_river_reach_sweep_matches_python(monkeypatch):
    rng = np.random.default_rng(1047)
    n = 300_000
    order = np.arange(n, dtype=np.int64)
    target = np.minimum(order + 1, n - 1)
    target[-1] = -1
    ocean = np.zeros(n, dtype=bool)
    stores = rng.random(n) < 0.5
    tau = rng.uniform(100, 1000, n)
    storage, runoff, melt, prior, runoff_loss, melt_loss = (rng.random(n) for _ in range(6))

    def run():
        inputs = [runoff.copy(), melt.copy(), prior.copy()]
        result = hydrology._sweep_river_reaches(
            order, target, ocean, stores, tau, 100.0, storage, *inputs, runoff_loss, melt_loss,
        )
        return list(result) + inputs

    actual = run()
    monkeypatch.setattr(_routing_rust, "RUST_AVAILABLE", False)
    expected = run()
    for a, b in zip(actual, expected, strict=True):
        np.testing.assert_array_equal(a, b)
