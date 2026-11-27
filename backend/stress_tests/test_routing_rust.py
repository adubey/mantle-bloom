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
