"""Realistic-grid parity for the native diagnostic climate kernels (#1046)."""

import numpy as np
import pytest

from app.hydroclimate import _climate_rust, climate


def test_diagnostic_climate_kernels_match_numpy_at_density_four(monkeypatch):
    exports = ("humidity_zonal_sweep_rs", "smooth_along_coast_rs", "offset_indices_rs", "moisture_divergence_rs", "smooth_field_rs")
    if not all(_climate_rust.has_kernel(name) for name in exports):
        pytest.skip("Rust climate kernels unavailable")
    rng = np.random.default_rng(1046)
    shape = (360, 720)
    lat = np.linspace(89.75, -89.75, shape[0])
    ocean = rng.random(shape) < 0.6
    elevation = rng.uniform(-4000, 6000, shape)
    humidity = rng.uniform(0, 18, shape)
    u = rng.uniform(-1, 1, shape)
    v = rng.uniform(-1, 1, shape)
    factor = rng.uniform(0, 1, shape)
    source = rng.uniform(0, 1, shape)

    def calculate():
        return (
            *climate._humidity_zonal_sweep(ocean, elevation, humidity, factor, lat, source),
            *climate._smooth_along_coast(u, v, ocean, lat, u, v),
            climate._sample_at_offset(elevation, u, v, 1.0, lat),
            climate.compute_moisture_flux_convergence(humidity, u, v, lat),
        )

    actual = calculate()
    monkeypatch.setattr(_climate_rust, "has_kernel", lambda _name: False)
    expected = calculate()
    for got, want in zip(actual, expected, strict=True):
        np.testing.assert_allclose(got, want, rtol=2e-14, atol=2e-14)
    np.testing.assert_array_equal(actual[4], expected[4])
