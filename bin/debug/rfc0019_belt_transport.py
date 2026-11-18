"""Broaden zonal convergence within fixed latitude windows for the G3 screen."""

from __future__ import annotations

import numpy as np

from rfc0019_transport_prototypes import EARTH_RADIUS_M

BELT_CENTERS_DEG = (-60.0, 0.0, 60.0)
HALF_WIDTH_DEG = 8.0
SMOOTHING_SIGMA_DEG = 4.0


def broaden_convergence(v: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    """Return a zonally uniform velocity correction on south-to-north latitude faces."""
    edges = np.linspace(-90.0, 90.0, v.shape[0] + 1)
    lengths = 2 * np.pi * EARTH_RADIUS_M * np.cos(np.deg2rad(edges[1:-1]))
    mean = v.mean(axis=1)
    original = np.r_[0.0, lengths * (mean[:-1] + mean[1:]) / 2, 0.0]
    broadened = original.copy()
    diagnostics = []
    for center in BELT_CENTERS_DEG:
        limits = np.array([center - HALF_WIDTH_DEG, center + HALF_WIDTH_DEG])
        first, last = np.searchsorted(edges, limits)
        if not np.allclose(edges[[first, last]], limits):
            raise ValueError("convergence windows must start and end on grid faces")
        window = original[first:last + 1]
        convergence = -np.diff(window)
        positive = np.maximum(convergence, 0)
        centers = (edges[first:last] + edges[first + 1:last + 1]) / 2
        separation = centers[:, None] - centers[None, :]
        kernel = np.exp(-0.5 * (separation / SMOOTHING_SIGMA_DEG)**2)
        # Each source row redistributes exactly its own positive convergence.
        adjusted = np.minimum(convergence, 0) + (kernel / kernel.sum(axis=0)) @ positive
        rebuilt = window[0] - np.cumsum(adjusted)
        roundoff = 64 * np.finfo(float).eps * np.max(np.abs(window))
        np.testing.assert_allclose(rebuilt[-1], window[-1], rtol=0, atol=roundoff)
        broadened[first + 1:last] = rebuilt[:-1]
        diagnostics.append({
            "latitude_bounds_deg": limits.tolist(),
            "peak_convergence_before_m2_s": float(convergence.max()),
            "peak_convergence_after_m2_s": float(adjusted.max()),
            "edge_transport_before_and_after_m2_s": [float(window[0]), float(window[-1])],
        })
    return (broadened[1:-1] - original[1:-1]) / lengths, diagnostics
