"""Equal-area map projections for rendering the sphere to a 2D canvas.

Both take/return radians and are vectorized over numpy arrays of any shape. Output is in
planar units for a unit-radius sphere (multiply by the desired pixels-per-radian scale when
rendering).
"""

from __future__ import annotations

import numpy as np

BEHRMANN_STANDARD_PARALLEL = np.radians(30.0)


def behrmann(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cylindrical equal-area projection, standard parallel 30 degrees."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    cos_phi0 = np.cos(BEHRMANN_STANDARD_PARALLEL)
    x = lon * cos_phi0
    y = np.sin(lat) / cos_phi0
    return x, y


_ECKERT4_C = 2.0 / np.sqrt(np.pi * (4.0 + np.pi))
_ECKERT4_RHS_SCALE = 2.0 + np.pi / 2.0
_NEWTON_ITERS = 30


def _eckert4_theta(lat: np.ndarray) -> np.ndarray:
    """Solve theta + sin(theta)cos(theta) + 2 sin(theta) = (2 + pi/2) sin(lat) for theta,
    via Newton's method, vectorized."""
    lat = np.asarray(lat, dtype=float)
    target = _ECKERT4_RHS_SCALE * np.sin(lat)
    theta = lat / 2.0  # standard starting guess
    for _ in range(_NEWTON_ITERS):
        sin_t, cos_t = np.sin(theta), np.cos(theta)
        f = theta + sin_t * cos_t + 2.0 * sin_t - target
        fprime = 2.0 * cos_t * (cos_t + 1.0)
        safe_fprime = np.where(np.abs(fprime) < 1e-9, 1e-9, fprime)
        theta = theta - f / safe_fprime
    # Exact at the poles, where the Newton denominator vanishes.
    theta = np.where(np.isclose(np.abs(lat), np.pi / 2), np.sign(lat) * (np.pi / 2), theta)
    return theta


def eckert4(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pseudocylindrical equal-area projection (Snyder 1987)."""
    lon = np.asarray(lon, dtype=float)
    theta = _eckert4_theta(lat)
    x = _ECKERT4_C * lon * (1.0 + np.cos(theta))
    y = _ECKERT4_C * np.pi * np.sin(theta)
    return x, y


PROJECTIONS = {
    "behrmann": behrmann,
    "eckert4": eckert4,
}


def project(name: str, lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        fn = PROJECTIONS[name]
    except KeyError as exc:
        raise ValueError(f"Unknown projection {name!r}; choices are {sorted(PROJECTIONS)}") from exc
    return fn(lat, lon)
