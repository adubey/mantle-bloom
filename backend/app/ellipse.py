"""Minimum-volume enclosing ellipse (MVEE) for a flat 2D point set -- Khachiyan's algorithm.

Deliberately sphere-agnostic: this only knows about plain 2D points, independent of how they
got there (see `plates.plate_bounding_ellipse`, which projects a plate's points into a local
flat-km coordinate system first and calls this). Kept as its own module so the numerics are
testable in isolation from any sphere/plate concerns.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Ellipse2D:
    center: np.ndarray  # (2,)
    semi_major: float
    semi_minor: float
    angle_rad: float  # major axis angle from +x, counter-clockwise


def _two_point_ellipse(points: np.ndarray) -> Ellipse2D:
    a, b = points[0], points[1]
    center = (a + b) / 2.0
    delta = b - a
    dist = float(np.hypot(delta[0], delta[1]))
    angle = float(np.arctan2(delta[1], delta[0])) if dist > 1e-12 else 0.0
    return Ellipse2D(center=center, semi_major=dist / 2.0, semi_minor=0.0, angle_rad=angle)


def min_enclosing_ellipse(points_xy: np.ndarray, tol: float = 1e-6, max_iter: int = 300) -> Ellipse2D:
    """The minimum-area ellipse containing every point in `points_xy` (N, 2) -- Khachiyan's
    algorithm for the d=2 minimum-volume enclosing ellipsoid, with a containment-safety
    post-step (see below) so the result reliably contains every input point rather than only
    doing so up to `tol`."""
    points_xy = np.asarray(points_xy, dtype=float)
    n = len(points_xy)
    if n == 0:
        raise ValueError("min_enclosing_ellipse requires at least one point")
    if n == 1:
        return Ellipse2D(center=points_xy[0].copy(), semi_major=0.0, semi_minor=0.0, angle_rad=0.0)
    if n == 2:
        return _two_point_ellipse(points_xy)
    # All points effectively coincide (relative to their own spread) -- handled explicitly,
    # like the N==1 case, rather than leaning on the covariance regularizer below: the
    # relative scale that regularizer is anchored to (the data's own spread) is itself ~0
    # here, so no fixed epsilon reliably rescues it from a singular/near-singular inverse.
    spread = float(np.ptp(points_xy, axis=0).max())
    scale_ref = float(np.abs(points_xy).max())
    if spread <= 1e-9 * max(scale_ref, 1.0):
        return Ellipse2D(center=points_xy.mean(axis=0), semi_major=0.0, semi_minor=0.0, angle_rad=0.0)

    d = 2
    q = np.vstack([points_xy.T, np.ones(n)])  # (3, N)
    u = np.full(n, 1.0 / n)

    for _ in range(max_iter):
        x = q @ (u[:, None] * q.T)  # (3, 3)
        x += 1e-12 * np.eye(3)  # regularize against collinear/near-degenerate point sets
        x_inv = np.linalg.inv(x)
        m = np.einsum("ij,ji->i", q.T, x_inv @ q)  # diag(Q^T X^-1 Q), length N
        j = int(np.argmax(m))
        max_m = float(m[j])
        if max_m - (d + 1) <= tol * (d + 1):
            break
        step = (max_m - (d + 1)) / ((d + 1) * (max_m - 1))
        u = (1.0 - step) * u
        u[j] += step

    center = points_xy.T @ u
    cov = (points_xy.T * u) @ points_xy - np.outer(center, center)
    # Regularized relative to the data's own scale (a fixed epsilon could be negligible at
    # large coordinate values, or unnecessarily large at tiny ones): perfectly (or nearly)
    # collinear points make `cov` singular -- rank 1 in exact arithmetic -- which would
    # otherwise crash the inverse below. This turns that into a very thin sliver ellipse
    # instead, which is the geometrically sensible degenerate answer.
    scale = max(float(np.trace(cov)), 1e-12)
    shape = np.linalg.inv(cov + 1e-9 * scale * np.eye(2)) / d  # (x-c)^T shape (x-c) <= 1

    # Containment safety: Khachiyan only converges to within `tol`, not exactly -- shrink
    # `shape` (grow the ellipse) just enough to cover whichever input point is worst, so the
    # result is *guaranteed* to contain every point rather than merely "usually does."
    centered = points_xy - center
    mahalanobis = np.einsum("ni,ij,nj->n", centered, shape, centered)
    worst = float(mahalanobis.max())
    if worst > 1.0:
        shape = shape / worst

    eigvals, eigvecs = np.linalg.eigh(shape)
    semi_axes = 1.0 / np.sqrt(eigvals)
    # eigh returns ascending eigenvalues -- the *smaller* eigenvalue gives the *larger* axis.
    minor_idx, major_idx = int(np.argmax(eigvals)), int(np.argmin(eigvals))
    major_vec = eigvecs[:, major_idx]
    return Ellipse2D(
        center=center,
        semi_major=float(semi_axes[major_idx]),
        semi_minor=float(semi_axes[minor_idx]),
        angle_rad=float(np.arctan2(major_vec[1], major_vec[0])),
    )
