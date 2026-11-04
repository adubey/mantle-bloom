"""Super arrows for the wind and ocean-current views (#805).

The views used to draw one short arrow per sampled cell, which shows the local field but not
its large-scale pattern: a trade-wind belt or a gyre reads as a field of ticks. Atlas maps
instead draw a few long, curved arrows that follow the prevailing flow. This module finds
those arrows in a vector field on the climate module's equirectangular grid. It only
summarises a field the climate model has already solved, so it changes no model state.

The method follows the issue's sketch, made concrete as evenly-spaced streamline placement
(Jobard & Lefer 1997):

1. **Cluster.** Smooth the field over a few degrees and score each place by speed times the
   *coherence* of its neighbourhood -- the length of the mean vector over the mean speed,
   1 where every vector points the same way, 0 where they cancel. The highest-scoring places
   not yet claimed by an arrow are the centres of coherent clusters, and seed the arrows.
2. **Grow.** Each seed grows into a streamline, forward and then backward, one small step at
   a time, until it comes within `stop_deg` of an existing arrow, runs off the valid region
   (land, for currents), slows below the speed floor, has turned through `max_turn_deg`, or
   reaches `max_length_deg`. The turn limit splits a gyre into a few arcs, as atlas maps do.
3. **Link.** Points every `control_spacing_deg` along the streamline become the control
   points of a Catmull-Rom spline, and the spline is the final arrow's path. Drawing (the
   tapered body and the head) is `render_image._draw_flow_arrows`.

Accepted arrows stamp a distance field on a one-degree grid, so later seeds and traces see
them in constant time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy import ndimage


@dataclass(frozen=True)
class FlowArrowParams:
    """Spacing and stopping rules for `trace_flow_arrows`. Angles are degrees of arc."""

    separation_deg: float = 9.0  # no new arrow starts closer than this to an existing one
    stop_deg: float = 4.0  # a growing arrow stops this close to an existing one
    step_deg: float = 0.75  # streamline integration step
    min_length_deg: float = 10.0  # shorter traces are dropped, not drawn
    max_length_deg: float = 42.0
    max_turn_deg: float = 130.0  # total turning along one arrow
    smoothing_deg: float = 3.0  # Gaussian sigma for the prevailing field
    min_speed_fraction: float = 0.2  # speed floor, as a fraction of the 90th-percentile speed
    min_coherence: float = 0.6  # seeds need at least this neighbourhood coherence
    seed_spacing_deg: float = 3.0  # candidate seed lattice
    control_spacing_deg: float = 4.0  # spline control point spacing
    max_abs_lat_deg: float = 82.0  # east/north degenerate at the poles; stay clear of them


# Wind belts are broad and nearly zonal, so the same spacing as currents draws a dense comb of
# parallel arrows; spacing them wider keeps a few arrows per belt, like the reference figure.
WIND_PARAMS = FlowArrowParams(separation_deg=18.0, stop_deg=6.0)
CURRENT_PARAMS = FlowArrowParams(separation_deg=10.0)


@dataclass(frozen=True)
class FlowArrow:
    """One super arrow, tail first."""

    points_xyz: np.ndarray  # (N, 3) unit vectors along the spline, tail to tip
    mean_speed: float  # smoothed field speed averaged along the arrow, in the field's units


_OCCUPANCY_RES_DEG = 1.0


def trace_flow_arrows(u: np.ndarray, v: np.ndarray, valid: np.ndarray, params: FlowArrowParams = FlowArrowParams()) -> list[FlowArrow]:
    """Super arrows for an eastward `u` / northward `v` field on an `(H, W)` equirectangular
    grid (row 0 at the north pole, column 0 at longitude -180, the climate module's layout).
    `valid` is where the field is defined (all True for wind, ocean for currents); arrows
    never cross an invalid cell. Deterministic: the same field always gives the same arrows."""
    u_s, v_s, coherence = _prevailing_field(u, v, valid, params.smoothing_deg)
    speed = np.hypot(u_s, v_s)
    if not np.any(valid):
        return []
    reference_speed = float(np.percentile(speed[valid], 90))
    if reference_speed <= 0.0:
        return []
    min_speed = params.min_speed_fraction * reference_speed

    mask = valid.astype(np.float64)
    occ_h, occ_w = int(round(180 / _OCCUPANCY_RES_DEG)), int(round(360 / _OCCUPANCY_RES_DEG))
    distance = np.full((occ_h, occ_w), np.inf)
    step = np.radians(params.step_deg)
    max_steps = int(np.ceil(params.max_length_deg / params.step_deg)) + 2
    stop_rad = np.radians(params.stop_deg)
    max_abs_lat = np.radians(params.max_abs_lat_deg)

    arrows: list[FlowArrow] = []
    for seed in _seed_candidates(speed, coherence, valid, min_speed, params):
        if _lookup_distance(distance, seed) < np.radians(params.separation_deg):
            continue
        pts_fwd = np.empty((max_steps + 1, 3))
        pts_bwd = np.empty((max_steps + 1, 3))
        n_fwd, length, turn = _trace(
            seed, 1.0, u_s, v_s, mask, distance, step, np.radians(params.max_length_deg), np.radians(params.max_turn_deg),
            stop_rad, min_speed, max_abs_lat, pts_fwd,
        )
        n_bwd, length_b, _ = _trace(
            seed, -1.0, u_s, v_s, mask, distance, step, np.radians(params.max_length_deg) - length,
            np.radians(params.max_turn_deg) - turn, stop_rad, min_speed, max_abs_lat, pts_bwd,
        )
        if np.degrees(length + length_b) < params.min_length_deg:
            continue
        path = np.concatenate([pts_bwd[:n_bwd][::-1], pts_fwd[1:n_fwd]])
        _stamp_distance(distance, path, np.radians(params.separation_deg))
        spline = _catmull_rom(path[_control_indices(len(path), params.control_spacing_deg / params.step_deg)])
        arrows.append(FlowArrow(points_xyz=spline, mean_speed=float(np.mean(_sample_speed(speed, path)))))
    return arrows


def _prevailing_field(u: np.ndarray, v: np.ndarray, valid: np.ndarray, sigma_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Smoothed (u, v) and their coherence, by normalised convolution over `valid` only, so a
    coastal current isn't dragged toward zero by the land beside it. Longitude wraps."""
    height, width = u.shape
    sigma = (sigma_deg / (180.0 / height), sigma_deg / (360.0 / width))
    mask = valid.astype(np.float64)

    def smooth(field: np.ndarray) -> np.ndarray:
        return ndimage.gaussian_filter(field * mask, sigma, mode=("nearest", "wrap"))

    weight = np.maximum(smooth(np.ones_like(mask)), 1e-12)
    u_s, v_s = smooth(u) / weight, smooth(v) / weight
    mean_speed = smooth(np.hypot(u, v)) / weight
    coherence = np.where(mean_speed > 0, np.hypot(u_s, v_s) / np.maximum(mean_speed, 1e-12), 0.0)
    return np.where(valid, u_s, 0.0), np.where(valid, v_s, 0.0), np.where(valid, coherence, 0.0)


def _seed_candidates(speed: np.ndarray, coherence: np.ndarray, valid: np.ndarray, min_speed: float, params: FlowArrowParams) -> np.ndarray:
    """(K, 3) seed points on a roughly equal-area lattice, best cluster centre first."""
    height, width = speed.shape
    lats = np.arange(-params.max_abs_lat_deg + params.seed_spacing_deg / 2, params.max_abs_lat_deg, params.seed_spacing_deg)
    lat_list, lon_list = [], []
    for lat in lats:
        n_lon = max(int(round(360.0 * np.cos(np.radians(lat)) / params.seed_spacing_deg)), 1)
        lat_list.append(np.full(n_lon, lat))
        lon_list.append(-180.0 + (np.arange(n_lon) + 0.5) * 360.0 / n_lon)
    lat_deg, lon_deg = np.concatenate(lat_list), np.concatenate(lon_list)
    rows = np.clip(((90.0 - lat_deg) / 180.0 * height).astype(int), 0, height - 1)
    cols = np.clip(((lon_deg + 180.0) / 360.0 * width).astype(int), 0, width - 1)
    keep = valid[rows, cols] & (speed[rows, cols] >= min_speed) & (coherence[rows, cols] >= params.min_coherence)
    score = (speed[rows, cols] * coherence[rows, cols] ** 2)[keep]
    # Stable sort, so ties resolve by lattice order and the result never depends on numpy's
    # choice of algorithm.
    order = np.argsort(-score, kind="stable")
    lat_r, lon_r = np.radians(lat_deg[keep][order]), np.radians(lon_deg[keep][order])
    return np.stack([np.cos(lat_r) * np.cos(lon_r), np.cos(lat_r) * np.sin(lon_r), np.sin(lat_r)], axis=-1)


def _control_indices(n: int, spacing_steps: float) -> np.ndarray:
    """Every `spacing_steps`-th index of an `n`-point path, spread evenly, both ends kept."""
    count = max(int(round((n - 1) / max(spacing_steps, 1.0))), 1) + 1
    return np.unique(np.round(np.linspace(0, n - 1, count)).astype(int))


def _catmull_rom(control: np.ndarray, samples_per_segment: int = 6) -> np.ndarray:
    """Uniform Catmull-Rom spline through `control` (K, 3), re-projected onto the sphere.
    The curve passes through every control point; the ends repeat their neighbours."""
    if len(control) < 3:
        return control.copy()
    padded = np.concatenate([control[:1] * 2 - control[1:2], control, control[-1:] * 2 - control[-2:-1]])
    t = np.linspace(0.0, 1.0, samples_per_segment, endpoint=False)[:, None]
    t2, t3 = t * t, t * t * t
    pieces = []
    for i in range(1, len(padded) - 2):
        p0, p1, p2, p3 = padded[i - 1], padded[i], padded[i + 1], padded[i + 2]
        pieces.append(0.5 * (2 * p1 + (p2 - p0) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 + (3 * p1 - p0 - 3 * p2 + p3) * t3))
    pieces.append(control[-1:])
    curve = np.concatenate(pieces)
    return curve / np.linalg.norm(curve, axis=1, keepdims=True)


def _sample_speed(speed: np.ndarray, points: np.ndarray) -> np.ndarray:
    height, width = speed.shape
    lat = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    rows = np.clip(((90.0 - lat) / 180.0 * height).astype(int), 0, height - 1)
    cols = np.clip(((lon + 180.0) / 360.0 * width).astype(int), 0, width - 1)
    return speed[rows, cols]


@njit(cache=True)
def _lookup_distance(distance: np.ndarray, p: np.ndarray) -> float:
    occ_h, occ_w = distance.shape
    lat = np.arcsin(min(max(p[2], -1.0), 1.0))
    lon = np.arctan2(p[1], p[0])
    row = min(max(int((np.pi / 2 - lat) / np.pi * occ_h), 0), occ_h - 1)
    col = min(max(int((lon + np.pi) / (2 * np.pi) * occ_w), 0), occ_w - 1)
    return distance[row, col]


@njit(cache=True)
def _sample_bilinear(u: np.ndarray, v: np.ndarray, mask: np.ndarray, lat: float, lon: float) -> tuple[float, float, float]:
    """(u, v, mask) at (lat, lon) in radians, bilinear, longitude wrapping."""
    height, width = u.shape
    fr = (np.pi / 2 - lat) / np.pi * height - 0.5
    fc = (lon + np.pi) / (2 * np.pi) * width - 0.5
    r0 = int(np.floor(fr))
    tr = fr - r0
    if r0 < 0:
        r0, tr = 0, 0.0
    if r0 >= height - 1:
        r0, tr = height - 1, 0.0
    r1 = min(r0 + 1, height - 1)
    c0f = np.floor(fc)
    tc = fc - c0f
    c0 = int(c0f) % width
    c1 = (c0 + 1) % width
    w00, w01, w10, w11 = (1 - tr) * (1 - tc), (1 - tr) * tc, tr * (1 - tc), tr * tc
    su = w00 * u[r0, c0] + w01 * u[r0, c1] + w10 * u[r1, c0] + w11 * u[r1, c1]
    sv = w00 * v[r0, c0] + w01 * v[r0, c1] + w10 * v[r1, c0] + w11 * v[r1, c1]
    sm = w00 * mask[r0, c0] + w01 * mask[r0, c1] + w10 * mask[r1, c0] + w11 * mask[r1, c1]
    return su, sv, sm


@njit(cache=True)
def _direction(u: np.ndarray, v: np.ndarray, mask: np.ndarray, p: np.ndarray, sign: float, min_speed: float):
    """Unit flow direction at `p` as a 3-vector (times `sign`), or zeros when the flow is
    invalid or slower than `min_speed` there."""
    lat = np.arcsin(min(max(p[2], -1.0), 1.0))
    lon = np.arctan2(p[1], p[0])
    su, sv, sm = _sample_bilinear(u, v, mask, lat, lon)
    out = np.zeros(3)
    speed = np.hypot(su, sv)
    if sm < 0.5 or speed < min_speed:
        return out
    sin_lat, cos_lat, sin_lon, cos_lon = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    # east = (-sin lon, cos lon, 0); north = (-sin lat cos lon, -sin lat sin lon, cos lat)
    out[0] = sign * (su * -sin_lon + sv * -sin_lat * cos_lon) / speed
    out[1] = sign * (su * cos_lon + sv * -sin_lat * sin_lon) / speed
    out[2] = sign * (sv * cos_lat) / speed
    return out


@njit(cache=True)
def _advance(p: np.ndarray, d: np.ndarray, h: float) -> np.ndarray:
    q = p + h * d
    return q / np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])


@njit(cache=True)
def _trace(
    seed, sign, u, v, mask, distance, step, max_length, max_turn, stop_rad, min_speed, max_abs_lat, out,
):
    """Grow a streamline from `seed` along `sign` times the flow (midpoint rule on the
    sphere), writing points into `out`. Returns (point count, arc length, turning used)."""
    out[0] = seed
    n = 1
    length = 0.0
    turn = 0.0
    p = seed.copy()
    prev_d = _direction(u, v, mask, p, sign, min_speed)
    if prev_d[0] == 0.0 and prev_d[1] == 0.0 and prev_d[2] == 0.0:
        return n, length, turn
    while n < out.shape[0] and length + step <= max_length + 1e-12:
        d = _direction(u, v, mask, p, sign, min_speed)
        if d[0] == 0.0 and d[1] == 0.0 and d[2] == 0.0:
            break
        mid_d = _direction(u, v, mask, _advance(p, d, 0.5 * step), sign, min_speed)
        if mid_d[0] == 0.0 and mid_d[1] == 0.0 and mid_d[2] == 0.0:
            break
        q = _advance(p, mid_d, step)
        cos_turn = min(max(prev_d[0] * mid_d[0] + prev_d[1] * mid_d[1] + prev_d[2] * mid_d[2], -1.0), 1.0)
        turn_next = turn + np.arccos(cos_turn)
        if turn_next > max_turn or abs(np.arcsin(min(max(q[2], -1.0), 1.0))) > max_abs_lat:
            break
        if _lookup_distance(distance, q) < stop_rad:
            break
        out[n] = q
        n += 1
        length += step
        turn = turn_next
        prev_d = mid_d
        p = q
    return n, length, turn


@njit(cache=True)
def _stamp_distance(distance: np.ndarray, path: np.ndarray, radius: float) -> None:
    """Lower `distance` (radians to the nearest arrow, per occupancy cell) to each cell's
    angular distance from `path`, for cells within `radius` of it."""
    occ_h, occ_w = distance.shape
    d_lat = np.pi / occ_h
    d_lon = 2 * np.pi / occ_w
    cos_radius = np.cos(radius)
    for k in range(path.shape[0]):
        p = path[k]
        lat = np.arcsin(min(max(p[2], -1.0), 1.0))
        lon = np.arctan2(p[1], p[0])
        row_lo = max(int((np.pi / 2 - lat - radius) / d_lat), 0)
        row_hi = min(int((np.pi / 2 - lat + radius) / d_lat), occ_h - 1)
        for row in range(row_lo, row_hi + 1):
            cell_lat = np.pi / 2 - (row + 0.5) * d_lat
            cos_cell_lat = np.cos(cell_lat)
            if cos_cell_lat * np.pi < radius:
                col_lo, col_hi = 0, occ_w - 1
            else:
                half = min(radius / cos_cell_lat, np.pi)
                col_lo = int(np.floor((lon - half + np.pi) / d_lon))
                col_hi = int(np.floor((lon + half + np.pi) / d_lon))
            for c in range(col_lo, col_hi + 1):
                col = c % occ_w
                cell_lon = -np.pi + (col + 0.5) * d_lon
                cx = cos_cell_lat * np.cos(cell_lon)
                cy = cos_cell_lat * np.sin(cell_lon)
                cz = np.sin(cell_lat)
                dot = cx * p[0] + cy * p[1] + cz * p[2]
                if dot < cos_radius:
                    continue
                angle = np.arccos(min(dot, 1.0))
                if angle < distance[row, col]:
                    distance[row, col] = angle
