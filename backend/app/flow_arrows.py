"""Super arrows for the wind and ocean-current views (#805).

The views used to draw one short arrow per sampled cell, which shows the local field but not
its large-scale pattern: a trade-wind belt or a gyre reads as a field of ticks. Atlas maps
instead draw a few long, curved arrows that follow the prevailing flow. This module finds
those arrows in a vector field on the climate module's equirectangular grid. It changes no
model state. The arrows are an interpretation, not a plot: smoothed for legibility, and near
coasts drawing the model's coastal rule wider than the model applies it, so there they can
disagree with the computed field (see docs/simulation-model/climate.md).

The method follows the issue's sketch, made concrete as evenly-spaced streamline placement
(Jobard & Lefer 1997):

1. **Cluster.** Smooth the field over a few degrees and score each place by speed times the
   *coherence* of its neighbourhood -- the length of the mean vector over the mean speed,
   1 where every vector points the same way, 0 where they cancel. The highest-scoring places
   not yet claimed by an arrow are the centres of coherent clusters, and seed the arrows.
2. **Grow.** Each seed grows into a streamline, forward and then backward, one small step at
   a time, until it comes within `stop_deg` of an existing arrow, runs off the valid region
   (land, for currents), slows below the speed floor, has turned `max_turn_deg` net (signed,
   so a loop is cut but an S-bend is not) or `max_total_turn_deg` in all (so a squiggle is
   cut), or reaches `max_length_deg`. The turn limit splits a gyre into a few arcs, as atlas
   maps do.
3. **Link.** Points every `control_spacing_deg` along the streamline become the control
   points of a B-spline that only approaches them, so the arrow is a smooth sweep like a
   hand-drawn diagram rather than a copy of every wobble in the trace. Where that would cut
   a corner over land, the spacing tightens until it doesn't. Drawing (the tapered body and
   the head) is `render_image._draw_flow_arrows`.

**Coasts.** A current that meets land has to turn. The model turns it with
`climate._deflect_into_coast`: the part of the flow aimed into land is cancelled and sent
along the coast, clockwise around land in the northern hemisphere and anticlockwise in the
southern, at `climate.DEFLECTION_SPEEDUP_FACTOR` times its speed. It does so only in the
ocean cells touching land, a band too narrow to see at map scale or for an arrow to follow
(on present-day Earth the flow one cell further out still points into the coast), and the
smoothed field blurs it away. Within `coast_band_deg` of an invalid cell the tracer
therefore applies the same rule to its smoothed field, scaled from nothing at the band's
edge to the full rule `coast_hold_off_deg` from land, so arrows curve before they reach the
coast and then run parallel to it. The coast normal is the gradient of a smoothed proximity
to land, steady where the cell outline steps.

In the band the trace steps at most `coast_step_deg`, and control points are
`coast_control_spacing_deg` apart. The turned flow is slow more
often near land (on present-day Earth, 14% of cells within a degree of land fall below the
open-water speed floor, against 3% in open water), mostly where the redirected flow meets
flow already running the other way along the coast, so the floor falls toward
`coast_speed_floor_fraction` of itself across the band. Arrows from both sides reach the
same stretch of coast, so the stop distance falls toward `coast_stop_fraction` of itself
too, letting an arrow run along the shore past its neighbours' ends.

An arrow seeded in open water often uses up its length before it reaches a coast, and the
spacing rule then keeps other arrows from starting there, so no arrow would show a turn.
Seeding therefore starts at *impact points*, where the coastal rule turns at least
`min_impact_fraction` of the flow, strongest first. Traced backward, such an arrow runs out
to sea along the incoming current; traced forward, it follows the turned flow along the
coast. Those coastal stretches show where and which way the model's rule turns flow, not the
computed current there, which only turns in the cells touching land.

Accepted arrows stamp a distance field on a one-degree grid, so later seeds and traces see
them in constant time.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy import interpolate, ndimage
from scipy.spatial import cKDTree

from .hydroclimate import climate


@dataclass(frozen=True)
class FlowArrowParams:
    """Spacing and stopping rules for `trace_flow_arrows`. Angles are degrees of arc."""

    separation_deg: float = 9.0  # no new arrow starts closer than this to an existing one
    stop_deg: float = 4.0  # a growing arrow stops this close to an existing one
    step_deg: float = 0.75  # streamline integration step
    min_length_deg: float = 10.0  # shorter traces are dropped, not drawn
    max_length_deg: float = 42.0
    max_turn_deg: float = 130.0  # net turning along one arrow, either way
    max_total_turn_deg: float = 200.0  # turning summed without sign: caps wiggles
    smoothing_deg: float = 3.0  # Gaussian sigma for the prevailing field
    min_speed_fraction: float = 0.2  # speed floor, as a fraction of the 90th-percentile speed
    min_coherence: float = 0.6  # seeds need at least this neighbourhood coherence
    seed_spacing_deg: float = 3.0  # candidate seed lattice
    control_spacing_deg: float = 6.0  # spline control point spacing
    max_abs_lat_deg: float = 82.0  # east/north degenerate at the poles; stay clear of them
    coast_band_deg: float = 3.0  # how far from an invalid cell coastal tracing applies
    coast_step_deg: float = 0.2  # step at the coast, capped at 0.4 grid cells
    coast_control_spacing_deg: float = 5.0  # control point spacing at the coast
    min_control_spacing_deg: float = 1.5  # the tightest spacing tried to keep a spline off land
    coast_hold_off_deg: float = 1.5  # arrows run along the coast this far out
    coast_normal_smoothing_deg: float = 1.5  # Gaussian sigma for the coast normal
    coast_speed_floor_fraction: float = 0.25  # the speed floor at the coast, as a fraction of the open-water one
    coast_stop_fraction: float = 0.4  # `stop_deg` at the coast, as a fraction of the open-water one
    min_impact_fraction: float = 0.5  # an impact seed turns at least this share of the flow


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
    never cross an invalid cell. Non-finite cells count as invalid. Deterministic: the same
    field always gives the same arrows."""
    finite = np.isfinite(u) & np.isfinite(v)
    valid = valid & finite
    u, v = np.where(finite, u, 0.0), np.where(finite, v, 0.0)
    u_s, v_s, coherence = _prevailing_field(u, v, valid, params.smoothing_deg)
    speed = np.hypot(u_s, v_s)
    if not np.any(valid):
        return []
    reference_speed = float(np.percentile(speed[valid], 90))
    if reference_speed <= 0.0:
        return []
    min_speed = params.min_speed_fraction * reference_speed
    speed_floors = np.array([min_speed, min_speed * params.coast_speed_floor_fraction])

    # Near a coast, turn the flow along it the way the model does; see the module docstring's
    # "Coasts".
    coast = _coast_weight(valid, params.coast_band_deg)
    u_t, v_t, impact = _turn_along_coast(u_s, v_s, coast, valid, params)
    fields = np.stack([u_t, v_t, valid.astype(np.float64), coast])
    occ_h, occ_w = int(round(180 / _OCCUPANCY_RES_DEG)), int(round(360 / _OCCUPANCY_RES_DEG))
    distance = np.full((occ_h, occ_w), np.inf)
    coast_step_deg = min(params.coast_step_deg, 0.4 * 180.0 / u.shape[0])
    steps = np.radians([params.step_deg, coast_step_deg])
    max_points = int(np.ceil(params.max_length_deg / coast_step_deg)) + 2
    limits = np.radians([params.stop_deg, params.stop_deg * params.coast_stop_fraction, params.max_abs_lat_deg])
    max_length = np.radians(params.max_length_deg)
    turn_limits = np.radians([params.max_turn_deg, params.max_total_turn_deg])

    arrows: list[FlowArrow] = []
    for seed in _seed_candidates(speed, coherence, valid, impact, min_speed, params):
        if _lookup_distance(distance, seed) < np.radians(params.separation_deg):
            continue
        fwd, fwd_coast = np.empty((max_points, 3)), np.empty(max_points)
        bwd, bwd_coast = np.empty((max_points, 3)), np.empty(max_points)
        n_fwd, length, net, total = _trace(seed, 1.0, fields, distance, steps, max_length, turn_limits, limits, speed_floors, fwd, fwd_coast)
        # The backward half gets what the forward half left of each budget. That is exact for
        # a bend that carries on through the seed. An S-bend's two halves turn opposite ways
        # and would partly cancel, but each is charged in full: conservative, so it can end
        # an S-bend early but never exceeds a limit.
        n_bwd, length_b, _, _ = _trace(
            seed, -1.0, fields, distance, steps, max_length - length, turn_limits - np.array([net, total]), limits, speed_floors, bwd, bwd_coast,
        )
        if np.degrees(length + length_b) < params.min_length_deg:
            continue
        path = np.concatenate([bwd[:n_bwd][::-1], fwd[1:n_fwd]])
        path_coast = np.concatenate([bwd_coast[:n_bwd][::-1], fwd_coast[1:n_fwd]])
        _stamp_distance(distance, path, np.radians(params.separation_deg))
        arrows.append(FlowArrow(points_xyz=_spline_off_land(path, path_coast, valid, params), mean_speed=float(np.mean(_sample_speed(speed, path)))))
    return arrows


def _turn_along_coast(u: np.ndarray, v: np.ndarray, coast: np.ndarray, valid: np.ndarray, params: FlowArrowParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`climate._deflect_into_coast`'s rule applied across the coastal band: the part of
    (u, v) aimed into land is cancelled and redirected along the coast in the model's
    hemisphere sense at its speed-up. The rule's share grows from 0 at the band's edge to all
    of it `coast_hold_off_deg` from land, so arrows run parallel to the coast that far out
    rather than tracing the grid's stepped shoreline. The outward normal is minus the gradient
    of proximity to land (`coast`, with land counted as 1), smoothed over
    `coast_normal_smoothing_deg` so a stepped coastline gives a steady normal. Returns the
    turned (u, v) and the part of the flow that was turned, in speed units."""
    height, width = coast.shape
    lat = np.radians(90.0 - (np.arange(height) + 0.5) * (180.0 / height))[:, None]
    sigma = (params.coast_normal_smoothing_deg / (180.0 / height), params.coast_normal_smoothing_deg / (360.0 / width))
    proximity = ndimage.gaussian_filter(np.where(valid, coast, 1.0), sigma, mode=("nearest", "wrap"))
    d_north = -np.gradient(proximity, axis=0) / np.radians(180.0 / height)
    d_east = (np.roll(proximity, -1, axis=1) - np.roll(proximity, 1, axis=1)) / (2 * np.radians(360.0 / width) * np.maximum(np.cos(lat), 1e-3))
    norm = np.hypot(d_east, d_north)
    has_normal = valid & (coast > 0) & (norm > 1e-9)
    nx = np.where(has_normal, -d_east / np.where(has_normal, norm, 1.0), 0.0)
    ny = np.where(has_normal, -d_north / np.where(has_normal, norm, 1.0), 0.0)
    hemisphere = np.where(lat >= 0, 1.0, -1.0)
    # climate._rotate_90(nx, ny, hemisphere): the along-coast direction the model uses.
    tx, ty = -hemisphere * ny, hemisphere * nx
    share = np.clip(coast * params.coast_band_deg / max(params.coast_band_deg - params.coast_hold_off_deg, 1e-9), 0.0, 1.0)
    into = np.clip(-(u * nx + v * ny), 0.0, None) * share
    redirected = into * climate.DEFLECTION_SPEEDUP_FACTOR
    return u + into * nx + redirected * tx, v + into * ny + redirected * ty, into


def _coast_weight(valid: np.ndarray, band_deg: float) -> np.ndarray:
    """Per cell, 1 at an invalid cell's edge falling linearly to 0 at `band_deg` away from
    every invalid cell (true great-circle distance); 0 everywhere when all cells are valid."""
    height, width = valid.shape
    weight = np.zeros((height, width))
    if valid.all() or not valid.any():
        return weight
    # Distance is measured cell centre to cell centre, and half a cell is the gap from a
    # centre to the shore, so the search reaches that much past the band.
    half_cell_deg = 0.5 * 180.0 / height
    reach_deg = band_deg + half_cell_deg
    # Only cells within reach can get a weight; find them cheaply on the grid first. The
    # column reach is widened for latitude, but only up to 85 degrees: closer to the poles a
    # cell may get too little weight, which is harmless because arrows stop at
    # `max_abs_lat_deg` (82 degrees).
    reach_rows = int(np.ceil(reach_deg / (180.0 / height))) + 1
    lat = np.radians(90.0 - (np.arange(height) + 0.5) * (180.0 / height))
    lon = np.radians(-180.0 + (np.arange(width) + 0.5) * (360.0 / width))
    cos_floor = max(np.cos(np.radians(90.0 - reach_deg - 180.0 / height)), np.cos(np.radians(85.0)))
    reach_cols = min(int(np.ceil(reach_deg / (360.0 / width) / cos_floor)) + 1, width // 2)
    near = ndimage.maximum_filter(~valid, size=(2 * reach_rows + 1, 2 * reach_cols + 1), mode=("nearest", "wrap")) & valid
    xyz = np.stack(np.broadcast_arrays(np.cos(lat)[:, None] * np.cos(lon)[None, :], np.cos(lat)[:, None] * np.sin(lon)[None, :], np.sin(lat)[:, None]), axis=-1)
    tree = cKDTree(xyz[~valid])
    chord, _ = tree.query(xyz[near], distance_upper_bound=2.0 * np.sin(np.radians(reach_deg) / 2.0))
    # Past the bound the query returns inf: those cells keep weight 0.
    angle_deg = np.degrees(2.0 * np.arcsin(np.clip(np.where(np.isfinite(chord), chord, 2.0) / 2.0, 0.0, 1.0)))
    edge_deg = np.maximum(angle_deg - half_cell_deg, 0.0)
    weight[near] = np.clip(1.0 - edge_deg / band_deg, 0.0, 1.0)
    return weight


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


def _seed_candidates(speed: np.ndarray, coherence: np.ndarray, valid: np.ndarray, impact: np.ndarray, min_speed: float, params: FlowArrowParams) -> np.ndarray:
    """(K, 3) seed points. Impact points come first, strongest first: grid cells where the
    coastal rule turns at least `min_impact_fraction` of the flow (`impact` is the turned
    part, in speed units). They lie within about a cell of land, finer than the open-water
    lattice, so they are taken from the grid itself. Then a roughly equal-area lattice of
    open-water points, best cluster centre first."""
    height, width = speed.shape
    lat_cells = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_cells = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    is_impact = valid & (speed >= min_speed) & (impact >= params.min_impact_fraction * speed) & (np.abs(lat_cells)[:, None] < params.max_abs_lat_deg)
    impact_rows, impact_cols = np.nonzero(is_impact)
    # Stable sorts, so ties resolve by grid or lattice order and the result never depends on
    # numpy's choice of algorithm.
    first = np.argsort(-impact[impact_rows, impact_cols], kind="stable")
    seed_lat, seed_lon = [lat_cells[impact_rows[first]]], [lon_cells[impact_cols[first]]]

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
    order = np.argsort(-(speed[rows, cols] * coherence[rows, cols] ** 2)[keep], kind="stable")
    seed_lat.append(lat_deg[keep][order])
    seed_lon.append(lon_deg[keep][order])

    lat_r, lon_r = np.radians(np.concatenate(seed_lat)), np.radians(np.concatenate(seed_lon))
    return np.stack([np.cos(lat_r) * np.cos(lon_r), np.cos(lat_r) * np.sin(lon_r), np.sin(lat_r)], axis=-1)


def _spline_off_land(path: np.ndarray, path_coast: np.ndarray, valid: np.ndarray, params: FlowArrowParams) -> np.ndarray:
    """The smoothest spline through `path` that stays on valid cells: wide control spacing
    first, tighter until no sample lands on an invalid cell (a smoothed arc can cut across a
    headland the trace went around), down to `min_control_spacing_deg`."""
    scale = 1.0
    while True:
        spacing = max(params.control_spacing_deg * scale, params.min_control_spacing_deg)
        coast_spacing = max(params.coast_control_spacing_deg * scale, params.min_control_spacing_deg)
        curve = _smooth_spline(path[_control_indices(path, path_coast, spacing, coast_spacing)])
        if valid.all() or valid[_cells(valid.shape, curve)].all() or max(spacing, coast_spacing) <= params.min_control_spacing_deg:
            return curve
        scale *= 0.6


def _cells(shape: tuple[int, int], points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(rows, cols) of the grid cells holding `points` (N, 3) on an equirectangular grid."""
    height, width = shape
    lat = np.degrees(np.arcsin(np.clip(points[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(points[:, 1], points[:, 0]))
    rows = np.clip(((90.0 - lat) / 180.0 * height).astype(int), 0, height - 1)
    cols = np.clip(((lon + 180.0) / 360.0 * width).astype(int), 0, width - 1)
    return rows, cols


def _control_indices(path: np.ndarray, coast: np.ndarray, spacing_deg: float, coast_spacing_deg: float) -> np.ndarray:
    """Indices of spline control points along `path`: one every `spacing_deg` of arc in open
    water, closing to `coast_spacing_deg` by coastal weight, so a turn at a coast survives the
    spline. Both ends are kept; a last gap shorter than half the spacing merges into the end."""
    seg = np.degrees(np.arccos(np.clip(np.sum(path[1:] * path[:-1], axis=1), -1.0, 1.0)))
    target = spacing_deg + (coast_spacing_deg - spacing_deg) * coast
    indices = [0]
    run = 0.0
    for i in range(1, len(path)):
        run += seg[i - 1]
        if run >= target[i]:
            indices.append(i)
            run = 0.0
    if indices[-1] != len(path) - 1:
        if len(indices) > 1 and run < 0.5 * target[-1]:
            indices[-1] = len(path) - 1
        else:
            indices.append(len(path) - 1)
    return np.array(indices)


def _smooth_spline(control: np.ndarray, sample_deg: float = 0.5) -> np.ndarray:
    """Clamped cubic B-spline over `control` (K, 3), sampled about every `sample_deg` and
    re-projected onto the sphere. It starts and ends at the end control points but only
    approaches the others, so a wobble in the trace becomes a gentle bend rather than a kink,
    as in hand-drawn flow diagrams (#805). The cost is that the drawn path can cut inside a
    tight bend of the traced flow."""
    count = len(control)
    if count < 3:
        return control.copy()
    degree = min(3, count - 1)
    interior = np.arange(1, count - degree, dtype=float)
    knots = np.concatenate([np.zeros(degree + 1), interior, np.full(degree + 1, count - degree)])
    spline = interpolate.BSpline(knots, control, degree)
    # Evaluate densely, then resample evenly by arc length: the spline's parameter runs
    # unevenly along the curve, fastest near the clamped ends.
    polygon_deg = np.degrees(np.sum(np.arccos(np.clip(np.sum(control[1:] * control[:-1], axis=1), -1.0, 1.0))))
    dense = spline(np.linspace(0.0, count - degree, 8 * max(int(np.ceil(polygon_deg / sample_deg)), 2) + 1))
    dense /= np.linalg.norm(dense, axis=1, keepdims=True)
    arc = np.concatenate([[0.0], np.cumsum(np.arccos(np.clip(np.sum(dense[1:] * dense[:-1], axis=1), -1.0, 1.0)))])
    targets = np.linspace(0.0, arc[-1], max(int(np.ceil(np.degrees(arc[-1]) / sample_deg)), 1) + 1)
    curve = np.stack([np.interp(targets, arc, dense[:, k]) for k in range(3)], axis=-1)
    return curve / np.linalg.norm(curve, axis=1, keepdims=True)


def _sample_speed(speed: np.ndarray, points: np.ndarray) -> np.ndarray:
    return speed[_cells(speed.shape, points)]


@njit(cache=True)
def _lookup_distance(distance: np.ndarray, p: np.ndarray) -> float:
    occ_h, occ_w = distance.shape
    lat = np.arcsin(min(max(p[2], -1.0), 1.0))
    lon = np.arctan2(p[1], p[0])
    row = min(max(int((np.pi / 2 - lat) / np.pi * occ_h), 0), occ_h - 1)
    col = min(max(int((lon + np.pi) / (2 * np.pi) * occ_w), 0), occ_w - 1)
    return distance[row, col]


@njit(cache=True)
def _bilinear(field: np.ndarray, lat: float, lon: float) -> float:
    """`field` (H, W) at (lat, lon) in radians, bilinear, longitude wrapping."""
    height, width = field.shape
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
    return (1 - tr) * ((1 - tc) * field[r0, c0] + tc * field[r0, c1]) + tr * ((1 - tc) * field[r1, c0] + tc * field[r1, c1])


@njit(cache=True)
def _latlon(p: np.ndarray) -> tuple[float, float]:
    return np.arcsin(min(max(p[2], -1.0), 1.0)), np.arctan2(p[1], p[0])


@njit(cache=True)
def _is_valid(fields: np.ndarray, p: np.ndarray) -> bool:
    lat, lon = _latlon(p)
    return _bilinear(fields[2], lat, lon) >= 0.5


@njit(cache=True)
def _direction(fields: np.ndarray, p: np.ndarray, sign: float, min_speed: float):
    """(unit flow direction at `p` as a 3-vector times `sign`, coastal weight there). The
    direction is zeros where the flow is invalid or slower than `min_speed`. `fields` is
    (u, v, valid, coast) stacked; `min_speed` is (open-water floor, floor at the coast)."""
    lat, lon = _latlon(p)
    su, sv = _bilinear(fields[0], lat, lon), _bilinear(fields[1], lat, lon)
    coast = _bilinear(fields[3], lat, lon)
    out = np.zeros(3)
    speed = np.hypot(su, sv)
    if _bilinear(fields[2], lat, lon) < 0.5 or speed < min_speed[0] + (min_speed[1] - min_speed[0]) * coast:
        return out, coast
    sin_lat, cos_lat, sin_lon, cos_lon = np.sin(lat), np.cos(lat), np.sin(lon), np.cos(lon)
    # east = (-sin lon, cos lon, 0); north = (-sin lat cos lon, -sin lat sin lon, cos lat)
    out[0] = sign * (su * -sin_lon + sv * -sin_lat * cos_lon) / speed
    out[1] = sign * (su * cos_lon + sv * -sin_lat * sin_lon) / speed
    out[2] = sign * (sv * cos_lat) / speed
    return out, coast


@njit(cache=True)
def _is_zero(d: np.ndarray) -> bool:
    return d[0] == 0.0 and d[1] == 0.0 and d[2] == 0.0


@njit(cache=True)
def _advance(p: np.ndarray, d: np.ndarray, h: float) -> np.ndarray:
    q = p + h * d
    return q / np.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])


@njit(cache=True)
def _trace(seed, sign, fields, distance, steps, max_length, turn_limits, limits, min_speed, out, out_coast):
    """Grow a streamline from `seed` along `sign` times the flow, writing points into `out`
    and each point's coastal weight into `out_coast`. `steps` is (open-water step, coastal
    step), `turn_limits` is (net turning, turning summed without sign) and `limits` is
    (open-water stop distance, coastal stop distance, maximum |latitude|), all radians.
    Returns (point count, arc length, |net turning|, turning summed without sign)."""
    out[0] = seed
    prev_d, coast = _direction(fields, seed, sign, min_speed)
    out_coast[0] = coast
    n, length, turn, total = 1, 0.0, 0.0, 0.0
    if _is_zero(prev_d):
        return n, length, turn, total
    p = seed.copy()
    while n < out.shape[0]:
        d, coast = _direction(fields, p, sign, min_speed)
        h = steps[0] + (steps[1] - steps[0]) * coast
        if _is_zero(d) or length + h > max_length + 1e-12:
            break
        used, _ = _direction(fields, _advance(p, d, 0.5 * h), sign, min_speed)
        if _is_zero(used):
            break
        q = _advance(p, used, h)
        if not _is_valid(fields, q):
            break
        # Signed about the local vertical, so a zigzag cancels and only net turning counts.
        cross = np.cross(prev_d, used)
        angle = np.arctan2(np.sum(cross * q), np.sum(prev_d * used))
        turn_next, total_next = turn + angle, total + abs(angle)
        stop = limits[0] + (limits[1] - limits[0]) * coast
        if abs(turn_next) > turn_limits[0] or total_next > turn_limits[1]:
            break
        if abs(_latlon(q)[0]) > limits[2] or _lookup_distance(distance, q) < stop:
            break
        out[n] = q
        out_coast[n] = coast
        n += 1
        length += h
        turn, total = turn_next, total_next
        prev_d = used
        p = q
    return n, length, abs(turn), total


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
        # The widest the cap around p reaches in longitude, on any row: arcsin(sin r / cos lat),
        # or all the way round once the cap covers a pole.
        if abs(lat) + radius >= np.pi / 2:
            half = np.pi
        else:
            half = np.arcsin(min(np.sin(radius) / np.cos(lat), 1.0))
        if half >= np.pi:
            col_lo, col_hi = 0, occ_w - 1
        else:
            col_lo = int(np.floor((lon - half + np.pi) / d_lon))
            col_hi = int(np.floor((lon + half + np.pi) / d_lon))
        for row in range(row_lo, row_hi + 1):
            cell_lat = np.pi / 2 - (row + 0.5) * d_lat
            cos_cell_lat = np.cos(cell_lat)
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
