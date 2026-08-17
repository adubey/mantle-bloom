"""Server-side rasterization for `GET /world/render`: bakes a requested "kind of map"
(elevation fill, plate-boundary/pole/velocity overlay, or raw per-plate node dots) directly
into a PNG at the client's requested resolution, rather than shipping the underlying
projected coordinates as JSON and letting the frontend draw them onto a `<canvas>` (the old
design -- see git history / docs/api-reference.md for what that wire format used to look
like). This module is a close port of what `frontend/src/MapCanvas.tsx` used to compute
client-side: same transform, same per-view drawing rules, just onto a Pillow image instead
of a 2D canvas context.

Moving rendering server-side also lifts a constraint the old JSON format was stuck with: the
render grid used to be deliberately coarsened (see GRID_SPACING_KM below) purely to keep the
wire payload small, since every grid point was serialized as JSON numbers. A PNG's size
depends on how compressible its *pixels* are, not on how many simulation samples went into
it, so there's no payload-size reason to keep the grid coarse -- it's a free-standing
"looks good rendered" choice now, independent of resolution requested.
"""

from __future__ import annotations

import base64
import io

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

from . import climate, geometry, mantle, plates, projections
from .world import World

# Climate views draw from climate.py's own fixed (H, W) grid, not the render grid below --
# see climate.py's module docstring for why. Handled by a separate code path
# (_render_climate_view) rather than threading a third data source through render_png's
# existing elevation/plates machinery.
CLIMATE_VIEWS = ("temperature", "wind", "oceanCurrents", "humidity", "precipitation")
VIEWS = ("elevation", "plates", "platesDetail") + CLIMATE_VIEWS

BACKGROUND_RGB = (11, 16, 32)  # #0b1020

# Visual constants below are all in pixel terms tuned at this reference width; render_png
# scales them by (requested width / REFERENCE_WIDTH_PX) so a higher-resolution request (e.g.
# doubling it for a sharper, same-*displayed*-size map) doesn't also make lines/dots/padding
# look proportionally thinner -- see render_png's `px` scale factor.
REFERENCE_WIDTH_PX = 1100
PADDING_PX = 20
POLE_RADIUS_PX = 5
ARROWHEAD_LENGTH_PX = 7
NODE_DOT_RADIUS_PX = 1.6
BOUNDARY_LINE_WIDTH_PX = 1.5
ARC_LINE_WIDTH_PX = 1.5
# Cells are drawn slightly larger than their measured half-extent so adjacent cells overlap
# a hair rather than risk a hairline gap from floating-point rounding. A ratio, not a pixel
# size, so it does not scale with requested resolution.
CELL_OVERLAP_FACTOR = 1.15
# The rotation-rate indicator (see _draw_rotation_arc) is a fixed-pixel-radius arc drawn
# around the pole marker itself, its length (degrees of arc) scaled between these two bounds
# by how fast the plate is moving relative to mantle.MAX_PLATE_RATE -- always at least a
# little visible, never more than most of a full circle.
ARC_RADIUS_PX = POLE_RADIUS_PX * 3.0
ARC_MIN_EXTENT_DEG = 40.0
ARC_MAX_EXTENT_DEG = 300.0
# How far (radians) from the displayed pole to sample a point for measuring the arc's true
# sweep direction once projected (see _draw_rotation_arc) -- small enough to stay a good
# local approximation of the tangent at the pole itself, far enough that its projection
# doesn't degenerate to a near-zero-length step.
ARC_DIRECTION_SAMPLE_RAD = 0.05

# Resolution of the render grid (see _render_grid_arrays), swept on a plate-independent
# global grid so the map's coverage never depends on how sparse any one plate's own line
# data looks once projected -- a fixed, display-oriented constant, deliberately *not* tied
# to plates.TARGET_LINE_SPACING_RAD (the simulation's physics resolution): the render grid
# only needs to look smooth once rasterized, which a resolution change in the physics has no
# bearing on. 100km (rather than a coarser value) is affordable now that this is a fully
# server-side numpy-slice fill (see _fill_rects) instead of one JSON number per point on the
# wire: ~117ms -> ~220ms per 2200x1222 render, for visibly smoother coastlines (confirmed
# side by side -- 250km left blocky, stair-stepped edges even at this canvas size).
GRID_SPACING_KM = 100.0
GRID_SPACING_RAD = GRID_SPACING_KM / plates.PLANET_RADIUS_KM

# A fixed categorical palette so each plate reads as a distinct region across
# generate/step calls (plate_id is stable within one world's lifetime).
PLATE_PALETTE = np.array(
    [
        (230, 25, 75), (60, 180, 75), (255, 225, 25), (67, 99, 216), (245, 130, 49),
        (66, 212, 244), (240, 50, 230), (188, 246, 12), (250, 190, 190), (70, 153, 144),
        (230, 190, 255), (154, 99, 36), (255, 250, 200), (128, 0, 0), (170, 255, 195),
        (128, 128, 0), (255, 216, 177), (0, 0, 117), (169, 169, 169), (255, 255, 255),
    ],
    dtype=np.uint8,
)

# A simple hypsometric tint: deep blue -> shallow blue -> coastal green -> tan -> brown ->
# white peaks. Kept in sync by hand with the palette this was ported from
# (frontend/src/elevationColor.ts, now deleted -- this is the sole copy).
_ELEVATION_STOP_E = np.array([-11000, -4000, -1500, -200, 0, 200, 1200, 3000, 6000, 9000], dtype=float)
_ELEVATION_STOP_RGB = np.array(
    [
        (10, 10, 40), (15, 40, 110), (40, 110, 190), (110, 170, 210), (200, 210, 150),
        (90, 150, 60), (170, 160, 90), (120, 90, 60), (230, 230, 230), (255, 255, 255),
    ],
    dtype=float,
)


def elevation_colors(elevations: np.ndarray) -> np.ndarray:
    """Vectorized elevation -> RGB: piecewise-linear interpolation through the hypsometric
    stops above, clamped at the ends (numpy.interp's default behavior already matches the
    old manual clamp)."""
    channels = [np.interp(elevations, _ELEVATION_STOP_E, _ELEVATION_STOP_RGB[:, c]) for c in range(3)]
    return np.clip(np.round(np.stack(channels, axis=-1)), 0, 255).astype(np.uint8)


# Diverging cold -> hot, degrees C. Covers climate.py's land/air temperature range
# comfortably (LAND_TEMP_MIN_C=-60 to LAND_TEMP_MIN_C+LAND_TEMP_RANGE_C=35).
_TEMPERATURE_STOP_C = np.array([-60, -30, -10, 0, 10, 20, 30, 40], dtype=float)
_TEMPERATURE_STOP_RGB = np.array(
    [
        (30, 20, 90), (40, 60, 170), (70, 130, 220), (150, 200, 230),
        (230, 230, 140), (240, 170, 60), (210, 80, 40), (140, 20, 20),
    ],
    dtype=float,
)

# Dry (tan) -> humid (deep teal/blue). climate.py's humidity is roughly [0, MAX_EVAPORATION_CEILING=1.4].
_HUMIDITY_STOP_V = np.array([0.0, 0.3, 0.6, 0.9, 1.2, 1.4], dtype=float)
_HUMIDITY_STOP_RGB = np.array(
    [
        (120, 100, 60), (170, 150, 90), (140, 170, 110), (70, 150, 140), (30, 100, 140), (15, 60, 110),
    ],
    dtype=float,
)

# Dry (tan) -> wet (dark blue), mm/year. climate.py's precipitation_mm typically runs 0-2000+.
_PRECIPITATION_STOP_MM = np.array([0, 250, 600, 1200, 2000, 3000], dtype=float)
_PRECIPITATION_STOP_RGB = np.array(
    [
        (180, 160, 100), (200, 190, 110), (140, 180, 90), (60, 140, 90), (40, 90, 150), (20, 40, 120),
    ],
    dtype=float,
)


def _interp_colors(values: np.ndarray, stops: np.ndarray, stop_rgb: np.ndarray) -> np.ndarray:
    channels = [np.interp(values, stops, stop_rgb[:, c]) for c in range(3)]
    return np.clip(np.round(np.stack(channels, axis=-1)), 0, 255).astype(np.uint8)


def temperature_colors(celsius: np.ndarray) -> np.ndarray:
    return _interp_colors(celsius, _TEMPERATURE_STOP_C, _TEMPERATURE_STOP_RGB)


def humidity_colors(humidity: np.ndarray) -> np.ndarray:
    return _interp_colors(humidity, _HUMIDITY_STOP_V, _HUMIDITY_STOP_RGB)


def precipitation_colors(precipitation_mm: np.ndarray) -> np.ndarray:
    return _interp_colors(precipitation_mm, _PRECIPITATION_STOP_MM, _PRECIPITATION_STOP_RGB)


def plate_colors(plate_ids: np.ndarray) -> np.ndarray:
    return PLATE_PALETTE[np.asarray(plate_ids, dtype=int) % len(PLATE_PALETTE)]


def _project_points(projection: str, world_pts: np.ndarray) -> np.ndarray:
    """World unit vectors (N, 3) -> projected (x, y), shape (N, 2)."""
    lat, lon = geometry.xyz_to_latlon(world_pts)
    x, y = projections.project(projection, lat, lon)
    return np.stack([x, y], axis=-1)


def _collect_all_points(world: World) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Every plate's current elevation-node positions, elevations, and owning plate_id,
    concatenated -- the source data the render grid resamples from."""
    points_list, elevation_list, owner_list = [], [], []
    for plate in world.plates:
        pts, elev = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        points_list.append(pts)
        elevation_list.append(elev)
        owner_list.append(np.full(len(pts), plate.plate_id))
    if not points_list:
        return None
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elevation_list, axis=0),
        np.concatenate(owner_list, axis=0),
    )


def _row_cell_half_extent(projection: str, phi: float, dtheta: float) -> tuple[float, float]:
    """How far apart (in *projected* units) this row's grid points end up, in each
    direction. See docs/simulation-model.md#render-image for the full derivation -- measured
    directly from the projection's local behavior (Behrmann, for instance, stretches
    longitude spacing near the poles) rather than assumed, so cells drawn at this size tile
    the map with no gaps anywhere, including at the poles."""
    origin = geometry.local_xyz(np.array([phi]), np.array([0.0]))
    theta_neighbor = geometry.local_xyz(np.array([phi]), np.array([dtheta]))
    phi_neighbor = geometry.local_xyz(np.array([phi + GRID_SPACING_RAD]), np.array([0.0]))
    (ox, oy), (tx, ty), (px, py) = (
        _project_points(projection, origin)[0],
        _project_points(projection, theta_neighbor)[0],
        _project_points(projection, phi_neighbor)[0],
    )
    return abs(tx - ox) / 2, abs(py - oy) / 2


def _render_grid_arrays(
    world: World, projection: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """A uniform lat/lon grid covering the whole sphere (GRID_SPACING_RAD, independent of
    any plate's own line spacing), each cell assigned its nearest elevation node's elevation
    and owning plate -- see docs/simulation-model.md#render-image. Returns flat concatenated
    (projected_xy, elevation, plate_id, cell_half_width, cell_half_height) arrays, or None
    for an empty world."""
    collected = _collect_all_points(world)
    if collected is None:
        return None
    all_points, all_elevation, all_owner = collected
    tree = cKDTree(all_points)

    xy_chunks, elev_chunks, owner_chunks, hw_chunks, hh_chunks = [], [], [], [], []
    for phi, theta_candidates, world_pts in plates.iter_local_lattice(np.eye(3), spacing_rad=GRID_SPACING_RAD):
        _, idx = tree.query(world_pts)
        xy_chunks.append(_project_points(projection, world_pts))
        elev_chunks.append(all_elevation[idx])
        owner_chunks.append(all_owner[idx])

        dtheta = GRID_SPACING_RAD / max(np.cos(phi), 1e-3)
        half_w, half_h = _row_cell_half_extent(projection, phi, dtheta)
        n = len(theta_candidates)
        hw_chunks.append(np.full(n, half_w))
        hh_chunks.append(np.full(n, half_h))

    return (
        np.concatenate(xy_chunks, axis=0),
        np.concatenate(elev_chunks, axis=0),
        np.concatenate(owner_chunks, axis=0),
        np.concatenate(hw_chunks, axis=0),
        np.concatenate(hh_chunks, axis=0),
    )


def _plate_tectonics(projection: str, plate) -> dict:
    """Pole marker, rotation arc, and boundary outline for a plate -- everything the
    "Plates"/"Plates (details)" views draw besides the elevation-fill/node dots."""
    speed = float(np.linalg.norm(plate.omega))

    pole = None
    rotation_arc = None
    if speed > 1e-15:
        # The true Euler pole -- real rotation axes are frequently nowhere near the plate
        # they belong to (this is physically normal, not a bug), so this can land anywhere
        # on the map. The pole marker is colored by plate (see render_png) specifically so
        # it still reads as "belonging to" the right plate even when it's drawn far away.
        pole_xyz = plate.omega / speed
        pole = _project_points(projection, pole_xyz[None, :])[0]

        intensity = np.clip(speed / mantle.MAX_PLATE_RATE, 0.3, 1.0)
        rotation_arc = {
            "pole_xyz": pole_xyz,
            "omega": plate.omega,
            "extent_deg": ARC_MIN_EXTENT_DEG + intensity * (ARC_MAX_EXTENT_DEG - ARC_MIN_EXTENT_DEG),
        }

    outline_world = plate.outline_world()
    boundary = _project_points(projection, outline_world) if len(outline_world) > 0 else np.zeros((0, 2))

    return {"pole": pole, "rotation_arc": rotation_arc, "boundary": boundary}


def _to_pixels(scale: float, offset_x: float, offset_y: float, xy: np.ndarray) -> np.ndarray:
    """(N, 2) projected coords -> (N, 2) pixel coords. y is flipped: projected y grows
    "up" (north), image y grows down, same as MapCanvas.tsx's old toPixel."""
    px = scale * xy[:, 0] + offset_x
    py = -scale * xy[:, 1] + offset_y
    return np.stack([px, py], axis=-1)


def _fill_rects(pixels: np.ndarray, centers: np.ndarray, half_w, half_h, colors: np.ndarray) -> None:
    """Fills each point's axis-aligned rectangle (centers[i] +/- (half_w, half_h)) directly
    via numpy slicing. Used for both the render-grid cells and (in "Plates (details)") the
    per-node dots -- at tens of thousands of points, one array slice per point is far
    cheaper than one Pillow draw call per point."""
    if len(centers) == 0:
        return
    height, width, _ = pixels.shape
    half_w = np.broadcast_to(np.asarray(half_w, dtype=float), (len(centers),))
    half_h = np.broadcast_to(np.asarray(half_h, dtype=float), (len(centers),))
    x0 = np.clip(np.round(centers[:, 0] - half_w).astype(int), 0, width)
    x1 = np.clip(np.round(centers[:, 0] + half_w).astype(int), 0, width)
    y0 = np.clip(np.round(centers[:, 1] - half_h).astype(int), 0, height)
    y1 = np.clip(np.round(centers[:, 1] + half_h).astype(int), 0, height)
    for i in range(len(centers)):
        if x1[i] > x0[i] and y1[i] > y0[i]:
            pixels[y0[i] : y1[i], x0[i] : x1[i]] = colors[i]


def _stroke_robust_loop(draw: ImageDraw.ImageDraw, pixel_pts: np.ndarray, color: tuple[int, int, int], width_px: float) -> None:
    """A plate's live outline (Plate.outline_world) is built from many short consecutive
    hops, but a plate whose local frame happens to put its own pole -- or the map
    projection's antimeridian -- inside its territory can still produce one or two segments
    that are wildly longer than the rest once projected. Rather than draw those as a
    straight line across the map, break the path there: find the typical (median) segment
    length for this loop and skip drawing anything far longer than that."""
    n = len(pixel_pts)
    if n < 2:
        return
    next_pts = np.roll(pixel_pts, -1, axis=0)
    seg_lengths = np.hypot(*(next_pts - pixel_pts).T)
    median = float(np.median(seg_lengths))
    break_threshold = max(median * 6, 20)
    line_width = max(int(round(width_px)), 1)
    for i in range(n):
        if seg_lengths[i] <= break_threshold:
            draw.line([tuple(pixel_pts[i]), tuple(next_pts[i])], fill=color, width=line_width)


def _draw_arrow_head(draw: ImageDraw.ImageDraw, tip: np.ndarray, direction: np.ndarray, color, length_px: float) -> None:
    """A small filled triangle at `tip`, pointing along the (already unit-length)
    `direction`, in pixel space."""
    angle = np.arctan2(direction[1], direction[0])
    p1 = tip - length_px * np.array([np.cos(angle - np.pi / 6), np.sin(angle - np.pi / 6)])
    p2 = tip - length_px * np.array([np.cos(angle + np.pi / 6), np.sin(angle + np.pi / 6)])
    draw.polygon([tuple(tip), tuple(p1), tuple(p2)], fill=color)


def _draw_rotation_arc(
    draw: ImageDraw.ImageDraw,
    projection: str,
    arc_info: dict,
    scale: float,
    offset_x: float,
    offset_y: float,
    pixel_scale: float,
) -> None:
    """A fixed-radius arc around the pole marker, swept in the plate's actual rotational
    direction (see below) by an angle representing its rotation rate (arc_info["extent_deg"],
    set in _plate_tectonics), with an arrowhead at the moving end -- replaces the old
    straight velocity arrow from the plate's seed point. Centering it on the pole (a single
    already-projected point, unlike the old arrow's separately-projected endpoint) also
    means it can no longer straddle a projection discontinuity the way the old arrow
    occasionally did near the antimeridian.

    **Sweep direction.** Image y grows down and projections aren't guaranteed to preserve
    on-screen handedness, so "clockwise" can't just be assumed from the sign of omega --
    it's measured directly: take a point near the pole, find its true tangential velocity
    (omega x point, the same formula boundary.py's closing_rate uses elsewhere), project
    both the point and a small step along that velocity into pixel space, and see whether
    the angle (around the pole, in PIL's arc-angle convention) increased or decreased.
    """
    pole_xyz = arc_info["pole_xyz"]
    omega = arc_info["omega"]
    extent_deg = arc_info["extent_deg"]

    pole_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, pole_xyz[None, :]))[0]

    helper = np.array([1.0, 0.0, 0.0]) if abs(pole_xyz[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    tangent_dir = geometry.normalize(np.cross(pole_xyz, helper))
    near_point = np.cos(ARC_DIRECTION_SAMPLE_RAD) * pole_xyz + np.sin(ARC_DIRECTION_SAMPLE_RAD) * tangent_dir
    velocity = np.cross(omega, near_point)
    if np.linalg.norm(velocity) < 1e-15:
        return
    step_point = geometry.normalize(near_point + 1e-4 * velocity)

    near_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, near_point[None, :]))[0]
    step_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, step_point[None, :]))[0]

    angle_near = np.degrees(np.arctan2(near_px[1] - pole_px[1], near_px[0] - pole_px[0]))
    angle_step = np.degrees(np.arctan2(step_px[1] - pole_px[1], step_px[0] - pole_px[0]))
    delta = ((angle_step - angle_near + 180) % 360) - 180  # signed shortest angular step

    radius = ARC_RADIUS_PX * pixel_scale
    if delta >= 0:
        start_deg, end_deg, head_angle, tangent_sign = angle_near, angle_near + extent_deg, angle_near + extent_deg, 1.0
    else:
        start_deg, end_deg, head_angle, tangent_sign = angle_near - extent_deg, angle_near, angle_near - extent_deg, -1.0

    color = (255, 255, 255)
    bbox = [pole_px[0] - radius, pole_px[1] - radius, pole_px[0] + radius, pole_px[1] + radius]
    draw.arc(bbox, start_deg, end_deg, fill=color, width=max(int(round(ARC_LINE_WIDTH_PX * pixel_scale)), 1))

    head_rad = np.radians(head_angle)
    head_px = pole_px + radius * np.array([np.cos(head_rad), np.sin(head_rad)])
    tangent_px = tangent_sign * np.array([-np.sin(head_rad), np.cos(head_rad)])
    _draw_arrow_head(draw, head_px, tangent_px, color, ARROWHEAD_LENGTH_PX * pixel_scale)


# Arrows are drawn at a coarser subsample of climate.py's own (90x180) computation grid --
# one arrow per cell would be unreadable clutter at these sizes.
ARROW_GRID_STRIDE = 6
ARROW_MAX_LENGTH_PX = 14.0
ARROW_LINE_WIDTH_PX = 1.3
SWELL_MARKER_RADIUS_PX = 4.0
WIND_ARROW_COLOR = (230, 230, 255)
CURRENT_ARROW_COLOR = (140, 210, 255)
CLIMATE_OCEAN_BACKDROP_RGB = np.array([18, 28, 55], dtype=np.uint8)
CLIMATE_LAND_BACKDROP_RGB = np.array([40, 46, 34], dtype=np.uint8)


def _climate_row_cell_half_extent(projection: str, lat_deg: float, dlat_deg: float, dlon_deg: float) -> tuple[float, float]:
    """Same measurement technique as _row_cell_half_extent, adapted to climate.py's plain
    fixed-dlat/dlon equirectangular grid rather than the render grid's cos(phi)-reduced one."""
    lat_r = np.radians(lat_deg)
    origin = geometry.latlon_to_xyz(np.array([lat_r]), np.array([0.0]))
    lon_neighbor = geometry.latlon_to_xyz(np.array([lat_r]), np.array([np.radians(dlon_deg)]))
    lat_neighbor = geometry.latlon_to_xyz(np.array([lat_r + np.radians(dlat_deg)]), np.array([0.0]))
    (ox, oy), (tx, ty), (px, py) = (
        _project_points(projection, origin)[0],
        _project_points(projection, lon_neighbor)[0],
        _project_points(projection, lat_neighbor)[0],
    )
    return abs(tx - ox) / 2, abs(py - oy) / 2


def _draw_climate_vectors(
    draw: ImageDraw.ImageDraw,
    fields: "climate.ClimateFields",
    u: np.ndarray,
    v: np.ndarray,
    projection: str,
    scale: float,
    offset_x: float,
    offset_y: float,
    pixel_scale: float,
    color: tuple[int, int, int],
) -> None:
    """Draws one arrow per subsampled grid cell, length scaled by speed relative to the max
    speed present. Direction is found the same way as the plate rotation arc's tangent: a
    small offset along the local (east, north) vector, projected, rather than assuming the
    projection preserves on-screen angles."""
    grid_h, grid_w = u.shape
    rows = np.arange(0, grid_h, ARROW_GRID_STRIDE)
    cols = np.arange(0, grid_w, ARROW_GRID_STRIDE)
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    rr, cc = rr.reshape(-1), cc.reshape(-1)

    base_xyz = fields.world_xyz[rr, cc]
    u_pts, v_pts = u[rr, cc], v[rr, cc]
    speed = np.hypot(u_pts, v_pts)
    keep = speed > 1e-6
    if not np.any(keep):
        return
    base_xyz, u_pts, v_pts, speed = base_xyz[keep], u_pts[keep], v_pts[keep], speed[keep]
    max_speed = float(speed.max())

    lon = np.arctan2(base_xyz[:, 1], base_xyz[:, 0])
    east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], axis=-1)
    north = np.cross(base_xyz, east)
    direction = geometry.normalize(u_pts[:, None] * east + v_pts[:, None] * north)
    offset_xyz = geometry.normalize(base_xyz + 0.02 * direction)

    base_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, base_xyz))
    offset_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, offset_xyz))

    arrow_vec = offset_px - base_px
    arrow_len = np.hypot(arrow_vec[:, 0], arrow_vec[:, 1])
    line_width = max(int(round(ARROW_LINE_WIDTH_PX * pixel_scale)), 1)
    head_len = ARROWHEAD_LENGTH_PX * 0.6 * pixel_scale

    for i in range(len(base_px)):
        if arrow_len[i] < 1e-6:
            continue
        dir_i = arrow_vec[i] / arrow_len[i]
        length_px = ARROW_MAX_LENGTH_PX * pixel_scale * (0.3 + 0.7 * speed[i] / max_speed)
        tip = base_px[i] + dir_i * length_px
        draw.line([tuple(base_px[i]), tuple(tip)], fill=color, width=line_width)
        _draw_arrow_head(draw, tip, dir_i, color, head_len)


def _draw_swell_markers(
    draw: ImageDraw.ImageDraw, fields: "climate.ClimateFields", projection: str, scale: float, offset_x: float, offset_y: float, pixel_scale: float,
) -> None:
    if len(fields.swell_rows) == 0:
        return
    xyz = fields.world_xyz[fields.swell_rows, fields.swell_cols]
    centers = _to_pixels(scale, offset_x, offset_y, _project_points(projection, xyz))
    r = SWELL_MARKER_RADIUS_PX * pixel_scale
    width_px = max(int(round(pixel_scale)), 1)
    for px, py in centers:
        draw.ellipse([px - r, py - r, px + r, py + r], outline=(255, 255, 255), width=width_px)


def _render_climate_view(world: World, projection: str, view: str, width: int, height: int) -> bytes:
    """Renders one of CLIMATE_VIEWS from climate.py's own fixed grid -- a separate path from
    the plate-tectonics views below since the data source (a real (H, W) array, always
    covering the whole sphere) is structurally different from the render grid's ragged
    lattice, so there's little to share beyond the pixel-space primitives."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    fields = climate.compute_climate(world)
    grid_h, grid_w = fields.elevation_m.shape
    flat_xyz = fields.world_xyz.reshape(-1, 3)
    flat_xy = _project_points(projection, flat_xyz)

    min_x, min_y = flat_xy.min(axis=0)
    max_x, max_y = flat_xy.max(axis=0)
    data_w = max(max_x - min_x, 1e-9)
    data_h = max(max_y - min_y, 1e-9)
    scale = min((width - 2 * padding_px) / data_w, (height - 2 * padding_px) / data_h)
    offset_x = width / 2 - scale * (min_x + max_x) / 2
    offset_y = height / 2 + scale * (min_y + max_y) / 2

    pixels = blank.copy()
    centers = _to_pixels(scale, offset_x, offset_y, flat_xy)

    dlat_deg = 180.0 / grid_h
    dlon_deg = 360.0 / grid_w
    half_w_by_row = np.empty(grid_h)
    half_h_by_row = np.empty(grid_h)
    for r in range(grid_h):
        half_w_by_row[r], half_h_by_row[r] = _climate_row_cell_half_extent(projection, float(fields.lat_deg[r]), dlat_deg, dlon_deg)
    half_w = np.repeat(half_w_by_row, grid_w) * scale * CELL_OVERLAP_FACTOR
    half_h = np.repeat(half_h_by_row, grid_w) * scale * CELL_OVERLAP_FACTOR

    if view == "temperature":
        # Whichever temperature is physically meaningful at that cell: ocean surface where
        # there's ocean, air temperature (already moderated toward nearby ocean) over land.
        display_temp = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
        colors = temperature_colors(display_temp.reshape(-1))
        _fill_rects(pixels, centers, half_w, half_h, colors)
    elif view == "humidity":
        colors = humidity_colors(fields.humidity.reshape(-1))
        _fill_rects(pixels, centers, half_w, half_h, colors)
    elif view == "precipitation":
        colors = precipitation_colors(fields.precipitation_mm.reshape(-1))
        _fill_rects(pixels, centers, half_w, half_h, colors)
    elif view in ("wind", "oceanCurrents"):
        backdrop = np.where(fields.is_ocean.reshape(-1)[:, None], CLIMATE_OCEAN_BACKDROP_RGB, CLIMATE_LAND_BACKDROP_RGB)
        _fill_rects(pixels, centers, half_w, half_h, backdrop)

    image = Image.fromarray(pixels, mode="RGB")

    if view == "wind":
        draw = ImageDraw.Draw(image)
        _draw_climate_vectors(draw, fields, fields.wind_u, fields.wind_v, projection, scale, offset_x, offset_y, pixel_scale, WIND_ARROW_COLOR)
    elif view == "oceanCurrents":
        draw = ImageDraw.Draw(image)
        _draw_climate_vectors(draw, fields, fields.current_u, fields.current_v, projection, scale, offset_x, offset_y, pixel_scale, CURRENT_ARROW_COLOR)
        _draw_swell_markers(draw, fields, projection, scale, offset_x, offset_y, pixel_scale)

    return _encode_image(image)


def render_png(world: World, projection: str, view: str, width: int, height: int) -> bytes:
    """Render `view` of `world` in `projection`, at `width`x`height` pixels, as PNG bytes.
    Mirrors what MapCanvas.tsx used to compute client-side from raw coordinate JSON -- this
    is now the only place that drawing logic lives."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    if view in CLIMATE_VIEWS:
        return _render_climate_view(world, projection, view, width, height)

    if not world.plates:
        return _encode_image(Image.fromarray(blank, mode="RGB"))

    grid = _render_grid_arrays(world, projection) if view in ("elevation", "plates") else None
    tectonics = {p.plate_id: _plate_tectonics(projection, p) for p in world.plates}

    detail_lines = []  # (projected_xy, elevation) per non-empty line, "platesDetail" only
    if view == "platesDetail":
        for plate in world.plates:
            for line in plate.lines:
                if len(line.theta) == 0:
                    continue
                detail_lines.append((_project_points(projection, line.world_xyz(plate.frame)), line.elevation))

    # Bounding box over every coordinate this view will draw (matches the old client-side
    # computation) so switching views never rescales or re-centers the map.
    bbox_chunks = []
    if grid is not None:
        bbox_chunks.append(grid[0])
    for info in tectonics.values():
        if len(info["boundary"]) > 0:
            bbox_chunks.append(info["boundary"])
        if info["pole"] is not None:
            bbox_chunks.append(info["pole"][None, :])
        # The rotation arc itself isn't included here: it's a fixed-pixel-radius decoration
        # drawn around the pole point after the transform is already fixed (see
        # _draw_rotation_arc), not additional world-space data to fit on screen.
    for xy, _ in detail_lines:
        bbox_chunks.append(xy)

    if not bbox_chunks:
        return _encode_image(Image.fromarray(blank, mode="RGB"))

    all_xy = np.concatenate(bbox_chunks, axis=0)
    min_x, min_y = all_xy.min(axis=0)
    max_x, max_y = all_xy.max(axis=0)
    data_w = max(max_x - min_x, 1e-9)
    data_h = max(max_y - min_y, 1e-9)
    scale = min((width - 2 * padding_px) / data_w, (height - 2 * padding_px) / data_h)
    offset_x = width / 2 - scale * (min_x + max_x) / 2
    offset_y = height / 2 + scale * (min_y + max_y) / 2

    pixels = blank.copy()

    if grid is not None:
        xy, elev, owner, half_w, half_h = grid
        centers = _to_pixels(scale, offset_x, offset_y, xy)
        hw_px = half_w * scale * CELL_OVERLAP_FACTOR
        hh_px = half_h * scale * CELL_OVERLAP_FACTOR
        colors = elevation_colors(elev) if view == "elevation" else plate_colors(owner)
        _fill_rects(pixels, centers, hw_px, hh_px, colors)

    if detail_lines:
        dot_radius = NODE_DOT_RADIUS_PX * pixel_scale
        for xy, elev in detail_lines:
            centers = _to_pixels(scale, offset_x, offset_y, xy)
            colors = elevation_colors(elev)
            _fill_rects(pixels, centers, dot_radius, dot_radius, colors)

    image = Image.fromarray(pixels, mode="RGB")

    if view in ("plates", "platesDetail"):
        draw = ImageDraw.Draw(image)
        for plate in world.plates:
            info = tectonics[plate.plate_id]
            color = tuple(int(c) for c in plate_colors(np.array([plate.plate_id]))[0])

            if len(info["boundary"]) > 0:
                boundary_px = _to_pixels(scale, offset_x, offset_y, info["boundary"])
                _stroke_robust_loop(draw, boundary_px, color, BOUNDARY_LINE_WIDTH_PX * pixel_scale)

            if view == "plates" and info["rotation_arc"] is not None:
                _draw_rotation_arc(draw, projection, info["rotation_arc"], scale, offset_x, offset_y, pixel_scale)

            if view == "plates" and info["pole"] is not None:
                px, py = _to_pixels(scale, offset_x, offset_y, info["pole"][None, :])[0]
                r = POLE_RADIUS_PX * pixel_scale
                # Filled with the plate's own color (matching its boundary/arc) rather than
                # a fixed color, with a white outline for contrast against the fill -- since
                # the pole can land anywhere on the map (see _plate_tectonics), color is what
                # ties a pole marker back to the plate it belongs to.
                draw.ellipse([px - r, py - r, px + r, py + r], fill=color, outline=(255, 255, 255), width=1)

    return _encode_image(image)


def _encode_image(image: Image.Image) -> bytes:
    buf = io.BytesIO()
    image.save(buf, format="PNG")
    return buf.getvalue()


def render_png_base64(world: World, projection: str, view: str, width: int, height: int) -> str:
    return base64.b64encode(render_png(world, projection, view, width, height)).decode("ascii")
