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

from . import geometry, mantle, plates, projections
from .world import World

VIEWS = ("elevation", "plates", "platesDetail")

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
ARROW_LINE_WIDTH_PX = 1.5
# Cells are drawn slightly larger than their measured half-extent so adjacent cells overlap
# a hair rather than risk a hairline gap from floating-point rounding. A ratio, not a pixel
# size, so it does not scale with requested resolution.
CELL_OVERLAP_FACTOR = 1.15
# Velocity arrows are always meant to be short, local indicators (ARROW_BASE_ANGULAR_LENGTH_RAD
# below caps them at ~0.15 rad on the sphere). If a plate's seed happens to sit near the map's
# antimeridian or a pole, projecting its (geodesically short) arrow endpoint can land it on the
# "other side" of the projection -- a tiny step in world space, but a huge jump in projected
# coordinates -- and drawing a straight line between the two would paint a stray line across
# most of the map. Any arrow this long on screen is that artifact, not a real velocity
# indicator, so it's skipped rather than drawn.
MAX_ARROW_FRACTION_OF_CANVAS = 0.15
ARROW_BASE_ANGULAR_LENGTH_RAD = 0.15

# Resolution of the render grid (see _render_grid_arrays), swept on a plate-independent
# global grid so the map's coverage never depends on how sparse any one plate's own line
# data looks once projected -- a fixed, display-oriented constant, deliberately *not* tied
# to plates.TARGET_LINE_SPACING_RAD (the simulation's physics resolution): the render grid
# only needs to look smooth once rasterized, which a resolution change in the physics has no
# bearing on.
GRID_SPACING_KM = 250.0
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
    """Pole marker, velocity arrow, and boundary outline for a plate -- everything the
    "Plates"/"Plates (details)" views draw besides the elevation-fill/node dots."""
    seed_xyz = plate.seed_world
    speed = float(np.linalg.norm(plate.omega))

    pole = None
    arrow = None
    if speed > 1e-15:
        pole_xyz = plate.omega / speed
        pole = _project_points(projection, pole_xyz[None, :])[0]

        direction = np.cross(plate.omega, seed_xyz) / speed
        intensity = np.clip(speed / mantle.MAX_PLATE_RATE, 0.3, 1.0)
        arrow_len = ARROW_BASE_ANGULAR_LENGTH_RAD * intensity
        end_xyz = np.cos(arrow_len) * seed_xyz + np.sin(arrow_len) * direction
        end_xyz = end_xyz / np.linalg.norm(end_xyz)
        start, end = _project_points(projection, np.stack([seed_xyz, end_xyz]))
        arrow = (start, end)

    outline_world = plate.outline_world()
    boundary = _project_points(projection, outline_world) if len(outline_world) > 0 else np.zeros((0, 2))

    return {"pole": pole, "velocity_arrow": arrow, "boundary": boundary}


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


def _draw_arrow(draw: ImageDraw.ImageDraw, x0, y0, x1, y1, color, width_px: float, head_len_px: float) -> None:
    line_width = max(int(round(width_px)), 1)
    draw.line([(x0, y0), (x1, y1)], fill=color, width=line_width)
    angle = np.arctan2(y1 - y0, x1 - x0)
    p1 = (x1 - head_len_px * np.cos(angle - np.pi / 6), y1 - head_len_px * np.sin(angle - np.pi / 6))
    p2 = (x1 - head_len_px * np.cos(angle + np.pi / 6), y1 - head_len_px * np.sin(angle + np.pi / 6))
    draw.polygon([(x1, y1), p1, p2], fill=color)


def render_png(world: World, projection: str, view: str, width: int, height: int) -> bytes:
    """Render `view` of `world` in `projection`, at `width`x`height` pixels, as PNG bytes.
    Mirrors what MapCanvas.tsx used to compute client-side from raw coordinate JSON -- this
    is now the only place that drawing logic lives."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

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
        if info["velocity_arrow"] is not None:
            bbox_chunks.append(np.stack(info["velocity_arrow"]))
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

            if view == "plates" and info["velocity_arrow"] is not None:
                start, end = info["velocity_arrow"]
                (sx, sy), (ex, ey) = _to_pixels(scale, offset_x, offset_y, np.stack([start, end]))
                max_arrow_px = min(width, height) * MAX_ARROW_FRACTION_OF_CANVAS
                if np.hypot(ex - sx, ey - sy) <= max_arrow_px:
                    _draw_arrow(
                        draw, sx, sy, ex, ey, (255, 255, 255),
                        ARROW_LINE_WIDTH_PX * pixel_scale, ARROWHEAD_LENGTH_PX * pixel_scale,
                    )

            if view == "plates" and info["pole"] is not None:
                px, py = _to_pixels(scale, offset_x, offset_y, info["pole"][None, :])[0]
                r = POLE_RADIUS_PX * pixel_scale
                draw.ellipse([px - r, py - r, px + r, py + r], fill=(255, 45, 85), outline=(255, 255, 255), width=1)

    return _encode_image(image)


def _encode_image(image: Image.Image) -> bytes:
    buf = io.BytesIO()
    image.save(buf, format="PNG")
    return buf.getvalue()


def render_png_base64(world: World, projection: str, view: str, width: int, height: int) -> str:
    return base64.b64encode(render_png(world, projection, view, width, height)).decode("ascii")
