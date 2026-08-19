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
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree

from . import biomes, climate, coastline, geometry, hydrology, mantle, plates, projections
from .world import World

# Climate views draw from climate.py's own fixed (H, W) grid, not the render grid below --
# see climate.py's module docstring for why. Handled by a separate code path
# (_render_climate_view) rather than threading a third data source through render_png's
# existing elevation/plates machinery. "biome" belongs here too, even though it isn't one of
# climate.py's own raw fields -- it's a pure classification (biomes.classify_biomes) derived
# entirely from two fields (temperature, precipitation) that path already has in hand.
CLIMATE_VIEWS = ("temperature", "wind", "oceanCurrents", "humidity", "precipitation", "biome")
VIEWS = ("elevation", "plates", "platesDetail") + CLIMATE_VIEWS

BACKGROUND_RGB = (11, 16, 32)  # #0b1020
# Muddier/less saturated than ocean blue (elevation_colors' own deep-water stop), matching
# plate-sim's own LAKE_COLOR_RGB choice for the same reason: a lake should read as visibly
# distinct from the open ocean, not just "more blue."
LAKE_COLOR_RGB = (58, 92, 122)
# Same color plate-sim's own frontend uses for its river overlay (#4dd8e6) -- a fixed color/
# width for every segment regardless of discharge, also matching plate-sim (only *which*
# segments get drawn varies with flow magnitude, via RIVER_FLOW_PERCENTILE below).
RIVER_COLOR_RGB = (77, 216, 230)
RIVER_LINE_WIDTH_PX = 1.1
# A second, independent cut on top of hydrology.py's own is_river classification
# (RIVER_FLOW_PERCENTILE): the main map views only draw a river segment if its own
# flow_accum also clears this absolute floor, so a large world's merely-top-decile trickles
# don't clutter every general-purpose view -- unlike the River Inspector (main.py's
# /world/rivers, RiverInspector.tsx), which deliberately keeps showing every is_river network
# regardless of flow_rate, since picking a small tributary out from the full list is exactly
# what that view is for.
RIVER_DRAW_MIN_FLOW = 10_000.0
# A pale icy blue-white -- deliberately distinct from both elevation_colors' own high-peak
# white/gray stops and LAKE_COLOR_RGB's darker muddy blue, so a glaciated node never reads
# as "just a tall mountain" or "just a lake" at a glance.
GLACIER_COLOR_RGB = (221, 240, 245)
# A hot, saturated red-orange -- distinct from elevation_colors' own high-peak browns/whites
# and from both LAKE_COLOR_RGB/GLACIER_COLOR_RGB's cool blues, so a volcanic node reads as
# "active/recent volcanism" at a glance rather than blending into ordinary high terrain.
VOLCANO_COLOR_RGB = (207, 63, 28)
# Coastline: drawn on views that have no other land/ocean cue at all (temperature/humidity/
# precipitation's color scales carry no land information on their own, unlike elevation's
# hypsometric coloring) -- see coastline.py. A single fixed color would vanish against parts
# of some of those scales (temperature's own white/black extremes in particular), so this is
# always drawn as a dark halo pass first, then a light line on top -- the same "halo" trick
# real maps use for a boundary that has to stay legible over an arbitrary backdrop color.
COASTLINE_COLOR_RGB = (235, 235, 235)
COASTLINE_HALO_RGB = (15, 15, 15)
COASTLINE_LINE_WIDTH_PX = 1.1
COASTLINE_HALO_WIDTH_PX = 2.6

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


# Spectral rainbow, degrees C: white at the cold end, violet -> indigo -> blue -> green ->
# yellow -> orange -> red through the middle, black at the hot end -- white/black are hard
# bounds (climate.py's land/air temperature can reach LAND_TEMP_MIN_C + LAND_TEMP_RANGE_C =
# 35, so the hottest cells do clamp to solid black, same as the coldest already clamp to
# solid white -- an intentional saturation effect, not a range that quietly clips real data
# without visual indication). green/yellow/red/black are pinned to specific requested
# values (0/10/20/30); violet/indigo/blue fill the -60..0 span evenly, orange sits at the
# midpoint of yellow and red since no exact value was requested for it.
_TEMPERATURE_STOP_C = np.array([-60.0, -45.0, -30.0, -15.0, 0.0, 10.0, 15.0, 20.0, 30.0])
_TEMPERATURE_STOP_RGB = np.array(
    [
        (255, 255, 255),  # white, -60
        (148, 0, 211),  # violet
        (75, 0, 130),  # indigo
        (0, 0, 255),  # blue
        (0, 200, 0),  # green
        (255, 255, 0),  # yellow
        (255, 140, 0),  # orange
        (255, 0, 0),  # red
        (0, 0, 0),  # black, +30
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
# Stops scaled 2x from an earlier 0-3000 range so the legend's max reads 6000, keeping the
# same relative color distribution rather than just extending the darkest stop's clamp range.
_PRECIPITATION_STOP_MM = np.array([0, 500, 1200, 2400, 4000, 6000], dtype=float)
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


def _project_offset(projection: str, world_pts: np.ndarray, reference_lon: np.ndarray) -> np.ndarray:
    """Like _project_points, but for a point known to be a *small* true angular step away
    from another point whose own longitude is `reference_lon` -- e.g. measuring a cell's
    extent by projecting a nearby offset. `xyz_to_latlon`'s atan2 jumps from +pi to -pi at
    the antimeridian, so two points a tiny 3D distance apart can still land on opposite sides
    of that jump (this becomes common, not a rare edge case, once the view can rotate
    arbitrarily -- the seam can fall anywhere). Longitude is unwrapped to the representative
    within pi of `reference_lon` before projecting, since a true small step should always be
    representable that way; ordinary data-point projection (_project_points) has no such
    reference and doesn't need this -- there, the seam simply *is* the map's left/right edge,
    correctly."""
    lat, lon = geometry.xyz_to_latlon(world_pts)
    lon = reference_lon + ((lon - reference_lon + np.pi) % (2 * np.pi) - np.pi)
    x, y = projections.project(projection, lat, lon)
    return np.stack([x, y], axis=-1)


def _rotate(world_pts: np.ndarray, view_rotation: np.ndarray) -> np.ndarray:
    """Applies the current view rotation to real world positions -- `view_rotation @ p` for
    each column-vector p, vectorized over row-shaped points as `world_pts @ view_rotation.T`
    (see geometry.py's to_world for the same row/column convention). Named `view_rotation`
    (not just `rotation`) throughout this module to keep it unambiguous next to
    `_plate_tectonics`'s unrelated `rotation_arc` (a plate's own Euler-pole tectonic
    rotation). This is a pure render-time transform: it only ever runs on positions
    immediately before they're projected to pixels, never on anything climate.py uses for
    actual physics (insolation, Coriolis, latitude-banded wind/currents all key off the
    true, un-rotated planetary frame -- see docs/simulation-model.md#rotating-the-view)."""
    return world_pts @ view_rotation.T


def _render_grid_arrays(
    world: World, projection: str, view_rotation: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None:
    """A uniform lat/lon grid covering the whole sphere (GRID_SPACING_RAD, independent of
    any plate's own line spacing), each cell assigned its nearest elevation node's elevation,
    owning plate, lake_depth, glacier_depth, and is_volcano -- see
    docs/simulation-model.md#render-image. Returns flat concatenated (projected_xy,
    elevation, plate_id, lake_depth, glacier_depth, is_volcano, cell_half_width,
    cell_half_height) arrays, or None for an empty world.

    Cell half-extents are measured per cell, not per row: at the identity rotation, a row of
    constant true latitude also has constant apparent latitude, so one measurement per row
    was enough (see git history for that older version). Once the view can be rotated
    arbitrarily that's no longer true -- a true-latitude row can span very different apparent
    positions depending on longitude, so cell size becomes a genuinely per-cell property.
    Measured via the cell's four *corners* (half a cell-step away on each diagonal), not its
    edge-midpoint neighbors: a rotation can turn what was an axis-aligned cell into a rotated
    or skewed quadrilateral, and a rotated square's axis-aligned bounding box is set by its
    corners, not by how far its edge midpoints sit from its center (which underestimates it,
    by a factor approaching sqrt(2) at a 45-degree rotation -- confirmed directly: this is
    what was producing widespread pinhole gaps under rotation before this was corner-based).
    Each corner is one extra projection call covering the whole row at once (four extra calls
    total per row, not per cell -- cells just no longer share a single measurement), via
    _project_offset rather than _project_points: once the view can rotate arbitrarily, the
    antimeridian seam (where xyz_to_latlon's atan2 jumps from +pi to -pi) can land anywhere
    within a row -- not a rare edge case, since the seam is a full line across the whole
    sphere that the grid's dense sweep is guaranteed to cross somewhere -- and a naive
    small-angle offset can still land on the far side of it, producing a wildly wrong extent
    (confirmed directly: this was the source of a handful of huge false cells smearing all
    the way across the render). _project_offset keeps the corner's longitude unwrapped
    relative to the cell's own center, which a true small step should always permit."""
    collected = plates.collect_all_points(world.plates)
    if collected is None:
        return None
    all_points, all_elevation, all_owner = collected
    all_lake_depth = plates.collect_all_lake_depth(world.plates)
    all_glacier_depth = plates.collect_all_glacier_depth(world.plates)
    all_is_volcano = plates.collect_all_is_volcano(world.plates)
    tree = cKDTree(all_points)

    xy_chunks, elev_chunks, owner_chunks, lake_chunks, glacier_chunks, volcano_chunks, hw_chunks, hh_chunks = [], [], [], [], [], [], [], []
    for phi, theta_candidates, world_pts in plates.iter_local_lattice(np.eye(3), spacing_rad=GRID_SPACING_RAD):
        _, idx = tree.query(world_pts)
        rotated = _rotate(world_pts, view_rotation)
        xy = _project_points(projection, rotated)
        xy_chunks.append(xy)
        elev_chunks.append(all_elevation[idx])
        owner_chunks.append(all_owner[idx])
        lake_chunks.append(all_lake_depth[idx])
        glacier_chunks.append(all_glacier_depth[idx])
        volcano_chunks.append(all_is_volcano[idx])

        _, center_lon = geometry.xyz_to_latlon(rotated)
        half_dtheta = GRID_SPACING_RAD / max(np.cos(phi), 1e-3) / 2
        half_dphi = GRID_SPACING_RAD / 2

        corner_steps = []
        for dphi_sign in (1.0, -1.0):
            for dtheta_sign in (1.0, -1.0):
                corner_local = geometry.local_xyz(
                    np.full_like(theta_candidates, phi + dphi_sign * half_dphi),
                    theta_candidates + dtheta_sign * half_dtheta,
                )
                corner_xy = _project_offset(projection, _rotate(corner_local, view_rotation), center_lon)
                corner_steps.append(corner_xy - xy)

        hw_chunks.append(np.max([np.abs(s[:, 0]) for s in corner_steps], axis=0))
        hh_chunks.append(np.max([np.abs(s[:, 1]) for s in corner_steps], axis=0))

    return (
        np.concatenate(xy_chunks, axis=0),
        np.concatenate(elev_chunks, axis=0),
        np.concatenate(owner_chunks, axis=0),
        np.concatenate(lake_chunks, axis=0),
        np.concatenate(glacier_chunks, axis=0),
        np.concatenate(volcano_chunks, axis=0),
        np.concatenate(hw_chunks, axis=0),
        np.concatenate(hh_chunks, axis=0),
    )


def _plate_tectonics(projection: str, plate, view_rotation: np.ndarray) -> dict:
    """Pole marker, rotation arc, and boundary outline for a plate -- everything the
    "Plates"/"Plates (details)" views draw besides the elevation-fill/node dots. `pole_xyz`
    inside `rotation_arc` is deliberately left in the *true* (un-rotated) frame -- it's
    combined with `omega` (also true-frame) for tangent/sweep-direction geometry in
    _draw_rotation_arc, which applies `view_rotation` itself at the point each of those
    intermediate points gets projected, same as everywhere else in this module."""
    speed = float(np.linalg.norm(plate.omega))

    pole = None
    rotation_arc = None
    if speed > 1e-15:
        # The true Euler pole -- real rotation axes are frequently nowhere near the plate
        # they belong to (this is physically normal, not a bug), so this can land anywhere
        # on the map. The pole marker is colored by plate (see render_png) specifically so
        # it still reads as "belonging to" the right plate even when it's drawn far away.
        pole_xyz = plate.omega / speed
        pole = _project_points(projection, _rotate(pole_xyz[None, :], view_rotation))[0]

        intensity = np.clip(speed / mantle.MAX_PLATE_RATE, 0.3, 1.0)
        rotation_arc = {
            "pole_xyz": pole_xyz,
            "omega": plate.omega,
            "extent_deg": ARC_MIN_EXTENT_DEG + intensity * (ARC_MAX_EXTENT_DEG - ARC_MIN_EXTENT_DEG),
        }

    outline_world = plate.outline_world()
    boundary = _project_points(projection, _rotate(outline_world, view_rotation)) if len(outline_world) > 0 else np.zeros((0, 2))

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
    view_rotation: np.ndarray,
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

    pole_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(pole_xyz[None, :], view_rotation)))[0]

    helper = np.array([1.0, 0.0, 0.0]) if abs(pole_xyz[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    tangent_dir = geometry.normalize(np.cross(pole_xyz, helper))
    near_point = np.cos(ARC_DIRECTION_SAMPLE_RAD) * pole_xyz + np.sin(ARC_DIRECTION_SAMPLE_RAD) * tangent_dir
    velocity = np.cross(omega, near_point)
    if np.linalg.norm(velocity) < 1e-15:
        return
    step_point = geometry.normalize(near_point + 1e-4 * velocity)

    near_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(near_point[None, :], view_rotation)))[0]
    step_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(step_point[None, :], view_rotation)))[0]

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
    view_rotation: np.ndarray,
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

    # Direction is computed from the *true* (un-rotated) local east/north basis, matching
    # how u/v were computed -- both the base and offset point are rotated together right
    # before projecting, same as everywhere else in this module.
    lon = np.arctan2(base_xyz[:, 1], base_xyz[:, 0])
    east = np.stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)], axis=-1)
    north = np.cross(base_xyz, east)
    direction = geometry.normalize(u_pts[:, None] * east + v_pts[:, None] * north)
    offset_xyz = geometry.normalize(base_xyz + 0.02 * direction)

    base_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(base_xyz, view_rotation)))
    offset_px = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(offset_xyz, view_rotation)))

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
    view_rotation: np.ndarray,
) -> None:
    if len(fields.swell_rows) == 0:
        return
    xyz = fields.world_xyz[fields.swell_rows, fields.swell_cols]
    centers = _to_pixels(scale, offset_x, offset_y, _project_points(projection, _rotate(xyz, view_rotation)))
    r = SWELL_MARKER_RADIUS_PX * pixel_scale
    width_px = max(int(round(pixel_scale)), 1)
    for px, py in centers:
        draw.ellipse([px - r, py - r, px + r, py + r], outline=(255, 255, 255), width=width_px)


# ---------------------------------------------------------------------------------------
# Legends
# ---------------------------------------------------------------------------------------

# All fixed pixel-space, scaled by pixel_scale like every other visual constant here (see
# module docstring / REFERENCE_WIDTH_PX). One shared panel layout (title, then an optional
# horizontal gradient bar with tick labels, then any number of symbol rows) covers every map
# view: a plain color scale for the heatmap views (elevation/temperature/humidity/
# precipitation), boundary/pole/rotation-arc symbols for the plate-tectonics views,
# arrow/backdrop/swell symbols for wind/oceanCurrents, and gradient + a boundary symbol for
# "Plates (details)" (its dots are colored by elevation, same scale as the elevation view).
LEGEND_MARGIN_PX = 16
LEGEND_PADDING_PX = 10
LEGEND_PANEL_WIDTH_PX = 172
LEGEND_TITLE_ROW_PX = 20
LEGEND_ROW_PX = 18
LEGEND_BAR_LENGTH_PX = 148
LEGEND_BAR_THICKNESS_PX = 14
LEGEND_TICK_ROW_PX = 16
LEGEND_SWATCH_PX = 14
LEGEND_TITLE_FONT_PX = 13
LEGEND_LABEL_FONT_PX = 11
LEGEND_BG_RGB = (16, 20, 34)
LEGEND_BORDER_RGB = (75, 80, 96)
LEGEND_TEXT_RGB = (222, 226, 235)
LEGEND_SYMBOL_RGB = (200, 205, 215)  # neutral color for symbols not tied to any one plate


def _legend_font(size_px: float) -> ImageFont.FreeTypeFont:
    return ImageFont.load_default(size=max(int(round(size_px)), 6))


def _draw_legend(
    image: Image.Image,
    draw: ImageDraw.ImageDraw,
    height: int,
    pixel_scale: float,
    title: str,
    gradient: tuple[object, float, float, list[tuple[float, str]]] | None = None,
    symbol_entries: list[tuple[object, str]] = (),
) -> None:
    """Draws a legend panel anchored at the image's bottom-left corner. `gradient`, if given,
    is `(colors_fn, value_min, value_max, [(tick_value, tick_label), ...])` -- colors_fn is
    one of this module's own vectorized value -> RGB functions (elevation_colors etc.),
    called here against a value ramp to build the bar image directly, so the legend is always
    in sync with whatever coloring the map itself used. `symbol_entries` is
    `[(draw_swatch(cx, cy, size), label), ...]`, stacked below the gradient (if any)."""
    pad = int(round(LEGEND_PADDING_PX * pixel_scale))
    title_h = int(round(LEGEND_TITLE_ROW_PX * pixel_scale))
    row_h = int(round(LEGEND_ROW_PX * pixel_scale))
    bar_w = int(round(LEGEND_BAR_LENGTH_PX * pixel_scale))
    bar_h = int(round(LEGEND_BAR_THICKNESS_PX * pixel_scale))
    tick_row_h = int(round(LEGEND_TICK_ROW_PX * pixel_scale))
    panel_w = max(int(round(LEGEND_PANEL_WIDTH_PX * pixel_scale)), bar_w + 2 * pad)

    content_h = (bar_h + tick_row_h if gradient is not None else 0) + row_h * len(symbol_entries)
    panel_h = title_h + content_h + 2 * pad

    margin = int(round(LEGEND_MARGIN_PX * pixel_scale))
    x0 = margin
    y0 = height - margin - panel_h
    draw.rectangle([x0, y0, x0 + panel_w, y0 + panel_h], fill=LEGEND_BG_RGB, outline=LEGEND_BORDER_RGB, width=1)

    title_font = _legend_font(LEGEND_TITLE_FONT_PX * pixel_scale)
    label_font = _legend_font(LEGEND_LABEL_FONT_PX * pixel_scale)
    inner_x = x0 + pad
    cy = y0 + pad
    draw.text((inner_x, cy), title, font=title_font, fill=LEGEND_TEXT_RGB)
    cy += title_h

    if gradient is not None:
        colors_fn, value_min, value_max, ticks = gradient
        values = np.linspace(value_min, value_max, bar_w)
        bar_row = colors_fn(values)
        bar_array = np.broadcast_to(bar_row[None, :, :], (bar_h, bar_w, 3))
        image.paste(Image.fromarray(bar_array, mode="RGB"), (inner_x, cy))
        tick_y = cy + bar_h
        span = value_max - value_min
        for value, label in ticks:
            frac = 0.0 if span <= 0 else min(max((value - value_min) / span, 0.0), 1.0)
            tick_x = inner_x + frac * (bar_w - 1)
            draw.line([(tick_x, cy), (tick_x, tick_y)], fill=(0, 0, 0), width=1)
            draw.text((tick_x, tick_y + 1), label, font=label_font, fill=LEGEND_TEXT_RGB, anchor="ma")
        cy += bar_h + tick_row_h

    swatch = int(round(LEGEND_SWATCH_PX * pixel_scale))
    label_x = inner_x + swatch + int(round(6 * pixel_scale))
    for draw_swatch, label in symbol_entries:
        row_mid = cy + row_h / 2
        draw_swatch(inner_x, row_mid, swatch)
        draw.text((label_x, row_mid), label, font=label_font, fill=LEGEND_TEXT_RGB, anchor="lm")
        cy += row_h


def _swatch_line(draw: ImageDraw.ImageDraw, color: tuple[int, int, int]):
    return lambda cx, cy, size: draw.line([(cx, cy), (cx + size, cy)], fill=color, width=max(int(round(size / 6)), 1))


def _swatch_square(draw: ImageDraw.ImageDraw, color) -> object:
    color = tuple(int(c) for c in color)
    return lambda cx, cy, size: draw.rectangle([cx, cy - size / 2, cx + size, cy + size / 2], fill=color)


def _swatch_circle(draw: ImageDraw.ImageDraw, fill_color: tuple[int, int, int], outline_color: tuple[int, int, int]):
    def paint(cx: float, cy: float, size: float) -> None:
        r = size / 2
        center = cx + r
        draw.ellipse([center - r, cy - r, center + r, cy + r], fill=fill_color, outline=outline_color, width=1)

    return paint


def _swatch_ring(draw: ImageDraw.ImageDraw, color: tuple[int, int, int]):
    def paint(cx: float, cy: float, size: float) -> None:
        r = size / 2
        center = cx + r
        draw.ellipse([center - r, cy - r, center + r, cy + r], outline=color, width=max(int(round(size / 8)), 1))

    return paint


def _swatch_arrow(draw: ImageDraw.ImageDraw, color: tuple[int, int, int]):
    def paint(cx: float, cy: float, size: float) -> None:
        draw.line([(cx, cy), (cx + size, cy)], fill=color, width=max(int(round(size / 7)), 1))
        _draw_arrow_head(draw, np.array([cx + size, cy]), np.array([1.0, 0.0]), color, size * 0.5)

    return paint


def _swatch_arc(draw: ImageDraw.ImageDraw, color: tuple[int, int, int]):
    def paint(cx: float, cy: float, size: float) -> None:
        r = size / 2
        center_x = cx + r
        draw.arc([center_x - r, cy - r, center_x + r, cy + r], 30, 300, fill=color, width=max(int(round(size / 8)), 1))
        head_rad = np.radians(300.0)
        head_px = np.array([center_x, cy]) + r * np.array([np.cos(head_rad), np.sin(head_rad)])
        tangent = np.array([-np.sin(head_rad), np.cos(head_rad)])
        _draw_arrow_head(draw, head_px, tangent, color, size * 0.35)

    return paint


def _elevation_legend_gradient() -> tuple[object, float, float, list[tuple[float, str]]]:
    lo, hi = float(_ELEVATION_STOP_E[0]), float(_ELEVATION_STOP_E[-1])
    ticks = [(-8000.0, "-8k"), (-4000.0, "-4k"), (0.0, "0"), (4000.0, "4k"), (8000.0, "8k")]
    return elevation_colors, lo, hi, ticks


def _render_climate_view(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """Renders one of CLIMATE_VIEWS from climate.py's own fixed grid -- a separate path from
    the plate-tectonics views below since the data source (a real (H, W) array, always
    covering the whole sphere) is structurally different from the render grid's ragged
    lattice, so there's little to share beyond the pixel-space primitives."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    # Reuses whatever erosion.py already computed this step (see World.climate_cache)
    # instead of triggering a second ~50ms recomputation, rather than calling
    # climate.compute_climate directly.
    fields = climate.compute_climate_cached(world)
    grid_h, grid_w = fields.elevation_m.shape
    flat_xyz = fields.world_xyz.reshape(-1, 3)

    # The bounding box is always the *whole* sphere's projected extent (climate's grid
    # always covers it), which is identical regardless of rotation -- rotation only permutes
    # which physical points land at which lat/lon, never changes the set of lat/lon values a
    # full sphere covers. So it's measured from the un-rotated grid (matching today's
    # behavior exactly at the identity rotation), while `flat_xy` below -- what's actually
    # drawn -- uses the rotated positions.
    bbox_xy = _project_points(projection, flat_xyz)
    flat_xy = _project_points(projection, _rotate(flat_xyz, view_rotation))

    min_x, min_y = bbox_xy.min(axis=0)
    max_x, max_y = bbox_xy.max(axis=0)
    data_w = max(max_x - min_x, 1e-9)
    data_h = max(max_y - min_y, 1e-9)
    scale = min((width - 2 * padding_px) / data_w, (height - 2 * padding_px) / data_h)
    offset_x = width / 2 - scale * (min_x + max_x) / 2
    offset_y = height / 2 + scale * (min_y + max_y) / 2

    pixels = blank.copy()
    centers = _to_pixels(scale, offset_x, offset_y, flat_xy)

    # Cell half-extents, per cell (see _render_grid_arrays's docstring for the full
    # reasoning, which applies identically here): measured via the cell's four *corners*
    # (half a cell-step away on each diagonal), not edge-midpoint neighbors, since a rotation
    # can skew an axis-aligned cell into a shape whose axis-aligned bounding box is set by its
    # corners; and via _project_offset (not _project_points), since a naive small-angle
    # corner offset can still land on the far side of the antimeridian seam once the view can
    # rotate arbitrarily -- the seam can fall anywhere, and climate's grid densely covers the
    # whole sphere, so some cell is always near it somewhere.
    xy_grid = flat_xy.reshape(grid_h, grid_w, 2)
    _, center_lon_grid = geometry.xyz_to_latlon(_rotate(flat_xyz, view_rotation))
    center_lon_grid = center_lon_grid.reshape(grid_h, grid_w)

    half_dlat_deg = 180.0 / grid_h / 2
    half_dlon_deg = 360.0 / grid_w / 2
    lat_grid_deg = np.repeat(fields.lat_deg[:, None], grid_w, axis=1)
    lon_grid_deg = np.repeat(fields.lon_deg[None, :], grid_h, axis=0)

    def _corner_xy(dlat: float, dlon: float) -> np.ndarray:
        xyz = geometry.latlon_to_xyz(np.radians(lat_grid_deg + dlat), np.radians(lon_grid_deg + dlon))
        rotated = _rotate(xyz.reshape(-1, 3), view_rotation)
        return _project_offset(projection, rotated, center_lon_grid.reshape(-1)).reshape(grid_h, grid_w, 2)

    corner_steps = [
        _corner_xy(dlat_sign * half_dlat_deg, dlon_sign * half_dlon_deg) - xy_grid
        for dlat_sign in (1.0, -1.0)
        for dlon_sign in (1.0, -1.0)
    ]
    half_w = np.max([np.abs(s[..., 0]) for s in corner_steps], axis=0).reshape(-1) * scale * CELL_OVERLAP_FACTOR
    half_h = np.max([np.abs(s[..., 1]) for s in corner_steps], axis=0).reshape(-1) * scale * CELL_OVERLAP_FACTOR

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
    elif view == "biome":
        # Same land/ocean temperature split the "temperature" view above uses -- ocean cells
        # always classify as Ocean regardless, so only the land branch actually matters here.
        display_temp = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
        biome_ids = biomes.classify_biomes(display_temp.reshape(-1), fields.precipitation_mm.reshape(-1), fields.is_ocean.reshape(-1))
        colors = biomes.BIOME_COLORS[biome_ids]
        _fill_rects(pixels, centers, half_w, half_h, colors)
    elif view in ("wind", "oceanCurrents"):
        backdrop = np.where(fields.is_ocean.reshape(-1)[:, None], CLIMATE_OCEAN_BACKDROP_RGB, CLIMATE_LAND_BACKDROP_RGB)
        _fill_rects(pixels, centers, half_w, half_h, backdrop)

    image = Image.fromarray(pixels, mode="RGB")
    draw = ImageDraw.Draw(image)

    if view in ("temperature", "humidity", "precipitation"):
        _draw_coastline(draw, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    if view == "wind":
        _draw_climate_vectors(draw, fields, fields.wind_u, fields.wind_v, projection, scale, offset_x, offset_y, pixel_scale, WIND_ARROW_COLOR, view_rotation)
    elif view == "oceanCurrents":
        _draw_climate_vectors(draw, fields, fields.current_u, fields.current_v, projection, scale, offset_x, offset_y, pixel_scale, CURRENT_ARROW_COLOR, view_rotation)
        _draw_swell_markers(draw, fields, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    coastline_entries = [(_swatch_line(draw, COASTLINE_COLOR_RGB), "Coastline")]
    if view == "temperature":
        ticks = [(-60.0, "-60°"), (0.0, "0°"), (10.0, "10°"), (20.0, "20°"), (30.0, "30°")]
        gradient = (temperature_colors, float(_TEMPERATURE_STOP_C[0]), float(_TEMPERATURE_STOP_C[-1]), ticks)
        _draw_legend(image, draw, height, pixel_scale, "Temperature (°C)", gradient=gradient, symbol_entries=coastline_entries)
    elif view == "humidity":
        lo, hi = float(_HUMIDITY_STOP_V[0]), float(_HUMIDITY_STOP_V[-1])
        gradient = (humidity_colors, lo, hi, [(lo, "Dry"), (hi, "Humid")])
        _draw_legend(image, draw, height, pixel_scale, "Humidity", gradient=gradient, symbol_entries=coastline_entries)
    elif view == "precipitation":
        lo, hi = float(_PRECIPITATION_STOP_MM[0]), float(_PRECIPITATION_STOP_MM[-1])
        ticks = [(lo, "0"), (hi / 2, f"{hi / 2:.0f}"), (hi, f"{hi:.0f}")]
        gradient = (precipitation_colors, lo, hi, ticks)
        _draw_legend(image, draw, height, pixel_scale, "Precipitation (mm/yr)", gradient=gradient, symbol_entries=coastline_entries)
    elif view == "wind":
        entries = [
            (_swatch_arrow(draw, WIND_ARROW_COLOR), "Speed (arrow length)"),
            (_swatch_square(draw, CLIMATE_OCEAN_BACKDROP_RGB), "Ocean"),
            (_swatch_square(draw, CLIMATE_LAND_BACKDROP_RGB), "Land"),
        ]
        _draw_legend(image, draw, height, pixel_scale, "Wind", symbol_entries=entries)
    elif view == "oceanCurrents":
        entries = [
            (_swatch_arrow(draw, CURRENT_ARROW_COLOR), "Speed (arrow length)"),
            (_swatch_square(draw, CLIMATE_OCEAN_BACKDROP_RGB), "Ocean"),
            (_swatch_square(draw, CLIMATE_LAND_BACKDROP_RGB), "Land"),
            (_swatch_ring(draw, (255, 255, 255)), "Ocean swell"),
        ]
        _draw_legend(image, draw, height, pixel_scale, "Ocean currents", symbol_entries=entries)
    elif view == "biome":
        entries = [(_swatch_square(draw, biomes.BIOME_COLORS[i]), name) for i, name in enumerate(biomes.BIOME_NAMES)]
        _draw_legend(image, draw, height, pixel_scale, "Biome", symbol_entries=entries)

    return _encode_image(image)


def _draw_rivers(
    draw: ImageDraw.ImageDraw, world: World, projection: str, scale: float, offset_x: float, offset_y: float, pixel_scale: float, view_rotation: np.ndarray
) -> None:
    """Draws each river node's edge to its own downstream flow target as a short line
    segment (see hydrology.py's is_river/flow_target) -- reads world.hydrology_cache
    directly, populated by erosion.py every step (None before the world has ever been
    stepped, in which case this draws nothing). Each segment is a real, short 3D hop between
    two adjacent-in-the-flow-graph nodes, so _project_offset (not two independent
    _project_points calls) keeps it from being wrongly split across the antimeridian seam --
    same technique _render_grid_arrays' own corner measurements already rely on. Also cut by
    RIVER_DRAW_MIN_FLOW -- see that constant's own comment for why this view is stricter than
    the River Inspector's own /world/rivers listing."""
    hydro = world.hydrology_cache
    if hydro is None:
        return
    river_idx = np.nonzero(hydro.is_river & (hydro.flow_target >= 0) & (hydro.flow_accum >= RIVER_DRAW_MIN_FLOW))[0]
    if len(river_idx) == 0:
        return
    target_idx = hydro.flow_target[river_idx]

    from_points = _rotate(hydro.points[river_idx], view_rotation)
    to_points = _rotate(hydro.points[target_idx], view_rotation)
    _, from_lon = geometry.xyz_to_latlon(from_points)
    from_xy = _project_points(projection, from_points)
    to_xy = _project_offset(projection, to_points, from_lon)

    from_px = _to_pixels(scale, offset_x, offset_y, from_xy)
    to_px = _to_pixels(scale, offset_x, offset_y, to_xy)
    width_px = max(int(round(RIVER_LINE_WIDTH_PX * pixel_scale)), 1)
    for (x1, y1), (x2, y2) in zip(from_px, to_px):
        draw.line([(x1, y1), (x2, y2)], fill=RIVER_COLOR_RGB, width=width_px)


def _draw_coastline(
    draw: ImageDraw.ImageDraw, world: World, projection: str, scale: float, offset_x: float, offset_y: float, pixel_scale: float, view_rotation: np.ndarray
) -> None:
    """Draws the land/lake-vs-ocean boundary (see coastline.compute_coastline_segments) as a
    dark-halo-plus-light-line stroke, so it stays legible over any backdrop color. Each
    segment is a real, short step between two adjacent grid-cell corners, so _project_offset
    (not two independent _project_points calls) keeps it from being wrongly split across the
    antimeridian seam -- same technique _draw_rivers already uses for the same reason."""
    point_a, point_b = coastline.compute_coastline_segments(world)
    if len(point_a) == 0:
        return

    from_points = _rotate(point_a, view_rotation)
    to_points = _rotate(point_b, view_rotation)
    _, from_lon = geometry.xyz_to_latlon(from_points)
    from_xy = _project_points(projection, from_points)
    to_xy = _project_offset(projection, to_points, from_lon)

    from_px = _to_pixels(scale, offset_x, offset_y, from_xy)
    to_px = _to_pixels(scale, offset_x, offset_y, to_xy)
    halo_width_px = max(int(round(COASTLINE_HALO_WIDTH_PX * pixel_scale)), 1)
    line_width_px = max(int(round(COASTLINE_LINE_WIDTH_PX * pixel_scale)), 1)
    segments = list(zip(from_px, to_px))
    for (x1, y1), (x2, y2) in segments:
        draw.line([(x1, y1), (x2, y2)], fill=COASTLINE_HALO_RGB, width=halo_width_px)
    for (x1, y1), (x2, y2) in segments:
        draw.line([(x1, y1), (x2, y2)], fill=COASTLINE_COLOR_RGB, width=line_width_px)


def render_png(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray | None = None) -> bytes:
    """Render `view` of `world` in `projection`, at `width`x`height` pixels, as PNG bytes.
    Mirrors what MapCanvas.tsx used to compute client-side from raw coordinate JSON -- this
    is now the only place that drawing logic lives. `view_rotation` (default identity, i.e.
    today's behavior exactly, center at lat=0/lon=0) is a pure render-time transform -- see
    _rotate's docstring and docs/simulation-model.md#rotating-the-view."""
    if view_rotation is None:
        view_rotation = np.eye(3)
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    if view in CLIMATE_VIEWS:
        return _render_climate_view(world, projection, view, width, height, view_rotation)

    if not world.plates:
        return _encode_image(Image.fromarray(blank, mode="RGB"))

    grid = _render_grid_arrays(world, projection, view_rotation) if view in ("elevation", "plates") else None
    tectonics = {p.plate_id: _plate_tectonics(projection, p, view_rotation) for p in world.plates}

    detail_lines = []  # (projected_xy, elevation) per non-empty line, "platesDetail" only
    if view == "platesDetail":
        for plate in world.plates:
            for line in plate.lines:
                if len(line.theta) == 0:
                    continue
                xy = _project_points(projection, _rotate(line.world_xyz(plate.frame), view_rotation))
                detail_lines.append((xy, line.elevation))

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
        xy, elev, owner, lake_depth, glacier_depth, is_volcano, half_w, half_h = grid
        centers = _to_pixels(scale, offset_x, offset_y, xy)
        hw_px = half_w * scale * CELL_OVERLAP_FACTOR
        hh_px = half_h * scale * CELL_OVERLAP_FACTOR
        colors = elevation_colors(elev) if view == "elevation" else plate_colors(owner)
        # Baked directly into the raster rather than a separate overlay/toggle -- same
        # reasoning plate-sim's own docs give for its own lake baking: always visible, no
        # separate overlay needed, and a lake (or a glacier) is meaningful on every view that
        # shows terrain at all (matches how ocean itself isn't specially toggle-able either).
        # Volcano drawn after lake (a volcanic vent is dry land in practice, but this keeps
        # the same precedence as everything else here), glacier drawn last so ice wins on the
        # rare cell where more than one would apply (a lake that just froze this same step
        # still shows lake_depth from the previous render's snapshot for one extra frame in
        # the worst case, and a volcano cold enough to glaciate should read as ice-covered,
        # not lava-red) -- ice is the more physically current/visually dominant state there.
        is_lake = lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
        if np.any(is_lake):
            colors = np.where(is_lake[:, None], np.array(LAKE_COLOR_RGB, dtype=np.uint8), colors)
        if np.any(is_volcano):
            colors = np.where(is_volcano[:, None], np.array(VOLCANO_COLOR_RGB, dtype=np.uint8), colors)
        is_glacier = glacier_depth > hydrology.GLACIER_VISIBLE_DEPTH_M
        if np.any(is_glacier):
            colors = np.where(is_glacier[:, None], np.array(GLACIER_COLOR_RGB, dtype=np.uint8), colors)
        _fill_rects(pixels, centers, hw_px, hh_px, colors)

    if detail_lines:
        dot_radius = NODE_DOT_RADIUS_PX * pixel_scale
        for xy, elev in detail_lines:
            centers = _to_pixels(scale, offset_x, offset_y, xy)
            colors = elevation_colors(elev)
            _fill_rects(pixels, centers, dot_radius, dot_radius, colors)

    image = Image.fromarray(pixels, mode="RGB")
    draw = ImageDraw.Draw(image)

    _draw_rivers(draw, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    if view in ("plates", "platesDetail"):
        for plate in world.plates:
            info = tectonics[plate.plate_id]
            color = tuple(int(c) for c in plate_colors(np.array([plate.plate_id]))[0])

            if len(info["boundary"]) > 0:
                boundary_px = _to_pixels(scale, offset_x, offset_y, info["boundary"])
                _stroke_robust_loop(draw, boundary_px, color, BOUNDARY_LINE_WIDTH_PX * pixel_scale)

            if view == "plates" and info["rotation_arc"] is not None:
                _draw_rotation_arc(draw, projection, info["rotation_arc"], scale, offset_x, offset_y, pixel_scale, view_rotation)

            if view == "plates" and info["pole"] is not None:
                px, py = _to_pixels(scale, offset_x, offset_y, info["pole"][None, :])[0]
                r = POLE_RADIUS_PX * pixel_scale
                # Filled with the plate's own color (matching its boundary/arc) rather than
                # a fixed color, with a white outline for contrast against the fill -- since
                # the pole can land anywhere on the map (see _plate_tectonics), color is what
                # ties a pole marker back to the plate it belongs to.
                draw.ellipse([px - r, py - r, px + r, py + r], fill=color, outline=(255, 255, 255), width=1)

    if view == "elevation":
        entries = [
            (_swatch_line(draw, RIVER_COLOR_RGB), "River"),
            (_swatch_square(draw, LAKE_COLOR_RGB), "Lake"),
            (_swatch_square(draw, VOLCANO_COLOR_RGB), "Volcano"),
            (_swatch_square(draw, GLACIER_COLOR_RGB), "Glacier"),
        ]
        _draw_legend(image, draw, height, pixel_scale, "Elevation (m)", gradient=_elevation_legend_gradient(), symbol_entries=entries)
    elif view == "plates":
        entries = [
            (_swatch_line(draw, LEGEND_SYMBOL_RGB), "Plate boundary"),
            (_swatch_circle(draw, LEGEND_SYMBOL_RGB, (255, 255, 255)), "Euler pole"),
            (_swatch_arc(draw, LEGEND_SYMBOL_RGB), "Rotation (rate & direction)"),
        ]
        _draw_legend(image, draw, height, pixel_scale, "Plates", symbol_entries=entries)
    elif view == "platesDetail":
        entries = [(_swatch_line(draw, LEGEND_SYMBOL_RGB), "Plate boundary")]
        _draw_legend(
            image, draw, height, pixel_scale, "Elevation (m)",
            gradient=_elevation_legend_gradient(), symbol_entries=entries,
        )

    return _encode_image(image)


def _encode_image(image: Image.Image) -> bytes:
    buf = io.BytesIO()
    image.save(buf, format="PNG")
    return buf.getvalue()


def render_png_base64(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray | None = None) -> str:
    return base64.b64encode(render_png(world, projection, view, width, height, view_rotation)).decode("ascii")
