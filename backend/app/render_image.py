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
from PIL import Image, ImageDraw, ImageFilter
from scipy.spatial import cKDTree

from . import biomes, climate, coastline, erosion, geology, geometry, hydrology, mantle, plates, projections, volcanism
from .world import World, step_world

# Climate views draw from climate.py's own fixed (H, W) grid, not the render grid below --
# see climate.py's module docstring for why. Handled by a separate code path
# (_render_climate_view) rather than threading a third data source through render_png's
# existing elevation/plates machinery. "biome" belongs here too, even though it isn't one of
# climate.py's own raw fields -- it's a pure classification (biomes.classify_biomes) derived
# entirely from two fields (temperature, precipitation) that path already has in hand.
CLIMATE_VIEWS = ("temperature", "wind", "oceanCurrents", "humidity", "precipitation", "biome")
# "combined" isn't a CLIMATE_VIEWS member (see _render_combined_view): it draws from the same
# fine per-node elevation data the elevation/plates views use, not climate.py's own grid, so
# it belongs with them structurally even though its coloring leans on biome classification.
# Resources/Soil Quality (see geology.py/volcanism.py) are node-cloud-derived like elevation/
# plates, not climate-grid-derived -- see _resource_fields -- but get their own tuple/dispatch
# branch (_render_resource_view) rather than folding into VIEWS' first group directly, since
# they share a render path with each other (one shared fine-grid resample) but not with
# elevation/plates' own render-grid machinery.
RESOURCE_VIEWS = ("resources", "soilQuality")
# Ocean/Atmospheric Fluid Dynamics's own map views -- see ocean_cfd.py/atmosphere_cfd.py and
# docs/simulation-model.md#ocean-atmospheric-fluid-dynamics. World.ocean_cfd_state/
# atmosphere_cfd_state are always populated now (see the former's own docstring), so these are
# always renderable, same as every other view here.
#
# TODO: oceanCfdVelocity/oceanCfdTemperature and atmosphereCfdVelocity/atmosphereCfdTemperature
# are now redundant with CLIMATE_VIEWS' own "wind"/"oceanCurrents"/"temperature" -- climate.py's
# compute_climate sources wind_u/wind_v/current_u/current_v straight from these same CFD
# states (see its own module docstring), just resampled onto climate_density's grid instead of
# fluid_density's. Worth consolidating (drop the duplicates, or make the plain views always use
# native fluid_density resolution) once it's clear which framing the UI wants long-term.
# oceanCfdSediment/oceanCfdDeposition stay unique -- nothing else produces sediment data.
OCEAN_CFD_VIEWS = ("oceanCfdVelocity", "oceanCfdTemperature", "oceanCfdSediment", "oceanCfdDeposition")
ATMOSPHERE_CFD_VIEWS = ("atmosphereCfdVelocity", "atmosphereCfdTemperature", "atmosphereCfdHumidity")
FLUID_VIEWS = OCEAN_CFD_VIEWS + ATMOSPHERE_CFD_VIEWS
VIEWS = ("elevation", "plates", "platesDetail", "combined") + CLIMATE_VIEWS + RESOURCE_VIEWS + FLUID_VIEWS

BACKGROUND_RGB = (11, 16, 32)  # #0b1020
# Muddier/less saturated than ocean blue (elevation_colors' own deep-water stop) -- a lake
# should read as visibly distinct from the open ocean, not just "more blue."
LAKE_COLOR_RGB = (58, 92, 122)
# Fixed river overlay color (#4dd8e6); line width is not fixed -- see _draw_rivers, which
# steps it directly off the segment's own flow_accum: exactly at river_draw_min_flow(world)
# (the minimum flow a segment can have and still be drawn at all) a river is a single pixel
# wide, and it only steps up to 2x/3x that width once flow_accum itself reaches 2x/3x the
# floor. Since flow_accum along an unbranched stretch of channel is constant (it only grows
# where a tributary's own flow merges in), this keeps a river's drawn width flat until a real
# confluence, and wider from then on toward the mouth -- never gradually thickening on its own
# -- only *which* segments get drawn at all varies with flow magnitude directly, via
# RIVER_FLOW_PERCENTILE below.
RIVER_COLOR_RGB = (77, 216, 230)
RIVER_LINE_WIDTH_PX = 1.1
# A light Gaussian blur of just the river-line mask before compositing the fixed river color
# in by the blurred mask's own value as a per-pixel alpha -- the same cheap-AA idea
# COMBINED_BLUR_RADIUS_PX already uses for cell edges (see that constant's own comment), just
# applied to a line mask instead of the filled-cell raster. Blurring only the (single-channel)
# mask and then blending a flat color by it, rather than Gaussian-blurring the drawn RGBA line
# directly, avoids the dark fringing a naive blur would produce against a transparent
# backdrop, and leaves whatever was already drawn underneath (coastline, plate boundaries)
# untouched -- unlike blurring the whole image, which would re-soften those too.
RIVER_BLUR_RADIUS_PX = 0.6
# Widest a river is ever drawn (at 3x river_draw_min_flow(world) or beyond) -- see
# RIVER_LINE_WIDTH_MAX_TIERS.
RIVER_LINE_WIDTH_MAX_TIERS = 3
# A second, independent cut on top of hydrology.py's own is_river classification
# (RIVER_FLOW_PERCENTILE): the main map views only draw a river segment if its own
# flow_accum also clears this absolute floor, so a large world's merely-top-decile trickles
# don't clutter every general-purpose view -- unlike the River Inspector (main.py's
# /world/rivers, RiverInspector.tsx), which deliberately keeps showing every is_river network
# regardless of flow_rate, since picking a small tributary out from the full list is exactly
# what that view is for. This is the floor at the *lowest* generation resolution (node_density
# == climate_density == 1.0) -- see river_draw_min_flow below for how a world generated at a
# higher elevation-point density and/or climate & biome resolution scales it up from here, not
# a value read directly by _draw_rivers.
RIVER_DRAW_MIN_FLOW = 5_000.0


def river_draw_min_flow(world: World) -> float:
    """The effective RIVER_DRAW_MIN_FLOW floor for `world`'s own generation resolution --
    higher at a higher elevation-point density (World.node_density) and/or climate & biome
    resolution (World.climate_density), since both mean flow_accum is being swept from a
    finer-grained network with more individual nodes/cells each contributing to it, so the
    same absolute floor would let proportionally more low-magnitude clutter through than it
    does at the reference (1.0, 1.0) resolution this constant was tuned at. Scales each
    density the same way its own module already relates it to node/cell *count* -- linearly
    for node_density (plates.line_spacing_rad's own docstring: node_density itself is the
    count-scaling factor relative to the reference), squared for climate_density
    (climate.grid_dimensions scales *both* grid dimensions by climate_density, so cell count
    scales by climate_density**2) -- rather than inventing a third, unrelated scaling rule
    here."""
    return RIVER_DRAW_MIN_FLOW * world.node_density * (world.climate_density ** 2)


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

# Resources view (see geology.py/volcanism.py): a muted, deliberately low-saturation
# land/ocean backdrop -- distinct from every other view's own colors -- so the three
# resource overlays (blended in on top by richness fraction, see _render_resources_view)
# read as the whole point of this view rather than competing with a busy backdrop.
RESOURCE_LAND_BACKDROP_RGB = np.array([82, 76, 68], dtype=float)
RESOURCE_OCEAN_BACKDROP_RGB = np.array([22, 32, 48], dtype=float)
COAL_COLOR_RGB = np.array([28, 24, 22], dtype=float)  # near-black charcoal
OIL_GAS_COLOR_RGB = np.array([92, 56, 14], dtype=float)  # dark amber/brown-black crude
MINERAL_COLOR_RGB = np.array([175, 62, 205], dtype=float)  # vivid purple -- ore veins

# Soil Quality view: barren bare rock -> rich, near-black fertile soil (real chernozem
# "black earth" is the visual reference). fertility (see _render_soil_view) is already
# [0, 1], so these stops are keyed on that directly rather than a physical unit.
_SOIL_STOP_V = np.array([0.0, 0.15, 0.35, 0.6, 1.0], dtype=float)
_SOIL_STOP_RGB = np.array(
    [
        (168, 150, 120),  # barren -- pale rocky tan
        (150, 120, 80),  # poor -- thin dry soil
        (120, 90, 55),  # moderate
        (80, 58, 36),  # fertile
        (35, 26, 18),  # rich -- dark, organic-and-mineral-rich earth
    ],
    dtype=float,
)
SOIL_OCEAN_BACKDROP_RGB = np.array([18, 28, 55], dtype=float)

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

# The Biome/Combined/Resources/Soil-Quality views' own render grid (see _biome_fields/
# _resource_fields) -- a fixed-shape equirectangular grid, like climate.py's native one, but
# swept at roughly GRID_SPACING_KM resolution (scaled by a world's own World.climate_density,
# see biome_grid_dimensions below) instead of climate's coarser native simulation grid, so a
# biome map reads with elevation-view-level detail even at climate_density's own default
# (1.0). This is purely a render-time upsample of climate.py's already-computed temperature/
# precipitation fields (see _bilinear_resample) plus a fresh nearest-node resample of the
# actual elevation data at this finer resolution for land/ocean/lake/glacier (see
# _biome_fields) -- the same technique _render_grid_arrays already uses for elevation.
BIOME_GRID_HEIGHT = round(np.pi / GRID_SPACING_RAD)
BIOME_GRID_WIDTH = round(2 * np.pi / GRID_SPACING_RAD)


def biome_grid_dimensions(climate_density: float) -> tuple[int, int]:
    """(height, width) for the Biome/Combined/Resources/Soil-Quality views' own render grid,
    scaled by a world's own World.climate_density the same direct-per-dimension way
    climate.grid_dimensions scales climate.py's native simulation grid -- density=2.0 means
    literally double BIOME_GRID_HEIGHT and double BIOME_GRID_WIDTH (4x the cells), matching
    the UI's own "double the density in each dimension" framing. BIOME_GRID_HEIGHT/WIDTH
    above stay as the density=1.0 reference values (also still the exact grid used by any
    caller not aware of climate_density -- there is none left after this change, but keeping
    them as plain names rather than folding them into this function preserves the "reference
    value the runtime one is scaled from" precedent plates.TARGET_LINE_SPACING_RAD and
    climate.GRID_HEIGHT/WIDTH both already set)."""
    spacing_rad = GRID_SPACING_RAD / climate_density
    return round(np.pi / spacing_rad), round(2 * np.pi / spacing_rad)

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
# pale rocky gray peaks. Kept in sync by hand with the palette this was ported from
# (frontend/src/elevationColor.ts, now deleted -- this is the sole copy). Deliberately no
# stop here is pure (255, 255, 255) white, even at the highest peaks -- white is reserved
# exclusively for actual ice cover (GLACIER_COLOR_RGB, applied as an overlay below) so it
# reads unambiguously as "glaciated," not "merely high/rocky."
_ELEVATION_STOP_E = np.array([-11000, -4000, -1500, -200, 0, 200, 1200, 3000, 6000, 9000], dtype=float)
_ELEVATION_STOP_RGB = np.array(
    [
        (10, 10, 40), (15, 40, 110), (40, 110, 190), (110, 170, 210), (200, 210, 150),
        (90, 150, 60), (170, 160, 90), (120, 90, 60), (195, 188, 178), (222, 217, 210),
    ],
    dtype=float,
)


def elevation_colors(elevations: np.ndarray, sea_level_m: float = 0.0) -> np.ndarray:
    """Vectorized elevation -> RGB: piecewise-linear interpolation through the hypsometric
    stops above, clamped at the ends (numpy.interp's default behavior already matches the
    old manual clamp). `sea_level_m` (World.sea_level_m, live-adjustable -- see
    World.sea_level_m/main.py's /world/controls) shifts every stop together, so the
    ocean/land color transition (the stop at 0) always tracks the *current* sea level rather
    than the fixed elevation value 0 -- the whole point of a real-time sea-level control is
    seeing the coastline visibly rise or recede on this same map."""
    channels = [np.interp(elevations - sea_level_m, _ELEVATION_STOP_E, _ELEVATION_STOP_RGB[:, c]) for c in range(3)]
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


def soil_fertility_colors(fertility: np.ndarray) -> np.ndarray:
    return _interp_colors(fertility, _SOIL_STOP_V, _SOIL_STOP_RGB)


# Clear (tan) -> silt-laden (murky brown), arbitrary concentration units (see ocean_cfd.py's
# own docstring -- sediment_concentration isn't calibrated against a real sediment density,
# only tuned so the pickup/settling balance produces a plausible-looking, slowly-varying
# field). Stops chosen empirically to keep a typical multi-step session's range legible.
_SEDIMENT_STOP_V = np.array([0.0, 0.5, 1.5, 3.0, 5.0], dtype=float)
_SEDIMENT_STOP_RGB = np.array([(150, 170, 160), (170, 150, 110), (150, 115, 70), (120, 85, 50), (90, 60, 35)], dtype=float)

# Bare seafloor (grey-teal) -> heavily deposited (dark ochre), meters -- ocean_cfd.py's
# sediment_deposited_m is a tracking-only field (see its own docstring: never mutates
# world.plates), so this is purely informational, not a hypsometric elevation scale.
_DEPOSITION_STOP_M = np.array([0.0, 0.1, 0.4, 1.0, 2.0], dtype=float)
_DEPOSITION_STOP_RGB = np.array([(60, 90, 95), (110, 95, 60), (140, 105, 45), (150, 90, 30), (120, 65, 20)], dtype=float)


def sediment_colors(concentration: np.ndarray) -> np.ndarray:
    return _interp_colors(concentration, _SEDIMENT_STOP_V, _SEDIMENT_STOP_RGB)


def sediment_deposition_colors(deposited_m: np.ndarray) -> np.ndarray:
    return _interp_colors(deposited_m, _DEPOSITION_STOP_M, _DEPOSITION_STOP_RGB)


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

    # At the default node_density, GRID_SPACING_RAD (100km) is already finer than the
    # physics resolution (plates.line_spacing_rad(1.0) = 125km), so it's the effective
    # resolution ceiling. At a higher node_density (e.g. 4x), the physics data itself gets
    # finer than 100km -- capping this render grid at the fixed 100km spacing regardless
    # made "Elevation" look visibly blockier than "Plates (details)" (which draws the raw
    # physics nodes directly, see render_png's `detail_lines`) even though both are drawing
    # the exact same underlying data. Taking the finer of the two keeps this view at least
    # as sharp as the data actually supports, without wasting resolution at the default
    # density where 100km was already the tighter bound.
    grid_spacing_rad = min(GRID_SPACING_RAD, plates.line_spacing_rad(world.node_density))

    xy_chunks, elev_chunks, owner_chunks, lake_chunks, glacier_chunks, volcano_chunks, hw_chunks, hh_chunks = [], [], [], [], [], [], [], []
    for phi, theta_candidates, world_pts in plates.iter_local_lattice(np.eye(3), spacing_rad=grid_spacing_rad):
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
        half_dtheta = grid_spacing_rad / max(np.cos(phi), 1e-3) / 2
        half_dphi = grid_spacing_rad / 2

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


def _biome_grid(height: int, width: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Same lat/lon/world_xyz construction as climate.py's own `_build_grid` (row 0 = north
    pole, increasing southward; column increasing eastward, wrapping) -- duplicated rather
    than imported since that's a private helper of climate.py's own native simulation grid,
    and this one is deliberately a different (finer) shape -- see BIOME_GRID_HEIGHT/WIDTH."""
    lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
    lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
    lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
    world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))
    return lat_deg, lon_deg, world_xyz


def _bilinear_resample(
    field: np.ndarray, src_lat_deg: np.ndarray, src_lon_deg: np.ndarray, dst_lat_deg: np.ndarray, dst_lon_deg: np.ndarray
) -> np.ndarray:
    """Upsamples a coarse (H, W) equirectangular field onto a finer (H2, W2) lat/lon grid by
    bilinear interpolation -- gives the Biome/Combined views' render grid (see
    BIOME_GRID_HEIGHT/WIDTH) a smoothly-varying temperature/precipitation field instead of
    climate.py's native grid's hard 2-degree cell boundaries, without re-running the
    simulation at a finer resolution. `src_lat_deg`/`src_lon_deg` are assumed uniformly
    spaced (true of climate.py's `_build_grid`, the only source this is ever called with).
    Latitude is clamped at the poles (not periodic); longitude wraps (periodic, matching the
    sphere's own topology) -- handled with explicit row/column index math rather than a
    generic scipy.interpolate call specifically so the two axes can differ that way."""
    src_h, src_w = field.shape
    row_step = 180.0 / src_h
    row_f = np.clip((90.0 - dst_lat_deg) / row_step - 0.5, 0, src_h - 1)
    row0 = np.floor(row_f).astype(int)
    row1 = np.clip(row0 + 1, 0, src_h - 1)
    row_t = (row_f - row0)[:, None]

    col_step = 360.0 / src_w
    lon_norm = ((dst_lon_deg + 180.0) % 360.0) - 180.0
    col_f = (lon_norm + 180.0) / col_step - 0.5
    col0 = np.floor(col_f).astype(int) % src_w
    col1 = (col0 + 1) % src_w
    col_t = (col_f - np.floor(col_f))[None, :]

    f00 = field[np.ix_(row0, col0)]
    f01 = field[np.ix_(row0, col1)]
    f10 = field[np.ix_(row1, col0)]
    f11 = field[np.ix_(row1, col1)]
    top = f00 * (1 - col_t) + f01 * col_t
    bot = f10 * (1 - col_t) + f11 * col_t
    return top * (1 - row_t) + bot * row_t


def _biome_fields(world: World, grid_h: int, grid_w: int):
    """The fine equirectangular grid shared by the Biome and Combined views: elevation/lake
    depth/glacier depth are a fresh nearest-node resample of the actual plate data at this
    grid's resolution (the same cKDTree technique _render_grid_arrays already uses for the
    Elevation/Plates views, via the same plates.collect_all_* helpers), so land/ocean/lake/
    glacier all line up with the Elevation view's own at matching detail; temperature/
    precipitation are bilinearly upsampled (_bilinear_resample) from climate.py's own
    coarser, fixed-shape simulation grid (see climate.compute_climate_cached) rather than
    resimulated at this resolution. Returns (lat_deg (H,), lon_deg (W,), world_xyz (H,W,3),
    elevation_m, is_ocean, air_temperature_c, ocean_temperature_c, precipitation_mm,
    lake_depth, glacier_depth), all (H, W) besides the first three."""
    lat_deg, lon_deg, world_xyz = _biome_grid(grid_h, grid_w)
    flat_xyz = world_xyz.reshape(-1, 3)
    shape = (grid_h, grid_w)

    collected = plates.collect_all_points(world.plates)
    if collected is None:
        elevation_m = np.zeros(shape)
        is_ocean = np.ones(shape, dtype=bool)
        lake_depth = np.zeros(shape)
        glacier_depth = np.zeros(shape)
    else:
        all_points, all_elevation, _ = collected
        all_lake_depth = plates.collect_all_lake_depth(world.plates)
        all_glacier_depth = plates.collect_all_glacier_depth(world.plates)
        tree = cKDTree(all_points)
        _, idx = tree.query(flat_xyz)
        elevation_m = all_elevation[idx].reshape(shape)
        lake_depth = all_lake_depth[idx].reshape(shape)
        glacier_depth = all_glacier_depth[idx].reshape(shape)
        is_ocean = elevation_m <= world.sea_level_m

    fields = climate.compute_climate_cached(world)
    air_temp = _bilinear_resample(fields.air_temperature_c, fields.lat_deg, fields.lon_deg, lat_deg, lon_deg)
    ocean_temp = _bilinear_resample(fields.ocean_temperature_c, fields.lat_deg, fields.lon_deg, lat_deg, lon_deg)
    precip = _bilinear_resample(fields.precipitation_mm, fields.lat_deg, fields.lon_deg, lat_deg, lon_deg)

    return lat_deg, lon_deg, world_xyz, elevation_m, is_ocean, air_temp, ocean_temp, precip, lake_depth, glacier_depth


def _resource_fields(world: World, grid_h: int, grid_w: int):
    """Fine equirectangular grid (same shape/construction as _biome_fields, see _biome_grid)
    for the Resources/Soil Quality views: coal/oil-gas/mineral deposit richness and soil
    depth/mineral-content/organic-content, each a fresh nearest-node resample of the actual
    plate data (plates.collect_all_* helpers) -- these views need none of climate.py's own
    fields, unlike Biome/Combined, so this is a separate, narrower resample rather than
    bloating _biome_fields' own return shape for views that don't need it. Both node-cloud
    fields (unlike climate.py's) are always defined, even before the first step -- soil/
    resources simply read as all-zero/barren then, same as any other freshly generated
    world's persistent fields (channel_depth, silt_depth, ...)."""
    lat_deg, lon_deg, world_xyz = _biome_grid(grid_h, grid_w)
    flat_xyz = world_xyz.reshape(-1, 3)
    shape = (grid_h, grid_w)

    collected = plates.collect_all_points(world.plates)
    if collected is None:
        z = np.zeros(shape)
        return lat_deg, lon_deg, world_xyz, np.ones(shape, dtype=bool), z, z, z, z, z, z

    all_points, all_elevation, _ = collected
    tree = cKDTree(all_points)
    _, idx = tree.query(flat_xyz)
    is_ocean = (all_elevation[idx].reshape(shape)) <= world.sea_level_m

    def resample(collector) -> np.ndarray:
        return collector(world.plates)[idx].reshape(shape)

    soil_depth = resample(plates.collect_all_soil_depth)
    soil_mineral = resample(plates.collect_all_soil_mineral_content)
    soil_organic = resample(plates.collect_all_soil_organic_content)
    coal = resample(plates.collect_all_coal_deposit)
    oil_gas = resample(plates.collect_all_oil_gas_deposit)
    mineral = resample(plates.collect_all_mineral_deposit)
    return lat_deg, lon_deg, world_xyz, is_ocean, soil_depth, soil_mineral, soil_organic, coal, oil_gas, mineral


def grid_slope(elevation_m: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """Dimensionless rise/run slope on the fine Biome/Combined grid (see _biome_fields) --
    real elevation difference to each cell's north/south or east/west neighbor (whichever is
    steeper), divided by that neighbor's real great-circle spacing in meters (longitude
    narrowed by cos(lat), same convention as everywhere else in this codebase -- see
    plates.iter_local_lattice). Feeds biomes.classify_wetland's own WETLAND_MAX_SLOPE cutoff --
    the same threshold erosion.compute_slope's own node-cloud slope (a different, finer
    discretization) is tuned against; see biomes.py's own module docstring for why an
    approximate, visually-tuned cutoff, not fit to any dataset, is this codebase's norm.
    np.roll wraps at the poles too (a minor, visually inconsequential artifact right at the
    map's own poles), the same "not worth special-casing" tradeoff this codebase already
    accepts elsewhere (e.g. the Plate Inspector's antipodal-projection limitation)."""
    grid_h, grid_w = elevation_m.shape
    dlat_km = (np.pi / grid_h) * plates.PLANET_RADIUS_KM
    dlon_km = np.maximum((2 * np.pi / grid_w) * plates.PLANET_RADIUS_KM * np.cos(np.radians(lat_deg))[:, None], 1.0)
    d_ns = np.abs(elevation_m - np.roll(elevation_m, 1, axis=0)) / (dlat_km * 1000.0)
    d_ew = np.abs(elevation_m - np.roll(elevation_m, 1, axis=1)) / (dlon_km * 1000.0)
    return np.maximum(d_ns, d_ew)


def _project_climate_grid(
    lat_deg: np.ndarray, lon_deg: np.ndarray, world_xyz: np.ndarray,
    projection: str, view_rotation: np.ndarray, width: int, height: int, padding_px: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float]:
    """Projects a full-sphere equirectangular (H, W, 3) grid to pixel space: fits the view's
    scale/offset to the whole sphere's projected extent (identical regardless of rotation --
    rotation only permutes which physical points land at which lat/lon, never changes the set
    of lat/lon values a full sphere covers), then returns each cell's pixel center and its
    per-cell axis-aligned half-extent in pixels (corner-based, not edge-midpoint-based --
    see _render_grid_arrays' own docstring for why that's required once the view can rotate).
    Shared by every CLIMATE_VIEWS renderer and the Biome/Combined views' own finer grid (see
    BIOME_GRID_HEIGHT/WIDTH) so they all place the sphere identically and never rescale or
    re-center relative to one another when the user switches views."""
    grid_h, grid_w = world_xyz.shape[:2]
    flat_xyz = world_xyz.reshape(-1, 3)

    bbox_xy = _project_points(projection, flat_xyz)
    flat_xy = _project_points(projection, _rotate(flat_xyz, view_rotation))

    min_x, min_y = bbox_xy.min(axis=0)
    max_x, max_y = bbox_xy.max(axis=0)
    data_w = max(max_x - min_x, 1e-9)
    data_h = max(max_y - min_y, 1e-9)
    scale = min((width - 2 * padding_px) / data_w, (height - 2 * padding_px) / data_h)
    offset_x = width / 2 - scale * (min_x + max_x) / 2
    offset_y = height / 2 + scale * (min_y + max_y) / 2

    centers = _to_pixels(scale, offset_x, offset_y, flat_xy)

    xy_grid = flat_xy.reshape(grid_h, grid_w, 2)
    _, center_lon_grid = geometry.xyz_to_latlon(_rotate(flat_xyz, view_rotation))
    center_lon_grid = center_lon_grid.reshape(grid_h, grid_w)

    half_dlat_deg = 180.0 / grid_h / 2
    half_dlon_deg = 360.0 / grid_w / 2
    lat_grid_deg = np.repeat(lat_deg[:, None], grid_w, axis=1)
    lon_grid_deg = np.repeat(lon_deg[None, :], grid_h, axis=0)

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

    return centers, half_w, half_h, scale, offset_x, offset_y


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

    outline_world = plate.get_bounding_polygon()
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


# Arrows are drawn at a coarser subsample of climate.py's own computation grid (90x180 at the
# default World.climate_density, finer if a world was generated at a higher one -- a fixed
# stride, so a finer grid draws proportionally more arrows, not fewer/thinner ones) -- one
# arrow per cell would be unreadable clutter at these sizes.
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


def _render_climate_view(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """Renders one of CLIMATE_VIEWS (besides "biome", which has its own finer-grid path --
    see _render_biome_view) from climate.py's own fixed grid -- a separate path from the
    plate-tectonics views below since the data source (a real (H, W) array, always covering
    the whole sphere) is structurally different from the render grid's ragged lattice, so
    there's little to share beyond the pixel-space primitives."""
    if view == "biome":
        return _render_biome_view(world, projection, width, height, view_rotation)

    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    blank = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    # Reuses whatever erosion.py already computed this step (see World.climate_cache)
    # instead of triggering a second ~50ms recomputation, rather than calling
    # climate.compute_climate directly.
    fields = climate.compute_climate_cached(world)

    pixels = blank.copy()
    centers, half_w, half_h, scale, offset_x, offset_y = _project_climate_grid(
        fields.lat_deg, fields.lon_deg, fields.world_xyz, projection, view_rotation, width, height, padding_px
    )

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
    draw = ImageDraw.Draw(image)

    if view in ("temperature", "humidity", "precipitation"):
        _draw_coastline(draw, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    if view == "wind":
        _draw_climate_vectors(draw, fields, fields.wind_u, fields.wind_v, projection, scale, offset_x, offset_y, pixel_scale, WIND_ARROW_COLOR, view_rotation)
    elif view == "oceanCurrents":
        _draw_climate_vectors(draw, fields, fields.current_u, fields.current_v, projection, scale, offset_x, offset_y, pixel_scale, CURRENT_ARROW_COLOR, view_rotation)
        _draw_swell_markers(draw, fields, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    # No server-side legend here -- see frontend/src/legendData.ts/Legend.tsx, which renders
    # it as a client-side HTML overlay instead (keyed purely on `view`, since none of this
    # module's legend content is actually data-dependent). Baking it into the PNG meant it
    # couldn't update mid-drag (the live-rotation preview draws a cheap client-side
    # graticule over the *last* rendered frame -- see MapCanvas.tsx) and meant every view's
    # legend text/gradient had to be hand-duplicated in Pillow drawing calls instead of
    # ordinary CSS/SVG.
    return _encode_image(image)


def _render_fluid_view(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """Renders one of OCEAN_CFD_VIEWS/ATMOSPHERE_CFD_VIEWS from the world's own permanent,
    always-on CFD state (World.ocean_cfd_state/atmosphere_cfd_state -- see ocean_cfd.py/
    atmosphere_cfd.py), reusing the same _project_climate_grid/_fill_rects/_draw_climate_vectors
    primitives _render_climate_view already uses for climate.py's own grid -- both state
    objects share the same (lat_deg, lon_deg, world_xyz) grid-geometry shape convention as
    climate.ClimateFields. Neither state is ever actually None once a world has been generated
    (see World.ocean_cfd_state's own docstring), but this function still degrades to a plain
    background-only image if it somehow is (same "always renders *something* standalone"
    contract every other view in VIEWS already has -- e.g. an elevation/plates render with no
    plates yet -- rather than raising)."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    pixels = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    is_ocean_view = view in OCEAN_CFD_VIEWS
    state = world.ocean_cfd_state if is_ocean_view else world.atmosphere_cfd_state
    if state is None:
        return _encode_image(Image.fromarray(pixels, mode="RGB"))

    centers, half_w, half_h, scale, offset_x, offset_y = _project_climate_grid(
        state.lat_deg, state.lon_deg, state.world_xyz, projection, view_rotation, width, height, padding_px
    )

    if view == "oceanCfdVelocity":
        backdrop = np.where(state.is_ocean.reshape(-1)[:, None], CLIMATE_OCEAN_BACKDROP_RGB, CLIMATE_LAND_BACKDROP_RGB)
        _fill_rects(pixels, centers, half_w, half_h, backdrop)
    elif view == "oceanCfdTemperature":
        _fill_rects(pixels, centers, half_w, half_h, temperature_colors(state.temperature_c.reshape(-1)))
    elif view == "oceanCfdSediment":
        _fill_rects(pixels, centers, half_w, half_h, sediment_colors(state.sediment_concentration.reshape(-1)))
    elif view == "oceanCfdDeposition":
        _fill_rects(pixels, centers, half_w, half_h, sediment_deposition_colors(state.sediment_deposited_m.reshape(-1)))
    elif view == "atmosphereCfdVelocity":
        backdrop = np.where(state.is_ocean.reshape(-1)[:, None], CLIMATE_OCEAN_BACKDROP_RGB, CLIMATE_LAND_BACKDROP_RGB)
        _fill_rects(pixels, centers, half_w, half_h, backdrop)
    elif view == "atmosphereCfdTemperature":
        _fill_rects(pixels, centers, half_w, half_h, temperature_colors(state.temperature_c.reshape(-1)))
    elif view == "atmosphereCfdHumidity":
        _fill_rects(pixels, centers, half_w, half_h, humidity_colors(state.humidity.reshape(-1)))

    image = Image.fromarray(pixels, mode="RGB")
    draw = ImageDraw.Draw(image)

    if view in ("oceanCfdTemperature", "oceanCfdSediment", "oceanCfdDeposition", "atmosphereCfdTemperature", "atmosphereCfdHumidity"):
        _draw_coastline(draw, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    if view == "oceanCfdVelocity":
        _draw_climate_vectors(draw, state, state.u, state.v, projection, scale, offset_x, offset_y, pixel_scale, CURRENT_ARROW_COLOR, view_rotation)
    elif view == "atmosphereCfdVelocity":
        _draw_climate_vectors(draw, state, state.u, state.v, projection, scale, offset_x, offset_y, pixel_scale, WIND_ARROW_COLOR, view_rotation)

    # No server-side legend here -- same reasoning as _render_climate_view's own trailing
    # comment (see frontend/src/legendData.ts/Legend.tsx).
    return _encode_image(image)


def _render_biome_view(world: World, projection: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """Biome, unlike the rest of CLIMATE_VIEWS, is rendered on its own much finer grid (see
    BIOME_GRID_HEIGHT/WIDTH and _biome_fields) rather than climate.py's native 90x180
    simulation grid directly, so its coastlines and color boundaries read at roughly
    Elevation-view resolution instead of climate's coarser simulation grid."""
    padding_px = PADDING_PX * (width / REFERENCE_WIDTH_PX)
    pixels = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    lat_deg, lon_deg, world_xyz, elevation_m, is_ocean, air_temp, ocean_temp, precip, _, _ = _biome_fields(
        world, *biome_grid_dimensions(world.climate_density)
    )
    display_temp = np.where(is_ocean, ocean_temp, air_temp)
    slope = grid_slope(elevation_m, lat_deg)
    biome_ids = biomes.classify_biomes(
        display_temp.reshape(-1), precip.reshape(-1), elevation_m.reshape(-1), slope.reshape(-1), is_ocean.reshape(-1), world.sea_level_m
    )
    colors = biomes.BIOME_COLORS[biome_ids]

    centers, half_w, half_h, _, _, _ = _project_climate_grid(
        lat_deg, lon_deg, world_xyz, projection, view_rotation, width, height, padding_px
    )
    _fill_rects(pixels, centers, half_w, half_h, colors)

    # No server-side legend/coastline overlay here -- same reasoning as _render_climate_view's
    # own trailing comment (see there); biome's colors already carry a land/ocean cue on
    # their own (Ocean is always the same fixed color), unlike temperature/humidity/
    # precipitation's scales.
    return _encode_image(Image.fromarray(pixels, mode="RGB"))


# Land colors blend toward the elevation-hypsometric shade by up to this fraction at the
# highest elevations -- a cheap "relief" cue (mountains read visibly lighter/rockier, matching
# real satellite natural-color imagery) without a full hillshade, which would need neighbor
# gradients the Biome/Combined grid's flat per-point sampling doesn't have on hand.
RELIEF_BLEND_MAX = 0.55
# Elevation (m above sea level) at which the relief blend reaches its max -- roughly
# ELEVATION_GRADIENT's own high-mountain stop, so full blend only kicks in near real peaks.
RELIEF_ELEVATION_RANGE_M = 6000.0

# Combined is the only view built from flat per-cell rectangles (_fill_rects) meant to read
# as a continuous true-color image rather than a legible data mosaic (Biome/Temperature/etc.
# lean into their own hard cell edges -- a viewer needs to tell one cell's exact color from
# its neighbor's). A light post-fill Gaussian blur softens those cell-edge jaggies into
# smooth coastlines/biome boundaries -- cheap anti-aliasing, in the same "approximate cue
# over an expensive exact one" spirit as biomes.biome_relative_shade_factor's own tiling
# above -- scaled by pixel_scale so it looks the same relative amount of soft at any requested
# resolution. Applied before rivers are drawn so their own lines stay crisp on top.
COMBINED_BLUR_RADIUS_PX = 1.0


def _render_combined_view(world: World, projection: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """"Combined": biome color for land, hypsometric ocean-depth shading for water (reusing
    elevation_colors, the same gradient the Elevation view itself uses) -- an approximation
    of what the planet would look like in true color from orbit, on the same fine grid the
    Biome view uses (see BIOME_GRID_HEIGHT/WIDTH and _biome_fields). Land color is first
    shaded by each cell's elevation *rank among other cells of the same biome* for a relief
    cue that still shows up even for biomes confined to a narrow absolute elevation band (see
    biomes.biome_relative_shade_factor), then blended toward that same hypsometric shade at
    high elevation for a further cue at real peaks (see RELIEF_BLEND_MAX); lakes/glaciers are
    overlaid the same way the Elevation view itself draws them, and rivers are drawn on top
    the same way too (see _draw_rivers), all at this grid's own resolution."""
    pixel_scale = width / REFERENCE_WIDTH_PX
    padding_px = PADDING_PX * pixel_scale
    pixels = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    lat_deg, lon_deg, world_xyz, elevation_m, is_ocean, air_temp, ocean_temp, precip, lake_depth, glacier_depth = _biome_fields(
        world, *biome_grid_dimensions(world.climate_density)
    )
    display_temp = np.where(is_ocean, ocean_temp, air_temp)
    slope = grid_slope(elevation_m, lat_deg)
    biome_ids = biomes.classify_biomes(
        display_temp.reshape(-1), precip.reshape(-1), elevation_m.reshape(-1), slope.reshape(-1), is_ocean.reshape(-1), world.sea_level_m
    )
    biome_rgb = biomes.BIOME_COLORS[biome_ids].astype(float)
    terrain_rgb = elevation_colors(elevation_m.reshape(-1), world.sea_level_m).astype(float)

    shade = biomes.biome_relative_shade_factor(biome_ids, elevation_m.reshape(-1))[:, None]
    shaded_biome_rgb = np.clip(biome_rgb * shade, 0, 255)

    relief_t = np.clip((elevation_m.reshape(-1) - world.sea_level_m) / RELIEF_ELEVATION_RANGE_M, 0.0, 1.0)
    blend = (relief_t * RELIEF_BLEND_MAX)[:, None]
    land_rgb = shaded_biome_rgb * (1 - blend) + terrain_rgb * blend

    colors = np.where(is_ocean.reshape(-1)[:, None], terrain_rgb, land_rgb)
    is_lake = lake_depth.reshape(-1) > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    if np.any(is_lake):
        colors = np.where(is_lake[:, None], np.array(LAKE_COLOR_RGB, dtype=float), colors)
    is_glacier = glacier_depth.reshape(-1) > hydrology.GLACIER_VISIBLE_DEPTH_M
    if np.any(is_glacier):
        colors = np.where(is_glacier[:, None], np.array(GLACIER_COLOR_RGB, dtype=float), colors)
    colors = np.clip(np.round(colors), 0, 255).astype(np.uint8)

    centers, half_w, half_h, scale, offset_x, offset_y = _project_climate_grid(
        lat_deg, lon_deg, world_xyz, projection, view_rotation, width, height, padding_px
    )
    _fill_rects(pixels, centers, half_w, half_h, colors)

    image = Image.fromarray(pixels, mode="RGB").filter(ImageFilter.GaussianBlur(radius=COMBINED_BLUR_RADIUS_PX * pixel_scale))
    image = _draw_rivers(image, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)

    return _encode_image(image)


def _blend(backdrop: np.ndarray, color: np.ndarray, fraction: np.ndarray) -> np.ndarray:
    """Per-cell linear blend of a flat `color` into `backdrop`, weighted by `fraction`
    ([0, 1], one value per cell) -- the same "clamped-fraction color blend" shape
    _land_shade_factor's own relief blend in _render_combined_view already uses, reused here
    for the Resources view's sequential coal/oil-gas-then-minerals overlay (see
    _render_resources_view)."""
    return backdrop * (1.0 - fraction[:, None]) + color[None, :] * fraction[:, None]


def _render_resources_view(is_ocean: np.ndarray, coal: np.ndarray, oil_gas: np.ndarray, mineral: np.ndarray) -> np.ndarray:
    """Categorical-ish overlay: a muted, low-saturation land/ocean backdrop, then coal (land)
    or oil & gas (ocean -- coal and oil & gas never spatially overlap, since one is strictly
    land-only and the other strictly ocean-only, see geology.py) blended in by richness
    fraction, then minerals blended on top last (can co-occur with either, since volcanism
    isn't restricted by crust type) -- drawing minerals last makes the rarer, more "exciting"
    deposit visually prominent wherever it does co-occur, the same "later layer wins where it
    applies" precedent render_image.py's own lake-before-glacier draw order already sets."""
    backdrop = np.where(is_ocean.reshape(-1)[:, None], RESOURCE_OCEAN_BACKDROP_RGB[None, :], RESOURCE_LAND_BACKDROP_RGB[None, :])
    coal_t = np.clip(coal.reshape(-1) / geology.MAX_COAL_DEPOSIT_M, 0.0, 1.0)
    oil_gas_t = np.clip(oil_gas.reshape(-1) / geology.MAX_OIL_GAS_DEPOSIT_M, 0.0, 1.0)
    mineral_t = np.clip(mineral.reshape(-1) / volcanism.MAX_MINERAL_DEPOSIT_M, 0.0, 1.0)

    colors = _blend(backdrop, COAL_COLOR_RGB, coal_t)
    colors = _blend(colors, OIL_GAS_COLOR_RGB, oil_gas_t)
    colors = _blend(colors, MINERAL_COLOR_RGB, mineral_t)
    return np.clip(np.round(colors), 0, 255).astype(np.uint8)


def _render_soil_view(is_ocean: np.ndarray, soil_depth: np.ndarray, soil_mineral: np.ndarray, soil_organic: np.ndarray) -> np.ndarray:
    """Continuous heatmap: fertility = sqrt(mineral * organic), the same "richest soil needs
    *both* high mineral and high organic content" scoring geology.py's own module docstring
    describes -- a geometric mean rewards having both far more than either alone, unlike a
    plain average. Zeroed wherever there's no soil at all (bare rock has nothing to hold
    either component, regardless of what the now-physically-meaningless relaxed fractions
    say) or over ocean (soil is a land-only concept)."""
    has_soil = soil_depth.reshape(-1) > 0.0
    fertility = np.sqrt(np.clip(soil_mineral.reshape(-1), 0.0, 1.0) * np.clip(soil_organic.reshape(-1), 0.0, 1.0))
    fertility = np.where(has_soil & ~is_ocean.reshape(-1), fertility, 0.0)
    colors = soil_fertility_colors(fertility)
    return np.where(is_ocean.reshape(-1)[:, None], SOIL_OCEAN_BACKDROP_RGB.astype(np.uint8)[None, :], colors)


def _render_resource_view(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray) -> bytes:
    """Renders "resources" or "soilQuality" (see RESOURCE_VIEWS) from the fine node-cloud
    resample _resource_fields provides -- structurally like _render_biome_view/
    _render_combined_view (same _biome_grid/_project_climate_grid machinery), but on data that
    exists independently of climate.py entirely, so (unlike Biome) this renders sensibly even
    before the first step has ever run."""
    padding_px = PADDING_PX * (width / REFERENCE_WIDTH_PX)
    pixels = np.full((height, width, 3), BACKGROUND_RGB, dtype=np.uint8)

    lat_deg, lon_deg, world_xyz, is_ocean, soil_depth, soil_mineral, soil_organic, coal, oil_gas, mineral = _resource_fields(
        world, *biome_grid_dimensions(world.climate_density)
    )
    if view == "resources":
        colors = _render_resources_view(is_ocean, coal, oil_gas, mineral)
    else:
        colors = _render_soil_view(is_ocean, soil_depth, soil_mineral, soil_organic)

    centers, half_w, half_h, scale, offset_x, offset_y = _project_climate_grid(
        lat_deg, lon_deg, world_xyz, projection, view_rotation, width, height, padding_px
    )
    _fill_rects(pixels, centers, half_w, half_h, colors)

    image = Image.fromarray(pixels, mode="RGB")
    if view == "soilQuality":
        # A continuous color scale carries no land/ocean cue on its own, same reasoning
        # temperature/humidity/precipitation already have -- resources' own backdrop already
        # distinguishes land/ocean directly, so it doesn't need this.
        draw = ImageDraw.Draw(image)
        _draw_coastline(draw, world, projection, scale, offset_x, offset_y, width / REFERENCE_WIDTH_PX, view_rotation)

    return _encode_image(image)


def _draw_rivers(
    image: Image.Image, world: World, projection: str, scale: float, offset_x: float, offset_y: float, pixel_scale: float, view_rotation: np.ndarray
) -> Image.Image:
    """Draws each river node's edge to its own downstream flow target as a short line
    segment (see hydrology.py's is_river/flow_target) -- reads world.hydrology_cache
    directly, populated by erosion.py every step (None before the world has ever been
    stepped, in which case this returns `image` unchanged). Each segment is a real, short 3D
    hop between two adjacent-in-the-flow-graph nodes, so _project_offset (not two independent
    _project_points calls) keeps it from being wrongly split across the antimeridian seam --
    same technique _render_grid_arrays' own corner measurements already rely on. Also cut by
    river_draw_min_flow(world) -- see that function's own docstring for why this view is
    stricter than the River Inspector's own /world/rivers listing, and why the floor itself
    isn't a single fixed constant. Line width steps directly off each segment's own
    flow_accum in whole multiples of that same floor -- exactly at the floor a segment is
    1 pixel wide, at 2x the floor it's 2 pixels, at 3x (RIVER_LINE_WIDTH_MAX_TIERS) or beyond
    it's 3 and no wider. flow_accum is constant along any unbranched stretch of channel (it
    only grows where a tributary's own flow actually merges in), so this widens a river only
    at real confluences -- never gradually along a single reach -- and only downstream of
    them, so it reads narrowest at the head and widest toward the mouth.

    Takes and returns a plain Image (rather than drawing onto a caller-owned
    ImageDraw.ImageDraw, like _draw_coastline and the other _draw_* helpers do) because
    antialiasing the lines here needs the mask-blur-then-composite done in pixel-array space --
    see RIVER_BLUR_RADIUS_PX. Callers must rebind their own ImageDraw to the returned image
    before drawing anything else on top."""
    hydro = world.hydrology_cache
    if hydro is None:
        return image
    min_flow = river_draw_min_flow(world)
    river_idx = np.nonzero(hydro.is_river & (hydro.flow_target >= 0) & (hydro.flow_accum >= min_flow))[0]
    if len(river_idx) == 0:
        return image
    target_idx = hydro.flow_target[river_idx]

    width_tier = np.clip(np.floor(hydro.flow_accum[river_idx] / min_flow), 1.0, RIVER_LINE_WIDTH_MAX_TIERS)
    width_px = np.maximum(np.round(width_tier * RIVER_LINE_WIDTH_PX * pixel_scale).astype(int), 1)

    from_points = _rotate(hydro.points[river_idx], view_rotation)
    to_points = _rotate(hydro.points[target_idx], view_rotation)
    _, from_lon = geometry.xyz_to_latlon(from_points)
    from_xy = _project_points(projection, from_points)
    to_xy = _project_offset(projection, to_points, from_lon)

    from_px = _to_pixels(scale, offset_x, offset_y, from_xy)
    to_px = _to_pixels(scale, offset_x, offset_y, to_xy)

    mask = Image.new("L", image.size, 0)
    mask_draw = ImageDraw.Draw(mask)
    for (x1, y1), (x2, y2), w in zip(from_px, to_px, width_px):
        mask_draw.line([(x1, y1), (x2, y2)], fill=255, width=int(w))
    mask = mask.filter(ImageFilter.GaussianBlur(radius=RIVER_BLUR_RADIUS_PX * pixel_scale))

    alpha = (np.asarray(mask, dtype=np.float32) / 255.0)[:, :, None]
    base_rgb = np.asarray(image, dtype=np.float32)
    blended = base_rgb * (1.0 - alpha) + np.array(RIVER_COLOR_RGB, dtype=np.float32) * alpha
    return Image.fromarray(np.clip(np.round(blended), 0, 255).astype(np.uint8), mode="RGB")


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
    if view in FLUID_VIEWS:
        return _render_fluid_view(world, projection, view, width, height, view_rotation)
    if view == "combined":
        return _render_combined_view(world, projection, width, height, view_rotation)
    if view in RESOURCE_VIEWS:
        return _render_resource_view(world, projection, view, width, height, view_rotation)

    if not world.plates:
        return _encode_image(Image.fromarray(blank, mode="RGB"))

    grid = _render_grid_arrays(world, projection, view_rotation) if view in ("elevation", "plates") else None
    tectonics = {p.plate_id: _plate_tectonics(projection, p, view_rotation) for p in world.plates}

    detail_lines = []  # (projected_xy, elevation) per non-empty line, "platesDetail" only
    if view == "platesDetail":
        for plate in world.plates:
            for line in plate.lines:
                if len(line) == 0:
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
        colors = elevation_colors(elev, world.sea_level_m) if view == "elevation" else plate_colors(owner)
        # Baked directly into the raster rather than a separate overlay/toggle: always
        # visible, no separate overlay needed, and a lake (or a glacier) is meaningful on
        # every view that shows terrain at all (matches how ocean itself isn't specially
        # toggle-able either).
        # Volcano is skipped on "elevation" specifically -- its saturated red used to
        # overwrite the actual hypsometric elevation color at every volcanic node, defeating
        # the point of an elevation view there (see VOLCANO_COLOR_RGB's own callers below);
        # it's still drawn on "plates", a categorical view with no elevation information to
        # lose. Glacier drawn last so ice wins on the rare cell where more than one overlay
        # would apply (a lake that just froze this same step still shows lake_depth from the
        # previous render's snapshot for one extra frame in the worst case, and a volcano
        # cold enough to glaciate should read as ice-covered, not lava-red) -- ice is the
        # more physically current/visually dominant state there.
        is_lake = lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
        if np.any(is_lake):
            colors = np.where(is_lake[:, None], np.array(LAKE_COLOR_RGB, dtype=np.uint8), colors)
        if view != "elevation" and np.any(is_volcano):
            colors = np.where(is_volcano[:, None], np.array(VOLCANO_COLOR_RGB, dtype=np.uint8), colors)
        is_glacier = glacier_depth > hydrology.GLACIER_VISIBLE_DEPTH_M
        if np.any(is_glacier):
            colors = np.where(is_glacier[:, None], np.array(GLACIER_COLOR_RGB, dtype=np.uint8), colors)
        _fill_rects(pixels, centers, hw_px, hh_px, colors)

    if detail_lines:
        dot_radius = NODE_DOT_RADIUS_PX * pixel_scale
        for xy, elev in detail_lines:
            centers = _to_pixels(scale, offset_x, offset_y, xy)
            colors = elevation_colors(elev, world.sea_level_m)
            _fill_rects(pixels, centers, dot_radius, dot_radius, colors)

    image = Image.fromarray(pixels, mode="RGB")
    image = _draw_rivers(image, world, projection, scale, offset_x, offset_y, pixel_scale, view_rotation)
    draw = ImageDraw.Draw(image)

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

    # No server-side legend here -- see the "No server-side legend" comment in
    # _render_climate_view above; the same client-side overlay covers elevation/plates/
    # platesDetail too.
    return _encode_image(image)


def _encode_image(image: Image.Image) -> bytes:
    buf = io.BytesIO()
    image.save(buf, format="PNG")
    return buf.getvalue()


def render_png_base64(world: World, projection: str, view: str, width: int, height: int, view_rotation: np.ndarray | None = None) -> str:
    return base64.b64encode(render_png(world, projection, view, width, height, view_rotation)).decode("ascii")


# How long each animation frame is shown, milliseconds -- matches frontend/src/App.tsx's own
# PLAY_INTERVAL_MS, so a saved GIF plays back at the same pace clicking Play already does.
ANIMATION_FRAME_DURATION_MS = 400


def render_animation_gif(
    world: World,
    projection: str,
    view: str,
    width: int,
    height: int,
    view_rotation: np.ndarray,
    years_per_frame: float,
    num_frames: int,
) -> bytes:
    """Renders an animated GIF of `world`'s progress in `view`/`projection`: frame 0 is the
    world's current state, and each of the `num_frames - 1` frames after it is
    `years_per_frame` further along -- calling step_world for real between frames, so this
    permanently advances `world` by `(num_frames - 1) * years_per_frame` years total (see
    main.py's `/world/animate` -- deliberately not a side-effect-free preview, same
    "the map really did move forward" semantics manually clicking Step that many times
    would have). Every frame is quantized against the *first* frame's own color palette
    rather than picking its own adaptive palette independently, which would otherwise make
    static regions (ocean, unchanged coastline) visibly flicker between playback frames --
    a well-known GIF-encoding pitfall, not the deliberately-changing regions this animation
    exists to show."""
    frames = []
    for i in range(num_frames):
        if i > 0:
            step_world(world, years_per_frame)
        png_bytes = render_png(world, projection, view, width, height, view_rotation)
        frames.append(Image.open(io.BytesIO(png_bytes)).convert("RGB"))

    reference_palette = frames[0].convert("P", palette=Image.ADAPTIVE, colors=256)
    quantized = [f.quantize(palette=reference_palette) for f in frames]

    buf = io.BytesIO()
    quantized[0].save(
        buf, format="GIF", save_all=True, append_images=quantized[1:], duration=ANIMATION_FRAME_DURATION_MS, loop=0
    )
    return buf.getvalue()


def render_animation_gif_base64(
    world: World,
    projection: str,
    view: str,
    width: int,
    height: int,
    view_rotation: np.ndarray,
    years_per_frame: float,
    num_frames: int,
) -> str:
    return base64.b64encode(
        render_animation_gif(world, projection, view, width, height, view_rotation, years_per_frame, num_frames)
    ).decode("ascii")
