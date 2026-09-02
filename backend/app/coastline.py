"""Coastline: the boundary between land-or-lake and open ocean, traced over climate.py's own
equirectangular grid.

**Why this needs its own module.** Several views (temperature/humidity/precipitation, and the
River Inspector) have no other visual cue for where land actually is: the elevation/plates
views show it implicitly through hypsometric coloring, and rivers draw their own segments
directly onto that same backdrop, but a climate-value color scale (or a blank canvas, for the
River Inspector) carries no land/ocean information on its own. This module computes exactly
one thing -- a set of line segments tracing that boundary -- that both a server-side PNG
render (render_image.py's `_draw_coastline`) and a client-rendered JSON view (main.py's
`GET /world/rivers`, consumed by RiverInspector.tsx) can each draw for themselves.

**Grid, not node cloud.** Traced over climate.py's own `(H, W)` grid (`world.climate_cache`,
reused rather than recomputed -- see climate.py's own caching philosophy), not the irregular
per-plate node cloud: a coastline is inherently about *land shape*, which requires deciding
"here" vs. "not here" at every point on the sphere, and climate's grid already covers the
whole sphere densely and uniformly, exactly what's needed. climate.py has no lake concept of
its own (see its module docstring's "deliberately not ported" list), so lake-ness is resampled
onto that same grid here, the same nearest-node `cKDTree` technique
`climate._sample_elevation_and_crust` already uses for elevation.

**A rectilinear trace, not a smoothed contour.** For every pair of horizontally- or
vertically-adjacent grid cells that disagree on land-or-lake vs. ocean, the shared edge
between them (half a cell-step in from each of their centers) is one coastline segment. This
deliberately matches the grid's own existing visual language everywhere else it's drawn
(render_image.py's `_fill_rects`, itself rectangular per cell) rather than introducing a new
marching-squares-style smoothing pass this codebase has no other use for.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, hydrology, plates

if TYPE_CHECKING:
    from .world import World


def _lake_mask_on_grid(world: "World", world_xyz: np.ndarray) -> np.ndarray:
    """Nearest-node resample of the current hydrology_cache's lake_depth onto `world_xyz`
    (H, W, 3) -- same cKDTree technique climate._sample_elevation_and_crust uses for
    elevation. All-False before the first step, or for a world with no lakes at all."""
    height, width, _ = world_xyz.shape
    hydro = world.hydrology_cache
    if hydro is None or len(hydro.points) == 0:
        return np.zeros((height, width), dtype=bool)
    tree = cKDTree(hydro.points)
    flat_xyz = world_xyz.reshape(-1, 3)
    _, idx = tree.query(flat_xyz, workers=plates.query_workers(len(flat_xyz)))
    lake_depth = hydro.lake_depth[idx].reshape(height, width)
    return lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M


def compute_coastline_segments(world: "World") -> tuple[np.ndarray, np.ndarray]:
    """Returns `(point_a, point_b)`, each `(K, 3)` world-space unit vectors -- one row per
    coastline edge, matching the segment-pair convention `hydrology`'s river edges already
    use (a flat edge list, not an ordered polyline: land/lake shapes aren't simple loops).
    Empty `(0, 3)` arrays before the first step (`climate_cache` is `None` until erosion.py
    runs once, same as `hydrology_cache` -- see `World.climate_cache`/`hydrology_cache`)."""
    fields = world.climate_cache
    if fields is None:
        empty = np.zeros((0, 3))
        return empty, empty

    height, width = fields.is_ocean.shape
    # Three mutually exclusive categories, not a two-way land-or-lake-vs-ocean split: a lake
    # entirely surrounded by land (the common case -- an inland lake, nowhere near open ocean)
    # still needs its own boundary drawn against that land, which a land-or-lake union would
    # never produce (both sides would read as the same "land-or-lake" value). `is_lake` is
    # masked to `~is_ocean` first so an ocean cell whose nearest hydrology node happens to be
    # a lake node (the two are separate nearest-neighbor resamples, so they can disagree right
    # at a coastline) never gets counted as a lake -- ocean always wins that conflict.
    is_lake = _lake_mask_on_grid(world, fields.world_xyz) & ~fields.is_ocean
    category = np.where(fields.is_ocean, 0, np.where(is_lake, 1, 2))  # 0=ocean, 1=lake, 2=land

    half_lat = 90.0 / height
    half_lon = 180.0 / width
    lat_deg, lon_deg = fields.lat_deg, fields.lon_deg

    a_latlon: list[np.ndarray] = []
    b_latlon: list[np.ndarray] = []

    # Vertical edges: column c vs. (c + 1) % width -- wraps at the antimeridian, since the
    # grid itself wraps there (see climate._build_grid's own convention).
    next_col = np.roll(category, -1, axis=1)
    rows, cols = np.nonzero(category != next_col)
    if len(rows) > 0:
        lon_edge = lon_deg[cols] + half_lon
        a_latlon.append(np.stack([lat_deg[rows] + half_lat, lon_edge], axis=1))
        b_latlon.append(np.stack([lat_deg[rows] - half_lat, lon_edge], axis=1))

    # Horizontal edges: row r vs. r + 1 -- no wraparound (the poles are the grid's own
    # top/bottom edge, not connected to each other).
    rows2, cols2 = np.nonzero(category[:-1, :] != category[1:, :])
    if len(rows2) > 0:
        lat_edge = lat_deg[rows2] - half_lat
        a_latlon.append(np.stack([lat_edge, lon_deg[cols2] - half_lon], axis=1))
        b_latlon.append(np.stack([lat_edge, lon_deg[cols2] + half_lon], axis=1))

    if not a_latlon:
        empty = np.zeros((0, 3))
        return empty, empty

    a = np.concatenate(a_latlon, axis=0)
    b = np.concatenate(b_latlon, axis=0)
    point_a = geometry.latlon_to_xyz(np.radians(a[:, 0]), np.radians(a[:, 1]))
    point_b = geometry.latlon_to_xyz(np.radians(b[:, 0]), np.radians(b[:, 1]))
    return point_a, point_b
