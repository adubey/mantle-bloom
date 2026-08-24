"""LatLongGrid: a fixed-shape equirectangular lat/long grid used as a representation-
agnostic bridge between a plate's own elevation-node storage (`PlateWithLines`' parallel
`ElevationLine`s or `PlateWithRTree`'s R-tree-indexed point cloud, see plates.py) and any
per-step pass that would rather read and modify terrain in plain grid space than deal with
either representation directly -- the same nearest-cell resample technique climate.py's own
grid already uses to sample plate elevation (`ClimateFields`/`_sample_elevation_and_crust`),
but bidirectional: `Plate.update_to_lat_long_grid` seeds a grid from a plate's current nodes,
arbitrary grid-space code reshapes it (`change_elevation`, below), and
`Plate.update_deltas_from_lat_long_grid` writes only the *net change* back onto that plate's
own nodes.

Only the delta round-trips, never the grid's absolute elevation: a grid cell coarser than a
plate's own node spacing resamples many nodes down to one value, so writing that value straight
back would flatten every node sharing a cell to the same elevation. Tracking the *change*
relative to each cell's own original (Plate-supplied) value instead means only genuinely new
grid-space modifications propagate back, at whatever resolution the plate's own nodes actually
have -- the same reasoning `plates.py`'s own `Plate.update_deltas_from_lat_long_grid` docstring
repeats from the plate side of this round trip.
"""

from __future__ import annotations

import numpy as np

from . import geometry
from .elevation_lines import MAX_ELEVATION_M, MIN_ELEVATION_M


class LatLongGrid:
    """Row 0 = north pole, row increases southward; column increases eastward, wraps -- the
    same layout convention climate.py's own grid uses (see that module's `_build_grid`), so a
    `LatLongGrid` and a `ClimateFields` built at the same (height, width) share cell indices."""

    def __init__(self, lat_deg: np.ndarray, lon_deg: np.ndarray, world_xyz: np.ndarray) -> None:
        self._lat_deg = lat_deg
        self._lon_deg = lon_deg
        self._world_xyz = world_xyz
        shape = (len(lat_deg), len(lon_deg))
        # Every cell starts at elevation 0 with no delta -- real values only arrive via
        # set_elevation (Plate.update_to_lat_long_grid), which stamps both arrays together.
        self._elevation = np.zeros(shape)
        self._original_elevation = np.zeros(shape)

    @classmethod
    def build(cls, height: int, width: int) -> "LatLongGrid":
        """A fresh (height, width) grid, every cell at elevation 0 -- populate it with real
        terrain via `Plate.update_to_lat_long_grid` (once per plate) before reading or
        changing it."""
        lat_deg = 90.0 - (np.arange(height) + 0.5) * (180.0 / height)
        lon_deg = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
        lat_grid = np.repeat(lat_deg[:, None], width, axis=1)
        lon_grid = np.repeat(lon_deg[None, :], height, axis=0)
        world_xyz = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))
        return cls(lat_deg, lon_deg, world_xyz)

    @property
    def height(self) -> int:
        return len(self._lat_deg)

    @property
    def width(self) -> int:
        return len(self._lon_deg)

    @property
    def lat_deg(self) -> np.ndarray:
        return self._lat_deg

    @property
    def lon_deg(self) -> np.ndarray:
        return self._lon_deg

    @property
    def world_xyz(self) -> np.ndarray:
        return self._world_xyz

    @property
    def elevation(self) -> np.ndarray:
        """Current (H, W) elevation: each cell's original Plate-supplied value plus whatever
        `change_elevation` has since added."""
        return self._elevation

    @property
    def elevation_delta(self) -> np.ndarray:
        """Net (H, W) change since this grid was last populated by a plate's own
        `update_to_lat_long_grid` -- exactly what `update_deltas_from_lat_long_grid` adds back
        onto that plate's nodes."""
        return self._elevation - self._original_elevation

    def row_col_for_world_xyz(self, world_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Nearest (row, col) for each world-space unit vector in `world_xyz` (..., 3) -- the
        same lat/lon -> cell-index rounding climate.py's own `_sample_at_offset` uses, so a
        point resamples to the same cell here as it would on a same-shape climate grid. Both
        `set_elevation`'s caller and `delta_at`'s caller go through this, so the same node
        always round-trips through the same cell."""
        lat, lon = geometry.xyz_to_latlon(world_xyz)
        lat_deg, lon_deg = np.degrees(lat), np.degrees(lon)
        height, width = self.height, self.width
        rows = np.clip(np.round((90.0 - lat_deg) / (180.0 / height) - 0.5).astype(np.int64), 0, height - 1)
        cols = np.round((lon_deg + 180.0) / (360.0 / width) - 0.5).astype(np.int64) % width
        return rows, cols

    def set_elevation(self, rows: np.ndarray, cols: np.ndarray, elevation: np.ndarray) -> None:
        """Seed cells (rows, cols) with `elevation`, a plate's own current node values --
        called by `Plate.update_to_lat_long_grid`. Also records the same value as each cell's
        *original* elevation, the baseline `elevation_delta` measures against; if two plates
        (or two nodes sharing one coarser cell) write to the same cell, the last write wins for
        both, same as any other nearest-neighbor resample onto a coarser grid (see climate.py's
        own `_sample_elevation_and_crust`)."""
        self._elevation[rows, cols] = elevation
        self._original_elevation[rows, cols] = elevation

    def change_elevation(self, delta: np.ndarray) -> None:
        """Add `delta` (H, W) to every cell's current elevation -- for grid-space code
        (erosion, orographic reshaping, ...) that wants to reshape terrain without knowing
        which plate, or which representation, actually owns whatever node(s) a given cell
        resamples. Each cell's *original* Plate-supplied elevation (`set_elevation`, above)
        is untouched, so `elevation_delta` always reflects the net change since this grid was
        last populated, no matter how many separate `change_elevation` calls contributed to
        it. Clipped to the same MIN_ELEVATION_M/MAX_ELEVATION_M bounds every other
        elevation-modifying pass in this codebase respects (see elevation_lines.py)."""
        self._elevation = np.clip(self._elevation + delta, MIN_ELEVATION_M, MAX_ELEVATION_M)

    def delta_at(self, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
        """`elevation_delta` at cells (rows, cols) -- what
        `Plate.update_deltas_from_lat_long_grid` adds onto each of a plate's own nodes."""
        return self.elevation_delta[rows, cols]
