"""Sparse plate-local quad patches -- the static quad surface of issue #228's Phase 2.

`PlateWithSparseQuadPatch` is a second implementation of the `PlateSurface` contract
(plates.py) alongside `PlateWithLines`. Its territory is a sparse set of active cells on an
equiangular cube-sphere lattice laid out in the plate's own local frame:

- Six faces, each an `n x n` grid of cells whose edges are great-circle arcs (lines of
  constant equiangular coordinate on a cube face are planes through the origin). Face 0 is
  centred on the plate's local (phi=0, theta=0) seed, so an ordinary plate lies mostly on
  one face; a plate of any size, including one owning a pole or more than a hemisphere, is
  still representable without a coordinate singularity. That answers the issue's "one patch
  or several charts?" question for Phase 2: one lattice, six charts, seams handled by
  construction because neighbouring faces' cells match edge for edge.
- Fields live on leaf cells (one node per active cell, at the cell's angular centre). Cells give
  every node an exact footprint (`SurfaceNodes.area_is_exact`), which is what Phase 3's
  conservative remapping needs, and cell-centred values are the finite-volume reading of
  the extensive fields in `surface_fields.SURFACE_FIELDS`.
- The cell set is authoritative. Adjacency, boundary loops, row/column intervals, node
  positions, and areas are derived from it and cached per topology revision; rigid rotation
  only changes `frame` and invalidates world-space caches (geometry revision), never the
  local ones. Phase 2 has no topology-changing operation at all -- remeshing is Phase 3 and
  tectonics Phase 4 -- so a static patch's `topology_revision` stays where it was built.

Cells are addressed by a packed `int64` key -- face, refinement level, row `j`, column `i`
-- whose sort order is the canonical node order, so the key doubles as the stable,
storage-independent node ID Phase 0a asked for. Refinement replaces a leaf by its four
children and coarsening reverses that operation. The key itself therefore records lineage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterator, Mapping

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
from .elevation_lines import (
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_OCEANIC,
    PLANET_RADIUS_KM,
    ElevationPoint,
    install_point_field_accessors,
)
from .plates import NEIGHBOUR_DISTANCE_RAD, Plate, SurfaceAdjacency, SurfaceNodes, _plates_within
from .surface_fields import SURFACE_FIELDS, RemapClass

# Bumped whenever the pickled state of a `PlateWithSparseQuadPatch` changes shape. A save
# carrying a newer version than this code understands is rejected on load rather than
# half-restored (see `__setstate__`).
QUAD_SURFACE_FORMAT_VERSION = 1

PLANET_RADIUS_M = PLANET_RADIUS_KM * 1000.0

# (normal, u-axis, v-axis) per face, in plate-local coordinates. Each triple is right-handed
# (u x v = normal), so counter-clockwise in a face's (u, v) grid is counter-clockwise seen
# from outside the sphere, and boundary loops come out with outer boundaries CCW and holes
# CW. Face 0 is centred on local +x (the plate's seed); its u/v axes are local east/north
# there, matching `PlateWithLines`' (theta, phi) orientation.
_FACE_AXES = np.array(
    [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
        [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
        [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
        [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
        [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    ],
    dtype=np.int64,
)
_FACE_NORMALS = _FACE_AXES[:, 0].astype(float)
_FACE_U = _FACE_AXES[:, 1].astype(float)
_FACE_V = _FACE_AXES[:, 2].astype(float)

# Key layout: face (3 bits) | level (8 bits) | j (24 bits) | i (24 bits). Row-major within
# a face, so sorting keys orders nodes face by face, row by row, column by column.
_I_BITS = 24
_J_SHIFT = _I_BITS
_LEVEL_SHIFT = 2 * _I_BITS
_FACE_SHIFT = _LEVEL_SHIFT + 8
_INDEX_MASK = (1 << _I_BITS) - 1

# How far past a shared cell edge (as a fraction of one cell) `_neighbour_keys` probes to
# find the cell on the other side. Small enough that crossing a cube-face seam, where the
# probe's straight-line path bends relative to the far face's grid, still lands in the
# matching row; large enough to be far above floating-point noise.
_NEIGHBOUR_PROBE_FRACTION = 0.1

# Direction order used by adjacency and boundary extraction: -v (below), +u (right),
# +v (above), -u (left) -- i.e. a cell's edges in counter-clockwise order.
_DIRECTION_OFFSETS = ((0.5, -_NEIGHBOUR_PROBE_FRACTION), (1.0 + _NEIGHBOUR_PROBE_FRACTION, 0.5), (0.5, 1.0 + _NEIGHBOUR_PROBE_FRACTION), (-_NEIGHBOUR_PROBE_FRACTION, 0.5))
# Each direction's exposed edge as (start corner, end corner) offsets from the cell's own
# (i, j) corner, traversed counter-clockwise.
_EDGE_CORNERS = (((0, 0), (1, 0)), ((1, 0), (1, 1)), ((1, 1), (0, 1)), ((0, 1), (0, 0)))


def cells_per_face_edge(spacing_rad: float) -> int:
    """Cells along one cube-face edge for a lattice whose *mean* cell area matches a line
    node's nominal `spacing_rad ** 2` footprint, so a quad world and a line world generated
    at the same `node_density` carry about the same number of nodes."""
    return max(1, int(round(np.sqrt(4.0 * np.pi / 6.0) / spacing_rad)))


def pack_cell_keys(face: np.ndarray, i: np.ndarray, j: np.ndarray, level: int | np.ndarray = 0) -> np.ndarray:
    face = np.asarray(face, dtype=np.int64)
    level = np.asarray(level, dtype=np.int64)
    return (
        (face << _FACE_SHIFT)
        | (level << _LEVEL_SHIFT)
        | (np.asarray(j, dtype=np.int64) << _J_SHIFT)
        | np.asarray(i, dtype=np.int64)
    )


def unpack_cell_keys(keys: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(face, level, i, j) for each packed key."""
    keys = np.asarray(keys, dtype=np.int64)
    return (
        keys >> _FACE_SHIFT,
        (keys >> _LEVEL_SHIFT) & 0xFF,
        keys & _INDEX_MASK,
        (keys >> _J_SHIFT) & _INDEX_MASK,
    )


def lattice_points(face: np.ndarray, a: np.ndarray, b: np.ndarray, n: int) -> np.ndarray:
    """Local unit vectors at fractional grid coordinates (`a` along u, `b` along v, each in
    cell units from the face's low corner) on `face`. Coordinates just outside `[0, n]` are
    allowed and land on the neighbouring face, which `_neighbour_keys` relies on."""
    step = (np.pi / 2.0) / n
    tan_a = np.tan(-np.pi / 4.0 + np.asarray(a, dtype=float) * step)
    tan_b = np.tan(-np.pi / 4.0 + np.asarray(b, dtype=float) * step)
    face = np.asarray(face, dtype=np.int64)
    vectors = _FACE_NORMALS[face] + tan_a[..., None] * _FACE_U[face] + tan_b[..., None] * _FACE_V[face]
    return vectors / np.linalg.norm(vectors, axis=-1, keepdims=True)


def locate_cells(local_xyz: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(face, i, j) of the level-0 cell containing each local unit vector."""
    local_xyz = np.asarray(local_xyz, dtype=float).reshape(-1, 3)
    axis = np.argmax(np.abs(local_xyz), axis=1)
    positive = local_xyz[np.arange(len(local_xyz)), axis] > 0.0
    # Axis/sign -> face index, inverting _FACE_AXES' normals.
    face_lookup = np.array([[2, 0], [3, 1], [5, 4]])
    face = face_lookup[axis, positive.astype(np.int64)]
    depth = np.einsum("ij,ij->i", local_xyz, _FACE_NORMALS[face])
    alpha = np.arctan(np.einsum("ij,ij->i", local_xyz, _FACE_U[face]) / depth)
    beta = np.arctan(np.einsum("ij,ij->i", local_xyz, _FACE_V[face]) / depth)
    step = (np.pi / 2.0) / n
    i = np.clip(np.floor((alpha + np.pi / 4.0) / step).astype(np.int64), 0, n - 1)
    j = np.clip(np.floor((beta + np.pi / 4.0) / step).astype(np.int64), 0, n - 1)
    return face, i, j


def cell_areas_sr(i: np.ndarray, j: np.ndarray, n: int) -> np.ndarray:
    """Exact solid angle (steradians) of each level-0 cell. On a face's gnomonic tangent
    plane, the solid angle of the rectangle [0, x] x [0, y] is atan(xy / sqrt(1 + x^2 + y^2)),
    so a cell is four evaluations of that by inclusion-exclusion. Face-independent: every
    face has the same grid."""
    step = (np.pi / 2.0) / n

    def tangent(k: np.ndarray) -> np.ndarray:
        return np.tan(-np.pi / 4.0 + np.asarray(k, dtype=float) * step)

    def corner(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        return np.arctan(x * y / np.sqrt(1.0 + x * x + y * y))

    x0, x1 = tangent(i), tangent(np.asarray(i) + 1)
    y0, y1 = tangent(j), tangent(np.asarray(j) + 1)
    return corner(x1, y1) - corner(x0, y1) - corner(x1, y0) + corner(x0, y0)


def parent_cell_keys(keys: np.ndarray) -> np.ndarray:
    """Immediate parent IDs, or ``-1`` for level-zero roots."""
    face, level, i, j = unpack_cell_keys(keys)
    return np.where(level == 0, -1, pack_cell_keys(face, i // 2, j // 2, level - 1))


def child_cell_keys(keys: np.ndarray) -> np.ndarray:
    """The four child IDs of each cell. Cell IDs themselves encode lineage."""
    face, level, i, j = unpack_cell_keys(keys)
    children = [
        pack_cell_keys(face, 2 * i + di, 2 * j + dj, level + 1)
        for di, dj in ((0, 0), (1, 0), (0, 1), (1, 1))
    ]
    return np.stack(children, axis=-1)


def _corner_keys(face: np.ndarray, ci: np.ndarray, cj: np.ndarray, n: int) -> np.ndarray:
    """A lattice-corner identity shared by every face that corner touches. Each corner maps
    to an integer point on the cube surface, `n * normal + (2ci - n) * u + (2cj - n) * v`:
    the equiangular parameterisation agrees on both sides of every cube edge, so the same
    corner seen from two or three faces gets the same integer point, with no float rounding
    involved."""
    face = np.asarray(face, dtype=np.int64)
    point = (
        n * _FACE_AXES[face, 0]
        + (2 * np.asarray(ci, dtype=np.int64) - n)[..., None] * _FACE_AXES[face, 1]
        + (2 * np.asarray(cj, dtype=np.int64) - n)[..., None] * _FACE_AXES[face, 2]
    )
    shifted = point + n
    return (shifted[..., 0] << 42) | (shifted[..., 1] << 21) | shifted[..., 2]


@dataclass(frozen=True)
class CellIntervals:
    """Runs of consecutive active cells along one lattice direction, sorted by
    (face, line, start). For row intervals `line` is the row `j` and `start`/`end` are
    inclusive column indices `i`; for column intervals, the reverse. A row cut by a notch,
    hole, or disconnected fragment simply contributes several runs -- nothing assumes one
    run per line. `node_run` maps each node (canonical order) to the run containing it."""

    face: np.ndarray
    line: np.ndarray
    start: np.ndarray
    end: np.ndarray
    node_run: np.ndarray
    level: np.ndarray


def _intervals(face: np.ndarray, line: np.ndarray, position: np.ndarray, level: np.ndarray | None = None) -> CellIntervals:
    level = np.zeros(len(face), dtype=np.int64) if level is None else np.asarray(level)
    order = np.lexsort((position, line, level, face))
    f, lev, l, p = face[order], level[order], line[order], position[order]
    breaks = np.ones(len(order), dtype=bool)
    if len(order) > 1:
        breaks[1:] = (f[1:] != f[:-1]) | (lev[1:] != lev[:-1]) | (l[1:] != l[:-1]) | (p[1:] != p[:-1] + 1)
    run_of_sorted = np.cumsum(breaks) - 1
    starts = np.flatnonzero(breaks)
    ends = np.append(starts[1:], len(order)) - 1
    node_run = np.empty(len(order), dtype=np.int64)
    node_run[order] = run_of_sorted
    return CellIntervals(face=f[starts], line=l[starts], start=p[starts], end=p[ends], node_run=node_run, level=lev[starts])


def _stitch_xyz_loops(loops: list[np.ndarray]) -> np.ndarray:
    """`plates._stitch_loops`' zero-width keyhole seams, over xyz loops and with a k-d tree
    nearest-pair search instead of an all-pairs scan (quad loops carry every boundary cell
    corner, so the quadratic scan would dominate outline construction on big plates)."""
    combined = loops[0]
    for extra in loops[1:]:
        dist, idx = cKDTree(combined).query(extra)
        j = int(np.argmin(dist))
        i = int(idx[j])
        rotated = np.concatenate([extra[j:], extra[:j]])
        combined = np.concatenate([combined[: i + 1], rotated, extra[j : j + 1], combined[i : i + 1], combined[i + 1 :]])
    return combined


@install_point_field_accessors
class ElevationPointInPatch:
    """One active cell of a `PlateWithSparseQuadPatch`, as a live `ElevationPoint` view --
    `set_*` writes straight into the plate's field arrays. Valid until the plate's topology
    changes, the same lifetime `ElevationPointOnLine` has against its line."""

    def __init__(self, plate: "PlateWithSparseQuadPatch", index: int) -> None:
        self._plate = plate
        self._index = index

    @property
    def phi(self) -> float:
        return float(self._plate._node_latlon()[0][self._index])

    def _field_array(self, name: str) -> np.ndarray:
        if name == "theta":
            return self._plate._node_latlon()[1]
        return self._plate._field_storage(name)


class PlateWithSparseQuadPatch(Plate):
    """A plate whose terrain is a sparse set of active cube-sphere cells -- see the module
    docstring. Static in Phase 2: it can be generated, queried, rendered, rigidly rotated,
    and saved, but not deformed, split, merged, or defragmented."""

    # Derived state rebuilt on demand from (`_n`, `_keys`, `_frame`); never pickled.
    _TOPOLOGY_CACHES = (
        "_local_cache",
        "_latlon_cache",
        "_area_cache",
        "_adjacency_cache",
        "_local_loops_cache",
        "_row_intervals_cache",
        "_column_intervals_cache",
    )
    _GEOMETRY_CACHES = (
        "_world_points_cache",
        "_bounding_polygon_cache",
        "_bounding_polygon_tree_cache",
        "_node_kdtree_cache",
    )

    def __init__(
        self,
        plate_id: int,
        frame: np.ndarray,
        crust_type: str,
        cells_per_edge: int,
        cell_keys: np.ndarray,
        fields: Mapping[str, np.ndarray] | None = None,
        omega: np.ndarray | None = None,
        age_steps: int = 0,
        internal_stress: float = 0.0,
    ) -> None:
        super().__init__(plate_id, frame, crust_type, omega=omega, age_steps=age_steps, internal_stress=internal_stress)
        if cells_per_edge < 1 or cells_per_edge > _INDEX_MASK:
            raise ValueError(f"cells_per_edge must be in [1, {_INDEX_MASK}], got {cells_per_edge}")
        keys = np.asarray(cell_keys, dtype=np.int64).reshape(-1)
        face, level, i, j = unpack_cell_keys(keys)
        max_level = int(np.floor(np.log2(_INDEX_MASK / cells_per_edge)))
        invalid_level = level > max_level
        resolution = np.array([cells_per_edge * (1 << int(value)) for value in level])
        if np.any((face < 0) | (face > 5) | invalid_level | (i >= resolution) | (j >= resolution)):
            raise ValueError("cell keys must address cells of this lattice")
        order = np.argsort(keys, kind="stable")
        keys = keys[order]
        if len(keys) > 1 and np.any(keys[1:] == keys[:-1]):
            raise ValueError("cell keys must be unique")
        self._n = int(cells_per_edge)
        self._keys = keys
        self._fields: dict[str, np.ndarray] = {}
        for name, values in (fields or {}).items():
            self._check_field(name, values, len(order))
            self._fields[name] = np.asarray(values, dtype=SURFACE_FIELDS[name].dtype)[order].copy()
        self._reset_caches()

    @classmethod
    def from_lattice(
        cls,
        plate_id: int,
        frame: np.ndarray,
        crust_type: str,
        spacing_rad: float,
        is_owned: Callable[[np.ndarray], np.ndarray],
        **kwargs,
    ) -> "PlateWithSparseQuadPatch":
        """Sweep this plate's whole local lattice and keep every cell whose centre
        `is_owned(world_pts)` claims -- the quad analogue of
        `elevation_lines.build_lines_from_lattice`, driven by the same ownership test."""
        n = cells_per_face_edge(spacing_rad)
        jj, ii = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
        ii, jj = ii.reshape(-1), jj.reshape(-1)
        kept: list[np.ndarray] = []
        for face in range(6):
            faces = np.full(len(ii), face)
            world_pts = geometry.to_world(frame, lattice_points(faces, ii + 0.5, jj + 0.5, n))
            owned = np.asarray(is_owned(world_pts), dtype=bool)
            kept.append(pack_cell_keys(faces[owned], ii[owned], jj[owned]))
        return cls(plate_id, frame, crust_type, n, np.concatenate(kept), **kwargs)

    # --- Storage-level accessors -------------------------------------------------------

    @property
    def cells_per_edge(self) -> int:
        return self._n

    @property
    def cell_keys(self) -> np.ndarray:
        """Packed keys of every active cell, in canonical node order (read-only view)."""
        view = self._keys.view()
        view.flags.writeable = False
        return view

    @staticmethod
    def _check_field(name: str, values: np.ndarray, expected: int) -> None:
        if name not in SURFACE_FIELDS:
            raise ValueError(f"unknown surface field(s): {name}")
        if np.asarray(values).shape != (expected,):
            raise ValueError(f"surface fields must have shape ({expected},); got {{{name!r}: {np.asarray(values).shape}}}")

    def _field_storage(self, name: str) -> np.ndarray:
        """The backing array for `name`, materialised at its registry default on first use --
        fields nothing has written yet cost no memory."""
        values = self._fields.get(name)
        if values is None:
            spec = SURFACE_FIELDS[name]
            values = np.full(len(self._keys), spec.default, dtype=spec.dtype)
            self._fields[name] = values
        return values

    # --- Caches ------------------------------------------------------------------------

    def _reset_caches(self) -> None:
        for name in self._TOPOLOGY_CACHES:
            setattr(self, name, None)
        self._invalidate_bounding_polygon()

    def _invalidate_bounding_polygon(self) -> None:
        super()._invalidate_bounding_polygon()
        self._world_points_cache = None

    def _unpacked(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return unpack_cell_keys(self._keys)

    def _local_centres(self) -> np.ndarray:
        if self._local_cache is None:
            face, level, i, j = self._unpacked()
            if len(face):
                scale = np.left_shift(1, level)
                self._local_cache = lattice_points(face, (i + 0.5) / scale, (j + 0.5) / scale, self._n)
            else:
                self._local_cache = np.zeros((0, 3))
        return self._local_cache

    def _node_latlon(self) -> tuple[np.ndarray, np.ndarray]:
        """Plate-local (phi, theta) of every node -- only for the `ElevationPoint` view,
        whose protocol still speaks in those terms."""
        if self._latlon_cache is None:
            self._latlon_cache = geometry.xyz_to_latlon(self._local_centres())
        return self._latlon_cache

    def _get_world_points(self) -> np.ndarray:
        if self._world_points_cache is None:
            self._world_points_cache = geometry.to_world(self._frame, self._local_centres())
        return self._world_points_cache

    def node_areas_m2(self) -> np.ndarray:
        if self._area_cache is None:
            _, level, i, j = self._unpacked()
            self._area_cache = np.array(
                [cell_areas_sr(np.array([ii]), np.array([jj]), self._n * (1 << int(ll)))[0] for ll, ii, jj in zip(level, i, j)]
            ) * PLANET_RADIUS_M**2
        return self._area_cache

    def row_intervals(self) -> CellIntervals:
        """Runs of active cells along each face row (constant `j`)."""
        if self._row_intervals_cache is None:
            face, level, i, j = self._unpacked()
            self._row_intervals_cache = _intervals(face, j, i, level)
        return self._row_intervals_cache

    def column_intervals(self) -> CellIntervals:
        """Runs of active cells along each face column (constant `i`)."""
        if self._column_intervals_cache is None:
            face, level, i, j = self._unpacked()
            self._column_intervals_cache = _intervals(face, i, j, level)
        return self._column_intervals_cache

    def _index_of_keys(self, keys: np.ndarray) -> np.ndarray:
        """Node index of each key, or -1 where the cell isn't active."""
        keys = np.asarray(keys, dtype=np.int64)
        if len(self._keys) == 0:
            return np.full(keys.shape, -1, dtype=np.int64)
        pos = np.clip(np.searchsorted(self._keys, keys), 0, len(self._keys) - 1)
        return np.where(self._keys[pos] == keys, pos, -1)

    def _leaf_key_at(self, points: np.ndarray) -> np.ndarray:
        """Active leaf containing each local point, or -1 outside this patch."""
        result = np.full(len(points), -1, dtype=np.int64)
        max_level = int(unpack_cell_keys(self._keys)[1].max(initial=0))
        for level in range(max_level + 1):
            face, i, j = locate_cells(points, self._n * (1 << level))
            keys = pack_cell_keys(face, i, j, level)
            found = self._index_of_keys(keys) >= 0
            result[found] = keys[found]
        return result

    def _neighbour_indices(self) -> np.ndarray:
        """(n, 4) node index of each edge neighbour, -1 where that side is exposed."""
        if self._adjacency_cache is None:
            face, level, i, j = self._unpacked()
            rows: list[list[int]] = []
            for f, lev, ii, jj in zip(face, level, i, j):
                scale = 1 << int(lev)
                neighbours: list[int] = []
                for direction in range(4):
                    for along in (0.25, 0.75):
                        if direction == 0:
                            a, b = ii + along, jj - _NEIGHBOUR_PROBE_FRACTION
                        elif direction == 1:
                            a, b = ii + 1 + _NEIGHBOUR_PROBE_FRACTION, jj + along
                        elif direction == 2:
                            a, b = ii + 1 - along, jj + 1 + _NEIGHBOUR_PROBE_FRACTION
                        else:
                            a, b = ii - _NEIGHBOUR_PROBE_FRACTION, jj + 1 - along
                        point = lattice_points(np.array([f]), np.array([a / scale]), np.array([b / scale]), self._n)
                        key = int(self._leaf_key_at(point)[0])
                        index = int(self._index_of_keys(np.array([key]))[0]) if key >= 0 else -1
                        if index >= 0 and index not in neighbours:
                            neighbours.append(index)
                rows.append(sorted(neighbours))
            width = max((len(row) for row in rows), default=0)
            self._adjacency_cache = np.full((len(rows), width), -1, dtype=np.int64)
            for row_index, row in enumerate(rows):
                self._adjacency_cache[row_index, : len(row)] = row
        return self._adjacency_cache

    # --- PlateSurface ------------------------------------------------------------------

    def node_count(self) -> int:
        return len(self._keys)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        return self._get_world_points(), self.collect("elevation")

    def surface_nodes(self, *field_names: str) -> SurfaceNodes:
        return SurfaceNodes(
            local_xyz=self._local_centres(),
            world_xyz=self._get_world_points(),
            node_ids=np.column_stack([self._keys.view(np.uint64), np.zeros(len(self._keys), dtype=np.uint64)]),
            area_m2=self.node_areas_m2(),
            area_is_exact=True,
            fields={name: self.collect(name) for name in field_names},
        )

    def node_index_for_id(self, node_id: np.ndarray | tuple[int, int]) -> int | None:
        target = np.asarray(node_id, dtype=np.uint64)
        if target.shape != (2,):
            raise ValueError("surface node ID must contain exactly two uint64 words")
        if target[1] != 0:
            return None
        index = int(self._index_of_keys(target[:1].view(np.int64))[0])
        return None if index < 0 else index

    def collect(self, field_name: str) -> np.ndarray:
        if field_name not in SURFACE_FIELDS:
            raise AttributeError(field_name)
        return self._field_storage(field_name).copy()

    def set_fields_on_plate(self, **fields: np.ndarray) -> None:
        expected = len(self._keys)
        invalid = sorted(name for name in fields if name not in SURFACE_FIELDS)
        if invalid:
            raise ValueError(f"unknown surface field(s): {', '.join(invalid)}")
        wrong = {name: np.asarray(values).shape for name, values in fields.items() if np.asarray(values).shape != (expected,)}
        if wrong:
            raise ValueError(f"surface fields must have shape ({expected},); got {wrong}")
        for name, values in fields.items():
            self._field_storage(name)[:] = values

    def adjacency(self) -> SurfaceAdjacency:
        """Edge-sharing (4-connected) cell adjacency, from the cells themselves rather than
        inferred from node spacing. Symmetric: cells match edge for edge across seams."""
        neighbours = self._neighbour_indices()
        valid = neighbours >= 0
        offsets = np.zeros(len(neighbours) + 1, dtype=np.int64)
        offsets[1:] = np.cumsum(valid.sum(axis=1))
        flat = np.sort(np.where(valid, neighbours, np.iinfo(np.int64).max), axis=1)
        return SurfaceAdjacency(offsets, flat[flat != np.iinfo(np.int64).max])

    def _local_boundary_loops(self) -> list[np.ndarray]:
        """Every closed loop of exposed cell edges, as local corner-point arrays -- outer
        boundaries counter-clockwise, holes clockwise. Where two loops touch at a single
        corner (diagonal-only contact), the traversal takes the left-most turn, which keeps
        each loop simple instead of figure-eighting through the pinch."""
        if self._local_loops_cache is not None:
            return self._local_loops_cache
        face, level, i, j = self._unpacked()
        max_level = int(level.max(initial=0))
        max_n = self._n * (1 << max_level)
        # Express every leaf edge as finest-level atomic segments. Shared edges then cancel
        # exactly, including a coarse edge facing two fine cells and cube-face seams.
        segments: dict[tuple[int, int], tuple[int, int, np.ndarray, np.ndarray]] = {}
        for f, lev, ii, jj in zip(face, level, i, j):
            scale = 1 << (max_level - int(lev))
            x0, y0 = int(ii) * scale, int(jj) * scale
            oriented = []
            for k in range(scale):
                oriented.extend(
                    [
                        ((x0 + k, y0), (x0 + k + 1, y0)),
                        ((x0 + scale, y0 + k), (x0 + scale, y0 + k + 1)),
                        ((x0 + scale - k, y0 + scale), (x0 + scale - k - 1, y0 + scale)),
                        ((x0, y0 + scale - k), (x0, y0 + scale - k - 1)),
                    ]
                )
            for (ai, aj), (bi, bj) in oriented:
                akey = int(_corner_keys(np.array([f]), np.array([ai]), np.array([aj]), max_n)[0])
                bkey = int(_corner_keys(np.array([f]), np.array([bi]), np.array([bj]), max_n)[0])
                identity = tuple(sorted((akey, bkey)))
                if identity in segments:
                    del segments[identity]
                else:
                    axyz = lattice_points(np.array([f]), np.array([ai]), np.array([aj]), max_n)[0]
                    bxyz = lattice_points(np.array([f]), np.array([bi]), np.array([bj]), max_n)[0]
                    segments[identity] = (akey, bkey, axyz, bxyz)
        if not segments:
            self._local_loops_cache = []
            return self._local_loops_cache
        start_key = np.array([segment[0] for segment in segments.values()])
        end_key = np.array([segment[1] for segment in segments.values()])
        start_xyz = np.array([segment[2] for segment in segments.values()])
        end_xyz = np.array([segment[3] for segment in segments.values()])

        by_start: dict[int, list[int]] = {}
        for e, key in enumerate(start_key.tolist()):
            by_start.setdefault(key, []).append(e)

        def successor(e: int) -> int:
            candidates = by_start[int(end_key[e])]
            if len(candidates) == 1:
                return candidates[0]
            vertex = end_xyz[e]
            incoming = vertex - start_xyz[e]
            return max(candidates, key=lambda c: float(np.dot(np.cross(incoming, end_xyz[c] - vertex), vertex)))

        used = np.zeros(len(segments), dtype=bool)
        loops: list[np.ndarray] = []
        for first in range(len(segments)):
            if used[first]:
                continue
            chain = []
            e = first
            while not used[e]:
                used[e] = True
                chain.append(e)
                e = successor(e)
            loops.append(start_xyz[chain])
        self._local_loops_cache = loops
        return loops

    def boundary_loops_world(self) -> tuple[np.ndarray, ...]:
        return tuple(geometry.to_world(self._frame, loop) for loop in self._local_boundary_loops())

    def outline_world(self) -> np.ndarray:
        """Every boundary loop joined into one keyholed vertex loop, the single-array form
        `get_bounding_polygon()` consumers expect (see `plates._stitch_loops`)."""
        loops = self._local_boundary_loops()
        if not loops:
            return np.zeros((0, 3))
        return geometry.to_world(self._frame, _stitch_xyz_loops(loops))

    def contains_batch(self, points_xyz: np.ndarray) -> np.ndarray:
        """Exact: a point is inside iff the lattice cell containing it is active. No
        polygon test and no boundary band -- the cells are the territory."""
        points_xyz = np.asarray(points_xyz, dtype=float).reshape(-1, 3)
        if len(points_xyz) == 0 or len(self._keys) == 0:
            return np.zeros(len(points_xyz), dtype=bool)
        return self._leaf_key_at(geometry.to_local(self._frame, points_xyz)) >= 0

    def contains(self, lat: float, lon: float) -> bool:
        return bool(self.contains_batch(geometry.latlon_to_xyz(np.asarray(lat), np.asarray(lon)))[0])

    def get_neighbours(self, all_plates: list[Plate], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list[Plate]:
        return _plates_within(self, all_plates, threshold_rad)

    # --- Per-node iteration --------------------------------------------------------------

    def __iter__(self) -> Iterator[ElevationPoint]:
        for index in range(len(self._keys)):
            yield ElevationPointInPatch(self, index)

    def map_world_points(self) -> Iterator[tuple[ElevationPoint, np.ndarray]]:
        world = self._get_world_points()
        for index in range(len(self._keys)):
            yield ElevationPointInPatch(self, index), world[index]

    def map_world_points_on_plate(self) -> Iterator[tuple[ElevationPoint, np.ndarray, float]]:
        """The "how far across the plate" fraction is measured along the node's own row run
        (0 and 1 at the run's end cells), the quad analogue of `PlateWithLines` measuring
        along its line."""
        world = self._get_world_points()
        rows = self.row_intervals()
        _, _, i, _ = self._unpacked()
        start = rows.start[rows.node_run]
        span = rows.end[rows.node_run] - start
        fraction = np.where(span == 0, 0.5, (i - start) / np.maximum(span, 1))
        for index in range(len(self._keys)):
            yield ElevationPointInPatch(self, index), world[index], float(fraction[index])

    # --- Adaptive remeshing ---------------------------------------------------------------

    def _replace_topology(self, keys: np.ndarray, fields: Mapping[str, np.ndarray]) -> None:
        order = np.argsort(keys, kind="stable")
        self._keys = np.asarray(keys, dtype=np.int64)[order]
        self._fields = {
            name: np.asarray(values, dtype=SURFACE_FIELDS[name].dtype)[order]
            for name, values in fields.items()
        }
        self._topology_revision += 1
        self._geometry_revision += 1
        self._reset_caches()

    def refine_cells(self, cell_ids: np.ndarray) -> dict[int, tuple[int, ...]]:
        """Subdivide selected leaves and any coarser edge neighbours needed for 2:1 balance.

        Returns explicit old-ID -> child-ID lineage. Unselected cell IDs remain stable.
        """
        requested = {int(key) for key in np.asarray(cell_ids, dtype=np.int64).reshape(-1)}
        active = set(map(int, self._keys))
        missing = requested - active
        if missing:
            raise ValueError(f"can only refine active leaf cells: {sorted(missing)}")
        max_level = int(np.floor(np.log2(_INDEX_MASK / self._n)))
        requested_levels = unpack_cell_keys(np.asarray(sorted(requested), dtype=np.int64))[1]
        if np.any(requested_levels >= max_level):
            raise ValueError(f"cannot refine past level {max_level}")
        selected = set(requested)
        changed = True
        while changed:
            changed = False
            graph = self.adjacency()
            _, levels, _, _ = unpack_cell_keys(self._keys)
            for index, key in enumerate(self._keys):
                if int(key) not in selected:
                    continue
                for neighbour in graph.neighbours[graph.offsets[index] : graph.offsets[index + 1]]:
                    if levels[neighbour] < levels[index] and int(self._keys[neighbour]) not in selected:
                        selected.add(int(self._keys[neighbour]))
                        changed = True
        lineage = {key: tuple(map(int, child_cell_keys(np.array([key]))[0])) for key in sorted(selected)}
        new_keys: list[int] = []
        sources: list[int] = []
        for index, key in enumerate(self._keys):
            children = lineage.get(int(key))
            if children is None:
                new_keys.append(int(key))
                sources.append(index)
            else:
                new_keys.extend(children)
                sources.extend([index] * 4)
        fields = {name: values[np.asarray(sources)] for name, values in self._fields.items()}
        self._replace_topology(np.asarray(new_keys), fields)
        return lineage

    def _coarsened_value(
        self,
        name: str,
        values: np.ndarray,
        areas: np.ndarray,
        all_fields: Mapping[str, np.ndarray],
    ) -> float | int | bool:
        spec = SURFACE_FIELDS[name]
        if name == "elevation" and all(
            field in all_fields
            for field in ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "crust_type_code")
        ):
            from . import lithosphere

            hc = all_fields["crustal_thickness_m"]
            hm = all_fields["mantle_lithosphere_thickness_m"]
            codes = all_fields["crust_type_code"]
            density = lithosphere.node_crust_density(codes, self.crust_type)
            residual = values - lithosphere.isostatic_elevation(hc, hm, density)
            new_hc = float(np.average(hc, weights=areas))
            new_hm = float(np.average(hm, weights=areas))
            new_code = self._coarsened_value("crust_type_code", codes, areas, all_fields)
            new_density = lithosphere.node_crust_density(np.array([new_code]), self.crust_type)
            return float(
                lithosphere.isostatic_elevation(np.array([new_hc]), np.array([new_hm]), new_density)[0]
                + np.average(residual, weights=areas)
            )
        if name == "channel_depth":
            return float(np.max(values))
        if name == "channel_width":
            return float(values[np.argmax(all_fields.get("channel_depth", np.zeros(len(values))))])
        weights = areas
        if name in ("soil_mineral_content", "soil_organic_content"):
            weights = weights * all_fields.get("soil_depth", np.ones(len(values)))
        if spec.remap_class in (RemapClass.EXTENSIVE, RemapClass.INTENSIVE, RemapClass.CLOCK):
            return float(np.average(values, weights=weights)) if weights.sum() else float(np.mean(values))
        if spec.remap_class == RemapClass.CATEGORICAL:
            if name == "crust_type_code":
                inherited = CRUST_TYPE_CONTINENTAL if self.crust_type == "continental" else CRUST_TYPE_OCEANIC
                values = np.where(values == 0, inherited, values)
            choices = np.unique(values)
            totals = np.array([areas[values == choice].sum() for choice in choices])
            tied = choices[totals == totals.max()]
            if name == "elev_change_reason":
                structural = tied[
                    ((tied >= 1) & (tied <= 8)) | ((tied >= 15) & (tied <= 17)) | (tied == 19)
                ]
                if len(structural):
                    tied = structural
            return tied[0]
        if spec.remap_class == RemapClass.BOOLEAN_PROVENANCE:
            return bool(np.any(values))
        if spec.remap_class == RemapClass.COUNTDOWN:
            return float(np.max(values))
        if spec.remap_class in (RemapClass.HISTORY, RemapClass.WRITE_ONCE_HISTORY):
            valid = values != spec.sentinel
            return spec.sentinel if not np.any(valid) else float(np.min(values[valid]))
        return float(np.average(values, weights=areas))

    def coarsen_cells(self, parent_ids: np.ndarray) -> dict[int, int]:
        """Replace complete sibling quartets by their parent, conservatively by exact area.

        A coarsen that would break the mesh's 2:1 level balance is rejected.
        Returns child-ID -> parent-ID lineage.
        """
        parents = {int(key) for key in np.asarray(parent_ids, dtype=np.int64).reshape(-1)}
        active = set(map(int, self._keys))
        groups = {parent: tuple(map(int, child_cell_keys(np.array([parent]))[0])) for parent in parents}
        incomplete = [parent for parent, children in groups.items() if not set(children) <= active]
        if incomplete:
            raise ValueError(f"coarsening requires four active children: {sorted(incomplete)}")
        graph = self.adjacency()
        _, levels, _, _ = unpack_cell_keys(self._keys)
        index = {int(key): k for k, key in enumerate(self._keys)}
        removing = {child for children in groups.values() for child in children}
        for parent, children in groups.items():
            parent_level = int(unpack_cell_keys(np.array([parent]))[1][0])
            for child in children:
                ci = index[child]
                for neighbour in graph.neighbours[graph.offsets[ci] : graph.offsets[ci + 1]]:
                    if int(self._keys[neighbour]) not in removing and int(levels[neighbour]) > parent_level + 1:
                        raise ValueError("coarsening would violate 2:1 balance")
        old_areas = self.node_areas_m2()
        new_keys = [int(key) for key in self._keys if int(key) not in removing] + sorted(parents)
        output = {name: [] for name in self._fields}
        retained = [k for k, key in enumerate(self._keys) if int(key) not in removing]
        for name, values in self._fields.items():
            output[name].extend(values[retained].tolist())
            for parent in sorted(parents):
                child_indices = np.array([index[child] for child in groups[parent]])
                child_fields = {field: stored[child_indices] for field, stored in self._fields.items()}
                output[name].append(self._coarsened_value(name, values[child_indices], old_areas[child_indices], child_fields))
        self._replace_topology(np.asarray(new_keys), {name: np.asarray(values) for name, values in output.items()})
        return {child: parent for parent, children in groups.items() for child in children}

    # --- Topology changes (later phases) ---------------------------------------------------

    def _plates_from_node_masks(self, masks: list[np.ndarray], ids: list[int]) -> list[Plate]:
        raise NotImplementedError("sparse quad patches are static until issue #228 Phase 4")

    # --- Persistence -----------------------------------------------------------------------

    def __getstate__(self) -> dict:
        state = {
            name: value
            for name, value in self.__dict__.items()
            if name not in self._TOPOLOGY_CACHES and name not in self._GEOMETRY_CACHES
        }
        state["_surface_format_version"] = QUAD_SURFACE_FORMAT_VERSION
        return state

    def __setstate__(self, state: dict) -> None:
        version = state.get("_surface_format_version")
        if version != QUAD_SURFACE_FORMAT_VERSION:
            raise ValueError(
                f"unsupported sparse quad surface format version {version!r}; "
                f"this build reads version {QUAD_SURFACE_FORMAT_VERSION}"
            )
        state = dict(state)
        del state["_surface_format_version"]
        self.__dict__.update(state)
        self._reset_caches()
