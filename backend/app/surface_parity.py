"""Plate-surface audit runs for issue #228 Phase 5 (issues #247/#249).

Steps one world per seed and records, at agreed ages, the measurements `surface_parity_gates`
judges. These runs began as line-vs-quad parity runs; since the line surface was retired
(#251), each run audits the quad surface on its own.

A run writes two documents:

- `seed<S>.json` -- deterministic metrics only (sorted keys, floats rounded), so a rerun of the
  same code on the same inputs reproduces it byte for byte.
- `seed<S>.timings.json` -- wall-clock timings: per-step phase buckets, derived-index build
  counts/times, index query timings and memory. Never byte-stable.

What a run records:

- **Audits** (every `audit_every` steps and at each checkpoint) -- hard invariants, each
  violation recorded with its step: finite, correctly sized, bounded fields; orthonormal
  frames; quad leaf topology (no overlapping leaves, 2:1 balance), neighbour validity
  (in range, symmetric, sharing an edge), folded cells; every *populated* derived cache on a
  plate or on `World` equal to a fresh rebuild; and revision tracking -- a plate whose local
  node set or frame changed must have bumped its topology or geometry revision, respectively.
  The atmospheric HEALPix grid must keep its resolution, independent of node count.
- **Checkpoints** -- area-weighted conservation totals (exact cell areas from
  `SurfaceNodes.area_m2`), land fraction, elevation percentiles, plate counts, whole-sphere
  coverage/overlap, neutral sample-cloud quality, quad lattice quality, climate and hydrology
  summaries, and derived-index parity: the cached world KD-tree against a fresh one (exact),
  and HEALPix nearest-value resampling against the KD-tree (statistical, broken out near
  coasts, plate boundaries, holes, poles and the antimeridian).
- **Load checks** at each checkpoint -- a save/load round trip must reproduce the
  authoritative state exactly and come back with every derived index dropped; at the final
  checkpoint, one step from the loaded copy must match one step from the in-memory world.

The rules that turn these into pass/warn/fail live in `surface_parity_gates`; see
docs/surface-parity.md for the command line and the tolerance table.
"""

from __future__ import annotations

import copy
import dataclasses
import functools
import hashlib
import inspect
import json
import time
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
from scipy.spatial import cKDTree

from . import (
    elevation_lines,
    erosion,
    eustasy,
    faults,
    gaps,
    geology,
    healpix_grid,
    hydrology,
    lithosphere,
    magma_transport,
    merge_split,
    persistence,
    plates as plates_mod,
    stranded_basins,
    volcanism,
    world as world_mod,
)
from .elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from .plates import Plate
from .sparse_quad_patch import (
    PlateWithSparseQuadPatch,
    _corner_coordinates,
    lattice_points,
    unpack_cell_keys,
)
from .surface_fields import SURFACE_FIELDS

# 2: runs no longer carry a `surface` (#251); every run is a quad run.
RESULT_SCHEMA_VERSION = 2

# Distances below are in multiples of the world's lattice spacing `s`, as in
# docs/plate-surface-baseline.md.
STACKED_SPACING = 0.5
VOID_SPACING = 1.5
BOUNDARY_SPACING = 1.5
CATEGORY_SPACING = 2.0
ANISOTROPY_NEIGHBOURS = 8
ANISOTROPY_CAP = 100.0
ANISOTROPIC_RATIO = 2.0
BOUNDARY_NEIGHBOURS = 12
BOUNDARY_SECTOR_RAD = np.deg2rad(135.0)
ASPECT_LIMIT = 4.0
SKEW_LIMIT_DEG = 45.0
POLE_LAT_DEG = 80.0
ANTIMERIDIAN_LON_DEG = 177.0
# Climate/hydrology series kept per step for the stability gate.
SERIES_STATS = ("sea_level_m", "land_fraction", "air_temperature_mean_c", "precipitation_mean_mm", "ocean_temperature_mean_c")
CHECKPOINT_STATS = SERIES_STATS + (
    "land_temperature_mean_c",
    "elevation_mean_m",
    "elevation_std_m",
    "ocean_depth_mean_m",
    "hc_mean_m",
    "hm_mean_m",
)


@dataclass(frozen=True)
class RunConfig:
    """Everything that determines a run's metrics. `checkpoints_myr` must be multiples of
    `step_years`; age 0 (the freshly generated world) is always a checkpoint."""

    name: str
    node_density: float
    step_years: float
    checkpoints_myr: tuple[float, ...]
    audit_every: int
    samples: int
    climate_density: float | None = None
    fluid_density: float = 1.0
    load_checks: bool = True
    render: bool = False
    render_size: tuple[int, int] = (720, 400)

    @property
    def steps(self) -> int:
        return self.checkpoint_steps[-1]

    @property
    def checkpoint_steps(self) -> tuple[int, ...]:
        steps = []
        for age in sorted(set((0.0,) + tuple(self.checkpoints_myr))):
            exact = age * 1e6 / self.step_years
            step = int(round(exact))
            if abs(step - exact) > 1e-6:
                raise ValueError(f"checkpoint {age} Myr is not a multiple of the {self.step_years:g}-year step")
            steps.append(step)
        return tuple(steps)

    def to_json(self) -> dict:
        data = dataclasses.asdict(self)
        data["checkpoints_myr"] = list(self.checkpoints_myr)
        data["render_size"] = list(self.render_size)
        return data


# `issue147` is the issue #147 profile world (seed 0 at density 4, climate 4, fluid 2, 100 kyr
# steps; analysis/issue147-profile-20260922), so its timings line up with that profile.
PRESETS: dict[str, RunConfig] = {
    "smoke": RunConfig("smoke", node_density=0.5, step_years=1e6, checkpoints_myr=(2.0, 4.0), audit_every=1, samples=20_000),
    "standard": RunConfig(
        "standard", node_density=1.0, step_years=1e6, checkpoints_myr=(30.0, 60.0, 120.0), audit_every=5, samples=200_000
    ),
    "long": RunConfig(
        "long", node_density=1.0, step_years=1e6, checkpoints_myr=(30.0, 60.0, 120.0, 240.0, 400.0), audit_every=10, samples=200_000
    ),
    "issue147": RunConfig(
        "issue147",
        node_density=4.0,
        step_years=1e5,
        checkpoints_myr=(10.0, 30.0, 60.0),
        audit_every=50,
        samples=200_000,
        climate_density=4.0,
        fluid_density=2.0,
    ),
}


# --- Small helpers --------------------------------------------------------------------------


def sig(value, digits: int = 6):
    """Round floats (recursively) to `digits` significant figures so reruns diff cleanly."""
    if isinstance(value, dict):
        return {str(k): sig(v, digits) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [sig(v, digits) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        value = float(value)
        if not np.isfinite(value):
            return repr(value)
        return float(f"{value:.{digits}g}")
    return value


def dumps(document: dict) -> str:
    return json.dumps(sig(document), indent=2, sort_keys=True) + "\n"


def fibonacci_sphere(n: int) -> np.ndarray:
    k = np.arange(n) + 0.5
    z = 1.0 - 2.0 * k / n
    r = np.sqrt(1.0 - z * z)
    lon = np.pi * (1.0 + 5**0.5) * k
    return np.column_stack([r * np.cos(lon), r * np.sin(lon), z])


def fraction(mask: np.ndarray) -> float | None:
    return float(np.mean(mask)) if len(mask) else None


def weighted_percentiles(values: np.ndarray, weights: np.ndarray, qs=(5, 50, 95)) -> dict:
    if len(values) == 0:
        return {f"p{q:02d}": None for q in qs}
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order])
    cumulative /= cumulative[-1]
    return {f"p{q:02d}": float(values[order][min(np.searchsorted(cumulative, q / 100.0), len(values) - 1)]) for q in qs}


def live_plates(world) -> list[Plate]:
    return [p for p in world.plates if p.node_count() > 0]


def state_hash(world) -> str:
    """SHA-256 of the world's authoritative plate state (frames, rotation, local node set and
    every surface field) plus sea level and age. Derived caches never contribute."""
    digest = hashlib.sha256()
    digest.update(np.float64([world.elapsed_years, world.sea_level_m]).tobytes())
    names = tuple(sorted(SURFACE_FIELDS))
    for plate in sorted(world.plates, key=lambda p: p.plate_id):
        digest.update(np.int64([plate.plate_id, plate.node_count()]).tobytes())
        digest.update(plate.crust_type.encode())
        digest.update(np.ascontiguousarray(plate.frame, dtype=np.float64).tobytes())
        digest.update(np.ascontiguousarray(plate.omega, dtype=np.float64).tobytes())
        if plate.node_count() == 0:
            continue
        nodes = plate.surface_nodes(*names)
        digest.update(np.ascontiguousarray(nodes.local_xyz).tobytes())
        for name in names:
            digest.update(np.ascontiguousarray(nodes.fields[name]).tobytes())
    return digest.hexdigest()


# --- Instrumentation ------------------------------------------------------------------------

# (owner, attribute, bucket). Top-level step phases only; nested timed calls are charged to
# the outermost bucket so buckets never double count. `deform` and `topology`, together with
# `gap_fill` and `overlap_tracking`, are the deformation/topology phases #228 calls out.
_PHASES = (
    (Plate, "shift", "shift"),
    (PlateWithSparseQuadPatch, "deform", "deform"),
    (faults, "update_faults", "faults"),
    (merge_split, "apply_topology_changes", "topology"),
    (faults, "reconcile_faults", "topology"),
    (merge_split, "update_overlap_tracking", "overlap_tracking"),
    (gaps, "fill_gaps_by_growing_neighbours", "gap_fill"),
    (gaps, "reconcile_gap_tracks", "gap_fill"),
    (magma_transport, "run_magma_transport", "magma_transport"),
    (world_mod, "_advance_fluid_dynamics", "fluid_dynamics"),
    (erosion, "apply_erosion", "climate_erosion_hydrology"),
    (volcanism, "apply_volcanic_activity", "volcanism"),
    (geology, "apply_resource_formation", "resource_formation"),
    (stranded_basins, "reconcile_world_tracks", "stranded_basins"),
    (eustasy, "update_sea_level", "sea_level"),
    (world_mod.World, "record_stats", "record_stats"),
)
DEFORM_TOPOLOGY_PHASES = ("deform", "topology", "gap_fill", "overlap_tracking")

# Derived-index builds: (owner, method, cache attribute checked before the call, index name).
_PLATE_INDEXES = (
    (Plate, "get_node_kdtree", "_node_kdtree_cache", "plate_node_kdtree"),
    (Plate, "get_bounding_polygon", "_bounding_polygon_cache", "plate_outline"),
    (Plate, "get_bounding_polygon_tree", "_bounding_polygon_tree_cache", "plate_outline_kdtree"),
    (PlateWithSparseQuadPatch, "_neighbour_indices", "_adjacency_cache", "quad_adjacency"),
    (PlateWithSparseQuadPatch, "_local_boundary_loops", "_local_loops_cache", "quad_boundary_loops"),
)
# Module-level world caches: modules that import the helper by name are patched too.
_WORLD_INDEXES = (
    ("cached_node_position_tree", "node_position_tree_cache", "world_node_kdtree", (plates_mod, hydrology)),
    ("cached_node_healpix_index", "node_healpix_index_cache", "world_node_healpix", (plates_mod,)),
)


class Instrumentation:
    """Per-step phase timers and derived-index build counters, installed by monkeypatching
    for the duration of `installed()` and removed afterwards."""

    def __init__(self) -> None:
        self.phase_seconds: dict[str, float] = defaultdict(float)
        self.index_builds: dict[str, int] = defaultdict(int)
        self.index_build_seconds: dict[str, float] = defaultdict(float)
        self.index_hits: dict[str, int] = defaultdict(int)
        self._depth = 0

    def reset(self) -> None:
        self.phase_seconds.clear()
        self.index_builds.clear()
        self.index_build_seconds.clear()
        self.index_hits.clear()

    def snapshot(self) -> dict:
        return {
            "phases": dict(self.phase_seconds),
            "index_builds": dict(self.index_builds),
            "index_build_seconds": dict(self.index_build_seconds),
            "index_hits": dict(self.index_hits),
        }

    def _timed(self, original, bucket: str):
        instrumentation = self

        if inspect.isgeneratorfunction(original):

            @functools.wraps(original)
            def generator(*args, **kwargs):
                iterator = original(*args, **kwargs)
                while True:
                    with instrumentation._phase(bucket):
                        try:
                            item = next(iterator)
                        except StopIteration as stop:
                            return stop.value
                    yield item

            return generator

        @functools.wraps(original)
        def timed(*args, **kwargs):
            with instrumentation._phase(bucket):
                return original(*args, **kwargs)

        return timed

    @contextmanager
    def _phase(self, bucket: str):
        self._depth += 1
        started = time.perf_counter()
        try:
            yield
        finally:
            self._depth -= 1
            if self._depth == 0:
                self.phase_seconds[bucket] += time.perf_counter() - started

    def _counted_method(self, original, cache_attr: str, name: str):
        instrumentation = self

        @functools.wraps(original)
        def counted(plate, *args, **kwargs):
            if getattr(plate, cache_attr, None) is not None:
                instrumentation.index_hits[name] += 1
                return original(plate, *args, **kwargs)
            started = time.perf_counter()
            try:
                return original(plate, *args, **kwargs)
            finally:
                instrumentation.index_builds[name] += 1
                instrumentation.index_build_seconds[name] += time.perf_counter() - started

        return counted

    def _counted_world(self, original, cache_attr: str, name: str):
        instrumentation = self

        @functools.wraps(original)
        def counted(world, points, *args, **kwargs):
            if world is not None and getattr(world, cache_attr, None) is not None:
                instrumentation.index_hits[name] += 1
                return original(world, points, *args, **kwargs)
            started = time.perf_counter()
            try:
                return original(world, points, *args, **kwargs)
            finally:
                instrumentation.index_builds[name] += 1
                instrumentation.index_build_seconds[name] += time.perf_counter() - started

        return counted

    @contextmanager
    def installed(self) -> Iterator["Instrumentation"]:
        patches: list[tuple[object, str, object]] = []

        def patch(owner, attribute: str, replacement) -> None:
            patches.append((owner, attribute, owner.__dict__[attribute] if isinstance(owner, type) else getattr(owner, attribute)))
            setattr(owner, attribute, replacement)

        try:
            for owner, attribute, bucket in _PHASES:
                patch(owner, attribute, self._timed(getattr(owner, attribute), bucket))
            for owner, method, cache_attr, name in _PLATE_INDEXES:
                patch(owner, method, self._counted_method(getattr(owner, method), cache_attr, name))
            for function, cache_attr, name, modules in _WORLD_INDEXES:
                counted = self._counted_world(getattr(plates_mod, function), cache_attr, name)
                for module in modules:
                    patch(module, function, counted)
            yield self
        finally:
            for owner, attribute, original in reversed(patches):
                setattr(owner, attribute, original)


# --- Audits ---------------------------------------------------------------------------------

# Cache attribute -> the method that (re)populates it, per representation.
_BASE_CACHE_BUILDERS = {
    "_bounding_polygon_cache": "get_bounding_polygon",
    "_bounding_polygon_tree_cache": "get_bounding_polygon_tree",
    "_node_kdtree_cache": "get_node_kdtree",
}
_QUAD_CACHE_BUILDERS = {
    "_local_cache": "_local_centres",
    "_latlon_cache": "_node_latlon",
    "_area_cache": "node_areas_m2",
    "_adjacency_cache": "_neighbour_indices",
    "_local_loops_cache": "_local_boundary_loops",
    "_row_intervals_cache": "row_intervals",
    "_column_intervals_cache": "column_intervals",
    "_probe_neighbours_cache": "_probe_neighbour_indices",
    "_world_points_cache": "_get_world_points",
}


def cache_builders(plate: Plate) -> dict[str, str]:
    builders = dict(_BASE_CACHE_BUILDERS)
    if isinstance(plate, PlateWithSparseQuadPatch):
        builders.update(_QUAD_CACHE_BUILDERS)
    return builders


def derived_equal(a, b) -> bool:
    """Structural equality for cached derived values: arrays, sequences, dataclasses, KD-trees
    (by indexed data), and plain objects (by attributes)."""
    if a is None or b is None:
        return a is b
    if isinstance(a, cKDTree):
        return isinstance(b, cKDTree) and np.array_equal(a.data, b.data)
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        a, b = np.asarray(a), np.asarray(b)
        return a.shape == b.shape and np.array_equal(a, b, equal_nan=a.dtype.kind == "f")
    if isinstance(a, (list, tuple)):
        return isinstance(b, (list, tuple)) and len(a) == len(b) and all(derived_equal(x, y) for x, y in zip(a, b))
    if dataclasses.is_dataclass(a):
        return type(a) is type(b) and all(derived_equal(getattr(a, f.name), getattr(b, f.name)) for f in dataclasses.fields(a))
    if hasattr(a, "__dict__"):
        return type(a) is type(b) and vars(a).keys() == vars(b).keys() and all(derived_equal(vars(a)[k], vars(b)[k]) for k in vars(a))
    return a == b


def fresh_copy(plate: Plate) -> Plate:
    """A shallow copy of `plate` with every derived cache dropped -- rebuilding a cache on it
    reads the same authoritative state without disturbing the original's caches."""
    fresh = copy.copy(plate)
    for attribute in cache_builders(plate):
        setattr(fresh, attribute, None)
    return fresh


def stale_plate_caches(plate: Plate) -> list[str]:
    """Names of this plate's populated derived caches that differ from a fresh rebuild."""
    fresh = fresh_copy(plate)
    stale = []
    for attribute, builder in cache_builders(plate).items():
        cached = getattr(plate, attribute, None)
        if cached is None:
            continue
        getattr(fresh, builder)()
        if not derived_equal(cached, getattr(fresh, attribute)):
            stale.append(attribute)
    return stale


def stale_world_caches(world) -> list[str]:
    """World-level node-cloud indexes that no longer describe the current node cloud."""
    points, _ = plates_mod.gather_node_positions(world.plates)
    stale = []
    position = world.node_position_tree_cache
    if position is not None and not (np.array_equal(position[0], points) and np.array_equal(position[1].data, points)):
        stale.append("node_position_tree_cache")
    render = world.node_kdtree_cache
    if render is not None:
        tree = render[3]
        if not np.array_equal(render[0], points):
            stale.append("node_kdtree_cache")
        elif isinstance(tree, cKDTree) and not np.array_equal(tree.data, points):
            stale.append("node_kdtree_cache")
    index = world.node_healpix_index_cache
    if index is not None and len(points):
        fresh = healpix_grid.build_node_pixel_index(index.grid, points)
        if index.grid.nside != healpix_grid.nside_for_node_count(len(points)) or not np.array_equal(index.pixel_to_node, fresh.pixel_to_node):
            stale.append("node_healpix_index_cache")
    return stale


def plate_signature(plate: Plate) -> dict:
    """What revision tracking must account for: a change in the local node set requires a
    topology revision bump, a change in the frame requires a geometry revision bump. Holds a
    reference to the plate so a later audit can tell a mutated plate from a replacement."""
    if isinstance(plate, PlateWithSparseQuadPatch):
        local = hashlib.sha1(np.int64([plate.cells_per_edge]).tobytes() + plate.cell_keys.tobytes()).hexdigest()
    elif plate.node_count():
        local = hashlib.sha1(np.ascontiguousarray(plate.surface_nodes().local_xyz).tobytes()).hexdigest()
    else:
        local = "empty"
    return {
        # Split/defragment build new plate objects that reuse a plate id with fresh counters;
        # revisions are only comparable on the same object. Never serialised.
        "object": plate,
        "topology_revision": int(plate.topology_revision),
        "geometry_revision": int(plate.geometry_revision),
        "nodes": int(plate.node_count()),
        "local": local,
        "frame": hashlib.sha1(np.ascontiguousarray(plate.frame, dtype=np.float64).tobytes()).hexdigest(),
    }


def revision_violations(before: dict | None, after: dict) -> list[str]:
    if before is None or before["object"] is not after["object"]:
        return []
    problems = []
    if (before["local"] != after["local"] or before["nodes"] != after["nodes"]) and before["topology_revision"] == after["topology_revision"]:
        problems.append("local node set changed without a topology revision bump")
    if before["frame"] != after["frame"] and before["geometry_revision"] == after["geometry_revision"]:
        problems.append("frame changed without a geometry revision bump")
    if after["topology_revision"] < before["topology_revision"] or after["geometry_revision"] < before["geometry_revision"]:
        problems.append("revision counter went backwards")
    return problems


def field_violations(plate: Plate) -> list[str]:
    """Hard field problems: wrong shape, non-finite values, non-positive areas."""
    return _field_problems(plate)[0]


def field_bound_violations(plate: Plate) -> list[str]:
    """Values outside elevation/Hc/Hm caps. Kept apart from `field_violations`: the engine
    clamps these caps in specific code paths only."""
    return _field_problems(plate)[1]


def _field_problems(plate: Plate) -> tuple[list[str], list[str]]:
    n = plate.node_count()
    if n == 0:
        return [], []
    nodes = plate.surface_nodes(*SURFACE_FIELDS)
    problems, out_of_bounds = [], []
    for name, values in nodes.fields.items():
        if values.shape != (n,):
            problems.append(f"{name}: shape {values.shape} != ({n},)")
        elif values.dtype.kind == "f" and not np.all(np.isfinite(values)):
            problems.append(f"{name}: {int(np.count_nonzero(~np.isfinite(values)))} non-finite values")
    for label, array in (("area_m2", nodes.area_m2), ("world_xyz", nodes.world_xyz), ("local_xyz", nodes.local_xyz)):
        if not np.all(np.isfinite(array)):
            problems.append(f"{label}: non-finite values")
    if not np.all(nodes.area_m2 > 0):
        problems.append(f"area_m2: {int(np.count_nonzero(nodes.area_m2 <= 0))} non-positive areas")
    bounds = (
        ("elevation", elevation_lines.MIN_ELEVATION_M, elevation_lines.MAX_ELEVATION_M),
        ("crustal_thickness_m", lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M),
        ("mantle_lithosphere_thickness_m", lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M),
    )
    for name, low, high in bounds:
        values = nodes.fields[name]
        # A relative slack of 1e-9 absorbs float rounding at a clamp.
        outside = (values < low - 1e-9 * abs(low)) | (values > high + 1e-9 * abs(high))
        if np.any(outside):
            out_of_bounds.append(f"{name}: {int(outside.sum())} values outside [{low:g}, {high:g}] (range {values.min():.6g}..{values.max():.6g})")
    return problems, out_of_bounds


def frame_violation(plate: Plate) -> str | None:
    frame = np.asarray(plate.frame, dtype=float)
    error = float(np.max(np.abs(frame @ frame.T - np.eye(3))))
    det = float(np.linalg.det(frame))
    if error > 1e-6 or abs(det - 1.0) > 1e-6:
        return f"frame not a proper rotation (|FF^T - I| = {error:.3g}, det = {det:.9f})"
    return None


def quad_cell_corners_local(plate: PlateWithSparseQuadPatch) -> np.ndarray:
    """(n, 4, 3) local unit vectors at every cell's corners, counter-clockwise seen from
    outside the sphere (each face's (u, v) axes are right-handed about its normal)."""
    face, level, i, j = unpack_cell_keys(plate.cell_keys)
    resolution = plate.cells_per_edge * (1 << level)
    corners = [lattice_points(face, i + di, j + dj, resolution) for di, dj in ((0, 0), (1, 0), (1, 1), (0, 1))]
    return np.stack(corners, axis=1)


def quad_cell_shape(plate: PlateWithSparseQuadPatch) -> dict[str, np.ndarray]:
    """Folded flag, aspect ratio and skew per cell, measured in the tangent plane at the cell
    centre. A cell is folded if any corner turns clockwise."""
    corners = quad_cell_corners_local(plate)
    centre = corners.sum(axis=1)
    centre /= np.linalg.norm(centre, axis=1, keepdims=True)
    planar = corners - np.einsum("nkc,nc->nk", corners, centre)[..., None] * centre[:, None, :]
    edges = np.roll(planar, -1, axis=1) - planar
    turn = np.einsum("nkc,nc->nk", np.cross(edges, np.roll(edges, -1, axis=1)), centre)
    lengths = np.linalg.norm(edges, axis=2)
    incoming = -np.roll(edges, 1, axis=1)
    cosine = np.einsum("nkc,nkc->nk", edges, incoming) / np.maximum(lengths * np.roll(lengths, 1, axis=1), 1e-300)
    angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return {
        "folded": np.any(turn <= 0, axis=1),
        "aspect": lengths.max(axis=1) / np.maximum(lengths.min(axis=1), 1e-300),
        "skew_deg": np.max(np.abs(angles - 90.0), axis=1),
    }


def quad_neighbour_violations(plate: PlateWithSparseQuadPatch) -> list[str]:
    """Adjacency must be in range, free of self-loops, symmetric, and join only cells that
    share an edge: two corners of one cell lie on the boundary of the other (exact integer
    cube-surface coordinates, which also covers 2:1 level jumps and cube-face seams)."""
    n = plate.node_count()
    graph = plate.adjacency()
    if len(graph.offsets) != n + 1 or graph.offsets[0] != 0 or graph.offsets[-1] != len(graph.neighbours):
        return ["adjacency offsets malformed"]
    source = np.repeat(np.arange(n), np.diff(graph.offsets))
    target = np.asarray(graph.neighbours, dtype=np.int64)
    problems = []
    if np.any((target < 0) | (target >= n)):
        return [f"{int(np.count_nonzero((target < 0) | (target >= n)))} neighbour indices out of range"]
    if np.any(source == target):
        problems.append(f"{int(np.count_nonzero(source == target))} self-neighbour entries")
    forward = np.sort(source * n + target)
    backward = np.sort(target * n + source)
    if not np.array_equal(forward, backward):
        problems.append("adjacency is not symmetric")
    if len(np.unique(forward)) != len(forward):
        problems.append("duplicate neighbour entries")
    if len(source) == 0:
        return problems
    face, level, i, j = unpack_cell_keys(plate.cell_keys)
    resolution = plate.cells_per_edge * (1 << level)
    common = int(resolution.max())
    offsets = np.array([(0, 0), (1, 0), (1, 1), (0, 1)])
    corners = np.stack(
        [_corner_coordinates(face, i + di, j + dj, resolution, common) for di, dj in offsets], axis=1
    )  # (n, 4, 3) int64
    shares_edge = np.maximum(_corners_on_boundary(corners[source], corners[target]), _corners_on_boundary(corners[target], corners[source])) >= 2
    if not np.all(shares_edge):
        problems.append(f"{int(np.count_nonzero(~shares_edge))} neighbour entries do not share an edge")
    return problems


def _corners_on_boundary(cells: np.ndarray, others: np.ndarray) -> np.ndarray:
    """For each pair, how many of `others`' corners lie on `cells`' boundary (integer exact)."""
    a = cells[:, :, None, :]  # edge start (m, 4, 1, 3)
    b = np.roll(cells, -1, axis=1)[:, :, None, :]
    p = others[:, None, :, :]  # (m, 1, 4, 3)
    edge = b - a
    rel = p - a
    collinear = np.all(np.cross(edge, rel) == 0, axis=-1)
    along = np.einsum("mekc,mekc->mek", np.broadcast_to(edge, rel.shape), rel)
    within = (along >= 0) & (along <= np.einsum("mekc,mekc->mek", edge, edge))
    return np.any(collinear & within, axis=1).sum(axis=1)


def quad_topology_violations(plate: PlateWithSparseQuadPatch) -> list[str]:
    problems = []
    try:
        fresh_copy(plate)._validate_leaf_topology()
    except ValueError as error:
        problems.append(f"leaf topology: {error}")
    keys = plate.cell_keys
    if len(keys) > 1 and np.any(np.diff(keys) <= 0):
        problems.append("cell keys not strictly increasing (duplicate or unsorted cells)")
    return problems


def atmosphere_signature(world) -> dict | None:
    state = world.atmosphere_cfd_state
    if state is None:
        return None
    return {"nside": int(state.grid.nside), "npix": int(state.grid.npix), "temperature_len": int(len(state.temperature_c))}


def audit_world(world, signatures: dict[int, dict]) -> tuple[list[dict], dict[int, dict]]:
    """Every hard invariant for the world's current state. Returns the violations and the new
    per-plate signatures (the caller passes them back in on the next audit)."""
    violations: list[dict] = []

    def add(kind: str, detail: str, plate: Plate | None = None) -> None:
        entry = {"kind": kind, "detail": detail}
        if plate is not None:
            entry["plate_id"] = int(plate.plate_id)
        violations.append(entry)

    new_signatures = {}
    for plate in world.plates:
        hard, bounds = _field_problems(plate)
        for detail in hard:
            add("fields", detail, plate)
        for detail in bounds:
            add("field_bounds", detail, plate)
        problem = frame_violation(plate)
        if problem:
            add("frame", problem, plate)
        stale = stale_plate_caches(plate)
        if stale:
            add("stale_plate_cache", ", ".join(stale), plate)
        signature = plate_signature(plate)
        for detail in revision_violations(signatures.get(plate.plate_id), signature):
            add("revision", detail, plate)
        new_signatures[plate.plate_id] = signature
        if isinstance(plate, PlateWithSparseQuadPatch) and plate.node_count():
            for detail in quad_topology_violations(plate):
                add("quad_topology", detail, plate)
            for detail in quad_neighbour_violations(plate):
                add("quad_neighbours", detail, plate)
            folded = int(quad_cell_shape(plate)["folded"].sum())
            if folded:
                add("quad_folded", f"{folded} folded cells", plate)
    stale = stale_world_caches(world)
    if stale:
        add("stale_world_cache", ", ".join(stale))
    return violations, new_signatures


# --- Checkpoint metrics ---------------------------------------------------------------------


def totals(world) -> dict:
    """Area-weighted totals from each node's actual area (`SurfaceNodes.area_m2`: exact cell
    areas)."""
    names = ("elevation", "crustal_thickness_m", "mantle_lithosphere_thickness_m", "crust_type_code")
    elevation, area, hc, hm, continental, owner_continental = [], [], [], [], [], []
    area_exact = True
    for plate in live_plates(world):
        nodes = plate.surface_nodes(*names)
        area_exact &= bool(nodes.area_is_exact)
        is_continental_plate = plate.crust_type == "continental"
        elevation.append(nodes.fields["elevation"])
        area.append(nodes.area_m2)
        hc.append(nodes.fields["crustal_thickness_m"])
        hm.append(nodes.fields["mantle_lithosphere_thickness_m"])
        continental.append(elevation_lines.effective_is_continental_from_codes(nodes.fields["crust_type_code"], is_continental_plate))
        owner_continental.append(is_continental_plate)
    if not area:
        return {"plates": len(world.plates), "nodes": 0}
    elevation, area, hc, hm, continental = map(np.concatenate, (elevation, area, hc, hm, continental))
    land = elevation > world.sea_level_m
    sphere_m2 = 4.0 * np.pi * (PLANET_RADIUS_KM * 1000.0) ** 2
    spacing_m = line_spacing_rad(world.node_density) * PLANET_RADIUS_KM * 1000.0
    return {
        "plates": len(world.plates),
        "live_plates": len(live_plates(world)),
        "continental_plates": int(sum(owner_continental)),
        "nodes": int(len(area)),
        "area_is_exact": area_exact,
        "area_over_sphere": float(area.sum() / sphere_m2),
        "nominal_area_over_sphere": float(len(area) * spacing_m**2 / sphere_m2),
        "hc_volume_km3": float(np.sum(hc * area) / 1e9),
        "hm_volume_km3": float(np.sum(hm * area) / 1e9),
        "continental_hc_volume_km3": float(np.sum((hc * area)[continental]) / 1e9),
        "continental_area_fraction": float(area[continental].sum() / area.sum()),
        "land_fraction": float(area[land].sum() / area.sum()),
        "sea_level_m": float(world.sea_level_m),
        "elevation": weighted_percentiles(elevation, area),
        "land_elevation": weighted_percentiles(elevation[land], area[land]),
    }


def coverage(world, sample: np.ndarray, node_tree: cKDTree | None) -> dict:
    """Whole-sphere coverage against each plate's own territory (`contains_batch`), voids
    (uncovered and farther than 1.5 s from any node), and nodes inside another plate."""
    spacing = line_spacing_rad(world.node_density)
    count = np.zeros(len(sample), dtype=np.int32)
    live = live_plates(world)
    for plate in live:
        count += plate.contains_batch(sample)
    uncovered = count == 0
    void = np.zeros(len(sample), dtype=bool)
    if node_tree is not None and np.any(uncovered):
        distance, _ = node_tree.query(sample[uncovered], workers=-1)
        void[uncovered] = distance > VOID_SPACING * spacing
    inside = total = 0
    for plate in live:
        points, _ = plate.all_points_and_elevation()
        others = plate.get_neighbours([p for p in live if p is not plate])
        mask = np.zeros(len(points), dtype=bool)
        for other in others:
            mask |= other.contains_batch(points)
        inside += int(mask.sum())
        total += len(points)
    return {
        "uncovered": fraction(uncovered),
        "multiply_covered": fraction(count >= 2),
        "void": fraction(void),
        "nodes_inside_other_plate": inside / total if total else None,
    }


def _tangent_coordinates(points: np.ndarray, neighbours: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Neighbour offsets in each node's local tangent plane, along plate-local east (the line
    row direction, `cross(z, p)`) and north. `points` (n, 3), `neighbours` (n, k, 3), both in
    plate-local coordinates."""
    east = np.cross(np.array([0.0, 0.0, 1.0]), points)
    norm = np.linalg.norm(east, axis=1, keepdims=True)
    east = np.where(norm > 1e-12, east / np.maximum(norm, 1e-300), np.array([1.0, 0.0, 0.0]))
    north = np.cross(points, east)
    offsets = neighbours - points[:, None, :]
    return np.einsum("nkc,nc->nk", offsets, east), np.einsum("nkc,nc->nk", offsets, north)


def sample_cloud(world) -> dict:
    """Representation-neutral lattice quality (docs/plate-surface-baseline.md section 1.1):
    stacked nodes, neighbourhood anisotropy and its alignment with plate-local rows, boundary
    and thin (tendril) nodes, and connectivity."""
    spacing = line_spacing_rad(world.node_density)
    counts = defaultdict(int)
    alignment = []
    for plate in live_plates(world):
        local = plate.surface_nodes().local_xyz
        n = len(local)
        counts["nodes"] += n
        components = plates_mod.node_components(local, elevation_lines.DEFRAG_CONNECT_RADIUS_MULT * spacing)
        sizes = np.bincount(components)
        counts["extra_components"] += max(0, len(sizes) - 1)
        counts["isolated"] += int(np.count_nonzero(sizes == 1))
        if n < 3:
            counts["boundary"] += n
            counts["thin"] += n
            continue
        tree = cKDTree(local)
        k = min(BOUNDARY_NEIGHBOURS + 1, n)
        distance, index = tree.query(local, k=k)
        distance, index = distance[:, 1:], index[:, 1:]
        stacked = distance < STACKED_SPACING * spacing
        counts["stacked"] += int(np.count_nonzero(stacked[:, 0]))
        a, b = _tangent_coordinates(local, local[index])

        near = ~stacked[:, :ANISOTROPY_NEIGHBOURS]
        aa, bb = a[:, :ANISOTROPY_NEIGHBOURS], b[:, :ANISOTROPY_NEIGHBOURS]
        weight = np.maximum(near.sum(axis=1), 1)
        caa = (aa * aa * near).sum(axis=1) / weight
        cbb = (bb * bb * near).sum(axis=1) / weight
        cab = (aa * bb * near).sum(axis=1) / weight
        spread = np.sqrt((caa - cbb) ** 2 + 4 * cab**2)
        major = (caa + cbb + spread) / 2
        minor = np.maximum((caa + cbb - spread) / 2, 0.0)
        ratio = np.minimum(np.sqrt(major / np.maximum(minor, major / ANISOTROPY_CAP**2 + 1e-300)), ANISOTROPY_CAP)
        anisotropic = ratio > ANISOTROPIC_RATIO
        counts["anisotropic"] += int(anisotropic.sum())
        alignment.append(((caa - cbb) / np.maximum(spread, 1e-300))[anisotropic])

        valid = ~stacked & (distance <= BOUNDARY_SPACING * spacing)
        angle = np.where(valid, np.arctan2(b, a), np.nan)
        angle = np.sort(angle, axis=1)  # NaNs last
        valid_count = valid.sum(axis=1)
        gaps_between = np.diff(angle, axis=1)
        gaps_between = np.where(np.arange(angle.shape[1] - 1)[None, :] < (valid_count - 1)[:, None], gaps_between, 0.0)
        first = angle[:, 0]
        last = np.take_along_axis(angle, np.maximum(valid_count - 1, 0)[:, None], axis=1)[:, 0]
        wrap = 2 * np.pi - (last - first)
        largest = np.maximum(np.nan_to_num(gaps_between, nan=0.0).max(axis=1), np.nan_to_num(wrap, nan=2 * np.pi))
        boundary = (valid_count < 2) | (largest > BOUNDARY_SECTOR_RAD)
        interior_neighbour = np.any(valid & ~boundary[index], axis=1)
        counts["boundary"] += int(boundary.sum())
        counts["thin"] += int(np.count_nonzero(boundary & ~interior_neighbour))
    total = max(counts["nodes"], 1)
    alignment = np.concatenate(alignment) if alignment else np.zeros(0)
    return {
        "stacked": counts["stacked"] / total,
        "anisotropic": counts["anisotropic"] / total,
        "anisotropic_row_alignment": float(alignment.mean()) if len(alignment) else None,
        "boundary": counts["boundary"] / total,
        "thin": counts["thin"] / total,
        "extra_components": counts["extra_components"],
        "isolated_nodes": counts["isolated"],
    }


def _loop_orientation(loop: np.ndarray) -> float:
    """+ for a counter-clockwise (outer) loop about its mean direction, - for a hole."""
    centre = loop.mean(axis=0)
    centre /= np.linalg.norm(centre)
    return float(np.sum(np.einsum("ij,j->i", np.cross(loop, np.roll(loop, -1, axis=0)), centre)))


def hole_vertices(world) -> np.ndarray:
    points = []
    for plate in live_plates(world):
        for loop in plate.boundary_loops_world():
            if len(loop) >= 3 and _loop_orientation(loop) < 0:
                points.append(loop)
    return np.concatenate(points) if points else np.zeros((0, 3))


def quad_lattice(world) -> dict | None:
    quads = [p for p in live_plates(world) if isinstance(p, PlateWithSparseQuadPatch)]
    if not quads:
        return None
    cells = refined = thin = spur = holes = folded = aspect = skew = 0
    max_level = 0
    for plate in quads:
        _, level, _, _ = unpack_cell_keys(plate.cell_keys)
        exposed = plate.exposed_sides()  # -v, +u, +v, -u
        degree = np.diff(plate.adjacency().offsets)
        shape = quad_cell_shape(plate)
        cells += plate.node_count()
        refined += int((level > 0).sum())
        max_level = max(max_level, int(level.max(initial=0)))
        thin += int(((exposed[:, 0] & exposed[:, 2]) | (exposed[:, 1] & exposed[:, 3])).sum())
        spur += int((degree <= 1).sum())
        holes += sum(1 for loop in plate._local_boundary_loops() if len(loop) >= 3 and _loop_orientation(loop) < 0)
        folded += int(shape["folded"].sum())
        aspect += int((shape["aspect"] > ASPECT_LIMIT).sum())
        skew += int((shape["skew_deg"] > SKEW_LIMIT_DEG).sum())
    return {
        "cells": cells,
        "refined": refined / cells,
        "max_level": max_level,
        "one_cell_thin": thin / cells,
        "degree_le_1": spur / cells,
        "hole_loops": holes,
        "folded": folded,
        "aspect_gt_4": aspect / cells,
        "skew_gt_45": skew / cells,
    }


def climate_hydrology(world) -> dict:
    stats = world.stats_history[-1] if world.stats_history else {}
    climate = {name: stats.get(name) for name in CHECKPOINT_STATS}
    biomes = stats.get("biome_land_fraction") or {}
    result = {"stats": climate, "biome_land_fraction": dict(sorted(biomes.items()))}
    fields = world.hydrology_cache
    if fields is not None and len(fields.points):
        land = ~np.asarray(fields.is_ocean, dtype=bool)
        land_count = max(int(land.sum()), 1)
        result["hydrology"] = {
            "river_fraction_of_land": float(np.count_nonzero(fields.is_river & land) / land_count),
            "lake_fraction_of_land": float(np.count_nonzero((fields.lake_depth > 0) & land) / land_count),
            "glacier_fraction_of_land": float(np.count_nonzero((fields.glacier_depth > 0) & land) / land_count),
            "land_sink_fraction": float(np.count_nonzero((fields.flow_target < 0) & land) / land_count),
            # `filled_elevation` is +inf by design for a basin with no spill rim
            # (lakes.compute_spill_routing); anything else non-finite is a real fault.
            "finite": bool(
                all(np.all(np.isfinite(np.asarray(getattr(fields, name), dtype=float))) for name in ("flow_accum", "lake_depth", "glacier_depth"))
                and not np.any(np.isnan(fields.filled_elevation) | np.isneginf(fields.filled_elevation))
            ),
        }
    return result


def index_parity(world, sample: np.ndarray) -> tuple[dict, dict]:
    """Exact KD-tree parity (the world's cached node tree, when present, answers exactly like a
    fresh one) and HEALPix nearest-value resampling against the KD-tree, overall and near
    coasts, plate boundaries, holes, poles and the antimeridian. Returns (metrics, timings)."""
    points, plates_in_order = plates_mod.gather_node_positions(world.plates)
    if len(points) == 0:
        return {}, {}
    spacing = line_spacing_rad(world.node_density)
    elevation = plates_mod.collect_all_elevation(plates_in_order)
    owner = np.repeat(np.arange(len(plates_in_order)), [p.node_count() for p in plates_in_order])
    timings = {}

    started = time.perf_counter()
    tree = cKDTree(points)
    timings["kdtree_build_s"] = time.perf_counter() - started
    started = time.perf_counter()
    distance, nearest = tree.query(sample, workers=-1)
    timings["kdtree_query_s"] = time.perf_counter() - started
    timings["kdtree_bytes"] = int(tree.data.nbytes + tree.indices.nbytes)

    metrics: dict = {"nodes": int(len(points)), "nearest_distance_spacing": {k: v / spacing for k, v in _percentiles(distance).items()}}
    cached = world.node_position_tree_cache
    if cached is not None and cached[0].shape == points.shape:
        cached_distance, cached_nearest = cached[1].query(sample, workers=-1)
        metrics["cached_kdtree_mismatches"] = int(np.count_nonzero((cached_nearest != nearest) & (cached_distance != distance)))
    else:
        metrics["cached_kdtree_mismatches"] = None

    started = time.perf_counter()
    grid = healpix_grid.build(healpix_grid.nside_for_node_count(len(points)))
    timings["healpix_grid_build_s"] = time.perf_counter() - started
    started = time.perf_counter()
    index = healpix_grid.build_node_pixel_index(grid, points)
    timings["healpix_index_build_s"] = time.perf_counter() - started
    started = time.perf_counter()
    _, approximate = index.query(sample)
    timings["healpix_query_s"] = time.perf_counter() - started
    timings["healpix_index_bytes"] = int(index.pixel_to_node.nbytes)
    timings["healpix_grid_bytes"] = int(sum(v.nbytes for v in vars(grid).values() if isinstance(v, np.ndarray)))

    approximate_distance = np.linalg.norm(points[approximate] - sample, axis=1)
    ratio = approximate_distance / np.maximum(distance, 1e-12)
    error = np.abs(elevation[approximate] - elevation[nearest])
    wrong_owner = owner[approximate] != owner[nearest]
    wrong_side = (elevation[approximate] > world.sea_level_m) != (elevation[nearest] > world.sea_level_m)

    _, near8 = tree.query(sample, k=min(8, len(points)), workers=-1)
    near8 = near8.reshape(len(sample), -1)
    land8 = elevation[near8] > world.sea_level_m
    owner8 = owner[near8]
    lat = np.degrees(np.arcsin(np.clip(sample[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(sample[:, 1], sample[:, 0]))
    holes = hole_vertices(world)
    near_hole = np.zeros(len(sample), dtype=bool)
    if len(holes):
        hole_distance, _ = cKDTree(holes).query(sample, workers=-1)
        near_hole = hole_distance <= CATEGORY_SPACING * spacing
    categories = {
        "all": np.ones(len(sample), dtype=bool),
        "coast": land8.any(axis=1) & ~land8.all(axis=1),
        "plate_boundary": np.any(owner8 != owner8[:, :1], axis=1),
        "hole": near_hole,
        "pole": np.abs(lat) >= POLE_LAT_DEG,
        "antimeridian": np.abs(lon) >= ANTIMERIDIAN_LON_DEG,
    }
    metrics["healpix"] = {
        "nside": int(grid.nside),
        "categories": {
            name: {
                "samples": int(mask.sum()),
                "same_node": fraction(approximate[mask] == nearest[mask]),
                "distance_ratio": _percentiles(ratio[mask], (50, 95, 99)),
                "distance_spacing": _percentiles(approximate_distance[mask] / spacing, (50, 95, 99)),
                "elevation_error_m": _percentiles(error[mask], (50, 95, 99)),
                "wrong_plate": fraction(wrong_owner[mask]),
                "wrong_land_sea": fraction(wrong_side[mask]),
            }
            for name, mask in categories.items()
        },
    }
    return metrics, timings


def _percentiles(values: np.ndarray, qs=(50, 95, 99, 100)) -> dict:
    if len(values) == 0:
        return {f"p{q}": None for q in qs}
    return {f"p{q}": float(np.percentile(values, q)) for q in qs}


def checkpoint_metrics(world, sample: np.ndarray) -> tuple[dict, dict]:
    points, _ = plates_mod.gather_node_positions(world.plates)
    tree = cKDTree(points) if len(points) else None
    parity, timings = index_parity(world, sample)
    metrics = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "totals": totals(world),
        "coverage": coverage(world, sample, tree),
        "sample_cloud": sample_cloud(world),
        "quad_lattice": quad_lattice(world),
        "climate_hydrology": climate_hydrology(world),
        "index_parity": parity,
        "atmosphere": atmosphere_signature(world),
    }
    return metrics, timings


# --- Load checks ----------------------------------------------------------------------------

_WORLD_DERIVED = ("node_kdtree_cache", "node_position_tree_cache", "node_healpix_index_cache", "node_kdtree_relief_cache", "node_hillshade_cache")


def load_check(world) -> tuple[dict, object]:
    """Save/load round trip: authoritative state identical, every derived index dropped.
    Returns (result, loaded world)."""
    before = state_hash(world)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    leftover = [name for name in _WORLD_DERIVED if getattr(loaded, name, None) is not None]
    for plate in loaded.plates:
        leftover += [f"plate {plate.plate_id}.{name}" for name in cache_builders(plate) if getattr(plate, name, None) is not None]
    return {"state_identical": state_hash(loaded) == before, "derived_state_after_load": leftover}, loaded


# --- Runs -----------------------------------------------------------------------------------


def run_name(seed: int) -> str:
    return f"seed{seed}"


def run_world(config: RunConfig, seed: int, out_dir: Path | None = None, log=print, initial_world: Path | None = None) -> tuple[dict, dict]:
    """Generate and step one world, returning (metrics document, timings document) and writing
    both under `out_dir` when given. With `initial_world`, continue a saved `.mbworld`
    instead of generating one: `seed` must match the save, checkpoint ages count from the
    save's own age, and the world's own densities apply (the config's are ignored). A
    line-backed save converts on load, like any other load."""
    checkpoint_steps = set(config.checkpoint_steps)
    sample = fibonacci_sphere(config.samples)
    instrumentation = Instrumentation()
    generation_kwargs = {"seed": seed, "node_density": config.node_density, "fluid_density": config.fluid_density}
    if config.climate_density is not None:
        generation_kwargs["climate_density"] = config.climate_density

    document: dict = {
        "schema_version": RESULT_SCHEMA_VERSION,
        "config": config.to_json(),
        "seed": seed,
        "checkpoints": [],
        "violations": [],
        "audits": 0,
        "series": [],
        "load_checks": [],
    }
    if initial_world is not None:
        data = Path(initial_world).read_bytes()
        document["initial_world"] = {"name": Path(initial_world).name, "sha256": hashlib.sha256(data).hexdigest()}
    timings: dict = {"schema_version": RESULT_SCHEMA_VERSION, "seed": seed, "steps": [], "checkpoints": []}
    renders_dir = out_dir / "renders" if out_dir is not None and config.render else None

    with instrumentation.installed():
        started = time.perf_counter()
        if initial_world is None:
            world = world_mod.generate_world(**generation_kwargs)
        else:
            world = persistence.load_world_bytes(data)
            del data
            if world.seed != seed:
                raise ValueError(f"{initial_world} is seed {world.seed}, not seed {seed}")
        timings["generation_s"] = time.perf_counter() - started
        signatures: dict[int, dict] = {}
        atmosphere = atmosphere_signature(world)
        for step in range(config.steps + 1):
            if step > 0:
                instrumentation.reset()
                revisions_before = {p.plate_id: (p.topology_revision, p.geometry_revision) for p in world.plates}
                started = time.perf_counter()
                world_mod.step_world(world, config.step_years)
                wall = time.perf_counter() - started
                bumps = [0, 0]
                for plate in world.plates:
                    old = revisions_before.get(plate.plate_id, (0, 0))
                    bumps[0] += max(0, plate.topology_revision - old[0])
                    bumps[1] += max(0, plate.geometry_revision - old[1])
                timings["steps"].append({"step": step, "wall_s": wall, "topology_revision_bumps": bumps[0], "geometry_revision_bumps": bumps[1], **instrumentation.snapshot()})
                stats = world.stats_history[-1] if world.stats_history else {}
                document["series"].append(
                    {"step": step, "nodes": sum(p.node_count() for p in world.plates), "plates": len(world.plates), **{name: stats.get(name) for name in SERIES_STATS}}
                )
            is_checkpoint = step in checkpoint_steps
            if step % config.audit_every == 0 or is_checkpoint:
                violations, signatures = audit_world(world, signatures)
                current = atmosphere_signature(world)
                if current != atmosphere:
                    violations.append({"kind": "atmosphere_grid", "detail": f"{atmosphere} -> {current}"})
                    atmosphere = current
                for violation in violations:
                    violation["step"] = step
                document["violations"].extend(violations)
                document["audits"] += 1
                if violations:
                    log(f"{run_name(seed)} step {step}: {len(violations)} violation(s): {violations[:3]}")
            if is_checkpoint:
                started = time.perf_counter()
                metrics, index_timings = checkpoint_metrics(world, sample)
                metrics["step"] = step
                document["checkpoints"].append(metrics)
                timing_row = {"step": step, "metrics_s": time.perf_counter() - started, "index": index_timings}
                totals_row = metrics["totals"]
                log(
                    f"{run_name(seed)} {world.elapsed_years / 1e6:g} Myr: {totals_row.get('plates')} plates, "
                    f"{totals_row.get('nodes')} nodes, land {totals_row.get('land_fraction', 0):.3f}, "
                    f"uncovered {metrics['coverage']['uncovered']:.4f}"
                )
                if renders_dir is not None:
                    timing_row["render_s"] = _render(world, renders_dir, run_name(seed), step, config.render_size)
                if config.load_checks:
                    result, loaded = load_check(world)
                    result["step"] = step
                    if step == config.steps:
                        # Steps both copies once more, so it runs after everything else here.
                        result["continuation_identical"] = _continuation_identical(world, loaded, config.step_years)
                    document["load_checks"].append(result)
                timings["checkpoints"].append(timing_row)
    timings["total_s"] = timings["generation_s"] + sum(row["wall_s"] for row in timings["steps"])
    if out_dir is not None:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / f"{run_name(seed)}.json").write_text(dumps(document))
        (out_dir / f"{run_name(seed)}.timings.json").write_text(dumps(timings))
    return document, timings


def _continuation_identical(world, loaded, years: float) -> bool:
    """One more step from the in-memory world (warm caches) and from its loaded copy (every
    derived index rebuilt lazily) must land on identical authoritative state."""
    world_mod.step_world(world, years)
    world_mod.step_world(loaded, years)
    return state_hash(world) == state_hash(loaded)


def _render(world, renders_dir: Path, name: str, step: int, size: tuple[int, int]) -> float:
    from . import render_image

    renders_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    for view in ("elevation", "platesDetail"):
        png = render_image.render_png(world, "eckert4", view, size[0], size[1])
        (renders_dir / f"{name}-{world.elapsed_years / 1e6:g}myr-{view}.png").write_bytes(png)
    return time.perf_counter() - started


def load_results(out_dir: Path) -> list[tuple[dict, dict | None]]:
    runs = []
    for path in sorted(out_dir.glob("seed*.json")):
        if path.name.endswith(".timings.json"):
            continue
        timings_path = path.with_name(path.stem + ".timings.json")
        runs.append((json.loads(path.read_text()), json.loads(timings_path.read_text()) if timings_path.exists() else None))
    return runs
