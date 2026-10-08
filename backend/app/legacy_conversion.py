"""One-way conversion of a line-backed world to sparse quads (issue #248, Phase 5b of #228).

docs/save-compatibility.md is the policy this implements; this docstring covers the method.

**Reading.** Line state is read structurally from each pickled object's `__dict__` -- a
plate's `_lines`, `_frame`, `_omega` and so on, and each line's `_phi`, `_theta` and
`_<field>` arrays -- never through line-class methods. A field a line predates (which the
retired `ElevationLine` backfilled lazily) reads as its `surface_fields.SURFACE_FIELDS`
default, which is the same value. State this build doesn't recognise is refused, not dropped.
The line classes themselves were retired in #251: `legacy_unpickler` loads them as inert
`LegacyRecord`s, and this module is the only code that still reads line state.

**Node areas.** Each node stands for the nominal footprint at the world's spacing
(`lithosphere.node_area_m2`, what the line engine's own budgets use), shared with the nodes of
its plate stacked within half a spacing of it -- duplicated row stubs, issue #230. These deduplicated areas weight every average and are the reference the
conversion's volume deltas are reported against.

**Territory.** Each plate keeps its id, frame and rigid motion, and gets a level-0
cube-sphere lattice in that frame at the world's spacing (`cells_per_face_edge`), so a
converted world carries about as many cells as a freshly generated quad world. A cell within
`CANDIDATE_DILATION_CELLS` of one of the plate's nodes becomes active on it when:

- its centre lies in some node's footprint (`FOOTPRINT_SPACINGS`), or in a pinhole enclosed by
  nodes on every side (`ENCLOSED_SPACINGS`). An edge facing open sphere therefore moves by
  about half a cell either way, and a real gap stays a gap; and
- this plate wins the area-weighted vote of the `OWNERSHIP_NEIGHBOURS` nearest nodes, of any
  plate, there (ties go to the plate of the nearest node).

The vote, rather than plain nearest-node ownership, keeps a territory overlap -- where two
line plates' rows interleave -- from turning into salt-and-pepper cells. Every plate is judged
against the same node cloud, so no cell centre is claimed twice; the plates' lattices are
rotated against each other, so their cell edges still overlap a little at boundaries, which
is the ordinary quad-world overlap (issue #255).

**Fields.** Each node goes to one target cell: the cell of its own plate's lattice it lies in
when that cell is active, otherwise the nearest active cell of any plate.

- A cell combines its own plate's nodes by each field's `RemapClass`, weighted by node area:
  the rules `PlateWithSparseQuadPatch.coarsen_cells` uses, applied to point samples.
  Extensive fields are thicknesses, so this is the area-weighted mean thickness, then scaled
  by one bounded ratio per plate and field so the cells hold exactly the volume of the nodes
  they took in (as `quad_merge` restores volume). Elevation is isostasy from the new column
  plus the carried non-isostatic residual. A cell no node targets
  copies its nearest same-plate node, except for volcanoes, which are point features.
- A node whose target belongs to another plate is crust the line world held twice, in an
  overlap. It stacks onto that cell as `quad_merge` stacks a suture, every field by its remap
  class: its extensive volume adds to the cell; volcano flags OR; countdowns and channel depth
  take the max; history takes the earliest set value; composition is voted by crust volume and
  other categories by area, ties keeping the cell's own; intensive and clock fields blend by
  area (soil contents by soil depth, channel width by channel depth); and the non-isostatic
  elevation residual blends by area.
- Thickness stacking pushes past the Hc/Hm caps spreads onto neighbouring cells with room,
  out to `OVERFLOW_SPREAD_RINGS`, so a suture thickens a belt rather than one column. Hc and
  Hm are then clamped into their caps (the suture-accretion limit), shifting elevation by the
  isostatic change. Clamping is the only place volume leaves, and it is reported.

**World.** Plate ids, `collision_progress`/`overlap_progress`, faults (plate id plus plate-local
traces), collision evidence and front records (plate id plus plate-local points), magma parcels
and the world-space logs carry over unchanged. Node-order-keyed caches
are dropped. The line-row tracker `_leading_row_retreat_years` has no quad meaning and is not
carried. The eustatic water budget is re-snapshotted at the current sea level against the
quad hypsometry, so the shoreline doesn't jump on the first step. `World.surface_conversion`
keeps the report's summary, and the event log records the conversion.
"""

from __future__ import annotations

import io
import pickle
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import geometry, lithosphere, rheology
from .elevation_lines import (
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_INHERIT,
    CRUST_TYPE_OCEANIC,
    DEFRAG_CONNECT_RADIUS_MULT,
    line_spacing_rad,
)
from .sparse_quad_patch import (
    PLANET_RADIUS_M,
    PlateWithSparseQuadPatch,
    cell_areas_sr,
    cells_per_face_edge,
    is_structural_reason,
    lattice_points,
    locate_cells,
    pack_cell_keys,
    unpack_cell_keys,
)
from .surface_fields import SURFACE_FIELDS, RemapClass

if TYPE_CHECKING:
    from .world import World

# A cell centre within this many line spacings of a node (of any plate) is in that node's
# footprint. 0.6 makes a straight edge of nodes one spacing apart claim, on average, the half
# spacing beyond its last row that the nodes' own footprints cover; interior points farther
# from every node than this are handled by the enclosure test below.
FOOTPRINT_SPACINGS = 0.6

# A centre beyond every footprint but within this many spacings of nodes on all sides -- no
# empty angular sector of `ENCLOSED_MAX_GAP_RAD` or more among its `OWNERSHIP_NEIGHBOURS` -- is
# an enclosed pinhole: the middle of a square of nodes, a stretched row's gap, or a sliver
# between two plates' outer rows. It is filled. A point beyond an open edge sees all its
# neighbours on one side, so it never passes. The reach is the void threshold of
# docs/plate-surface-baseline.md §1.3.
ENCLOSED_SPACINGS = 1.5
ENCLOSED_MAX_GAP_RAD = np.pi

# How many nearest nodes vote on a cell's owner. Enough to smooth a territory overlap's
# interleaved rows into a contiguous boundary, few enough to stay local (about one ring).
OWNERSHIP_NEIGHBOURS = 8

# Nodes of one plate closer than this many spacings stand for the same patch of crust (the
# stacked fraction of docs/plate-surface-baseline.md §1.1) and share one footprint.
STACKED_SPACINGS = 0.5

# Bounds on the per-plate ratio that restores each extensive field's volume after the mean
# transfer. Healthy plates need well under 5%; the bound only stops a degenerate plate (a
# handful of nodes, mostly stacked away) being inflated or drained.
CONSERVATION_RATIO_BOUNDS = (0.8, 1.25)

# How many rings of neighbouring cells a stacked column's excess over its cap may spread
# across before the remainder is clamped (removed, as the suture-accretion cap does). Twelve
# rings is ~750 km at the default density: a collision plateau's width, which a converted
# collision between two overlapping line plates can need (measured on the issue #248 saves).
OVERFLOW_SPREAD_RINGS = 12

# How far (in cells) around each node's own cell to look for candidate cells. Two rings reach
# every cell whose centre can be within ENCLOSED_SPACINGS of a node, since cells are at least
# ~0.75 spacings across.
CANDIDATE_DILATION_CELLS = 2

# Coverage in the report is measured on this many Fibonacci-sphere samples.
COVERAGE_SAMPLES = 200_000

# The legacy classes a line-backed save pickles. Retiring them (#251) leaves the converter
# reading their state through `legacy_unpickler`. `_RowLookup` is a line plate's cached row
# index (`_row_lookup_cache`, dropped on conversion), pickled by any plate whose containment
# fast path had run -- 20 of the 51 real saves in docs/save-compatibility.md hold one.
LEGACY_LINE_CLASSES = frozenset(
    {
        ("app.plates", "PlateWithLines"),
        ("app.plates", "_RowLookup"),
        ("app.lithosphere_plate", "LithospherePlate"),
        ("app.elevation_lines", "ElevationLine"),
    }
)

# Every attribute a pickled line plate may carry. Revisions and caches are derived state, and
# `_leading_row_retreat_years` tracks retreat at a row's extreme ends, which has no quad
# meaning; none of them carry over. Anything else is state this build doesn't know how to
# convert, so conversion refuses rather than drop it.
_LINE_PLATE_STATE = frozenset({"_plate_id", "_frame", "_crust_type", "_omega", "_age_steps", "_internal_stress", "_lines"})
_LINE_PLATE_DROPPED = frozenset(
    {
        "_topology_revision",
        "_geometry_revision",
        "_bounding_polygon_cache",
        "_bounding_polygon_tree_cache",
        "_node_kdtree_cache",
        "_row_lookup_cache",
        "_world_points_cache",
        "_leading_row_retreat_years",
    }
)
_LINE_STATE = frozenset({"_phi", "_theta"} | {f"_{name}" for name in SURFACE_FIELDS})


class LegacyRecord:
    """Inert stand-in for a retired legacy class: unpickling restores its `__dict__` and
    nothing else, which is all the structural reader needs."""

    legacy_class: str = ""

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)


def legacy_unpickler(data: bytes, retired: frozenset[tuple[str, str]] = LEGACY_LINE_CLASSES) -> pickle.Unpickler:
    """An `Unpickler` that loads the `retired` (module, class) pairs as `LegacyRecord`s."""
    stand_ins: dict[tuple[str, str], type] = {}

    class _Unpickler(pickle.Unpickler):
        def find_class(self, module: str, name: str):
            if (module, name) in retired:
                if (module, name) not in stand_ins:
                    stand_ins[(module, name)] = type(name, (LegacyRecord,), {"legacy_class": f"{module}.{name}"})
                return stand_ins[(module, name)]
            return super().find_class(module, name)

    return _Unpickler(io.BytesIO(data))


def is_line_plate(plate: object) -> bool:
    """A line-backed plate, live or a `LegacyRecord`: one that stores `_lines`."""
    return "_lines" in getattr(plate, "__dict__", {})


def is_line_world(world: "World") -> bool:
    return any(is_line_plate(p) for p in world.plates)


@dataclass
class LegacyPlate:
    """One line plate's persistent state, read from its pickled attributes."""

    plate_id: int
    frame: np.ndarray
    crust_type: str
    omega: np.ndarray
    age_steps: int
    internal_stress: float
    local_xyz: np.ndarray
    fields: dict[str, np.ndarray]
    line_count: int
    one_node_lines: int

    @property
    def node_count(self) -> int:
        return len(self.local_xyz)


def read_line_plate(plate: object) -> LegacyPlate:
    """The persistent state of a line plate, from `__dict__` alone (see the module docstring).
    Raises `ValueError` on state it doesn't recognise, mismatched field lengths, non-finite
    fields or coordinates, or malformed plate motion (frame, omega, age, stress, crust type)."""
    state = plate.__dict__
    unknown = sorted(set(state) - _LINE_PLATE_STATE - _LINE_PLATE_DROPPED)
    if unknown:
        raise ValueError(f"plate {state.get('_plate_id')!r} carries state this build can't convert: {', '.join(unknown)}")
    lines = list(state["_lines"])
    phi, theta = [], []
    columns: dict[str, list[np.ndarray]] = {name: [] for name in SURFACE_FIELDS}
    for line in lines:
        line_state = line.__dict__
        unknown = sorted(set(line_state) - _LINE_STATE)
        if unknown:
            raise ValueError(f"plate {state.get('_plate_id')!r}: a line carries fields this build can't convert: {', '.join(unknown)}")
        line_theta = np.asarray(line_state["_theta"], dtype=float).reshape(-1)
        phi.append(np.full(len(line_theta), float(line_state["_phi"])))
        theta.append(line_theta)
        for name, spec in SURFACE_FIELDS.items():
            values = line_state.get(f"_{name}")
            if values is None:
                values = np.full(len(line_theta), spec.default, dtype=spec.dtype)
            values = np.asarray(values, dtype=spec.dtype).reshape(-1)
            if len(values) != len(line_theta):
                raise ValueError(f"plate {state.get('_plate_id')!r}: line field {name!r} has {len(values)} values for {len(line_theta)} nodes")
            columns[name].append(values)
    fields = {
        name: np.concatenate(chunks) if chunks else np.zeros(0, dtype=SURFACE_FIELDS[name].dtype)
        for name, chunks in columns.items()
    }
    for name, values in fields.items():
        if values.dtype.kind == "f" and not np.all(np.isfinite(values)):
            raise ValueError(f"plate {state.get('_plate_id')!r}: field {name!r} holds non-finite values")
    phi = np.concatenate(phi) if phi else np.zeros(0)
    theta = np.concatenate(theta) if theta else np.zeros(0)
    if not (np.all(np.isfinite(phi)) and np.all(np.isfinite(theta))):
        raise ValueError(f"plate {state.get('_plate_id')!r}: line coordinates hold non-finite values")
    local = geometry.local_xyz(phi, theta) if len(phi) else np.zeros((0, 3))
    frame, omega, age_steps, internal_stress, crust_type = _plate_motion_state(state)
    return LegacyPlate(
        plate_id=int(state["_plate_id"]),
        frame=frame,
        crust_type=crust_type,
        omega=omega,
        age_steps=age_steps,
        internal_stress=internal_stress,
        local_xyz=np.asarray(local, dtype=float).reshape(-1, 3),
        fields=fields,
        line_count=len(lines),
        one_node_lines=sum(len(line.__dict__["_theta"]) == 1 for line in lines),
    )


# How far a pickled frame may be from a proper rotation (the stress tests' own tolerance).
FRAME_TOLERANCE = 1e-6


def _plate_motion_state(state: dict) -> tuple[np.ndarray, np.ndarray, int, float, str]:
    """(frame, omega, age_steps, internal_stress, crust_type) of a pickled plate, checked:
    a converted world would carry anything malformed here straight into its simulation."""
    plate_id = state.get("_plate_id")
    frame = np.asarray(state["_frame"], dtype=float)
    if frame.shape != (3, 3) or not np.all(np.isfinite(frame)):
        raise ValueError(f"plate {plate_id!r}: frame must be a finite 3x3 matrix")
    if not np.allclose(frame @ frame.T, np.eye(3), atol=FRAME_TOLERANCE) or np.linalg.det(frame) <= 0.0:
        raise ValueError(f"plate {plate_id!r}: frame is not a proper rotation")
    omega = np.asarray(state.get("_omega", np.zeros(3)), dtype=float)
    if omega.shape != (3,) or not np.all(np.isfinite(omega)):
        raise ValueError(f"plate {plate_id!r}: omega must be a finite 3-vector")
    internal_stress = float(state.get("_internal_stress", 0.0))
    if not np.isfinite(internal_stress):
        raise ValueError(f"plate {plate_id!r}: internal stress is not finite")
    age_steps = state.get("_age_steps", 0)
    if isinstance(age_steps, bool) or not isinstance(age_steps, (int, np.integer)) or age_steps < 0:
        raise ValueError(f"plate {plate_id!r}: age must be a non-negative step count")
    crust_type = state["_crust_type"]
    if crust_type not in ("continental", "oceanic"):
        raise ValueError(f"plate {plate_id!r}: unknown crust type {crust_type!r}")
    return frame.copy(), omega.copy(), int(age_steps), internal_stress, crust_type


@dataclass
class ConversionReport:
    """What a conversion kept and changed. Volumes are m^3 (field value x area); fractions
    are of the whole sphere unless noted. docs/save-compatibility.md has the tolerances."""

    cells_per_edge: int
    spacing_rad: float
    source: dict[str, Any] = field(default_factory=dict)
    result: dict[str, Any] = field(default_factory=dict)
    # Per extensive field: the line total by deduplicated and by nominal node area, the quad
    # total by exact cell area, and the relative delta of the quad total from each.
    extensive: dict[str, dict[str, float]] = field(default_factory=dict)
    # Per plate: node/cell counts, nodes stacked onto or from other plates, the relative area
    # and crust-volume deltas (crust against the line plate clipped to the caps), and components
    # before (nodes linked within the defragment radius) and after (edge-adjacent cells).
    plates: list[dict[str, Any]] = field(default_factory=list)
    # Nodes whose own cell wasn't active on their plate: `strays` joined a nearby cell of their
    # own plate; `stacked` went onto another plate's cell (an overlap).
    targets: dict[str, float] = field(default_factory=dict)
    # Hc/Hm volume removed by clamping into the caps, and the volume the line world itself
    # already held past them (saves from before issue #256; the line engine's own end-of-step
    # clamp removes that too), each as a fraction of the line total.
    clamped: dict[str, float] = field(default_factory=dict)
    source_over_cap: dict[str, float] = field(default_factory=dict)
    # Shares of the covered sphere, before and after, that are land at the current sea level
    # and effectively continental crust, sampled the same way for both (`_sampled_fraction`).
    area_fractions: dict[str, dict[str, float]] = field(default_factory=dict)
    # Provenance violations (all must be 0) and field ranges.
    preserved: dict[str, Any] = field(default_factory=dict)
    coverage: dict[str, Any] = field(default_factory=dict)
    sea_level: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return _plain(self.__dict__)

    def summary(self) -> dict[str, Any]:
        """The compact form kept on `World.surface_conversion`."""
        hc = self.extensive.get("crustal_thickness_m", {})
        return _plain(
            {
                "from": "lines",
                "to": "quad",
                "cells_per_edge": self.cells_per_edge,
                "line_nodes": self.source.get("nodes"),
                "quad_cells": self.result.get("cells"),
                "crust_volume_delta": hc.get("delta_vs_capped"),
                "crust_volume_delta_vs_nominal": hc.get("delta_vs_nominal"),
                "crust_volume_clamped": self.clamped.get("crustal_thickness_m"),
                "stacked_node_fraction": self.targets.get("stacked_fraction"),
                "uncovered_fraction": self.coverage.get("after", {}).get("uncovered"),
            }
        )


def _plain(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _chord(angle_rad: float) -> float:
    return 2.0 * np.sin(0.5 * angle_rad)


def convert_world_to_quad(world: "World", *, coverage_samples: int = COVERAGE_SAMPLES) -> ConversionReport:
    """Convert every line plate of `world` to a `PlateWithSparseQuadPatch` in place and return
    the report. One-way: nothing is kept that could restore the lines. Raises `ValueError` on
    a world with no line plates, one that mixes in other plates, or line state this build
    can't convert faithfully (see `read_line_plate`)."""
    if not is_line_world(world):
        raise ValueError("world has no line-backed plates to convert")
    if any(not is_line_plate(p) for p in world.plates):
        raise ValueError("world mixes line-backed and other plates; refusing to convert part of it")

    spacing = line_spacing_rad(world.node_density)
    n = cells_per_face_edge(spacing)
    report = ConversionReport(cells_per_edge=n, spacing_rad=spacing)
    legacy = [read_line_plate(p) for p in world.plates]
    sea_level = float(world.sea_level_m)
    budget_before = getattr(world, "ocean_water_volume_m3", None)

    world_xyz = [geometry.to_world(p.frame, p.local_xyz) for p in legacy]
    areas = [_node_areas(xyz, spacing) for xyz in world_xyz]
    node_tree = cKDTree(np.concatenate(world_xyz))
    node_plate = np.concatenate([np.full(p.node_count, k) for k, p in enumerate(legacy)]).astype(np.int64)
    node_area = np.concatenate(areas)
    keys = [
        _territory(plate, k, n, node_tree, node_plate, node_area, _chord(FOOTPRINT_SPACINGS * spacing), _chord(ENCLOSED_SPACINGS * spacing))
        for k, plate in enumerate(legacy)
    ]
    targets = _targets(legacy, keys, n)

    converted = []
    conservation_ratios: list[dict[str, float]] = []
    for k, plate in enumerate(legacy):
        t_plate, t_cell = targets[k]
        own = t_plate == k
        fields, ratios = _cell_fields(plate, keys[k], n, t_cell, own, areas[k])
        conservation_ratios.append(ratios)
        converted.append(
            PlateWithSparseQuadPatch(
                plate.plate_id,
                plate.frame.copy(),
                plate.crust_type,
                n,
                keys[k],
                fields=fields,
                omega=plate.omega.copy(),
                age_steps=plate.age_steps,
                internal_stress=plate.internal_stress,
            )
        )
    _stack_foreign(legacy, converted, targets, areas)
    clamped = _clamp_columns(converted)

    nodes = int(node_plate.size)
    own_cell = [_own_cell_mask(plate, keys[k], n) for k, plate in enumerate(legacy)]
    stacked = sum(int(np.count_nonzero(targets[k][0] != k)) for k in range(len(legacy)))
    strays = sum(int(np.count_nonzero((targets[k][0] == k) & ~own_cell[k])) for k in range(len(legacy)))
    report.targets = {
        "strays": strays,
        "stray_fraction": strays / max(nodes, 1),
        "stacked": stacked,
        "stacked_fraction": stacked / max(nodes, 1),
    }
    report.source = {
        "plates": len(legacy),
        "lines": sum(p.line_count for p in legacy),
        "nodes": nodes,
        "one_node_lines": sum(p.one_node_lines for p in legacy),
        "deduplicated_area_m2": float(node_area.sum()),
        "nominal_area_m2": float(lithosphere.node_area_m2(spacing) * nodes),
    }
    for name, spec in SURFACE_FIELDS.items():
        if spec.remap_class != RemapClass.EXTENSIVE:
            continue
        deduplicated = float(sum(np.sum(p.fields[name] * a) for p, a in zip(legacy, areas)))
        nominal = float(sum(np.sum(p.fields[name]) for p in legacy) * lithosphere.node_area_m2(spacing))
        quad = float(sum(np.sum(q.collect(name) * q.node_areas_m2()) for q in converted if q.node_count()))
        report.extensive[name] = {
            "line_deduplicated": deduplicated,
            "line_nominal": nominal,
            "quad": quad,
            "delta_vs_deduplicated": _relative(quad, deduplicated),
            "delta_vs_nominal": _relative(quad, nominal),
        }
        if name in _CAPS:
            # What the conversion itself changed: the line world clipped to the caps it would
            # be clamped to on its next step anyway.
            capped = float(sum(np.sum(np.clip(p.fields[name], *_CAPS[name]) * a) for p, a in zip(legacy, areas)))
            report.extensive[name].update(line_capped=capped, delta_vs_capped=_relative(quad, capped))
    report.clamped = {
        name: volume / report.extensive[name]["line_deduplicated"] if report.extensive[name]["line_deduplicated"] else 0.0
        for name, volume in clamped.items()
    }
    report.source_over_cap = {
        name: float(sum(np.sum((p.fields[name] - np.clip(p.fields[name], *_CAPS[name])) * a) for p, a in zip(legacy, areas)))
        / report.extensive[name]["line_deduplicated"]
        if report.extensive[name]["line_deduplicated"]
        else 0.0
        for name in _CAPS
    }
    if coverage_samples:
        sample = _fibonacci_sphere(coverage_samples)
        reach = _chord(ENCLOSED_SPACINGS * spacing)
        report.area_fractions = {
            name: {
                "before": _sampled_fraction(sample, reach, world_xyz, [p.fields["elevation"] for p in legacy], [_effective_continental(p.fields["crust_type_code"], p.crust_type) for p in legacy], sea_level, name),
                "after": _sampled_fraction(
                    sample,
                    reach,
                    [q.all_points_and_elevation()[0] for q in converted if q.node_count()],
                    [q.collect("elevation") for q in converted if q.node_count()],
                    [_effective_continental(q.collect("crust_type_code"), q.crust_type) for q in converted if q.node_count()],
                    sea_level,
                    name,
                ),
            }
            for name in ("land", "continental")
        }
    for k, (plate, quad) in enumerate(zip(legacy, converted)):
        before = float(np.sum(np.clip(plate.fields["crustal_thickness_m"], *_CAPS["crustal_thickness_m"]) * areas[k]))
        after = float(np.sum(quad.collect("crustal_thickness_m") * quad.node_areas_m2())) if quad.node_count() else 0.0
        report.plates.append(
            {
                "plate_id": plate.plate_id,
                "line_nodes": plate.node_count,
                "lines": plate.line_count,
                "quad_cells": quad.node_count(),
                "crust_conservation_ratio": conservation_ratios[k].get("crustal_thickness_m", 1.0),
                "stacked_nodes": int(np.count_nonzero(targets[k][0] != k)),
                "stacked_onto": int(sum(np.count_nonzero(targets[m][0] == k) for m in range(len(legacy)) if m != k)),
                "area_delta": _relative(float(quad.node_areas_m2().sum()) if quad.node_count() else 0.0, float(areas[k].sum())),
                "crust_volume_delta": _relative(after, before),
                "line_components": _node_components(world_xyz[k], DEFRAG_CONNECT_RADIUS_MULT * spacing),
                "quad_components": _cell_components(quad),
            }
        )
    report.preserved = _preserved(legacy, converted, targets)

    world.plates = converted
    _reset_world_caches(world)
    from . import eustasy

    eustasy.initialize_water_budget(world)
    report.result = {
        "plates": len(converted),
        "cells": sum(q.node_count() for q in converted),
        "area_m2": float(sum(q.node_areas_m2().sum() for q in converted if q.node_count())),
    }
    # "before" would need the retired line plates' own containment test; kept for the
    # report's shape.
    report.coverage = {"before": None, "after": _quad_coverage(converted, coverage_samples)}
    report.sea_level = {
        "before_m": sea_level,
        "after_m": float(world.sea_level_m),
        "water_budget_before_m3": float(budget_before) if budget_before is not None else float("nan"),
        "water_budget_after_m3": float(world.ocean_water_volume_m3),
    }
    world.surface_conversion = report.summary()
    world.log_event(
        f"Converted line-backed world to sparse quads: {nodes} line nodes -> {report.result['cells']} cells "
        f"(crust volume {100 * report.extensive['crustal_thickness_m']['delta_vs_capped']:+.2f}%)."
    )
    return report


def _relative(after: float, before: float) -> float:
    return (after - before) / before if before else 0.0


_CAPS = {
    "crustal_thickness_m": (lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M),
    "mantle_lithosphere_thickness_m": (lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M),
}


def _sampled_fraction(
    sample: np.ndarray, reach: float, points: list[np.ndarray], elevation: list[np.ndarray], continental: list[np.ndarray], sea_level: float, name: str
) -> float:
    """Share of the covered sphere that is land (or continental crust), reading each sample
    point from its nearest node or cell centre within `reach` -- what a render shows. Line
    plates that overlap count the overlap once, as cells do, which per-node area sums can't."""
    dist, nearest = cKDTree(np.concatenate(points)).query(sample, distance_upper_bound=reach)
    covered = np.isfinite(dist)
    flag = np.concatenate(elevation) > sea_level if name == "land" else np.concatenate(continental)
    return float(np.mean(flag[nearest[covered]])) if np.any(covered) else 0.0


def _effective_continental(codes: np.ndarray, crust_type: str) -> np.ndarray:
    return np.where(codes == CRUST_TYPE_INHERIT, crust_type == "continental", codes == CRUST_TYPE_CONTINENTAL)


def _node_areas(world_xyz: np.ndarray, spacing: float) -> np.ndarray:
    """Nominal footprint per node, divided by how many nodes of the plate (itself included)
    lie within `STACKED_SPACINGS`: a local density, so two stacked nodes take half each and a
    row squeezed to 0.4 spacings takes about a third each, without the chaining a connected-
    component grouping has along a stretched or stacked row."""
    nominal = lithosphere.node_area_m2(spacing)
    if len(world_xyz) < 2:
        return np.full(len(world_xyz), nominal)
    counts = cKDTree(world_xyz).query_ball_point(world_xyz, _chord(STACKED_SPACINGS * spacing), return_length=True)
    return nominal / np.maximum(counts, 1)


def _territory(
    plate: LegacyPlate,
    index: int,
    n: int,
    node_tree: cKDTree,
    node_plate: np.ndarray,
    node_area: np.ndarray,
    footprint_chord: float,
    enclosed_chord: float,
) -> np.ndarray:
    """Sorted level-0 cell keys this plate keeps (see the module docstring)."""
    if plate.node_count == 0:
        return np.zeros(0, dtype=np.int64)
    candidates = np.unique(pack_cell_keys(*locate_cells(plate.local_xyz, n)))
    di, dj = (offset.ravel() for offset in np.meshgrid([-1, 0, 1], [-1, 0, 1], indexing="ij"))
    for _ in range(CANDIDATE_DILATION_CELLS):
        face, _, i, j = unpack_cell_keys(candidates)
        probes = lattice_points(np.repeat(face, 9), np.repeat(i, 9) + np.tile(di, len(i)) + 0.5, np.repeat(j, 9) + np.tile(dj, len(j)) + 0.5, n)
        candidates = np.unique(pack_cell_keys(*locate_cells(probes, n)))
    face, _, i, j = unpack_cell_keys(candidates)
    centres = geometry.to_world(plate.frame, lattice_points(face, i + 0.5, j + 0.5, n))
    return candidates[_vote_owners(node_tree, node_plate, node_area, centres, footprint_chord, enclosed_chord) == index]


def _vote_owners(
    tree: cKDTree, node_plate: np.ndarray, node_area: np.ndarray, centres: np.ndarray, footprint_chord: float, enclosed_chord: float
) -> np.ndarray:
    """Owning plate index of each cell centre, or -1 outside every footprint and not enclosed."""
    dist, idx = tree.query(centres, k=OWNERSHIP_NEIGHBOURS, distance_upper_bound=enclosed_chord)
    dist, idx = dist.reshape(len(centres), -1), idx.reshape(len(centres), -1)
    valid = np.isfinite(dist)
    inside = (valid[:, 0] & (dist[:, 0] <= footprint_chord)) | _enclosed(tree.data, centres, idx, valid)
    valid &= inside[:, None]
    rows = np.broadcast_to(np.arange(len(centres))[:, None], valid.shape)[valid]
    plates = node_plate[idx[valid]]
    # A plate's every node has positive area; the epsilon only guards a degenerate tie at 0.
    weights = node_area[idx[valid]] + 1e-9
    nearest = np.where(valid[:, 0], node_plate[np.where(valid[:, 0], idx[:, 0], 0)], -1)
    width = int(node_plate.max(initial=0)) + 1
    pair, inverse = np.unique(rows * width + plates, return_inverse=True)
    totals = np.bincount(inverse, weights=weights)
    row, plate = pair // width, pair % width
    order = np.lexsort((plate != nearest[row], -totals, row))
    first = order[np.r_[True, row[order][1:] != row[order][:-1]]] if len(order) else order
    owner = np.full(len(centres), -1, dtype=np.int64)
    owner[row[first]] = plate[first]
    return owner


def _enclosed(points: np.ndarray, centres: np.ndarray, idx: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Whether each centre's valid neighbours surround it: every empty angular sector between
    them, seen in the centre's tangent plane, is narrower than `ENCLOSED_MAX_GAP_RAD`."""
    if not len(centres):
        return np.zeros(0, dtype=bool)
    reference = np.where(np.abs(centres[:, 2:3]) < 0.9, [[0.0, 0.0, 1.0]], [[1.0, 0.0, 0.0]])
    east = np.cross(reference, centres)
    east /= np.linalg.norm(east, axis=1, keepdims=True)
    north = np.cross(centres, east)
    offsets = points[np.where(valid, idx, 0)] - centres[:, None, :]
    angles = np.arctan2(np.einsum("nkj,nj->nk", offsets, north), np.einsum("nkj,nj->nk", offsets, east))
    count = valid.sum(axis=1)
    angles = np.sort(np.where(valid, angles, np.inf), axis=1)
    # Pad each row's missing neighbours with its own largest angle, so they add no gap.
    last = np.take_along_axis(angles, np.maximum(count - 1, 0)[:, None], axis=1)
    last = np.where(count[:, None] > 0, last, 0.0)
    angles = np.where(np.isfinite(angles), angles, last)
    widest = np.maximum(np.diff(angles, axis=1).max(axis=1), angles[:, 0] + 2.0 * np.pi - last[:, 0])
    return (count >= 3) & (widest < ENCLOSED_MAX_GAP_RAD)


def _own_cell_mask(plate: LegacyPlate, keys: np.ndarray, n: int) -> np.ndarray:
    """Which nodes lie in an active cell of their own plate's lattice."""
    if plate.node_count == 0 or len(keys) == 0:
        return np.zeros(plate.node_count, dtype=bool)
    node_key = pack_cell_keys(*locate_cells(plate.local_xyz, n))
    position = np.minimum(np.searchsorted(keys, node_key), len(keys) - 1)
    return keys[position] == node_key


def _targets(legacy: list[LegacyPlate], keys: list[np.ndarray], n: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """(target plate index, target cell index) per node of each plate: its own cell when
    active, otherwise the nearest active cell of any plate."""
    live = [k for k in range(len(legacy)) if len(keys[k])]
    if not live:
        raise ValueError("conversion found no territory for any plate")
    centres = []
    for k in live:
        face, _, i, j = unpack_cell_keys(keys[k])
        centres.append(geometry.to_world(legacy[k].frame, lattice_points(face, i + 0.5, j + 0.5, n)))
    tree = cKDTree(np.concatenate(centres))
    centre_plate = np.concatenate([np.full(len(keys[k]), k) for k in live])
    centre_cell = np.concatenate([np.arange(len(keys[k])) for k in live])

    targets = []
    for k, plate in enumerate(legacy):
        t_plate = np.full(plate.node_count, -1, dtype=np.int64)
        t_cell = np.full(plate.node_count, -1, dtype=np.int64)
        own = _own_cell_mask(plate, keys[k], n)
        if np.any(own):
            node_key = pack_cell_keys(*locate_cells(plate.local_xyz[own], n))
            t_plate[own], t_cell[own] = k, np.searchsorted(keys[k], node_key)
        if not np.all(own):
            _, hit = tree.query(geometry.to_world(plate.frame, plate.local_xyz[~own]))
            t_plate[~own], t_cell[~own] = centre_plate[hit], centre_cell[hit]
        targets.append((t_plate, t_cell))
    return targets


def _cell_fields(
    plate: LegacyPlate,
    keys: np.ndarray,
    n: int,
    cell_of_node: np.ndarray,
    contributes: np.ndarray,
    node_area: np.ndarray,
) -> dict[str, np.ndarray]:
    """Every surface field on `keys`, combined from the `contributes` nodes by remap class,
    and the per-field conservation ratio applied to each extensive field."""
    count = len(keys)
    if count == 0:
        return {}, {}
    cell = cell_of_node[contributes]
    weight = node_area[contributes]
    src = {name: values[contributes] for name, values in plate.fields.items()}
    has_node = np.bincount(cell, minlength=count) > 0
    total_weight = np.bincount(cell, weights=weight, minlength=count)

    def mean(values: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
        if w is None:
            return np.bincount(cell, weights=weight * values, minlength=count) / np.where(has_node, total_weight, 1.0)
        total = np.bincount(cell, weights=w, minlength=count)
        return np.where(total > 0.0, np.bincount(cell, weights=w * values, minlength=count) / np.where(total > 0.0, total, 1.0), mean(values))

    def vote(values: np.ndarray, w: np.ndarray, prefer=None) -> np.ndarray:
        """Largest total weight wins; ties go to a preferred value, then the lowest."""
        choices, rank = np.unique(values, return_inverse=True)
        pair, inverse = np.unique(cell * len(choices) + rank, return_inverse=True)
        totals = np.bincount(inverse, weights=w)
        pair_cell, pair_rank = pair // len(choices), pair % len(choices)
        not_preferred = np.zeros(len(pair), dtype=bool) if prefer is None else ~prefer(choices)[pair_rank]
        order = np.lexsort((pair_rank, not_preferred, -totals, pair_cell))
        first = order[np.r_[True, pair_cell[order][1:] != pair_cell[order][:-1]]]
        out = np.zeros(count, dtype=values.dtype)
        out[pair_cell[first]] = choices[pair_rank[first]]
        return out

    own_code = CRUST_TYPE_CONTINENTAL if plate.crust_type == "continental" else CRUST_TYPE_OCEANIC
    out: dict[str, np.ndarray] = {}
    for name, spec in SURFACE_FIELDS.items():
        values = src[name]
        if name == "elevation":
            continue
        if name == "crust_type_code":
            # Composition votes by crust volume, as in quad_merge: a column's density follows
            # the crust that makes it up.
            explicit = np.where(values == CRUST_TYPE_INHERIT, own_code, values)
            won = vote(explicit, weight * np.maximum(src["crustal_thickness_m"], 0.0) + 1e-9)
            out[name] = np.where(won == own_code, CRUST_TYPE_INHERIT, won).astype(spec.dtype)
        elif name == "elev_change_reason":
            out[name] = vote(values, weight, prefer=is_structural_reason)
        elif name == "channel_depth" or spec.remap_class == RemapClass.COUNTDOWN:
            peak = np.full(count, -np.inf)
            np.maximum.at(peak, cell, values)
            out[name] = np.where(has_node, peak, spec.default)
        elif name == "channel_width":
            order = np.lexsort((-src["channel_depth"], cell))
            first = order[np.r_[True, cell[order][1:] != cell[order][:-1]]] if len(order) else order
            width = np.full(count, spec.default, dtype=float)
            width[cell[first]] = values[first]
            out[name] = width
        elif spec.remap_class in (RemapClass.HISTORY, RemapClass.WRITE_ONCE_HISTORY):
            earliest = np.full(count, np.inf)
            set_ = values != spec.sentinel
            np.minimum.at(earliest, cell[set_], values[set_])
            out[name] = np.where(np.isfinite(earliest), earliest, spec.sentinel)
        elif spec.remap_class == RemapClass.BOOLEAN_PROVENANCE:
            out[name] = np.bincount(cell, weights=values.astype(float), minlength=count) > 0.0
        elif spec.remap_class == RemapClass.CATEGORICAL:
            out[name] = vote(values, weight)
        elif name in ("soil_mineral_content", "soil_organic_content"):
            out[name] = mean(values, weight * src["soil_depth"])
        else:
            out[name] = mean(values)

    # Elevation: the isostatic column the new Hc/Hm/composition floats at, plus the carried
    # non-isostatic residual (erosion, texture, flexure), as coarsen_cells and quad_merge do.
    density = lithosphere.node_crust_density(src["crust_type_code"], plate.crust_type)
    residual = src["elevation"] - lithosphere.isostatic_elevation(src["crustal_thickness_m"], src["mantle_lithosphere_thickness_m"], density)
    cell_density = lithosphere.node_crust_density(out["crust_type_code"], plate.crust_type)
    out["elevation"] = rheology.clip_elevation_bounds(
        lithosphere.isostatic_elevation(out["crustal_thickness_m"], out["mantle_lithosphere_thickness_m"], cell_density) + mean(residual)
    )

    # A cell no node targets (between rows, or an enclosed pinhole) copies its nearest node,
    # elevation included: a copied column floats where it did. Volcanoes are point features,
    # so their provenance and countdown aren't copied -- that would multiply them.
    empty = ~has_node
    if np.any(empty):
        face, _, i, j = unpack_cell_keys(keys[empty])
        _, nearest = cKDTree(plate.local_xyz).query(lattice_points(face, i + 0.5, j + 0.5, n))
        for name, values in plate.fields.items():
            spec = SURFACE_FIELDS[name]
            point_feature = spec.remap_class in (RemapClass.BOOLEAN_PROVENANCE, RemapClass.COUNTDOWN)
            out[name][empty] = spec.default if point_feature else values[nearest]

    # Conserve every extensive field per plate, as quad_merge does: one ratio restores the
    # volume of the nodes the cells took in (the mean transfer alone gains or loses where an
    # edge moves by part of a cell). Bounded, so a degenerate plate can't be inflated; the
    # report shows the ratios. Elevation follows Hc/Hm isostatically.
    face, _, i, j = unpack_cell_keys(keys)
    cell_area = cell_areas_sr(i, j, n) * PLANET_RADIUS_M**2
    cell_density = lithosphere.node_crust_density(out["crust_type_code"], plate.crust_type)
    ratios = {}
    hc0, hm0 = out["crustal_thickness_m"].copy(), out["mantle_lithosphere_thickness_m"].copy()
    for name, spec in SURFACE_FIELDS.items():
        if spec.remap_class != RemapClass.EXTENSIVE:
            continue
        placed = float(np.sum(out[name] * cell_area))
        owed = float(np.sum(src[name] * weight))
        ratio = float(np.clip(owed / placed, *CONSERVATION_RATIO_BOUNDS)) if placed > 0.0 else 1.0
        out[name] = out[name] * ratio
        ratios[name] = ratio
    shift = lithosphere.isostatic_elevation(out["crustal_thickness_m"], out["mantle_lithosphere_thickness_m"], cell_density) - lithosphere.isostatic_elevation(
        hc0, hm0, cell_density
    )
    out["elevation"] = rheology.clip_elevation_bounds(out["elevation"] + shift)
    return out, ratios


def _stack_foreign(
    legacy: list[LegacyPlate],
    converted: list[PlateWithSparseQuadPatch],
    targets: list[tuple[np.ndarray, np.ndarray]],
    areas: list[np.ndarray],
) -> None:
    """Stack every node whose target cell is on another plate onto that cell, the way
    `quad_merge` stacks a suture (see the module docstring). Caps are left to
    `_clamp_columns`."""
    records = []
    for k, (t_plate, t_cell) in enumerate(targets):
        foreign = np.flatnonzero(t_plate != k)
        if len(foreign):
            records.append((t_plate[foreign], t_cell[foreign], np.full(len(foreign), k), foreign))
    if not records:
        return
    t_plate, t_cell, s_plate, s_node = (np.concatenate(parts) for parts in zip(*records))

    for t in np.unique(t_plate):
        quad = converted[t]
        at = t_plate == t
        cell, src_plate, src_node = t_cell[at], s_plate[at], s_node[at]
        count = quad.node_count()
        area = quad.node_areas_m2()

        def source(name: str) -> np.ndarray:
            values = np.empty(len(cell), dtype=SURFACE_FIELDS[name].dtype)
            for k in np.unique(src_plate):
                mine = src_plate == k
                values[mine] = legacy[k].fields[name][src_node[mine]]
            return values

        node_area = np.empty(len(cell))
        src_explicit = np.empty(len(cell), dtype=np.int8)
        for k in np.unique(src_plate):
            mine = src_plate == k
            node_area[mine] = areas[k][src_node[mine]]
            codes = legacy[k].fields["crust_type_code"][src_node[mine]]
            inherited = CRUST_TYPE_CONTINENTAL if legacy[k].crust_type == "continental" else CRUST_TYPE_OCEANIC
            src_explicit[mine] = np.where(codes == CRUST_TYPE_INHERIT, inherited, codes)
        receiving = np.zeros(count, dtype=bool)
        receiving[cell] = True

        fields = {name: quad.collect(name) for name in SURFACE_FIELDS}
        own = {name: values.copy() for name, values in fields.items()}
        old_hc, old_hm = own["crustal_thickness_m"], own["mantle_lithosphere_thickness_m"]

        # Composition: explicit codes on both sides, voted by crust volume, then re-encoded
        # relative to the receiving plate.
        own_code = CRUST_TYPE_CONTINENTAL if quad.crust_type == "continental" else CRUST_TYPE_OCEANIC
        cell_explicit = np.where(own["crust_type_code"] == CRUST_TYPE_INHERIT, own_code, own["crust_type_code"])
        src_volume = node_area * np.maximum(source("crustal_thickness_m"), 0.0)
        continental = area * old_hc * (cell_explicit == CRUST_TYPE_CONTINENTAL)
        oceanic = area * old_hc * (cell_explicit != CRUST_TYPE_CONTINENTAL)
        continental += np.bincount(cell, weights=src_volume * (src_explicit == CRUST_TYPE_CONTINENTAL), minlength=count)
        oceanic += np.bincount(cell, weights=src_volume * (src_explicit != CRUST_TYPE_CONTINENTAL), minlength=count)
        winner = np.where(continental > oceanic, CRUST_TYPE_CONTINENTAL, np.where(oceanic > continental, CRUST_TYPE_OCEANIC, cell_explicit))
        fields["crust_type_code"] = np.where(receiving, np.where(winner == own_code, CRUST_TYPE_INHERIT, winner), own["crust_type_code"]).astype(np.int8)

        # Every other field combines the cell's own value (weighted by its area) with the
        # incoming nodes' (weighted by theirs) by remap class -- quad_merge's suture rules.
        for name, spec in SURFACE_FIELDS.items():
            if name in ("elevation", "crust_type_code"):
                continue
            values = source(name)
            if spec.remap_class == RemapClass.EXTENSIVE:
                fields[name] = own[name] + np.bincount(cell, weights=values * node_area, minlength=count) / area
                if name in _CAPS:
                    fields[name] = _spread_overflow(quad, fields[name], own[name], _CAPS[name][1])
            elif spec.remap_class == RemapClass.BOOLEAN_PROVENANCE:
                fields[name] = own[name] | (np.bincount(cell, weights=values.astype(float), minlength=count) > 0.0)
            elif name == "channel_depth" or spec.remap_class == RemapClass.COUNTDOWN:
                np.maximum.at(fields[name], cell, values)
            elif spec.remap_class in (RemapClass.HISTORY, RemapClass.WRITE_ONCE_HISTORY):
                earliest = np.where(own[name] != spec.sentinel, own[name], np.inf)
                set_ = values != spec.sentinel
                np.minimum.at(earliest, cell[set_], values[set_])
                fields[name] = np.where(np.isfinite(earliest), earliest, spec.sentinel)
            elif spec.remap_class == RemapClass.CATEGORICAL:
                fields[name] = _vote_with_own(cell, values, node_area, own[name], area, receiving)
            else:
                # Intensive and clock fields: an area-weighted blend; a concentration or a
                # channel width is weighted by the column or channel it describes.
                coupled = _COUPLED_WEIGHTS.get(name)
                own_weight = area * (own[coupled] if coupled else 1.0)
                incoming_weight = node_area * (source(coupled) if coupled else 1.0)
                total = own_weight + np.bincount(cell, weights=incoming_weight, minlength=count)
                blended = own_weight * own[name] + np.bincount(cell, weights=incoming_weight * values, minlength=count)
                mixed = receiving & (total > 0.0)
                fields[name] = np.where(mixed, blended / np.where(total > 0.0, total, 1.0), own[name])

        # Elevation: isostasy from the stacked column, plus both sides' non-isostatic residual
        # blended by area, as quad_merge does at a suture.
        own_residual = own["elevation"] - lithosphere.isostatic_elevation(old_hc, old_hm, lithosphere.node_crust_density(own["crust_type_code"], quad.crust_type))
        incoming_residual = np.empty(len(cell))
        for k in np.unique(src_plate):
            mine = src_plate == k
            plate_fields = legacy[k].fields
            nodes = src_node[mine]
            density = lithosphere.node_crust_density(plate_fields["crust_type_code"][nodes], legacy[k].crust_type)
            incoming_residual[mine] = plate_fields["elevation"][nodes] - lithosphere.isostatic_elevation(
                plate_fields["crustal_thickness_m"][nodes], plate_fields["mantle_lithosphere_thickness_m"][nodes], density
            )
        residual = (area * own_residual + np.bincount(cell, weights=node_area * incoming_residual, minlength=count)) / (
            area + np.bincount(cell, weights=node_area, minlength=count)
        )
        density = lithosphere.node_crust_density(fields["crust_type_code"], quad.crust_type)
        rebuilt = lithosphere.isostatic_elevation(fields["crustal_thickness_m"], fields["mantle_lithosphere_thickness_m"], density) + residual
        changed = receiving | (fields["crustal_thickness_m"] != old_hc) | (fields["mantle_lithosphere_thickness_m"] != old_hm)
        fields["elevation"] = np.where(changed, rheology.clip_elevation_bounds(rebuilt), own["elevation"])
        quad.set_fields_on_plate(**fields)


# Intensive fields whose blend is weighted by another field: a concentration in the soil by
# the soil's depth, a channel's width by its depth (as quad_merge weights them).
_COUPLED_WEIGHTS = {"soil_mineral_content": "soil_depth", "soil_organic_content": "soil_depth", "channel_width": "channel_depth"}


def _vote_with_own(cell: np.ndarray, values: np.ndarray, weight: np.ndarray, own: np.ndarray, own_weight: np.ndarray, receiving: np.ndarray) -> np.ndarray:
    """A categorical field on cells receiving incoming nodes: the area-weighted vote between
    the cell's own value (weighted by its area) and the incoming nodes'. Ties keep the cell's
    own value, as quad_merge does."""
    choices = np.unique(np.concatenate([values, own[receiving]]))
    count = len(own)
    votes = np.stack([np.bincount(cell, weights=weight * (values == c), minlength=count) + own_weight * (own == c) for c in choices])
    best = votes.max(axis=0)
    own_votes = votes[np.minimum(np.searchsorted(choices, own), len(choices) - 1), np.arange(count)]
    own_votes = np.where(choices[np.minimum(np.searchsorted(choices, own), len(choices) - 1)] == own, own_votes, -np.inf)
    winner = np.where(own_votes >= best * (1.0 - 1e-9), own, choices[np.argmax(votes, axis=0)])
    return np.where(receiving, winner, own).astype(own.dtype)


def _spread_overflow(quad: PlateWithSparseQuadPatch, values: np.ndarray, before: np.ndarray, cap: float) -> np.ndarray:
    """Move the thickness stacking pushed past `cap` (above the cap or the cell's own
    pre-stack value, whichever is higher) outward across the plate's cells, conserving volume:
    each round, every over-full cell hands its excess to its edge neighbours in equal volume
    shares, so it passes through full cells and travels one ring per round, out to
    `OVERFLOW_SPREAD_RINGS`. A stacked suture thickens a belt rather than one column. What is
    still over the cap afterwards is left for `_clamp_columns`."""
    values = values.copy()
    area = quad.node_areas_m2()
    graph = quad.adjacency()
    degree = np.diff(graph.offsets)
    rows = np.repeat(np.arange(len(values)), degree)
    ceiling = np.maximum(cap, before)
    for _ in range(OVERFLOW_SPREAD_RINGS):
        excess = np.where(degree > 0, np.maximum(values - ceiling, 0.0) * area, 0.0)
        if not np.any(excess > 0.0):
            break
        moved = excess[rows] / degree[rows]
        values += (np.bincount(graph.neighbours, weights=moved, minlength=len(values)) - excess) / area
    return values


def _clamp_columns(converted: list[PlateWithSparseQuadPatch]) -> dict[str, float]:
    """Clamp every cell's Hc/Hm into their caps, shifting elevation by the isostatic change so
    the non-isostatic residual survives. Returns the net volume (m^3) removed per field."""
    bounds = _CAPS
    removed = dict.fromkeys(bounds, 0.0)
    for quad in converted:
        if not quad.node_count():
            continue
        area = quad.node_areas_m2()
        old = {name: quad.collect(name) for name in bounds}
        new = {name: np.clip(old[name], *bounds[name]) for name in bounds}
        changed = np.zeros(quad.node_count(), dtype=bool)
        for name in bounds:
            removed[name] += float(np.sum((old[name] - new[name]) * area))
            changed |= new[name] != old[name]
        if not np.any(changed):
            continue
        density = lithosphere.node_crust_density(quad.collect("crust_type_code"), quad.crust_type)
        shift = lithosphere.isostatic_elevation(new["crustal_thickness_m"], new["mantle_lithosphere_thickness_m"], density) - lithosphere.isostatic_elevation(
            old["crustal_thickness_m"], old["mantle_lithosphere_thickness_m"], density
        )
        elevation = quad.collect("elevation")
        elevation[changed] = rheology.clip_elevation_bounds(elevation[changed] + shift[changed])
        quad.set_fields_on_plate(elevation=elevation, **new)
    return removed


def _preserved(legacy: list[LegacyPlate], converted: list[PlateWithSparseQuadPatch], targets: list[tuple[np.ndarray, np.ndarray]]) -> dict[str, Any]:
    """Every node's provenance must be represented in the cell it went to: a volcano node's
    cell is a volcano, an active one's countdown is at least as long, and a node's set
    creation time is no earlier than its cell's. Counts of violations (all must be 0), plus
    totals for context."""
    violations = {"is_volcano": 0, "volcano_active_years_remaining": 0, "node_created_years": 0, "untargeted": 0}
    for k, plate in enumerate(legacy):
        t_plate, t_cell = targets[k]
        violations["untargeted"] += int(np.count_nonzero(t_plate < 0))
        for t in np.unique(t_plate[t_plate >= 0]):
            at = t_plate == t
            cell = t_cell[at]
            quad = converted[t]
            volcano = plate.fields["is_volcano"][at]
            violations["is_volcano"] += int(np.count_nonzero(volcano & ~quad.collect("is_volcano")[cell]))
            countdown = plate.fields["volcano_active_years_remaining"][at]
            violations["volcano_active_years_remaining"] += int(np.count_nonzero(quad.collect("volcano_active_years_remaining")[cell] < countdown))
            created = plate.fields["node_created_years"][at]
            cell_created = quad.collect("node_created_years")[cell]
            set_ = created >= 0.0
            violations["node_created_years"] += int(np.count_nonzero(set_ & ((cell_created < 0.0) | (cell_created > created))))

    def before(name: str) -> np.ndarray:
        return np.concatenate([p.fields[name] for p in legacy]) if legacy else np.zeros(0)

    def after(name: str) -> np.ndarray:
        chunks = [q.collect(name) for q in converted if q.node_count()]
        return np.concatenate(chunks) if chunks else np.zeros(0)

    return {
        "provenance_violations": violations,
        "volcano_nodes_before": int(np.count_nonzero(before("is_volcano"))),
        "volcano_cells_after": int(np.count_nonzero(after("is_volcano"))),
        "elevation_range_before": [float(before("elevation").min(initial=0.0)), float(before("elevation").max(initial=0.0))],
        "elevation_range_after": [float(after("elevation").min(initial=0.0)), float(after("elevation").max(initial=0.0))],
    }


def _node_components(world_xyz: np.ndarray, connect_rad: float) -> int:
    if len(world_xyz) == 0:
        return 0
    pairs = cKDTree(world_xyz).query_pairs(_chord(connect_rad), output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(world_xyz),) * 2)
    return int(connected_components(graph, directed=False)[0])


def _cell_components(plate: PlateWithSparseQuadPatch) -> int:
    count = plate.node_count()
    if count == 0:
        return 0
    graph = plate.adjacency()
    rows = np.repeat(np.arange(count), np.diff(graph.offsets))
    matrix = coo_matrix((np.ones(len(rows)), (rows, graph.neighbours)), shape=(count, count))
    return int(connected_components(matrix, directed=False)[0])


def _fibonacci_sphere(count: int) -> np.ndarray:
    k = np.arange(count) + 0.5
    z = 1.0 - 2.0 * k / count
    r = np.sqrt(1.0 - z * z)
    lon = np.pi * (1.0 + 5**0.5) * k
    return np.column_stack([r * np.cos(lon), r * np.sin(lon), z])


def _coverage(plates: list, count: int) -> dict[str, float]:
    samples = _fibonacci_sphere(count)
    claims = np.zeros(count, dtype=np.int64)
    for plate in plates:
        if plate.node_count():
            claims += np.asarray(plate.contains_batch(samples), dtype=bool)
    return {"uncovered": float(np.mean(claims == 0)), "multiply_covered": float(np.mean(claims > 1))}


def _quad_coverage(plates: list[PlateWithSparseQuadPatch], count: int) -> dict[str, float]:
    return _coverage(plates, count) if count else {}


def _reset_world_caches(world: "World") -> None:
    """Drop every world cache keyed by node order or node count. Climate is an Eulerian grid,
    but its snapshot was sampled from the line nodes; recomputing it is cheap."""
    for name in (
        "climate_cache",
        "hydrology_cache",
        "hydrology_cache_step",
        "erosion_cache",
        "land_kdtree_cache",
        "node_kdtree_cache",
        "node_position_tree_cache",
        "node_healpix_grid_cache",
        "node_healpix_index_cache",
        "node_kdtree_relief_cache",
        "node_hillshade_cache",
        "collision_polarity_frame",
        "boundary_search_cache",
    ):
        if hasattr(world, name):
            setattr(world, name, None)
