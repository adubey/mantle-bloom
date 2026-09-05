"""Whole-sphere coverage maintenance: spawn new crust into any region no plate currently
covers -- oceanic almost everywhere, continental only where a gap point genuinely borders a
still-standing continental coastline (see GAP_LAND_ADOPTION_RADIUS_MULT).

`LithospherePlate.deform()`'s per-step boundary growth only ever extends a line from an
*existing* node -- a plate can spread into space right next to its own current edge, but
nothing grows crust somewhere no plate has any nearby line at all. Most of the time that's
fine: a newly-opened divergent gap is one or two nodes wide and next step's ordinary growth
closes it. But once every oceanic plate bordering a stretch of open ocean has been ground
down by subduction and fully removed (see `merge_split.remove_defunct_plates`), the sphere
area it used to occupy has no plate left anywhere near it -- there is nothing there to grow,
so it just stays empty forever. Confirmed on a 399 My / node_density=4 save (seed 920135003,
see docs/TODO.md's "Very-long-run collapse" section): ~42% of the sphere had zero elevation
nodes, all of it sphere area no live plate's lines reached.

This module finds those genuinely-uncovered regions periodically (same cadence as
`merge_split.defragment_plates` -- a whole-world k-d-tree pass, cheap but not free) and fills
each big-enough one with a brand-new plate. It deliberately does *not* try to instead grow an
existing neighbouring plate into the gap -- besides needing a partition/absorption scheme of
its own (the pre-refactor `gaps.py` this replaces did that too, "if bordered mainly by one
plate, absorb it into that plate," never ported to this engine -- see docs/TODO.md), handing
a large freshly-vacated region to whichever plate happens to be nearest would feed exactly
the continental-growth ratchet already tracked there. A brand-new plate is neutral: it can
still merge, subduct, or get absorbed by ordinary boundary growth like any other plate once
it has a real neighbour again.

The new plate's own composition is decided per node, not blanket-oceanic: real new crust in
open water is oceanic (the same crust type any mid-ocean ridge produces), but a gap point
right at a still-standing continental coastline -- e.g. a fully-subducted marginal sea
landlocked by continent -- comes back continental instead (see `_spawn_plate_from_gap`'s own
`node_is_continental`). The spawned plate's own `crust_type` label is the majority of what it
actually ended up with (`elevation_lines.majority_crust_type`), so it is oceanic in practice
for all but that rare landlocked case.

Known stopgap, not the real fix -- see docs/TODO.md ("`gaps.py`'s plate-spawn is a stopgap,
not the real fix"): conjuring a whole fully-formed plate into existence after the fact isn't
how new ocean floor actually forms (continuous mid-ocean-ridge spreading off an existing
plate's own divergent edge is). The real fix is upstream, in `deform()`'s own per-step
growth; this module should shrink back to a rare fallback once that exists, not stay the
primary mechanism.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import geometry, mantle
from .elevation_lines import (
    COVERAGE_RADIUS_MULT as _SHARED_COVERAGE_RADIUS_MULT,
    effective_is_continental_from_codes,
    iter_local_lattice,
    line_spacing_rad,
)
from .lithosphere_plate import LithospherePlate, new_plate

if TYPE_CHECKING:
    from .world import World

# Cadence: `fill_gaps` is a whole-sphere lattice sweep (O(nodes) at full density), cheap but
# not free, and coverage doesn't collapse fast -- same reasoning and same cadence as
# merge_split.DEFRAG_INTERVAL_STEPS, which world.step_world calls it alongside.
GAP_FILL_INTERVAL_STEPS = 4

# "Covered" (see elevation_lines.COVERAGE_RADIUS_MULT, shared with LithospherePlate's own
# local divergent-boundary growth so both agree on the same tolerance) -- re-exported under
# this module's own name since callers/tests already refer to it as gaps.COVERAGE_RADIUS_MULT.
COVERAGE_RADIUS_MULT = _SHARED_COVERAGE_RADIUS_MULT
# Two uncovered lattice points within this of each other belong to the same gap cluster --
# a bit looser than COVERAGE_RADIUS_MULT so one contiguous void isn't sliced into several
# clusters at its own lattice resolution.
CLUSTER_RADIUS_MULT = 2.0
# A node count, not a distance -- scales with node_density directly, same reasoning as
# merge_split.SPLIT_MIN_NODES/DEFRAG_FRAGMENT_MIN_NODES. Deliberately in the same range as
# SPLIT_MIN_NODES (a split's own minimum daughter size): anything smaller than "big enough to
# be its own plate" is left alone as ordinary boundary-growth catch-up lag rather than
# spawning a sliver plate at every busy divergent boundary every interval.
MIN_GAP_NODES = 500

# Sweeping in the identity frame reuses iter_local_lattice as a plain global lat/lon
# lattice -- fine for a one-off detection query even though it has the usual pole bias,
# since nothing here is carried forward as persistent state.
_GLOBAL_FRAME = np.eye(3)

# A gap point adopts the *continental* type only if the nearest pre-existing node is itself
# continental, still above sea level, and within this many line-spacings -- hugging a real
# coastline (e.g. a fully-subducted marginal sea landlocked by continent), not reaching all
# the way across an ocean basin to a far-off continent. Every other gap point (the
# overwhelming majority -- gaps are, per this module's own docstring, almost always open
# ocean a fully-subducted plate vacated) stays oceanic, exactly as before this field existed.
GAP_LAND_ADOPTION_RADIUS_MULT = 3.0


class _ExistingNodeContext:
    """Every currently-live node's position plus enough context to decide a newly-upwelled
    gap point's own crust type by what actually borders it -- see `_spawn_plate_from_gap`."""

    def __init__(self, tree: cKDTree, is_continental: np.ndarray, elevation: np.ndarray) -> None:
        self.tree = tree
        self.is_continental = is_continental
        self.elevation = elevation


def _existing_node_tree(world: "World") -> _ExistingNodeContext | None:
    point_chunks = []
    continental_chunks = []
    elevation_chunks = []
    for plate in world.plates:
        pts, elev = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        point_chunks.append(pts)
        elevation_chunks.append(elev)
        continental_chunks.append(effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental"))
    if not point_chunks:
        return None
    return _ExistingNodeContext(
        cKDTree(np.concatenate(point_chunks, axis=0)),
        np.concatenate(continental_chunks, axis=0),
        np.concatenate(elevation_chunks, axis=0),
    )


def _find_gap_points(existing_tree: _ExistingNodeContext, spacing_rad: float) -> np.ndarray:
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad
    chunks = []
    for _, _, world_pts in iter_local_lattice(_GLOBAL_FRAME, spacing_rad=spacing_rad):
        dist, _ = existing_tree.tree.query(world_pts)
        uncovered = dist > coverage_radius_rad
        if np.any(uncovered):
            chunks.append(world_pts[uncovered])
    return np.concatenate(chunks, axis=0) if chunks else np.zeros((0, 3))


def _cluster(points: np.ndarray, radius_rad: float) -> np.ndarray:
    """Connected-component label per point: points within `radius_rad` of each other
    (directly or transitively) share a label -- same technique plates.node_components uses."""
    n = len(points)
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = cKDTree(points).query_pairs(r=radius_rad, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def _spawn_plate_from_gap(
    world: "World", cluster_points: np.ndarray, spacing_rad: float, existing_context: _ExistingNodeContext
) -> LithospherePlate:
    centroid = geometry.normalize(cluster_points.mean(axis=0))
    frame = geometry.plate_frame_from_seed(centroid)
    cluster_tree = cKDTree(cluster_points)
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad
    land_adoption_radius_rad = GAP_LAND_ADOPTION_RADIUS_MULT * spacing_rad

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = cluster_tree.query(world_pts)
        return dist < coverage_radius_rad

    def node_is_continental(world_pts: np.ndarray) -> np.ndarray:
        # A new gap node adopts continental type only where it's genuinely hugging a real,
        # still-standing coastline -- the nearest *pre-existing* node (not the nearest gap
        # point) is itself continental, still above sea level, and close by. Everywhere else
        # in the gap -- the overwhelming majority of one, per this module's own docstring --
        # stays oceanic, matching the behaviour before this rule existed.
        dist, idx = existing_context.tree.query(world_pts)
        borders_land = existing_context.is_continental[idx] & (existing_context.elevation[idx] > 0.0)
        return borders_land & (dist <= land_adoption_radius_rad)

    plate = new_plate(
        world.next_plate_id, frame, "oceanic", spacing_rad, world.seed, is_owned=is_owned, node_is_continental=node_is_continental
    )
    world.next_plate_id += 1

    points, _ = plate.all_points_and_elevation()
    if len(points) > 0:
        velocities = mantle.flow_at(points, world.mantle_centers)
        plate.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(points, velocities)))
    return plate


def fill_gaps(world: "World") -> list[str]:
    """Find every sphere region no live plate currently covers and, for each one at least
    `MIN_GAP_NODES` (scaled by `world.node_density`) large, spawn a new plate to cover it --
    oceanic almost everywhere (real gaps are overwhelmingly open water a fully-subducted
    plate vacated), except nodes genuinely hugging a still-standing continental coastline
    (see GAP_LAND_ADOPTION_RADIUS_MULT), which come back continental. Mutates
    `world.plates`/`world.next_plate_id` in place; returns event strings for the UI's
    console."""
    existing_context = _existing_node_tree(world)
    if existing_context is None:
        return []

    spacing_rad = line_spacing_rad(world.node_density)
    gap_points = _find_gap_points(existing_context, spacing_rad)
    if len(gap_points) == 0:
        return []

    labels = _cluster(gap_points, CLUSTER_RADIUS_MULT * spacing_rad)
    min_gap_nodes = max(1, round(MIN_GAP_NODES * world.node_density))

    events: list[str] = []
    for label in np.unique(labels):
        cluster_points = gap_points[labels == label]
        if len(cluster_points) < min_gap_nodes:
            continue
        plate = _spawn_plate_from_gap(world, cluster_points, spacing_rad, existing_context)
        if plate.node_count() == 0:
            continue
        world.plates.append(plate)
        where = "in open water no plate had reached in a long time" if plate.crust_type == "oceanic" else "over a long-vacated, landlocked gap"
        events.append(
            f"New {plate.crust_type} crust formed as plate {plate.plate_id} ({plate.node_count()} nodes) {where}."
        )
    return events
