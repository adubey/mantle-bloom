"""Whole-sphere coverage maintenance: spawn new oceanic crust into any region no plate
currently covers.

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
each big-enough one with a brand-new oceanic plate (new crust forming in open water, the
same crust type any mid-ocean ridge produces). It deliberately does *not* try to instead grow
an existing neighbouring plate into the gap -- besides needing a partition/absorption scheme
of its own (the pre-refactor `gaps.py` this replaces did that too, "if bordered mainly by one
plate, absorb it into that plate," never ported to this engine -- see docs/TODO.md), handing
a large freshly-vacated region to whichever plate happens to be nearest would feed exactly
the continental-growth ratchet already tracked there. A brand-new plate is neutral: it can
still merge, subduct, or get absorbed by ordinary boundary growth like any other plate once
it has a real neighbour again.

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
from .elevation_lines import iter_local_lattice, line_spacing_rad
from .lithosphere_plate import LithospherePlate, new_plate

if TYPE_CHECKING:
    from .world import World

# Cadence: `fill_gaps` is a whole-sphere lattice sweep (O(nodes) at full density), cheap but
# not free, and coverage doesn't collapse fast -- same reasoning and same cadence as
# merge_split.DEFRAG_INTERVAL_STEPS, which world.step_world calls it alongside.
GAP_FILL_INTERVAL_STEPS = 4

# A node is "covered" if some existing plate has a node within this of it -- comfortably
# more than one line spacing so ordinary per-step catch-up growth (deform() closing a
# freshly-opened one-node-wide divergent gap) is never mistaken for a genuine void.
COVERAGE_RADIUS_MULT = 1.5
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


def _existing_node_tree(world: "World") -> cKDTree | None:
    chunks = []
    for plate in world.plates:
        pts, _ = plate.all_points_and_elevation()
        if len(pts) > 0:
            chunks.append(pts)
    if not chunks:
        return None
    return cKDTree(np.concatenate(chunks, axis=0))


def _find_gap_points(existing_tree: cKDTree, spacing_rad: float) -> np.ndarray:
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad
    chunks = []
    for _, _, world_pts in iter_local_lattice(_GLOBAL_FRAME, spacing_rad=spacing_rad):
        dist, _ = existing_tree.query(world_pts)
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


def _spawn_plate_from_gap(world: "World", cluster_points: np.ndarray, spacing_rad: float) -> LithospherePlate:
    centroid = geometry.normalize(cluster_points.mean(axis=0))
    frame = geometry.plate_frame_from_seed(centroid)
    cluster_tree = cKDTree(cluster_points)
    coverage_radius_rad = COVERAGE_RADIUS_MULT * spacing_rad

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = cluster_tree.query(world_pts)
        return dist < coverage_radius_rad

    plate = new_plate(world.next_plate_id, frame, "oceanic", spacing_rad, world.seed, is_owned=is_owned)
    world.next_plate_id += 1

    points, _ = plate.all_points_and_elevation()
    if len(points) > 0:
        velocities = mantle.flow_at(points, world.mantle_centers)
        plate.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(points, velocities)))
    return plate


def fill_gaps(world: "World") -> list[str]:
    """Find every sphere region no live plate currently covers and, for each one at least
    `MIN_GAP_NODES` (scaled by `world.node_density`) large, spawn a new oceanic plate to
    cover it. Mutates `world.plates`/`world.next_plate_id` in place; returns event strings
    for the UI's console."""
    existing_tree = _existing_node_tree(world)
    if existing_tree is None:
        return []

    spacing_rad = line_spacing_rad(world.node_density)
    gap_points = _find_gap_points(existing_tree, spacing_rad)
    if len(gap_points) == 0:
        return []

    labels = _cluster(gap_points, CLUSTER_RADIUS_MULT * spacing_rad)
    min_gap_nodes = max(1, round(MIN_GAP_NODES * world.node_density))

    events: list[str] = []
    for label in np.unique(labels):
        cluster_points = gap_points[labels == label]
        if len(cluster_points) < min_gap_nodes:
            continue
        plate = _spawn_plate_from_gap(world, cluster_points, spacing_rad)
        if plate.node_count() == 0:
            continue
        world.plates.append(plate)
        events.append(
            f"New oceanic crust formed as plate {plate.plate_id} ({plate.node_count()} nodes) "
            "in open water no plate had reached in a long time."
        )
    return events
