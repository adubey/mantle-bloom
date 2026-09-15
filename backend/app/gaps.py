"""Whole-sphere coverage maintenance: grow the plate(s) adjacent to any region no plate
currently covers into it, falling back to spawning new crust -- oceanic almost everywhere,
continental only where a gap point genuinely borders a still-standing continental coastline
(see GAP_LAND_ADOPTION_RADIUS_MULT) -- only when nothing is adjacent.

`LithospherePlate.deform()`'s per-step boundary growth only ever extends a line from an
*existing* node -- a plate can spread into space right next to its own current edge, but
nothing grows crust somewhere no plate has any nearby line at all. Most of the time that's
fine: a newly-opened divergent gap is one or two nodes wide and next step's ordinary growth
closes it. But once every oceanic plate bordering a stretch of open ocean has been ground
down by subduction and fully removed (see `merge_split.remove_defunct_plates`), the sphere
area it used to occupy has no plate left anywhere near it -- there is nothing there to grow,
so it just stays empty forever. Confirmed on a 399 My / node_density=4 save (seed 920135003,
see GitHub issue #126's "Very-long-run collapse" section): ~42% of the sphere had zero elevation
nodes, all of it sphere area no live plate's lines reached.

This module finds those genuinely-uncovered regions periodically (same cadence as
`merge_split.defragment_plates` -- a whole-world k-d-tree pass, cheap but not free) and, for
each big-enough one, grows the plate(s) genuinely adjacent to it into the gap node by node
(`fill_gaps_by_growing_neighbours`, via `gap_fill_frontier.fill_gap_by_growing_plates`) --
falling back to spawning a brand-new neutral plate (`_spawn_plate_from_gap`) only when nothing
is adjacent at all (a fully-vacated region with no live plate left nearby to grow).

The new plate's own composition (spawn fallback) is decided per node, not blanket-oceanic: real
new crust in open water is oceanic (the same crust type any mid-ocean ridge produces), but a gap
point right at a still-standing continental coastline -- e.g. a fully-subducted marginal sea
landlocked by continent -- comes back continental instead (see `_spawn_plate_from_gap`'s own
`node_is_continental`). The spawned plate's own `crust_type` label is the majority of what it
actually ended up with (`elevation_lines.majority_crust_type`), so it is oceanic in practice
for all but that rare landlocked case.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import gap_fill_frontier, geometry, mantle
from .elevation_lines import (
    COVERAGE_RADIUS_MULT as _SHARED_COVERAGE_RADIUS_MULT,
    effective_is_continental_from_codes,
    iter_local_lattice,
    line_spacing_rad,
)
from .lithosphere_plate import LithospherePlate, new_plate

if TYPE_CHECKING:
    from .world import World

# Cadence: `fill_gaps_by_growing_neighbours` is a whole-sphere lattice sweep (O(nodes) at full
# density), cheap but not free, and coverage doesn't collapse fast -- same reasoning and same
# cadence as merge_split.DEFRAG_INTERVAL_STEPS, which world.step_world calls it alongside.
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

# `fill_gaps_by_growing_neighbours`'s own "detect adjacent plates" gate: a plate with at least
# one node within this multiple of spacing of a gap cluster is genuinely adjacent to it, not
# just the nearest thing on an otherwise-empty sphere. Wider than CLUSTER_RADIUS_MULT (which
# only has to bridge one lattice step to the *next* gap point) since this has to reach a real
# plate's actual edge, which -- unlike the gap lattice itself -- was never guaranteed to sit
# right at the cluster boundary; narrower than CORNER_NOTCH_NEIGHBOUR_REACH_MARGIN_MULT-scale
# reaches (lithosphere_plate.py), which are about a single step's worst-case motion rather than
# "is there a plate here at all."
ADJACENT_PLATE_REACH_MULT = 6.0

# Sweeping in the identity frame reuses iter_local_lattice as a plain global lat/lon
# lattice -- fine for a one-off detection query even though it has the usual pole bias,
# since nothing here is carried forward as persistent state.
_GLOBAL_FRAME = np.eye(3)

# Gap-age tracking (see reconcile_gap_tracks/GapTrack): a gap cluster has no stable identity
# across steps the way a plate id does, so -- same problem stranded_basins.py already solved
# for endorheic basins -- persistence comes from matching this step's cluster centroids
# against last step's by proximity, not by any persistent key. Tighter than
# stranded_basins.MATCH_DISTANCE_RAD (0.08 rad): a gap cluster is typically a thin sliver along
# a boundary rather than a basin's own compact catchment, so a looser gate risks fusing two
# genuinely distinct nearby gaps (e.g. both sides of a triple junction) into one track.
GAP_MATCH_DISTANCE_RAD = 0.05
# Deliberately much lower than MIN_GAP_NODES (which gates a real plate *spawn*, an expensive
# and disruptive act) -- age-tracking's whole point is to surface exactly the small, easy-to-
# miss persistent notches (a stuck triple-junction corner, a thin sliver along a rift) that
# never reach MIN_GAP_NODES and so never trigger a spawn. 1 is the floor, not a tuned value:
# tracking must never hide a gap a spawn-gate would still (eventually) act on.
GAP_AGE_MIN_CLUSTER_NODES = 1

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


def _adjacent_plates_to_cluster(world: "World", cluster_points: np.ndarray, spacing_rad: float) -> list[LithospherePlate]:
    """Live plates with at least one node within `ADJACENT_PLATE_REACH_MULT * spacing_rad` of
    `cluster_points` -- "detect adjacent plates" for `fill_gaps_by_growing_neighbours`. Can
    come back empty (a fully-vacated region with no live plate left nearby), which its caller
    falls back to spawning a plate for."""
    reach_rad = ADJACENT_PLATE_REACH_MULT * spacing_rad
    adjacent = []
    for plate in world.plates:
        points, _ = plate.all_points_and_elevation()
        if len(points) == 0:
            continue
        dist, _ = cKDTree(points).query(cluster_points, k=1, distance_upper_bound=reach_rad)
        if np.any(np.isfinite(dist)):
            adjacent.append(plate)
    return adjacent


def fill_gaps_by_growing_neighbours(world: "World") -> list[str]:
    """Find every sphere region no live plate currently covers and, for each one at least
    `MIN_GAP_NODES` (scaled by `world.node_density`) large, first look for plate(s) actually
    adjacent to it (`_adjacent_plates_to_cluster`) and, when there are any, grow those
    *existing* plates into the cluster node by node (`gap_fill_frontier.fill_gap_by_growing_
    plates` -- see that module's own docstring). Falls back to `_spawn_plate_from_gap` only
    when a cluster has no adjacent plate at all -- a fully-vacated region with nothing nearby
    to grow. Mutates `world.plates`/`world.next_plate_id` (spawn fallback) or existing plates'
    own lines (grow path) in place; returns event strings for the UI's console."""
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
        adjacent = _adjacent_plates_to_cluster(world, cluster_points, spacing_rad)
        if not adjacent:
            plate = _spawn_plate_from_gap(world, cluster_points, spacing_rad, existing_context)
            if plate.node_count() == 0:
                continue
            world.plates.append(plate)
            where = "in open water no plate had reached in a long time" if plate.crust_type == "oceanic" else "over a long-vacated, landlocked gap"
            events.append(
                f"New {plate.crust_type} crust formed as plate {plate.plate_id} ({plate.node_count()} nodes) {where} (no adjacent plate to grow)."
            )
            continue
        added = gap_fill_frontier.fill_gap_by_growing_plates(world, cluster_points, adjacent, spacing_rad)
        for plate_id, n in added.items():
            events.append(f"Plate {plate_id} grew by {n} nodes into a long-vacated gap.")
    return events


@dataclass
class GapTrack:
    """One entry in `world.gap_tracks`: the minimal cross-step memory a still-uncovered
    lattice cluster needs, reconciled by centroid proximity each step -- same shape and same
    reasoning as `stranded_basins.StrandedBasinTrack` (a gap cluster has no persistent identity
    across steps any more than a stranded basin does). Diagnostic only, surfaced by the
    `overlapAge` debug render view's gap-age layer (see docs/debugging.md); nothing in the
    physics reads it back."""

    centroid_xyz: np.ndarray  # (3,) unit vector
    first_seen_years: float
    last_seen_years: float
    steps_seen: int
    node_count: int


def _nearest_gap_track(
    prev_centroids: np.ndarray | None, used: np.ndarray, centroid: np.ndarray
) -> int | None:
    """Index of the closest not-yet-claimed previous track within `GAP_MATCH_DISTANCE_RAD` of
    `centroid`, or `None` -- identical shape to `stranded_basins._nearest_track`."""
    if prev_centroids is None or len(prev_centroids) == 0:
        return None
    dist = geometry.angular_distance(prev_centroids, centroid)
    dist = np.where(used, np.inf, dist)
    j = int(np.argmin(dist))
    return j if dist[j] <= GAP_MATCH_DISTANCE_RAD else None


def reconcile_gap_tracks(world: "World") -> None:
    """Recompute this step's uncovered-lattice clusters (the same whole-sphere sweep
    `fill_gaps_by_growing_neighbours` uses, but with no `MIN_GAP_NODES` floor -- see
    `GAP_AGE_MIN_CLUSTER_NODES`'s own comment on why age-tracking needs a much lower one) and
    reconcile `world.gap_tracks` by centroid proximity: a cluster matching a previous track
    keeps its `first_seen_years` and bumps `steps_seen`; an unmatched cluster starts a fresh
    track; a track with no matching cluster this step is dropped by omission. Same "replace
    wholesale" pattern as `stranded_basins.reconcile_world_tracks`. Called at the same cadence
    as `fill_gaps_by_growing_neighbours` (`GAP_FILL_INTERVAL_STEPS`) from `world.step_world`,
    right alongside it."""
    existing_context = _existing_node_tree(world)
    if existing_context is None:
        world.gap_tracks = []
        return

    spacing_rad = line_spacing_rad(world.node_density)
    gap_points = _find_gap_points(existing_context, spacing_rad)
    if len(gap_points) == 0:
        world.gap_tracks = []
        return

    labels = _cluster(gap_points, CLUSTER_RADIUS_MULT * spacing_rad)
    prev_tracks = world.gap_tracks
    prev_centroids = np.array([t.centroid_xyz for t in prev_tracks]) if prev_tracks else None
    used = np.zeros(len(prev_tracks), dtype=bool)

    new_tracks: list[GapTrack] = []
    for label in np.unique(labels):
        cluster_points = gap_points[labels == label]
        if len(cluster_points) < GAP_AGE_MIN_CLUSTER_NODES:
            continue
        centroid = geometry.normalize(cluster_points.mean(axis=0))
        match = _nearest_gap_track(prev_centroids, used, centroid)
        if match is None:
            new_tracks.append(GapTrack(centroid, world.elapsed_years, world.elapsed_years, 1, len(cluster_points)))
        else:
            used[match] = True
            old = prev_tracks[match]
            new_tracks.append(
                GapTrack(centroid, old.first_seen_years, world.elapsed_years, old.steps_seen + 1, len(cluster_points))
            )
    world.gap_tracks = new_tracks
