"""Whole-sphere coverage maintenance.

Per-step boundary evolution (boundary.py) only ever grows/shrinks a line's two theta-
direction ends, so a plate can only spread in the direction its existing rows already run
-- it never gains a whole new row toward its own local pole, and territory vacated by a
fully-subducted plate isn't automatically reclaimed by its neighbors (see plates.py's
docstring and docs/simulation-model.md for the underlying per-line design). Both leave
literal gaps: sphere area no plate currently covers.

This module finds those gaps periodically (same cadence as line_regrid's line
regularization, see world.step_world) and resolves each one of two ways depending on whether
one plate dominates its border:

- One plate accounts for most of the gap's border (DOMINANT_BORDER_FRACTION), or, failing
  that, a *young* plate (age_steps <= YOUNG_PLATE_AGE_STEPS) has a meaningful share of it
  (YOUNG_PLATE_MIN_BORDER_FRACTION): absorbed into that plate. The plate's line set is
  rebuilt from its own local lattice (plates.build_lines_from_lattice), preserving elevation
  wherever old data exists (nearest-neighbor lookup, the same one-time-resample technique
  merge_split.merge_plates uses for merging two plates) and giving the newly-claimed area
  fresh ridge/rift elevation. This is what actually lets a plate grow toward its own pole,
  or reclaim a neighbor's vacated territory -- i.e. "add points to the polygon as it gets
  bigger." Both the ring a single pass can claim (GROWTH_RING_RAD) and the total a plate can
  absorb in one call (MAX_ABSORB_NODES_PER_PLATE_PER_CALL) are capped, so this stays gradual
  like every other change in the model rather than able to jump a plate's size in one step.
- No plate dominates and no young plate has a meaningful share (or the border is empty):
  spawned as a brand new plate -- new crust forming in the open space between separating
  plates, rather than being arbitrarily assigned to one side.

Gap clusters smaller than MIN_GAP_POINTS are left alone -- ordinary one-step growth lag
that boundary.py's per-line extension will already be closing on its own; reacting to those
here would just mean needlessly rebuilding a plate's lines every pass, or worse, repeatedly
spawning a new sliver plate at the same busy boundary (see MIN_GAP_POINTS' comment below).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import geometry, line_regrid, mantle
from .noise import SphereNoise
from .plates import (
    TARGET_LINE_SPACING_RAD,
    Plate,
    base_elevation,
    build_lines_from_lattice,
    iter_local_lattice,
    noise_amplitude,
)

if TYPE_CHECKING:
    from .world import World

COVERAGE_RADIUS_RAD = 1.2 * TARGET_LINE_SPACING_RAD
BORDER_RADIUS_RAD = 2.5 * TARGET_LINE_SPACING_RAD
CLUSTER_RADIUS_RAD = 1.5 * TARGET_LINE_SPACING_RAD
# How far into a gap a single absorb pass reaches -- see _absorb_gap_into_plate.
GROWTH_RING_RAD = 2.0 * TARGET_LINE_SPACING_RAD
# Ordinary per-line theta growth (boundary.py) already closes most divergent-boundary gaps
# on its own, one node per line per step; what's left over after a GC interval is normally
# just catch-up lag, not a genuinely new region. MIN_GAP_POINTS needs to sit well above that
# routine lag or this reacts to it every pass -- an early version used 6 and, at a busy
# divergent boundary, kept spawning a new sliver plate next to the last one it had just
# spawned, every GC interval, fragmenting a boundary into dozens of near-parallel micro-
# plates within a few hundred Myr (the same failure shape merge_split.py's split threshold
# had before it got a similar tightening -- see the note there).
# A point count over a 2D area, not a distance -- doesn't scale automatically with
# TARGET_LINE_SPACING_RAD the way the distance-based radii above do. Point count for a
# given physical area scales with the *square* of resolution, so this is scaled by ~4x
# alongside plates.py's most recent resolution doubling, to keep meaning "the same physical
# gap size" rather than silently reacting to 4x smaller gaps than before.
MIN_GAP_POINTS = 160
# A gap cluster's border is rarely split cleanly between plates -- usually it's mostly one
# plate's edge with a few incidentally-close points from a neighbor. Only treat it as a
# genuine new rift (spawn a new plate) when no single bordering plate dominates; otherwise
# absorb it into whichever plate borders most of it.
DOMINANT_BORDER_FRACTION = 0.75
# A plate with an already-large perimeter tends to dominate the border of *every* nearby
# gap (see _dominant_border_plate), so GROWTH_RING_RAD alone doesn't stop it steadily
# outgrowing everyone else -- it just means the growth arrives over more calls instead of
# one. This caps how many nodes any single plate can absorb in one fill_gaps call in total,
# across every gap it dominates that call, so gap-filling stays a modest coverage-cleanup
# mechanism rather than becoming the dominant growth channel for whichever plate is already
# biggest (ordinary per-step boundary evolution, boundary.py, is still where most growth
# should come from).
# Also an area-based count (see MIN_GAP_POINTS above for why that means ~4x, not the same
# factor as the 1D distance-based constants).
MAX_ABSORB_NODES_PER_PLATE_PER_CALL = 160
# A wide, long-lived rift with no dominant side would, left alone, spawn a brand new sliver
# plate at *every* gap-fill pass -- a fresh spawn never dominates its own border any better
# than the last one did, since neither side of a genuinely symmetric spread is "winning".
# Giving a recently-created plate first claim on its own neighborhood (even without an
# outright majority) breaks that chain: it gets a few passes to consolidate the rift it was
# born into before being treated as just another equal competitor.
YOUNG_PLATE_AGE_STEPS = 4 * line_regrid.REGULARIZE_INTERVAL_STEPS
YOUNG_PLATE_MIN_BORDER_FRACTION = 0.3

# Sweeping in the identity frame reuses plates.iter_local_lattice as a plain global lat/lon
# lattice -- fine for a one-off detection query even though it has the usual pole bias,
# since nothing here is carried forward as persistent state.
_GLOBAL_FRAME = np.eye(3)


def _all_existing_points(world: "World") -> tuple[np.ndarray, np.ndarray]:
    """Every plate's current elevation-node positions, concatenated, alongside a
    same-length array of which plate_id owns each one."""
    points_list = []
    owner_list = []
    for plate in world.plates:
        pts, _ = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        points_list.append(pts)
        owner_list.append(np.full(len(pts), plate.plate_id))
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    return np.concatenate(points_list, axis=0), np.concatenate(owner_list, axis=0)


def _find_gap_points(existing_tree: cKDTree) -> np.ndarray:
    gap_chunks = []
    for _, _, world_pts in iter_local_lattice(_GLOBAL_FRAME):
        dist, _ = existing_tree.query(world_pts)
        uncovered = dist > COVERAGE_RADIUS_RAD
        if np.any(uncovered):
            gap_chunks.append(world_pts[uncovered])
    if not gap_chunks:
        return np.zeros((0, 3))
    return np.concatenate(gap_chunks, axis=0)


def _cluster(points: np.ndarray, radius: float) -> np.ndarray:
    """Connected-component label per point: points within `radius` of each other (directly
    or transitively) share a label."""
    n = len(points)
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = cKDTree(points).query_pairs(r=radius, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def _preferred_border_plate(
    cluster_points: np.ndarray,
    existing_tree: cKDTree,
    existing_owner: np.ndarray,
    plate_by_id: dict[int, Plate],
) -> int | None:
    """The plate id this gap cluster should be absorbed into, if any -- otherwise None,
    meaning the gap is genuinely shared and should become a new plate.

    First choice: a plate accounting for at least DOMINANT_BORDER_FRACTION of the nearby
    border. Failing that, a *young* plate (age_steps <= YOUNG_PLATE_AGE_STEPS) with at
    least YOUNG_PLATE_MIN_BORDER_FRACTION -- see YOUNG_PLATE_AGE_STEPS for why that
    exception exists."""
    dist, idx = existing_tree.query(cluster_points)
    close = dist < BORDER_RADIUS_RAD
    owners = existing_owner[idx[close]]
    if len(owners) == 0:
        return None
    ids, counts = np.unique(owners, return_counts=True)
    total = len(owners)
    dominant = np.argmax(counts)
    if counts[dominant] / total >= DOMINANT_BORDER_FRACTION:
        return int(ids[dominant])

    for i in np.argsort(-counts):
        pid = int(ids[i])
        fraction = counts[i] / total
        if fraction < YOUNG_PLATE_MIN_BORDER_FRACTION:
            break
        plate = plate_by_id.get(pid)
        if plate is not None and plate.age_steps <= YOUNG_PLATE_AGE_STEPS:
            return pid
    return None


def _absorb_gap_into_plate(
    plate: Plate, gap_points: np.ndarray, rng: np.random.Generator, max_new_points: int
) -> int:
    """Grow `plate` into (some of) `gap_points`, claiming at most `max_new_points` of them
    (the ones closest to the plate's existing territory first) -- see
    MAX_ABSORB_NODES_PER_PLATE_PER_CALL for why that cap exists. Also restricted to
    GROWTH_RING_RAD of the plate's *existing* nodes -- one ring of new territory, not the
    whole (possibly huge, e.g. an entire uncovered polar cap nobody's lines reach yet)
    cluster at once. Whatever isn't claimed this pass -- because of either limit -- is
    picked up over subsequent gap-fill passes instead. Returns how many points were
    actually claimed, so the caller can track the per-call budget across multiple gaps the
    same plate dominates."""
    old_points, old_elevation = plate.all_points_and_elevation()
    if len(old_points) == 0:
        claimed_points = gap_points[:max_new_points]
    else:
        dist_to_old, _ = cKDTree(old_points).query(gap_points)
        within_ring = dist_to_old < GROWTH_RING_RAD
        candidates = gap_points[within_ring]
        order = np.argsort(dist_to_old[within_ring])
        claimed_points = candidates[order][:max_new_points]
    if len(claimed_points) == 0:
        return 0

    base = base_elevation(plate.crust_type)
    amp = noise_amplitude(plate.crust_type)
    noise = SphereNoise(rng, octaves=3, base_freq=2.5)
    gap_elevation = base + amp * noise.sample(claimed_points)

    combined_points = np.concatenate([old_points, claimed_points], axis=0)
    combined_elevation = np.concatenate([old_elevation, gap_elevation], axis=0)
    tree = cKDTree(combined_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < COVERAGE_RADIUS_RAD

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        _, idx = tree.query(world_pts)
        return combined_elevation[idx]

    plate.lines = build_lines_from_lattice(plate.frame, is_owned, elevation_at)
    return len(claimed_points)


def _spawn_plate_from_gap(world: "World", gap_points: np.ndarray, rng: np.random.Generator) -> Plate:
    centroid = geometry.normalize(gap_points.mean(axis=0))
    frame = geometry.plate_frame_from_seed(centroid)
    crust_type = "oceanic"  # new crust between separating plates

    base = base_elevation(crust_type)
    amp = noise_amplitude(crust_type)
    noise = SphereNoise(rng, octaves=3, base_freq=2.5)
    tree = cKDTree(gap_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < COVERAGE_RADIUS_RAD

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return base + amp * noise.sample(world_pts)

    lines = build_lines_from_lattice(frame, is_owned, elevation_at)
    plate = Plate(plate_id=world.next_plate_id, frame=frame, crust_type=crust_type, lines=lines)
    world.next_plate_id += 1

    points, _ = plate.all_points_and_elevation()
    if len(points) > 0:
        velocities = mantle.flow_at(points, world.mantle_centers)
        plate.omega = mantle.clamp_rate(mantle.fit_euler_pole(points, velocities))
    return plate


def fill_gaps(world: "World") -> list[str]:
    """Find every sphere region no plate currently covers and resolve it (absorb into a
    single bordering plate, or spawn a new one). Mutates world.plates/world.next_plate_id
    in place. Returns one human-readable event message per newly spawned plate, for the UI's
    event console -- absorption isn't logged, since it grows a plate that already exists
    rather than adding or removing one (see world.step_world for where these get logged)."""
    existing_points, existing_owner = _all_existing_points(world)
    if len(existing_points) == 0:
        return []
    existing_tree = cKDTree(existing_points)

    gap_points = _find_gap_points(existing_tree)
    if len(gap_points) == 0:
        return []

    labels = _cluster(gap_points, CLUSTER_RADIUS_RAD)
    rng = np.random.default_rng((world.seed, world.gap_fill_calls))
    world.gap_fill_calls += 1

    plate_by_id = {p.plate_id: p for p in world.plates}
    absorb_budget: dict[int, int] = {}
    new_plates: list[Plate] = []
    for label in np.unique(labels):
        cluster_points = gap_points[labels == label]
        if len(cluster_points) < MIN_GAP_POINTS:
            continue

        dominant = _preferred_border_plate(cluster_points, existing_tree, existing_owner, plate_by_id)
        if dominant is not None:
            remaining = absorb_budget.setdefault(dominant, MAX_ABSORB_NODES_PER_PLATE_PER_CALL)
            if remaining <= 0:
                continue
            claimed = _absorb_gap_into_plate(plate_by_id[dominant], cluster_points, rng, remaining)
            absorb_budget[dominant] -= claimed
        else:
            new_plates.append(_spawn_plate_from_gap(world, cluster_points, rng))

    world.plates.extend(new_plates)
    return [
        f"Plate {plate.plate_id} ({plate.crust_type}) formed from an uncovered gap."
        for plate in new_plates
    ]
