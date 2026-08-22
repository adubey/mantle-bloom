"""Plate topology changes: consumption, continental-collision merging, and mantle-flow-
driven splitting.

Unlike rotation and boundary evolution, these are rare, discrete, topology-changing events,
so a one-time resample onto a fresh local lattice (merge) or a plane cut through existing
nodes (split) is an acceptable cost here -- the exact, no-resampling guarantee only matters
for routine per-step motion.

Merging in particular is deliberately slow: a pair of continental plates has to stay
continuously close and converging (see find_continental_collision_pairs) for a sustained,
randomized 50-100 Myr (COLLISION_MERGE_MIN/MAX_YEARS, tracked in
World.collision_progress) before they actually fuse, and at most one merge happens per
step -- see update_collision_progress and apply_topology_changes. apply_topology_changes
returns a list of human-readable event strings for whatever happened, which world.step_world
threads through to the API for the frontend's event console.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial import cKDTree

from . import geometry, mantle
from .boundary import MERGE_THRESHOLD_RAD, TRANSFORM_RATE_THRESHOLD, closing_rate
from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate, build_lines_from_lattice, line_spacing_rad

if TYPE_CHECKING:
    from .world import World

# Reference (World.node_density == 1.0) values -- kept as bare module constants (still the
# default any direct caller/test without a World gets), scaled at the point of use inside
# find_continental_collision_pairs/merge_plates by this world's own node_density instead of
# read bare, same reasoning as boundary.py's own _far_threshold_rad and friends.
MERGE_CONTACT_DISTANCE_RAD = MERGE_THRESHOLD_RAD
MERGE_MIN_CONTACT_NODES = 4
MERGE_COVERAGE_RADIUS_RAD = 1.2 * TARGET_LINE_SPACING_RAD
# A single step can already show a real (not just proximity-driven) closing rate over part
# of a shared boundary without the two plates being in anything like a genuine, sustained
# collision -- a curving boundary is often locally convergent in one stretch even while the
# plates as a whole are just sliding past each other. Checking the instantaneous rate alone
# still merged plates within a step or two of first contact. Real continental collisions
# play out over tens of millions of years, so a pair only actually merges once they've been
# continuously close-and-converging (see find_continental_collision_pairs) for a sustained
# duration, randomized per pair within this range so unrelated collisions don't all resolve
# in lockstep.
COLLISION_MERGE_MIN_YEARS = 50_000_000
COLLISION_MERGE_MAX_YEARS = 100_000_000

# A node count, not a distance -- doesn't scale automatically with TARGET_LINE_SPACING_RAD
# the way the distance-based thresholds elsewhere do. Node count for a given physical plate
# area scales with the *square* of resolution (more rows and more samples per row), so this
# is scaled by ~4x alongside plates.py's most recent resolution doubling, to keep meaning
# "the same physical plate size" rather than silently becoming 4x more split-happy.
SPLIT_MIN_NODES = 1200
# A single rigid rotation essentially never fits a wide footprint's flow samples exactly
# -- any spatially-varying field sampled over a large angular extent has *some* residual,
# even for a plate that has no business splitting. These thresholds need to sit well above
# that everyday baseline, or the split check fires on ordinary plates every time it runs.
SPLIT_RMS_RESIDUAL_THRESHOLD = mantle.cm_per_yr_to_rad_per_yr(9.0)
SPLIT_MIN_POLE_SEPARATION = mantle.cm_per_yr_to_rad_per_yr(6.0)
# Cooldown after a plate is created (by generation, split, or merge) before it's eligible
# to split again. Without this, a freshly-split daughter plate -- still not a perfect rigid-
# rotation fit, since the field it was cut from was continuous, not truly bimodal -- would
# often clear the thresholds above again on the very next step, and the next, recursively
# slicing the original plate into many thin near-parallel slivers within a handful of
# steps (this happened during development: what looked like elevation "banding" turned out
# to be dozens of sliver plates, each rotating almost identically to its neighbors).
SPLIT_MIN_AGE_STEPS = 15


def remove_defunct_plates(world: "World") -> None:
    """A plate whose every elevation node was deleted (fully subducted, see boundary.py), or
    that's been eroded down to a single line (or none) -- no real remaining territory, just
    a sliver along one latitude -- simply vanishes. No special-cased merge algorithm needed
    either way; see apply_topology_changes for the distinct log messages for each case."""
    world.plates = [p for p in world.plates if len(p.lines) > 1]


def find_continental_collision_pairs(world: "World") -> list[tuple[int, int]]:
    """Continental plate pairs that are both close *and* actively converging.

    Distance alone isn't enough: plates.py's tiling has no gaps, so every pair of
    neighboring plates is already touching along their shared boundary the moment they're
    generated, whether that boundary is convergent, divergent, or transform. Checking
    distance only would flag ordinary neighbors as "colliding" before they've moved at all
    -- confirmed directly: for a meaningful fraction of random seeds, this fired on the very
    first step regardless of how small `years` was, because the pre-existing generation-time
    proximity was already enough on its own. Requiring a real closing rate (the same check
    boundary.py uses to classify convergent boundaries) is what actually distinguishes a
    genuine collision from any other pair of neighbors.

    Most continental pairs are nowhere near each other (separated by oceanic plates, or just
    on opposite sides of the sphere), and building/querying a full point-cloud k-d tree for
    every one of the O(n^2) pairs dominated step time once plates carry thousands of nodes
    each (see docs/simulation-model.md's resolution note). Each plate's cheap bounding
    sphere (see geometry.bounding_sphere) rules most pairs out with one arccos instead of a
    full tree query, and each plate's tree is built once and reused across every pair it
    appears in, rather than rebuilt per pair."""
    # Scaled to this world's own node_density (see the module constants' own comment) --
    # ratio against TARGET_LINE_SPACING_RAD rather than a fresh 0.4 multiplier, so this
    # doesn't duplicate boundary.py's own MERGE_THRESHOLD_RAD derivation.
    merge_contact_distance_rad = MERGE_THRESHOLD_RAD * (line_spacing_rad(world.node_density) / TARGET_LINE_SPACING_RAD)

    continental = [p for p in world.plates if p.crust_type == "continental" and p.node_count() > 0]
    points = {p.plate_id: p.all_points_and_elevation()[0] for p in continental}
    spheres = {pid: geometry.bounding_sphere(pts) for pid, pts in points.items()}
    trees: dict[int, cKDTree] = {}

    pairs = []
    for i, a in enumerate(continental):
        for b in continental[i + 1 :]:
            pa, pb = points[a.plate_id], points[b.plate_id]
            if len(pa) == 0 or len(pb) == 0:
                continue

            ca, ra = spheres[a.plate_id]
            cb, rb = spheres[b.plate_id]
            centroid_dist = float(geometry.angular_distance(ca, cb))
            if centroid_dist - ra - rb > merge_contact_distance_rad:
                continue  # no point in either cloud can possibly be within contact distance

            if b.plate_id not in trees:
                trees[b.plate_id] = cKDTree(pb)

            # Same triangle-inequality reasoning as the pair-level bounding-sphere prune
            # above, just applied per point: a point of `pa` farther than rb +
            # merge_contact_distance_rad from cb cannot possibly be within contact distance
            # of any point in pb, so it can't affect the result below. Without this, the
            # k-d tree query -- the actual cost once a pair survives the prune above -- ran
            # over a continental plate's *entire* node cloud (thousands of interior nodes
            # nowhere near the other plate) instead of just the handful of nodes near the
            # shared edge that could ever come back "close".
            candidate_mask = geometry.angular_distance(pa, cb) <= rb + merge_contact_distance_rad
            if not np.any(candidate_mask):
                continue
            candidate_idx = np.nonzero(candidate_mask)[0]
            dist, idx = trees[b.plate_id].query(pa[candidate_idx])
            close = dist < merge_contact_distance_rad
            if np.sum(close) < MERGE_MIN_CONTACT_NODES:
                continue

            close_points = pa[candidate_idx[close]]
            neighbor_points = pb[idx[close]]
            closing = closing_rate(close_points, a.omega, b.omega, neighbor_points)
            if np.sum(closing > TRANSFORM_RATE_THRESHOLD) >= MERGE_MIN_CONTACT_NODES:
                # Sorted so a pair's identity as a dict key (collision-duration tracking,
                # see update_collision_progress) doesn't depend on world.plates' current
                # iteration order, which can change step to step.
                pairs.append(tuple(sorted((a.plate_id, b.plate_id))))
    return pairs


def _collision_threshold_years(seed: int, pair: tuple[int, int]) -> float:
    """How long `pair` needs to stay continuously close-and-converging before it actually
    merges -- deterministic per (world seed, pair) so it doesn't need to be stored, just
    recomputed whenever it's needed."""
    rng = np.random.default_rng((seed, pair[0], pair[1]))
    return float(rng.uniform(COLLISION_MERGE_MIN_YEARS, COLLISION_MERGE_MAX_YEARS))


def update_collision_progress(world: "World", years: float) -> list[tuple[int, int]]:
    """Advance sustained-collision tracking by `years` for every pair currently close and
    converging (world.collision_progress: pair -> accumulated convergent years). Returns
    pairs that have now accumulated enough to merge. A pair that stops being close-and-
    converging before reaching its threshold has its progress dropped entirely -- the
    collision didn't sustain, so it doesn't get partial credit toward a future one.

    A collision *starting* isn't logged to the UI's event console (see
    apply_topology_changes) -- plates.py's tiling has every neighbor pair already touching
    at generation, so a real fraction of pairs pass this check at some point without ever
    actually merging; only the outcomes that actually change the world (a completed merge,
    a plate disappearing) are worth surfacing there."""
    current_pairs = find_continental_collision_pairs(world)
    current_set = set(current_pairs)

    ready = []
    for pair in current_pairs:
        accumulated = world.collision_progress.get(pair, 0.0) + years
        world.collision_progress[pair] = accumulated
        if accumulated >= _collision_threshold_years(world.seed, pair):
            ready.append(pair)

    for pair in list(world.collision_progress):
        if pair not in current_set:
            del world.collision_progress[pair]

    return ready


def merge_plates(world: "World", id_keep: int, id_absorb: int) -> None:
    """Fuse `id_absorb` into `id_keep`: keep `id_keep`'s frame, resample the union
    footprint from scratch (a one-time nearest-neighbor lookup into the pre-merge combined
    point cloud), and drop `id_absorb`."""
    keep = next(p for p in world.plates if p.plate_id == id_keep)
    absorb = next(p for p in world.plates if p.plate_id == id_absorb)

    spacing_rad = line_spacing_rad(world.node_density)
    coverage_radius_rad = MERGE_COVERAGE_RADIUS_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)

    keep_pts, keep_elev = keep.all_points_and_elevation()
    absorb_pts, absorb_elev = absorb.all_points_and_elevation()
    old_points = np.concatenate([keep_pts, absorb_pts], axis=0)
    old_elevation = np.concatenate([keep_elev, absorb_elev], axis=0)
    tree = cKDTree(old_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < coverage_radius_rad

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        _, idx = tree.query(world_pts)
        return old_elevation[idx]

    keep.lines = build_lines_from_lattice(keep.frame, is_owned, elevation_at, spacing_rad=spacing_rad)
    keep.omega = mantle.clamp_rate((keep.omega + absorb.omega) / 2.0)
    keep.age_steps = 0
    world.plates = [p for p in world.plates if p.plate_id != id_absorb]


def _fit_residual_rms(points: np.ndarray, velocities: np.ndarray, omega: np.ndarray) -> float:
    predicted = np.cross(omega, points)
    return float(np.sqrt(np.mean(np.sum((predicted - velocities) ** 2, axis=-1))))


def maybe_split_plate(world: "World", plate: Plate) -> tuple[Plate, Plate] | None:
    """If a single rigid rotation poorly explains the mantle flow across this plate's
    footprint, cluster the flow into two regimes and, if they're genuinely different, cut
    the plate along the great circle separating them."""
    # An area-based count -- scales with node_density directly (see SPLIT_MIN_NODES' own
    # comment), not with spacing_rad the way the distance thresholds elsewhere do.
    min_nodes = max(1, round(SPLIT_MIN_NODES * world.node_density))
    if plate.node_count() < min_nodes:
        return None
    if plate.age_steps < SPLIT_MIN_AGE_STEPS:
        return None

    points, _ = plate.all_points_and_elevation()
    velocities = mantle.flow_at(points, world.mantle_centers)
    if _fit_residual_rms(points, velocities, plate.omega) < SPLIT_RMS_RESIDUAL_THRESHOLD:
        return None

    _, labels = kmeans2(velocities, k=2, minit="++", seed=plate.plate_id)
    if len(np.unique(labels)) < 2:
        return None

    mask_a = labels == 0
    mask_b = labels == 1
    pole_a = mantle.fit_euler_pole(points[mask_a], velocities[mask_a])
    pole_b = mantle.fit_euler_pole(points[mask_b], velocities[mask_b])
    if np.linalg.norm(pole_a - pole_b) < SPLIT_MIN_POLE_SEPARATION:
        return None

    centroid_a = geometry.normalize(points[mask_a].mean(axis=0))
    centroid_b = geometry.normalize(points[mask_b].mean(axis=0))
    # The great circle equidistant from two points has normal (a - b): P.a == P.b iff
    # P.(a-b) == 0.
    cut_normal = geometry.normalize(centroid_a - centroid_b)

    lines_a: list[ElevationLine] = []
    lines_b: list[ElevationLine] = []
    for line in plate.lines:
        world_pts = line.world_xyz(plate.frame)
        side = np.sum(world_pts * cut_normal, axis=-1) > 0
        if np.any(side):
            lines_a.append(
                ElevationLine(
                    phi=line.phi,
                    theta=line.theta[side],
                    elevation=line.elevation[side],
                    channel_depth=line.channel_depth[side],
                    channel_width=line.channel_width[side],
                    lake_depth=line.lake_depth[side],
                    glacier_depth=line.glacier_depth[side],
                    silt_depth=line.silt_depth[side],
                    is_volcano=line.is_volcano[side],
                    volcano_active_years_remaining=line.volcano_active_years_remaining[side],
                    soil_depth=line.soil_depth[side],
                    soil_mineral_content=line.soil_mineral_content[side],
                    soil_organic_content=line.soil_organic_content[side],
                    coal_deposit_m=line.coal_deposit_m[side],
                    oil_gas_deposit_m=line.oil_gas_deposit_m[side],
                    mineral_deposit_m=line.mineral_deposit_m[side],
                )
            )
        if np.any(~side):
            lines_b.append(
                ElevationLine(
                    phi=line.phi,
                    theta=line.theta[~side],
                    elevation=line.elevation[~side],
                    channel_depth=line.channel_depth[~side],
                    channel_width=line.channel_width[~side],
                    lake_depth=line.lake_depth[~side],
                    glacier_depth=line.glacier_depth[~side],
                    silt_depth=line.silt_depth[~side],
                    is_volcano=line.is_volcano[~side],
                    volcano_active_years_remaining=line.volcano_active_years_remaining[~side],
                    soil_depth=line.soil_depth[~side],
                    soil_mineral_content=line.soil_mineral_content[~side],
                    soil_organic_content=line.soil_organic_content[~side],
                    coal_deposit_m=line.coal_deposit_m[~side],
                    oil_gas_deposit_m=line.oil_gas_deposit_m[~side],
                    mineral_deposit_m=line.mineral_deposit_m[~side],
                )
            )

    if sum(len(l.theta) for l in lines_a) < min_nodes or sum(len(l.theta) for l in lines_b) < min_nodes:
        return None

    new_id = world.next_plate_id
    world.next_plate_id += 1
    plate_a = Plate(
        plate_id=plate.plate_id,
        frame=plate.frame.copy(),
        crust_type=plate.crust_type,
        omega=mantle.clamp_rate(pole_a),
        lines=lines_a,
    )
    plate_b = Plate(
        plate_id=new_id,
        frame=plate.frame.copy(),
        crust_type=plate.crust_type,
        omega=mantle.clamp_rate(pole_b),
        lines=lines_b,
    )
    return plate_a, plate_b


def apply_topology_changes(world: "World", years: float) -> list[str]:
    """Consumption, then at most one collision merge, then splits. Returns human-readable
    event messages for anything that happened, for the UI's event console -- a plate
    disappearing (fully subducted, or eroded down to no real land), two plates merging, or
    a new plate being created by a split. A collision merely *starting* is deliberately not
    logged here -- see update_collision_progress."""
    events: list[str] = []
    for plate in world.plates:
        plate.age_steps += 1

    consumed = [p for p in world.plates if p.node_count() == 0]
    for p in consumed:
        events.append(f"Plate {p.plate_id} ({p.crust_type}) was fully subducted and disappeared.")

    # "No land" here means no real remaining territory -- reduced to a single line (or, if
    # it also has zero nodes, already counted as consumed above) -- not the crust_type sense
    # of dry-land-above-sea-level; see remove_defunct_plates.
    no_land = [p for p in world.plates if p.node_count() > 0 and len(p.lines) <= 1]
    for p in no_land:
        events.append(f"Plate {p.plate_id} ({p.crust_type}) had no land left and disappeared.")

    remove_defunct_plates(world)

    ready_pairs = update_collision_progress(world, years)
    if ready_pairs:
        # Real continental collisions don't resolve all at once, and merging every ready
        # pair in the same call could still cascade through a whole chain of plates in one
        # step -- fuse at most one collision per step, same as every other change here.
        id_keep, id_absorb = ready_pairs[0]
        elapsed_years = world.collision_progress.pop((id_keep, id_absorb), 0.0)
        merge_plates(world, id_keep, id_absorb)
        events.append(
            f"Plates {id_keep} and {id_absorb} collided and merged into plate {id_keep} "
            f"after {elapsed_years / 1e6:.0f} million years."
        )

    new_plates: list[Plate] = []
    for plate in world.plates:
        split_result = maybe_split_plate(world, plate)
        if split_result is None:
            new_plates.append(plate)
        else:
            new_plates.extend(split_result)
            events.append(
                f"Plate {plate.plate_id} split into plates {split_result[0].plate_id} "
                f"and {split_result[1].plate_id}."
            )
    world.plates = new_plates
    return events
