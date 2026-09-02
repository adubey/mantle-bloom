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

from . import geometry, mantle, plates as plates_mod
from .boundary import MERGE_THRESHOLD_RAD, TRANSFORM_RATE_THRESHOLD, closing_rate
from .elevation_lines import TARGET_LINE_SPACING_RAD, line_spacing_rad
from .plates import Plate, query_workers

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

# Even a pair that's collided long enough shouldn't merge for certain if the result would be
# a large chunk of the whole world -- real supercontinents form slowly and are the exception,
# not the norm every sustained collision should produce. Combined *node-count* share of the
# world is resolution-independent (a fraction, not a raw count) and directly targets what a
# merge actually produces: how much of the sphere the surviving plate would cover.
MERGE_SIZE_UNLIKELY_FRACTION = 0.25
MERGE_PROBABILITY_FLOOR = 0.02  # never impossible -- real supercontinents do form eventually

# A node count, not a distance -- doesn't scale automatically with TARGET_LINE_SPACING_RAD
# the way the distance-based thresholds elsewhere do. Node count for a given physical plate
# area scales with the *square* of resolution (more rows and more samples per row), so this
# is scaled by node_density directly at the point of use.
#
# Halved 1200 -> 600 (2026): the great-circle cut between the two k-means flow centroids is
# frequently lopsided, so at 1200 a plate needed well over 2*1200 well-distributed nodes for
# *both* halves to clear the floor -- instrumenting a real run showed plates repeatedly
# clearing both physics gates only to have plate.split() reject the cut on size, and true
# rift-splits essentially never happening while merges/consumption ran unchecked (reported:
# "the number of plates is decreasing over time"). 600 still means each half is a genuinely
# plate-sized fragment, not a near-continent; SPLIT_MIN_AGE_STEPS is what actually blocks the
# freshly-split-daughter sliver cascade.
SPLIT_MIN_NODES = 700
# A single rigid rotation essentially never fits a wide footprint's flow samples exactly
# -- any spatially-varying field sampled over a large angular extent has *some* residual,
# even for a plate that has no business splitting. These thresholds need to sit well above
# that everyday baseline, or the split check fires on ordinary plates every time it runs.
# Lowered from 9.0/6.0 cm/yr (2026) -- with the torque engine a large continental plate gets
# a *good* rigid-rotation fit (measured residual ~0.3-0.7x the old threshold), so ordinary
# supercontinent-scale plates never tripped the residual gate and never rifted; 6.0/4.0 is
# still ~3-4x the everyday oceanic-plate baseline the instrumented run showed.
SPLIT_RMS_RESIDUAL_THRESHOLD = mantle.cm_per_yr_to_rad_per_yr(6.0)
SPLIT_MIN_POLE_SEPARATION = mantle.cm_per_yr_to_rad_per_yr(4.0)
# Cooldown after a plate is created (by generation, split, or merge) before it's eligible
# to split again. Without this, a freshly-split daughter plate -- still not a perfect rigid-
# rotation fit, since the field it was cut from was continuous, not truly bimodal -- would
# often clear the thresholds above again on the very next step, and the next, recursively
# slicing the original plate into many thin near-parallel slivers within a handful of
# steps (this happened during development: what looked like elevation "banding" turned out
# to be dozens of sliver plates, each rotating almost identically to its neighbors).
# Raised 15 -> 20 (2026) alongside the more permissive split gates below -- it also spaces
# out the otherwise-synchronized first rift of every freshly-generated plate (they all clear
# their cooldown on the same step).
SPLIT_MIN_AGE_STEPS = 20

# A plate's own sheer size independently raises its odds of rifting, on top of (not instead
# of) the mantle-flow-fit criteria above -- a bigger footprint is both more likely to
# straddle genuinely different mantle flow regimes and, mechanically, more likely to run its
# local (phi, theta) parametrization into trouble the longer it's left uncut (a boundary
# line that ends up spiraling many times around a plate's own local-frame pole as the plate
# grows past it is exactly what overwhelmed elevation_lines.py's periodic regularization pass
# once). Modeled as linearly relaxing the two split gates (residual-fit and pole-separation)
# toward zero as the plate's own angular radius (`geometry.bounding_sphere`, centroid to
# farthest node) approaches SPLIT_SIZE_CERTAIN_RIFT_RAD -- at that size essentially any
# mantle-flow variation at all, however small, is enough to both trigger and pass a split,
# i.e. rifting probability is effectively 100%. At radius 0 this changes nothing (both gates
# at their normal full strength); the plate still needs kmeans to find two genuinely
# separate clusters either way, since there's no sensible cut without one.
#
# Lowered pi -> 2.2 rad (2026): at pi the relaxation only bit for near-hemisphere plates, so
# real supercontinent-scale plates (radius ~1.3-2.0 rad) barely felt it and effectively
# never rifted. 2.2 rad still forces a genuinely huge plate to rift eventually while leaving
# ordinary large plates to the flow-fit gates -- 1.8 was tried first and shattered every
# freshly-generated plate within the first ~30 Myr (they all start near radius 1.3-1.5).
SPLIT_SIZE_CERTAIN_RIFT_RAD = 2.2

# Defragmentation (see defragment_plates / Plate.defragment). Ordinary deform() never
# deletes a line's last node and only ever shrinks a line's ends, so subduction/transform
# can sever a plate's node cloud into two disconnected landmasses (still carried as one
# Plate) or leave a comb of stranded one-node rows behind; maybe_split_plate only cuts on
# mantle-flow disagreement, not geometry, so neither case is caught. This pass finds them.
#
# Two nodes count as connected if they're within this multiple of the world's own line
# spacing -- one row-step is ~1x spacing in phi and >= 1x in theta, a diagonal neighbour
# ~1.4x, so 2.5x comfortably links a genuinely contiguous patch while still separating two
# lobes across a real (>~300km) subduction gap. Validated against real saved worlds: every
# healthy plate comes back as a single component at this radius.
DEFRAG_CONNECT_RADIUS_MULT = 2.5
# A component smaller than this becomes stranded crust and is dropped rather than promoted
# to its own plate. A node count, not a distance -- scales with node_density directly (an
# area), same reasoning as SPLIT_MIN_NODES.
DEFRAG_FRAGMENT_MIN_NODES = 50
# Cadence: this is a whole-world O(nodes) k-d-tree pass, cheap but not free, and plate
# topology doesn't fragment fast. cf. the removed reassign.py's REASSIGN_INTERVAL_STEPS = 5.
DEFRAG_INTERVAL_STEPS = 4


def remove_defunct_plates(world: "World") -> None:
    """A plate whose every elevation node was deleted (fully subducted, see boundary.py), or
    that's been eroded down to a single line (or none) -- no real remaining territory, just
    a sliver along one latitude -- simply vanishes. No special-cased merge algorithm needed
    either way; see apply_topology_changes for the distinct log messages for each case.

    `node_count() > 0` is its own check, not implied by `not p.has_negligible_territory()`: a
    plate can have two or more lines that have each individually shrunk to zero nodes (see
    `_grow_or_shrink_line`) without ever dropping below the line-count threshold, which
    would otherwise leave an empty-but-not-removed plate sitting in world.plates -- e.g.
    still counted by /world/summary's num_plates, still iterated by every other per-step
    pass -- indefinitely. Called every step via apply_topology_changes, so this never lingers
    more than one step. `has_negligible_territory` is representation-generic (`Plate`'s own
    method, see plates.py) -- `PlateWithLines` still means "at most one line left" by it,
    just expressed through the abstract interface now instead of reaching into `.lines`
    directly, so this works for any `Plate` subclass, not just that one."""
    world.plates = [p for p in world.plates if p.node_count() > 0 and not p.has_negligible_territory()]


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
            dist, idx = trees[b.plate_id].query(pa[candidate_idx], workers=query_workers(len(candidate_idx)))
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


def _merge_probability(world: "World", pair: tuple[int, int]) -> float:
    """How likely `pair` is to actually merge this step, given it's already met the sustained-
    duration requirement -- 1.0 for two small plates, relaxing linearly toward
    MERGE_PROBABILITY_FLOOR as their combined share of the world's total nodes approaches
    MERGE_SIZE_UNLIKELY_FRACTION. Duration alone isn't enough for a pair that would produce a
    large chunk of the whole world; this is an independent, size-only gate on top of it."""
    total_nodes = sum(p.node_count() for p in world.plates)
    if total_nodes == 0:
        return 1.0
    a = next(p for p in world.plates if p.plate_id == pair[0])
    b = next(p for p in world.plates if p.plate_id == pair[1])
    combined_frac = (a.node_count() + b.node_count()) / total_nodes
    size_frac = min(1.0, combined_frac / MERGE_SIZE_UNLIKELY_FRACTION)
    return 1.0 - size_frac * (1.0 - MERGE_PROBABILITY_FLOOR)


def update_collision_progress(world: "World", years: float) -> list[tuple[int, int]]:
    """Advance sustained-collision tracking by `years` for every pair currently close and
    converging (world.collision_progress: pair -> accumulated convergent years). Returns
    pairs that have now accumulated enough to merge *and* cleared this step's size-based
    probability roll (see _merge_probability) -- a pair that clears the duration threshold but
    fails the roll keeps its accumulated progress untouched (it's still >= threshold, so it's
    re-rolled fresh, with a freshly recomputed size, on every subsequent step it stays close-
    and-converging) rather than losing credit toward a future attempt. A pair that stops being
    close-and-converging before reaching its threshold has its progress dropped entirely --
    the collision didn't sustain, so it doesn't get partial credit toward a future one.

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
            # Keyed by elapsed_years (rounded, since np.random.default_rng's seed tuple only
            # accepts ints) alongside the pair, so a pair that fails this roll gets an
            # independent fresh one next step rather than being stuck forever.
            rng = np.random.default_rng((world.seed, pair[0], pair[1], round(world.elapsed_years)))
            if rng.random() < _merge_probability(world, pair):
                ready.append(pair)

    for pair in list(world.collision_progress):
        if pair not in current_set:
            del world.collision_progress[pair]

    return ready


def merge_plates(world: "World", id_keep: int, id_absorb: int) -> None:
    """Fuse `id_absorb` into `id_keep`: keep `id_keep`'s frame, absorb `id_absorb`'s territory
    (see Plate.merge_with for how each representation actually folds the two node sets
    together), and drop `id_absorb`.

    Exclusivity against every *other* live plate is enforced the same way initial generation
    guarantees it (see plates.generate_plates' nearest-seed Voronoi query) -- see
    Plate.merge_with's own docstring for the cross-plate overlap bug this guards against."""
    keep = next(p for p in world.plates if p.plate_id == id_keep)
    absorb = next(p for p in world.plates if p.plate_id == id_absorb)

    spacing_rad = line_spacing_rad(world.node_density)
    coverage_radius_rad = MERGE_COVERAGE_RADIUS_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)

    other_points_list = [p.all_points_and_elevation()[0] for p in world.plates if p.plate_id not in (id_keep, id_absorb)]
    other_points = np.concatenate(other_points_list, axis=0) if other_points_list else np.zeros((0, 3))

    keep.merge_with(absorb, spacing_rad, coverage_radius_rad, other_points)
    world.plates = [p for p in world.plates if p.plate_id != id_absorb]


def _fit_residual_rms(points: np.ndarray, velocities: np.ndarray, omega: np.ndarray) -> float:
    predicted = np.cross(omega, points)
    return float(np.sqrt(np.mean(np.sum((predicted - velocities) ** 2, axis=-1))))


def maybe_split_plate(world: "World", plate: Plate) -> tuple[Plate, Plate] | None:
    """If a single rigid rotation poorly explains the mantle flow across this plate's
    footprint, cluster the flow into two regimes and, if they're genuinely different, cut
    the plate along the great circle separating them (see Plate.split for the actual
    representation-specific partition)."""
    # An area-based count -- scales with node_density directly (see SPLIT_MIN_NODES' own
    # comment), not with spacing_rad the way the distance thresholds elsewhere do.
    min_nodes = max(1, round(SPLIT_MIN_NODES * world.node_density))
    if plate.node_count() < min_nodes:
        return None
    if plate.age_steps < SPLIT_MIN_AGE_STEPS:
        return None

    points, _ = plate.all_points_and_elevation()
    velocities = mantle.flow_at(points, world.mantle_centers)

    # See SPLIT_SIZE_CERTAIN_RIFT_RAD's own comment: the bigger this plate already is, the
    # less additional mantle-flow evidence it should take to justify cutting it.
    _, radius_rad = geometry.bounding_sphere(points)
    size_frac = min(1.0, radius_rad / SPLIT_SIZE_CERTAIN_RIFT_RAD)
    residual_threshold = SPLIT_RMS_RESIDUAL_THRESHOLD * (1.0 - size_frac)
    pole_separation_threshold = SPLIT_MIN_POLE_SEPARATION * (1.0 - size_frac)

    if _fit_residual_rms(points, velocities, plate.omega) < residual_threshold:
        return None

    _, labels = kmeans2(velocities, k=2, minit="++", seed=plate.plate_id)
    if len(np.unique(labels)) < 2:
        return None

    mask_a = labels == 0
    mask_b = labels == 1
    pole_a = mantle.fit_euler_pole(points[mask_a], velocities[mask_a])
    pole_b = mantle.fit_euler_pole(points[mask_b], velocities[mask_b])
    if np.linalg.norm(pole_a - pole_b) < pole_separation_threshold:
        return None

    centroid_a = geometry.normalize(points[mask_a].mean(axis=0))
    centroid_b = geometry.normalize(points[mask_b].mean(axis=0))
    # The great circle equidistant from two points has normal (a - b): P.a == P.b iff
    # P.(a-b) == 0.
    cut_normal = geometry.normalize(centroid_a - centroid_b)

    # Peek at next_plate_id without consuming it yet -- plate.split returns None (no id
    # actually used) if either resulting half would be too small, and a rejected split
    # shouldn't burn an id.
    split_result = plate.split(world.next_plate_id, cut_normal, min_nodes)
    if split_result is None:
        return None
    world.next_plate_id += 1

    plate_a, plate_b = split_result
    plate_a.set_omega(mantle.clamp_rate(pole_a))
    plate_b.set_omega(mantle.clamp_rate(pole_b))
    return plate_a, plate_b


def defragment_plates(world: "World") -> list[str]:
    """Split any plate whose nodes form more than one disconnected landmass into that many
    plates, and drop stranded sub-fragments -- the geometric cleanup ordinary deform()/
    maybe_split_plate structurally can't do (see Plate.defragment). Mutates world.plates and
    world.next_plate_id in place; returns event strings for whatever actually changed."""
    connect_radius_rad = DEFRAG_CONNECT_RADIUS_MULT * line_spacing_rad(world.node_density)
    min_fragment_nodes = max(1, round(DEFRAG_FRAGMENT_MIN_NODES * world.node_density))

    events: list[str] = []
    new_plates: list[Plate] = []
    for plate in world.plates:
        before = plate.node_count()
        result = plate.defragment(world.next_plate_id, connect_radius_rad, min_fragment_nodes)
        if result is None:
            new_plates.append(plate)
            continue

        replacements, ids_consumed = result
        world.next_plate_id += ids_consumed
        new_plates.extend(replacements)

        shed = before - sum(p.node_count() for p in replacements)
        if len(replacements) > 1:
            spawned = ", ".join(str(p.plate_id) for p in replacements[1:])
            events.append(
                f"Plate {plate.plate_id} fragmented into disconnected landmasses; spawned plate(s) {spawned}."
            )
        if shed > 0:
            events.append(f"Plate {plate.plate_id} shed {shed} stranded nodes.")

    world.plates = new_plates

    # A fragmented/removed plate can leave stale collision-progress keys behind -- same
    # cleanup update_collision_progress does for pairs that stop being close-and-converging.
    live_ids = {p.plate_id for p in world.plates}
    for pair in list(world.collision_progress):
        if pair[0] not in live_ids or pair[1] not in live_ids:
            del world.collision_progress[pair]

    return events


def update_overlap_tracking(world: "World") -> None:
    """Diagnostic per-node bookkeeping (nothing in the physics reads it back): stamp
    `ElevationLine.overlap_onset_years` with `world.elapsed_years` on every node that has
    *just* started sitting on top of another plate's territory, and clear it (back to 0.0) on
    every node no longer overlapping anything. Run once per step from `world.step_world`
    right after `apply_topology_changes`, so it sees this step's final geometry.

    Uses the same `plates.compute_node_overlap` the Plate Inspector / `plate_diagnostics.py`
    overlap view goes through, so "15% of plate 21 is on plate 0" and "...since 178 My" are
    the same underlying node set. Mirrors `stranded_basins.reconcile_world_tracks` /
    `World.collision_progress` -- a lightweight first-seen tracker, persisted in the save."""
    tol = plates_mod.OVERLAP_TOLERANCE_MULT * line_spacing_rad(world.node_density)
    overlap = plates_mod.compute_node_overlap(world.plates, tol)
    for plate in world.plates:
        n = plate.node_count()
        if n == 0:
            continue
        info = overlap.get(plate.plate_id)
        mask = info["overlap_mask"] if info is not None else np.zeros(n, dtype=bool)
        onset = np.asarray(plate.collect("overlap_onset_years"), dtype=float).copy()
        onset[mask & (onset == 0.0)] = world.elapsed_years
        onset[~mask] = 0.0
        plate.set_fields_on_plate(overlap_onset_years=onset)


def apply_topology_changes(world: "World", years: float) -> list[str]:
    """Consumption, then at most one collision merge, then splits. Returns human-readable
    event messages for anything that happened, for the UI's event console -- a plate
    disappearing (fully subducted, or eroded down to no real land), two plates merging, or
    a new plate being created by a split. A collision merely *starting* is deliberately not
    logged here -- see update_collision_progress."""
    events: list[str] = []
    for plate in world.plates:
        plate.age_one_step()

    consumed = [p for p in world.plates if p.node_count() == 0]
    for p in consumed:
        events.append(f"Plate {p.plate_id} ({p.crust_type}) was fully subducted and disappeared.")

    # "No land" here means no real remaining territory (or, if it also has zero nodes,
    # already counted as consumed above) -- not the crust_type sense of dry-land-above-sea-
    # level; see remove_defunct_plates and Plate.has_negligible_territory.
    no_land = [p for p in world.plates if p.node_count() > 0 and p.has_negligible_territory()]
    for p in no_land:
        events.append(f"Plate {p.plate_id} ({p.crust_type}) had no land left and disappeared.")

    remove_defunct_plates(world)

    # Geometric cleanup before collision/split so a severed lobe or a ghost comb of stranded
    # nodes stops polluting neighbour polygons and collision detection. Gated to every
    # DEFRAG_INTERVAL_STEPS steps -- it's a whole-world pass and topology doesn't fragment
    # fast (see the constant's own comment).
    if world.steps_taken % DEFRAG_INTERVAL_STEPS == 0:
        events.extend(defragment_plates(world))

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

    # At most one rift per step, same incremental-change rule as the merge above -- real
    # continental breakup plays out over tens of Myr, and splitting every eligible plate in
    # one call shatters a freshly-generated world (every plate clears its identical cooldown
    # on the same step) instead of staggering the rifts over time.
    new_plates: list[Plate] = []
    split_done = False
    for plate in world.plates:
        split_result = None if split_done else maybe_split_plate(world, plate)
        if split_result is None:
            new_plates.append(plate)
        else:
            split_done = True
            new_plates.extend(split_result)
            events.append(
                f"Plate {plate.plate_id} split into plates {split_result[0].plate_id} "
                f"and {split_result[1].plate_id}."
            )
    world.plates = new_plates
    return events
