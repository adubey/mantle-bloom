"""Volcanic fields: new continental crust forming where plates are separating.

**Detection, run in the clean-up phase (world.py's `steps_since_regularize`-gated block,
alongside gaps.py/elevation_lines.py's own regularization), not every step**, in two passes:

1. Every elevation point's nearest neighbor, whole-world, unrestricted by plate -- the
   median of all these distances is "how far apart elevation points normally sit."
2. Only *boundary* points -- each plate's own line endpoints (`Plate.outline_world()`'s own
   definition of a plate's territory, reused directly rather than re-derived) -- get checked
   against each other: for each boundary point, the nearest boundary point on a *different*
   plate. If that's more than `GAP_OUTLIER_FACTOR` times the pass-1 median, it fires.

Restricting pass 2 to boundary points (rather than every point, whole-world, as an earlier
version of this did) turned out to be essential, not a stylistic choice: a boundary point
sitting right next to a genuinely wide inter-plate gap still has plenty of *same-plate*
interior neighbors much closer than that gap, so an unrestricted whole-world nearest-neighbor
query never sees the gap at all -- confirmed directly (dozens of seeds, several step sizes,
even with boundary.py's own line-growth completely disabled): real cross-plate gaps up to 5x
the typical spacing existed the entire time, invisible to that first design. Restricting the
search to boundary-vs-boundary removes the same-plate interior points that were masking the
signal, without changing what "normal spacing" means (pass 1 stays whole-world, since that's
a stable reference regardless).

Every qualifying pair contributes one new volcano point (the great-circle midpoint between
the two boundary points -- the actual empty space between the separating plates, not either
plate's own territory). i < j (candidate index order) dedupes a pair found from both
directions at once, when both sides independently qualify. All of this pass's new volcano
points are then clustered by proximity (`gaps.cluster_points`, the same connected-components
technique gaps.py's own gap-clustering already uses) and each cluster becomes one brand-new
`Plate` with `crust_type="continental"` ("volcanic fields ... result in continental plates,"
per spec -- the rock has continental physical properties regardless of whether the two
separating plates were themselves oceanic or continental). Every node of a freshly-spawned
field starts as an active volcano (see below).

A boundary point belonging to a plate that's *currently* a tracked volcanic field
(`World.volcanic_field_plate_ids`, see below) is excluded as a pass-2 *source* candidate --
without this, a field's own still-forming edge could immediately re-fire against the very
neighbor it just separated from, spawning another field on top of the last one every single
clean-up pass (the exact failure shape gaps.py's own docstring documents having fixed once
already, for its own gap-spawning). It's still a valid *target* for some other plate's own
check, so a genuinely separate rift on the field's far side isn't blocked.

**A plate stops being tracked as a volcanic field once fewer than
`VOLCANO_FRACTION_DORMANT_THRESHOLD` (5%) of its own nodes are still `is_volcano` --** not a
fixed elapsed-time countdown. `is_volcano` never reverts to False once set, so this ratio can
only ever fall, and only by dilution: as the field's own edges grow via ordinary boundary
evolution (or absorb gap territory via gaps.py), each newly-added node starts non-volcanic,
so a field that keeps growing eventually reads as "just an ordinary continental plate that
happens to have a few old volcanoes embedded in it" -- checked every step (`apply_volcanic_
activity`, alongside the eruption roll below), so the transition is caught within one step of
crossing the threshold, not lagged to the next clean-up interval.

**Eruption, run every step** (`apply_volcanic_activity`, alongside erosion.py/bathymetry.py):
each individual volcano point has its own `volcano_active_years_remaining`
(`VOLCANO_ACTIVE_MIN/MAX_YEARS`, drawn once at creation), decremented every step. While
active, it rolls a per-step eruption chance
(`1 - exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1e6)`, the same
exponential-arrival-rate shape used elsewhere in this codebase, e.g. lake evaporation's own
retention factor) and, if it erupts, adds `ERUPTION_ELEVATION_M` of new land. Deterministic
per `(seed, elapsed_years, plate_id, line_index)`, the same reproducibility precedent
merge_split.py's own per-pair collision threshold sets. `active_years_this_step` is clamped
to the volcano's own *remaining* life, not the full step size -- a large step (the UI offers
up to 10 Myr) shouldn't roll eruption chances for years past when a short-lived volcano
actually went dormant.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import gaps, geometry, mantle
from .elevation_lines import (
    MAX_ELEVATION_M,
    MIN_ELEVATION_M,
    TARGET_LINE_SPACING_RAD,
    ElevationLine,
    build_lines_from_lattice,
    iter_local_lattice,
    line_spacing_rad,
)
from .noise import SphereNoise
from .plates import Plate, PlateWithLines, base_elevation, noise_amplitude, query_workers

if TYPE_CHECKING:
    from .world import World

GAP_OUTLIER_FACTOR = 3.0
VOLCANO_ACTIVE_MIN_YEARS = 100_000
VOLCANO_ACTIVE_MAX_YEARS = 1_000_000
# Once a volcanic-field plate's own is_volcano fraction dilutes below this, it's no longer
# tracked as a field -- see module docstring.
VOLCANO_FRACTION_DORMANT_THRESHOLD = 0.05

# Deliberately much wider than gaps.py's own CLUSTER_RADIUS_RAD (~1.5x spacing, ~187km):
# gap-clustering there groups points from one dense, contiguous coverage scan, but pass 2's
# candidate points are individual boundary-point pairs spread out along a whole divergent
# boundary's length -- at gaps.py's own radius, a single long rift's many independently-
# qualifying points never merge into one cluster at all, each spawning its own tiny field
# instead (confirmed directly: 11-40+ new plates from one clean-up pass on a 10-plate world).
# 15x spacing (~1875km) merges same-rift detections into one field while still keeping
# genuinely separate rifts elsewhere on the sphere apart. Reference (World.node_density ==
# 1.0) value -- scaled at the point of use in detect_and_spawn_volcanic_fields, same
# reasoning as boundary.py's own _far_threshold_rad and friends.
VOLCANIC_FIELD_CLUSTER_RADIUS_RAD = 15.0 * TARGET_LINE_SPACING_RAD
# How far from a detected volcano point the new plate's own lattice extends -- same question
# gaps._spawn_plate_from_gap answers for its own new-plate spawning, at its own (tighter)
# scale, since this is about "the new plate's own coverage," not "how far apart do
# same-rift detections still count as one field" (VOLCANIC_FIELD_CLUSTER_RADIUS_RAD above).
# Also a reference value -- see gaps.COVERAGE_RADIUS_RAD's own comment.
VOLCANIC_FIELD_COVERAGE_RADIUS_RAD = gaps.COVERAGE_RADIUS_RAD

# Expected number of eruption events over a volcano's full active life is
# ERUPTION_RATE_PER_MYR * (active life in Myr) -- e.g. at the low end of VOLCANO_ACTIVE
# (0.1 Myr), that's 0.3 expected events (p(>=1) ~= 26%, so most short-lived volcanoes erupt
# zero or one time); at the high end (1 Myr), 3 expected events. "Occasionally," not "every
# step" or "constantly."
ERUPTION_RATE_PER_MYR = 3.0
# A single eruption's land contribution -- comparable order of magnitude to
# boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr) applied over a fraction of a Myr,
# consistent with "a discrete volcanic event" rather than a smooth continuous uplift rate.
ERUPTION_ELEVATION_M = 100.0

# Two currently-tracked volcanic-field plates whose nearest points come within this distance
# are treated as the same underlying province and get fused (see merge_close_volcanic_fields)
# -- e.g. two separate rifts that were independently detected/spawned nearby, or two fields
# that have drifted close together since. Distinct from VOLCANIC_FIELD_CLUSTER_RADIUS_RAD
# above (which clusters brand-new volcano *points* within a single detection pass, before any
# plate exists yet) -- this instead re-examines plates that already exist, on a later pass, so
# it needs its own constant even though the underlying idea ("this close together, at rift
# scale, is really one field") is the same. Reference (World.node_density == 1.0) value.
VOLCANIC_FIELD_MERGE_DISTANCE_RAD = 6.0 * TARGET_LINE_SPACING_RAD
# Which lattice nodes count as "covered" by the merged point cloud when resampling -- same
# role/value as merge_split.MERGE_COVERAGE_RADIUS_RAD/gaps.COVERAGE_RADIUS_RAD.
MERGE_FIELD_COVERAGE_RADIUS_RAD = 1.2 * TARGET_LINE_SPACING_RAD
# A handful of the closest cross-plate point pairs is enough to bridge a rift-width gap --
# every qualifying pair would produce near-duplicate, nearly-parallel bridge strips otherwise,
# without meaningfully improving coverage (the final coverage-radius union resample already
# fills in the rest around whichever bridge strips get built).
BRIDGE_MAX_PAIRS = 6

# np.random.default_rng's seed tuple only accepts ints (unlike Python's own hash-based
# tuple hashing) -- these are arbitrary distinguishing tags so merge/isolated-growth's own
# noise texture RNG streams don't collide with any other (seed, ...)-keyed stream elsewhere
# in this codebase (eruption, gap-fill, field spawning, etc. each use their own distinct
# tuple shape already).
_MERGE_RNG_TAG = 3_500_017
_ISOLATED_GROWTH_RNG_TAG = 3_500_021

# Mineral deposits: real hydrothermal circulation around an active volcanic vent precipitates
# metal-rich ore (porphyry-copper/VMS-style deposits), so mineral_deposit_m is grown right
# here, at the same eruption roll that already adds ERUPTION_ELEVATION_M -- "an eruption
# deposits mineral-rich material" is exactly what that mask already means, no separate
# detection pass needed. This only ever fires where this codebase's own volcanism model
# already places is_volcano nodes (rift-spawned volcanic fields) -- a fair collapse of real
# geology's two ore-forming settings (mid-ocean-ridge VMS deposits, arc porphyry deposits)
# into the one volcanism mechanism this simulator has, not a new detection concept.
# Monotonically non-decreasing (see plates.ElevationLine), same self-reinforcing convention
# silt_depth/coal_deposit_m/oil_gas_deposit_m already use.
MINERAL_DEPOSIT_PER_ERUPTION_M = 0.5
MAX_MINERAL_DEPOSIT_M = 20.0


def _whole_world_median_spacing(world: "World") -> float:
    """Pass 1: every elevation point's own nearest neighbor, whole-world, unrestricted by
    plate -- the median of all these distances, "how far apart elevation points normally
    sit." Returns 0.0 for a world too small to have a meaningful median."""
    points_list = [line.world_xyz(plate.frame) for plate in world.plates for line in plate.lines if len(line) > 0]
    if not points_list:
        return 0.0
    points = np.concatenate(points_list, axis=0)
    if len(points) <= 1:
        return 0.0
    # Same build-once/query-once tradeoff as boundary.step_boundaries' own per-plate tree
    # (see its comment) -- this tree is built fresh every call and queried exactly once.
    tree = cKDTree(points, balanced_tree=False, compact_nodes=False)
    dist, _ = tree.query(points, k=2, workers=query_workers(len(points)))  # column 0 is always the point itself, at distance 0
    return float(np.median(dist[:, 1]))


def _boundary_points(world: "World") -> tuple[np.ndarray, np.ndarray]:
    """Every plate's own boundary points (Plate.outline_world()'s line-endpoint definition of
    a plate's territory) and owning plate_id, concatenated -- pass 2's candidate/search
    population, deliberately much smaller than every node in the world (see module
    docstring)."""
    points_list, owner_list = [], []
    for plate in world.plates:
        outline = plate.outline_world()
        if len(outline) == 0:
            continue
        points_list.append(outline)
        owner_list.append(np.full(len(outline), plate.plate_id))
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0, dtype=int)
    return np.concatenate(points_list, axis=0), np.concatenate(owner_list, axis=0)


def _nearest_other_plate_boundary(points: np.ndarray, owner: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For each boundary point, the distance to and index of the nearest *other* boundary
    point belonging to a different plate -- one cKDTree per plate (built from every other
    plate's own boundary points), the same "query against everyone else" shape boundary.py's
    own per-step adjacency detection uses, just over the much smaller boundary-point
    population here."""
    n = len(points)
    nearest_dist = np.full(n, np.inf)
    nearest_idx = np.full(n, -1, dtype=int)
    for plate_id in np.unique(owner):
        own_mask = owner == plate_id
        other_mask = ~own_mask
        if not np.any(other_mask):
            continue
        other_indices = np.nonzero(other_mask)[0]
        # Same build-once/query-once tradeoff as boundary.step_boundaries' own per-plate
        # tree (see its comment) -- built fresh per plate_id here too, queried exactly once.
        own_points = points[own_mask]
        tree = cKDTree(points[other_mask], balanced_tree=False, compact_nodes=False)
        dist, idx = tree.query(own_points, k=1, workers=query_workers(len(own_points)))
        own_indices = np.nonzero(own_mask)[0]
        nearest_dist[own_indices] = dist
        nearest_idx[own_indices] = other_indices[idx]
    return nearest_dist, nearest_idx


def _spawn_volcanic_field_plate(world: "World", cluster_points: np.ndarray, rng: np.random.Generator) -> PlateWithLines:
    """Builds a brand-new continental Plate covering the lattice around `cluster_points`
    (this pass's own newly-detected volcano points, already clustered by proximity) --
    same `build_lines_from_lattice`-around-a-fresh-seed-frame construction
    gaps._spawn_plate_from_gap uses for an ordinary gap-spawned plate, but continental crust
    (volcanic fields "result in continental plates," regardless of what kind of plates were
    separating) and with every resulting node marked as its own active volcano."""
    spacing_rad = line_spacing_rad(world.node_density)
    coverage_radius_rad = VOLCANIC_FIELD_COVERAGE_RADIUS_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)

    centroid = geometry.normalize(cluster_points.mean(axis=0))
    frame = geometry.plate_frame_from_seed(centroid)
    crust_type = "continental"

    base = base_elevation(crust_type)
    amp = noise_amplitude(crust_type)
    noise = SphereNoise(rng, octaves=3, base_freq=2.5)
    tree = cKDTree(cluster_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < coverage_radius_rad

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return base + amp * noise.sample(world_pts)

    lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
    volcanic_lines = [
        line.replace(
            is_volcano=np.ones(len(line), dtype=bool),
            volcano_active_years_remaining=rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=len(line)),
        )
        for line in lines
    ]

    plate = PlateWithLines(plate_id=world.next_plate_id, frame=frame, crust_type=crust_type, lines=volcanic_lines)
    world.next_plate_id += 1

    points, _ = plate.all_points_and_elevation()
    if len(points) > 0:
        velocities = mantle.flow_at(points, world.mantle_centers)
        plate.set_omega(mantle.clamp_rate(mantle.fit_euler_pole(points, velocities)))
    return plate


def detect_and_spawn_volcanic_fields(world: "World") -> list[str]:
    """Finds divergent gaps (see module docstring) and spawns one new continental Plate per
    cluster of newly-detected volcano points. Mutates world.plates/world.next_plate_id/
    world.volcanic_field_plate_ids in place. Returns one human-readable event message per
    newly spawned field, for the UI's event console."""
    median_spacing = _whole_world_median_spacing(world)
    if median_spacing <= 0:
        return []

    points, owner = _boundary_points(world)
    if len(points) == 0:
        return []
    nearest_dist, nearest_idx = _nearest_other_plate_boundary(points, owner)

    is_gap_point = np.isfinite(nearest_dist) & (nearest_dist > GAP_OUTLIER_FACTOR * median_spacing)
    # A field's own still-forming edge doesn't re-fire against the neighbor it just
    # separated from -- see module docstring.
    is_gap_point &= ~np.isin(owner, list(world.volcanic_field_plate_ids))
    candidate_idx = np.nonzero(is_gap_point)[0]
    if len(candidate_idx) == 0:
        return []

    partner = nearest_idx[candidate_idx]
    keep = candidate_idx < partner  # dedupe a pair found from both directions at once
    pair_i = candidate_idx[keep]
    pair_j = partner[keep]
    if len(pair_i) == 0:
        return []

    new_volcano_points = geometry.normalize((points[pair_i] + points[pair_j]) / 2.0)

    cluster_radius_rad = VOLCANIC_FIELD_CLUSTER_RADIUS_RAD * (line_spacing_rad(world.node_density) / TARGET_LINE_SPACING_RAD)
    labels = gaps.cluster_points(new_volcano_points, cluster_radius_rad)
    rng = np.random.default_rng((world.seed, world.volcanic_field_calls))
    world.volcanic_field_calls += 1

    events = []
    for label in np.unique(labels):
        cluster = new_volcano_points[labels == label]
        plate = _spawn_volcanic_field_plate(world, cluster, rng)
        if plate.node_count() == 0:
            continue
        world.plates.append(plate)
        world.volcanic_field_plate_ids.add(plate.plate_id)
        events.append(
            f"A volcanic field emerged where plates are separating, forming new continental plate {plate.plate_id} ({plate.node_count()} nodes)."
        )
    return events


def _plate_volcano_state(plate: Plate) -> tuple[np.ndarray, np.ndarray]:
    """This plate's own is_volcano/volcano_active_years_remaining, concatenated in the exact
    same order Plate.all_points_and_elevation uses -- so the two can be indexed together
    (see merge_close_volcanic_fields)."""
    return plate.collect("is_volcano"), plate.collect("volcano_active_years_remaining")


def find_close_volcanic_field_pairs(world: "World") -> list[tuple[int, int]]:
    """Currently-tracked volcanic-field plate pairs whose outlines are within
    VOLCANIC_FIELD_MERGE_DISTANCE_RAD of each other. Unlike merge_split.
    find_continental_collision_pairs, proximity alone is enough here -- no closing-rate check
    -- because these are already newly-formed rift crust of the same kind, not two mature
    plates that happen to be neighbors purely by plates.py's generation-time tiling (see that
    function's own docstring for why *it* needs the extra check and this doesn't). Delegates
    the actual proximity test to Plate.get_neighbours -- the field count involved is
    typically small, so its per-pair bounding-sphere prescreen + outline k-d tree query is
    plenty, without needing find_continental_collision_pairs' own cross-pair tree cache."""
    merge_distance_rad = VOLCANIC_FIELD_MERGE_DISTANCE_RAD * (line_spacing_rad(world.node_density) / TARGET_LINE_SPACING_RAD)

    fields = [p for p in world.plates if p.plate_id in world.volcanic_field_plate_ids and p.node_count() > 0]
    pairs = {
        tuple(sorted((a.plate_id, b.plate_id)))
        for a in fields
        for b in a.get_neighbours(fields, threshold_rad=merge_distance_rad)
    }
    return sorted(pairs)


def _cluster_plate_pairs(pairs: list[tuple[int, int]]) -> list[set[int]]:
    """Connected components over the pair graph -- a chain of three or more mutually-close
    fields (A close to B, B close to C, even if A and C aren't directly close) merges into one
    plate in a single call, rather than needing multiple passes to fully settle one cluster.
    Only clusters of 2+ plates are returned (a plate with no close partner isn't a cluster)."""
    if not pairs:
        return []
    ids = sorted({pid for pair in pairs for pid in pair})
    index = {pid: i for i, pid in enumerate(ids)}
    n = len(ids)
    rows = [index[a] for a, b in pairs]
    cols = [index[b] for a, b in pairs]
    graph = coo_matrix((np.ones(len(pairs)), (rows, cols)), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    clusters: dict[int, set[int]] = {}
    for pid, label in zip(ids, labels):
        clusters.setdefault(int(label), set()).add(pid)
    return [cluster for cluster in clusters.values() if len(cluster) > 1]


def _slerp(a: np.ndarray, b: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Points a t-fraction of the way from unit vector `a` to unit vector `b` along their
    shared great circle (spherical, not linear, interpolation -- a bridge point needs to sit
    exactly on the sphere's surface like every other node here, not slightly inside it the way
    a naive linear blend plus renormalization would leave a wide gap's midpoint)."""
    omega = max(float(geometry.angular_distance(a, b)), 1e-9)
    sin_omega = np.sin(omega)
    t = np.asarray(t)[:, None]
    return (np.sin((1.0 - t) * omega) * a + np.sin(t * omega) * b) / sin_omega


def _bridge_points(points_a: np.ndarray, points_b: np.ndarray, spacing_rad: float, max_gap_rad: float) -> np.ndarray:
    """New intermediate points spanning the gap between `points_a` and `points_b`'s closest
    pairs -- real interpolated land connecting two merging fields, rather than leaving a hole
    once their footprints are simply unioned (build_lines_from_lattice-style resampling only
    keeps lattice nodes near an *existing* point, so without this the strip between two fields
    that were close but not yet touching would stay empty after the merge). Only pairs within
    `max_gap_rad` are bridged, and only the closest BRIDGE_MAX_PAIRS of those -- see its own
    comment."""
    if len(points_a) == 0 or len(points_b) == 0:
        return np.zeros((0, 3))
    # Same build-once/query-once tradeoff as boundary.step_boundaries' own per-plate tree
    # (see its comment) -- built fresh per call, queried exactly once.
    tree_b = cKDTree(points_b, balanced_tree=False, compact_nodes=False)
    dist, idx = tree_b.query(points_a, workers=query_workers(len(points_a)))
    close = dist < max_gap_rad
    if not np.any(close):
        return np.zeros((0, 3))

    order = np.argsort(dist[close])[:BRIDGE_MAX_PAIRS]
    a_close = points_a[close][order]
    b_close = points_b[idx[close]][order]
    d_close = dist[close][order]

    chunks = []
    for a, b, d in zip(a_close, b_close, d_close):
        n_steps = max(int(d / spacing_rad), 1)
        t = np.arange(1, n_steps) / n_steps  # exclude both endpoints -- already real nodes
        if len(t) == 0:
            continue
        chunks.append(_slerp(a, b, t))
    if not chunks:
        return np.zeros((0, 3))
    return np.concatenate(chunks, axis=0)


def _build_field_lines_with_volcano_state(
    frame: np.ndarray,
    points: np.ndarray,
    elevation: np.ndarray,
    is_volcano: np.ndarray,
    volcano_active_years_remaining: np.ndarray,
    coverage_radius_rad: float,
    spacing_rad: float,
) -> list[ElevationLine]:
    """Like plates.build_lines_from_lattice, but also carries is_volcano/volcano_active_years_
    remaining through the same nearest-neighbor lookup elevation already uses. A plain
    build_lines_from_lattice call only threads elevation, which would otherwise silently reset
    every merged node's volcanic status to False -- see plates.ElevationLine's own docstring
    for the exact bug class this avoids (the same one erosion.py/bathymetry.py's own
    reconstruction sites hit before is_volcano existed)."""
    tree = cKDTree(points)
    lines: list[ElevationLine] = []
    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        dist, idx = tree.query(world_pts)
        owned = dist < coverage_radius_rad
        if not np.any(owned):
            continue
        nearest = idx[owned]
        lines.append(
            ElevationLine(
                phi=phi,
                theta=theta_candidates[owned],
                elevation=elevation[nearest],
                is_volcano=is_volcano[nearest],
                volcano_active_years_remaining=volcano_active_years_remaining[nearest],
            )
        )
    return lines


def merge_close_volcanic_fields(world: "World") -> list[str]:
    """Finds clusters of currently-tracked volcanic-field plates that have ended up within
    VOLCANIC_FIELD_MERGE_DISTANCE_RAD of each other (independently spawned nearby, or having
    drifted close since) and fuses each cluster into a single plate -- real volcanic provinces
    this close together read as the same underlying feature, not separate fields that happen
    to be neighbors. Bridges any remaining physical gap between them with new interpolated
    land/fresh volcanoes (see _bridge_points) before the final union resample, so merging
    doesn't leave a hole where two fields were close but not yet touching. Mutates
    world.plates/world.volcanic_field_plate_ids in place. Returns one human-readable event
    message per merge, for the UI's event console."""
    clusters = _cluster_plate_pairs(find_close_volcanic_field_pairs(world))
    if not clusters:
        return []

    spacing_rad = line_spacing_rad(world.node_density)
    coverage_radius_rad = MERGE_FIELD_COVERAGE_RADIUS_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)
    bridge_gap_rad = VOLCANIC_FIELD_MERGE_DISTANCE_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)
    rng = np.random.default_rng((world.seed, _MERGE_RNG_TAG, world.volcanic_field_calls))
    world.volcanic_field_calls += 1

    plate_by_id = {p.plate_id: p for p in world.plates}
    events: list[str] = []
    removed_ids: set[int] = set()

    for cluster in clusters:
        cluster_plates = [plate_by_id[pid] for pid in cluster if pid in plate_by_id]
        if len(cluster_plates) < 2:
            continue
        # Keep whichever plate already has the most territory -- an arbitrary but stable
        # choice, the same "keep one side's frame" pattern merge_split.merge_plates uses.
        primary = max(cluster_plates, key=lambda p: p.node_count())
        others = [p for p in cluster_plates if p.plate_id != primary.plate_id]

        points_list, elevation_list, volcano_list, remaining_list = [], [], [], []
        for p in cluster_plates:
            pts, elev = p.all_points_and_elevation()
            is_volcano, remaining = _plate_volcano_state(p)
            points_list.append(pts)
            elevation_list.append(elev)
            volcano_list.append(is_volcano)
            remaining_list.append(remaining)

        bridge_chunks = [
            _bridge_points(a.all_points_and_elevation()[0], b.all_points_and_elevation()[0], spacing_rad, bridge_gap_rad)
            for i, a in enumerate(cluster_plates)
            for b in cluster_plates[i + 1 :]
        ]
        bridge_chunks = [chunk for chunk in bridge_chunks if len(chunk) > 0]
        if bridge_chunks:
            bridge_points = np.concatenate(bridge_chunks, axis=0)
            base = base_elevation("continental")
            amp = noise_amplitude("continental")
            noise = SphereNoise(rng, octaves=3, base_freq=2.5)
            points_list.append(bridge_points)
            elevation_list.append(base + amp * noise.sample(bridge_points))
            volcano_list.append(np.ones(len(bridge_points), dtype=bool))
            remaining_list.append(rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=len(bridge_points)))

        primary.set_lines(
            _build_field_lines_with_volcano_state(
                primary.frame,
                np.concatenate(points_list, axis=0),
                np.concatenate(elevation_list, axis=0),
                np.concatenate(volcano_list, axis=0),
                np.concatenate(remaining_list, axis=0),
                coverage_radius_rad,
                spacing_rad,
            )
        )
        primary.set_omega(mantle.clamp_rate(np.mean(np.stack([p.omega for p in cluster_plates], axis=0), axis=0)))
        primary.reset_age()

        absorbed_ids = [p.plate_id for p in others]
        removed_ids.update(absorbed_ids)
        world.volcanic_field_plate_ids.difference_update(absorbed_ids)
        world.volcanic_field_plate_ids.add(primary.plate_id)
        events.append(
            f"Nearby volcanic fields merged into plate {primary.plate_id} "
            f"(plates {', '.join(str(pid) for pid in absorbed_ids)} absorbed, {primary.node_count()} nodes)."
        )

    if removed_ids:
        world.plates = [p for p in world.plates if p.plate_id not in removed_ids]
    return events


# "Isolated" means no other plate's own elevation node is anywhere near this line-end -- a
# field drifting (or spawned) into open water with no neighbor at all. Deliberately much wider
# than boundary.EXTEND_THRESHOLD_RAD, which reacts to an ordinary, still-*adjacent* divergent
# neighbor that boundary.step_boundaries already grows on its own every step --
# ISOLATED_GROWTH_CLEARANCE_RAD is about there being no neighbor there at all, not about how
# wide an ordinary gap has opened. Reference (World.node_density == 1.0) value.
ISOLATED_GROWTH_CLEARANCE_RAD = 6.0 * TARGET_LINE_SPACING_RAD
# One small ring of new nodes per periodic pass -- same "grow gradually, not all at once"
# philosophy as gaps.GROWTH_RING_RAD/_absorb_gap_into_plate. There's no adjacent plate here to
# measure a gap size against (open space has no far edge to close against), so this is a fixed
# small count rather than a distance-derived one.
ISOLATED_GROWTH_NODES_PER_END = 3

# Persistent per-node fields that a brand-new node (isolated growth, or a merge bridge) starts
# at zero -- no history to carry, same reasoning boundary.py's own line-growth already uses.
_ZERO_FILL_FIELDS = (
    "channel_depth",
    "channel_width",
    "lake_depth",
    "glacier_depth",
    "silt_depth",
    "soil_depth",
    "soil_mineral_content",
    "soil_organic_content",
    "coal_deposit_m",
    "oil_gas_deposit_m",
    "mineral_deposit_m",
)


def _new_isolated_nodes(
    frame: np.ndarray, phi: float, new_theta: np.ndarray, base: float, amp: float, noise: SphereNoise, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """elevation / is_volcano (all True) / volcano_active_years_remaining for brand-new
    isolated-growth nodes -- new land forming at the field's own edge in open space, the same
    "new active volcano" convention _spawn_volcanic_field_plate uses for a field's very first
    nodes."""
    new_pts = geometry.to_world(frame, geometry.local_xyz(np.full_like(new_theta, phi), new_theta))
    elevation = base + amp * noise.sample(new_pts)
    is_volcano = np.ones(len(new_theta), dtype=bool)
    remaining = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=len(new_theta))
    return elevation, is_volcano, remaining


def grow_isolated_volcanic_fields(world: "World") -> None:
    """Extends a volcanic-field plate's own line ends outward when there's genuinely no other
    plate nearby to bound it. boundary.step_boundaries only ever grows a line relative to a
    neighbor's own motion (see its own docstring) -- a field that has drifted, or was spawned,
    into open space with no neighbor within ISOLATED_GROWTH_CLEARANCE_RAD never gets a
    boundary classification there at all, so without this it would sit frozen at that edge
    forever no matter how long the world runs. Mutates world.plates in place. Not logged to
    World.events -- ordinary incremental growth, the same reasoning gaps.py's own gap
    absorption isn't logged either."""
    if not world.volcanic_field_plate_ids:
        return

    spacing_rad = line_spacing_rad(world.node_density)
    clearance_rad = ISOLATED_GROWTH_CLEARANCE_RAD * (spacing_rad / TARGET_LINE_SPACING_RAD)

    plate_by_id = {p.plate_id: p for p in world.plates}
    for plate_id in list(world.volcanic_field_plate_ids):
        plate = plate_by_id.get(plate_id)
        if plate is None:
            continue
        lines_with_nodes = [line for line in plate.lines if len(line) > 0]
        if not lines_with_nodes:
            continue

        # Plate.get_neighbours' own bounding-sphere-plus-outline prescreen narrows the tree
        # down to plates whose territory could plausibly have a point within clearance_rad,
        # before paying for their full point clouds -- the same prefilter
        # boundary.step_boundaries uses, just via the shared helper instead of re-derived
        # here.
        nearby = plate.get_neighbours(world.plates, threshold_rad=clearance_rad)
        other_points_list = [other.all_points_and_elevation()[0] for other in nearby if other.node_count() > 0]
        # Same build-once/query-once tradeoff as boundary.step_boundaries' own per-plate
        # tree (see its comment) -- built fresh per volcanic-field plate_id here, then
        # queried only twice (high_world/low_world below, each a handful of line-ends) --
        # build cost dominates, so favoring a faster build is worth it even more here than
        # in step_boundaries.
        other_tree = (
            cKDTree(np.concatenate(other_points_list, axis=0), balanced_tree=False, compact_nodes=False)
            if other_points_list
            else None
        )

        high_phi = np.array([line.phi for line in lines_with_nodes])
        high_theta = np.array([line[-1].get_theta() for line in lines_with_nodes])
        low_theta = np.array([line[0].get_theta() for line in lines_with_nodes])
        high_world = geometry.to_world(plate.frame, geometry.local_xyz(high_phi, high_theta))
        low_world = geometry.to_world(plate.frame, geometry.local_xyz(high_phi, low_theta))

        if other_tree is None:
            high_isolated = np.ones(len(lines_with_nodes), dtype=bool)
            low_isolated = np.ones(len(lines_with_nodes), dtype=bool)
        else:
            high_dist, _ = other_tree.query(high_world, workers=query_workers(len(high_world)))
            low_dist, _ = other_tree.query(low_world, workers=query_workers(len(low_world)))
            high_isolated = high_dist > clearance_rad
            low_isolated = low_dist > clearance_rad

        if not np.any(high_isolated) and not np.any(low_isolated):
            continue

        base = base_elevation(plate.crust_type)
        amp = noise_amplitude(plate.crust_type)
        rng = np.random.default_rng((world.seed, _ISOLATED_GROWTH_RNG_TAG, plate_id, round(world.elapsed_years)))
        noise = SphereNoise(rng, octaves=3, base_freq=2.5)

        new_lines: list[ElevationLine] = []
        for i, line in enumerate(lines_with_nodes):
            if not high_isolated[i] and not low_isolated[i]:
                new_lines.append(line)
                continue

            theta = line.theta.copy()
            elevation = line.elevation.copy()
            is_volcano = line.is_volcano.copy()
            remaining = line.volcano_active_years_remaining.copy()
            fields = {name: getattr(line, name).copy() for name in _ZERO_FILL_FIELDS}
            dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)

            # High end first so the low-end index (0) is unaffected by any change made here --
            # same ordering reasoning as boundary._grow_or_shrink_line.
            if high_isolated[i]:
                new_theta = theta[-1] + dtheta * np.arange(1, ISOLATED_GROWTH_NODES_PER_END + 1)
                new_elev, new_volcano, new_remaining = _new_isolated_nodes(plate.frame, line.phi, new_theta, base, amp, noise, rng)
                theta = np.append(theta, new_theta)
                elevation = np.append(elevation, new_elev)
                is_volcano = np.append(is_volcano, new_volcano)
                remaining = np.append(remaining, new_remaining)
                for name, values in fields.items():
                    fields[name] = np.append(values, np.zeros(len(new_theta), dtype=values.dtype))

            if low_isolated[i]:
                new_theta = theta[0] - dtheta * np.arange(ISOLATED_GROWTH_NODES_PER_END, 0, -1)
                new_elev, new_volcano, new_remaining = _new_isolated_nodes(plate.frame, line.phi, new_theta, base, amp, noise, rng)
                theta = np.insert(theta, 0, new_theta)
                elevation = np.insert(elevation, 0, new_elev)
                is_volcano = np.insert(is_volcano, 0, new_volcano)
                remaining = np.insert(remaining, 0, new_remaining)
                for name, values in fields.items():
                    fields[name] = np.insert(values, 0, np.zeros(len(new_theta), dtype=values.dtype))

            new_lines.append(
                ElevationLine(
                    phi=line.phi,
                    theta=theta,
                    elevation=elevation,
                    is_volcano=is_volcano,
                    volcano_active_years_remaining=remaining,
                    **fields,
                )
            )
        plate.set_lines(new_lines)


def _volcano_fraction(plate: Plate) -> float:
    total = plate.node_count()
    if total == 0:
        return 0.0
    return int(plate.collect("is_volcano").sum()) / total


def apply_volcanic_activity(world: "World", years: float) -> list[str]:
    """Every step: relabels any volcanic-field plate whose is_volcano fraction has diluted
    below VOLCANO_FRACTION_DORMANT_THRESHOLD as an ordinary continental plate (see module
    docstring), and rolls each individual active volcano's own eruption chance, adding
    ERUPTION_ELEVATION_M wherever it erupts. Mutates world.plates/world.volcanic_field_
    plate_ids in place. Returns one human-readable event message per field that just cooled
    into ordinary continental crust this step, for the UI's event console."""
    events = []
    plate_by_id = {plate.plate_id: plate for plate in world.plates}
    for plate_id in list(world.volcanic_field_plate_ids):
        plate = plate_by_id.get(plate_id)
        if plate is None or _volcano_fraction(plate) < VOLCANO_FRACTION_DORMANT_THRESHOLD:
            world.volcanic_field_plate_ids.discard(plate_id)
            if plate is not None:
                events.append(f"The volcanic field on plate {plate_id} has cooled into ordinary continental crust.")

    for plate in world.plates:
        for line_index, line in enumerate(plate.lines):
            if len(line) == 0 or not np.any(line.is_volcano):
                continue
            active_mask = line.is_volcano & (line.volcano_active_years_remaining > 0)
            if not np.any(active_mask):
                continue

            active_years_this_step = np.minimum(years, line.volcano_active_years_remaining)
            p_erupt = 1.0 - np.exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1_000_000.0)
            rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate.plate_id, line_index))
            erupts = active_mask & (rng.random(len(line)) < p_erupt)

            new_elevation = line.elevation.copy()
            new_elevation[erupts] += ERUPTION_ELEVATION_M
            new_elevation = np.clip(new_elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            new_remaining = np.clip(line.volcano_active_years_remaining - years, 0.0, None)
            new_mineral_deposit = np.clip(
                line.mineral_deposit_m + np.where(erupts, MINERAL_DEPOSIT_PER_ERUPTION_M, 0.0), 0.0, MAX_MINERAL_DEPOSIT_M
            )

            # theta unchanged -- line.replace copies every other field (including
            # channel_width) from the existing line automatically. See plates.ElevationLine's
            # own docstring for why this pattern replaced explicit field-by-field
            # reconstruction here.
            plate.replace_line(
                line_index,
                line.replace(
                    elevation=new_elevation, volcano_active_years_remaining=new_remaining, mineral_deposit_m=new_mineral_deposit
                ),
            )
    return events
