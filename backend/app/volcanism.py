"""Volcanic fields: new continental crust forming where plates are separating.

**Detection, run in the clean-up phase (world.py's `steps_since_regularize`-gated block,
alongside gaps.py/line_regrid.py), not every step**, in two passes:

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
from scipy.spatial import cKDTree

from . import gaps, geometry, mantle
from .boundary import MAX_ELEVATION_M, MIN_ELEVATION_M
from .noise import SphereNoise
from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate, base_elevation, build_lines_from_lattice, noise_amplitude

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
# genuinely separate rifts elsewhere on the sphere apart.
VOLCANIC_FIELD_CLUSTER_RADIUS_RAD = 15.0 * TARGET_LINE_SPACING_RAD
# How far from a detected volcano point the new plate's own lattice extends -- same question
# gaps._spawn_plate_from_gap answers for its own new-plate spawning, at its own (tighter)
# scale, since this is about "the new plate's own coverage," not "how far apart do
# same-rift detections still count as one field" (VOLCANIC_FIELD_CLUSTER_RADIUS_RAD above).
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


def _whole_world_median_spacing(world: "World") -> float:
    """Pass 1: every elevation point's own nearest neighbor, whole-world, unrestricted by
    plate -- the median of all these distances, "how far apart elevation points normally
    sit." Returns 0.0 for a world too small to have a meaningful median."""
    points_list = [line.world_xyz(plate.frame) for plate in world.plates for line in plate.lines if len(line.theta) > 0]
    if not points_list:
        return 0.0
    points = np.concatenate(points_list, axis=0)
    if len(points) <= 1:
        return 0.0
    dist, _ = cKDTree(points).query(points, k=2)  # column 0 is always the point itself, at distance 0
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
        tree = cKDTree(points[other_mask])
        dist, idx = tree.query(points[own_mask], k=1)
        own_indices = np.nonzero(own_mask)[0]
        nearest_dist[own_indices] = dist
        nearest_idx[own_indices] = other_indices[idx]
    return nearest_dist, nearest_idx


def _spawn_volcanic_field_plate(world: "World", cluster_points: np.ndarray, rng: np.random.Generator) -> Plate:
    """Builds a brand-new continental Plate covering the lattice around `cluster_points`
    (this pass's own newly-detected volcano points, already clustered by proximity) --
    same `build_lines_from_lattice`-around-a-fresh-seed-frame construction
    gaps._spawn_plate_from_gap uses for an ordinary gap-spawned plate, but continental crust
    (volcanic fields "result in continental plates," regardless of what kind of plates were
    separating) and with every resulting node marked as its own active volcano."""
    centroid = geometry.normalize(cluster_points.mean(axis=0))
    frame = geometry.plate_frame_from_seed(centroid)
    crust_type = "continental"

    base = base_elevation(crust_type)
    amp = noise_amplitude(crust_type)
    noise = SphereNoise(rng, octaves=3, base_freq=2.5)
    tree = cKDTree(cluster_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < VOLCANIC_FIELD_COVERAGE_RADIUS_RAD

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return base + amp * noise.sample(world_pts)

    lines = build_lines_from_lattice(frame, is_owned, elevation_at)
    volcanic_lines = [
        ElevationLine(
            phi=line.phi,
            theta=line.theta,
            elevation=line.elevation,
            is_volcano=np.ones(len(line.theta), dtype=bool),
            volcano_active_years_remaining=rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=len(line.theta)),
        )
        for line in lines
    ]

    plate = Plate(plate_id=world.next_plate_id, frame=frame, crust_type=crust_type, lines=volcanic_lines)
    world.next_plate_id += 1

    points, _ = plate.all_points_and_elevation()
    if len(points) > 0:
        velocities = mantle.flow_at(points, world.mantle_centers)
        plate.omega = mantle.clamp_rate(mantle.fit_euler_pole(points, velocities))
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

    labels = gaps.cluster_points(new_volcano_points, VOLCANIC_FIELD_CLUSTER_RADIUS_RAD)
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


def _volcano_fraction(plate: Plate) -> float:
    total = plate.node_count()
    if total == 0:
        return 0.0
    return sum(int(line.is_volcano.sum()) for line in plate.lines) / total


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
            if len(line.theta) == 0 or not np.any(line.is_volcano):
                continue
            active_mask = line.is_volcano & (line.volcano_active_years_remaining > 0)
            if not np.any(active_mask):
                continue

            active_years_this_step = np.minimum(years, line.volcano_active_years_remaining)
            p_erupt = 1.0 - np.exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1_000_000.0)
            rng = np.random.default_rng((world.seed, round(world.elapsed_years), plate.plate_id, line_index))
            erupts = active_mask & (rng.random(len(line.theta)) < p_erupt)

            new_elevation = line.elevation.copy()
            new_elevation[erupts] += ERUPTION_ELEVATION_M
            new_elevation = np.clip(new_elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            new_remaining = np.clip(line.volcano_active_years_remaining - years, 0.0, None)

            plate.lines[line_index] = ElevationLine(
                phi=line.phi,
                theta=line.theta,
                elevation=new_elevation,
                channel_depth=line.channel_depth,
                lake_depth=line.lake_depth,
                glacier_depth=line.glacier_depth,
                is_volcano=line.is_volcano,
                volcano_active_years_remaining=new_remaining,
            )
    return events
