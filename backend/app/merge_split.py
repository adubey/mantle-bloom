"""Plate topology changes: consumption, continental-collision merging, and mantle-flow-
driven splitting.

Unlike rotation and boundary evolution, these are rare, discrete, topology-changing events,
so a one-time resample onto a fresh local lattice (merge) or a plane cut through existing
nodes (split) is an acceptable cost here -- the exact, no-resampling guarantee only matters
for routine per-step motion.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial import cKDTree

from . import geometry, mantle
from .boundary import MERGE_THRESHOLD_RAD
from .plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate, build_lines_from_lattice

if TYPE_CHECKING:
    from .world import World

MERGE_CONTACT_DISTANCE_RAD = MERGE_THRESHOLD_RAD
MERGE_MIN_CONTACT_NODES = 4
MERGE_COVERAGE_RADIUS_RAD = 1.2 * TARGET_LINE_SPACING_RAD

SPLIT_MIN_NODES = 300
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


def remove_consumed_plates(world: "World") -> None:
    """A plate whose every elevation node was deleted (fully subducted, see boundary.py)
    simply vanishes -- no special-cased merge algorithm needed."""
    world.plates = [p for p in world.plates if p.node_count() > 0]


def find_continental_collision_pairs(world: "World") -> list[tuple[int, int]]:
    continental = [p for p in world.plates if p.crust_type == "continental" and p.node_count() > 0]
    points = {p.plate_id: p.all_points_and_elevation()[0] for p in continental}

    pairs = []
    for i, a in enumerate(continental):
        for b in continental[i + 1 :]:
            pa, pb = points[a.plate_id], points[b.plate_id]
            if len(pa) == 0 or len(pb) == 0:
                continue
            dist, _ = cKDTree(pb).query(pa)
            if np.sum(dist < MERGE_CONTACT_DISTANCE_RAD) >= MERGE_MIN_CONTACT_NODES:
                pairs.append((a.plate_id, b.plate_id))
    return pairs


def merge_plates(world: "World", id_keep: int, id_absorb: int) -> None:
    """Fuse `id_absorb` into `id_keep`: keep `id_keep`'s frame, resample the union
    footprint from scratch (a one-time nearest-neighbor lookup into the pre-merge combined
    point cloud), and drop `id_absorb`."""
    keep = next(p for p in world.plates if p.plate_id == id_keep)
    absorb = next(p for p in world.plates if p.plate_id == id_absorb)

    keep_pts, keep_elev = keep.all_points_and_elevation()
    absorb_pts, absorb_elev = absorb.all_points_and_elevation()
    old_points = np.concatenate([keep_pts, absorb_pts], axis=0)
    old_elevation = np.concatenate([keep_elev, absorb_elev], axis=0)
    tree = cKDTree(old_points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(world_pts)
        return dist < MERGE_COVERAGE_RADIUS_RAD

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        _, idx = tree.query(world_pts)
        return old_elevation[idx]

    keep.lines = build_lines_from_lattice(keep.frame, is_owned, elevation_at)
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
    if plate.node_count() < SPLIT_MIN_NODES:
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
            lines_a.append(ElevationLine(phi=line.phi, theta=line.theta[side], elevation=line.elevation[side]))
        if np.any(~side):
            lines_b.append(
                ElevationLine(phi=line.phi, theta=line.theta[~side], elevation=line.elevation[~side])
            )

    if sum(len(l.theta) for l in lines_a) < SPLIT_MIN_NODES or sum(len(l.theta) for l in lines_b) < SPLIT_MIN_NODES:
        return None

    new_id = world.next_plate_id
    world.next_plate_id += 1
    plate_a = Plate(
        plate_id=plate.plate_id,
        frame=plate.frame.copy(),
        crust_type=plate.crust_type,
        omega=mantle.clamp_rate(pole_a),
        boundary_local=plate.boundary_local.copy(),
        lines=lines_a,
    )
    plate_b = Plate(
        plate_id=new_id,
        frame=plate.frame.copy(),
        crust_type=plate.crust_type,
        omega=mantle.clamp_rate(pole_b),
        boundary_local=plate.boundary_local.copy(),
        lines=lines_b,
    )
    return plate_a, plate_b


def apply_topology_changes(world: "World") -> None:
    for plate in world.plates:
        plate.age_steps += 1

    remove_consumed_plates(world)

    merged = True
    while merged:
        merged = False
        pairs = find_continental_collision_pairs(world)
        if pairs:
            merge_plates(world, *pairs[0])
            merged = True

    new_plates: list[Plate] = []
    for plate in world.plates:
        split_result = maybe_split_plate(world, plate)
        if split_result is None:
            new_plates.append(plate)
        else:
            new_plates.extend(split_result)
    world.plates = new_plates
