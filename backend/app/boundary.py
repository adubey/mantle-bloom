"""Boundary evolution: where a plate's elevation-line nodes are close to another plate's,
classify the local relative motion (convergent/divergent/transform) and apply the
consequences -- uplift, trench deepening, ridge/rift relaxation, and, at each line's two
ends, inserting new nodes (crust created at a divergent boundary) or deleting them (crust
destroyed/folded at a convergent one).

Adjacency is *not* a maintained topological structure -- every step, each plate's nodes are
matched against a fresh k-d tree of every other plate's current nodes. This is what lets
plates evolve independently (see plates.py) while boundaries still behave sensibly: it's
self-healing every step rather than requiring an always-consistent shared-edge structure.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import mantle
from .plates import TARGET_LINE_SPACING_RAD, ElevationLine

if TYPE_CHECKING:
    from .world import World

FAR_THRESHOLD_RAD = 1.6 * TARGET_LINE_SPACING_RAD
EXTEND_THRESHOLD_RAD = 1.3 * TARGET_LINE_SPACING_RAD
MERGE_THRESHOLD_RAD = 0.4 * TARGET_LINE_SPACING_RAD

TRANSFORM_RATE_THRESHOLD = mantle.cm_per_yr_to_rad_per_yr(1.0)

CONVERGENT_MOUNTAIN_RATE_M_PER_MYR = 800.0
CONVERGENT_TRENCH_RATE_M_PER_MYR = 700.0
DIVERGENT_RIDGE_TARGET_M = -1500.0  # new oceanic crust at a mid-ocean ridge
DIVERGENT_RIFT_TARGET_M = -200.0  # new continental crust in a rift valley
DIVERGENT_RELAX_RATE_PER_MYR = 0.5

MIN_ELEVATION_M = -11000.0
MAX_ELEVATION_M = 9000.0

# Safety cap on how many nodes a single step can insert at one line end. Not meant to bind
# in practice -- even at MAX_PLATE_RATE with the largest step size the UI offers, the real
# gap is only ever a handful of spacing units (see the comment on _grow_or_shrink_line for
# why a fixed one-node-per-step used to fall far short of that and why it matters).
MAX_EXTEND_NODES_PER_STEP = 200


def _divergent_target(crust_type: str) -> float:
    return DIVERGENT_RIDGE_TARGET_M if crust_type == "oceanic" else DIVERGENT_RIFT_TARGET_M


def _closing_rate(
    points: np.ndarray, self_omega: np.ndarray, neighbor_omega: np.ndarray, neighbor_points: np.ndarray
) -> np.ndarray:
    """Positive = this plate's material is moving toward the neighbor's (convergent) at
    this point; negative = moving apart (divergent)."""
    v_self = np.cross(self_omega, points)
    v_neighbor = np.cross(neighbor_omega, points)
    normal_dir = neighbor_points - points
    norm = np.linalg.norm(normal_dir, axis=-1, keepdims=True)
    safe_norm = np.where(norm < 1e-12, 1.0, norm)
    normal_dir = normal_dir / safe_norm
    return np.sum((v_self - v_neighbor) * normal_dir, axis=-1)


def _grow_or_shrink_line(
    line: ElevationLine,
    dist: np.ndarray,
    closing: np.ndarray,
    crust_type: str,
) -> ElevationLine:
    """Grow or shrink a line's two ends based on this step's boundary classification there.

    Growth inserts as many nodes as it takes to actually close the gap (`dist`), not just
    one -- at a small step size the two are the same thing, but at a large `years` step (the
    UI offers up to 10 Myr/step) a fast-diverging boundary can open by many spacing units in
    a single step, and inserting only one node per step falls further and further behind.
    That leftover, perpetually-reopening gap looked to gaps.py's periodic gap-filling like a
    genuinely new, unclosable gap and kept spawning fresh micro-plates at the same busy
    boundary every interval -- fixing the actual growth rate here is what stops that at the
    source, rather than only reacting to its symptom in gaps.py."""
    theta = line.theta.copy()
    elevation = line.elevation.copy()
    if len(theta) == 0:
        return ElevationLine(phi=line.phi, theta=theta, elevation=elevation)

    dtheta = TARGET_LINE_SPACING_RAD / max(np.cos(line.phi), 1e-3)
    target = _divergent_target(crust_type)

    # High end first so the low-end index (0) is unaffected by any change made here.
    if dist[-1] > EXTEND_THRESHOLD_RAD and closing[-1] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[-1] / TARGET_LINE_SPACING_RAD), 1), MAX_EXTEND_NODES_PER_STEP)
        new_theta = theta[-1] + dtheta * np.arange(1, n_new + 1)
        theta = np.append(theta, new_theta)
        elevation = np.append(elevation, np.full(n_new, target))
    elif dist[-1] < MERGE_THRESHOLD_RAD and closing[-1] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        theta = theta[:-1]
        elevation = elevation[:-1]

    if len(theta) == 0:
        return ElevationLine(phi=line.phi, theta=theta, elevation=elevation)

    if dist[0] > EXTEND_THRESHOLD_RAD and closing[0] < -TRANSFORM_RATE_THRESHOLD:
        n_new = min(max(int(dist[0] / TARGET_LINE_SPACING_RAD), 1), MAX_EXTEND_NODES_PER_STEP)
        new_theta = theta[0] - dtheta * np.arange(n_new, 0, -1)
        theta = np.insert(theta, 0, new_theta)
        elevation = np.insert(elevation, 0, np.full(n_new, target))
    elif dist[0] < MERGE_THRESHOLD_RAD and closing[0] > TRANSFORM_RATE_THRESHOLD and len(theta) > 1:
        theta = theta[1:]
        elevation = elevation[1:]

    return ElevationLine(phi=line.phi, theta=theta, elevation=elevation)


def step_boundaries(world: World, years: float) -> None:
    """Mutates world.plates in place: applies boundary elevation deltas and grows/shrinks
    each line's two ends. Call after rotation, using the plates' new positions."""
    years_myr = years / 1e6

    plate_by_id = {plate.plate_id: plate for plate in world.plates}

    plate_points: dict[int, np.ndarray] = {}
    for plate in world.plates:
        pieces = [line.world_xyz(plate.frame) for line in plate.lines]
        plate_points[plate.plate_id] = (
            np.concatenate(pieces, axis=0) if pieces else np.zeros((0, 3))
        )

    for plate in world.plates:
        other_points = []
        other_owner = []
        for other in world.plates:
            if other.plate_id == plate.plate_id or len(plate_points[other.plate_id]) == 0:
                continue
            other_points.append(plate_points[other.plate_id])
            other_owner.append(np.full(len(plate_points[other.plate_id]), other.plate_id))
        if not other_points or not plate.lines:
            continue
        other_points = np.concatenate(other_points, axis=0)
        other_owner = np.concatenate(other_owner, axis=0)
        tree = cKDTree(other_points)

        new_lines = []
        for line in plate.lines:
            pts = line.world_xyz(plate.frame)
            dist, idx = tree.query(pts)
            neighbor_owner = other_owner[idx]
            neighbor_points = other_points[idx]

            neighbor_omega = np.array([plate_by_id[o].omega for o in neighbor_owner])
            closing = _closing_rate(pts, plate.omega, neighbor_omega, neighbor_points)

            intensity = np.clip(1.0 - dist / FAR_THRESHOLD_RAD, 0.0, 1.0)
            near_boundary = dist < FAR_THRESHOLD_RAD
            convergent = near_boundary & (closing > TRANSFORM_RATE_THRESHOLD)
            divergent = near_boundary & (closing < -TRANSFORM_RATE_THRESHOLD)

            elevation = line.elevation.copy()
            if plate.crust_type == "continental":
                elevation[convergent] += CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * years_myr * intensity[convergent]
            else:
                elevation[convergent] -= CONVERGENT_TRENCH_RATE_M_PER_MYR * years_myr * intensity[convergent]

            target = _divergent_target(plate.crust_type)
            relax_factor = 1.0 - np.exp(-DIVERGENT_RELAX_RATE_PER_MYR * years_myr)
            elevation[divergent] += (target - elevation[divergent]) * relax_factor * intensity[divergent]

            elevation = np.clip(elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            updated_line = ElevationLine(phi=line.phi, theta=line.theta, elevation=elevation)
            grown_line = _grow_or_shrink_line(updated_line, dist, closing, plate.crust_type)
            if len(grown_line.theta) > 0:
                new_lines.append(grown_line)

        plate.lines = new_lines
