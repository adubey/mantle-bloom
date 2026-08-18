"""Boundary point reassignment.

Ordinary per-step boundary evolution (boundary.py) only ever grows or shrinks a line's two
theta-direction ends -- it never looks at whether an *interior* node still actually belongs
to the plate that carries it. After enough shearing along a transform boundary, or a slow
rotational drift, a node can end up geometrically embedded in a neighboring plate's own
territory while still being carried by its original plate's line. This pass finds those and
hands them over.

For every node P (across every plate, not just its own), find its NEIGHBOR_COUNT nearest
neighbors globally. If at least MIN_FOREIGN_NEIGHBORS of them belong to the same other plate
L, P is misplaced: it's moved onto whichever of L's lines is closest to it (see
_find_target_line -- almost always one its own matching neighbors already sit on, so no full
per-plate line search is needed in the common case), its local phi snapped to exactly that
line's phi (an ElevationLine is only ever defined at one fixed phi, so "moved onto a line"
and "phi equals that line's phi" are the same thing -- a small shift, since P was already
established to be close to L by the neighbor-majority test), and its elevation interpolated
from its own prior value and its new nearest neighbor on that line.

Deliberately staggered against gaps.fill_gaps (see world.step_world) rather than sharing its
cadence -- both are whole-sphere, cross-plate structural passes; running them the same step
would mean each one's "which plate owns what" picture could already be stale by the time it
acts on it (fill_gaps can spawn or grow plates out from under this pass's neighbor query, and
vice versa).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
from .plates import ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

REASSIGN_INTERVAL_STEPS = 5
NEIGHBOR_COUNT = 4
MIN_FOREIGN_NEIGHBORS = 3


def _gather_nodes(world: "World") -> tuple[np.ndarray, np.ndarray, list[tuple[int, int, int]]]:
    """Every node's world position and elevation, concatenated, alongside same-length
    (plate_id, line_index, node_index) provenance -- unlike plates.collect_all_points, this
    keeps enough detail to actually move a node from one line to another, not just identify
    its owning plate."""
    points_list, elevation_list = [], []
    owners: list[tuple[int, int, int]] = []
    for plate in world.plates:
        for line_index, line in enumerate(plate.lines):
            if len(line.theta) == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            elevation_list.append(line.elevation)
            owners.extend((plate.plate_id, line_index, node_index) for node_index in range(len(line.theta)))
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0), []
    return np.concatenate(points_list, axis=0), np.concatenate(elevation_list, axis=0), owners


def _find_target_line(plate: Plate, candidate_line_indices: set[int], point_world: np.ndarray) -> int | None:
    """Which of `plate`'s lines is closest to `point_world` -- searching only
    `candidate_line_indices` when given (see the optimization note in this module's
    docstring), falling back to every line in the plate only if that's empty."""
    indices = candidate_line_indices if candidate_line_indices else range(len(plate.lines))
    best_index, best_dist = None, np.inf
    for line_index in indices:
        line = plate.lines[line_index]
        if len(line.theta) == 0:
            continue
        dist = float(np.min(geometry.angular_distance(line.world_xyz(plate.frame), point_world)))
        if dist < best_dist:
            best_dist, best_index = dist, line_index
    return best_index


def reassign_misplaced_points(world: "World") -> None:
    """Find every node whose NEIGHBOR_COUNT nearest neighbors are mostly owned by a
    different plate and hand it over -- see this module's docstring for the full algorithm.
    Mutates world.plates' lines in place. Not logged to the UI's event console -- this runs
    every REASSIGN_INTERVAL_STEPS steps and, on a busy boundary, can touch many plates and
    hundreds of points per pass, which flooded the console with little useful signal."""
    points, elevations, owners = _gather_nodes(world)
    n = len(points)
    if n <= NEIGHBOR_COUNT:
        return

    tree = cKDTree(points)
    _, neighbor_idx = tree.query(points, k=NEIGHBOR_COUNT + 1)
    neighbor_idx = neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0

    owner_plate_ids = np.array([o[0] for o in owners])
    neighbor_owner_ids = owner_plate_ids[neighbor_idx]  # (n, NEIGHBOR_COUNT)
    foreign_mask_all = neighbor_owner_ids != owner_plate_ids[:, None]
    # Cheap vectorized pre-filter: most nodes sit deep in their own plate's interior, with
    # every neighbor sharing their own owner, so this keeps the more expensive per-node work
    # below limited to the (usually small) set of actual boundary candidates.
    candidate_indices = np.nonzero(foreign_mask_all.sum(axis=1) >= MIN_FOREIGN_NEIGHBORS)[0]

    plate_by_id = {p.plate_id: p for p in world.plates}

    remove_by_source: dict[tuple[int, int], set[int]] = {}
    add_by_target: dict[tuple[int, int], list[tuple[float, float]]] = {}
    transfer_counts: dict[tuple[int, int], int] = {}

    for i in candidate_indices:
        own_plate_id, own_line_index, own_node_index = owners[i]
        neigh_indices = neighbor_idx[i]
        neigh_owner_ids = neighbor_owner_ids[i]
        foreign_ids = neigh_owner_ids[foreign_mask_all[i]]

        candidate_ids, counts = np.unique(foreign_ids, return_counts=True)
        best = int(np.argmax(counts))
        if counts[best] < MIN_FOREIGN_NEIGHBORS:
            continue
        target_plate_id = int(candidate_ids[best])
        target_plate = plate_by_id.get(target_plate_id)
        if target_plate is None:
            continue

        # Optimization hint: line A is almost always one of these matching neighbors' own
        # lines -- check just those instead of every line in the target plate.
        candidate_line_indices = {
            owners[neigh_indices[k]][1]
            for k in range(len(neigh_indices))
            if neigh_owner_ids[k] == target_plate_id
        }
        target_line_index = _find_target_line(target_plate, candidate_line_indices, points[i])
        if target_line_index is None:
            continue
        target_line = target_plate.lines[target_line_index]

        local_xyz = geometry.to_local(target_plate.frame, points[i])
        _, new_theta = geometry.xyz_to_latlon(local_xyz)
        new_theta = float(new_theta)

        insert_pos = int(np.searchsorted(target_line.theta, new_theta))
        neighbor_candidates = [c for c in (insert_pos - 1, insert_pos) if 0 <= c < len(target_line.theta)]
        nearest_on_line = min(neighbor_candidates, key=lambda c: abs(target_line.theta[c] - new_theta))
        new_elevation = (float(elevations[i]) + float(target_line.elevation[nearest_on_line])) / 2.0

        remove_by_source.setdefault((own_plate_id, own_line_index), set()).add(own_node_index)
        add_by_target.setdefault((target_plate_id, target_line_index), []).append((new_theta, new_elevation))
        pair_key = (own_plate_id, target_plate_id)
        transfer_counts[pair_key] = transfer_counts.get(pair_key, 0) + 1

    if not transfer_counts:
        return

    for (plate_id, line_index), remove_indices in remove_by_source.items():
        plate = plate_by_id[plate_id]
        line = plate.lines[line_index]
        keep_mask = np.ones(len(line.theta), dtype=bool)
        keep_mask[list(remove_indices)] = False
        plate.lines[line_index] = ElevationLine(
            phi=line.phi,
            theta=line.theta[keep_mask],
            elevation=line.elevation[keep_mask],
            channel_depth=line.channel_depth[keep_mask],
            lake_depth=line.lake_depth[keep_mask],
        )

    for (plate_id, line_index), new_nodes in add_by_target.items():
        plate = plate_by_id[plate_id]
        line = plate.lines[line_index]
        new_theta = np.array([t for t, _ in new_nodes])
        new_elevation = np.array([e for _, e in new_nodes])
        combined_theta = np.concatenate([line.theta, new_theta])
        combined_elevation = np.concatenate([line.elevation, new_elevation])
        # A reassigned point's own channel_depth/lake_depth resets to 0 -- this pass only
        # ever touches a small number of boundary-adjacent nodes at a time, so losing a
        # carved channel right at the moment its owning plate changes is an acceptable
        # simplification (unlike line_regrid.py's regularize_line, which runs on nearly
        # every line, every regularize interval -- resetting there would erase rivers
        # constantly, not rarely).
        combined_channel_depth = np.concatenate([line.channel_depth, np.zeros(len(new_nodes))])
        combined_lake_depth = np.concatenate([line.lake_depth, np.zeros(len(new_nodes))])
        order = np.argsort(combined_theta)
        plate.lines[line_index] = ElevationLine(
            phi=line.phi,
            theta=combined_theta[order],
            elevation=combined_elevation[order],
            channel_depth=combined_channel_depth[order],
            lake_depth=combined_lake_depth[order],
        )
