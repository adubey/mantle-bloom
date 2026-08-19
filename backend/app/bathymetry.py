"""Submerged continental crust: continental shelf vs deep water.

Newly generated continental crust, and continental crust pulled underwater by rifting (see
boundary.py) or born underwater by generation-time noise, has nothing else in this codebase
pulling it toward any particular depth once it's submerged. In reality, submerged continental
crust close to a coastline sits on the continental shelf (shallow), while farther out it's
genuinely deep water -- distinct from, though shallower than, oceanic crust's own abyssal
depth. This module is a slow background relaxation (the same exponential-toward-a-target
style boundary.py's divergent relaxation already uses) pulling every submerged continental
node toward whichever of those two targets its distance to the nearest land node implies.

Deliberately continental-only: oceanic crust's own average depth already comes from its
generation-time baseline (plates.BASE_OCEANIC_M = -3800) and nothing erodes or otherwise
drifts it away from that on its own (erosion.py explicitly excludes ocean nodes from both its
sources), so it doesn't need a parallel correction here.

"Land" is computed globally (elevation > 0, any plate) -- this is a geographic
distance-to-coastline question, not a plate-boundary one, so a submerged node's own plate
membership is irrelevant to how close it is to *some* coastline.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from .plates import PLANET_RADIUS_KM, Plate

if TYPE_CHECKING:
    from .world import World

# Real continental shelves are shallow and comparatively narrow before the "shelf break"
# drops off toward deep water -- SHELF_TARGET_M has no exact figure to port (none was
# given), picked as a reasonable shelf-like depth, clearly shallower than open ocean.
SHELF_RANGE_KM = 200.0
SHELF_RANGE_RAD = SHELF_RANGE_KM / PLANET_RADIUS_KM
SHELF_TARGET_M = -100.0
DEEP_CONTINENTAL_TARGET_M = -3000.0
# Slower than boundary.DIVERGENT_RELAX_RATE_PER_MYR (0.5) -- this is a passive background
# equilibration of already-submerged, non-actively-deforming crust, not an active tectonic
# process, so it's given a gentler pace. A starting point, not a derived constant.
BATHYMETRY_RELAX_RATE_PER_MYR = 0.3

MIN_ELEVATION_M = -11000.0
MAX_ELEVATION_M = 9000.0


def _gather_nodes(world: "World") -> tuple[np.ndarray, np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position, elevation, and whether its plate is continental,
    concatenated, alongside (plate, line_index, start, end) references -- same shape as
    erosion.py's own _gather_nodes, for the same reason (slicing a flat per-node result
    straight back onto each line's own contiguous elevation range)."""
    points_list, elev_list, continental_list = [], [], []
    line_refs: list[tuple[Plate, int, int, int]] = []
    offset = 0
    for plate in world.plates:
        is_continental = plate.crust_type == "continental"
        for line_index, line in enumerate(plate.lines):
            n = len(line.theta)
            if n == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            elev_list.append(line.elevation)
            continental_list.append(np.full(n, is_continental))
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0), np.zeros(0, dtype=bool), []
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elev_list, axis=0),
        np.concatenate(continental_list, axis=0),
        line_refs,
    )


def apply_bathymetry(world: "World", years: float) -> None:
    """Relaxes every submerged (elevation <= 0) continental node toward SHELF_TARGET_M if
    within SHELF_RANGE_RAD of the nearest land node (elevation > 0, any plate), or
    DEEP_CONTINENTAL_TARGET_M otherwise. Mutates world.plates' line elevations in place;
    never touches node positions or line topology, so this can't interact with line
    regularization or point reassignment at all."""
    points, elevation, is_continental, line_refs = _gather_nodes(world)
    if len(points) == 0:
        return

    is_land = elevation > 0.0
    if not np.any(is_land):
        return  # no coastline anywhere to measure distance from

    submerged_indices = np.nonzero(is_continental & ~is_land)[0]
    if len(submerged_indices) == 0:
        return

    land_tree = cKDTree(points[is_land])
    dist_to_land, _ = land_tree.query(points[submerged_indices])
    target = np.where(dist_to_land <= SHELF_RANGE_RAD, SHELF_TARGET_M, DEEP_CONTINENTAL_TARGET_M)

    years_myr = years / 1_000_000.0
    relax_factor = 1.0 - np.exp(-BATHYMETRY_RELAX_RATE_PER_MYR * years_myr)

    new_elevation = elevation.copy()
    new_elevation[submerged_indices] += (target - elevation[submerged_indices]) * relax_factor
    new_elevation = np.clip(new_elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)

    # theta (and therefore every other parallel array's shape) is never touched here -- only
    # elevation changes -- so dataclasses.replace copies every other field (channel_depth/
    # channel_width/lake_depth/glacier_depth/is_volcano/volcano_active_years_remaining, and
    # whatever else gets added later) from the existing line automatically, rather than
    # needing to name each one explicitly here. See plates.ElevationLine's own docstring for
    # the bug this exact pattern replaces (this site used to silently wipe is_volcano to
    # False every step).
    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        plate.lines[line_index] = dataclasses.replace(line, elevation=new_elevation[start:end])
