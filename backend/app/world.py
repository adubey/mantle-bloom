"""Top-level world state: a collection of plates, a mantle-flow field, and elapsed time."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import boundary, geometry, mantle
from .plates import Plate, generate_plates

DEFAULT_NUM_PLATES = 12
DEFAULT_MANTLE_CENTERS = 8


@dataclass
class World:
    seed: int
    plates: list[Plate] = field(default_factory=list)
    mantle_centers: list[mantle.ConvectionCenter] = field(default_factory=list)
    elapsed_years: float = 0.0


def _plate_sample_points(plate: Plate) -> np.ndarray:
    """World positions used to sample the mantle flow field for fitting this plate's Euler
    pole: every elevation-line node plus the boundary loop -- its current footprint, for
    free (no separate sampling grid needed)."""
    pieces = [line.world_xyz(plate.frame) for line in plate.lines]
    pieces.append(plate.boundary_world())
    if not pieces:
        return np.zeros((0, 3))
    return np.concatenate(pieces, axis=0)


def _update_plate_omega(
    plate: Plate, centers: list[mantle.ConvectionCenter], damping: float | None
) -> None:
    sample_points = _plate_sample_points(plate)
    if len(sample_points) == 0:
        return
    velocities = mantle.flow_at(sample_points, centers)
    target_omega = mantle.fit_euler_pole(sample_points, velocities)
    if damping is None:
        new_omega = target_omega
    else:
        new_omega = plate.omega + damping * (target_omega - plate.omega)
    plate.omega = mantle.clamp_rate(new_omega)


def generate_world(
    seed: int,
    num_plates: int = DEFAULT_NUM_PLATES,
    num_mantle_centers: int = DEFAULT_MANTLE_CENTERS,
) -> World:
    plates = generate_plates(seed, num_plates=num_plates)
    # Separate RNG stream so changing num_mantle_centers doesn't reshuffle plate layout.
    mantle_rng = np.random.default_rng(seed + 1)
    mantle_centers = mantle.generate_convection_centers(mantle_rng, n_centers=num_mantle_centers)

    world = World(seed=seed, plates=plates, mantle_centers=mantle_centers, elapsed_years=0.0)
    for plate in world.plates:
        _update_plate_omega(plate, world.mantle_centers, damping=None)
    return world


def step_world(world: World, years: float) -> None:
    """Advance the world by `years`: refit each plate's Euler pole from the mantle flow
    field (damped toward the new target), rotate the plate rigidly by that pole for
    `years` (exact for every carried point, no resampling), then let boundaries evolve --
    uplift/trench/ridge/rift elevation deltas and line growth/shrinkage where plates are
    now close to each other."""
    for plate in world.plates:
        _update_plate_omega(plate, world.mantle_centers, damping=mantle.VELOCITY_DAMPING)
        increment = geometry.rotation_matrix_from_omega(plate.omega, years)
        plate.frame = increment @ plate.frame
    boundary.step_boundaries(world, years)
    world.elapsed_years += years
