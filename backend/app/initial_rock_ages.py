"""Starting rock ages, read from each world's layout without using any random stream."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import elevation_lines, torque, trace
from .elevation_lines import line_spacing_rad
from .surface_fields import SURFACE_FIELDS

if TYPE_CHECKING:
    from .world import World

OCEANIC_MEAN_AGE_YEARS = 64.2e6
OCEANIC_MAX_AGE_YEARS = 3.0 * OCEANIC_MEAN_AGE_YEARS

# Cumulative area fraction, youngest first, against age. Geological boundaries are ages
# before present, each era is uniform within itself as a deliberately simple proxy, and the
# era shares are the regional Laurentia split in #744.
CONTINENTAL_AREA_FRACTIONS = (0.0, 0.24, 0.67, 1.0)
CONTINENTAL_AGES_YEARS = (0.0, 541e6, 2.5e9, 4.0e9)


@trace.traced("Dating the initial rocks", trace.SEEDING, user=True)
def assign_initial_rock_ages(
    world: World, omegas: dict[int, np.ndarray], *, undated_only: bool = False
) -> dict[int, dict[str, float | str]]:
    """Stamp starting cells with ages that vary across each plate (#747).

    Ocean ages with spreading time from its own plate's ridges and continental crust with
    distance from its cratons. Area-weighted rank then maps each crust type onto its Earth-like distribution,
    so the world's age shape is fixed while its pattern follows this layout. Stored ages are
    negative simulation years: ``-1`` remains the pre-tracking sentinel for legacy worlds.

    ``omegas`` is each plate's angular velocity, read for the ridges. ``undated_only`` is for
    dating an existing save (#747): only sentinel cells change and take part in the ranking,
    and a plate with none is left out of the record.
    """
    undated = SURFACE_FIELDS["node_created_years"].sentinel
    cells = [_Cells.of(plate) for plate in world.plates]
    stamped = [
        c.created == undated if undated_only else np.ones(len(c.created), dtype=bool) for c in cells
    ]
    points = np.concatenate([c.points for c in cells])
    areas = np.concatenate([c.areas for c in cells])
    continental = np.concatenate([c.continental for c in cells])
    target = np.concatenate(stamped)

    ages = np.full(len(points), np.nan)
    ocean = target & ~continental
    if np.any(ocean):
        spreading = np.concatenate(_spreading_times(world, omegas, cells))
        ages[ocean] = _oceanic_age_years(_rank(spreading[ocean], areas[ocean]))
    land = target & continental
    if np.any(land):
        craton = np.concatenate([c.craton for c in cells])
        ages[land] = _continental_age_years(_rank(_youth(points, continental, craton, land), areas[land]))
    ages[ages == 1.0] = np.nextafter(1.0, np.inf)  # Stored as -age; -1.0 is the sentinel.

    records: dict[int, dict[str, float | str]] = {}
    start = 0
    for plate, plate_cells, mask in zip(world.plates, cells, stamped):
        stop = start + len(mask)
        if np.any(mask):
            plate_ages = ages[start:stop][mask]
            plate.set_fields_on_plate(node_created_years=np.where(mask, -ages[start:stop], plate_cells.created))
            records[plate.plate_id] = {
                "primary_crust_type": plate.crust_type,
                "youngest_age_years": float(plate_ages.min()),
                "mean_age_years": float(plate_ages.mean()),
                "oldest_age_years": float(plate_ages.max()),
            }
        start = stop
    return records


class _Cells:
    """One plate's cells in world coordinates, as the ranking reads them."""

    def __init__(self, points, areas, created, continental, craton):
        self.points, self.areas, self.created = points, areas, created
        self.continental, self.craton = continental, craton

    @classmethod
    def of(cls, plate) -> _Cells:
        points, _ = plate.all_points_and_elevation()
        continental = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        return cls(points, plate.node_areas_m2(), plate.collect("node_created_years"), continental,
                   plate.collect("craton_crust_m") > 0.0)


def _ridges(world: World, omegas: dict[int, np.ndarray]) -> dict[int, tuple[np.ndarray, float]]:
    """Each plate's oceanic cells at opening boundaries, with its mean half-spreading rate."""
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    ridges = {}
    for plate in world.plates:
        neighbours = plate.get_neighbours(world.plates, threshold_rad=reach)
        if plate.node_count() == 0 or not neighbours:
            continue
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach)
        closing = torque.closing_rate_under(omegas, plate, neighbours, inputs)
        oceanic = ~elevation_lines.effective_is_continental_from_codes(
            inputs.own_crust_type_codes, plate.crust_type == "continental"
        )
        ridge = (inputs.dist_to_neighbor <= reach) & oceanic & (closing < -torque.BOUNDARY_TRANSFORM_RATE_M_PER_S)
        if np.any(ridge):
            ridges[plate.plate_id] = (inputs.own_points[ridge], float(np.mean(-closing[ridge])) / 2.0)
    return ridges


def _spreading_times(world: World, omegas: dict[int, np.ndarray], cells: list[_Cells]) -> list[np.ndarray]:
    """Each plate's cells' distance from its own ridges over its own spreading rate.

    Measured per plate because the age field must jump at a trench, where two plates with
    different histories meet: one nearest-ridge distance across the world is continuous, and
    tied converging contacts at a median 3 Myr (#747). A plate with no ridge of its own reads
    the world's ridges at their median rate, and a world with none ages from one arbitrary
    point.
    """
    ridges = _ridges(world, omegas)
    shared = np.concatenate([points for points, _rate in ridges.values()]) if ridges else np.zeros((1, 3)) + (1.0, 0.0, 0.0)
    shared_rate = float(np.median([rate for _points, rate in ridges.values()])) if ridges else 1.0
    times = []
    for plate, plate_cells in zip(world.plates, cells):
        points, rate = ridges.get(plate.plate_id, (shared, shared_rate))
        times.append(cKDTree(points).query(plate_cells.points)[0] / rate)
    return times


def _youth(points: np.ndarray, continental: np.ndarray, craton: np.ndarray, land: np.ndarray) -> np.ndarray:
    """How young each `land` cell should rank: far from cratons, or nearer the coast if none."""
    if np.any(craton):
        return -cKDTree(points[craton]).query(points[land])[0]
    if np.all(continental):
        return np.zeros(int(np.count_nonzero(land)))
    return -cKDTree(points[~continental]).query(points[land])[0]


def _rank(key: np.ndarray, areas: np.ndarray) -> np.ndarray:
    """Each cell's area-weighted quantile by `key`, strictly between 0 and 1."""
    order = np.argsort(key, kind="stable")
    ordered_areas = areas[order]
    ranks = np.empty(len(key))
    ranks[order] = (np.cumsum(ordered_areas) - 0.5 * ordered_areas) / ordered_areas.sum()
    return ranks


def _oceanic_age_years(rank: np.ndarray) -> np.ndarray:
    """The inverse CDF of a triangular age distribution from 0 to the maximum, mode zero."""
    return OCEANIC_MAX_AGE_YEARS * (1.0 - np.sqrt(1.0 - rank))


def _continental_age_years(rank: np.ndarray) -> np.ndarray:
    return np.interp(rank, CONTINENTAL_AREA_FRACTIONS, CONTINENTAL_AGES_YEARS)
