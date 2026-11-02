import copy
import functools

import numpy as np
from scipy.spatial import cKDTree

from app import elevation_lines, initial_rock_ages
from app.world import generate_world


@functools.cache
def _world():
    return generate_world(seed=744, num_plates=6, node_density=0.5)


def _oceanic_cells(world):
    rows = []
    for plate in world.plates:
        continental = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        points, _ = plate.all_points_and_elevation()
        rows.append((points[~continental], plate.node_areas_m2()[~continental], -plate.collect("node_created_years")[~continental]))
    return tuple(np.concatenate(columns) for columns in zip(*rows))


def _omegas(world):
    from app import torque

    return torque.solve_plates_quasistatic(world)


def test_rank_is_area_weighted_and_stays_inside_the_unit_interval():
    ranks = initial_rock_ages._rank(np.array([3.0, 1.0, 2.0]), np.array([1.0, 2.0, 1.0]))
    np.testing.assert_allclose(ranks, [0.875, 0.25, 0.625])


def test_oceanic_ranks_reproduce_the_triangular_distribution():
    rank = (np.arange(100000) + 0.5) / 100000
    ages = initial_rock_ages._oceanic_age_years(rank)
    np.testing.assert_allclose(ages.mean(), initial_rock_ages.OCEANIC_MEAN_AGE_YEARS, rtol=1e-3)
    assert ages.max() < initial_rock_ages.OCEANIC_MAX_AGE_YEARS
    assert np.all(np.diff(ages) > 0.0)


def test_continental_ranks_reproduce_the_era_shares():
    ages = initial_rock_ages._continental_age_years(np.array([0.12, 0.45, 0.8]))
    assert ages[0] < 541e6 <= ages[1] < 2.5e9 <= ages[2] < 4.0e9


def test_every_starting_cell_is_dated_and_ages_vary_within_each_plate():
    world = _world()
    assert set(world.initial_rock_ages) == {plate.plate_id for plate in world.plates}
    for plate in world.plates:
        created = plate.collect("node_created_years")
        record = world.initial_rock_ages[plate.plate_id]
        assert np.all(created < 0.0)
        assert len(np.unique(created)) > 1
        assert record["primary_crust_type"] == plate.crust_type
        assert record["youngest_age_years"] == -created.max()
        assert record["oldest_age_years"] == -created.min()


def test_the_worlds_ocean_has_the_earth_mean_age_by_area():
    _points, areas, ages = _oceanic_cells(_world())
    np.testing.assert_allclose(np.average(ages, weights=areas), initial_rock_ages.OCEANIC_MEAN_AGE_YEARS, rtol=0.01)


def test_ocean_ages_rise_away_from_each_plates_own_ridges():
    world = copy.deepcopy(_world())
    ridges = initial_rock_ages._ridges(world, _omegas(world))
    assert ridges
    for plate in world.plates:
        if plate.plate_id not in ridges or plate.crust_type != "oceanic":
            continue
        points, _elevation = plate.all_points_and_elevation()
        ages = -plate.collect("node_created_years")
        distance = cKDTree(ridges[plate.plate_id][0]).query(points)[0]
        order = np.argsort(distance, kind="stable")
        assert np.corrcoef(distance[order], ages[order])[0, 1] > 0.9


def test_spreading_time_is_distance_over_the_plates_own_rate():
    # Spreading time, not distance, is what keeps two plates from tying at a trench (#747).
    world = copy.deepcopy(_world())
    omegas = _omegas(world)
    cells = [initial_rock_ages._Cells.of(plate) for plate in world.plates]
    ridges = initial_rock_ages._ridges(world, omegas)
    plate_id, (points, rate) = next(iter(ridges.items()))
    index = [plate.plate_id for plate in world.plates].index(plate_id)
    times = initial_rock_ages._spreading_times(world, omegas, cells)[index]
    distance = cKDTree(points).query(cells[index].points)[0]
    np.testing.assert_allclose(times, distance / rate)


def test_continental_crust_is_oldest_at_its_cratons():
    points = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    continental = np.array([True, True, False])
    craton = np.array([True, False, False])
    youth = initial_rock_ages._youth(points, continental, craton, continental)
    assert youth[0] > youth[1]


def test_without_cratons_the_interior_is_oldest():
    points = np.array([[1.0, 0.0, 0.0], [0.8, 0.6, 0.0], [0.0, 1.0, 0.0]])
    continental = np.array([True, True, False])
    youth = initial_rock_ages._youth(points, continental, np.zeros(3, dtype=bool), continental)
    assert youth[0] < youth[1]


def test_a_stamped_age_never_lands_on_the_sentinel():
    world = copy.deepcopy(_world())
    initial_rock_ages.assign_initial_rock_ages(world, _omegas(world))
    for plate in world.plates:
        assert not np.any(plate.collect("node_created_years") == -1.0)


def test_the_same_seed_keeps_the_same_ages():
    second = generate_world(seed=744, num_plates=6, node_density=0.5)
    assert second.initial_rock_ages == _world().initial_rock_ages
