import numpy as np

from app import elevation_lines, initial_rock_ages
from app.world import generate_world


def test_plate_age_samples_match_the_configured_distributions():
    rng = np.random.default_rng(744)
    samples = [initial_rock_ages.sample_plate_ages(rng) for _ in range(5000)]
    oceanic = np.array([pair[0] for pair in samples])
    continental = np.array([pair[1] for pair in samples])

    assert np.all((oceanic >= 0.0) & (oceanic <= initial_rock_ages.OCEANIC_MAX_AGE_YEARS))
    np.testing.assert_allclose(
        np.mean(oceanic), initial_rock_ages.OCEANIC_MEAN_AGE_YEARS, atol=2e6
    )

    in_archean = (continental >= 2.5e9) & (continental <= 4.0e9)
    in_proterozoic = (continental >= 541e6) & (continental < 2.5e9)
    in_phanerozoic = continental < 541e6
    shares = np.array([in_archean.mean(), in_proterozoic.mean(), in_phanerozoic.mean()])
    np.testing.assert_allclose(shares, initial_rock_ages.CONTINENTAL_ERA_PROBABILITIES, atol=0.03)


def test_exact_zero_and_one_year_draws_stay_distinct_from_time_zero_and_sentinel():
    ages = [initial_rock_ages._non_sentinel_age(age) for age in (0.0, 1.0)]

    assert all(age > 0.0 for age in ages)
    assert all(-age < 0.0 for age in ages)
    assert all(-age != -1.0 for age in ages)


def test_generation_stamps_each_crust_type_and_records_both_plate_ages():
    world = generate_world(seed=744, num_plates=6, node_density=0.5)

    assert set(world.initial_rock_ages) == {plate.plate_id for plate in world.plates}
    for plate in world.plates:
        record = world.initial_rock_ages[plate.plate_id]
        assert record["primary_crust_type"] == plate.crust_type
        assert 0.0 <= record["oceanic_age_years"] <= initial_rock_ages.OCEANIC_MAX_AGE_YEARS
        assert 0.0 <= record["continental_age_years"] <= 4.0e9

        codes = plate.collect("crust_type_code")
        is_continental = elevation_lines.effective_is_continental_from_codes(
            codes, plate.crust_type == "continental"
        )
        created = plate.collect("node_created_years")
        np.testing.assert_array_equal(
            created[~is_continental], -record["oceanic_age_years"]
        )
        np.testing.assert_array_equal(
            created[is_continental], -record["continental_age_years"]
        )
        assert np.all(created != -1.0)


def test_same_seed_keeps_the_same_initial_rock_age_records():
    first = generate_world(seed=744, num_plates=6, node_density=0.5)
    second = generate_world(seed=744, num_plates=6, node_density=0.5)

    assert first.initial_rock_ages == second.initial_rock_ages
