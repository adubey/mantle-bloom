"""Whole-world checks of generation's synthetic-history relief (#218), too slow for unit_tests."""

import numpy as np
import pytest

from app import lithosphere, plates
from app import world as world_mod
from app.elevation_lines import effective_is_continental_from_codes

SEED = 447276134


def _generate(**options):
    return world_mod.generate_world(SEED, node_density=0.5, land_fraction=0.29, **options)


@pytest.fixture(scope="module")
def worlds():
    return _generate(cratons_first=True, synthetic_relief=False), _generate(cratons_first=True, synthetic_relief=True)


def _continental(world) -> np.ndarray:
    return np.concatenate(
        [effective_is_continental_from_codes(p.collect("crust_type_code"), p.crust_type == "continental") for p in world.plates]
    )


def test_every_ice_free_column_floats_at_its_elevation(worlds):
    _, world = worlds
    hc = plates.collect_all_crustal_thickness(world.plates)
    hm = plates.collect_all_mantle_lithosphere_thickness(world.plates)
    rho = np.concatenate([lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in world.plates])
    ice = plates.collect_all_glacier_depth(world.plates) > 0.0
    residual = np.abs(plates.collect_all_elevation(world.plates) - lithosphere.isostatic_elevation(hc, hm, rho))

    assert residual[~ice].max() < 1e-6
    assert np.all((hc >= lithosphere.MIN_CRUSTAL_THICKNESS_M) & (hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M))


def test_land_stays_on_continental_crust_at_the_requested_share(worlds):
    control, world = worlds
    land = plates.collect_all_elevation(world.plates) > world.sea_level_m
    control_land = plates.collect_all_elevation(control.plates) > control.sea_level_m

    assert not np.any(land & ~_continental(world))
    assert abs(land.mean() - control_land.mean()) <= 0.01


def test_oceanic_crust_and_plates_are_left_as_generated(worlds):
    control, world = worlds
    oceanic = ~_continental(world)

    assert [len(p.cell_keys) for p in world.plates] == [len(p.cell_keys) for p in control.plates]
    assert np.array_equal(
        plates.collect_all_crustal_thickness(world.plates)[oceanic],
        plates.collect_all_crustal_thickness(control.plates)[oceanic],
    )


def test_cratonic_cores_stay_cratonic_through_their_whole_new_crust(worlds):
    _, world = worlds
    for plate in world.plates:
        craton = plate.collect("craton_crust_m")
        core = craton > 0.0
        assert np.allclose(craton[core], plate.collect("crustal_thickness_m")[core])


def test_synthetic_relief_needs_the_cratons_first_layout():
    plain = _generate(cratons_first=False, synthetic_relief=False)
    flagged = _generate(cratons_first=False, synthetic_relief=True)

    assert np.array_equal(plates.collect_all_elevation(plain.plates), plates.collect_all_elevation(flagged.plates))
