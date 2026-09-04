"""Eustatic sea level (`app.eustasy`): `World.sea_level_m` re-solved each step from a
conserved ocean water volume against the world's hypsometry."""

import numpy as np
import pytest

from app import eustasy
from app.world import generate_world, step_world


def test_total_water_column_is_monotonic_and_exact():
    z = np.array([-4000.0, -1000.0, -200.0, 50.0, 800.0])
    # Only the nodes below the datum contribute, each its own depth.
    assert eustasy.total_water_column_m(z, 0.0) == 4000.0 + 1000.0 + 200.0
    assert eustasy.total_water_column_m(z, -300.0) == 3700.0 + 700.0
    assert eustasy.total_water_column_m(z, -5000.0) == 0.0
    # Strictly increasing in sea level.
    levels = np.linspace(-4000.0, 1000.0, 50)
    w = [eustasy.total_water_column_m(z, h) for h in levels]
    assert np.all(np.diff(w) >= 0.0)


def test_solve_sea_level_inverts_total_water_column():
    rng = np.random.default_rng(0)
    z = np.sort(rng.uniform(-6000.0, 3000.0, 500))
    for target_level in (-3000.0, -500.0, 0.0, 250.0):
        w = eustasy.total_water_column_m(z, target_level)
        solved = eustasy.solve_sea_level(z, w)
        assert abs(solved - target_level) < 1.0


def test_deepening_the_ocean_basins_drops_sea_level():
    """The whole point: move crust down (a spreading basin) and the same water volume no
    longer reaches as high -- sea level falls, handing dry land back as freeboard."""
    world = generate_world(seed=555, num_plates=8, continental_fraction=0.5, node_density=0.3)
    eustasy.update_sea_level(world)
    sea_before = world.sea_level_m

    for plate in world.plates:
        if plate.crust_type != "oceanic":
            continue
        for i, line in enumerate(plate.lines):
            plate.replace_line(i, line.replace(elevation=line.elevation - 800.0))
    eustasy.update_sea_level(world)

    assert world.sea_level_m < sea_before - 50.0
    # Water volume itself is unchanged -- only where it sits.
    assert eustasy.water_column_for_sea_level(world, world.sea_level_m) == \
        pytest.approx(world.ocean_water_column_m, abs=5.0)


def test_generate_world_snapshots_the_budget_and_keeps_sea_level_at_zero():
    world = generate_world(seed=7, num_plates=8, continental_fraction=0.5, node_density=0.3)
    assert world.ocean_water_column_m is not None
    assert world.sea_level_m == 0.0  # the flat starting datum
    # Re-solving an unchanged world reproduces the same level (the solve is the exact inverse).
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m) < 1.0


def test_set_sea_level_via_water_budget_is_conserved_across_a_step():
    world = generate_world(seed=9, num_plates=8, continental_fraction=0.5, node_density=0.3)
    world.simulate_plate_movement = False
    world.simulate_climate_biomes = False

    eustasy.set_sea_level_via_water_budget(world, 120.0)
    assert world.sea_level_m == 120.0
    budget = world.ocean_water_column_m

    step_world(world, years=1_000_000)
    # Nothing moved (both sims off), so the level and the budget both hold.
    assert abs(world.sea_level_m - 120.0) < 1.0
    assert abs(world.ocean_water_column_m - budget) < 1e-6


def _pile_ice_everywhere(world, depth_m):
    for plate in world.plates:
        for i, line in enumerate(plate.lines):
            gd = line.glacier_depth.copy()
            gd[:] = depth_m
            plate.replace_line(i, line.replace(glacier_depth=gd))


def test_trapped_ice_lowers_sea_level_and_total_budget_is_conserved():
    """Water frozen into ice caps/glaciers is debited from the ocean's share of the conserved
    total budget -- so the shoreline falls as the ice grows (glacio-eustasy), while the total
    surface-water budget itself doesn't change."""
    world = generate_world(seed=321, num_plates=8, continental_fraction=0.5, node_density=0.3)
    eustasy.update_sea_level(world)
    sea_before = world.sea_level_m
    total_before = world.ocean_water_column_m
    assert eustasy.trapped_water_column_m(world) == 0.0

    _pile_ice_everywhere(world, 400.0)
    eustasy.update_sea_level(world)

    assert eustasy.trapped_water_column_m(world) > 0.0
    assert world.sea_level_m < sea_before - 50.0  # ice locked up ocean water -> lower stand
    assert world.ocean_water_column_m == total_before  # the *total* budget is untouched
    # The ocean's own column plus what's trapped still sums back to the conserved total.
    assert eustasy.water_column_for_sea_level(world, world.sea_level_m) + eustasy.trapped_water_column_m(world) == \
        pytest.approx(world.ocean_water_column_m, abs=5.0)

    # Melt it all back off -> sea level returns to where it started.
    _pile_ice_everywhere(world, 0.0)
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m - sea_before) < 1.0


def test_set_sea_level_via_water_budget_folds_in_trapped_water():
    world = generate_world(seed=654, num_plates=8, continental_fraction=0.5, node_density=0.3)
    _pile_ice_everywhere(world, 300.0)

    eustasy.set_sea_level_via_water_budget(world, 0.0)
    trapped = eustasy.trapped_water_column_m(world)
    assert trapped > 0.0
    # The stored budget is ocean-at-0 *plus* the currently-trapped water.
    assert world.ocean_water_column_m == pytest.approx(
        eustasy.water_column_for_sea_level(world, 0.0) + trapped, abs=1.0
    )
    # Re-solving against that budget reproduces the requested level (ice unchanged).
    eustasy.update_sea_level(world)
    assert abs(world.sea_level_m) < 1.0
