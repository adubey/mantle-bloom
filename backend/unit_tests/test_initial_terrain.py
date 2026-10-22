import numpy as np
import pytest

from app import initial_terrain
from app.world import World

from .quad_fixtures import quad_plate


def test_lowland_spherical_relief_tapers_at_coasts_and_highlands_without_flooding_land():
    elevation = np.array([0.0, 1.0, 200.0, 400.0, 1_400.0, 3_200.0, 5_000.0, 6_000.0])
    noise_sample = np.array([-1.0, -1.0, -1.0, -10.0, 1.0, 1.0, 1.0, 1.0])
    eligible = np.ones(len(elevation), dtype=bool)
    eligible[0] = eligible[-1] = False

    gain = initial_terrain._lowland_relief_gain(elevation, 0.0, noise_sample, eligible)

    assert gain[0] == 0.0  # Sea-level cells are not eligible.
    assert elevation[1] + gain[1] > 0.0  # Negative noise cannot flood the cell.
    assert gain[2] == pytest.approx(-LOWLAND_COAST_EXPECTED_M)
    assert gain[3] == pytest.approx(-400.0 + 1e-6)  # The hard land-mask bound is enforced.
    assert gain[4] == pytest.approx(initial_terrain.LOWLAND_RELIEF_AMPLITUDE_M)
    assert 0.0 < gain[5] < initial_terrain.LOWLAND_RELIEF_AMPLITUDE_M
    assert gain[6] == 0.0
    assert gain[7] == 0.0
    assert np.all(elevation[eligible] + gain[eligible] > 0.0)


LOWLAND_COAST_EXPECTED_M = initial_terrain.LOWLAND_RELIEF_AMPLITUDE_M * (0.5**2 * (3.0 - 2.0 * 0.5))


def _lowland_world(seed: int) -> World:
    plate = quad_plate(
        0,
        "continental",
        columns=range(8),
        elevation=np.array([100.0, 500.0, 1_000.0, 1_400.0, 2_500.0, 5_000.0, -200.0, 600.0]),
    )
    return World(seed=seed, plates=[plate])


def test_lowland_spherical_relief_is_seed_stable_bounded_and_only_changes_lowland_land():
    first = _lowland_world(31)
    second = _lowland_world(31)
    other_seed = _lowland_world(32)
    before = first.plates[0].collect("elevation").copy()

    initial_terrain.apply_lowland_spherical_relief(first)
    initial_terrain.apply_lowland_spherical_relief(second)
    initial_terrain.apply_lowland_spherical_relief(other_seed)

    changed = first.plates[0].collect("elevation")
    same_seed = second.plates[0].collect("elevation")
    different_seed = other_seed.plates[0].collect("elevation")
    assert np.array_equal(changed, same_seed)
    assert not np.array_equal(changed, different_seed)
    assert changed[6] == before[6]  # Ocean is untouched.
    assert changed[5] == before[5]  # The taper is zero at 5 km.
    assert np.any(changed != before)
    assert np.array_equal(changed > 0.0, before > 0.0)  # The initial coastline is fixed.
    assert np.max(np.abs(changed[:6] - before[:6])) <= initial_terrain.LOWLAND_RELIEF_AMPLITUDE_M



def test_river_bootstrap_cuts_through_a_wet_pit_and_keeps_the_notch(monkeypatch):
    # The zero-duration bootstrap had no time to carve or fill, so every pit ended a river and
    # most generated land drained to closed sinks (#218). With the drainage history it breaches
    # a shallow pit in a wet climate, and the notch persists for the first step.
    from app.hydroclimate import climate, hydrology

    monkeypatch.setattr(hydrology, "BREACH_DEPRESSIONS_ENABLED", True)
    elevation = np.full(40, 100.0)
    elevation[0:15] = -100.0
    elevation[28:31] = 80.0
    plate = quad_plate(0, "continental", columns=np.arange(40) - 20, elevation=elevation)
    world = World(seed=0, plates=[plate])
    real = climate.compute_climate

    def wet(*args, **kwargs):
        fields = real(*args, **kwargs)
        fields.precipitation_mm[:] = 2000.0
        fields.air_temperature_c[:] = 10.0
        fields.ocean_temperature_c[:] = 10.0
        return fields

    monkeypatch.setattr(climate, "compute_climate", wet)
    initial_terrain.initialize_river_channels(world)

    notch = plate.collect("breach_notch_depth_m")
    assert notch.max() == pytest.approx(20.0)
    assert plate.collect("lake_depth")[28:31].max() == 0.0
    flow = world.hydrology_cache.flow_target
    node = 30
    for _ in range(40):
        if flow[node] < 0:
            break
        node = flow[node]
    assert world.hydrology_cache.is_ocean[node]
