import numpy as np

from app import stats
from app.plates import ElevationLine, Plate
from app.world import World


def _all_ocean_world() -> World:
    line = ElevationLine(phi=0.0, theta=np.linspace(-np.pi, np.pi, 20, endpoint=False), elevation=np.full(20, -3800.0))
    plate = Plate(plate_id=0, frame=np.eye(3), crust_type="oceanic", lines=[line])
    return World(seed=0, plates=[plate])


def test_compute_stats_land_and_ocean_fractions_sum_to_one():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["land_fraction"] == 0.0
    assert result["ocean_fraction"] == 1.0


def test_compute_stats_land_temperature_none_for_all_ocean_world():
    # No land grid cells at all -- land/air temperature stats must not divide by zero,
    # they should report None instead (same defensive pattern plate-sim's own stats
    # endpoint uses for an empty mask).
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["land_temperature_mean_c"] is None
    assert result["air_temperature_mean_c"] is None
    assert result["ocean_temperature_mean_c"] is not None


def test_compute_stats_elevation_bounds_are_consistent():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["elevation_min_m"] <= result["elevation_mean_m"] <= result["elevation_max_m"]
