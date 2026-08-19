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
    # they should report None instead (the defensive pattern used here for an empty mask).
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["land_temperature_mean_c"] is None
    assert result["air_temperature_mean_c"] is None
    assert result["ocean_temperature_mean_c"] is not None


def test_compute_stats_elevation_is_land_only():
    # No land cells at all -- elevation_* (now land-only, see stats.py's module docstring)
    # must report None the same way land_temperature_mean_c already does, not divide by zero
    # or return a bogus range built from ocean cells.
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["elevation_min_m"] is None
    assert result["elevation_max_m"] is None
    assert result["elevation_mean_m"] is None


def test_compute_stats_ocean_depth_bounds_are_consistent():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["ocean_depth_min_m"] <= result["ocean_depth_mean_m"] <= result["ocean_depth_max_m"]
    # Every cell here is a uniform -3800m ocean floor at the default sea_level_m=0.0, so
    # depth should come out to a uniform +3800m.
    assert result["ocean_depth_min_m"] == result["ocean_depth_max_m"] == 3800.0


def test_compute_stats_biome_land_fraction_excludes_ocean_and_sums_to_one():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    # No land cells at all -- every biome fraction should be omitted, not reported as a
    # permanent (and misleading) 0.0/0.
    assert result["biome_land_fraction"] == {}


def _land_and_ocean_world() -> World:
    theta = np.linspace(-np.pi, np.pi, 40, endpoint=False)
    elevation = np.where(np.abs(theta) < np.pi / 2, 500.0, -3800.0)  # half land, half ocean
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)
    plate = Plate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    return World(seed=0, plates=[plate])


def test_compute_stats_biome_land_fraction_sums_to_one_with_land():
    world = _land_and_ocean_world()
    result = stats.compute_stats(world)
    assert result["biome_land_fraction"] != {}
    assert np.isclose(sum(result["biome_land_fraction"].values()), 1.0)
    assert "Ocean" not in result["biome_land_fraction"]
