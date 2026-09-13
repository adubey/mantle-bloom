import numpy as np
import pytest

from app import lithosphere, stats
from app.elevation_lines import line_spacing_rad
from app.plates import ElevationLine, PlateWithLines
from app.world import World


def _all_ocean_world() -> World:
    line = ElevationLine(phi=0.0, theta=np.linspace(-np.pi, np.pi, 20, endpoint=False), elevation=np.full(20, -3800.0))
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="oceanic", lines=[line])
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


def test_compute_stats_total_land_area_and_continental_volume_are_zero_for_all_ocean_world():
    # Read straight off world.plates (see stats.py's own docstring), not the climate grid --
    # an all-ocean world has no land nodes and no continental plates at all, so both must be
    # exactly zero (a running total, unlike the None-for-empty-domain convention the
    # climate-grid stats above use).
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["total_land_area_km2"] == 0.0
    assert result["total_continental_crust_volume_km3"] == 0.0


def test_compute_stats_total_continental_crust_volume_matches_hand_computed_sum():
    # One continental plate, every node holding a known Hc -- total volume should be exactly
    # node_area_m2 * sum(Hc), converted to km^3, regardless of how that Hc happens to be
    # distributed across nodes (see _total_land_area_and_continental_volume's own docstring).
    n = 20
    hc = np.linspace(20_000.0, 40_000.0, n)
    line = ElevationLine(
        phi=0.0, theta=np.linspace(-np.pi, np.pi, n, endpoint=False),
        elevation=np.full(n, 500.0), crustal_thickness_m=hc,
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])
    result = stats.compute_stats(world)

    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    expected_km3 = float(hc.sum()) * area_m2 / 1.0e9
    assert result["total_continental_crust_volume_km3"] == pytest.approx(expected_km3)
    # Every node here sits above the default sea_level_m=0.0, so the whole plate is land.
    assert result["total_land_area_km2"] == pytest.approx(n * area_m2 / 1.0e6)


def test_compute_stats_biome_land_fraction_excludes_ocean_and_sums_to_one():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    # No land cells at all -- every biome fraction should be omitted, not reported as a
    # permanent (and misleading) 0.0/0.
    assert result["biome_land_fraction"] == {}


def test_compute_stats_biome_ocean_fraction_excludes_land_and_sums_to_one():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    # Every cell is ocean here -- the pelagic classes should cover 100% of it, and no Köppen
    # land class should appear as a key.
    assert result["biome_ocean_fraction"] != {}
    assert np.isclose(sum(result["biome_ocean_fraction"].values()), 1.0)
    from app import biomes

    assert all(name in biomes.PELAGIC_NAMES for name in result["biome_ocean_fraction"])


def _land_and_ocean_world() -> World:
    theta = np.linspace(-np.pi, np.pi, 40, endpoint=False)
    elevation = np.where(np.abs(theta) < np.pi / 2, 500.0, -3800.0)  # half land, half ocean
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    return World(seed=0, plates=[plate])


def test_compute_stats_biome_land_fraction_sums_to_one_with_land():
    world = _land_and_ocean_world()
    result = stats.compute_stats(world)
    assert result["biome_land_fraction"] != {}
    assert np.isclose(sum(result["biome_land_fraction"].values()), 1.0)
    assert "Ocean" not in result["biome_land_fraction"]
    # The ocean half is classified into pelagic provinces the same way.
    assert np.isclose(sum(result["biome_ocean_fraction"].values()), 1.0)
    assert set(result["biome_land_fraction"]).isdisjoint(result["biome_ocean_fraction"])


def test_compute_stats_biome_land_fraction_reads_the_stored_climate_cache_biome_ids():
    # stats.py no longer runs its own classify_biomes -- it reads ClimateFields.biome_ids,
    # the same stored field compute_climate now computes once and every other biome-consuming
    # caller shares (see climate.py). Recomputing the expected fractions directly from
    # world.climate_cache.biome_ids (populated as a side effect of compute_stats calling
    # compute_climate_cached) should match stats.py's own result exactly.
    from app import biomes

    world = _land_and_ocean_world()
    result = stats.compute_stats(world)

    biome_ids = world.climate_cache.biome_ids
    is_land = ~world.climate_cache.is_ocean
    land_biome_ids = biome_ids[is_land]
    n_land = int(is_land.sum())
    expected = {
        name: float(np.count_nonzero(land_biome_ids == i)) / n_land
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i not in biomes.OCEAN_IDS and n_land > 0
    }
    assert result["biome_land_fraction"] == expected
