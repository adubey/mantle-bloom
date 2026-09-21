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
    assert result["land_near_max_elevation_fraction"] is None


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


def test_compute_stats_land_and_ocean_fractions_count_a_lake_as_water():
    # A landlocked lake sits on land elevation-wise (never below world.sea_level_m, and never
    # part of the connected ocean -- that's exactly what makes it a lake), so `is_ocean` alone
    # would miss it entirely. land_fraction/ocean_fraction ("Land"/"Water" in the frontend's
    # Stats panel) must still count it as water, not land -- see stats.py's own docstring.
    n = 40
    theta = np.linspace(-np.pi, np.pi, n, endpoint=False)
    elevation = np.full(n, 500.0)  # all land, well above the default sea_level_m=0.0
    lake_depth = np.where(np.abs(theta) < np.pi / 4, 20.0, 0.0)  # a lake over a quarter of it
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation, lake_depth=lake_depth)
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    result = stats.compute_stats(world)
    lake_fraction = float(np.count_nonzero(lake_depth > 0.0)) / n
    assert result["ocean_fraction"] == pytest.approx(lake_fraction, abs=0.05)
    assert result["land_fraction"] == pytest.approx(1.0 - lake_fraction, abs=0.05)
    assert result["land_fraction"] + result["ocean_fraction"] == pytest.approx(1.0)


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


def test_compute_stats_land_fraction_node_matches_raw_elevation_count():
    # land_fraction_node (see stats.py's own docstring, GitHub issue #121) is a plain
    # elevation > sea_level_m node count/total -- unlike land_fraction, not connectivity-aware
    # and not resampled from any hydrology cache. Compare against an independent recount
    # straight off world.plates (line regularization means node count isn't necessarily the
    # ElevationLine's original theta length, so don't assume exactly half).
    world = _land_and_ocean_world()
    result = stats.compute_stats(world)

    land_nodes = 0
    total_nodes = 0
    for plate in world.plates:
        _, elevation = plate.all_points_and_elevation()
        land_nodes += int(np.count_nonzero(elevation > world.sea_level_m))
        total_nodes += len(elevation)
    assert result["land_fraction_node"] == pytest.approx(land_nodes / total_nodes)
    assert 0.0 < result["land_fraction_node"] < 1.0  # genuinely a mix, not all-one-or-the-other


def test_compute_stats_land_fraction_node_zero_for_all_ocean_world():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["land_fraction_node"] == 0.0


def test_compute_stats_land_fraction_stale_false_when_hydrology_cache_never_populated():
    # A freshly constructed World has hydrology_cache=None -- land_fraction is computed live
    # off a fresh climate solve in that case (see hydrology.sample_is_ocean's elevation-only
    # fallback), not resampled from a frozen cache, so it must never read as stale.
    world = _land_and_ocean_world()
    assert world.hydrology_cache is None
    result = stats.compute_stats(world)
    assert result["land_fraction_stale"] is False


def test_compute_stats_land_fraction_stale_true_once_climate_is_toggled_off():
    # End-to-end version of test_world_stepping.py's
    # test_hydrology_cache_step_freezes_when_climate_toggled_off: once simulate_climate_biomes
    # is off, world.hydrology_cache stops advancing with world.steps_taken, and
    # compute_stats's land_fraction_stale must flag exactly that (GitHub issue #121).
    from app.world import generate_world, step_world

    world = generate_world(seed=10, num_plates=8)
    step_world(world, years=1_000_000)
    assert stats.compute_stats(world)["land_fraction_stale"] is False

    world.simulate_climate_biomes = False
    step_world(world, years=1_000_000)
    assert stats.compute_stats(world)["land_fraction_stale"] is True


@pytest.mark.parametrize(
    "heights, sea_level, expected",
    [
        ([1000, 950, 949, 100, -1000], 0, 0.5),
        ([1000, 950, 949, 100, -1000], 200, 0.5),
        ([500, 500, 500], 0, 1.0),
        ([0, 0, 0], 0, 1.0),
        ([-100, -105, -106], 0, 2 / 3),
    ],
)
def test_land_near_max_elevation_fraction(monkeypatch, heights, sea_level, expected):
    from types import SimpleNamespace

    height = np.array(heights, dtype=float)
    # Stub the climate grid so the exact 95% boundary is not blurred by resampling.
    # Deliberately stale ocean labels on peaks must be reconciled before counting.
    fields = SimpleNamespace(
        elevation_m=height + sea_level,
        is_ocean=height > 0,
        lake_depth_m=np.zeros(height.size),
        land_temperature_c=np.zeros(height.size),
        air_temperature_c=np.zeros(height.size),
        ocean_temperature_c=np.zeros(height.size),
        precipitation_mm=np.zeros(height.size),
        biome_ids=np.zeros(height.size, dtype=int),
    )
    monkeypatch.setattr(stats.climate, "compute_climate_cached", lambda world: fields)
    world = _all_ocean_world()
    world.sea_level_m = sea_level
    assert stats.compute_stats(world)["land_near_max_elevation_fraction"] == pytest.approx(expected)
