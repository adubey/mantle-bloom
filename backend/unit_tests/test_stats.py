import numpy as np
import pytest

from app import geometry, lithosphere, stats
from app.hydroclimate import climate
from app.elevation_lines import line_spacing_rad
from app.world import World, generate_world

from .quad_fixtures import globe_plate, quad_plate

_GLOBE_DENSITY = 0.5


def _globe_world(crust_type: str, **fields_at) -> World:
    return World(seed=0, plates=[globe_plate(crust_type, _GLOBE_DENSITY, **fields_at)], node_density=_GLOBE_DENSITY)


def _longitude(points: np.ndarray) -> np.ndarray:
    return geometry.xyz_to_latlon(points)[1]


def _all_ocean_world() -> World:
    return _globe_world("oceanic", elevation=lambda pts: np.full(len(pts), -3800.0))


@pytest.mark.parametrize("climate_density", [0.5, 2.0])
def test_spherical_area_weights_match_analytic_latitude_band(climate_density):
    height, width = climate.grid_dimensions(climate_density)
    lat_deg, _, _ = climate._build_grid(height, width)
    weights = stats._spherical_area_weights(lat_deg, width)

    # Choose boundaries that coincide exactly with this grid at both tested resolutions.
    band = (lat_deg <= 30.0) & (lat_deg >= -30.0)
    actual = weights[band].sum() / weights.sum()
    expected = (np.sin(np.radians(30.0)) - np.sin(np.radians(-30.0))) / 2.0
    assert actual == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("climate_density", [0.5, 2.0])
def test_equal_cell_counts_at_equator_and_pole_use_physical_area(climate_density):
    height, width = climate.grid_dimensions(climate_density)
    lat_deg, _, _ = climate._build_grid(height, width)
    weights = stats._spherical_area_weights(lat_deg, width)
    cell_count = width // 4

    equatorial = np.zeros((height, width), dtype=bool)
    polar = np.zeros_like(equatorial)
    equatorial[np.argmin(np.abs(lat_deg)), :cell_count] = True
    polar[0, :cell_count] = True

    assert equatorial.sum() == polar.sum()
    assert weights[equatorial].sum() > weights[polar].sum()


def test_weighted_summary_preserves_extrema_and_weights_mean_and_std():
    values = np.array([0.0, 10.0])
    weights = np.array([3.0, 1.0])
    summary = stats.Stat4.of(values, weights)
    assert summary.min == 0.0
    assert summary.max == 10.0
    assert summary.mean == pytest.approx(2.5)
    assert summary.std == pytest.approx(np.sqrt(18.75))


def test_compute_stats_land_and_ocean_fractions_sum_to_one():
    world = _all_ocean_world()
    result = stats.compute_stats(world)
    assert result["land_fraction"] == 0.0
    assert result["ocean_fraction"] == 1.0


def test_compute_stats_all_land_is_exactly_one():
    world = _all_ocean_world()
    for plate in world.plates:
        plate.set_fields_on_plate(elevation=np.full(plate.node_count(), 500.0))
    result = stats.compute_stats(world)
    assert result["land_fraction"] == 1.0
    assert result["ocean_fraction"] == 0.0


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
    # One continental plate, every cell holding a known Hc -- total volume should be exactly
    # sum(Hc * cell area), converted to km^3 (see Plate.accounting_areas_m2).
    n = 20
    hc = np.linspace(20_000.0, 40_000.0, n)
    plate = quad_plate(0, "continental", columns=range(n), elevation=500.0, crustal_thickness_m=hc)
    world = World(seed=0, plates=[plate])
    result = stats.compute_stats(world)

    area_m2 = plate.node_areas_m2()
    assert result["total_continental_crust_volume_km3"] == pytest.approx(float(hc @ area_m2) / 1.0e9)
    # Every cell here sits above the default sea_level_m=0.0, so the whole plate is land.
    assert result["total_land_area_km2"] == pytest.approx(float(area_m2.sum()) / 1.0e6)


def test_overlap_once_totals_count_stacked_ground_once():
    # Issue #289: two continental plates whose 10-cell rows sit exactly on top of each other.
    # The plain totals count both plates in full (as the ledger does); the dedup totals count
    # the shared ground once.
    n = 10
    hc = 35_000.0

    def plate(plate_id):
        return quad_plate(plate_id, "continental", columns=range(n), elevation=500.0, crustal_thickness_m=hc)

    world = World(seed=0, plates=[plate(0), plate(1)], next_plate_id=2)
    result = stats.compute_stats(world)

    area_m2 = float(world.plates[0].node_areas_m2().sum())
    assert result["total_continental_crust_volume_km3"] == pytest.approx(2 * hc * area_m2 / 1.0e9)
    assert result["total_continental_crust_volume_dedup_km3"] == pytest.approx(hc * area_m2 / 1.0e9)
    assert result["total_land_area_km2"] == pytest.approx(2 * area_m2 / 1.0e6)
    assert result["total_land_area_dedup_km2"] == pytest.approx(area_m2 / 1.0e6)
    assert result["plate_overlap_area_fraction"] == pytest.approx(2 * area_m2 / stats.SPHERE_AREA_M2)
    assert result["represented_area_fraction"] == pytest.approx(2 * area_m2 / stats.SPHERE_AREA_M2)


def test_continental_dedup_volume_ignores_an_oceanic_plate_on_top():
    # PR #295 review: an oceanic plate overlapping a continental one duplicates no continental
    # crust, so the continental plate keeps its full volume; land area still counts once.
    n = 10
    hc = 35_000.0

    def plate(plate_id, crust_type):
        return quad_plate(plate_id, crust_type, columns=range(n), elevation=500.0, crustal_thickness_m=hc)

    world = World(seed=0, plates=[plate(0, "continental"), plate(1, "oceanic")], next_plate_id=2)
    result = stats.compute_stats(world)

    area_m2 = float(world.plates[0].node_areas_m2().sum())
    assert result["total_continental_crust_volume_dedup_km3"] == pytest.approx(hc * area_m2 / 1.0e9)
    assert result["total_land_area_dedup_km2"] == pytest.approx(area_m2 / 1.0e6)


def test_overlap_once_totals_match_plain_totals_without_overlap():
    world = World(seed=0, plates=[quad_plate(0, "continental", columns=range(20), elevation=500.0, crustal_thickness_m=35_000.0)])
    result = stats.compute_stats(world)
    assert result["total_continental_crust_volume_dedup_km3"] == result["total_continental_crust_volume_km3"]
    assert result["plate_overlap_area_fraction"] == 0.0


def test_total_land_area_uses_each_quad_cells_own_area():
    # Issue #257: quad cells are not equal-area, so land area is the sum of the land cells'
    # own areas, not land-node count times the nominal area.
    world = generate_world(seed=5, num_plates=5)
    land_area = 0.0
    for plate in world.plates:
        nodes = plate.surface_nodes("elevation")
        land_area += float(nodes.area_m2[nodes.fields["elevation"] > world.sea_level_m].sum())
    land_area_m2, _, land_nodes = stats._total_land_area_and_continental_volume(world)
    assert land_nodes > 0
    assert land_area_m2 == pytest.approx(land_area)
    nominal = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    assert land_area_m2 != pytest.approx(land_nodes * nominal, rel=1e-3)


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
    from app.hydroclimate import biomes

    assert all(name in biomes.PELAGIC_NAMES for name in result["biome_ocean_fraction"])


def _land_and_ocean_world() -> World:
    # Half land (the hemisphere around longitude 0), half ocean.
    return _globe_world("continental", elevation=lambda pts: np.where(np.abs(_longitude(pts)) < np.pi / 2, 500.0, -3800.0))


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
    world = _globe_world(
        "continental",
        elevation=lambda pts: np.full(len(pts), 500.0),  # all land, well above the default sea_level_m=0.0
        lake_depth=lambda pts: np.where(np.abs(_longitude(pts)) < np.pi / 4, 20.0, 0.0),  # a lake over a quarter of it
    )
    plate = world.plates[0]

    result = stats.compute_stats(world)
    areas = plate.node_areas_m2()
    lake_fraction = float(areas[plate.collect("lake_depth") > 0.0].sum() / areas.sum())
    assert result["ocean_fraction"] == pytest.approx(lake_fraction, abs=0.05)
    assert result["land_fraction"] == pytest.approx(1.0 - lake_fraction, abs=0.05)
    assert result["land_fraction"] + result["ocean_fraction"] == pytest.approx(1.0)


def test_compute_stats_biome_land_fraction_reads_the_stored_climate_cache_biome_ids():
    # stats.py no longer runs its own classify_biomes -- it reads ClimateFields.biome_ids,
    # the same stored field compute_climate now computes once and every other biome-consuming
    # caller shares (see climate.py). Recomputing the expected fractions directly from
    # world.climate_cache.biome_ids (populated as a side effect of compute_stats calling
    # compute_climate_cached) should match stats.py's own result exactly.
    from app.hydroclimate import biomes

    world = _land_and_ocean_world()
    result = stats.compute_stats(world)

    biome_ids = world.climate_cache.biome_ids
    is_land = ~world.climate_cache.is_ocean
    weights = stats._spherical_area_weights(world.climate_cache.lat_deg, biome_ids.shape[1])
    land_biome_ids = biome_ids[is_land]
    land_weights = weights[is_land]
    land_weight = float(land_weights.sum())
    expected = {
        name: float(land_weights[land_biome_ids == i].sum()) / land_weight
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i not in biomes.OCEAN_IDS and land_weight > 0.0
    }
    assert result["biome_land_fraction"] == expected


def test_compute_stats_land_fraction_node_matches_raw_elevation_count():
    # land_fraction_node (see stats.py's own docstring, GitHub issue #121) is a plain
    # elevation > sea_level_m node count/total -- unlike land_fraction, not connectivity-aware
    # and not resampled from any hydrology cache. Compare against an independent recount
    # straight off world.plates.
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
        lat_deg=np.array([0.0]),
        elevation_m=(height + sea_level)[None, :],
        is_ocean=(height > 0)[None, :],
        lake_depth_m=np.zeros((1, height.size)),
        land_temperature_c=np.zeros((1, height.size)),
        air_temperature_c=np.zeros((1, height.size)),
        ocean_temperature_c=np.zeros((1, height.size)),
        precipitation_mm=np.zeros((1, height.size)),
        biome_ids=np.zeros((1, height.size), dtype=int),
    )
    monkeypatch.setattr(stats.climate, "compute_climate_cached", lambda world: fields)
    world = _all_ocean_world()
    world.sea_level_m = sea_level
    assert stats.compute_stats(world)["land_near_max_elevation_fraction"] == pytest.approx(expected)
