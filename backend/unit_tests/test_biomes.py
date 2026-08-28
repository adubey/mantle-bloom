import numpy as np
import pytest

from app import biomes, hydrology

# A deliberately "not flat, not low-lying" elevation/slope pair -- used by every test below
# that only cares about the ordinary temperature/precipitation bands, so a stray match against
# classify_wetland's own flat-and-low condition can't silently change which biome a cell lands
# in. See test_classify_wetland_* below for tests that actually exercise that axis.
_NOT_WETLAND_ELEVATION = np.array([1000.0])
_NOT_WETLAND_SLOPE = np.array([0.01])


def test_classify_biomes_ocean_wins_regardless_of_temperature_and_precipitation():
    # Even a "hot and wet" cell classifies as Ocean if is_ocean says so -- temperature/
    # precipitation are land-surface concepts, is_ocean already settles the question.
    temp = np.array([30.0])
    precip = np.array([3000.0])
    is_ocean = np.array([True])
    elevation = np.array([-5000.0])  # deep, so this doesn't also read as Intertidal
    slope = np.array([0.0])
    assert biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)[0] == biomes.OCEAN


def test_classify_biomes_ice_below_the_glacier_accumulation_threshold():
    # Reuses hydrology.GLACIER_ACCUMULATION_TEMP_C directly -- a biome map's Ice region should
    # line up with where the simulation would actually grow a glacier.
    temp = np.array([hydrology.GLACIER_ACCUMULATION_TEMP_C - 1.0])
    precip = np.array([500.0])
    is_ocean = np.array([False])
    assert biomes.classify_biomes(temp, precip, _NOT_WETLAND_ELEVATION, _NOT_WETLAND_SLOPE, is_ocean)[0] == biomes.ICE


def test_classify_biomes_covers_every_ordinary_temperature_precipitation_biome():
    # A coarse sweep across plausible temperature/precipitation combinations (elevation/slope
    # held fixed at a "not wetland" value) should produce every ordinary temp/precip-only
    # biome at least once -- confirms no band is unreachable/dead code. Wetland/Carboniferous
    # Forest/Intertidal Zone need their own axes to vary (see the dedicated tests below), so
    # they're deliberately excluded from "expected" here rather than left to fail.
    temps = np.linspace(-30.0, 35.0, 40)
    precips = np.linspace(0.0, 4000.0, 40)
    temp_grid, precip_grid = np.meshgrid(temps, precips)
    is_ocean = np.zeros_like(temp_grid, dtype=bool)
    elevation = np.full_like(temp_grid, 1000.0)
    slope = np.full_like(temp_grid, 0.01)

    result = biomes.classify_biomes(temp_grid, precip_grid, elevation, slope, is_ocean)
    seen = set(np.unique(result).tolist())
    expected = set(range(len(biomes.BIOME_NAMES))) - {biomes.OCEAN, biomes.WETLAND, biomes.CARBONIFEROUS_FOREST, biomes.INTERTIDAL}
    assert expected <= seen


def test_classify_biomes_matches_shape_of_inputs():
    temp = np.zeros((5, 7))
    precip = np.full((5, 7), 800.0)
    elevation = np.full((5, 7), 1000.0)
    slope = np.full((5, 7), 0.01)
    is_ocean = np.zeros((5, 7), dtype=bool)
    result = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)
    assert result.shape == (5, 7)


def test_classify_biomes_boundary_values_go_to_the_wetter_warmer_side():
    # np.select's first-matching-condition semantics mean a boundary value (e.g. exactly
    # COLD_TEMP_C) belongs to the band that starts there (>=), not the band below it -- pin
    # this down explicitly so a future refactor can't silently flip which side a boundary
    # falls on.
    is_ocean = np.array([False])
    args = (_NOT_WETLAND_ELEVATION, _NOT_WETLAND_SLOPE, is_ocean)
    assert biomes.classify_biomes(np.array([biomes.COLD_TEMP_C]), np.array([1000.0]), *args)[0] != biomes.TUNDRA
    assert biomes.classify_biomes(np.array([biomes.ICE_TEMP_C]), np.array([1000.0]), *args)[0] != biomes.ICE
    assert biomes.classify_biomes(np.array([20.0]), np.array([biomes.SEMI_ARID_MM]), *args)[0] != biomes.SUBTROPICAL_DESERT


def test_biome_colors_index_aligned_with_biome_names():
    assert len(biomes.BIOME_COLORS) == len(biomes.BIOME_NAMES)
    assert biomes.BIOME_COLORS.shape == (len(biomes.BIOME_NAMES), 3)
    assert biomes.BIOME_COLORS.dtype == np.uint8


def test_classify_wetland_requires_flat_low_land():
    # Warm and wet, but steep -- real relief drains rather than waterlogs, so this shouldn't
    # read as wetland even though temperature/precipitation alone would qualify.
    temp = np.array([25.0])
    precip = np.array([2500.0])
    is_ocean = np.array([False])
    is_wetland, is_carboniferous = biomes.classify_wetland(temp, precip, np.array([10.0]), np.array([0.05]), is_ocean)
    assert not is_wetland[0] and not is_carboniferous[0]

    # Same temperature/precipitation, but flat and low -- now it should qualify, and (being
    # warm + very wet) as the Carboniferous Forest subtype specifically.
    is_wetland, is_carboniferous = biomes.classify_wetland(temp, precip, np.array([10.0]), np.array([0.0001]), is_ocean)
    assert is_carboniferous[0] and not is_wetland[0]


def test_classify_wetland_cooler_flat_low_land_is_plain_wetland_not_carboniferous():
    temp = np.array([10.0])  # below CARBONIFEROUS_MIN_TEMP_C
    precip = np.array([1500.0])
    is_ocean = np.array([False])
    is_wetland, is_carboniferous = biomes.classify_wetland(temp, precip, np.array([5.0]), np.array([0.0001]), is_ocean)
    assert is_wetland[0] and not is_carboniferous[0]


def test_classify_wetland_excludes_ocean():
    is_wetland, is_carboniferous = biomes.classify_wetland(
        np.array([25.0]), np.array([2500.0]), np.array([-5.0]), np.array([0.0001]), np.array([True])
    )
    assert not is_wetland[0] and not is_carboniferous[0]


def test_classify_biomes_wetland_and_carboniferous_forest_reachable():
    is_ocean = np.array([False, False])
    elevation = np.array([5.0, 5.0])
    slope = np.array([0.0001, 0.0001])
    temp = np.array([10.0, 25.0])
    precip = np.array([1500.0, 2500.0])
    result = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)
    assert result[0] == biomes.WETLAND
    assert result[1] == biomes.CARBONIFEROUS_FOREST


def test_biome_relative_shade_factor_matches_shape_of_biome_ids():
    biome_ids = np.array([biomes.TUNDRA, biomes.TUNDRA, biomes.OCEAN])
    elevation = np.array([100.0, 2000.0, -500.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert result.shape == (3,)


def test_biome_relative_shade_factor_stays_within_the_amplitude_bounds():
    biome_ids = np.array([biomes.TUNDRA, biomes.TUNDRA, biomes.TUNDRA, biomes.OCEAN])
    elevation = np.array([100.0, 2000.0, 900.0, -500.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert np.all(result >= biomes.BIOME_SHADE_MIN - 1e-9)
    assert np.all(result <= biomes.BIOME_SHADE_MAX + 1e-9)


def test_biome_relative_shade_factor_spans_the_full_range_end_to_end():
    # A biome's lowest cell hits BIOME_SHADE_MIN and its highest hits BIOME_SHADE_MAX.
    n = 200
    biome_ids = np.full(n, biomes.TEMPERATE_GRASSLAND)
    elevation = np.linspace(0.0, 1000.0, n)
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert result.min() == pytest.approx(biomes.BIOME_SHADE_MIN)
    assert result.max() == pytest.approx(biomes.BIOME_SHADE_MAX)


def test_biome_relative_shade_factor_is_continuous_not_a_few_discrete_tiers():
    # The whole point of the rewrite: a large biome varies smoothly, so within-biome relief
    # doesn't read as visible bands. Many distinct factor values, not ~3.
    n = 300
    biome_ids = np.full(n, biomes.TEMPERATE_GRASSLAND)
    elevation = np.linspace(0.0, 1000.0, n)
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert len(np.unique(result)) == n


def test_biome_relative_shade_factor_increases_monotonically_with_elevation_rank():
    n = 50
    biome_ids = np.full(n, biomes.SAVANNA)
    elevation = np.linspace(0.0, 3000.0, n)
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert np.all(np.diff(result) > 0)


def test_biome_relative_shade_factor_ranks_higher_elevation_cells_brighter_within_a_biome():
    biome_ids = np.array([biomes.BOREAL_FOREST, biomes.BOREAL_FOREST])
    elevation = np.array([0.0, 5000.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert result[1] > result[0]


def test_biome_relative_shade_factor_ranks_are_relative_to_each_biome_separately():
    # A biome confined to a narrow absolute elevation band still spans the full brightness
    # range, since ranking is relative to that biome's own cells, not a shared absolute scale.
    n = 60
    biome_ids = np.concatenate([np.full(n, biomes.WETLAND), np.full(n, biomes.TUNDRA)])
    elevation = np.concatenate([np.linspace(0.0, 10.0, n), np.linspace(0.0, 3000.0, n)])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    for segment in (result[:n], result[n:]):
        assert segment.min() == pytest.approx(biomes.BIOME_SHADE_MIN)
        assert segment.max() == pytest.approx(biomes.BIOME_SHADE_MAX)


def test_classify_biomes_intertidal_is_shallow_ocean_only():
    is_ocean = np.array([True, True])
    elevation = np.array([-5.0, -5000.0])  # shallow vs. deep
    slope = np.array([0.0, 0.0])
    temp = np.array([20.0, 20.0])
    precip = np.array([1000.0, 1000.0])
    result = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean, sea_level_m=0.0)
    assert result[0] == biomes.INTERTIDAL
    assert result[1] == biomes.OCEAN


def test_grid_slope_is_zero_on_a_flat_grid():
    lat_deg = np.array([30.0, 0.0, -30.0])
    elevation_m = np.full((3, 4), 500.0)
    assert np.all(biomes.grid_slope(elevation_m, lat_deg) == 0.0)


def test_grid_slope_matches_a_known_north_south_step():
    # A single elevation step between two rows, everything else flat: the rise/run slope at
    # the step should equal that rise divided by the real great-circle spacing between rows
    # (dlat_km, in meters) -- np.roll compares each row to the row above it (axis=0), so row 1
    # (the step) reads a nonzero d_ns against row 0, and (via wraparound) row 0 reads that same
    # step against the last row.
    lat_deg = np.array([30.0, 0.0, -30.0])
    elevation_m = np.zeros((3, 4))
    elevation_m[1, :] = 1000.0
    slope = biomes.grid_slope(elevation_m, lat_deg)
    dlat_km = (np.pi / 3) * biomes.PLANET_RADIUS_KM
    expected = 1000.0 / (dlat_km * 1000.0)
    assert np.allclose(slope[1, :], expected)


# --- smooth_biome_field: the stateless boundary-cleanup pass -----------------------------

# Temperate, sub-humid, upland: WOODLAND_SHRUBLAND, comfortably away from every band edge.
_WOODLAND_TEMP_C = 10.0
_WOODLAND_PRECIP_MM = 750.0


def _uniform(shape, temp, precip, elevation=1000.0, slope=0.01, is_ocean=False):
    return (
        np.full(shape, float(temp)),
        np.full(shape, float(precip)),
        np.full(shape, float(elevation)),
        np.full(shape, float(slope)),
        np.full(shape, bool(is_ocean)),
    )


def test_smooth_biome_field_rejects_a_flat_array():
    temp, precip, elevation, slope, is_ocean = (a.reshape(-1) for a in _uniform((4, 4), 10.0, 750.0))
    with pytest.raises(ValueError):
        biomes.smooth_biome_field(temp, precip, elevation, slope, is_ocean)


def test_smooth_biome_field_preserves_shape():
    args = _uniform((5, 8), _WOODLAND_TEMP_C, _WOODLAND_PRECIP_MM)
    assert biomes.smooth_biome_field(*args).shape == (5, 8)


def test_smooth_biome_field_is_a_noop_on_a_spatially_uniform_climate():
    # Nothing sits near a band edge and every cell agrees with its neighbours, so the cleanup
    # pass has nothing to do -- the result is exactly classify_biomes.
    args = _uniform((6, 9), _WOODLAND_TEMP_C, _WOODLAND_PRECIP_MM)
    raw = biomes.classify_biomes(*args)
    assert np.array_equal(biomes.smooth_biome_field(*args), raw)


def test_smooth_biome_field_outvotes_a_lone_cusp_speckle():
    # One cell nudged just across the humid band edge (within BIOME_CUSP_MARGIN_PRECIP_MM of
    # it) in an otherwise uniform WOODLAND_SHRUBLAND field -> reverts to its neighbours' biome.
    temp, precip, elevation, slope, is_ocean = _uniform((7, 7), _WOODLAND_TEMP_C, _WOODLAND_PRECIP_MM)
    precip[3, 3] = biomes.SUB_HUMID_MM + 40.0  # -> TEMPERATE_SEASONAL_FOREST, 40 mm over the edge
    raw = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)
    assert raw[3, 3] == biomes.TEMPERATE_SEASONAL_FOREST

    smoothed = biomes.smooth_biome_field(temp, precip, elevation, slope, is_ocean)
    assert smoothed[3, 3] == biomes.WOODLAND_SHRUBLAND


def test_smooth_biome_field_leaves_a_clean_two_biome_interface_in_place():
    # A straight interface: each side's cells see only ~5/8 same-biome neighbours, below the
    # BIOME_VOTE_MIN_NEIGHBOUR_FRACTION supermajority, so neither side is eroded.
    temp, precip, elevation, slope, is_ocean = _uniform((8, 8), _WOODLAND_TEMP_C, _WOODLAND_PRECIP_MM)
    precip[:, 4:] = biomes.SUB_HUMID_MM + 40.0  # right half just over the humid edge
    raw = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)
    assert np.array_equal(biomes.smooth_biome_field(temp, precip, elevation, slope, is_ocean), raw)


def test_smooth_biome_field_keeps_a_small_solid_region_far_from_any_cutoff():
    # A 3x3 block of a different biome, its climate nowhere near a band edge: no cell is a
    # cusp cell and no cell is outnumbered enough to be "speckle", so the block is untouched.
    temp, precip, elevation, slope, is_ocean = _uniform((7, 7), _WOODLAND_TEMP_C, _WOODLAND_PRECIP_MM)
    precip[2:5, 2:5] = 1600.0  # deep in the TEMPERATE_SEASONAL_FOREST band (edges at 1000 / 2000)
    smoothed = biomes.smooth_biome_field(temp, precip, elevation, slope, is_ocean)
    assert np.all(smoothed[2:5, 2:5] == biomes.TEMPERATE_SEASONAL_FOREST)
    assert np.all(smoothed[0, :] == biomes.WOODLAND_SHRUBLAND)


def test_smooth_biome_field_revotes_a_single_cell_spike_inside_a_wetland():
    # A lone elevation spike disqualifies one cell from Wetland (elevation > WETLAND_MAX_
    # ELEVATION_M) inside an otherwise solid flat/low/wet region; the neighbour vote, which
    # always reconsiders a fully-outnumbered land cell, puts it back.
    shape = (7, 7)
    temp, precip, elevation, slope, is_ocean = _uniform(
        shape, 12.0, 1200.0, elevation=20.0, slope=0.0005, is_ocean=False
    )
    elevation[3, 3] = 200.0
    raw = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean)
    assert raw[3, 3] != biomes.WETLAND and np.all(raw[0, :] == biomes.WETLAND)

    smoothed = biomes.smooth_biome_field(temp, precip, elevation, slope, is_ocean)
    assert smoothed[3, 3] == biomes.WETLAND


def test_smooth_for_bands_attenuates_a_spike_and_spreads_it():
    field = np.zeros((40, 40))
    field[20, 20] = 100.0
    out = biomes._smooth_for_bands(field, cell_km=30.0)  # 30 km cells -> a multi-cell window
    assert out[20, 20] < 20.0  # spike knocked down
    assert out[20, 21] > 0.0 and out[19, 20] > 0.0  # mass spread to neighbours
    assert out.max() <= field.max() + 1e-9
