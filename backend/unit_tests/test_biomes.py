import numpy as np

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


def test_biome_relative_shade_factor_only_ever_returns_a_known_factor():
    biome_ids = np.array([biomes.TUNDRA, biomes.TUNDRA, biomes.OCEAN])
    elevation = np.array([100.0, 2000.0, -500.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert set(result.tolist()) <= set(biomes.BIOME_SHADE_FACTORS.tolist())


def test_biome_relative_shade_factor_splits_a_large_uniform_biome_into_multiple_tiers():
    # A single biome across many cells with a real elevation spread should use more than one
    # shade tier -- otherwise a large biome area would still render as one flat color.
    n = 300
    biome_ids = np.full(n, biomes.TEMPERATE_GRASSLAND)
    elevation = np.linspace(0.0, 1000.0, n)
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert len(np.unique(result)) >= biomes.BIOME_SHADE_TIERS


def test_biome_relative_shade_factor_shaded_colors_stay_close_to_the_flat_biome_color():
    # Every shaded variant should still be near BIOME_COLORS' own flat entry -- close enough
    # that frontend/src/legendData.ts's own click-to-highlight tolerance can match it without
    # a second legend swatch per biome (see BIOME_SHADE_FACTORS' own docstring).
    n = 50
    biome_ids = np.full(n, biomes.SAVANNA)
    elevation = np.linspace(0.0, 3000.0, n)
    factor = biomes.biome_relative_shade_factor(biome_ids, elevation)
    base = biomes.BIOME_COLORS[biomes.SAVANNA].astype(float)
    shaded = np.clip(base[None, :] * factor[:, None], 0, 255)
    dist = np.linalg.norm(shaded - base[None, :], axis=1)
    assert np.all(dist <= 35)


def test_biome_relative_shade_factor_ranks_higher_elevation_cells_brighter_within_a_biome():
    biome_ids = np.array([biomes.BOREAL_FOREST, biomes.BOREAL_FOREST])
    elevation = np.array([0.0, 5000.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert result[1] > result[0]


def test_biome_relative_shade_factor_ranks_are_relative_to_each_biome_separately():
    # A biome confined to a narrow absolute elevation band should still get real tier spread,
    # since ranking is relative to that biome's own cells, not a shared absolute scale.
    n = 60
    biome_ids = np.concatenate([np.full(n, biomes.WETLAND), np.full(n, biomes.TUNDRA)])
    elevation = np.concatenate([np.linspace(0.0, 10.0, n), np.linspace(0.0, 3000.0, n)])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    wetland_tiers = len(np.unique(result[:n]))
    tundra_tiers = len(np.unique(result[n:]))
    assert wetland_tiers == biomes.BIOME_SHADE_TIERS
    assert tundra_tiers == biomes.BIOME_SHADE_TIERS


def test_classify_biomes_intertidal_is_shallow_ocean_only():
    is_ocean = np.array([True, True])
    elevation = np.array([-5.0, -5000.0])  # shallow vs. deep
    slope = np.array([0.0, 0.0])
    temp = np.array([20.0, 20.0])
    precip = np.array([1000.0, 1000.0])
    result = biomes.classify_biomes(temp, precip, elevation, slope, is_ocean, sea_level_m=0.0)
    assert result[0] == biomes.INTERTIDAL
    assert result[1] == biomes.OCEAN
