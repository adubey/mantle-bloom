import numpy as np

from app import biomes, hydrology


def test_classify_biomes_ocean_wins_regardless_of_temperature_and_precipitation():
    # Even a "hot and wet" cell classifies as Ocean if is_ocean says so -- temperature/
    # precipitation are land-surface concepts, is_ocean already settles the question.
    temp = np.array([30.0])
    precip = np.array([3000.0])
    is_ocean = np.array([True])
    assert biomes.classify_biomes(temp, precip, is_ocean)[0] == biomes.OCEAN


def test_classify_biomes_ice_below_the_glacier_accumulation_threshold():
    # Reuses hydrology.GLACIER_ACCUMULATION_TEMP_C directly -- a biome map's Ice region should
    # line up with where the simulation would actually grow a glacier.
    temp = np.array([hydrology.GLACIER_ACCUMULATION_TEMP_C - 1.0])
    precip = np.array([500.0])
    is_ocean = np.array([False])
    assert biomes.classify_biomes(temp, precip, is_ocean)[0] == biomes.ICE


def test_classify_biomes_covers_every_named_biome_across_the_temperature_precipitation_grid():
    # A coarse sweep across plausible temperature/precipitation combinations should produce
    # every non-Ocean biome at least once -- confirms no band is unreachable/dead code.
    temps = np.linspace(-30.0, 35.0, 40)
    precips = np.linspace(0.0, 4000.0, 40)
    temp_grid, precip_grid = np.meshgrid(temps, precips)
    is_ocean = np.zeros_like(temp_grid, dtype=bool)

    result = biomes.classify_biomes(temp_grid, precip_grid, is_ocean)
    seen = set(np.unique(result).tolist())
    expected = set(range(len(biomes.BIOME_NAMES))) - {biomes.OCEAN}
    assert expected <= seen


def test_classify_biomes_matches_shape_of_inputs():
    temp = np.zeros((5, 7))
    precip = np.full((5, 7), 800.0)
    is_ocean = np.zeros((5, 7), dtype=bool)
    result = biomes.classify_biomes(temp, precip, is_ocean)
    assert result.shape == (5, 7)


def test_classify_biomes_boundary_values_go_to_the_wetter_warmer_side():
    # np.select's first-matching-condition semantics mean a boundary value (e.g. exactly
    # COLD_TEMP_C) belongs to the band that starts there (>=), not the band below it -- pin
    # this down explicitly so a future refactor can't silently flip which side a boundary
    # falls on.
    is_ocean = np.array([False])
    assert biomes.classify_biomes(np.array([biomes.COLD_TEMP_C]), np.array([1000.0]), is_ocean)[0] != biomes.TUNDRA
    assert biomes.classify_biomes(np.array([biomes.ICE_TEMP_C]), np.array([1000.0]), is_ocean)[0] != biomes.ICE
    assert biomes.classify_biomes(np.array([20.0]), np.array([biomes.SEMI_ARID_MM]), is_ocean)[0] != biomes.SUBTROPICAL_DESERT


def test_biome_colors_index_aligned_with_biome_names():
    assert len(biomes.BIOME_COLORS) == len(biomes.BIOME_NAMES)
    assert biomes.BIOME_COLORS.shape == (len(biomes.BIOME_NAMES), 3)
    assert biomes.BIOME_COLORS.dtype == np.uint8
