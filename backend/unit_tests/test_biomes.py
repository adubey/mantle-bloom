import numpy as np
import pytest

from app import biomes, hydrology

_UPLAND_ELEVATION = np.array([1000.0])
_UPLAND_SLOPE = np.array([0.01])


def _koppen(t, p, lat=0.0, cont=0.0, tilt=biomes.DEFAULT_AXIAL_TILT_DEG):
    arr = lambda v: np.array([float(v)])
    i = int(biomes.classify_koppen(arr(t), arr(p), arr(lat), arr(cont), tilt)[0])
    return biomes.KOPPEN_CODES[i]


def _classify(temp, precip, is_ocean, lat=0.0, elevation=1000.0, slope=0.01, **kw):
    shape = np.shape(temp)
    return biomes.classify_biomes(
        np.asarray(temp, dtype=float),
        np.asarray(precip, dtype=float),
        np.full(shape, float(elevation)) if np.ndim(elevation) == 0 else np.asarray(elevation, dtype=float),
        np.full(shape, float(slope)) if np.ndim(slope) == 0 else np.asarray(slope, dtype=float),
        np.asarray(is_ocean),
        lat_deg=np.full(shape, float(lat)) if np.ndim(lat) == 0 else np.asarray(lat, dtype=float),
        **kw,
    )


# --- Köppen land classification --------------------------------------------------------


def test_classify_koppen_main_classes_reachable():
    # One representative (temperature, precipitation, latitude, continentality) per Köppen
    # main class -- confirms A/B/C/D/E are all reachable.
    assert _koppen(26, 2200, lat=-3, cont=0.3)[0] == "A"
    assert _koppen(25, 20, lat=22, cont=0.6)[0] == "B"
    assert _koppen(11, 900, lat=48, cont=0.1)[0] == "C"
    assert _koppen(-4, 400, lat=60, cont=0.9)[0] == "D"
    assert _koppen(-25, 180, lat=75, cont=0.6)[0] == "E"


def test_classify_koppen_desert_vs_steppe_and_hot_vs_cold():
    assert _koppen(25, 20, lat=22, cont=0.6) == "BWh"
    assert _koppen(25, 320, lat=18, cont=0.4) in ("BSh", "BWh")
    assert _koppen(7, 90, lat=44, cont=1.0) == "BWk"


def test_classify_koppen_tropical_rainforest_needs_no_dry_month():
    assert _koppen(26, 2400, lat=-2, cont=0.2) == "Af"
    # Same warmth, far less rain, strong dry season -> savanna, not rainforest.
    assert _koppen(27, 900, lat=14, cont=0.3) == "Aw"


def test_classify_koppen_ice_cap_uses_the_glacier_threshold_region():
    # A cell whose warmest month stays below freezing is an ice cap.
    assert _koppen(-30, 150, lat=80, cont=0.5) == "EF"
    # Warmest month above 0 but below 10 -> tundra.
    assert _koppen(-8, 250, lat=72, cont=0.3) == "ET"


def test_classify_koppen_tilt_zero_world_has_no_seasonal_subtypes():
    # With no axial tilt there is no seasonal temperature swing, so `s` (dry-summer) and `w`
    # (dry-winter) subtypes and the `d` extreme-winter third letter can never occur.
    rng = np.random.default_rng(0)
    lat = rng.uniform(-85, 85, 4000)
    t = 28 - 0.55 * np.abs(lat) + rng.normal(0, 4, lat.shape)
    p = rng.uniform(20, 3000, lat.shape)
    cont = rng.uniform(0, 1, lat.shape)
    ids = biomes.classify_koppen(t, p, lat, cont, axial_tilt_deg=0.0)
    codes = {biomes.KOPPEN_CODES[i] for i in np.unique(ids)}
    # No temperate/continental dry-summer or dry-winter subtypes, and no `d` (extreme cold
    # winter) third letter -- all of those need a seasonal temperature cycle.
    assert not any(c[0] in ("C", "D") and c[1] in ("s", "w") for c in codes)
    assert not any(c[0] in ("C", "D") and c.endswith("d") for c in codes)


def test_classify_koppen_mediterranean_appears_for_a_maritime_subtropical_world():
    # Warm dry-summer subtropical, near a coast (low continentality), modest precipitation.
    assert _koppen(16, 500, lat=38, cont=0.15) in ("Csa", "Csb")


# --- Pelagic ocean classification -----------------------------------------------------


def test_classify_pelagic_realm_from_sst():
    lat = np.array([15.0, 40.0, 55.0, 75.0])
    sst = np.array([27.0, 15.0, 6.0, -1.5])
    dist = np.full(4, 1.0)  # far from land
    ids = biomes.classify_pelagic(sst, lat, dist)
    names = [biomes.BIOME_NAMES[i] for i in ids]
    assert names == [
        "Tropical Open Ocean",
        "Temperate Open Ocean",
        "Cold-Temperate Open Ocean",
        "Polar Sea Ice",
    ]


def test_classify_pelagic_shelf_vs_open_ocean():
    sst = np.array([26.0, 26.0])
    lat = np.array([20.0, 20.0])
    dist = np.array([biomes.SHELF_RANGE_RAD * 0.5, biomes.SHELF_RANGE_RAD * 3.0])
    ids = biomes.classify_pelagic(sst, lat, dist)
    assert biomes.BIOME_NAMES[ids[0]] == "Tropical Coastal Waters"
    assert biomes.BIOME_NAMES[ids[1]] in ("Tropical Open Ocean", "Subtropical Gyre", "Equatorial Divergence")


def test_classify_pelagic_equatorial_divergence_band():
    ids = biomes.classify_pelagic(np.array([27.0]), np.array([2.0]), np.array([1.0]))
    assert biomes.BIOME_NAMES[int(ids[0])] == "Equatorial Divergence"


# --- classify_biomes: land/ocean split ----------------------------------------------


def test_classify_biomes_ocean_wins_regardless_of_temperature_and_precipitation():
    result = _classify([30.0], [3000.0], [True], lat=0.0, elevation=-5000.0)
    assert int(result[0]) in biomes.OCEAN_IDS


def test_classify_biomes_land_cell_is_a_koppen_class():
    result = _classify([10.0], [800.0], [False], lat=48.0)
    assert int(result[0]) not in biomes.OCEAN_IDS
    assert 0 <= int(result[0]) < len(biomes.KOPPEN_CODES)


def test_classify_biomes_matches_shape_of_inputs():
    temp = np.zeros((5, 7))
    precip = np.full((5, 7), 800.0)
    is_ocean = np.zeros((5, 7), dtype=bool)
    result = _classify(temp, precip, is_ocean, lat=np.full((5, 7), 40.0))
    assert result.shape == (5, 7)


def test_classify_biomes_covers_a_broad_climate_sweep():
    temps = np.linspace(-35.0, 35.0, 50)
    precips = np.linspace(0.0, 4000.0, 50)
    lats = np.linspace(-80.0, 80.0, 50)
    t, p = np.meshgrid(temps, precips)
    lat = np.repeat(lats[:, None], 50, axis=1)
    is_ocean = np.zeros_like(t, dtype=bool)
    result = biomes.classify_koppen(t, p, lat, np.full_like(t, 0.5))
    # At least half of all 31 Köppen classes should show up in a sweep this wide.
    assert len(np.unique(result)) >= len(biomes.KOPPEN_CODES) // 2


def test_biome_colors_index_aligned_with_biome_names():
    assert len(biomes.BIOME_COLORS) == len(biomes.BIOME_NAMES)
    assert biomes.BIOME_COLORS.shape == (len(biomes.BIOME_NAMES), 3)
    assert biomes.BIOME_COLORS.dtype == np.uint8


def test_ocean_ids_are_exactly_the_pelagic_classes():
    assert biomes.OCEAN_IDS == frozenset(range(len(biomes.KOPPEN_CODES), len(biomes.BIOME_NAMES)))
    assert len(biomes.OCEAN_IDS) == len(biomes.PELAGIC_NAMES)


# --- classify_wetland (unchanged -- still used by geology.py) ------------------------


def test_classify_wetland_requires_flat_low_land():
    temp = np.array([25.0])
    precip = np.array([2500.0])
    is_ocean = np.array([False])
    is_wetland, is_carboniferous = biomes.classify_wetland(temp, precip, np.array([10.0]), np.array([0.05]), is_ocean)
    assert not is_wetland[0] and not is_carboniferous[0]

    is_wetland, is_carboniferous = biomes.classify_wetland(temp, precip, np.array([10.0]), np.array([0.0001]), is_ocean)
    assert is_carboniferous[0] and not is_wetland[0]


def test_classify_wetland_cooler_flat_low_land_is_plain_wetland_not_carboniferous():
    is_wetland, is_carboniferous = biomes.classify_wetland(
        np.array([10.0]), np.array([1500.0]), np.array([5.0]), np.array([0.0001]), np.array([False])
    )
    assert is_wetland[0] and not is_carboniferous[0]


def test_classify_wetland_excludes_ocean():
    is_wetland, is_carboniferous = biomes.classify_wetland(
        np.array([25.0]), np.array([2500.0]), np.array([-5.0]), np.array([0.0001]), np.array([True])
    )
    assert not is_wetland[0] and not is_carboniferous[0]


# --- biome_relative_shade_factor ---------------------------------------------------


def test_biome_relative_shade_factor_matches_shape_of_biome_ids():
    biome_ids = np.array([biomes.TUNDRA, biomes.TUNDRA, biomes.POLAR_SEA_ICE])
    elevation = np.array([100.0, 2000.0, -500.0])
    assert biomes.biome_relative_shade_factor(biome_ids, elevation).shape == (3,)


def test_biome_relative_shade_factor_stays_within_the_amplitude_bounds():
    biome_ids = np.full(4, biomes.TUNDRA)
    elevation = np.array([100.0, 2000.0, 900.0, 50.0])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert np.all(result >= biomes.BIOME_SHADE_MIN - 1e-9)
    assert np.all(result <= biomes.BIOME_SHADE_MAX + 1e-9)


def test_biome_relative_shade_factor_spans_the_full_range_and_is_continuous():
    n = 300
    biome_ids = np.full(n, biomes.koppen_index("BSk"))
    elevation = np.linspace(0.0, 1000.0, n)
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert result.min() == pytest.approx(biomes.BIOME_SHADE_MIN)
    assert result.max() == pytest.approx(biomes.BIOME_SHADE_MAX)
    assert len(np.unique(result)) == n
    assert np.all(np.diff(result) > 0)


def test_biome_relative_shade_factor_gives_equal_elevation_equal_shade():
    # Two big flat plateaus (many exactly-equal values) plus one peak. Cells that share an
    # elevation must share a shade -- no tie-break by array position (the horizontal-banding bug).
    biome_ids = np.full(201, biomes.koppen_index("BWh"))
    elevation = np.concatenate([np.full(100, 200.0), np.full(100, 201.0), [9000.0]])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    assert len(np.unique(result[:100])) == 1
    assert len(np.unique(result[100:200])) == 1
    assert result[0] < result[100] < result[200]


def test_biome_relative_shade_factor_ranks_are_relative_to_each_class_separately():
    n = 60
    biome_ids = np.concatenate([np.full(n, biomes.koppen_index("Af")), np.full(n, biomes.TUNDRA)])
    elevation = np.concatenate([np.linspace(0.0, 10.0, n), np.linspace(0.0, 3000.0, n)])
    result = biomes.biome_relative_shade_factor(biome_ids, elevation)
    for segment in (result[:n], result[n:]):
        assert segment.min() == pytest.approx(biomes.BIOME_SHADE_MIN)
        assert segment.max() == pytest.approx(biomes.BIOME_SHADE_MAX)


# --- grid_slope (unchanged) --------------------------------------------------------


def test_grid_slope_is_zero_on_a_flat_grid():
    lat_deg = np.array([30.0, 0.0, -30.0])
    assert np.all(biomes.grid_slope(np.full((3, 4), 500.0), lat_deg) == 0.0)


def test_grid_slope_matches_a_known_north_south_step():
    lat_deg = np.array([30.0, 0.0, -30.0])
    elevation_m = np.zeros((3, 4))
    elevation_m[1, :] = 1000.0
    slope = biomes.grid_slope(elevation_m, lat_deg)
    dlat_km = (np.pi / 3) * biomes.PLANET_RADIUS_KM
    assert np.allclose(slope[1, :], 1000.0 / (dlat_km * 1000.0))


# --- smooth_biome_field -----------------------------------------------------------


def _uniform_grid(shape, temp, precip, elevation=1000.0, slope=0.01, is_ocean=False):
    return dict(
        temperature_c=np.full(shape, float(temp)),
        precipitation_mm=np.full(shape, float(precip)),
        elevation_m=np.full(shape, float(elevation)),
        slope=np.full(shape, float(slope)),
        is_ocean=np.full(shape, bool(is_ocean)),
    )


_SMOOTH_LAT = 9.0  # tropical, so precipitation contrast alone drives distinct classes (Af/Aw/BS)


def _smooth(d, lat=_SMOOTH_LAT):
    shape = d["temperature_c"].shape
    return biomes.smooth_biome_field(
        d["temperature_c"], d["precipitation_mm"], d["elevation_m"], d["slope"], d["is_ocean"],
        lat_deg=np.full(shape, float(lat)),
    )


def _raw(d, lat=_SMOOTH_LAT):
    # Mirror smooth_biome_field's own continentality/coast-distance derivation so a noop
    # comparison isolates the vote pass, not the geometry inputs.
    shape = d["temperature_c"].shape
    return biomes.classify_biomes(
        d["temperature_c"], d["precipitation_mm"], d["elevation_m"], d["slope"], d["is_ocean"],
        lat_deg=np.full(shape, float(lat)),
        continentality=biomes.grid_continentality(d["is_ocean"]),
        dist_to_land_rad=biomes.grid_dist_to_land_rad(d["is_ocean"]),
    )


def test_smooth_biome_field_rejects_a_flat_array():
    with pytest.raises(ValueError):
        biomes.smooth_biome_field(
            np.zeros(16), np.full(16, 750.0), np.full(16, 1000.0), np.zeros(16), np.zeros(16, dtype=bool),
            lat_deg=np.zeros(16),
        )


def test_smooth_biome_field_preserves_shape():
    d = _uniform_grid((5, 8), 25.0, 1500.0)
    assert _smooth(d).shape == (5, 8)


def test_smooth_biome_field_is_a_noop_on_a_spatially_uniform_climate():
    d = _uniform_grid((6, 9), 25.0, 1500.0)
    assert np.array_equal(_smooth(d), _raw(d))


def test_smooth_biome_field_outvotes_a_lone_speckle():
    d = _uniform_grid((7, 7), 25.0, 1800.0)
    d["precipitation_mm"][3, 3] = 200.0  # one cell shoved into a much drier class
    raw = _raw(d)
    assert raw[3, 3] != raw[0, 0]
    smoothed = _smooth(d)
    assert smoothed[3, 3] == smoothed[0, 0]


def test_smooth_biome_field_leaves_a_clean_two_class_interface_in_place():
    d = _uniform_grid((8, 8), 25.0, 2400.0)
    d["precipitation_mm"][:, 4:] = 500.0
    raw = _raw(d)
    assert np.array_equal(_smooth(d), raw)


def test_smooth_biome_field_keeps_a_small_solid_region():
    d = _uniform_grid((7, 7), 25.0, 2400.0)
    d["precipitation_mm"][2:5, 2:5] = 400.0
    smoothed = _smooth(d)
    assert len(np.unique(smoothed[2:5, 2:5])) == 1
    assert smoothed[2, 2] != smoothed[0, 0]
