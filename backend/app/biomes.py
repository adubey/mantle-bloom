"""Köppen-Geiger climate classification for land, plus a Pelagic-Provinces-of-the-World-style
classification for ocean -- buckets already-computed fields (annual mean temperature, annual
precipitation, elevation/slope, latitude, ocean-surface temperature) into one of a fixed set
of named classes. Like render_image.py's own temperature_colors/humidity_colors this is a
stateless function of values other modules already produce every step -- no new simulation
state, nothing cached here.

The real Köppen scheme keys off *sub-annual* quantities this model never produces -- coldest-
and warmest-month temperature, and the seasonal (summer vs. winter) split of precipitation.
We synthesize those from what we do have: latitude, continentality (distance inland), and the
world's axial tilt drive a seasonal temperature amplitude (`_seasonal_temp_amplitude`);
latitude alone drives a summer precipitation share and a seasonality concentration
(`_precip_season`). A tilt-0 world gets amplitude 0 -- no seasons -- so its `s`/`w`/`d`
classes simply never occur, which is the physically right behavior. The boundary values below
are the standard Köppen thresholds where they apply directly; the synthesis constants are
picked (like erosion.py's own RAIN_EROSION_COEFFICIENT, "tuned by rough order-of-magnitude
reasoning") for a visually sensible spread of climate zones on a generated world, not fit
against a real-world dataset.

Ocean cells can't carry the literal PPOW province names (those are Earth geography -- the
Gulf Stream, the Agulhas Current); instead `classify_pelagic` reproduces PPOW's *abiotic*
hierarchy: a thermal realm (polar / cold-temperate / temperate / subtropical / tropical) from
sea-surface temperature, crossed with a structural zone (open ocean, coastal shelf,
equatorial divergence, sea ice).

`BIOME_NAMES` / `BIOME_COLORS` stay the names the rest of the codebase reads (climate.py's
`VEGETATION_TRANSPIRATION_BY_BIOME`, stats.py's `biome_land_fraction`, render_image.py, the
frontend legend): the Köppen classes first, then the pelagic classes, then -- for
`smooth_biome_field`'s callers only -- the classification is a plain per-cell band lookup.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import distance_transform_edt

from .elevation_lines import PLANET_RADIUS_KM
from .hydrology import GLACIER_ACCUMULATION_TEMP_C

# Earth's real tilt -- the reference the seasonal-amplitude synthesis is calibrated against,
# and the default when a caller has no World handy (kept in sync with world.DEFAULT_AXIAL_TILT_DEG
# by hand, same one-line-constant precedent as ICE_TEMP_C below; not imported, to avoid a
# world -> biomes import edge).
DEFAULT_AXIAL_TILT_DEG = 23.5

# Same threshold hydrology.py's own glacier accumulation logic already uses for "cold enough
# to permanently freeze" -- reused rather than inventing a second, potentially-inconsistent
# cold cutoff, so the Ice Cap class lines up with where the simulation would actually grow a
# glacier.
ICE_TEMP_C = GLACIER_ACCUMULATION_TEMP_C
COLD_TEMP_C = 5.0
# Köppen's tropical / arid `h`-vs-`k` isotherm, and the warmth at which soil biomass
# productivity saturates in geology.py -- one shared constant, same precedent as ICE_TEMP_C.
TROPICAL_TEMP_C = 18.0

# Precipitation band references, still imported by geology.py (soil productivity saturates at
# HUMID_MM) and reused inside classify_wetland below.
ARID_MM = 250.0
SEMI_ARID_MM = 500.0
SUB_HUMID_MM = 1000.0
HUMID_MM = 2000.0

# Wetland/Carboniferous Forest are no longer a displayed climate class (Köppen has no such
# category), but classify_wetland and these constants stay: geology.py's own coal formation
# calls classify_wetland directly, and the map's relief/erosion code references the same
# flat-and-low idea. See classify_wetland's docstring.
WETLAND_MAX_ELEVATION_M = 50.0
WETLAND_MAX_SLOPE = 0.001
WETLAND_MIN_PRECIP_MM = SEMI_ARID_MM
CARBONIFEROUS_MIN_TEMP_C = TROPICAL_TEMP_C
CARBONIFEROUS_MIN_PRECIP_MM = HUMID_MM

# Ocean shallow enough to read as "coastal shelf" rather than open water -- an arc distance
# from the nearest land, matching geology.SHELF_RANGE_RAD's own shelf definition (kept as a
# literal here to avoid a geology <-> biomes import edge; ~2 degrees of arc).
SHELF_RANGE_RAD = np.radians(2.0)

# ---------------------------------------------------------------------------------------
# Class tables
# ---------------------------------------------------------------------------------------
# Every standard Köppen-Geiger 3rd-level code, in the conventional A/B/C/D/E order, each with
# a descriptive display name (the user-facing label -- "Humid Subtropical", not "Cfa"). The
# frontend legend (frontend/src/legendData.ts) groups these into 1st/2nd-level rows for
# display; the classification itself always resolves to one of these.
KOPPEN_CODES = [
    "Af", "Am", "Aw", "As",
    "BWh", "BWk", "BSh", "BSk",
    "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Cfa", "Cfb", "Cfc",
    "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "Dfa", "Dfb", "Dfc", "Dfd",
    "ET", "EF",
]
KOPPEN_NAMES = [
    "Tropical Rainforest", "Tropical Monsoon", "Tropical Savanna", "Tropical Savanna (Dry Summer)",
    "Hot Desert", "Cold Desert", "Hot Semi-Arid", "Cold Semi-Arid",
    "Hot-Summer Mediterranean", "Warm-Summer Mediterranean", "Cold-Summer Mediterranean",
    "Humid Subtropical (Dry Winter)", "Subtropical Highland", "Cold Subtropical Highland",
    "Humid Subtropical", "Oceanic", "Subpolar Oceanic",
    "Mediterranean Continental (Hot Summer)", "Mediterranean Continental (Warm Summer)",
    "Mediterranean Subarctic", "Extremely Cold Mediterranean Subarctic",
    "Monsoon Continental (Hot Summer)", "Monsoon Continental (Warm Summer)",
    "Monsoon Subarctic", "Extremely Cold Monsoon Subarctic",
    "Hot-Summer Humid Continental", "Warm-Summer Humid Continental",
    "Subarctic (Boreal)", "Extremely Cold Subarctic",
    "Tundra", "Ice Cap",
]

# The pelagic (ocean) classes -- PPOW's abiotic realm x zone hierarchy, descriptive names.
PELAGIC_NAMES = [
    "Tropical Open Ocean", "Subtropical Gyre", "Equatorial Divergence", "Tropical Coastal Waters",
    "Temperate Open Ocean", "Temperate Shelf",
    "Cold-Temperate Open Ocean", "Cold-Temperate Shelf",
    "Polar Ocean", "Polar Sea Ice",
]

BIOME_NAMES = KOPPEN_NAMES + PELAGIC_NAMES

_N_KOPPEN = len(KOPPEN_NAMES)


def koppen_index(code: str) -> int:
    """Id (index into BIOME_NAMES/BIOME_COLORS) of a Köppen 3rd-level code, e.g. "Cfb"."""
    return KOPPEN_CODES.index(code)


def pelagic_index(name: str) -> int:
    """Id of a pelagic class by its display name, e.g. "Polar Sea Ice"."""
    return _N_KOPPEN + PELAGIC_NAMES.index(name)


# A few named ids still referenced by name elsewhere / in tests.
TUNDRA = koppen_index("ET")
ICE_CAP = koppen_index("EF")
POLAR_SEA_ICE = pelagic_index("Polar Sea Ice")
# The set of ids that are ocean, not land -- stats.py excludes these from biome_land_fraction,
# render_image.py's Combined view colors them by pelagic province rather than Köppen.
OCEAN_IDS = frozenset(range(_N_KOPPEN, len(BIOME_NAMES)))

# Index-aligned with BIOME_NAMES. Natural "true colors" -- the muted, vegetation-and-season-
# averaged gamut of the reference image the user supplied (Earth rendered with each Köppen
# class painted its own average true color): dark greens through the wet tropics, tans across
# the deserts, olive/khaki through the mid-latitude grasslands and the snow-averaged boreal
# belt, near-white at the ice caps; the oceans run deep navy in the oligotrophic gyres,
# brightening to blue over the shelves and toward the poles, near-white where sea ice covers
# them. Anchor values (Sahara/Amazon/Congo/Siberia/deep-Pacific/Antarctica) were sampled from
# that image; the finer 3rd-level variants are interpolated within the same gamut and are the
# obvious place to hand-tune against a render.
_KOPPEN_COLORS = [
    (26, 51, 21),    # Af  Tropical Rainforest
    (34, 66, 29),    # Am  Tropical Monsoon
    (74, 92, 44),    # Aw  Tropical Savanna
    (92, 103, 52),   # As  Tropical Savanna (Dry Summer)
    (176, 141, 99),  # BWh Hot Desert
    (150, 136, 104), # BWk Cold Desert
    (156, 134, 78),  # BSh Hot Semi-Arid
    (124, 128, 82),  # BSk Cold Semi-Arid
    (150, 138, 66),  # Csa Hot-Summer Mediterranean
    (128, 136, 76),  # Csb Warm-Summer Mediterranean
    (112, 126, 86),  # Csc Cold-Summer Mediterranean
    (94, 116, 54),   # Cwa Humid Subtropical (Dry Winter)
    (100, 122, 66),  # Cwb Subtropical Highland
    (106, 124, 80),  # Cwc Cold Subtropical Highland
    (58, 96, 44),    # Cfa Humid Subtropical
    (72, 104, 56),   # Cfb Oceanic
    (92, 112, 76),   # Cfc Subpolar Oceanic
    (132, 126, 80),  # Dsa Mediterranean Continental (Hot Summer)
    (124, 124, 84),  # Dsb Mediterranean Continental (Warm Summer)
    (120, 122, 92),  # Dsc Mediterranean Subarctic
    (128, 126, 104), # Dsd Extremely Cold Mediterranean Subarctic
    (96, 110, 56),   # Dwa Monsoon Continental (Hot Summer)
    (104, 116, 68),  # Dwb Monsoon Continental (Warm Summer)
    (116, 118, 84),  # Dwc Monsoon Subarctic
    (128, 122, 100), # Dwd Extremely Cold Monsoon Subarctic
    (84, 104, 52),   # Dfa Hot-Summer Humid Continental
    (100, 114, 66),  # Dfb Warm-Summer Humid Continental
    (122, 112, 80),  # Dfc Subarctic (Boreal)
    (140, 128, 100), # Dfd Extremely Cold Subarctic
    (150, 152, 132), # ET  Tundra
    (234, 238, 240), # EF  Ice Cap
]
_PELAGIC_COLORS = [
    (7, 17, 43),     # Tropical Open Ocean
    (4, 12, 36),     # Subtropical Gyre -- darkest, oligotrophic
    (14, 40, 78),    # Equatorial Divergence -- upwelling, slightly brighter
    (34, 82, 122),   # Tropical Coastal Waters
    (10, 26, 58),    # Temperate Open Ocean
    (30, 70, 106),   # Temperate Shelf
    (18, 42, 72),    # Cold-Temperate Open Ocean
    (38, 74, 104),   # Cold-Temperate Shelf
    (32, 62, 94),    # Polar Ocean
    (226, 233, 237), # Polar Sea Ice
]

BIOME_COLORS = np.array(_KOPPEN_COLORS + _PELAGIC_COLORS, dtype=np.uint8)
assert BIOME_COLORS.shape == (len(BIOME_NAMES), 3)


# ---------------------------------------------------------------------------------------
# Relief shading (unchanged behavior -- see render_image.py's _render_combined_view)
# ---------------------------------------------------------------------------------------
# Peak-to-trough brightness swing biome_relative_shade_factor spreads each class's flat
# BIOME_COLORS entry across, keyed to a cell's elevation *rank among other cells of the same
# class* -- a wide +-25% shaded-relief look. Combined writes the class id into the render's
# alpha channel (see render_image.COMBINED_*_ID_CODE), so a wide color spread never confuses
# the frontend's click-to-highlight.
BIOME_SHADE_AMPLITUDE = 0.25
BIOME_SHADE_MIN = 1.0 - BIOME_SHADE_AMPLITUDE
BIOME_SHADE_MAX = 1.0 + BIOME_SHADE_AMPLITUDE


def biome_relative_shade_factor(biome_ids: np.ndarray, elevation_m: np.ndarray) -> np.ndarray:
    """Per-cell brightness multiplier, same shape as `biome_ids`: a continuous ramp from
    BIOME_SHADE_MIN at a class's lowest-elevation cell to BIOME_SHADE_MAX at its highest,
    linear in that cell's elevation *rank among same-class cells only*. Rank-based (an even
    spread over the sorted order) rather than a quantile-of-value ramp, so a class with one
    outlier peak and otherwise flat terrain still spans the full range smoothly. Multiply this
    elementwise into a class's flat BIOME_COLORS entry to shade it."""
    biome_ids = np.asarray(biome_ids)
    elevation_m = np.asarray(elevation_m)
    flat_ids = biome_ids.reshape(-1)
    flat_elevation = elevation_m.reshape(-1)

    factor = np.ones(flat_ids.shape, dtype=float)
    for biome_id in np.unique(flat_ids):
        mask = flat_ids == biome_id
        n = int(mask.sum())
        if n <= 1:
            continue
        order = np.argsort(flat_elevation[mask])
        ranks = np.empty(n, dtype=np.int64)
        ranks[order] = np.arange(n)
        frac = ranks / (n - 1)
        factor[mask] = BIOME_SHADE_MIN + (BIOME_SHADE_MAX - BIOME_SHADE_MIN) * frac

    return factor.reshape(biome_ids.shape)


def grid_slope(elevation_m: np.ndarray, lat_deg: np.ndarray) -> np.ndarray:
    """Dimensionless rise/run slope on a fixed (H, W) lat/lon grid -- real elevation difference
    to each cell's north/south or east/west neighbor (whichever is steeper), divided by that
    neighbor's real great-circle spacing in meters (longitude narrowed by cos(lat)). Feeds
    classify_wetland's own WETLAND_MAX_SLOPE cutoff. np.roll wraps at the poles too (a minor,
    visually inconsequential artifact right at the map's own poles)."""
    grid_h, grid_w = elevation_m.shape
    dlat_km = (np.pi / grid_h) * PLANET_RADIUS_KM
    dlon_km = np.maximum((2 * np.pi / grid_w) * PLANET_RADIUS_KM * np.cos(np.radians(lat_deg))[:, None], 1.0)
    d_ns = np.abs(elevation_m - np.roll(elevation_m, 1, axis=0)) / (dlat_km * 1000.0)
    d_ew = np.abs(elevation_m - np.roll(elevation_m, 1, axis=1)) / (dlon_km * 1000.0)
    return np.maximum(d_ns, d_ew)


def classify_wetland(
    temperature_c: np.ndarray, precipitation_mm: np.ndarray, elevation_m: np.ndarray, slope: np.ndarray, is_ocean: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(is_wetland, is_carboniferous_forest) boolean arrays, same shape as the inputs -- used
    by geology.py's per-node coal-formation accumulation (warm wet swamp forest forms coal
    fastest). No longer a displayed climate class, but the flat-and-low predicate stays here
    as the one shared definition. Both require flat, low-lying land -- a floodplain, delta, or
    coastal marsh, not an upland bog; Carboniferous Forest is the warm
    (>= CARBONIFEROUS_MIN_TEMP_C), very wet (>= CARBONIFEROUS_MIN_PRECIP_MM) subtype, Wetland
    the cooler/drier remainder (>= ICE_TEMP_C, >= WETLAND_MIN_PRECIP_MM)."""
    elevation_m = np.asarray(elevation_m)
    slope = np.asarray(slope)
    temp = np.asarray(temperature_c)
    precip = np.asarray(precipitation_mm)
    is_ocean = np.asarray(is_ocean)

    flat_low = ~is_ocean & (elevation_m > 0) & (elevation_m <= WETLAND_MAX_ELEVATION_M) & (slope <= WETLAND_MAX_SLOPE)
    is_carboniferous = flat_low & (temp >= CARBONIFEROUS_MIN_TEMP_C) & (precip >= CARBONIFEROUS_MIN_PRECIP_MM)
    is_wetland = flat_low & (temp >= ICE_TEMP_C) & (precip >= WETLAND_MIN_PRECIP_MM) & ~is_carboniferous
    return is_wetland, is_carboniferous


# ---------------------------------------------------------------------------------------
# Synthesized seasonality
# ---------------------------------------------------------------------------------------
# The model has annual-mean temperature and annual precipitation only. Köppen needs the
# coldest/warmest month and the summer-vs-winter precipitation split; we synthesize both from
# latitude + continentality + axial tilt. All constants here are visually tuned (module
# docstring), calibrated so an Earth-like world lands its major zones roughly where Earth's
# are.

# Seasonal temperature half-amplitude (mean-to-peak): grows with |latitude|, is amplified
# deep in a continental interior (a maritime cell sees a much smaller swing at the same
# latitude), and scales with how tilted the world is. A tilt-0 world has amplitude 0
# everywhere -- so its `s`/`w`/`d` classes never occur, which is correct.
_SEASON_LAT_AMPLITUDE_C = 15.5
_SEASON_LAT_EXPONENT = 1.3
_SEASON_CONTINENTAL_GAIN = 0.9   # continental interior multiplies the latitudinal swing by up to this
_SEASON_CONTINENTAL_OFFSET_C = 8.0  # plus a flat interior boost (winter radiative cooling)
_SEASON_TILT_REFERENCE_DEG = 23.5
# Earth's largest real mean-to-peak seasonal swing (interior NE Siberia) is ~33 C -- cap the
# synthesized value near there so an extreme cold + fully-continental cell doesn't get an
# unphysical +30 C summer bump lifting it out of the polar classes.
_SEASON_MAX_AMPLITUDE_C = 33.0


def _seasonal_temp_amplitude(lat_deg: np.ndarray, continentality: np.ndarray, axial_tilt_deg: float) -> np.ndarray:
    """Half of (warmest-month mean - coldest-month mean), same shape as `lat_deg`.
    `continentality` is 0 at the coast, 1 deep in a landmass' interior."""
    lat_deg = np.asarray(lat_deg, dtype=float)
    cont = np.clip(np.asarray(continentality, dtype=float), 0.0, 1.0)
    tilt_factor = np.clip(axial_tilt_deg / _SEASON_TILT_REFERENCE_DEG, 0.0, 1.7)
    lat_term = _SEASON_LAT_AMPLITUDE_C * (np.abs(lat_deg) / 90.0) ** _SEASON_LAT_EXPONENT
    amplitude = lat_term * (1.0 + _SEASON_CONTINENTAL_GAIN * cont) + _SEASON_CONTINENTAL_OFFSET_C * cont
    return tilt_factor * np.minimum(amplitude, _SEASON_MAX_AMPLITUDE_C)


def _precip_season(lat_deg: np.ndarray, season_strength: float | np.ndarray = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """(summer_share, concentration), both same shape as `lat_deg`:
    `summer_share` is the fraction of annual precipitation falling in the warm half-year
    (>= 0.7 => Köppen treats the regime as summer-wet, <= 0.3 as winter-wet); `concentration`
    is how peaked precipitation is into one season (0 = even every month, up toward 1 = one
    short wet season). Latitude-driven -- deep tropics monsoonal (summer-wet, peaked),
    subtropics dry-summer Mediterranean, mid-latitudes even -- and scaled by `season_strength`
    (the axial-tilt factor: with no tilt the ITCZ doesn't migrate, so there is no wet/dry
    season at all and every regime collapses toward 'f')."""
    lat = np.abs(np.asarray(lat_deg, dtype=float))
    # The monsoon belt fades to ~0 right at the equator (the ITCZ crosses twice a year there
    # -- wettest, least seasonal) and peaks in the outer tropics / low subtropics.
    equator_ramp = np.clip(lat / 8.0, 0.0, 1.0)
    monsoon_wet = np.exp(-((lat - 14.0) / 13.0) ** 2) * equator_ramp
    # A narrow subtropical dry-summer dip. Latitude alone can't tell a west coast
    # (Mediterranean) from an east coast (humid) at the same latitude -- this dip is
    # deliberately shallow, and classify_koppen gates `s` further on maritime continentality
    # and modest precipitation.
    subtropical_dry_summer = np.exp(-((lat - 35.0) / 7.0) ** 2)
    summer_share = 0.5 + season_strength * (0.36 * monsoon_wet - 0.24 * subtropical_dry_summer)
    summer_share = np.clip(summer_share, 0.12, 0.9)

    monsoon = np.exp(-((lat - 14.0) / 15.0) ** 2) * equator_ramp
    subtropical = np.exp(-((lat - 35.0) / 10.0) ** 2)
    concentration = np.clip(0.10 + season_strength * (0.66 * monsoon + 0.28 * subtropical), 0.05, 0.9)
    return summer_share, concentration


def _months_above(t_annual: np.ndarray, amplitude: np.ndarray, threshold_c: float) -> np.ndarray:
    """Number of months (0..12) whose mean temperature exceeds `threshold_c`, for a sinusoidal
    annual cycle of the given mean and half-amplitude. Closed form: the fraction of the cycle
    above the threshold is arccos(x)/pi with x = (threshold - mean) / amplitude."""
    amp = np.maximum(np.asarray(amplitude, dtype=float), 1e-6)
    x = np.clip((threshold_c - np.asarray(t_annual, dtype=float)) / amp, -1.0, 1.0)
    return 12.0 * np.arccos(x) / np.pi


def _seasonal_extremes(
    t_annual: np.ndarray, precip_mm: np.ndarray, lat_deg: np.ndarray, continentality: np.ndarray, axial_tilt_deg: float
):
    """Everything the Köppen decision needs beyond the raw annual pair, all same-shape:
    (t_cold, t_warm, amplitude, driest_month_mm, wettest_summer_month_mm,
    driest_summer_month_mm, wettest_winter_month_mm, summer_share)."""
    amplitude = _seasonal_temp_amplitude(lat_deg, continentality, axial_tilt_deg)
    t_annual = np.asarray(t_annual, dtype=float)
    precip_mm = np.asarray(precip_mm, dtype=float)
    t_cold = t_annual - amplitude
    t_warm = t_annual + amplitude

    season_strength = np.clip(axial_tilt_deg / _SEASON_TILT_REFERENCE_DEG, 0.0, 1.3)
    summer_share, conc = _precip_season(lat_deg, season_strength)
    monthly_mean = precip_mm / 12.0
    # A season's monthly rate = annual_share / 6 months, then spread within the year by the
    # concentration: the wet season's peak month runs well above its mean, the dry season's
    # low month well below.
    summer_month = precip_mm * summer_share / 6.0
    winter_month = precip_mm * (1.0 - summer_share) / 6.0
    wet_is_summer = summer_share >= 0.5
    wettest_month = np.where(wet_is_summer, summer_month, winter_month) * (1.0 + 4.0 * conc)
    driest_other = np.where(wet_is_summer, winter_month, summer_month) * (1.0 - conc)
    driest_month = np.minimum(driest_other, monthly_mean * (1.0 - conc))

    wettest_summer_month = summer_month * (1.0 + 4.0 * conc * wet_is_summer)
    driest_summer_month = summer_month * (1.0 - conc)
    wettest_winter_month = winter_month * (1.0 + 4.0 * conc * ~wet_is_summer)
    return (
        t_cold, t_warm, amplitude,
        driest_month, wettest_summer_month, driest_summer_month, wettest_winter_month, summer_share,
    )


# ---------------------------------------------------------------------------------------
# Köppen classification
# ---------------------------------------------------------------------------------------


def classify_koppen(
    t_annual: np.ndarray,
    precip_mm: np.ndarray,
    lat_deg: np.ndarray,
    continentality: np.ndarray | None = None,
    axial_tilt_deg: float = DEFAULT_AXIAL_TILT_DEG,
) -> np.ndarray:
    """Per-cell Köppen id (index into BIOME_NAMES/KOPPEN_CODES). Standard A/B/C/D/E thresholds,
    with the coldest/warmest-month temperature and seasonal precipitation split synthesized
    (see the block above). All array inputs must share a shape; `continentality` defaults to
    0 (everywhere coastal -- no continental seasonal boost) when a caller has no landmass
    geometry on hand."""
    t_annual = np.asarray(t_annual, dtype=float)
    precip_mm = np.asarray(precip_mm, dtype=float)
    lat_deg = np.asarray(lat_deg, dtype=float)
    if continentality is None:
        continentality = np.zeros_like(t_annual)

    (
        t_cold, t_warm, amplitude,
        driest_month, wettest_summer_month, driest_summer_month, wettest_winter_month, summer_share,
    ) = _seasonal_extremes(t_annual, precip_mm, lat_deg, continentality, axial_tilt_deg)

    # --- Arid (B) dryness threshold: 20*T + {280 summer-wet | 140 even | 0 winter-wet} ---
    pth_offset = np.where(summer_share >= 0.7, 280.0, np.where(summer_share <= 0.3, 0.0, 140.0))
    pth = 20.0 * t_annual + pth_offset
    is_b = precip_mm < pth
    is_bw = is_b & (precip_mm < 0.5 * pth)
    is_hot = t_annual >= TROPICAL_TEMP_C

    # --- Polar (E) ---
    is_e = t_warm < 10.0
    is_ef = is_e & (t_warm < 0.0)

    # --- Tropical (A): coldest month >= 18 C ---
    is_a = (t_cold >= 18.0) & ~is_b
    af = is_a & (driest_month >= 60.0)
    am = is_a & ~af & (driest_month >= 100.0 - precip_mm / 25.0)
    a_dry = is_a & ~af & ~am
    aw = a_dry & (summer_share >= 0.5)  # dry season in the low-sun half -> "dry winter"
    a_s = a_dry & ~aw

    # --- Second letter (s / w / f) for C and D ---
    # Köppen `s`: driest summer month < 40 mm and < 1/3 the wettest winter month. Further
    # gated to a maritime (low-continentality), modest-precipitation subtropical band -- a
    # latitude-only model otherwise scatters spurious Mediterranean across every dry
    # subtropical interior and highland (module docstring).
    mediterranean_band = (np.abs(lat_deg) >= 27.0) & (np.abs(lat_deg) <= 46.0) & (continentality < 0.4)
    summer_dry = (
        mediterranean_band
        & (precip_mm < 900.0)
        & (summer_share < 0.45)
        & (driest_summer_month < 40.0)
        & (driest_summer_month < wettest_winter_month / 3.0)
    )
    # Köppen `w`: driest month (a winter month once the wet season is summer) < 1/10 the
    # wettest summer month.
    winter_dry = (driest_month < wettest_summer_month / 10.0) & (summer_share >= 0.62)

    # --- Third letter (a / b / c / d) for C and D ---
    warm_a = t_warm >= 22.0
    months_10 = _months_above(t_annual, amplitude, 10.0)
    warm_b = ~warm_a & (months_10 >= 4.0)
    very_cold_d = t_cold < -38.0

    is_cd = ~is_a & ~is_b & ~is_e & (t_warm >= 10.0)
    is_d = is_cd & (t_cold < 0.0)
    is_c = is_cd & ~is_d

    def third(is_grp):
        return {
            "a": is_grp & warm_a,
            "b": is_grp & warm_b,
            "c": is_grp & ~warm_a & ~warm_b,
        }

    c_s, c_w = is_c & summer_dry, is_c & winter_dry & ~summer_dry
    c_f = is_c & ~c_s & ~c_w
    d_s, d_w = is_d & summer_dry, is_d & winter_dry & ~summer_dry
    d_f = is_d & ~d_s & ~d_w

    cs, cw, cf = third(c_s), third(c_w), third(c_f)
    ds, dw, df = third(d_s), third(d_w), third(d_f)

    # B (aridity) is assigned before E/A/C/D -- a cold-enough desert is still BWk, not EF
    # (Peel et al. 2007's own precedence).
    condlist = [
        is_bw & is_hot, is_bw, is_b & is_hot, is_b,
        is_ef, is_e,
        af, am, aw, a_s,
        cs["a"], cs["b"], cs["c"],
        cw["a"], cw["b"], cw["c"],
        cf["a"], cf["b"], cf["c"],
        ds["a"], ds["b"], ds["c"] & ~very_cold_d, ds["c"] & very_cold_d,
        dw["a"], dw["b"], dw["c"] & ~very_cold_d, dw["c"] & very_cold_d,
        df["a"], df["b"], df["c"] & ~very_cold_d, df["c"] & very_cold_d,
    ]
    choicelist = [
        koppen_index("BWh"), koppen_index("BWk"), koppen_index("BSh"), koppen_index("BSk"),
        koppen_index("EF"), koppen_index("ET"),
        koppen_index("Af"), koppen_index("Am"), koppen_index("Aw"), koppen_index("As"),
        koppen_index("Csa"), koppen_index("Csb"), koppen_index("Csc"),
        koppen_index("Cwa"), koppen_index("Cwb"), koppen_index("Cwc"),
        koppen_index("Cfa"), koppen_index("Cfb"), koppen_index("Cfc"),
        koppen_index("Dsa"), koppen_index("Dsb"), koppen_index("Dsc"), koppen_index("Dsd"),
        koppen_index("Dwa"), koppen_index("Dwb"), koppen_index("Dwc"), koppen_index("Dwd"),
        koppen_index("Dfa"), koppen_index("Dfb"), koppen_index("Dfc"), koppen_index("Dfd"),
    ]
    # Default: anything the bands above didn't catch is a cold maritime remainder -> Cfc.
    return np.select(condlist, choicelist, default=koppen_index("Cfc")).astype(np.int64)


# ---------------------------------------------------------------------------------------
# Pelagic (ocean) classification
# ---------------------------------------------------------------------------------------

# Sea-surface-temperature realm cutoffs (deg C): polar / cold-temperate / temperate /
# subtropical / tropical.
_PELAGIC_SST_EDGES = (2.0, 12.0, 19.0, 24.0)
_EQUATORIAL_LAT_DEG = 6.0


def classify_pelagic(
    ocean_temp_c: np.ndarray,
    lat_deg: np.ndarray,
    dist_to_land_rad: np.ndarray | None = None,
    has_sea_ice: np.ndarray | None = None,
) -> np.ndarray:
    """Per-cell pelagic id (index into BIOME_NAMES) for ocean cells. Thermal realm from SST x
    structural zone (coastal shelf / equatorial divergence / open ocean / sea ice). Non-ocean
    cells get an arbitrary id here -- callers mask with `is_ocean` (see classify_biomes)."""
    sst = np.asarray(ocean_temp_c, dtype=float)
    lat = np.abs(np.asarray(lat_deg, dtype=float))
    if dist_to_land_rad is None:
        dist_to_land_rad = np.full(sst.shape, np.inf)
    else:
        dist_to_land_rad = np.asarray(dist_to_land_rad, dtype=float)
    if has_sea_ice is None:
        has_sea_ice = np.zeros(sst.shape, dtype=bool)
    else:
        has_sea_ice = np.asarray(has_sea_ice, dtype=bool)

    is_ice = has_sea_ice | (sst < -1.0)
    is_shelf = dist_to_land_rad <= SHELF_RANGE_RAD
    e0, e1, e2, e3 = _PELAGIC_SST_EDGES
    polar = sst < e0
    cold_temp = (sst >= e0) & (sst < e1)
    temp = (sst >= e1) & (sst < e2)
    subtropical = (sst >= e2) & (sst < e3)
    tropical = sst >= e2  # subtropical + tropical share the "warm" pelagic classes below

    p = pelagic_index
    condlist = [
        is_ice,
        polar,
        cold_temp & is_shelf,
        cold_temp,
        temp & is_shelf,
        temp,
        tropical & is_shelf,
        tropical & (lat <= _EQUATORIAL_LAT_DEG),
        subtropical,
        tropical,
    ]
    choicelist = [
        p("Polar Sea Ice"),
        p("Polar Ocean"),
        p("Cold-Temperate Shelf"),
        p("Cold-Temperate Open Ocean"),
        p("Temperate Shelf"),
        p("Temperate Open Ocean"),
        p("Tropical Coastal Waters"),
        p("Equatorial Divergence"),
        p("Subtropical Gyre"),
        p("Tropical Open Ocean"),
    ]
    return np.select(condlist, choicelist, default=p("Temperate Open Ocean")).astype(np.int64)


# ---------------------------------------------------------------------------------------
# Combined entry point
# ---------------------------------------------------------------------------------------


def classify_biomes(
    temperature_c: np.ndarray,
    precipitation_mm: np.ndarray,
    elevation_m: np.ndarray,
    slope: np.ndarray,
    is_ocean: np.ndarray,
    sea_level_m: float = 0.0,
    *,
    lat_deg: np.ndarray,
    axial_tilt_deg: float = DEFAULT_AXIAL_TILT_DEG,
    continentality: np.ndarray | None = None,
    dist_to_land_rad: np.ndarray | None = None,
    has_sea_ice: np.ndarray | None = None,
) -> np.ndarray:
    """Per-cell class id: land cells classified by classify_koppen, ocean cells by
    classify_pelagic, with `is_ocean` settling the split. Stays a pure any-shape band lookup
    (its unit tests and geology.py's classify_wetland call depend on that); the stateful
    spatial cleanup lives in smooth_biome_field below. `lat_deg` must broadcast against the
    other inputs; `continentality` / `dist_to_land_rad` / `has_sea_ice` are optional (the
    grid-only smoothing pass computes them -- a flat-array caller may omit them).

    `elevation_m` / `slope` / `sea_level_m` are accepted for signature stability with the
    former Whittaker classifier and its callers; they no longer feed the classification
    directly (Köppen is a pure climate scheme) but are still threaded through so
    smooth_biome_field's blurred-elevation substitution has somewhere to land if reintroduced."""
    temp = np.asarray(temperature_c, dtype=float)
    precip = np.asarray(precipitation_mm, dtype=float)
    is_ocean = np.asarray(is_ocean)
    lat_deg = np.broadcast_to(np.asarray(lat_deg, dtype=float), temp.shape)

    land_ids = classify_koppen(temp, precip, lat_deg, continentality, axial_tilt_deg)
    ocean_ids = classify_pelagic(temp, lat_deg, dist_to_land_rad, has_sea_ice)
    return np.where(is_ocean, ocean_ids, land_ids).astype(np.int64)


# ---------------------------------------------------------------------------------------
# Boundary cleanup: smooth_biome_field
# ---------------------------------------------------------------------------------------
# classify_biomes buckets continuous fields into hard-edged bands. Wherever the underlying
# field runs nearly parallel to a band cutoff, adjacent render cells flip back and forth
# across it and the un-blurred Biome view reads as a dithered smear along every interface plus
# scattered single-cell speckle. `smooth_biome_field` is a stateless cleanup pass the
# render/stats callers apply on top of classify_biomes: a cell that already disagrees with a
# strong majority of its 8 neighbours (an edge or speckle cell) can be out-voted by that
# majority. A cell comfortably inside a region is never touched.

# Adopt the most common non-ocean neighbour class only if at least this fraction of a cell's
# valid (non-ocean) neighbours agree -- 0.65 means 6+ of 8. Above 0.5 so a genuine two-class
# interface (roughly half and half) is left in place.
BIOME_VOTE_MIN_NEIGHBOUR_FRACTION = 0.65
BIOME_VOTE_ITERATIONS = 2

# A land cell disagreeing with at least this many of its 8 land neighbours is an outnumbered
# speckle/edge cell eligible to be out-voted. High (a clean straight interface has ~3/8
# disagreeing, a diagonal ~2) so region edges and necks stay intact.
BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS = 6

# Continentality (0 at the coast, 1 in the interior) saturates at this arc distance inland.
CONTINENTALITY_SATURATION_RAD = np.radians(12.0)


def _equator_cell_km(grid_w: int) -> float:
    return 2.0 * np.pi * PLANET_RADIUS_KM / grid_w


def grid_continentality(is_ocean: np.ndarray) -> np.ndarray:
    """(H, W) continentality: 0 over ocean and at every coast, ramping to 1 by
    CONTINENTALITY_SATURATION_RAD inland. Distance transform of the land mask, in cells,
    converted to arc via the grid's own resolution (longitude compression at high latitude is
    ignored -- a visually inconsequential over-estimate near the poles, same tradeoff the rest
    of this module already accepts)."""
    land = ~np.asarray(is_ocean)
    if not land.any():
        return np.zeros(land.shape, dtype=float)
    dist_rad = distance_transform_edt(land) * (_equator_cell_km(land.shape[1]) / PLANET_RADIUS_KM)
    return np.clip(dist_rad / CONTINENTALITY_SATURATION_RAD, 0.0, 1.0)


def grid_dist_to_land_rad(is_ocean: np.ndarray) -> np.ndarray:
    """(H, W) arc distance from each ocean cell to the nearest land cell (0 over land)."""
    ocean = np.asarray(is_ocean)
    if not ocean.any() or ocean.all():
        return np.zeros(ocean.shape, dtype=float)
    return distance_transform_edt(ocean) * (_equator_cell_km(ocean.shape[1]) / PLANET_RADIUS_KM)


_NEIGHBOUR_OFFSETS = tuple((dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0))


def _disagreeing_neighbour_count(biome_ids: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """Per-cell count of the 8 neighbours that are land and a different class."""
    count = np.zeros(biome_ids.shape, dtype=np.int16)
    for dy, dx in _NEIGHBOUR_OFFSETS:
        nb = np.roll(np.roll(biome_ids, dy, axis=0), dx, axis=1)
        nb_ocean = np.roll(np.roll(is_ocean, dy, axis=0), dx, axis=1)
        count += ((nb != biome_ids) & ~nb_ocean).astype(np.int16)
    return count


def _neighbour_vote(biome_ids: np.ndarray, eligible: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """One majority-vote pass: each `eligible` land cell adopts the most common land class
    among its 8 neighbours when at least BIOME_VOTE_MIN_NEIGHBOUR_FRACTION of its valid
    (non-ocean) neighbours agree and it differs from the cell's current class."""
    stack = np.stack([np.roll(np.roll(biome_ids, dy, axis=0), dx, axis=1) for dy, dx in _NEIGHBOUR_OFFSETS])
    ocean_stack = np.stack([np.roll(np.roll(is_ocean, dy, axis=0), dx, axis=1) for dy, dx in _NEIGHBOUR_OFFSETS])
    stack = np.where(ocean_stack, np.int16(-1), stack.astype(np.int16))
    valid = (~ocean_stack).sum(axis=0).astype(np.int16)

    stack.sort(axis=0)
    best_count = np.zeros(biome_ids.shape, dtype=np.int16)
    best_biome = np.full(biome_ids.shape, -1, dtype=np.int16)
    run = np.zeros(biome_ids.shape, dtype=np.int16)
    prev = np.full(biome_ids.shape, -2, dtype=np.int16)
    for k in range(stack.shape[0]):
        cur = stack[k]
        run = np.where(cur == prev, run + 1, 1)
        better = (cur >= 0) & (run > best_count)
        best_count = np.where(better, run, best_count)
        best_biome = np.where(better, cur, best_biome)
        prev = cur

    threshold = np.ceil(BIOME_VOTE_MIN_NEIGHBOUR_FRACTION * valid).astype(np.int16)
    take = eligible & ~is_ocean & (valid > 0) & (best_count >= threshold) & (best_biome != biome_ids)
    return np.where(take, best_biome.astype(biome_ids.dtype), biome_ids)


def smooth_biome_field(
    temperature_c: np.ndarray,
    precipitation_mm: np.ndarray,
    elevation_m: np.ndarray,
    slope: np.ndarray,
    is_ocean: np.ndarray,
    sea_level_m: float = 0.0,
    *,
    lat_deg: np.ndarray,
    axial_tilt_deg: float = DEFAULT_AXIAL_TILT_DEG,
    glacier_depth_m: np.ndarray | None = None,
) -> np.ndarray:
    """classify_biomes plus the stateless boundary-cleanup pass -- this is the classification
    the Biome/Combined map views and /world/stats use, while `classify_biomes` stays the raw
    per-cell primitive. Inputs must all be the same 2D (H, W) grid; `lat_deg` may be a (H,)
    row vector or a full grid. Computes continentality and coast distance from `is_ocean`
    itself, and sea-ice cover from `glacier_depth_m` if given."""
    temperature_c = np.asarray(temperature_c)
    precipitation_mm = np.asarray(precipitation_mm)
    elevation_m = np.asarray(elevation_m)
    slope = np.asarray(slope)
    is_ocean = np.asarray(is_ocean)
    if elevation_m.ndim != 2:
        raise ValueError("smooth_biome_field needs 2D grid inputs; use classify_biomes for a flat array")

    grid_h, grid_w = elevation_m.shape
    lat_arr = np.asarray(lat_deg, dtype=float)
    lat_grid = np.broadcast_to(lat_arr.reshape(-1, 1) if lat_arr.ndim == 1 else lat_arr, (grid_h, grid_w))

    continentality = grid_continentality(is_ocean)
    dist_to_land_rad = grid_dist_to_land_rad(is_ocean)
    has_sea_ice = None
    if glacier_depth_m is not None:
        has_sea_ice = is_ocean & (np.asarray(glacier_depth_m) > 0.0)

    biome_ids = classify_biomes(
        temperature_c, precipitation_mm, elevation_m, slope, is_ocean, sea_level_m,
        lat_deg=lat_grid, axial_tilt_deg=axial_tilt_deg, continentality=continentality,
        dist_to_land_rad=dist_to_land_rad, has_sea_ice=has_sea_ice,
    )

    land = ~is_ocean
    eligible = land & (_disagreeing_neighbour_count(biome_ids, is_ocean) >= BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS)
    for _ in range(BIOME_VOTE_ITERATIONS):
        biome_ids = _neighbour_vote(biome_ids, eligible, is_ocean)
    return biome_ids
