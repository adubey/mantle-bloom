"""Whittaker-diagram-inspired biome classification: buckets a (temperature, precipitation)
pair -- plus elevation/slope for three biomes that need more than climate alone (Wetland,
Carboniferous Forest, Intertidal Zone) -- into one of a fixed set of named biomes, purely from
already-computed fields -- no new simulation state, nothing cached anywhere, since (like
render_image.py's own temperature_colors/humidity_colors) this is a stateless function of
values other modules already produce every step. The boundary values below are a
simplification, same spirit as this codebase's other openly-approximate constants (e.g.
erosion.py's own RAIN_EROSION_COEFFICIENT, "tuned by rough order-of-magnitude reasoning" not
lifted from a literal scientific reference) -- picked for a visually sensible spread of named
regions on a generated world, not fit against any specific real-world dataset.
"""

from __future__ import annotations

import numpy as np

from .elevation_lines import PLANET_RADIUS_KM
from .hydrology import GLACIER_ACCUMULATION_TEMP_C

# Same threshold hydrology.py's own glacier accumulation logic already uses for "cold enough
# to permanently freeze" -- reused rather than inventing a second, potentially-inconsistent
# cold cutoff, so this view's Ice region lines up with where the simulation would actually
# grow a glacier.
ICE_TEMP_C = GLACIER_ACCUMULATION_TEMP_C
COLD_TEMP_C = 5.0
TROPICAL_TEMP_C = 18.0

ARID_MM = 250.0
SEMI_ARID_MM = 500.0
SUB_HUMID_MM = 1000.0
HUMID_MM = 2000.0

# Wetland/Carboniferous Forest need two axes the rest of this module's classification doesn't:
# elevation (a wetland sits right near sea level -- a floodplain, delta, or coastal marsh, not
# an upland bog) and slope (flat enough for water to pool/waterlog rather than drain away).
# Both are real quantities this codebase already computes elsewhere (erosion.py's own
# node-cloud slope, or grid_slope below's grid analogue for the map view) -- not new
# simulation state. WETLAND_MAX_SLOPE is a dimensionless rise/run, the same convention
# erosion.py's own slope uses, picked (like every other cutoff in this module) for a visually
# sensible result rather than fit against a dataset.
WETLAND_MAX_ELEVATION_M = 50.0
WETLAND_MAX_SLOPE = 0.001
WETLAND_MIN_PRECIP_MM = SEMI_ARID_MM
# Carboniferous Forest is the warm, very-wet subtype of wetland -- the real geological analog
# of Carboniferous/Permian coal swamps, which were predominantly tropical lowland swamp
# forest, not the cooler bogs/marshes "Wetland" alone covers. geology.py's own coal formation
# reuses this exact split (see classify_wetland below) so the map's Carboniferous Forest
# region always lines up with where the simulation is actually forming coal fastest -- the
# same "one shared cutoff" precedent ICE_TEMP_C above already sets.
CARBONIFEROUS_MIN_TEMP_C = TROPICAL_TEMP_C
CARBONIFEROUS_MIN_PRECIP_MM = HUMID_MM

# Intertidal Zone: a pure elevation-band split of the ocean category (real tides aren't
# modeled) -- ocean water shallow enough to read as "coastal shallows" rather than open ocean.
INTERTIDAL_MAX_DEPTH_M = 30.0

BIOME_NAMES = [
    "Ocean",
    "Ice",
    "Tundra",
    "Boreal Forest",
    "Temperate Desert",
    "Temperate Grassland",
    "Woodland/Shrubland",
    "Temperate Seasonal Forest",
    "Temperate Rainforest",
    "Subtropical Desert",
    "Savanna",
    "Tropical Seasonal Forest",
    "Tropical Rainforest",
    # Appended, not interleaved -- keeps every existing index stable for any external
    # reference to BIOME_NAMES/BIOME_COLORS by position.
    "Wetland",
    "Carboniferous Forest",
    "Intertidal Zone",
]

(
    OCEAN,
    ICE,
    TUNDRA,
    BOREAL_FOREST,
    TEMPERATE_DESERT,
    TEMPERATE_GRASSLAND,
    WOODLAND_SHRUBLAND,
    TEMPERATE_SEASONAL_FOREST,
    TEMPERATE_RAINFOREST,
    SUBTROPICAL_DESERT,
    SAVANNA,
    TROPICAL_SEASONAL_FOREST,
    TROPICAL_RAINFOREST,
    WETLAND,
    CARBONIFEROUS_FOREST,
    INTERTIDAL,
) = range(len(BIOME_NAMES))

# Index-aligned with BIOME_NAMES. Ocean reuses render_image.CLIMATE_OCEAN_BACKDROP_RGB's own
# value directly (kept as a literal here, not an import, to avoid a render_image -> biomes ->
# render_image cycle -- render_image.py is the one that imports this module, not the reverse)
# so a biome map's ocean reads the same as every other climate view's own ocean backdrop.
# Ice is close to render_image.GLACIER_COLOR_RGB for the same "read consistently across
# views" reason. Everything else is a hand-picked, ecologically-suggestive color (browns/tans
# toward dry biomes, greens deepening with both warmth and moisture toward rainforest).
BIOME_COLORS = np.array(
    [
        [18, 28, 55],  # Ocean
        [223, 235, 240],  # Ice
        [156, 171, 158],  # Tundra
        [61, 96, 74],  # Boreal Forest
        [176, 152, 116],  # Temperate Desert
        [168, 178, 107],  # Temperate Grassland
        [126, 143, 90],  # Woodland/Shrubland
        [79, 121, 66],  # Temperate Seasonal Forest
        [42, 94, 68],  # Temperate Rainforest
        [214, 178, 115],  # Subtropical Desert
        [196, 178, 92],  # Savanna
        [58, 122, 66],  # Tropical Seasonal Forest
        [26, 84, 46],  # Tropical Rainforest
        [101, 111, 66],  # Wetland -- muddy marsh/bog green-brown
        [20, 66, 40],  # Carboniferous Forest -- dense, dark primeval swamp-forest green
        [70, 120, 130],  # Intertidal Zone -- lighter, sandier turquoise than open Ocean
    ],
    dtype=np.uint8,
)


# Peak-to-trough brightness swing biome_relative_shade_factor spreads each biome's flat
# BIOME_COLORS entry across, used by render_image.py's _render_combined_view as the "layer
# tint" its land cells are shaded by before the further hypsometric relief blend near real
# peaks (see that module's RELIEF_BLEND_MAX). The shading is deliberately keyed to each cell's
# elevation *rank among other cells of the same biome*, not an absolute elevation scale --
# many biomes occupy a much narrower absolute elevation band by definition (Wetland is capped
# at WETLAND_MAX_ELEVATION_M, Ice is a temperature band that can sit at sea level or on a
# peak) and would render as one flat, undifferentiated color under a fixed absolute scale
# even though real terrain still varies within them.
#
# The swing is a wide +-25% (a proper shaded-relief look, not the old subtle +-8%): a biome's
# lowest and highest terrain can read as visibly light/dark variants of its color. That's only
# tolerable because _render_combined_view no longer relies on a cell's *color* to say which
# biome it is -- it writes the biome id into the render's alpha channel instead (see
# COMBINED_*_ID_CODE there), so frontend/src/legendData.ts's click-to-highlight reads the id
# directly rather than reverse-matching a small set of shaded RGB variants within a tolerance.
BIOME_SHADE_AMPLITUDE = 0.25
BIOME_SHADE_MIN = 1.0 - BIOME_SHADE_AMPLITUDE
BIOME_SHADE_MAX = 1.0 + BIOME_SHADE_AMPLITUDE


def biome_relative_shade_factor(biome_ids: np.ndarray, elevation_m: np.ndarray) -> np.ndarray:
    """Per-cell brightness multiplier, same shape as `biome_ids`: a *continuous* ramp from
    BIOME_SHADE_MIN at a biome's lowest-elevation cell to BIOME_SHADE_MAX at its highest,
    linear in that cell's elevation *rank among same-biome cells only* (see
    BIOME_SHADE_AMPLITUDE's own docstring for why relative-to-biome rather than an absolute
    elevation scale). Rank-based (an even spread over the sorted order) rather than a
    quantile-of-value ramp, so a biome with one outlier peak and otherwise dead-flat terrain
    still spans the full range smoothly instead of piling almost every cell at one end. A
    continuous ramp (rather than the old 3 discrete tiers) is what keeps within-biome relief
    from reading as visible bands/contour lines. `biome_ids`/`elevation_m` must be the same
    shape; multiply this elementwise into a biome's flat BIOME_COLORS entry to shade it (see
    render_image.py's _render_combined_view)."""
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
    neighbor's real great-circle spacing in meters (longitude narrowed by cos(lat), same
    convention as everywhere else in this codebase -- see plates.iter_local_lattice). Feeds
    classify_wetland's own WETLAND_MAX_SLOPE cutoff -- the same threshold erosion.compute_slope's
    own node-cloud slope (a different, finer discretization) is tuned against; see this module's
    own docstring for why an approximate, visually-tuned cutoff, not fit to any dataset, is this
    codebase's norm. np.roll wraps at the poles too (a minor, visually inconsequential artifact
    right at the map's own poles), the same "not worth special-casing" tradeoff this codebase
    already accepts elsewhere (e.g. the Plate Inspector's antipodal-projection limitation)."""
    grid_h, grid_w = elevation_m.shape
    dlat_km = (np.pi / grid_h) * PLANET_RADIUS_KM
    dlon_km = np.maximum((2 * np.pi / grid_w) * PLANET_RADIUS_KM * np.cos(np.radians(lat_deg))[:, None], 1.0)
    d_ns = np.abs(elevation_m - np.roll(elevation_m, 1, axis=0)) / (dlat_km * 1000.0)
    d_ew = np.abs(elevation_m - np.roll(elevation_m, 1, axis=1)) / (dlon_km * 1000.0)
    return np.maximum(d_ns, d_ew)


def classify_wetland(
    temperature_c: np.ndarray, precipitation_mm: np.ndarray, elevation_m: np.ndarray, slope: np.ndarray, is_ocean: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """(is_wetland, is_carboniferous_forest) boolean arrays, same shape as the inputs -- shared
    by classify_biomes (the map view, below) and geology.py's own per-node coal-formation
    accumulation, so the map's Carboniferous Forest region always lines up with where the
    simulation is actually forming coal fastest (see WETLAND_MAX_ELEVATION_M's own comment for
    the ICE_TEMP_C precedent this follows). Both require flat, low-lying land -- a floodplain,
    delta, or coastal marsh, not an upland bog perched on a slope; Carboniferous Forest is the
    warm (>= CARBONIFEROUS_MIN_TEMP_C), very wet (>= CARBONIFEROUS_MIN_PRECIP_MM) subtype,
    Wetland the cooler/drier remainder (>= ICE_TEMP_C so it excludes frozen ground, >=
    WETLAND_MIN_PRECIP_MM)."""
    elevation_m = np.asarray(elevation_m)
    slope = np.asarray(slope)
    temp = np.asarray(temperature_c)
    precip = np.asarray(precipitation_mm)
    is_ocean = np.asarray(is_ocean)

    flat_low = ~is_ocean & (elevation_m > 0) & (elevation_m <= WETLAND_MAX_ELEVATION_M) & (slope <= WETLAND_MAX_SLOPE)
    is_carboniferous = flat_low & (temp >= CARBONIFEROUS_MIN_TEMP_C) & (precip >= CARBONIFEROUS_MIN_PRECIP_MM)
    is_wetland = flat_low & (temp >= ICE_TEMP_C) & (precip >= WETLAND_MIN_PRECIP_MM) & ~is_carboniferous
    return is_wetland, is_carboniferous


def classify_biomes(
    temperature_c: np.ndarray,
    precipitation_mm: np.ndarray,
    elevation_m: np.ndarray,
    slope: np.ndarray,
    is_ocean: np.ndarray,
    sea_level_m: float = 0.0,
) -> np.ndarray:
    """Per-cell biome id (index into BIOME_NAMES/BIOME_COLORS) from mean temperature and
    annual precipitation, plus elevation/slope for Wetland/Carboniferous Forest/Intertidal
    Zone (see classify_wetland and INTERTIDAL_MAX_DEPTH_M) -- the core temp/precip axes are
    the same two the real Whittaker diagram uses, in the same broad relative order
    (cold-to-hot, dry-to-wet), though with this module's own boundary values (see module
    docstring for why). `is_ocean` still wins over temperature/precipitation for the ordinary
    land bands -- temperature/precipitation/wetland-detection are land-surface concepts, and
    is_ocean already settles the question for a water cell -- Intertidal Zone is purely a
    further elevation-band split *within* the ocean category, not an exception to that rule.
    `temperature_c`/`precipitation_mm`/`elevation_m`/`slope`/`is_ocean` must all be the same
    shape; the result is that same shape. Uses np.select (first-matching condition wins)
    rather than chained np.where overwrites, so each band's cutoffs stay a self-contained,
    independently checkable list instead of relying on write order to get right-of-boundary
    cells correct."""
    temp = np.asarray(temperature_c)
    precip = np.asarray(precipitation_mm)
    elevation_m = np.asarray(elevation_m)
    slope = np.asarray(slope)
    is_ocean = np.asarray(is_ocean)

    tropical = temp >= TROPICAL_TEMP_C
    temperate = (temp >= COLD_TEMP_C) & (temp < TROPICAL_TEMP_C)
    cold = (temp >= ICE_TEMP_C) & (temp < COLD_TEMP_C)
    ice = temp < ICE_TEMP_C
    is_intertidal = is_ocean & ((sea_level_m - elevation_m) <= INTERTIDAL_MAX_DEPTH_M)
    is_wetland, is_carboniferous = classify_wetland(temp, precip, elevation_m, slope, is_ocean)

    condlist = [
        is_ocean & is_intertidal,
        is_ocean,
        ice,
        is_carboniferous,
        is_wetland,
        cold & (precip < SEMI_ARID_MM),
        cold,
        temperate & (precip < ARID_MM),
        temperate & (precip < SEMI_ARID_MM),
        temperate & (precip < SUB_HUMID_MM),
        temperate & (precip < HUMID_MM),
        temperate,
        tropical & (precip < SEMI_ARID_MM),
        tropical & (precip < SUB_HUMID_MM),
        tropical & (precip < HUMID_MM),
        tropical,
    ]
    choicelist = [
        INTERTIDAL,
        OCEAN,
        ICE,
        CARBONIFEROUS_FOREST,
        WETLAND,
        TUNDRA,
        BOREAL_FOREST,
        TEMPERATE_DESERT,
        TEMPERATE_GRASSLAND,
        WOODLAND_SHRUBLAND,
        TEMPERATE_SEASONAL_FOREST,
        TEMPERATE_RAINFOREST,
        SUBTROPICAL_DESERT,
        SAVANNA,
        TROPICAL_SEASONAL_FOREST,
        TROPICAL_RAINFOREST,
    ]
    return np.select(condlist, choicelist, default=OCEAN).astype(np.int64)
