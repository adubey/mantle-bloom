"""Whittaker-diagram-inspired biome classification: buckets a (temperature, precipitation)
pair into one of a fixed set of named terrestrial biomes, purely from climate.py's own
already-computed grid fields -- no new simulation state, nothing cached anywhere, since (like
render_image.py's own temperature_colors/humidity_colors) this is a stateless function of
values climate.py already produces every step. The boundary values below are a
simplification, same spirit as this codebase's other openly-approximate constants (e.g.
erosion.py's own RAIN_EROSION_COEFFICIENT, "tuned by rough order-of-magnitude reasoning" not
lifted from a literal scientific reference) -- picked for a visually sensible spread of named
regions on a generated world, not fit against any specific real-world dataset.
"""

from __future__ import annotations

import numpy as np

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
    ],
    dtype=np.uint8,
)


def classify_biomes(temperature_c: np.ndarray, precipitation_mm: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """Per-cell biome id (index into BIOME_NAMES/BIOME_COLORS) from mean temperature and
    annual precipitation alone -- the same two axes the real Whittaker diagram uses, in the
    same broad relative order (cold-to-hot, dry-to-wet), though with this module's own
    boundary values (see module docstring for why). `is_ocean` wins over everything else --
    temperature/precipitation are land-surface climate concepts, and is_ocean already settles
    the question for a water cell. All three arrays must be the same shape; the result is
    that same shape. Uses np.select (first-matching condition wins) rather than chained
    np.where overwrites, so each band's cutoffs stay a self-contained, independently
    checkable list instead of relying on write order to get right-of-boundary cells correct."""
    temp = np.asarray(temperature_c)
    precip = np.asarray(precipitation_mm)
    is_ocean = np.asarray(is_ocean)

    tropical = temp >= TROPICAL_TEMP_C
    temperate = (temp >= COLD_TEMP_C) & (temp < TROPICAL_TEMP_C)
    cold = (temp >= ICE_TEMP_C) & (temp < COLD_TEMP_C)
    ice = temp < ICE_TEMP_C

    condlist = [
        is_ocean,
        ice,
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
    ]
    return np.select(condlist, choicelist, default=OCEAN).astype(np.int64)
