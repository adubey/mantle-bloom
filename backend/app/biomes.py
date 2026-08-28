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
from scipy.ndimage import uniform_filter1d

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


# ---------------------------------------------------------------------------------------
# Boundary cleanup: smooth_biome_field
# ---------------------------------------------------------------------------------------
# classify_biomes buckets continuous (temp, precip, elevation, slope) into hard-edged bands.
# Wherever the underlying field runs nearly parallel to a band cutoff -- or the nearest-node
# resample of elevation feeding the Wetland/Carboniferous/Intertidal tests carries
# stair-stepping -- adjacent render cells flip back and forth across the cutoff, and the
# Biome view (un-blurred by design, unlike Combined) reads as a dithered smear along every
# biome interface plus scattered single-cell speckle. `smooth_biome_field` is a *stateless*
# cleanup pass the render/stats callers apply on top of classify_biomes; classify_biomes
# itself stays a pure, any-shape band lookup (its unit tests and geology.py's own
# classify_wetland call both depend on that). Two mechanisms, both keyed to "this cell is
# only marginally its biome":
#   1. the Wetland/Carboniferous/Intertidal elevation tests (the `<= WETLAND_MAX_ELEVATION_M`
#      band and the Intertidal depth split) run against a lightly box-blurred copy of
#      `elevation` -- never the elevation used for relief shading -- so a nearest-node
#      stair-step near a coast or a 50 m contour doesn't dither. `slope` is left raw: it
#      already gates hard (WETLAND_MAX_SLOPE is 1 m/km), and blurring a noisy positive
#      magnitude only biases it upward and erases these biomes wholesale;
#   2. a cell that is either within a small margin of the nearest temp/precip band edge --
#      the band the coherent climate noise alone (see climate.HUMIDITY_NOISE_STD / the
#      MFC_NOISE_* block) can shove it across -- or already surrounded by mostly other biomes
#      (an edge or speckle cell, regardless of which cutoff put it there) can be out-voted by
#      a strong majority of its 8 neighbours. A cell comfortably inside a biome, with mostly
#      same-biome neighbours, is never touched -- so genuine regions keep their shape (bar a
#      one-cell edge shave) and only outnumbered speckle / single ragged steps straighten.
#
# All picked (like every other cutoff in this module) for a visually sensible result on a
# generated world, not fit against a dataset.

# A cell within this of the nearest temperature/precipitation band edge is "on the cusp" and
# eligible to be out-voted; further than this it is solidly its biome and is left alone.
BIOME_CUSP_MARGIN_TEMP_C = 1.5
BIOME_CUSP_MARGIN_PRECIP_MM = 150.0

# Adopt the most common non-ocean neighbour biome only if at least this fraction of a cell's
# valid (non-ocean) neighbours agree on it -- 0.65 means 6+ of 8. Kept well above 0.5 so a
# genuine two-biome interface (roughly half and half) is left in place.
BIOME_VOTE_MIN_NEIGHBOUR_FRACTION = 0.65
BIOME_VOTE_ITERATIONS = 2

# A land cell that already disagrees with at least this many of its 8 non-ocean neighbours is
# treated as a cusp cell too, regardless of the temp/precip margin -- catches speckle driven
# by the elevation/slope resample rather than a climate band. Deliberately high (a genuinely
# outnumbered cell, not just a boundary cell -- a clean straight interface has ~3/8
# disagreeing, a diagonal ~2) so region edges and necks are left intact.
BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS = 6

# Real-world radius of the box blur applied to the `elevation` copy used by the
# Wetland/Carboniferous/Intertidal elevation/depth predicates -- about one default node
# spacing (elevation_lines.TARGET_LINE_SPACING_KM / DEFAULT_NODE_DENSITY ~= 62 km), enough to
# bridge the nearest-node stair-stepping the resample introduces without moving a real
# coastline or floodplain edge appreciably. Expressed in km, converted to a cell count from
# the grid's own resolution, so it covers the same distance regardless of World.climate_density.
ELEVATION_BAND_SMOOTH_RADIUS_KM = 75.0


def _smooth_for_bands(field: np.ndarray, cell_km: float) -> np.ndarray:
    """Separable box blur at ELEVATION_BAND_SMOOTH_RADIUS_KM, longitude-wrapping (axis 1) and
    edge-replicating at the poles (axis 0). One uniform_filter1d per axis -- O(n) in the
    window size, unlike an iterated Jacobi sweep."""
    size = max(1, int(round(2.0 * ELEVATION_BAND_SMOOTH_RADIUS_KM / cell_km)))
    out = uniform_filter1d(np.asarray(field, dtype=float), size, axis=1, mode="wrap")
    return uniform_filter1d(out, size, axis=0, mode="nearest")


_NEIGHBOUR_OFFSETS = tuple((dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0))


def _band_cusp_mask(temperature_c: np.ndarray, precipitation_mm: np.ndarray) -> np.ndarray:
    """Boolean, same shape as the inputs: True where the cell sits within
    BIOME_CUSP_MARGIN_TEMP_C / _PRECIP_MM of the nearest temperature or precipitation band
    edge classify_biomes splits on."""
    temp_edges = np.array([ICE_TEMP_C, COLD_TEMP_C, TROPICAL_TEMP_C])
    precip_edges = np.array([ARID_MM, SEMI_ARID_MM, SUB_HUMID_MM, HUMID_MM])
    temp = np.asarray(temperature_c)[..., None]
    precip = np.asarray(precipitation_mm)[..., None]
    near_temp = np.min(np.abs(temp - temp_edges), axis=-1) <= BIOME_CUSP_MARGIN_TEMP_C
    near_precip = np.min(np.abs(precip - precip_edges), axis=-1) <= BIOME_CUSP_MARGIN_PRECIP_MM
    return near_temp | near_precip


def _disagreeing_neighbour_count(biome_ids: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """Per-cell count of the 8 neighbours that are land and a different biome."""
    count = np.zeros(biome_ids.shape, dtype=np.int16)
    for dy, dx in _NEIGHBOUR_OFFSETS:
        nb = np.roll(np.roll(biome_ids, dy, axis=0), dx, axis=1)
        nb_ocean = np.roll(np.roll(is_ocean, dy, axis=0), dx, axis=1)
        count += ((nb != biome_ids) & ~nb_ocean).astype(np.int16)
    return count


def _neighbour_vote(biome_ids: np.ndarray, eligible: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """One majority-vote pass: each `eligible` land cell adopts the most common land biome
    among its 8 neighbours when at least BIOME_VOTE_MIN_NEIGHBOUR_FRACTION of its valid
    (non-ocean) neighbours agree on it and it differs from the cell's current biome. The
    modal neighbour is found by sorting the 8 neighbour ids per cell and taking the longest
    equal run -- O(1) array ops in the biome count, not a per-biome loop."""
    stack = np.stack([np.roll(np.roll(biome_ids, dy, axis=0), dx, axis=1) for dy, dx in _NEIGHBOUR_OFFSETS])
    ocean_stack = np.stack([np.roll(np.roll(is_ocean, dy, axis=0), dx, axis=1) for dy, dx in _NEIGHBOUR_OFFSETS])
    # Ocean neighbours -> -1 so they sort to the front and can never form the winning run.
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
) -> np.ndarray:
    """classify_biomes plus the stateless boundary-cleanup pass described in the block above
    -- this is the classification the Biome/Combined map views and /world/stats use, while
    `classify_biomes` stays the raw per-cell primitive. Inputs must all be the same 2D
    (H, W) grid (the vote needs spatial neighbours); the result is that same shape."""
    temperature_c = np.asarray(temperature_c)
    precipitation_mm = np.asarray(precipitation_mm)
    elevation_m = np.asarray(elevation_m)
    slope = np.asarray(slope)
    is_ocean = np.asarray(is_ocean)
    if elevation_m.ndim != 2:
        raise ValueError("smooth_biome_field needs 2D grid inputs; use classify_biomes for a flat array")

    equator_cell_km = 2.0 * np.pi * PLANET_RADIUS_KM / elevation_m.shape[1]
    elev_smoothed = _smooth_for_bands(elevation_m, equator_cell_km)

    # Substituting the blurred elevation for the elevation-band / Intertidal-depth tests is
    # safe: inside classify_biomes elevation feeds only those -- every ordinary
    # temperature/precipitation band ignores it. `slope` stays raw (see mechanism 1 above).
    biome_ids = classify_biomes(temperature_c, precipitation_mm, elev_smoothed, slope, is_ocean, sea_level_m)

    # Eligibility is fixed from the initial classification, not recomputed each iteration --
    # the band-cusp part can't change (temp/precip don't move) and freezing the ragged part
    # keeps a second pass from chasing its own tail across a wide boundary.
    land = ~is_ocean
    eligible = land & (
        _band_cusp_mask(temperature_c, precipitation_mm)
        | (_disagreeing_neighbour_count(biome_ids, is_ocean) >= BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS)
    )
    for _ in range(BIOME_VOTE_ITERATIONS):
        biome_ids = _neighbour_vote(biome_ids, eligible, is_ocean)
    return biome_ids
