"""Resource deposits (coal, oil & gas, mineral-rich ore) and soil formation -- all persistent
per-node fields (see plates.ElevationLine) that accumulate over sustained geologic conditions,
the same "rides along for free with plate rotation" pattern channel_depth/silt_depth/
glacier_depth already use.

**Coal** forms from peat accumulating in waterlogged, low-relief wetland/swamp terrain (see
biomes.classify_wetland) -- warm tropical swamp forest (Carboniferous Forest) accumulates peat
fastest, real Carboniferous/Permian coal being predominantly of exactly that tropical
lowland-swamp origin, but any sufficiently flat, low, wet land accumulates some (a cooler bog
or marsh, plain "Wetland"). Once formed, coal_deposit_m never decreases -- buried peat isn't
un-buried by a later climate shift, the same self-reinforcing, monotonic convention
silt_depth/channel_depth already establish.

**Oil and gas** form the same way, in shallow marine shelf water instead of on land --
organic-rich (mostly planktonic) sediment settling in a quiet, sunlit, nutrient-fed shelf sea,
real petroleum source rocks. Nutrient input is boosted near a river mouth: reuses
hydrology.py's own `water_deposited` directly -- an ocean node a river empties into already
tells this module "land-derived nutrient runoff reaches the sea here," no separate river-mouth
search needed. Also monotonic, gated to shelf water via the same shelf-distance concept
bathymetry.py already established (`bathymetry.SHELF_RANGE_RAD`).

**Minerals** come from volcanism -- grown directly inside volcanism.apply_volcanic_activity's
own eruption roll (see that module), not here. This module only *reads* mineral_deposit_m back,
for the soil mineral-content term below; it never writes it.

**Soil** forms from weathered rock (a fraction of erosion.py's own weathering term, already
computed this step, becomes fresh regolith in place) plus organic matter contributed by
"pioneering species and plant roots" (approximated as a warmth/moisture-driven productivity
term, since this codebase has no explicit vegetation field), is stripped by fast erosion
(rain/river erosion carries topsoil away), and gets an extra deposit wherever a slow, big river
drops its sediment load (erosion.py's own floodplain deposition term, reused directly rather
than re-derived -- soil "carried down river ... and deposited in flood plains" is literally the
same mechanism). Mineral content relaxes toward a target driven by the node's own accumulated
mineral_deposit_m (weathered/hydrothermal rock feeding the soil above it) plus a small baseline
from ordinary rock weathering; organic content relaxes toward the same productivity term soil
depth's own organic contribution uses. Both land-only, [0, 1]-clamped, and -- unlike
coal/oil-gas/minerals -- can fall as well as rise (real soil erodes away, unlike buried ore).
The richest soil needs *both* high mineral and high organic content at once: the Soil Quality
map view (render_image.py) scores fertility as `sqrt(mineral * organic)`, which rewards having
both far more than either alone.

Runs every step, right after erosion.py/bathymetry.py/volcanism.py in world.step_world, reusing
World.hydrology_cache (populated by erosion.apply_erosion, called first) and the ErosionResult
that same call returns -- rather than recomputing climate or flow routing a second time. Its own
node ordering comes directly from `World.hydrology_cache.line_refs`, the exact (plate,
line_index, start, end) list that produced hydro's own flat arrays, so there's no separate
gather to keep in sync with it.
"""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import biomes
from .bathymetry import SHELF_RANGE_RAD
from .noise import SphereNoise
from .plates import Plate

if TYPE_CHECKING:
    from .erosion import ErosionResult
    from .world import World

# Coal: Carboniferous Forest (warm tropical swamp) accumulates several times faster than plain
# Wetland (a cooler bog/marsh) -- real Carboniferous/Permian coal is predominantly of that
# tropical-swamp origin. Starting points, not derived from any dataset -- same "rough
# order-of-magnitude reasoning" spirit as this codebase's other from-scratch rates (e.g.
# erosion.py's own RAIN_EROSION_COEFFICIENT), tuned so a sustained multi-hundred-Myr wetland
# produces a visually meaningful deposit without saturating within a handful of steps.
COAL_BASE_RATE_M_PER_MYR = 0.05
COAL_CARBONIFEROUS_RATE_M_PER_MYR = 0.4
MAX_COAL_DEPOSIT_M = 50.0

# Oil & gas: a baseline shelf-sea productivity rate, boosted near a river mouth (real deltas --
# the Gulf of Mexico, Niger Delta -- are major petroleum provinces for exactly this reason).
# OIL_GAS_RIVER_REFERENCE_M is the water_deposited scale (hydrology.py's own stylized units, see
# that module) at which the river-mouth boost saturates.
OIL_GAS_BASE_RATE_M_PER_MYR = 0.02
OIL_GAS_RIVER_BOOST_RATE_M_PER_MYR = 0.3
OIL_GAS_RIVER_REFERENCE_M = 5.0
MAX_OIL_GAS_DEPOSIT_M = 30.0

# Soil formation/loss. WEATHERING_TO_REGOLITH is the fraction of this step's weathering erosion
# (erosion.py's own term) that becomes fresh soil in place rather than bare-rock loss/dust;
# ORGANIC_RATE is organic matter's own peak contribution to soil depth, scaled by `productivity`
# below; STRIP_FRACTION is the fraction of rain+river erosion that also carries away existing
# topsoil (a real, distinct effect from weathering forming new soil); DEPOSITION_FRACTION reuses
# erosion.py's own floodplain-deposit amount directly, at its own (smaller) share, since not all
# deposited sediment is soil-grade.
SOIL_WEATHERING_TO_REGOLITH = 0.4
SOIL_ORGANIC_RATE_M_PER_MYR = 0.15
SOIL_STRIP_FRACTION = 0.5
SOIL_DEPOSITION_FRACTION = 0.3
MAX_SOIL_DEPTH_M = 15.0

# Organic/mineral content relax toward a target (same exponential-toward-target style
# bathymetry.py's own shelf relaxation already uses) rather than integrate directly -- a
# saturating [0, 1] fraction, not an unbounded accumulating depth.
SOIL_ORGANIC_RELAX_RATE_PER_MYR = 0.3
SOIL_MINERAL_RELAX_RATE_PER_MYR = 0.2
# Ordinary rock weathering alone still releases some minerals even with no nearby volcanism;
# MINERAL_CONTENT_REFERENCE_M is the mineral_deposit_m scale (volcanism.py's own units) at
# which the volcanic-enrichment boost on top of that baseline saturates.
MINERAL_CONTENT_BASELINE = 0.2
MINERAL_CONTENT_REFERENCE_M = 5.0

# Generation-time "initial soil maturity" slider (see world.generate_world) -- a mature
# (maturity=1.0) world's average starting soil depth, and how much per-node noise texture (see
# SphereNoise) varies it, the same "texture, not a flat value" treatment plates.py's own initial
# elevation noise already uses.
INITIAL_SOIL_MAX_DEPTH_M = 8.0
INITIAL_SOIL_NOISE_AMPLITUDE = 0.3


def apply_resource_formation(world: "World", years: float, erosion_result: "ErosionResult") -> None:
    """Grows/relaxes soil_depth/soil_mineral_content/soil_organic_content and accumulates
    coal_deposit_m/oil_gas_deposit_m for every plate's elevation nodes, from this step's already
    -computed erosion terms (`erosion_result`) and flow routing (`world.hydrology_cache`) --
    see module docstring for the formulas and reasoning. Mutates world.plates' line soil/
    resource fields in place; never touches node positions, line topology, or elevation, so
    this can't interact with line regularization, point reassignment, or erosion/bathymetry at
    all. A no-op if `world.hydrology_cache` hasn't been populated yet (mirrors erosion.py's own
    empty-world early return -- see World.hydrology_cache)."""
    hydro = world.hydrology_cache
    if hydro is None:
        return
    line_refs = hydro.line_refs
    n = len(erosion_result.points)
    if n == 0:
        return

    prior_soil_depth = np.zeros(n)
    prior_soil_mineral_content = np.zeros(n)
    prior_soil_organic_content = np.zeros(n)
    prior_coal_deposit_m = np.zeros(n)
    prior_oil_gas_deposit_m = np.zeros(n)
    prior_mineral_deposit_m = np.zeros(n)
    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        prior_soil_depth[start:end] = line.soil_depth
        prior_soil_mineral_content[start:end] = line.soil_mineral_content
        prior_soil_organic_content[start:end] = line.soil_organic_content
        prior_coal_deposit_m[start:end] = line.coal_deposit_m
        prior_oil_gas_deposit_m[start:end] = line.oil_gas_deposit_m
        prior_mineral_deposit_m[start:end] = line.mineral_deposit_m

    is_ocean = hydro.is_ocean
    is_land = ~is_ocean
    dt_myr = years / 1_000_000.0

    is_wetland, is_carboniferous = biomes.classify_wetland(
        erosion_result.temperature_c, erosion_result.precipitation_mm, erosion_result.elevation, erosion_result.slope, is_ocean
    )

    # --- Coal (land, monotonic) ---
    coal_rate = np.where(is_carboniferous, COAL_CARBONIFEROUS_RATE_M_PER_MYR, np.where(is_wetland, COAL_BASE_RATE_M_PER_MYR, 0.0))
    new_coal_deposit_m = np.clip(prior_coal_deposit_m + coal_rate * dt_myr, 0.0, MAX_COAL_DEPOSIT_M)

    # --- Oil & gas (ocean, monotonic), gated to shelf water and boosted near a river mouth ---
    if np.any(is_land):
        land_tree = cKDTree(erosion_result.points[is_land])
        dist_to_land, _ = land_tree.query(erosion_result.points)
    else:
        dist_to_land = np.full(n, np.inf)
    is_shelf = is_ocean & (dist_to_land <= SHELF_RANGE_RAD)
    river_boost = np.clip(hydro.water_deposited / OIL_GAS_RIVER_REFERENCE_M, 0.0, 1.0)
    oil_gas_rate = np.where(is_shelf, OIL_GAS_BASE_RATE_M_PER_MYR + OIL_GAS_RIVER_BOOST_RATE_M_PER_MYR * river_boost, 0.0)
    new_oil_gas_deposit_m = np.clip(prior_oil_gas_deposit_m + oil_gas_rate * dt_myr, 0.0, MAX_OIL_GAS_DEPOSIT_M)

    # --- Soil ---
    # Warmer + wetter land supports more biomass -- 0 at/below freezing or bone dry, saturating
    # at HUMID_MM precipitation and TROPICAL_TEMP_C warmth. Also drives organic_content's own
    # relaxation target below, so the two stay consistent with each other.
    productivity = np.clip(erosion_result.precipitation_mm / biomes.HUMID_MM, 0.0, 1.0) * np.clip(
        erosion_result.temperature_c / biomes.TROPICAL_TEMP_C, 0.0, 1.0
    )
    soil_gain = SOIL_WEATHERING_TO_REGOLITH * erosion_result.weathering + SOIL_ORGANIC_RATE_M_PER_MYR * productivity * dt_myr
    deposition_gain = SOIL_DEPOSITION_FRACTION * erosion_result.sediment_deposited
    soil_loss = SOIL_STRIP_FRACTION * (erosion_result.rain + erosion_result.river)
    new_soil_depth = np.clip(prior_soil_depth + soil_gain + deposition_gain - soil_loss, 0.0, MAX_SOIL_DEPTH_M)
    new_soil_depth = np.where(is_ocean, 0.0, new_soil_depth)

    organic_relax = 1.0 - np.exp(-SOIL_ORGANIC_RELAX_RATE_PER_MYR * dt_myr)
    new_soil_organic_content = prior_soil_organic_content + (productivity - prior_soil_organic_content) * organic_relax

    mineral_target = np.clip(MINERAL_CONTENT_BASELINE + (1.0 - MINERAL_CONTENT_BASELINE) * (prior_mineral_deposit_m / MINERAL_CONTENT_REFERENCE_M), 0.0, 1.0)
    mineral_relax = 1.0 - np.exp(-SOIL_MINERAL_RELAX_RATE_PER_MYR * dt_myr)
    new_soil_mineral_content = prior_soil_mineral_content + (mineral_target - prior_soil_mineral_content) * mineral_relax

    has_soil = new_soil_depth > 0.0
    new_soil_organic_content = np.clip(np.where(is_land & has_soil, new_soil_organic_content, 0.0), 0.0, 1.0)
    new_soil_mineral_content = np.clip(np.where(is_land & has_soil, new_soil_mineral_content, 0.0), 0.0, 1.0)

    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        plate.lines[line_index] = dataclasses.replace(
            line,
            soil_depth=new_soil_depth[start:end],
            soil_mineral_content=new_soil_mineral_content[start:end],
            soil_organic_content=new_soil_organic_content[start:end],
            coal_deposit_m=new_coal_deposit_m[start:end],
            oil_gas_deposit_m=new_oil_gas_deposit_m[start:end],
        )


def seed_initial_soil(plate_list: list[Plate], seed: int, initial_soil_maturity: float) -> None:
    """Seeds soil_depth/soil_organic_content/soil_mineral_content on every land node (elevation
    > 0), scaled by `initial_soil_maturity` (the UI's "initial soil maturity" slider, 0 to 1) --
    called once from world.generate_world, right after plates.generate_plates. At 0 (the
    default), every land node starts at exactly zero soil -- ElevationLine's own zero defaults
    already give this, so this function is a no-op rather than special-casing it. Deliberately
    *not* climate-informed at seed time (no biome differentiation yet at generation, unlike the
    organic-content relaxation apply_resource_formation drives every step thereafter) -- just a
    uniform-ish starting maturity with a noise texture for visual variation, the same role
    plates.py's own initial elevation noise plays. Uses a fresh RNG stream (`seed + 2`) so it
    doesn't disturb plates.generate_plates' own stream or world.py's `seed + 1` mantle-center
    stream."""
    if initial_soil_maturity <= 0.0:
        return
    maturity = min(initial_soil_maturity, 1.0)
    rng = np.random.default_rng(seed + 2)
    noise = SphereNoise(rng, octaves=3, base_freq=3.0)
    for plate in plate_list:
        for line_index, line in enumerate(plate.lines):
            if len(line.theta) == 0:
                continue
            world_pts = line.world_xyz(plate.frame)
            is_land = line.elevation > 0.0
            texture = np.clip(1.0 + INITIAL_SOIL_NOISE_AMPLITUDE * noise.sample(world_pts), 0.0, None)
            depth = np.where(is_land, maturity * INITIAL_SOIL_MAX_DEPTH_M * texture, 0.0)
            fraction = np.where(is_land, np.clip(maturity * texture, 0.0, 1.0), 0.0)
            plate.lines[line_index] = dataclasses.replace(
                line,
                soil_depth=depth,
                soil_organic_content=fraction,
                soil_mineral_content=fraction,
            )
