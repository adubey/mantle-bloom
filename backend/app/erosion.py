"""Weather-driven erosion: rain/sheet erosion, river-channelized erosion, and weathering,
all directly reducing elevation every step; slow/big rivers deposit part of what they carry
downstream instead of losing it all to the coast. The other direction of the weather<->
geology coupling -- terrain influencing weather (lapse-rate cooling, mountain wind
deflection, orographic rain shadow) -- already exists in climate.py, as part of its climate
pipeline. This module is the new half: climate.py's fixed grid, and hydrology.py's flow
routing over the geology node cloud, feeding back onto elevation.

**The mapping problem, and why it's easier in this direction.** climate.py already solves
node-cloud -> grid (`_sample_elevation_and_crust`'s cKDTree nearest-neighbor resample) to
build the grid in the first place. This module needs the reverse, grid -> node-cloud, which
turns out to be simpler: the climate grid is a plain regular lat/lon lattice, so a node's
own world position converts straight to (lat, lon) (geometry.xyz_to_latlon) and then to a
grid (row, col) by direct arithmetic -- no tree, no resampling, just array indexing (see
`climate_grid_indices`, whose row/col convention mirrors climate._build_grid exactly -- public,
not private, since geology.py also needs it, see that module).

**Slope is the one genuinely new piece of math climate.py's grid can't hand over for free**
(it gets slope from neighbor-index differences; an irregular node cloud has no such
structure) -- this reuses reassign.py's whole-world cKDTree pattern (build once, query k
nearest neighbors) instead: for each node, the elevation drop to the *lowest* of its nearest
neighbors (the "slope to lowest neighbor" definition used throughout this module), divided
by the real great-circle distance to that neighbor. This is a genuine dimensionless rise/run
-- not elevation-drop-per-grid-step, which would conflate slope with grid resolution -- so
RAIN_EROSION_COEFFICIENT below is tuned for this rise/run scale (order 0.001-0.1), not for
meters-of-drop-per-grid-cell (order 10-100s of meters); a coefficient tuned against that
coarser, resolution-dependent scale would make rain erosion negligible here.
WEATHERING_COEFFICIENT needs no such rescaling, since it depends only on wind speed, whose
scale (MERIDIONAL_BASE_SPEED = 6.0) is fixed by this module's own climate pipeline rather
than by any grid convention. River erosion needs the same reasoning as rain erosion, one
level further removed: it depends on `water_accum` (see hydrology.py), which is itself
downstream-accumulated *precipitation* (not a raw grid-cell count, which would implicitly
scale with resolution) -- RIVER_EROSION_COEFFICIENT is derived the same order-of-magnitude
way RAIN_EROSION_COEFFICIENT is.

**Erosion sources implemented here, and what's still out of scope.** Submarine erosion (deep
bottom currents plus pressure-driven mass wasting slumping a freshly-uplifted submerged range
back down -- see SUBMARINE_EROSION_* below) and coastal erosion (wave attack plus freeze-thaw
frost shattering wherever the climate cycles through freezing -- see COASTAL_EROSION_* below)
were both previously out of scope ("a distinct source never implemented here") and are now
implemented: they are the sea floor's and the shoreline's counterparts to the subaerial
sources, and their eroded rock sheds onto the surrounding sea floor as marine sediment
(_spread_marine_sediment) rather than into any river's flow graph. Together they are also what
keeps a range built by two *submerged* plates colliding growing far more slowly than a
subaerial one -- the sea floor has denudation the dry interior doesn't. Weathering's vegetation
boost is out of scope (no vegetation field). River-channelized erosion, downstream deposition,
and glacier erosion -- previously out of scope for the same reason ("needs flow routing over a
rotating, irregular per-plate lattice, a separate, harder problem") -- are now implemented, see
hydrology.py for how that flow-routing graph is built and how glacier_depth itself is
grown/melted/flowed. Glacier-driven **flattening** (broad terrain smoothing under an ice
sheet, distinct from the directional erosion term) and **seismic erosion** (earthquake-
triggered landsliding, scaled by elevation as a stand-in for how tectonically active a range
is -- this model has no separate fault/stress field) are both mantle-bloom-original additions
-- see `_flatten` and SEISMIC_EROSION_* below. Glacially-eroded material that isn't deposited
as immediate subglacial till also travels along the glacier's own real flow path rather than
water's, dumping out only once it clears the ice -- a terminal moraine/outwash deposit, see
the comment above GLACIER_TILL_FRACTION.

**No lag.** This module always erodes against the *current* step's climate and flow routing,
never a stale previous-step snapshot: climate.py here is fully stateless and cheap enough to
call fresh every step. Flow routing, though, is comparatively expensive (no JIT), so rather
than computing it twice, this module computes it once (via hydrology.compute_hydrology) and
reuses the result for both erosion and the world's cached river/lake fields
(World.hydrology_cache).
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import biomes, climate, continental_ledger, cratons, faults, geometry, hydrology, lithosphere, mobile_cover
from .elevation_lines import (
    ELEV_CHANGE_COASTAL_LEVELING,
    ELEV_CHANGE_COLLISION,
    ELEV_CHANGE_DEPOSITION,
    ELEV_CHANGE_EROSION,
    ELEV_CHANGE_FAULT_NORMAL,
    ELEV_CHANGE_FAULT_STRIKE_SLIP,
    ELEV_CHANGE_GLACIAL_FLATTEN,
    ELEV_CHANGE_LAKE_SILT,
    ELEV_CHANGE_LATERAL_MAGMA,
    ELEV_CHANGE_MARINE,
    ELEV_CHANGE_MIN_DELTA_M,
    ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR,
    ELEV_CHANGE_VOLCANO,
    MAX_ELEVATION_M,
    MIN_ELEVATION_M,
    PLANET_RADIUS_KM,
    line_spacing_rad,
)
from .plates import (
    Plate,
    cached_node_position_tree,
    collect_all_channel_depth,
    collect_all_channel_width,
    collect_all_crustal_thickness,
    collect_all_elev_change_reason,
    collect_all_elevation,
    collect_all_glacier_depth,
    collect_all_ice_load_deflection,
    collect_all_mantle_lithosphere_thickness,
    gather_node_positions,
    query_workers,
)

if TYPE_CHECKING:
    from .world import World

SLOPE_NEIGHBOR_COUNT = 4

# Starting points, not final -- see module docstring. Tuned by rough order-of-magnitude
# reasoning against boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr): at a
# moderately steep slope (~0.05) and moderate precipitation (~1000 mm/yr), this coefficient
# gives rain erosion the same order of magnitude as mountain-building uplift -- the design
# goal being that erosion should roughly balance typical uplift rates, not swamp them or be
# swamped by them (the raw number can't be picked from a grid-cell-based slope convention;
# see above for why).
RAIN_EROSION_COEFFICIENT = 6000.0
# Needs no rescaling the way the rain coefficient does, since wind speed's scale
# (MERIDIONAL_BASE_SPEED = 6.0) is fixed by this module's own climate pipeline, not tied to a
# grid-resolution-dependent slope convention.
WEATHERING_COEFFICIENT = 3.0
# Humidity level at which the weathering-humidity factor saturates to 1.0.
HUMIDITY_REFERENCE = 1.0
# Unlike rain/river erosion (both already multiply by `slope`), weathering had no relief
# dependence at all until this was added -- confirmed as the main reason flat, low-lying
# coastal land was drowning within a step or two of being created: wind/humidity-driven
# weathering ate a real coastal plain (whose slope is close to 0) exactly as fast as a
# mountain flank, even though its whole elevation buffer is often just a few meters. Real
# denudation is strongly relief-correlated -- flat lowlands erode far slower than steep
# terrain -- so `relief_factor` reintroduces that, same normalize-and-saturate idiom as
# `humidity_norm` above. WEATHERING_RELIEF_REFERENCE_SLOPE is the slope at which weathering
# reaches full strength; picked against a real generated world's own slope distribution
# (median land slope ~2.7e-4, far below this -- so typical flat terrain gets a small
# fraction of full weathering -- while land already at the ~90th percentile of steepness,
# ~3e-3, reaches full strength).
WEATHERING_RELIEF_REFERENCE_SLOPE = 0.002

# Stream-power river erosion: coefficient*channel_boost*water_accum^FLOW_EXPONENT*
# slope^SLOPE_EXPONENT. Coefficient re-derived by the same order-of-magnitude reasoning as
# RAIN_EROSION_COEFFICIENT (see module docstring); the two exponents need no such
# re-derivation, being dimensionless and not tied to any grid/unit convention.
RIVER_EROSION_COEFFICIENT = 100.0
RIVER_FLOW_EXPONENT = 0.5
RIVER_SLOPE_EXPONENT = 1.0
# A river preferentially re-carves its own established channel (real rivers meander within,
# not across, their valley). Channel depth is a real length (meters), so these values are
# plain physical quantities, not tied to any grid/unit convention that would need
# re-derivation.
CHANNEL_EROSION_BOOST = 0.6
CHANNEL_BOOST_REFERENCE_M = 200.0
MAX_CHANNEL_DEPTH_M = 2000.0

# Channel width is a mantle-bloom-original addition. "Larger flows make a channel larger":
# standard hydraulic-geometry scaling (width ~ discharge^b, b commonly observed near 0.5 for
# real rivers) -- same discharge exponent as RIVER_FLOW_EXPONENT above, but width is purely a
# function of how much water
# passes through, not of slope (a wide, lazy lowland river and a narrow, steep mountain
# torrent can carry comparable discharge, but width is driven by the water, not the
# gradient). Grows the same way channel_depth does -- persistent, monotonically
# non-decreasing, capped -- rather than a stateless function of this step's flow alone, so a
# river that temporarily dries up doesn't instantly narrow either.
WIDTH_GROWTH_COEFFICIENT = 50.0
WIDTH_FLOW_EXPONENT = 0.5
MAX_CHANNEL_WIDTH_M = 5000.0

# A big, slow river drops part of its sediment load locally (floodplain/delta) instead of
# carrying all of it to the coast. river_speed here is a stylized, unitless quantity (see
# hydrology.compute_river_speed) rather than a literal speed -- so its threshold isn't a
# physical speed either, just derived against this codebase's own river_speed scale.
# DEPOSITION_MIN_FLOW_M and DEPOSITION_FRACTION need no such derivation (both are already
# dimensionless/fractional, not tied to any particular scale convention).
DEPOSITION_SPEED_THRESHOLD = 2.0
DEPOSITION_MIN_FLOW_M = 0.05
DEPOSITION_FRACTION = 0.15

# Glacier erosion: scales with slope and actual accumulated ice depth (hydrology.py's
# glacier_depth, a real persistent field, not a stateless cold proxy) -- depth*slope
# approximates basal shear stress, a standard real glacial-erosion proxy: a flat-bottomed
# accumulation bowl still correctly erodes near zero regardless of ice depth.
# Temperature/precipitation-driven ice depth is a real length in consistent units (unlike the
# slope-based rain/river coefficients above, which needed re-derivation for a rise/run
# convention), so this coefficient doesn't need the same treatment. Coefficient and max factor
# both raised from an earlier, more timid starting point (0.05/2.0): basal shear stress under
# real ice scales with the ice's own *weight* (depth), and a thin valley glacier's few hundred
# meters versus a real continental ice sheet's kilometers is a large enough range that capping
# the erosive multiplier at 2x was leaving thick ice under-erosive relative to thin ice -- the
# higher cap lets genuinely deep, heavy ice keep scaling up over a wider depth range before
# saturating.
#
# GLACIER_EROSION_COEFFICIENT re-derived (see GitHub issue #145's investigation) against a
# published global compilation of measured erosion rates (Nature Geoscience 2025's "Drivers of
# global glacial erosion rates", building on Koppes & Montgomery 2009): glacial erosion's
# log-mean rate (~0.51 mm/yr, 99% of glaciers between 0.02-2.68 mm/yr, alpine tidewater
# glaciers -- the fastest-eroding type -- averaging ~2.2 mm/yr) sits within the *same order of
# magnitude* as fluvial erosion measured specifically in actively uplifting mountain terrain
# (1-10 mm/yr) -- not a full order of magnitude above it the way glacial-vs-*all*-rivers
# (including ordinary lowland ones, log-mean ~0.067 mm/yr) misleadingly suggests. Since
# RAIN_EROSION_COEFFICIENT above is itself calibrated for exactly that terrain (steep, wet,
# actively uplifting), the fair comparison is glacier-vs-rain at matched reference conditions,
# not glacier-vs-global-river-average -- so the coefficient is set equal to
# RAIN_EROSION_COEFFICIENT (the same "put it in the same bucket as the thing it's being
# compared to" precedent SEISMIC_EROSION_COEFFICIENT below already uses), rather than the
# ~650,000x jump a naive "glacial is 10x fluvial" headline would have implied. At their
# respective reference points (ice_factor=1, i.e. GLACIER_EROSION_REFERENCE_DEPTH_M of ice;
# 1000 mm/yr precipitation) the two terms now erode at the same rate for a given slope, with
# GLACIER_EROSION_MAX_FACTOR still letting genuinely thick ice (a real continental sheet) scale
# up to 4x that before saturating -- the ceiling a modest valley glacier's ~100 m couldn't
# plausibly reach.
GLACIER_EROSION_COEFFICIENT = 6000.0
GLACIER_EROSION_REFERENCE_DEPTH_M = 100.0
GLACIER_EROSION_MAX_FACTOR = 4.0

# Glacier flattening is a mantle-bloom-original addition modeling how real continental ice
# sheets grind down local relief over broad areas (e.g. the Canadian Shield/Fennoscandia read
# as glacially smoothed bedrock, not just eroded-lower). Implemented as a relaxation of each
# node's elevation toward the mean of its hydrology.py flow-graph neighbors (a genuine local
# blur, not a directional erosion/deposition), scaled by the same ice_factor glacier erosion
# uses -- reuses hydrology's own k=FLOW_NEIGHBOR_COUNT neighbor graph rather than a separate
# query. GLACIER_FLATTEN_RATE_PER_MYR has no established real-world value to draw on; picked
# as a starting point, checked against a live run the same way every other from-scratch rate
# in this codebase was. Raised alongside the erosion coefficient above, same "heavier ice grinds
# harder" reasoning.
GLACIER_FLATTEN_RATE_PER_MYR = 0.3

# Seismic erosion: earthquake-triggered landsliding, a real, well-documented contributor in
# young, actively-uplifting ranges (the Himalaya, the Andes) that this model didn't have a
# source for at all. mantle-bloom has no explicit fault/stress field to drive this from
# directly, so elevation itself stands in for "how much active convergent uplift has piled up
# here" -- in this model, sustained height above a couple thousand meters *is* the signature of
# an actively colliding, actively seismic belt (see plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR),
# there being no other elevation source that reaches these heights. Needs a slope to fail down,
# same rise/run definition rain/river erosion already use -- an earthquake can't landslide a
# perfectly flat plateau regardless of how high it sits. The height term is deliberately
# superlinear (SEISMIC_EROSION_ELEVATION_EXPONENT): real seismicity in young orogens grows
# faster than linearly with height/ongoing convergence, and this is also this model's main
# mechanism for capping just how tall an actively colliding range can get before its own
# erosion catches up with uplift -- a real "geomorphic ceiling" effect, not just cosmetic
# variation. SEISMIC_EROSION_ELEVATION_REFERENCE_M is picked at real Himalaya/Andes-base scale;
# SEISMIC_EROSION_COEFFICIENT by the same order-of-magnitude reasoning as
# RAIN_EROSION_COEFFICIENT (against CONVERGENT_MOUNTAIN_RATE_M_PER_MYR, 800 m/Myr): at a
# moderately steep mountain slope (~0.05) right at the reference elevation, this gives seismic
# erosion a meaningful fraction of the uplift rate -- a real but not dominant contributor at
# ordinary mountain heights, growing sharply toward and past truly extreme (Everest-scale)
# peaks.
SEISMIC_EROSION_COEFFICIENT = 6000.0
SEISMIC_EROSION_ELEVATION_REFERENCE_M = 3000.0
SEISMIC_EROSION_ELEVATION_EXPONENT = 2.0
SEISMIC_EROSION_MAX_HEIGHT_FACTOR = 3.0

# Earthquake-triggered landsliding burst: when faults.py records recent ruptures
# (World.earthquakes -- transient located events, see faults.Earthquake), the seismic-erosion
# term above is locally multiplied by (1 + burst) around each epicentre, where a single
# epicentre contributes EARTHQUAKE_EROSION_PEAK_BOOST * 10**(Mw - EARTHQUAKE_EROSION_MW_REF)
# at the epicentre, tapering linearly to zero at EARTHQUAKE_EROSION_REACH_KM_PER_MW * Mw.
# This is the direct-fault-driven counterpart to the elevation stand-in above -- now that
# there *is* a fault/stress signal, use it where it exists. Recency-weighted so a quake's
# landsliding pulse fades over the retention window rather than cutting off abruptly.
EARTHQUAKE_EROSION_PEAK_BOOST = 4.0
EARTHQUAKE_EROSION_MW_REF = 6.5
EARTHQUAKE_EROSION_REACH_KM_PER_MW = 12.0

# Mass-wasting runout (issue #275 phase 4). Landslide and rockfall debris moves by gravity, not
# with liquid water, so seismic erosion no longer joins the water-routed pool. That pool stops
# dead at a frozen node (hydrology gives frozen land no flow target), and nearly every column
# at the Hc cap is a frozen summit: its debris was "deposited" straight back where it came
# from, so the highest belts never lost any net rock. `_route_mass_wasting` instead carries the
# debris down the bare-rock steepest-descent path, frozen or not, and settles it where the
# ground flattens out -- the foreland. On ground at or above MASS_WASTING_RUNOUT_SLOPE (a
# mountain flank) none of it settles; on flat ground it all does, linearly in between. A
# column within MASS_WASTING_HEADROOM_TAPER_M of the Hc cap takes proportionally less, and none
# at the cap, so debris runs on past a saturated plateau rather than overflowing it. Debris
# that reaches the sea spreads onto the shelf and down into the ocean basin like the other
# marine sediment (_spread_marine_sediment) -- turbidity currents off a collision margin.
# MASS_WASTING_RUNOUT_SLOPE sits between a mountain flank (~0.05) and a foreland plain (median
# land slope ~3e-4, see WEATHERING_RELIEF_REFERENCE_SLOPE), a starting point like the others.
MASS_WASTING_RUNOUT_SLOPE = 0.005
MASS_WASTING_HEADROOM_TAPER_M = 2000.0
# A land pit fills only to its Hc room; the rest spills over the basin rim (hydrology's
# spill_target) and runs on. This many downhill sweeps, then any debris still spilling stays
# in its pit.
MASS_WASTING_SPILL_PASSES = 4

# Submarine erosion (mantle-bloom-original): plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR uplifts a
# colliding node regardless of whether it sits above or below sea level, so two *submerged*
# plates colliding will, unchecked, raise a range to the surface at the same rate a subaerial
# collision does. Real submerged relief grows much more slowly, because the sea floor has its
# own denudation the land lacks: deep bottom currents sweeping the flanks, and gravity-driven
# slope failure that the surrounding water does nothing to buttress (a steep, freshly-uplifted
# submarine scarp slumps under its own weight and the water pressure on it). Modeled as two
# slope-driven terms -- a current-driven baseline (SUBMARINE_EROSION_COEFFICIENT, slope alone)
# plus a pressure/depth term (SUBMARINE_PRESSURE_COEFFICIENT, slope x how deep the node still
# sits, saturating at SUBMARINE_PRESSURE_REFERENCE_M) -- so the brake is hardest on a deep
# abyssal range and eases continuously as a crest climbs toward the surface, handing off to
# ordinary subaerial erosion the moment it breaches. Slope-driven (same dimensionless rise/run
# as rain/river erosion) so a flat abyssal plain still erodes near zero; capped at the drop to
# the lowest neighbor like every other term here. Coefficients picked by the same
# order-of-magnitude reasoning as RAIN_EROSION_COEFFICIENT against
# CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr): at a moderate submarine-ridge slope (~0.02)
# and mid-depth (~half the reference), the two terms sum to roughly half the uplift rate -- a
# real, sustained drag on submarine orogeny that still lets a persistent collision eventually
# breach, not a hard ceiling.
SUBMARINE_EROSION_COEFFICIENT = 8000.0
SUBMARINE_PRESSURE_COEFFICIENT = 24000.0
SUBMARINE_PRESSURE_REFERENCE_M = 4000.0

# Coastal erosion (mantle-bloom-original): the shoreline itself -- and the crest of a mid-ocean
# range that rises into the same near-surface zone -- takes erosion both the open sea floor and
# the dry interior escape. Two mechanisms, both integrated per year (dt-scaled like every rate
# here): wave attack (COASTAL_EROSION_WAVE_RATE_M_PER_MYR, a flat rate across the band -- swell
# energy reaching the coast doesn't depend on the node's own relief), and frost shattering --
# water seeping into rock, freezing, expanding, prying it apart, then melting -- which needs the
# climate to actually *cycle* through freezing: it peaks at COASTAL_FROST_PEAK_C (just below 0C,
# where a seasonal/diurnal climate crosses the freezing point most often) and falls off both
# toward permanently-frozen and toward never-freezing climates, a Gaussian of width
# COASTAL_FROST_WIDTH_C in the node's mean temperature. COASTAL_EROSION_BAND_M is how far above
# and below sea level still counts as "coast" (a real wave-cut platform plus intertidal/spray
# zone is tens to low hundreds of meters of vertical reach); the rate tapers linearly to zero at
# the band edge. Eroded rock sheds seaward as marine sediment, same as submarine erosion.
COASTAL_EROSION_BAND_M = 200.0
COASTAL_EROSION_WAVE_RATE_M_PER_MYR = 400.0
COASTAL_FROST_MAX_RATE_M_PER_MYR = 500.0
COASTAL_FROST_PEAK_C = -2.0
COASTAL_FROST_WIDTH_C = 6.0

# Marine sediment spreading: submarine + coastal erosion both shed rock onto the surrounding sea
# floor (underwater slumping / longshore drift move it a short distance downslope, not along
# route_downstream's precipitation-fed river graph). A submerged source keeps
# MARINE_SEDIMENT_LOCAL_FRACTION where it eroded; a subaerial sea-cliff source sheds all of it
# into the water. The rest spreads, inverse-distance weighted, across up to
# MARINE_SPREAD_NEIGHBOR_COUNT of the nearest *lower* ocean nodes within MARINE_SPREAD_RANGE_KM
# (sediment runs downhill and settles into basins). A source with no lower ocean node in range
# keeps the full amount locally rather than losing it. Exactly conserves the eroded total.
MARINE_SEDIMENT_LOCAL_FRACTION = 0.35
MARINE_SPREAD_RANGE_KM = 120.0
MARINE_SPREAD_RANGE_RAD = MARINE_SPREAD_RANGE_KM / PLANET_RADIUS_KM
MARINE_SPREAD_NEIGHBOR_COUNT = 8

# Deposition, not just erosion -- eroded material has to go somewhere, and "wherever
# route_downstream's single water-flow graph happens to carry it" is only right for rain/
# river erosion. Four more pathways, alongside the existing river/runoff floodplain
# deposition above (DEPOSITION_*, already shared by every source through that one routed
# pool): a glacier's own debris load mostly drops close to where the ice picked it up
# (subglacial till), not carried far by meltwater; wind-eroded material mostly resettles a
# short real distance downwind (dust/sand -- genuinely different physics from water-driven
# transport, so it needs its own transport step, not the water flow_target graph); material
# that reaches the ocean spreads along nearby shallow coast as beach/nearshore sediment
# instead of piling entirely onto the single river-mouth node route_downstream happened to
# route it to; and material that reaches an enclosed lake/sea instead spreads across that
# whole standing body of water (see LAKE_SEDIMENT_UNIFORM_FRACTION/_spread_lake_sediment
# below) rather than piling onto the single sink node its basin's own catchment happens to
# drain to. All four conserve mass exactly, same as route_downstream's own retain_fraction --
# nothing here is lost, only moved (deep-ocean beach spreading is the one partial exception,
# see BEACH_SHELF_DEPTH_M below, and even that falls back to full local deposit rather than
# vanishing when no shallow water is in range).

# A glacier's scoured load splits two ways: GLACIER_TILL_FRACTION settles immediately,
# subglacial till dropped right where the ice picked it up; the rest travels *with the ice
# itself* rather than joining the water-routed pool -- apply_erosion reuses hydrology.
# route_downstream directly (the same "elevation-descending sweep along a flow_target graph,
# retain_fraction settles material locally" engine the ordinary river-deposition pool already
# uses below), but along hydrology.HydrologyFields.ice_flow_target (the ice's own real downhill
# flow path, not water's -- see that field's own comment for the frozen-node routing this
# depends on) instead of flow_target, retaining in full the moment a hop reaches a node with
# less than hydrology.GLACIER_VISIBLE_DEPTH_M of its own ice -- genuinely outside the glacier --
# rather than continuing to travel once there's no more ice there to carry it. A real terminal
# moraine/outwash deposit built beyond the ice margin, not debris stranded throughout the
# glacier's interior.
GLACIER_TILL_FRACTION = 0.5

# Aeolian transport: wind-eroded material moves with the wind, not downhill with water, so it
# needs its own single-hop transport rather than joining route_downstream's flow_target graph.
# WIND_DEPOSITION_FRACTION is how much of weathering's own erosion counts as this genuinely
# wind-carried fraction (the complement stays in the ordinary water-routed pool -- weathered
# material that rain then washes downhill); WIND_TRANSPORT_DISTANCE_KM is how far downwind it
# typically resettles before the next step's wind moves it again. Both starting points.
WIND_DEPOSITION_FRACTION = 0.4
WIND_TRANSPORT_DISTANCE_KM = 40.0
WIND_TRANSPORT_DISTANCE_RAD = WIND_TRANSPORT_DISTANCE_KM / PLANET_RADIUS_KM

# Marine/beach sediment: BEACH_DEPOSITION_LOCAL_FRACTION of whatever reaches the ocean at a
# given node stays right there (the river mouth/delta itself); the rest spreads across nearby
# shallow water (elevation above BEACH_SHELF_DEPTH_M -- a real continental-shelf-like depth,
# not the open abyssal seafloor) within BEACH_DEPOSITION_RANGE_KM, inverse-distance weighted.
# Deep-water coastlines with no shallow neighbour in range keep the full amount locally rather
# than losing it -- there's nowhere shallower to spread it to, not a reason to destroy mass.
BEACH_DEPOSITION_LOCAL_FRACTION = 0.4
BEACH_DEPOSITION_RANGE_KM = 150.0
BEACH_DEPOSITION_RANGE_RAD = BEACH_DEPOSITION_RANGE_KM / PLANET_RADIUS_KM
BEACH_SHELF_DEPTH_M = -200.0
BEACH_SPREAD_NEIGHBOR_COUNT = 8

# Coastal leveling feedback (mantle-bloom-original -- see GitHub issue #122, "Speckled low-relief
# coastlines"). Every source above is either purely subaerial or purely submarine, and none
# of them look at coastal *connectivity*: a marginally-submerged flat continental shelf
# sitting right on the waterline is a stable fixed point that just dithers land<->ocean
# node-by-node forever (per-node elevation noise > the surface's height above/below sea level),
# and the water-routed deposition pool makes it worse -- route_downstream funnels a slow
# river's whole load down one discretised channel and drops 50-250 m on a single band node
# while its neighbour gets ~0 (real rivers on a delta plain avulse into distributaries and
# spread that load laterally across the whole flat fan). This pass makes a coherent coast the
# stable state instead. It is *one symmetric operation* over a band that straddles sea level
# (COASTAL_LEVELING_BAND_M): every band node has a local target datum, and the pass grinds
# the nodes standing above it (wave-cut planation) and silts up the ones sitting below it
# (sheltered-shelf / interdistributary infill) toward that common datum -- conservatively
# (np.add.at, same as the _spread_* helpers), fed from the ground-off rock, a redirected
# share of the submarine/coastal erosion pool, and DELTA_REDIRECT_FRACTION of the clumped
# near-sea-level river-deposition lump (the distributary spread). Barrier islands and their
# back-barrier lagoons/marshes then emerge across steps from the same shelter field, with no
# explicit lagoon detection (see _spread_coastal_leveling's docstring).

# "Wave exposure" proxy: the fraction of nodes within COASTAL_OPENNESS_RANGE_KM that are open
# ocean (hydrology's connectivity-aware is_ocean, so an inland-lake shore reads as fully
# enclosed and is never planed). ~150 km is a coarse fetch scale -- far enough that a straight
# open coast sits near 0.5, a bay interior well below, a headland/peninsula tip above. Density-
# independent (a radius count, not a fixed k).
COASTAL_OPENNESS_RANGE_KM = 150.0
COASTAL_OPENNESS_RANGE_RAD = COASTAL_OPENNESS_RANGE_KM / PLANET_RADIUS_KM

# The single near-sea-level band the leveling pass operates on: |elevation - sea_level| within
# COASTAL_LEVELING_BAND_M, both sides. Deliberately far shallower than COASTAL_EROSION_BAND_M
# (200 m) -- this targets the low-relief drowned shelf / coastal plain specifically, not every
# sea cliff. Subsumes the old separate PLANATION_BAND_M (land, +60 m) and INFILL_DEPTH_M
# (ocean, -70 m) one-sided gates, which each only saw half the checkerboard.
COASTAL_LEVELING_BAND_M = 45.0

# Each band node's target datum (see leveling_datum_m) is a single continuous function of
# wave exposure: sea_level - LEVELING_PLATFORM_UNDERCUT_M * exposure + LEVELING_MARSH_CREST_M
# * shelter. A genuinely wave-exposed low sheet is cut down into open water (real shore
# platforms are planed to roughly low-tide level and a little below), a sheltered embayment
# silts up to a marsh sitting just proud of the waterline (biomes.classify_wetland), and a
# straight open coast lands near sea level. Cutting/filling to a datum that is *off* the
# waterline (rather than balanced exactly on it) is what lets the feedback *resolve* a
# checkerboard into a coast that follows the shelter field, not just damp its amplitude.
# Grinding rate = LEVELING_RATE_M_PER_MYR * coastal * proximity * prominence (coastal
# saturating once LEVELING_EXPOSURE_REF of the neighbourhood is open water, so a landlocked
# lowland is untouched; proximity tapers linearly to 0 at the band edge; prominence -- see
# PROMINENCE_* below -- boosts a protruding speck and nearly spares a flat sheet). Also caps
# the fill side's per-step per-node receipt (see _spread_coastal_leveling), so a hollow and
# its neighbouring bump converge on the datum at the same pace instead of the hollow slamming
# shut in one step. This is a from-scratch coefficient, order-of-magnitude-checked against
# boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr) and swept on real saves (round-1's
# 250 over-planed the drowned shelf and *raised* the node-flip rate; 60 both cuts the dither
# and holds |dElev|/step in the near-sea-level band well down -- flip fraction ~0.40 -> ~0.29,
# band |dElev| p90 ~45 m -> ~25 m, across two seeds, an apply_erosion-only loop).
LEVELING_RATE_M_PER_MYR = 60.0
LEVELING_EXPOSURE_REF = 0.3
LEVELING_PLATFORM_UNDERCUT_M = 6.0
LEVELING_MARSH_CREST_M = 4.0
# A fill sink must have at least this fraction of open water in its neighbourhood -- an
# inland-lake shore reads ~0 on the connectivity-aware openness field and must never silt up
# through this coastal pass (it has its own lake siltation). Just above 0, not a real gate on
# genuine coast.
LEVELING_MIN_OPENNESS = 0.03

# The fill side is a capped water-fill (see _spread_coastal_leveling): each below-datum band
# node has a hard per-step capacity -- its own metres of room to its datum, so nothing is
# raised past the datum in one step -- and LEVELING_FILL_ITERS passes distribute each source's
# remaining load across its still-open sinks (up to LEVELING_SPREAD_NEIGHBOR_COUNT of them
# within INFILL_RANGE_KM), weighted by shelter x hollow x room-left / distance, capping and
# carrying the overflow to the next pass. That is what makes a 200 m lump genuinely spread
# across the flat plain -- each near sink takes only its ~30 m of room and the rest flows on
# -- instead of piling onto the single most-weighted node (the failure mode the round-1
# _spread_coastal_infill still had). LEVELING_LOCAL_FRACTION stays put at a source that is
# itself a below-datum node (its own hollow to fill -- a river emptying into a genuine basin
# isn't fully swept clean); a pure source spreads in full. Whatever no sink had room for after
# the last pass stays on the source (net-zero, no mass loss).
INFILL_RANGE_KM = 120.0
INFILL_RANGE_RAD = INFILL_RANGE_KM / PLANET_RADIUS_KM
INFILL_SHELTER_REF = 0.5
LEVELING_SPREAD_NEIGHBOR_COUNT = 48
LEVELING_FILL_ITERS = 3
LEVELING_LOCAL_FRACTION = 0.25
# How much of each step's submarine + coastal erosion is redirected from _spread_marine_sediment
# (which only runs spoil *downhill to deeper* water -- actively counterproductive in a shallow
# embayment) into the leveling fill sink, which gets first call on it. The rest still spreads
# to deep water as before. Feeding the band matters: on real saves, routing more of this pool
# (and more of the distributary redirect) into the capped water-fill is what actually pulls
# the near-sea-level flip rate and |dElev|/step down -- the fill can only ever raise a node to
# its datum, so there is no runaway, and whatever the thin band can't hold bounces back to the
# eroding node.
COASTAL_INFILL_MARINE_FRACTION = 0.5

# Distributary redirect: a big, slow river on the near-sea-level plain (the existing
# DEPOSITION_* "depositing" test -- slow river_speed, real flow) has route_downstream drop its
# retained load on whichever single node its discretised flow path runs through. A real delta
# splits into distributaries and spreads that load across the whole flat fan. This fraction of
# such a band node's water-routed deposit is pulled back out and fed into the leveling fill
# spread (which scatters it across the band's hollows -- the emergent distributary fan); the
# rest stays put as the active channel / natural levee. Swept high on real saves (0.8 beats
# 0.6 beats 0) -- the capped water-fill trickles it out over several steps anyway, so an
# aggressive redirect is safe and de-clumps faster.
DELTA_REDIRECT_FRACTION = 0.8

# Barrier islands: a shallow near-sea-level sink that has land within BARRIER_LANDWARD_KM (a
# shore to parallel) but still faces open water (openness >= BARRIER_MIN_OPENNESS -- it's on
# the outer edge of the shelf, not deep in a bay) is a barrier candidate. It gets a priority
# multiplier on the infill weight (longshore drift piles sediment onto the bar first) and its
# fill cap is raised BARRIER_CREST_M above sea level, so the band can just breach the surface.
# The water it then encloses loses open-ocean neighbours, so next step its own openness drops,
# its shelter rises, and the back-barrier lagoon silts up to a flat near-sea-level sheet that
# biomes.classify_wetland reads as marsh -- all emergent, no lagoon flag.
BARRIER_LANDWARD_KM = 45.0
BARRIER_LANDWARD_RAD = BARRIER_LANDWARD_KM / PLANET_RADIUS_KM
BARRIER_MIN_OPENNESS = 0.3
BARRIER_PRIORITY = 3.0
BARRIER_CREST_M = 3.0

# Prominence: waves plane protrusions and fill hollows, so a near-sea-level node that stands
# `PROMINENCE_REF_M` above its own neighbourhood mean is planed ~PROMINENCE_MAX times harder
# than a flat one (and one sitting below its neighbours ~0), while a band node that sits
# *below* its neighbourhood mean is the preferred fill sink. This is what actually collapses a
# pixel-by-pixel land<->ocean checkerboard into a coherent shoreline -- openness alone,
# measured at a ~150 km fetch scale, is far too smooth to make that call at ~60 km node
# spacing. It only reweights the existing grind / fill (no new mass term): the rock a
# prominence-boosted node loses still goes through the conservative leveling pool.
PROMINENCE_REF_M = 18.0
PROMINENCE_MAX = 3.0
# A barrier candidate whose openness is above INFILL_SHELTER_REF has shelter == 0; without a
# floor it would attract no sediment and never form. This is the minimum attractiveness a
# barrier candidate keeps regardless of shelter.
BARRIER_ATTRACT_FLOOR = 0.35

# Lake/sea deposition: an enclosed standing body of water has none of the wave/current energy
# that keeps ocean sediment moving toward deeper water (see SUBMARINE_EROSION_*/
# _spread_marine_sediment above) -- it's quiet enough that essentially all the sediment a river
# carries into it settles somewhere inside the basin, rather than being carried on through (there
# *is* nowhere further to go -- a lake/sea's own sink node has no lower neighbor by construction,
# unlike an ordinary reach along a river's course). Left alone, route_downstream's ordinary
# dead-end rule (see hydrology.route_downstream) piles the *entire* accumulated sediment load
# onto that single sink node every step -- a single growing spike right at the inflow point, not
# a real lake-bed blanket. `_spread_lake_sediment` instead redistributes whatever lands on a lake
# member node across that lake's whole currently-flooded extent (`hydrology.lake_components`, the
# same connected-body grouping hydrology.py's own breach-erosion term uses), weighted
# LAKE_SEDIMENT_UNIFORM_FRACTION uniform across every member (fine suspended sediment settling
# out of quiet water genuinely blankets a whole lake floor fairly evenly) plus the remainder
# weighted by each member's own current `lake_depth` -- a real lake fills its deepest holes
# fastest, since there's more standing water (and so more settling sediment reaching the bed)
# above a deep spot than a shallow one. Depth-weighting is what actually flattens the basin over
# repeated steps rather than just relocating one pile to another: it's a relaxation toward
# uniform depth, i.e. the discrete analogue of "valleys fill in faster than hills, converging on
# a flat floor" -- the enclosed-basin counterpart to _spread_marine_sediment's open-ocean
# spreading. Conserves the redistributed total exactly, via np.add.at, same as every other
# _spread_* helper here; a sink with no real standing water yet (a dry playa/closed basin) is
# left untouched, keeping this pathway's old, un-lake-aware single-point-pile behavior there.
LAKE_SEDIMENT_UNIFORM_FRACTION = 0.5

# Depositional Hc-cap overflow (issue #288). Every pathway above settles its load without
# regard to the receiving column's room under lithosphere.MAX_CRUSTAL_THICKNESS_M (only
# landslide runout looks at it), so a foreland, lake floor or shelf node that tectonics has
# already driven to the cap used to have the excess clipped off at the Hc write-back -- and its
# continental share booked as numerical loss. `_carry_overflow` instead carries what a full
# receiver can't hold on to the next receiver with room, by the transport its location implies:
# downhill (steepest descent, spilling over a basin rim) on land, a capped spread across the
# whole lake in a flooded basin, a capped spread onto the lower shelf and basin at sea
# (MARINE_SPREAD_*), and along the ice's own flow path under a glacier. The load travels mixed,
# so every receiver gets its share of each pathway and of the continental tracer. Each name in
# OVERFLOW_PATHWAYS is one diagnostic column in the erosion budget.
OVERFLOW_PATHWAYS = ("river", "lake", "lake_silt", "glacial", "aeolian", "mass_wasting", "coastal_leveling", "marine")
OVERFLOW_GLACIAL = OVERFLOW_PATHWAYS.index("glacial")
OVERFLOW_LAKE_SILT = OVERFLOW_PATHWAYS.index("lake_silt")
# Every column but lake silt: material some erosion source gave up this step.
OVERFLOW_SEDIMENT = np.arange(len(OVERFLOW_PATHWAYS)) != OVERFLOW_LAKE_SILT
# Hops the carried load may take before whatever is still moving counts as stuck -- far more
# than any real descent path, so it only binds on a routing cycle.
OVERFLOW_MAX_HOPS = 1024
# Load still stuck after the transport above (a saturated pit with no outlet, a full lake with
# no spill, a deep-sea node with nothing lower in range, a cycle) fills outward from there in
# rising elevation order -- the full basin overflows its lowest saddle -- onto the first nodes
# with room. What that can't place goes to the nearest receivers with room: the
# OVERFLOW_SEARCH_NEIGHBOR_COUNT nearest nodes that still have room. Both stay within
# OVERFLOW_SEARCH_RANGE_KM -- far enough to reach past the edge of a Tibet-scale saturated
# plateau (~1000 km across; quad cells are ~125 km apart at the default density) to its
# flanks. The search hands load out inverse-distance weighted, in capped passes that each fill
# at least one oversubscribed receiver (OVERFLOW_FILL_PASSES, enough to reach every candidate,
# so a load many receivers' room thick still finds them all). What even that can't place is
# the terminal remainder: a column at the cap
# that keeps being loaded sheds the excess from its root instead (load-driven delamination of
# lower crust pushed past the eclogite transition), booked as the continental ledger's
# `overloaded_root_delaminated_m3` -- distinct from suture-accretion delamination.
OVERFLOW_STUCK_REASONS = ("pit", "closed_lake", "sea_floor", "hop_limit")
OVERFLOW_SEARCH_RANGE_KM = 1000.0
OVERFLOW_SEARCH_RANGE_RAD = OVERFLOW_SEARCH_RANGE_KM / PLANET_RADIUS_KM
OVERFLOW_SEARCH_NEIGHBOR_COUNT = 128
OVERFLOW_FILL_PASSES = OVERFLOW_SEARCH_NEIGHBOR_COUNT
# Overflow carried under a glacier rides the ice rather than settling under it (the same rule
# glacial transport and landslide debris follow), and the debris-laden ice scours its bed as it
# goes: each node it crosses gives up OVERFLOW_ICE_SCOUR_FRACTION of the passing load's
# thickness-equivalent, at most OVERFLOW_ICE_SCOUR_MAX_M_PER_MYR, within the node's remaining
# Hc and neighbour-drop headroom. The scoured rock joins the load (conserved like everything
# else here) and deepens channel_depth -- a glacial trough that rivers inherit once the ice is
# gone, the same way glacial abrasion's own channel growth does. The per-node cap keeps the
# load's growth along a long ice path additive rather than compounding.
OVERFLOW_ICE_SCOUR_FRACTION = 0.1
OVERFLOW_ICE_SCOUR_MAX_M_PER_MYR = 100.0

# Mobile cover (issue #297 phase 2). Loose sediment and regolith is not intact rock: every
# column carries `mobile_cover_m`, the unconsolidated share of its Hc at the top of the
# column, and `mobile_cover_continental_m`, the continental-derived share of that cover (a
# share of `continental_material_m`, so the existing ledger is unchanged). Everything any
# pathway settles -- river, lake, till, moraine, wind, landslide, beach, marine, leveling fill,
# flattening, overflow and lake silt -- enters the cover, so weathered rock joins it where it
# settles, and every removal takes cover before substrate. The erosion laws above are
# substrate laws: a rate law (subaerial, submarine and coastal erosion) entrains cover
# MOBILE_COVER_ERODIBILITY_FACTOR times faster than the same forcing detaches substrate, and
# reaches substrate only once the cover is gone (cover shields its bed; no tool effect is
# modelled), so a bare column erodes exactly as it did before the cover existed. Craton
# resistance is a property of coherent substrate, so it scales only the substrate share.
# SPACE 1.0 (Shobe et al. 2017) makes the same split; sediment/bedrock erodibility ratios of
# several-fold are typical, so this is a starting point. There is no in-place regolith
# production yet: the laws above are calibrated against real denudation, and a production
# term would let every low-relief column entrain up to this factor times its calibrated rate
# -- a global retune, left to the coupled solve.
MOBILE_COVER_ERODIBILITY_FACTOR = 5.0
# Burial consolidates cover: what lies deeper than MOBILE_COVER_CONSOLIDATION_DEPTH_M below the
# cover's top lithifies back into substrate with an e-folding time of
# MOBILE_COVER_CONSOLIDATION_TIMESCALE_MYR. Hc does not change; the continental share goes with
# it. Both starting points (diagenesis to rock takes km of burial and Myr-to-tens-of-Myr).
MOBILE_COVER_CONSOLIDATION_DEPTH_M = 1000.0
MOBILE_COVER_CONSOLIDATION_TIMESCALE_MYR = 20.0


@dataclass
class ErosionResult:
    """This step's per-node erosion breakdown, returned by apply_erosion for geology.py's own
    soil/coal/oil-gas formation (see that module) to reuse directly rather than re-deriving --
    the same "don't recompute, reuse what a prior step already produced" precedent World.
    climate_cache/World.hydrology_cache already set, except this one is threaded as a plain
    function return/argument rather than cached on World, since (unlike climate/hydrology) it
    has exactly one consumer, called from the same step_world function that produced it -- no
    later same-turn caller (rendering, /world/stats) ever needs it. Index-aligned with
    World.hydrology_cache's own arrays (points/plates_in_order/is_ocean in particular -- this
    dataclass deliberately doesn't repeat those, geology.py reads them off
    World.hydrology_cache directly, including its own sea_level_m-aware is_ocean rather than
    this module's own internal elevation<=0.0 shorthand): both come from an identical
    per-plate gather over the same (unreordered) world.plates within the same step_world
    call.

    `sediment_deposited` is every deposition pathway's combined total at each node -- ordinary
    river/runoff floodplain deposit, glacier till, glacier-transported moraine/outwash material,
    landslide runout (see MASS_WASTING_*), wind-blown resettling, beach/nearshore spreading, and marine sediment shed onto the sea
    floor by submarine and coastal erosion, all summed together (see apply_erosion's own comment
    for how they're split and redistributed) -- not just the water-routed share alone.

    `net_elevation_change_m` is this step's total geomorphic elevation delta per node
    (post-erosion `new_elevation` minus the pre-erosion `elevation` above): erosion removed,
    every deposition pathway added, plus the glacial-flattening exchange and the small
    non-conservative lake-siltation term -- but *not* tectonic deform, isostasy, or
    volcanism, which move elevation outside this module. Signed (negative where the step net-lowered a node, positive where it
    net-raised it). Retained on `World.erosion_cache` for the Geomorph Rate debug view
    (`render_image._render_geomorph_view`) -- the lumpiness of near-sea-level deposition (a
    +200 m spike on one node, ~0 on its neighbour) is invisible in every other view but is
    the whole coastal-speckle mechanism (see GitHub issue #122).

    `is_river_depositing` is exactly the mask DEPOSITION_SPEED_THRESHOLD/DEPOSITION_MIN_FLOW_M
    already select internally -- "a big, slow river is actively settling its sediment load
    here this step" -- exposed so geology.py can reuse the same, already-tuned "big and slow"
    test for its own riparian-vegetation boost (see that module's RIPARIAN_* constants)
    instead of re-deriving river_speed/water_accum_m a second time."""

    points: np.ndarray
    elevation: np.ndarray  # this step's *pre*-erosion elevation (same array hydro.elevation holds)
    slope: np.ndarray
    rain: np.ndarray
    river: np.ndarray
    weathering: np.ndarray
    sediment_deposited: np.ndarray
    is_river_depositing: np.ndarray  # bool -- a big, slow river settling sediment here this step
    net_elevation_change_m: np.ndarray
    temperature_c: np.ndarray
    precipitation_mm: np.ndarray
    # Issue #275: this step's area-weighted material budget (m^3) -- see apply_erosion.
    budget: dict[str, float] | None = None
    # Issue #275 phase 3: the vertical stress (Pa) this step's grounded ice puts on each
    # column, its change since last step (negative = unloading), and the resulting
    # isostatic elevation change (m) -- also *not* part of net_elevation_change_m.
    ice_load_pa: np.ndarray | None = None
    ice_load_change_pa: np.ndarray | None = None
    ice_deflection_change_m: np.ndarray | None = None


def _strip_mobile_cover(cover_m: np.ndarray, capacity_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """One step of a rate law acting on a column with `cover_m` of mobile cover on top -- see
    MOBILE_COVER_*. `capacity_m` is the substrate the step's forcing would detach from bare
    rock. The cover goes first, MOBILE_COVER_ERODIBILITY_FACTOR times faster; whatever part of
    the step is left once it is gone detaches substrate at the bare rate. Returns (cover
    entrained, substrate detached), meters. Exact for constant forcing, so one long step and
    several shorter ones agree."""
    capacity = np.asarray(capacity_m, dtype=float)
    entrained = np.minimum(cover_m, MOBILE_COVER_ERODIBILITY_FACTOR * capacity)
    detached = np.clip(capacity - entrained / MOBILE_COVER_ERODIBILITY_FACTOR, 0.0, None)
    return entrained, detached


def _gather_nodes(
    world: "World",
    node_cloud: tuple[np.ndarray, list[Plate]] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[Plate]]:
    """Every node's world position, elevation, and prior channel_depth/channel_width/
    glacier_depth, concatenated, alongside the ordered list of plates that contributed them --
    unlike reassign.py's own _gather_nodes (which needs per-node plate/line identity, since
    nodes there can move between lines), this only needs each field in the same order
    `points` and the write-back loop below (`plate.set_fields_on_plate`) already agree
    on, since erosion never moves a node or changes line topology. `node_cloud`, when passed
    (see apply_erosion), reuses an already-gathered (points, plates_in_order) pair instead of
    re-deriving every node's world position from scratch -- see plates.gather_node_positions's
    own docstring for why."""
    points, plates_in_order = node_cloud if node_cloud is not None else gather_node_positions(world.plates)
    if not plates_in_order:
        return np.zeros((0, 3)), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), []
    return (
        points,
        collect_all_elevation(plates_in_order),
        collect_all_channel_depth(plates_in_order),
        collect_all_channel_width(plates_in_order),
        collect_all_glacier_depth(plates_in_order),
        plates_in_order,
    )


def _gather_areas(world: "World", plates_in_order: list[Plate]) -> np.ndarray:
    """Every node's accounting area (m^2), in `_gather_nodes` order -- the same per-plate
    `accounting_areas_m2` continental_ledger.py integrates the tracer over, so erosion's volume
    transfers and the ledger agree exactly (exact cell areas on quad surfaces)."""
    spacing = line_spacing_rad(world.node_density)
    if not plates_in_order:
        return np.zeros(0)
    return np.concatenate([p.accounting_areas_m2(spacing) for p in plates_in_order])


def _earthquake_erosion_multiplier(world: "World", points: np.ndarray) -> np.ndarray:
    """Per-node `(1 + burst)` factor for the seismic-erosion term: a recency-weighted,
    magnitude-scaled bump around every recent epicentre in `World.earthquakes`, tapering
    linearly to 1.0 (no boost) beyond a per-quake reach. All-ones when nothing has ruptured
    recently -- so this is exactly a no-op in `"boundary"` mode with no active fast faults."""
    quakes = getattr(world, "earthquakes", None)
    if not quakes or len(points) == 0:
        return np.ones(len(points))
    boost = np.zeros(len(points))
    retain_years = faults.EARTHQUAKE_RETAIN_MYR * 1_000_000.0
    epicentres = np.array([q.epicenter_world for q in quakes])
    tree = cached_node_position_tree(world, points)
    for q, epi in zip(quakes, epicentres):
        recency = np.clip(1.0 - (world.elapsed_years - q.birth_years) / max(retain_years, 1.0), 0.0, 1.0)
        if recency <= 0.0:
            continue
        reach_rad = EARTHQUAKE_EROSION_REACH_KM_PER_MW * q.magnitude / PLANET_RADIUS_KM
        near = tree.query_ball_point(epi, reach_rad)
        if not near:
            continue
        near = np.asarray(near)
        d = geometry.angular_distance(points[near], epi)
        taper = np.clip(1.0 - d / reach_rad, 0.0, 1.0)
        peak = EARTHQUAKE_EROSION_PEAK_BOOST * 10.0 ** (q.magnitude - EARTHQUAKE_EROSION_MW_REF)
        boost[near] += peak * recency * taper
    return 1.0 + boost


def compute_slope(points: np.ndarray, elevation: np.ndarray, world: "World | None" = None) -> tuple[np.ndarray, np.ndarray]:
    """Per-node dimensionless rise/run -- elevation drop to the *lowest* of each node's
    SLOPE_NEIGHBOR_COUNT nearest neighbors (0 if this node is already a local minimum -- the
    "slope to lowest neighbor" definition used throughout this module), divided by the real
    great-circle distance to that specific neighbor -- alongside the raw drop in meters,
    unnormalized. Returns (slope, drop_m): `slope` drives the erosion rate formula, `drop_m`
    caps it (see apply_erosion) so a single step can't erode a node *past* its lowest
    neighbor's own elevation, which would carve a new pit lower than the valley it drains
    into. The two are tracked as separate return values because `slope` here is normalized
    (dimensionless rise/run) while the cap needs the raw, unnormalized drop -- capping against
    the normalized value would bound elevation change in the wrong units entirely. `world`,
    when passed (apply_erosion's own call does), shares this tree with every other per-step
    full-node-cloud query via plates.cached_node_position_tree instead of rebuilding an
    equivalent one from scratch -- see that function's own docstring."""
    n = len(points)
    if n <= SLOPE_NEIGHBOR_COUNT:
        return np.zeros(n), np.zeros(n)

    tree = cached_node_position_tree(world, points)
    _, neighbor_idx = tree.query(points, k=SLOPE_NEIGHBOR_COUNT + 1, workers=query_workers(n))
    neighbor_idx = neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0

    neighbor_elevation = elevation[neighbor_idx]
    rows = np.arange(n)
    lowest_col = np.argmin(neighbor_elevation, axis=1)
    lowest_elevation = neighbor_elevation[rows, lowest_col]
    lowest_idx = neighbor_idx[rows, lowest_col]

    drop_m = np.clip(elevation - lowest_elevation, 0.0, None)
    run_m = geometry.angular_distance(points, points[lowest_idx]) * PLANET_RADIUS_KM * 1000.0
    run_m = np.maximum(run_m, 1.0)  # avoid a divide-by-zero for (near-)coincident points
    return drop_m / run_m, drop_m


def submarine_erosion_amount(elevation: np.ndarray, slope: np.ndarray, is_ocean: np.ndarray, dt_myr: float) -> np.ndarray:
    """Meters eroded this step from submerged terrain by bottom currents and pressure-driven
    slope failure -- see SUBMARINE_EROSION_* constants for the reasoning. A current-driven
    baseline (slope alone) plus a term that scales with how deep the node still sits (water
    pressure / column height above it, saturating at SUBMARINE_PRESSURE_REFERENCE_M), the whole
    thing multiplied by `slope` so a flat abyssal plain erodes near zero and a steep,
    freshly-uplifted submarine scarp erodes fast. Zero on subaerial nodes (they get the ordinary
    subaerial sources instead). Not yet capped at the drop to the lowest neighbor -- apply_erosion
    does that, jointly with coastal erosion, against whatever drop subaerial erosion left."""
    depth_below_sea_m = np.clip(-elevation, 0.0, None)
    rate = SUBMARINE_EROSION_COEFFICIENT + SUBMARINE_PRESSURE_COEFFICIENT * np.clip(
        depth_below_sea_m / SUBMARINE_PRESSURE_REFERENCE_M, 0.0, 1.0
    )
    return np.where(is_ocean, rate * slope * dt_myr, 0.0)


def coastal_erosion_amount(elevation: np.ndarray, temperature_c: np.ndarray, dt_myr: float) -> np.ndarray:
    """Meters eroded this step from the near-sea-level band (|elevation| within
    COASTAL_EROSION_BAND_M, tapering linearly to zero at the band edge) by wave attack and
    freeze-thaw frost shattering -- see COASTAL_EROSION_* constants. Wave attack is a flat rate
    across the band; the frost term is a Gaussian in temperature peaked at COASTAL_FROST_PEAK_C
    (a climate has to cycle through freezing for frost wedging to do anything). Applies on both
    sides of the shoreline, so it also gnaws at the crest of a mid-ocean range that rises into
    the band. Not relief-gated (unlike weathering) -- a flat wave-cut bench erodes as readily as
    a cliff. Not yet capped -- apply_erosion caps it jointly with submarine erosion."""
    proximity = np.clip(1.0 - np.abs(elevation) / COASTAL_EROSION_BAND_M, 0.0, 1.0)
    frost = np.exp(-(((temperature_c - COASTAL_FROST_PEAK_C) / COASTAL_FROST_WIDTH_C) ** 2))
    rate = COASTAL_EROSION_WAVE_RATE_M_PER_MYR + COASTAL_FROST_MAX_RATE_M_PER_MYR * frost
    return rate * proximity * dt_myr


def climate_grid_indices(world_xyz: np.ndarray, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Nearest climate-grid (row, col) for each node's world position -- direct array
    indexing, not a tree lookup, since (unlike the geology side) the climate grid is
    already a plain regular lat/lon lattice. Mirrors climate._build_grid's own convention
    exactly: row 0 = north pole, row increases southward; column increases eastward,
    wrapping at the antimeridian."""
    lat, lon = geometry.xyz_to_latlon(world_xyz)
    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    row = np.clip(np.floor((90.0 - lat_deg) / (180.0 / height)).astype(int), 0, height - 1)
    col = np.floor((lon_deg + 180.0) / (360.0 / width)).astype(int) % width
    return row, col


def _flatten(
    hydro: "hydrology.HydrologyFields",
    ice_factor: np.ndarray,
    years: float,
    removable_m: np.ndarray | None = None,
    multiplier: float = 1.0,
) -> np.ndarray:
    """Glacier flattening (mantle-bloom-original, see module docstring): ice relaxes each
    node toward its hydrology.py flow-graph neighbors, scaled by GLACIER_FLATTEN_RATE_PER_MYR
    and the same ice_factor glacier erosion uses -- a genuine local blur (can raise a valley
    or lower a peak).

    Issue #275: done as a pairwise exchange rather than a per-node relaxation toward the
    neighbour mean, so it moves rock rather than creating/destroying it. A node passes
    `relax * (its height - lower neighbour's height) / k` thickness to each lower neighbour
    (k = neighbour count) -- the same per-edge rate the old relaxation applied, only the
    downhill half of it, since the lower node receives what the higher one sheds. Returns the
    (n, k) thickness each node sends along each `hydro.neighbor_idx` edge, scaled down per
    node so it never sends more than `removable_m` (its Hc headroom above the floor).
    `multiplier` is World.glacier_erosion_multiplier, applied before that cap."""
    years_myr = years / 1_000_000.0
    relax = 1.0 - np.exp(-GLACIER_FLATTEN_RATE_PER_MYR * ice_factor * years_myr)
    k = hydro.neighbor_idx.shape[1]
    drop = hydro.elevation[:, None] - hydro.elevation[hydro.neighbor_idx]
    send = multiplier * relax[:, None] * np.clip(drop, 0.0, None) / max(k, 1)
    if removable_m is not None:
        total = send.sum(axis=1)
        scale = np.divide(removable_m, total, out=np.ones_like(total), where=total > removable_m)
        send = send * scale[:, None]
    return send


def _apply_ice_load(
    elevation: np.ndarray,
    prior_deflection: np.ndarray,
    glacier_depth: np.ndarray,
    hc: np.ndarray,
    hm: np.ndarray,
    rho_c: np.ndarray,
    sea_level_m: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Move `elevation` (which carries `prior_deflection` of ice depression) to the
    depression this step's `glacier_depth` implies. Returns (new elevation, new deflection,
    grounded load in kg/m^2) -- issue #275 phase 3.

    The load is evaluated on the bed with the old deflection taken out, so the depression
    can't feed back into its own flotation test. Elevation moves by only the change in
    deflection, so melting the ice rebounds the surface by exactly what loading took. The
    stored deflection is what the elevation bounds actually let through, not the raw Airy
    response: storing more than was applied would rebuild too high an unloaded bed on melt."""
    unloaded_bed = elevation - prior_deflection
    load = lithosphere.grounded_ice_load_kg_m2(glacier_depth, unloaded_bed, sea_level_m)
    target = lithosphere.ice_load_deflection(hc, hm, rho_c, load)
    new_elevation = np.clip(unloaded_bed + target, MIN_ELEVATION_M, MAX_ELEVATION_M)
    return new_elevation, new_elevation - unloaded_bed, load


def _route_mass_wasting(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    slope: np.ndarray,
    source_vol: np.ndarray,
    source_tagged: np.ndarray,
    capacity_vol: np.ndarray,
    headroom_m: np.ndarray,
    spill_target: np.ndarray | None = None,
    on_ice: np.ndarray | None = None,
    inflow_vol: np.ndarray | None = None,
    inflow_tagged: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Gravity runout of landslide debris (see MASS_WASTING_*). Each node hands its debris to
    its lowest strictly-lower `neighbor_idx` neighbour -- bare rock, frozen or not, with no
    lake-spill or ice edges, so elevation-descending order visits every node before its target.
    Debris leaves its source; from the first node downslope on, a share settles at each node it
    passes (MASS_WASTING_RUNOUT_SLOPE's linear taper, scaled by the column's Hc headroom), never
    more than `capacity_vol` there. Debris reaching an ocean node stops there.

    A land pit fills up to its capacity, and the rest spills over the basin rim to
    `spill_target` (hydrology's one-hop basin escape) and runs on from there in a further sweep,
    up to MASS_WASTING_SPILL_PASSES sweeps. Whatever is still left -- a pit with no spill
    target, or debris still spilling after the last sweep -- settles in its pit, room or not.

    Debris that starts on, or runs onto, a land node marked `on_ice` stops there and is handed
    back as "ice arrival": a glacier carries rockfall on and in the ice to its margin (lateral
    and medial moraine), and out through a basin's outlet once the ice overtops the rim, so the
    caller carries it along `ice_flow_target` and hands what the ice drops back here as
    `inflow_vol`: debris arriving at a node rather than leaving one, so it can settle right
    where it is put down.

    Returns (land deposit, ocean arrival, ice arrival, and the tagged share of each), all
    volumes. `source_tagged` is the continental share of `source_vol`; it travels mixed with the
    rest, so every deposit carries the tagged fraction of the flux it came from. Conserves both
    totals exactly."""
    n = len(elevation)
    deposit = np.zeros(n)
    deposit_tagged = np.zeros(n)
    arrival = np.zeros(n)
    arrival_tagged = np.zeros(n)
    ice = np.zeros(n)
    ice_tagged = np.zeros(n)
    inflow = inflow_vol if inflow_vol is not None else np.zeros(n)
    if n == 0 or not (np.any(source_vol > 0) or np.any(inflow > 0)):
        return deposit, arrival, ice, deposit_tagged, arrival_tagged, ice_tagged

    rows = np.arange(n)
    neighbor_elevation = elevation[neighbor_idx]
    lowest = np.argmin(neighbor_elevation, axis=1)
    target = np.where(neighbor_elevation[rows, lowest] < elevation, neighbor_idx[rows, lowest], -1)
    settle = np.clip(1.0 - slope / MASS_WASTING_RUNOUT_SLOPE, 0.0, 1.0)
    settle = settle * np.clip(headroom_m / MASS_WASTING_HEADROOM_TAPER_M, 0.0, 1.0)
    spill = spill_target if spill_target is not None else np.full(n, -1)
    glaciated = (on_ice if on_ice is not None else np.zeros(n, dtype=bool)) & ~is_ocean

    # Debris shed onto the ice rides it; a source with nowhere lower to send its debris keeps it.
    rides = (source_vol > 0) & glaciated
    ice[rides] += source_vol[rides]
    ice_tagged[rides] += source_tagged[rides]
    stays = (source_vol > 0) & (target < 0) & ~glaciated
    deposit[stays] += source_vol[stays]
    deposit_tagged[stays] += source_tagged[stays]
    moves = (source_vol > 0) & (target >= 0) & ~glaciated
    flux = np.zeros(n)
    flux_tagged = np.zeros(n)
    np.add.at(flux, target[moves], source_vol[moves])
    np.add.at(flux_tagged, target[moves], source_tagged[moves])
    flux += inflow
    if inflow_tagged is not None:
        flux_tagged += inflow_tagged

    order = np.argsort(-elevation).tolist()
    target_list = target.tolist()
    spill_list = spill.tolist()
    settle_list = settle.tolist()
    room = np.clip(capacity_vol, 0.0, None).tolist()
    is_ocean_list = is_ocean.tolist()
    glaciated_list = glaciated.tolist()
    deposit_list = [0.0] * n
    deposit_tagged_list = [0.0] * n
    flux_list = flux.tolist()
    tagged_list = flux_tagged.tolist()
    for sweep in range(MASS_WASTING_SPILL_PASSES):
        last = sweep == MASS_WASTING_SPILL_PASSES - 1
        spilled = [0.0] * n
        spilled_tagged = [0.0] * n
        for i in order:
            f = flux_list[i]
            if f <= 0.0:
                continue
            ft = tagged_list[i]
            if is_ocean_list[i]:
                arrival[i] += f
                arrival_tagged[i] += ft
                continue
            if glaciated_list[i]:
                ice[i] += f
                ice_tagged[i] += ft
                continue
            j = target_list[i]
            if j >= 0:
                dropped = min(f * settle_list[i], room[i])
            elif last or spill_list[i] < 0:
                dropped = f
            else:
                dropped = min(f, room[i])
            if dropped > 0.0:
                dropped_tagged = ft * (dropped / f)
                deposit_list[i] += dropped
                deposit_tagged_list[i] += dropped_tagged
                room[i] -= dropped
                f -= dropped
                ft -= dropped_tagged
            if f <= 0.0:
                continue
            if j >= 0:
                flux_list[j] += f
                tagged_list[j] += ft
            else:
                spilled[spill_list[i]] += f
                spilled_tagged[spill_list[i]] += ft
        flux_list = spilled
        tagged_list = spilled_tagged
        if not any(spilled):
            break
    deposit += np.array(deposit_list)
    deposit_tagged += np.array(deposit_tagged_list)
    return deposit, arrival, ice, deposit_tagged, arrival_tagged, ice_tagged


def _route_wind_deposit(
    points: np.ndarray, wind_u: np.ndarray, wind_v: np.ndarray, source_amount: np.ndarray, world: "World | None" = None
) -> np.ndarray:
    """Single-hop aeolian transport: every node with `source_amount > 0` moves that amount to
    whichever real node sits nearest a point WIND_TRANSPORT_DISTANCE_RAD further downwind on
    the sphere (the exact small-circle geodesic step `point*cos(d) + tangent*sin(d)`, `tangent`
    resolved from (wind_u, wind_v) via `geometry.local_tangent_frame_batch`'s own (east, north)
    convention -- the same convention climate.py's wind field itself uses). A node with
    negligible wind (`speed` near 0) has no real direction to carry it, so it just redeposits
    in place. Exactly conserves `source_amount`'s total -- every unit that leaves a source node
    lands on exactly one target node, via `np.add.at`. `world`, when passed, shares the source
    tree via plates.cached_node_position_tree instead of rebuilding it (see that function's
    own docstring) -- this is the *source* tree only; `target_points` still needs its own query
    since it's a different point set (each source node projected downwind), not this step's
    shared node cloud."""
    n = len(points)
    result = np.zeros(n)
    active = source_amount > 0
    if not np.any(active):
        return result

    east, north = geometry.local_tangent_frame_batch(points)
    speed = np.hypot(wind_u, wind_v)
    has_wind = speed > 1e-9
    safe_speed = np.where(has_wind, speed, 1.0)
    direction = east * np.where(has_wind, wind_u / safe_speed, 0.0)[:, None] + north * np.where(has_wind, wind_v / safe_speed, 0.0)[:, None]

    target_points = np.where(
        has_wind[:, None],
        points * np.cos(WIND_TRANSPORT_DISTANCE_RAD) + direction * np.sin(WIND_TRANSPORT_DISTANCE_RAD),
        points,
    )

    tree = cached_node_position_tree(world, points)
    active_idx = np.nonzero(active)[0]
    _, nearest = tree.query(target_points[active_idx], k=1, workers=query_workers(len(active_idx)))
    np.add.at(result, nearest, source_amount[active_idx])
    return result


def _spread_beach_sediment(points: np.ndarray, elevation: np.ndarray, is_ocean: np.ndarray, ocean_terminal_deposit: np.ndarray) -> np.ndarray:
    """Redistributes whatever `route_downstream` piled onto a single ocean node (the river
    mouth it happened to route to) across nearby shallow water instead -- BEACH_DEPOSITION_
    LOCAL_FRACTION stays right at that node (the delta/mouth itself), the rest spreads,
    inverse-distance weighted, across up to BEACH_SPREAD_NEIGHBOR_COUNT of the nearest shallow
    (shallower than BEACH_SHELF_DEPTH_M) ocean nodes within BEACH_DEPOSITION_RANGE_RAD. A
    source node with no shallow neighbour in range (an open, deep-water coastline) keeps its
    full amount locally rather than losing the difference -- see module comment. Exactly
    conserves `ocean_terminal_deposit`'s total."""
    n = len(points)
    result = np.zeros(n)
    source_idx = np.nonzero(ocean_terminal_deposit > 0)[0]
    if len(source_idx) == 0:
        return result

    shallow_idx = np.nonzero(is_ocean & (elevation > BEACH_SHELF_DEPTH_M))[0]
    if len(shallow_idx) == 0:
        result[source_idx] = ocean_terminal_deposit[source_idx]
        return result

    tree = cKDTree(points[shallow_idx], balanced_tree=False, compact_nodes=False)
    k = min(BEACH_SPREAD_NEIGHBOR_COUNT, len(shallow_idx))
    dist, nearby = tree.query(points[source_idx], k=k, workers=query_workers(len(source_idx)))
    if k == 1:
        dist = dist[:, None]
        nearby = nearby[:, None]

    within_range = dist <= BEACH_DEPOSITION_RANGE_RAD
    weight = np.where(within_range, 1.0 / np.maximum(dist, 1e-9), 0.0)
    weight_sum = weight.sum(axis=1)
    has_shallow_neighbour = weight_sum > 0

    total = ocean_terminal_deposit[source_idx]
    local_share = np.where(has_shallow_neighbour, total * BEACH_DEPOSITION_LOCAL_FRACTION, total)
    spread_share = total - local_share
    result[source_idx] += local_share

    normalized_weight = np.divide(weight, weight_sum[:, None], out=np.zeros_like(weight), where=weight_sum[:, None] > 0)
    contribution = normalized_weight * spread_share[:, None]
    target_idx = shallow_idx[nearby]
    np.add.at(result, target_idx.ravel(), contribution.ravel())
    return result


def _spread_marine_sediment(
    points: np.ndarray, elevation: np.ndarray, is_ocean: np.ndarray, source_amount: np.ndarray
) -> np.ndarray:
    """Redistributes submarine + coastal erosion (see SUBMARINE_EROSION_* / COASTAL_EROSION_*)
    onto the sea floor around each source node -- underwater currents and slope failure carry it
    a short distance downslope, not along the water-flow graph the river-routed pool uses. A
    submerged source keeps MARINE_SEDIMENT_LOCAL_FRACTION locally; a subaerial sea-cliff source
    (is_ocean False) sheds the whole amount seaward. The remainder spreads, inverse-distance
    weighted, across up to MARINE_SPREAD_NEIGHBOR_COUNT of the nearest ocean nodes *lower* than
    the source within MARINE_SPREAD_RANGE_RAD (sediment runs downhill, it doesn't climb the far
    flank of the range it came off). A source with no lower ocean node in range keeps its full
    amount locally rather than losing it. Exactly conserves source_amount's total, via
    np.add.at -- same shape as _spread_beach_sediment."""
    n = len(points)
    result = np.zeros(n)
    source_idx = np.nonzero(source_amount > 0)[0]
    if len(source_idx) == 0:
        return result

    ocean_idx = np.nonzero(is_ocean)[0]
    if len(ocean_idx) == 0:
        result[source_idx] = source_amount[source_idx]
        return result

    tree = cKDTree(points[ocean_idx], balanced_tree=False, compact_nodes=False)
    k = min(MARINE_SPREAD_NEIGHBOR_COUNT, len(ocean_idx))
    dist, nearby = tree.query(points[source_idx], k=k, workers=query_workers(len(source_idx)))
    if k == 1:
        dist = dist[:, None]
        nearby = nearby[:, None]

    target_idx = ocean_idx[nearby]
    is_lower = elevation[target_idx] < elevation[source_idx][:, None]
    within_range = dist <= MARINE_SPREAD_RANGE_RAD
    weight = np.where(is_lower & within_range, 1.0 / np.maximum(dist, 1e-9), 0.0)
    weight_sum = weight.sum(axis=1)
    has_target = weight_sum > 0

    source_is_ocean = is_ocean[source_idx]
    total = source_amount[source_idx]
    local_share = np.where(
        has_target, np.where(source_is_ocean, total * MARINE_SEDIMENT_LOCAL_FRACTION, 0.0), total
    )
    spread_share = total - local_share
    result[source_idx] += local_share

    normalized_weight = np.divide(weight, weight_sum[:, None], out=np.zeros_like(weight), where=weight_sum[:, None] > 0)
    contribution = normalized_weight * spread_share[:, None]
    np.add.at(result, target_idx.ravel(), contribution.ravel())
    return result


def _spread_lake_sediment(lake_depth: np.ndarray, neighbor_idx: np.ndarray, source_amount: np.ndarray) -> np.ndarray:
    """Redistributes whatever `source_amount` (route_downstream's own per-node deposit, from the
    ordinary water-routed sediment pool) landed on a lake/sea member node across that lake's
    whole currently-flooded extent -- see LAKE_SEDIMENT_UNIFORM_FRACTION's own comment for why. A
    source node not currently part of a real lake (`lake_depth` at or below
    `hydrology.LAKE_MIN_VISIBLE_DEPTH_M` -- an ordinary dry floodplain sink, or a closed basin
    that hasn't filled with standing water yet) keeps its full amount right there, same as
    route_downstream's own un-lake-aware dead-end rule.

    Grouped by `hydrology.lake_components` (the same connected-body grouping hydrology.py's own
    breach-erosion term uses) so two separate catchments that happen to have merged into one
    lake this step spread their combined sediment across the *whole* merged body, not just each
    catchment's own original half. A component with zero total lake_depth across every member
    (shouldn't happen for a node `hydrology.lake_components` itself classified as lake, but
    guarded the same defensive way every other weighted spread here is) falls back to a plain
    uniform split. Exactly conserves `source_amount`'s total, via np.add.at, same as
    _spread_beach_sediment/_spread_marine_sediment above."""
    n = len(lake_depth)
    result = np.zeros(n)
    source_idx = np.nonzero(source_amount > 0)[0]
    if len(source_idx) == 0:
        return result

    is_lake = lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    off_lake = source_idx[~is_lake[source_idx]]
    result[off_lake] = source_amount[off_lake]
    if not np.any(is_lake[source_idx]):
        return result

    for members in hydrology.lake_components(is_lake, neighbor_idx):
        total = source_amount[members].sum()
        if total <= 0.0:
            continue
        result[members] += total * _lake_member_weights(lake_depth[members])
    return result


def _lake_member_weights(member_depth: np.ndarray) -> np.ndarray:
    """Each lake member's share of its lake's sediment: LAKE_SEDIMENT_UNIFORM_FRACTION spread
    evenly, the rest by depth (evenly too if the lake has no depth at all)."""
    uniform_share = np.full(len(member_depth), 1.0 / len(member_depth))
    depth_sum = member_depth.sum()
    depth_share = (member_depth / depth_sum) if depth_sum > 0.0 else uniform_share
    return LAKE_SEDIMENT_UNIFORM_FRACTION * uniform_share + (1.0 - LAKE_SEDIMENT_UNIFORM_FRACTION) * depth_share


def _spread_lake_sediment_capped(
    lake_depth: np.ndarray,
    neighbor_idx: np.ndarray,
    source_amount: np.ndarray,
    tagged_amount: np.ndarray,
    capacity: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """`_spread_lake_sediment` for a deposit that already respects per-node `capacity` (landslide
    debris, see `_route_mass_wasting`), keeping it that way: a lake member takes its weighted
    share only up to its capacity, and the remainder goes to the members with room left, by the
    same weights. Whatever the lake as a whole has no room for stays where it was deposited.
    The lake's load is mixed, so every member's deposit carries the lake's tagged fraction.
    Returns (spread amount, spread tagged amount); conserves both totals exactly."""
    result = np.where(source_amount > 0, source_amount, 0.0)
    tagged = np.where(source_amount > 0, tagged_amount, 0.0)
    is_lake = lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    if not np.any(is_lake & (source_amount > 0)):
        return result, tagged

    for members in hydrology.lake_components(is_lake, neighbor_idx):
        total = source_amount[members].sum()
        if total <= 0.0:
            continue
        tagged_fraction = tagged_amount[members].sum() / total
        weight = _lake_member_weights(lake_depth[members])
        room = np.clip(capacity[members], 0.0, None)
        placed = np.zeros(len(members))
        remaining = total
        # Each pass fills at least one member to capacity or places everything, so this ends.
        for _ in range(len(members)):
            open_weight = np.where(placed < room, weight, 0.0)
            if remaining <= 0.0 or open_weight.sum() <= 0.0:
                break
            give = np.minimum(remaining * open_weight / open_weight.sum(), room - placed)
            placed += give
            remaining -= give.sum()
        if remaining > 0.0:
            # No room left in the lake: the unplaced share stays on the nodes it landed on.
            placed += remaining * source_amount[members] / total
        result[members] = placed
        tagged[members] = placed * tagged_fraction
    return result, tagged


@dataclass
class OverflowCarry:
    """What `_carry_overflow` did with the load full receivers turned away. Load matrices have
    one volume column per OVERFLOW_PATHWAYS entry, then the continental tracer's share of it."""

    placed: np.ndarray  # (n, C) settled at each onward receiver
    terminal: np.ndarray  # (C,) no receiver in reach: overloaded-root delamination
    scour_m: np.ndarray  # (n,) rock the debris-laden ice scoured off each node it crossed
    scour_tagged_m: np.ndarray  # (n,) continental share of scour_m
    # Volume the transport left with nowhere to go, by why (before the nearest-receiver search):
    # "pit" (land or ice with no way on), "closed_lake", "sea_floor" (nothing lower with room in
    # range), "hop_limit".
    stuck_m3: dict[str, float] = field(default_factory=lambda: dict.fromkeys(OVERFLOW_STUCK_REASONS, 0.0))
    # How that stuck load was placed: by filling outward over the basin's saddle, or by the
    # nearest-receiver search (the rest is `terminal`).
    basin_fill_m3: float = 0.0
    search_m3: float = 0.0
    # (n,) the part of scour_m that was mobile cover (taken before substrate); none if omitted.
    scour_cover_m: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.scour_cover_m is None:
            self.scour_cover_m = np.zeros_like(self.scour_m)


def _capped_fill(
    amount: np.ndarray, candidates: np.ndarray, weight: np.ndarray, room: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Hand each source's `amount` to its (m, k) `candidates` (-1 = none) in proportion to
    `weight`, never past a receiver's `room` (a volume, decremented in place): OVERFLOW_FILL_
    PASSES passes, each scaling back the gives to an oversubscribed receiver and offering the
    rest to the receivers still open -- the capped water-fill `_spread_coastal_leveling` uses.
    Returns (given (m, k), what each source could not place)."""
    valid = candidates >= 0
    safe = np.where(valid, candidates, 0)
    base = np.where(valid, weight, 0.0)
    given = np.zeros(candidates.shape)
    left = np.array(amount, dtype=float)
    for _ in range(OVERFLOW_FILL_PASSES):
        if not np.any(left > 0.0):
            break
        open_weight = np.where(room[safe] > 0.0, base, 0.0)
        weight_sum = open_weight.sum(axis=1)
        want = np.divide(open_weight, weight_sum[:, None], out=np.zeros_like(open_weight), where=weight_sum[:, None] > 0) * left[:, None]
        tentative = np.zeros(len(room))
        np.add.at(tentative, safe.ravel(), want.ravel())
        scale = np.divide(room, tentative, out=np.ones_like(tentative), where=tentative > room)
        give = want * scale[safe]
        np.add.at(room, safe.ravel(), -give.ravel())
        np.maximum(room, 0.0, out=room)  # round-off
        given += give
        left = np.clip(left - give.sum(axis=1), 0.0, None)
    return given, left


def _carry_overflow(
    points: np.ndarray,
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    on_ice: np.ndarray,
    lake_depth: np.ndarray,
    neighbor_idx: np.ndarray,
    spill_target: np.ndarray,
    ice_flow_target: np.ndarray,
    area: np.ndarray,
    room_vol: np.ndarray,
    excess: np.ndarray,
    scour_limit_m: np.ndarray,
    scour_material_m: np.ndarray,
    world: "World | None" = None,
    scour_cover_m: np.ndarray | None = None,
    scour_cover_material_m: np.ndarray | None = None,
) -> OverflowCarry:
    """Carry the load full receivers turned away (`excess`, (n, C): a volume per
    OVERFLOW_PATHWAYS entry, then its continental share) on to receivers with room
    (`room_vol`, per-node volume left under the Hc cap; inf for a node with no Hc column) --
    issue #288, see OVERFLOW_* for the transport rules.

    The load moves in hops, mixed, so every receiver takes the same pathway and tracer shares
    as the load it came from. Each hop, a node holding load:

    - at sea, settles what room it has, and spreads the rest onto up to
      MARINE_SPREAD_NEIGHBOR_COUNT lower ocean nodes within MARINE_SPREAD_RANGE_RAD;
    - under a glacier (`on_ice`), settles nothing and rides `ice_flow_target`, scouring its
      bed by OVERFLOW_ICE_SCOUR_* within `scour_limit_m` (thickness). It takes the mobile
      cover first (`scour_cover_m`, continental share `scour_cover_material_m`, mixed through
      it), then substrate, whose continental tracer (`scour_material_m`) is given up first;
    - in a lake, spreads across the whole lake up to each member's room, by
      `_lake_member_weights`, and leaves over the lake's spill point once the lake is full;
    - on other land, settles what room it has, and passes the rest to its lowest lower
      neighbour, or over its basin rim (`spill_target`) from a pit.

    Load arriving back at a node it already passed goes on by bare rock (downhill, or over the
    rim), and strictly downhill after that, so an ice or spill loop can't trap it.

    Load with nowhere to go after that -- a full basin or lake with no way out, a sea floor with
    nothing lower, or load still moving after OVERFLOW_MAX_HOPS -- fills outward from where it
    stuck in rising elevation order (a priority flood: the basin fills and overflows its lowest
    saddle), settling on the first nodes with room; then whatever is left is offered to the
    nearest nodes with room. Both stay within OVERFLOW_SEARCH_RANGE_KM; what neither can place
    is `terminal`.
    Conserves every column: excess.sum(0) + scour == placed.sum(0) + terminal."""
    n, c = excess.shape
    tag_col = c - 1
    placed = np.zeros((n, c))
    terminal = np.zeros(c)
    scour_m = np.zeros(n)
    scour_tagged_m = np.zeros(n)
    scour_cover_taken = np.zeros(n)
    idx = np.nonzero(excess[:, :tag_col].sum(axis=1) > 0.0)[0]
    if len(idx) == 0:
        return OverflowCarry(placed, terminal, scour_m, scour_tagged_m, scour_cover_m=scour_cover_taken)
    load = excess[idx].copy()
    stuck_m3 = dict.fromkeys(OVERFLOW_STUCK_REASONS, 0.0)
    basin_fill_m3 = 0.0
    search_m3 = 0.0

    room = np.array(room_vol, dtype=float)
    scour_left = np.clip(np.array(scour_limit_m, dtype=float), 0.0, None)
    material_left = np.clip(np.array(scour_material_m, dtype=float), 0.0, None)
    cover_left = np.zeros(n) if scour_cover_m is None else np.clip(np.array(scour_cover_m, dtype=float), 0.0, None)
    cover_material = np.zeros(n) if scour_cover_material_m is None else np.asarray(scour_cover_material_m, dtype=float)
    cover_fraction = np.clip(np.divide(cover_material, cover_left, out=np.zeros(n), where=cover_left > 0.0), 0.0, 1.0)
    rows = np.arange(n)
    neighbor_elevation = elevation[neighbor_idx]
    lowest = np.argmin(neighbor_elevation, axis=1)
    descent = np.where(neighbor_elevation[rows, lowest] < elevation, neighbor_idx[rows, lowest], -1)
    land_next = np.where(descent >= 0, descent, spill_target)
    riding = on_ice & ~is_ocean
    next_hop = np.where(riding & (ice_flow_target >= 0), ice_flow_target, land_next)
    is_lake = (lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M) & ~is_ocean & ~riding

    def place(at: np.ndarray, volume: np.ndarray, share: np.ndarray) -> None:
        np.add.at(placed, at, volume[:, None] * share)

    def place_given(candidates: np.ndarray, given: np.ndarray, share: np.ndarray) -> None:
        hit = (candidates >= 0) & (given > 0.0)
        place(candidates[hit], given[hit], np.broadcast_to(share[:, None, :], (*candidates.shape, c))[hit])

    lake_of: np.ndarray | None = None
    lakes: list[np.ndarray] = []
    ocean_idx: np.ndarray | None = None
    ocean_tree: cKDTree | None = None
    moved_idx: list[np.ndarray] = []
    moved_load: list[np.ndarray] = []
    stuck_idx: list[np.ndarray] = []
    stuck_load: list[np.ndarray] = []

    def get_stuck(at: np.ndarray, rest: np.ndarray, reason: str) -> None:
        stuck_idx.append(at)
        stuck_load.append(rest)
        stuck_m3[reason] += float(rest[:, :tag_col].sum())

    # Ice and spill edges can point uphill and close a loop (see hydrology._routing_order), so
    # a node the load comes back to sends it on by bare rock instead -- downhill or over its
    # basin rim -- and from a third visit on strictly downhill, which can't loop: debris
    # circling on an ice cap still drains off it by gravity.
    visits = np.zeros(n, dtype=np.int64)

    def move_on(at: np.ndarray, rest: np.ndarray) -> None:
        seen = visits[at]
        target = np.where(seen == 0, next_hop[at], np.where(seen == 1, land_next[at], descent[at]))
        visits[at] += 1
        moved_idx.append(target[target >= 0])
        moved_load.append(rest[target >= 0])
        get_stuck(at[target < 0], rest[target < 0], "pit")

    def settle_here(at: np.ndarray, part: np.ndarray, volume: np.ndarray) -> np.ndarray:
        take = np.minimum(volume, room[at])
        room[at] -= take
        share = part / volume[:, None]
        place(at, take, share)
        return (volume - take)[:, None] * share

    for _ in range(OVERFLOW_MAX_HOPS):
        if len(idx) == 0:
            break
        volume = load[:, :tag_col].sum(axis=1)
        live = volume > 0.0
        idx, load, volume = idx[live], load[live], volume[live]
        moved_idx, moved_load = [], []

        sea = is_ocean[idx]
        if np.any(sea):
            at = idx[sea]
            rest = settle_here(at, load[sea], volume[sea])
            rest_volume = rest[:, :tag_col].sum(axis=1)
            go = rest_volume > 0.0
            if np.any(go):
                if ocean_tree is None:
                    ocean_idx = np.nonzero(is_ocean)[0]
                    ocean_tree = cKDTree(points[ocean_idx], balanced_tree=False, compact_nodes=False)
                at, rest, rest_volume = at[go], rest[go], rest_volume[go]
                k = min(MARINE_SPREAD_NEIGHBOR_COUNT + 1, len(ocean_idx))
                dist, near = ocean_tree.query(points[at], k=k, workers=query_workers(len(at)))
                dist, near = dist.reshape(len(at), k), near.reshape(len(at), k)
                candidates = ocean_idx[near]
                downslope = (dist <= MARINE_SPREAD_RANGE_RAD) & (elevation[candidates] < elevation[at][:, None])
                candidates = np.where(downslope, candidates, -1)
                share = rest / rest_volume[:, None]
                given, left = _capped_fill(rest_volume, candidates, 1.0 / np.maximum(dist, 1e-9), room)
                place_given(candidates, given, share)
                get_stuck(at[left > 0.0], left[left > 0.0][:, None] * share[left > 0.0], "sea_floor")

        ice = riding[idx]
        if np.any(ice):
            at, rest = idx[ice], load[ice].copy()
            scour = np.minimum(OVERFLOW_ICE_SCOUR_FRACTION * volume[ice] / area[at], scour_left[at])
            from_cover = np.minimum(scour, cover_left[at])
            from_substrate = np.minimum(scour - from_cover, material_left[at])
            scour_tagged = from_cover * cover_fraction[at] + from_substrate
            scour_left[at] -= scour
            cover_left[at] -= from_cover
            material_left[at] -= from_substrate
            scour_cover_taken[at] += from_cover
            scour_m[at] += scour
            scour_tagged_m[at] += scour_tagged
            room[at] += scour * area[at]
            rest[:, OVERFLOW_GLACIAL] += scour * area[at]
            rest[:, tag_col] += scour_tagged * area[at]
            move_on(at, rest)

        lake = is_lake[idx]
        if np.any(lake):
            if lake_of is None:
                lakes = hydrology.lake_components(is_lake, neighbor_idx)
                lake_of = np.full(n, -1)
                for component, members in enumerate(lakes):
                    lake_of[members] = component
            at, arriving = idx[lake], load[lake]
            component_of = lake_of[at]
            for component in np.unique(component_of):
                members = lakes[component]
                total = arriving[component_of == component].sum(axis=0)
                total_volume = total[:tag_col].sum()
                share = total / total_volume
                weight = _lake_member_weights(lake_depth[members])
                member_room = np.clip(room[members], 0.0, None)
                taken = np.zeros(len(members))
                remaining = total_volume
                # Each pass fills at least one member to its room or places everything.
                for _ in range(len(members)):
                    open_weight = np.where(taken < member_room, weight, 0.0)
                    if remaining <= 0.0 or open_weight.sum() <= 0.0:
                        break
                    give = np.minimum(remaining * open_weight / open_weight.sum(), member_room - taken)
                    taken += give
                    remaining -= give.sum()
                room[members] -= taken
                placed[members] += taken[:, None] * share
                if remaining <= 0.0:
                    continue
                outlet = spill_target[members]
                outlet = outlet[(outlet >= 0) & (lake_of[np.maximum(outlet, 0)] != component)]
                if len(outlet):
                    moved_idx.append(outlet[:1])
                    moved_load.append(remaining * share[None, :])
                else:
                    get_stuck(at[component_of == component][:1], remaining * share[None, :], "closed_lake")

        land = ~sea & ~ice & ~lake
        if np.any(land):
            at = idx[land]
            rest = settle_here(at, load[land], volume[land])
            go = rest[:, :tag_col].sum(axis=1) > 0.0
            move_on(at[go], rest[go])

        if not moved_idx:
            idx = np.zeros(0, dtype=int)
            break
        all_idx = np.concatenate(moved_idx)
        idx, inverse = np.unique(all_idx, return_inverse=True)
        load = np.zeros((len(idx), c))
        np.add.at(load, inverse, np.concatenate(moved_load))
    if len(idx):
        get_stuck(idx, load, "hop_limit")

    if stuck_idx:
        at, inverse = np.unique(np.concatenate(stuck_idx), return_inverse=True)
        stuck = np.zeros((len(at), c))
        np.add.at(stuck, inverse, np.concatenate(stuck_load))
        stuck_volume = stuck[:, :tag_col].sum(axis=1)
        live = stuck_volume > 0.0
        at, stuck, stuck_volume = at[live], stuck[live], stuck_volume[live]
        share = stuck / stuck_volume[:, None]
        # Fill the full basin and overflow its lowest saddle: a priority flood out from where
        # the load stuck, in rising elevation order, settling on the first nodes with room it
        # reaches (never under the ice), within OVERFLOW_SEARCH_RANGE_KM.
        neighbor_lists = neighbor_idx.tolist()
        elevation_list = elevation.tolist()
        can_settle = (~riding).tolist()
        min_dot = float(np.cos(OVERFLOW_SEARCH_RANGE_RAD))
        for row, start in enumerate(at.tolist()):
            amount = float(stuck_volume[row])
            heap = [(elevation_list[start], start)]
            reached = {start}
            to: list[int] = []
            took: list[float] = []
            while heap and amount > 0.0:
                _, i = heapq.heappop(heap)
                if can_settle[i] and room[i] > 0.0:
                    take = min(float(room[i]), amount)
                    room[i] -= take
                    amount -= take
                    to.append(i)
                    took.append(take)
                for j in neighbor_lists[i]:
                    if j not in reached:
                        reached.add(j)
                        if float(points[j] @ points[start]) >= min_dot:
                            heapq.heappush(heap, (elevation_list[j], j))
            if to:
                place(np.array(to), np.array(took), np.broadcast_to(share[row], (len(to), c)))
                basin_fill_m3 += float(stuck_volume[row]) - amount
            stuck_volume[row] = amount
        stuck = stuck_volume[:, None] * share
        live = stuck_volume > 0.0
        at, stuck, stuck_volume = at[live], stuck[live], stuck_volume[live]
        open_idx = np.nonzero(room > 0.0)[0]
        if len(at) and len(open_idx) == 0:
            terminal += stuck.sum(axis=0)
        elif len(at):
            share = stuck / stuck_volume[:, None]
            # Only nodes that still have room are worth finding: inside a saturated plateau the
            # nearest nodes of any kind are all full. This also reaches room under the ice.
            tree = cKDTree(points[open_idx], balanced_tree=False, compact_nodes=False)
            k = min(OVERFLOW_SEARCH_NEIGHBOR_COUNT, len(open_idx))
            dist, near = tree.query(points[at], k=k, workers=query_workers(len(at)))
            dist, near = dist.reshape(len(at), k), open_idx[near.reshape(len(at), k)]
            candidates = np.where(dist <= OVERFLOW_SEARCH_RANGE_RAD, near, -1)
            given, left = _capped_fill(stuck_volume, candidates, 1.0 / np.maximum(dist, 1e-9), room)
            place_given(candidates, given, share)
            search_m3 += float(given.sum())
            terminal += (left[:, None] * share).sum(axis=0)
    return OverflowCarry(placed, terminal, scour_m, scour_tagged_m, stuck_m3, basin_fill_m3, search_m3, scour_cover_taken)


def _coastal_openness(points: np.ndarray, is_ocean: np.ndarray) -> np.ndarray:
    """Per-node "wave exposure" proxy in [0, 1] -- the fraction of nodes within
    COASTAL_OPENNESS_RANGE_RAD that are open ocean (`is_ocean` here is hydrology's
    connectivity-aware mask, so an enclosed interior pit / inland lake counts as *not* open,
    and a node ringed by such water still reads as sheltered). ~0.5 on a straight open coast,
    well below that inside a bay or behind a barrier, above it on a headland or peninsula
    tip. Two radius counts (`query_ball_point(..., return_length=True)`), density-independent
    -- not a fixed-k neighbourhood, whose real radius would shrink as node_density rises. 0
    everywhere for a world too small to have a meaningful neighbourhood.

    The second count is of each node's *non*-open neighbours (land + connectivity-enclosed
    water), not its open-ocean ones, with `ocean_count` recovered as `total - non_open_count`:
    every node is in exactly one of the two sets (self included), so the subtraction is exact
    and the returned field is bit-identical, but on a mostly-ocean world the non-open set is
    several times smaller, making both that tree build and its radius query much cheaper."""
    n = len(points)
    if n <= SLOPE_NEIGHBOR_COUNT:
        return np.zeros(n)
    workers = query_workers(n)
    tree = cKDTree(points, balanced_tree=False, compact_nodes=False)
    total = tree.query_ball_point(points, COASTAL_OPENNESS_RANGE_RAD, workers=workers, return_length=True)
    if not np.any(is_ocean):
        return np.zeros(n)
    non_open_idx = np.nonzero(~np.asarray(is_ocean, dtype=bool))[0]
    if len(non_open_idx) == 0:
        return total / np.maximum(total, 1)
    non_open_tree = cKDTree(points[non_open_idx], balanced_tree=False, compact_nodes=False)
    non_open_count = non_open_tree.query_ball_point(points, COASTAL_OPENNESS_RANGE_RAD, workers=workers, return_length=True)
    ocean_count = total - non_open_count
    return ocean_count / np.maximum(total, 1)


def leveling_datum_m(
    sea_level_m: float,
    openness: np.ndarray,
    dist_to_land: np.ndarray,
) -> np.ndarray:
    """Each node's local target elevation for the symmetric coastal-leveling pass, a single
    continuous function of wave exposure: `sea_level - LEVELING_PLATFORM_UNDERCUT_M * exposure
    + LEVELING_MARSH_CREST_M * shelter`. An exposed node (openness >= LEVELING_EXPOSURE_REF)
    targets a wave-cut platform a few metres *below* the waterline (real shore platforms are
    planed to roughly low-tide level and a little below); a sheltered node (openness <
    INFILL_SHELTER_REF) targets a marsh crest a few metres *above* it; a straight open coast
    lands near sea level. Cutting/filling to a datum that is *off* the waterline (rather than
    balanced exactly on it) is what lets the pass *resolve* a checkerboard into a coast that
    follows the shelter field, not just damp its amplitude. A barrier candidate (land within
    BARRIER_LANDWARD_RAD but still facing open water, openness >= BARRIER_MIN_OPENNESS)
    instead targets `sea_level + BARRIER_CREST_M`, so a shore-parallel bar can just breach the
    surface. Both the grind side and the fill side aim at this one number, so they can never
    fight each other over a node."""
    exposure = np.clip(openness / LEVELING_EXPOSURE_REF, 0.0, 1.0)
    shelter = np.clip(1.0 - openness / INFILL_SHELTER_REF, 0.0, 1.0)
    datum = sea_level_m - LEVELING_PLATFORM_UNDERCUT_M * exposure + LEVELING_MARSH_CREST_M * shelter
    is_barrier = (dist_to_land <= BARRIER_LANDWARD_RAD) & (openness >= BARRIER_MIN_OPENNESS)
    return np.where(is_barrier, sea_level_m + BARRIER_CREST_M, datum)


def coastal_leveling_grind(
    elevation: np.ndarray,
    sea_level_m: float,
    openness: np.ndarray,
    datum_m: np.ndarray,
    dt_myr: float,
    local_relief_m: np.ndarray | None = None,
) -> np.ndarray:
    """Meters ground off each band node this step (the grind / source half of the symmetric
    leveling pass): a node within COASTAL_LEVELING_BAND_M of sea level that stands *above* its
    local `datum_m` (see `leveling_datum_m`) is planed down toward it at `LEVELING_RATE_M_PER_MYR
    * coastal * proximity * prominence` -- `coastal` saturating once LEVELING_EXPOSURE_REF of
    the neighbourhood is open water (a landlocked lowland, openness ~ 0, is untouched),
    `proximity` tapering linearly to 0 at the band edge, and `prominence` (from
    `local_relief_m`, this node's height above its own neighbourhood mean -- 1.0 when omitted)
    making a protruding speck plane several times faster than a flat sheet and a hollow one
    barely at all. Applies to a just-submerged shoal as readily as to low land -- the old
    planation gate ignored everything below sea level. Returned uncapped against the neighbour
    drop -- a wave-cut platform genuinely planes a bench below adjacent terrain -- but never
    past `datum_m`; apply_erosion caps it further against the erosion already taken this
    step."""
    height_above_sea_m = elevation - sea_level_m
    in_band = np.abs(height_above_sea_m) <= COASTAL_LEVELING_BAND_M
    proximity = np.clip(1.0 - np.abs(height_above_sea_m) / COASTAL_LEVELING_BAND_M, 0.0, 1.0)
    coastal = np.clip(openness / LEVELING_EXPOSURE_REF, 0.0, 1.0)
    if local_relief_m is None:
        prominence = 1.0
    else:
        prominence = np.clip(1.0 + local_relief_m / PROMINENCE_REF_M, 0.0, PROMINENCE_MAX)
    rate_m = LEVELING_RATE_M_PER_MYR * coastal * proximity * prominence * dt_myr
    room_to_datum_m = np.clip(elevation - datum_m, 0.0, None)
    return np.where(in_band, np.minimum(rate_m, room_to_datum_m), 0.0)


def _spread_coastal_leveling(
    points: np.ndarray,
    elevation: np.ndarray,
    openness: np.ndarray,
    dist_to_land: np.ndarray,
    sea_level_m: float,
    datum_m: np.ndarray,
    source_amount: np.ndarray,
    dt_myr: float,
    local_relief_m: np.ndarray | None = None,
    area_m2: np.ndarray | None = None,
    tagged_amount: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Redistributes `source_amount` (ground-off rock + the redirected share of submarine/
    coastal erosion + the redirected distributary share of the river-deposition lump) onto
    band nodes sitting *below* their local `datum_m` (see `leveling_datum_m`) -- the fill /
    sink half of the symmetric leveling pass, land or ocean alike, so a dry interdistributary
    low just above sea level fills as readily as sheltered shallow water. A sink must have at
    least LEVELING_MIN_OPENNESS of open water in its neighbourhood (an inland-lake shore,
    connectivity-openness ~ 0, never silts up here).

    This is a **capped water-fill**, not a one-shot weighted scatter: no sink can be raised
    past its own datum in a step (`below_datum_m` is a hard per-node capacity), and the fill
    runs LEVELING_FILL_ITERS passes -- each pass distributes what's left of every source
    across its still-open sinks in proportion to `priority * attract * room_left / dist`, then
    caps every sink at its capacity and carries the overflow into the next pass -- so a big
    lump genuinely spreads across the whole flat plain (each near sink takes only its ~30 m of
    room, the rest flows on) instead of piling 200 m onto the single most-weighted node.
    `attract = shelter * hollow` (shelter = how enclosed -- the most sheltered embayment silts
    up fastest, that is where a marsh belongs; hollow = how far below its own neighbourhood
    mean, so an isolated pond fills faster than a broad shelf -- 1.0 when `local_relief_m` is
    omitted). A barrier candidate (see `leveling_datum_m`) gets a BARRIER_PRIORITY boost and a
    BARRIER_ATTRACT_FLOOR. A source that is itself a below-datum node keeps
    LEVELING_LOCAL_FRACTION (its own hollow to fill); a pure source spreads in full. The
    spread never returns to its own source node; whatever no sink had room for after the last
    pass stays on the source. Exactly conserves source_amount's total via np.add.at.

    `area_m2` (issue #275), when passed, makes `source_amount` a *volume* per node rather than
    a thickness: each sink's per-step capacity becomes its thickness capacity times its own
    area, so the result is a volume too and conserves area-weighted material on unequal
    cells. `tagged_amount` is a sub-share of `source_amount` (continental-derived material,
    never more than it) carried along exactly the same transfers; returns `(result,
    tagged_result)`, the second all-zero when `tagged_amount` is omitted."""
    n = len(points)
    result = np.zeros(n)
    tagged_result = np.zeros(n)
    area = np.ones(n) if area_m2 is None else area_m2
    tagged = np.zeros(n) if tagged_amount is None else tagged_amount
    tag_fraction = np.divide(tagged, source_amount, out=np.zeros(n), where=source_amount > 0)
    source_idx = np.nonzero(source_amount > 0)[0]
    if len(source_idx) == 0:
        return result, tagged_result

    height_above_sea_m = elevation - sea_level_m
    below_datum_m = datum_m - elevation
    in_band = np.abs(height_above_sea_m) <= COASTAL_LEVELING_BAND_M
    is_sink = in_band & (below_datum_m > 0.0) & (openness > LEVELING_MIN_OPENNESS)
    sink_idx = np.nonzero(is_sink)[0]
    if len(sink_idx) == 0:
        result[source_idx] = source_amount[source_idx]
        return result, tagged.copy()

    # Only sources within INFILL_RANGE_RAD of some sink can place anything. On a mostly-ocean
    # world `source_idx` is ~every ocean node (each carries a sliver of redirected submarine/
    # coastal spoil) while the sink band is a thin coastal ribbon, so this usually cuts the
    # k-NN query below from ~10^5 source points to a few thousand. Every excluded source keeps
    # its full amount -- bit-identical to the `any_reachable is False` branch that would
    # otherwise handle it once its k nearest sinks all came back beyond the range.
    sink_tree = cKDTree(points[sink_idx], balanced_tree=False, compact_nodes=False)
    near_sink_count = sink_tree.query_ball_point(
        points[source_idx], INFILL_RANGE_RAD, workers=query_workers(len(source_idx)), return_length=True
    )
    out_of_range = source_idx[near_sink_count == 0]
    result[out_of_range] = source_amount[out_of_range]
    tagged_result[out_of_range] = tagged[out_of_range]
    source_idx = source_idx[near_sink_count > 0]
    if len(source_idx) == 0:
        return result, tagged_result

    # Per-step capacity: room to the datum, but never more than the grind side can take off in
    # the same step -- so a checkerboard hollow and its neighbouring bump both move ~one
    # LEVELING_RATE step toward the datum, damping the dither instead of slamming the hollow
    # shut in a single step and setting up a new oscillation.
    # In the units of `source_amount` -- a volume when `area_m2` is passed.
    capacity = np.minimum(below_datum_m[sink_idx], LEVELING_RATE_M_PER_MYR * dt_myr) * area[sink_idx]
    shelter = np.clip(1.0 - openness[sink_idx] / INFILL_SHELTER_REF, 0.0, 1.0)
    if local_relief_m is None:
        hollow = 1.0
    else:
        hollow = np.clip(1.0 - local_relief_m[sink_idx] / PROMINENCE_REF_M, 0.0, PROMINENCE_MAX)
    is_barrier = (dist_to_land[sink_idx] <= BARRIER_LANDWARD_RAD) & (openness[sink_idx] >= BARRIER_MIN_OPENNESS)
    base_attract = shelter * hollow
    attract = np.where(is_barrier, np.maximum(base_attract, BARRIER_ATTRACT_FLOOR), base_attract)
    priority = np.where(is_barrier, BARRIER_PRIORITY, 1.0)
    sink_pref = priority * attract  # distance- and room-independent part of the weight

    k = min(LEVELING_SPREAD_NEIGHBOR_COUNT, len(sink_idx))
    dist, nearby = sink_tree.query(points[source_idx], k=k, workers=query_workers(len(source_idx)))
    if k == 1:
        dist = dist[:, None]
        nearby = nearby[:, None]

    # A source that is itself a sink must not spread onto itself.
    is_self = sink_idx[nearby] == source_idx[:, None]
    reachable = (dist <= INFILL_RANGE_RAD) & ~is_self
    inv_dist = np.where(reachable, 1.0 / np.maximum(dist, 1e-9), 0.0)
    pref = sink_pref[nearby] * inv_dist  # (n_source, k), distance-weighted preference

    total = source_amount[source_idx]
    source_is_sink = is_sink[source_idx]
    any_reachable = reachable.any(axis=1)
    local_share = np.where(
        any_reachable, np.where(source_is_sink, total * LEVELING_LOCAL_FRACTION, 0.0), total
    )
    amt = total - local_share  # per-source amount still to place
    result[source_idx] += local_share
    source_tag = tag_fraction[source_idx]
    tagged_result[source_idx] += local_share * source_tag
    tagged_received = np.zeros(len(sink_idx))

    received = np.zeros(len(sink_idx))
    for _ in range(LEVELING_FILL_ITERS):
        if amt.sum() <= 1e-6:
            break
        room_left = np.clip(capacity - received, 0.0, None)
        weight = pref * room_left[nearby]
        weight_sum = weight.sum(axis=1)
        frac = np.divide(weight, weight_sum[:, None], out=np.zeros_like(weight), where=weight_sum[:, None] > 0)
        want = frac * amt[:, None]  # (n_source, k) tentative give this pass

        tentative = np.zeros(len(sink_idx))
        np.add.at(tentative, nearby.ravel(), want.ravel())
        # Scale back every source's give to whichever sinks oversubscribed this pass.
        scale = np.where(tentative > room_left, np.divide(room_left, tentative, out=np.ones_like(tentative), where=tentative > 0), 1.0)
        given = want * scale[nearby]
        np.add.at(received, nearby.ravel(), given.ravel())
        np.add.at(tagged_received, nearby.ravel(), (given * source_tag[:, None]).ravel())
        amt = amt - given.sum(axis=1)

    result[source_idx] += amt  # whatever no sink had room for stays on the source
    tagged_result[source_idx] += amt * source_tag
    np.add.at(result, sink_idx, received)
    np.add.at(tagged_result, sink_idx, tagged_received)
    return result, tagged_result


def apply_erosion(
    world: "World",
    years: float,
    node_cloud: tuple[np.ndarray, list[Plate]] | None = None,
) -> ErosionResult | None:
    """Erodes every plate's elevation nodes based on the world's current climate and flow
    routing -- rain/sheet erosion (precipitation x slope), river-channelized erosion
    (accumulated flow x slope, boosted by the node's own established channel), weathering
    (wind speed x humidity), glacier erosion (accumulated ice depth x slope, see
    GLACIER_EROSION_* for how ice's own weight drives this), and seismic erosion
    (earthquake-triggered landsliding, scaled by elevation as a stand-in for how tectonically
    active a range is -- see SEISMIC_EROSION_* constants; its debris runs out by gravity onto
    the foreland, see MASS_WASTING_*) -- then routes the rest of the eroded material downstream, redepositing part of it wherever a big, slow river drops its load (a
    floodplain/delta) instead of losing everything to the coast. Glacially-eroded material
    (net of GLACIER_TILL_FRACTION's own immediate local deposit) is routed separately, along
    the ice's own real flow path rather than water's, settling only once it reaches the
    glacier's actual melting margin -- a terminal moraine/outwash deposit pushed outside the
    ice by the glacier's own flow, see the comment above GLACIER_TILL_FRACTION. Separately
    relaxes elevation under thick ice toward its local neighborhood mean (glacial flattening,
    see `_flatten`). Also grows channel_depth (from this step's river- and glacier-erosion
    terms, plus the scour of overflow-laden ice -- see OVERFLOW_ICE_SCOUR_*) and
    channel_width (from discharge alone -- larger flows carve a wider channel); lake_depth/
    glacier_depth/silt_depth are hydrology.py's own state transitions, read directly from
    World.hydrology_cache. Keeps each column's mobile cover (see MOBILE_COVER_*): removals take
    it first and everything that settles adds to it. All persistent, see plates.ElevationLine.
    Mutates world.plates'
    line elevations in place; never touches node positions or line topology, so this can't
    interact with line regularization or point reassignment at all (both of those are
    purely about node density/position/ownership).

    Always computes climate fresh (never reuses World.climate_cache itself -- this runs
    right after this step's own tectonic/topology changes, so a cache from a previous step
    would already be stale for erosion's own purposes) and stores the result back onto
    World.climate_cache/World.hydrology_cache, so /world/stats and a map render don't each
    also trigger their own recomputation this same turn -- see
    climate.compute_climate_cached.

    Returns this step's ErosionResult (or None for an empty world, mirroring the
    World.hydrology_cache = None branch below) so geology.py's own soil/coal/oil-gas formation
    (called right after this, from world.step_world) can reuse these same per-node terms
    instead of re-deriving them -- see ErosionResult's own docstring.

    `node_cloud`, when passed (world.py's step_world computes it once, right after this
    step's own rotation/boundary evolution/topology changes settle), is an already-gathered
    (points, plates_in_order) pair -- see plates.gather_node_positions -- reused here, in
    climate.compute_climate, and in hydrology.compute_hydrology, all three of which would
    otherwise independently redo the identical world_xyz rotation for every node this same
    step (node positions don't move again until the next step's rotation, so one gather
    upfront is enough for all three)."""
    node_cloud = node_cloud if node_cloud is not None else gather_node_positions(world.plates)
    # Seeds the continental-material tracer on a world that never had one (a hand-built test
    # world, or a save predating it) before this step starts moving it.
    continental_ledger.ensure_initialized(world)
    mobile_cover.ensure_ledger(world)
    fields = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density), node_cloud=node_cloud)
    world.climate_cache = fields

    points, elevation, prior_channel_depth, prior_channel_width, prior_glacier_depth, plates_in_order = _gather_nodes(world, node_cloud=node_cloud)
    n = len(points)
    if n == 0:
        world.hydrology_cache = None
        world.hydrology_cache_step = None
        return None

    height, width = fields.precipitation_mm.shape
    row, col = climate_grid_indices(points, height, width)
    precipitation_mm = fields.precipitation_mm[row, col]
    wind_u_at_nodes = fields.wind_u[row, col]
    wind_v_at_nodes = fields.wind_v[row, col]
    wind_speed = np.hypot(wind_u_at_nodes, wind_v_at_nodes)
    humidity = fields.humidity[row, col]
    # Preliminary, elevation-only ocean test -- just for the ocean-vs-air temperature pick
    # feeding the flow solve below. The authoritative, connectivity-aware mask
    # (`hydro.is_ocean`, which also excludes an enclosed interior sub-sea-level pit) is only
    # available after compute_hydrology; `is_ocean_node` is rebound to it right after.
    is_ocean_node = elevation <= world.sea_level_m
    # The same real temperature a node actually experiences that render_image.py's own
    # temperature view displays -- ocean surface over water, moderated air over land.
    temperature = np.where(is_ocean_node, fields.ocean_temperature_c[row, col], fields.air_temperature_c[row, col])

    slope, drop_to_lowest_neighbor_m = compute_slope(points, elevation, world=world)
    dt_myr = years / 1_000_000.0

    # Seasons aren't stepped through, but glacier melt imputes them: each node's seasonal
    # half-amplitude around its annual mean (see hydrology.GLACIER_MELT_THRESHOLD_C).
    seasonal_amplitude = biomes.seasonal_temp_amplitude(
        fields.lat_deg[row],
        biomes.grid_continentality(fields.is_ocean)[row, col],
        world.axial_tilt_deg,
        relief_m=elevation - world.sea_level_m,
    )
    hydro = hydrology.compute_hydrology(
        world, precipitation_mm, temperature, years, node_cloud=node_cloud, seasonal_amplitude_at_nodes=seasonal_amplitude
    )
    world.hydrology_cache = hydro
    world.hydrology_cache_step = world.steps_taken
    # From here on use hydrology's connectivity-aware mask: an interior pit that dipped below
    # sea level without connecting to open ocean now gets the *subaerial* erosion/deposition
    # pathways (and its lake silt, folded into elevation below), not the marine ones.
    is_ocean_node = hydro.is_ocean
    # Lake merge/split transitions (hydro.lake_events, via lakes.summarize_lake_events) are
    # deliberately *not* routed to World.log_event: even after the near-sea-level aggregation
    # they churn constantly along a dithering coastline and drown out the genuine tectonic
    # events in the console. The structured events are still computed and available for tests.
    water_accum_m = hydro.flow_accum / 1000.0

    # `world.*_erosion_multiplier` (and the deposition/leveling knobs further down) are the
    # user's live geomorphic-budget tuning knobs -- 1.0 everywhere == untuned behaviour (see
    # World's field group). Applied right where each term is formed so every downstream cap/
    # split/isostasy step just sees a scaled amount and stays self-consistent.
    rain = RAIN_EROSION_COEFFICIENT * world.rain_erosion_multiplier * slope * (precipitation_mm / 1000.0) * dt_myr
    channel_boost = 1.0 + CHANNEL_EROSION_BOOST * np.clip(prior_channel_depth / CHANNEL_BOOST_REFERENCE_M, 0.0, 1.0)
    river = (
        RIVER_EROSION_COEFFICIENT
        * world.river_erosion_multiplier
        * channel_boost
        * np.power(np.clip(water_accum_m, 0.0, None), RIVER_FLOW_EXPONENT)
        * np.power(slope, RIVER_SLOPE_EXPONENT)
        * dt_myr
    )
    humidity_norm = np.clip(humidity / HUMIDITY_REFERENCE, 0.0, 1.0)
    relief_factor = np.clip(slope / WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    weathering = WEATHERING_COEFFICIENT * world.wind_erosion_multiplier * wind_speed * humidity_norm * relief_factor * dt_myr
    ice_factor = np.clip(prior_glacier_depth / GLACIER_EROSION_REFERENCE_DEPTH_M, 0.0, GLACIER_EROSION_MAX_FACTOR)
    glacier = GLACIER_EROSION_COEFFICIENT * world.glacier_erosion_multiplier * slope * ice_factor * dt_myr
    # See SEISMIC_EROSION_* constants' own comment: elevation (clipped/normalized against
    # SEISMIC_EROSION_ELEVATION_REFERENCE_M, then raised to a superlinear power) stands in for
    # how tectonically active/seismic a mountain range is, this model having no separate
    # fault/stress field to drive it from directly. np.clip(elevation, 0, None) first since a
    # negative elevation would otherwise flip sign under the exponent -- ocean nodes are zeroed
    # out below regardless, but this keeps the intermediate factor itself well-defined.
    mountain_height_factor = np.clip(np.clip(elevation, 0.0, None) / SEISMIC_EROSION_ELEVATION_REFERENCE_M, 0.0, SEISMIC_EROSION_MAX_HEIGHT_FACTOR)
    seismic = (
        SEISMIC_EROSION_COEFFICIENT
        * world.seismic_erosion_multiplier
        * slope
        * np.power(mountain_height_factor, SEISMIC_EROSION_ELEVATION_EXPONENT)
        * dt_myr
    )
    # Direct earthquake-driven landsliding burst around each recent epicentre (see
    # EARTHQUAKE_EROSION_* and faults.Earthquake). No-op when nothing has ruptured recently.
    seismic = seismic * _earthquake_erosion_multiplier(world, points)
    # Zeroed over ocean nodes (elevation <= sea level, the same convention climate.py/plates.py
    # use everywhere else): every source here is a subaerial process. The sea floor and the
    # shoreline get their own erosion separately, below (submarine + coastal erosion -- see
    # SUBMARINE_EROSION_* / COASTAL_EROSION_*).
    raw_erosion_total = rain + river + weathering + glacier + seismic
    detach_capacity = np.where(is_ocean_node, 0.0, np.clip(raw_erosion_total, 0.0, None))
    # Issue #275: every removal below is also capped at its column's Hc headroom above
    # lithosphere.MIN_CRUSTAL_THICKNESS_M *before* anything is routed. The Hc floor clip at the
    # write-back used to be the only guard, so a column eroded to the floor still handed its
    # full computed load downstream -- rock its source never gave up. Nodes with no column (v1
    # PlateWithLines, Hc 0) keep the bare elevation-only caps.
    prior_hc = collect_all_crustal_thickness(plates_in_order)
    prior_hm = collect_all_mantle_lithosphere_thickness(plates_in_order)
    has_column = prior_hc > 0.0
    removable_m = np.where(has_column, np.clip(prior_hc - lithosphere.MIN_CRUSTAL_THICKNESS_M, 0.0, None), np.inf)
    # Cratons (cratons.py) resist erosional unroofing: at full strength a craton column gives up
    # only (1 - CRATON_EROSION_RESISTANCE) of the substrate that would otherwise be removed.
    # Applied here, at the source, so everything routed downstream stays exactly what the
    # sources gave up.
    craton_keep = 1.0 - cratons.CRATON_EROSION_RESISTANCE * cratons.strength(
        np.concatenate([p.collect("craton_crust_m") for p in plates_in_order])
    )
    # Mobile cover (issue #297 phase 2, see MOBILE_COVER_*) and the continental tracer
    # (`continental_material_m`, see "Volume and provenance" below). A stored tracer or cover
    # above its column's Hc is material some *other* step already removed without updating it
    # (tectonic thinning, not yet instrumented -- see #272). The continental ledger's
    # inventories never count tracer excess, so clipping it is budget-neutral and
    # `stale_tracer_excess_m3` only reports it. The cover ledger's live inventory does count the
    # cover's, so `stale_mobile_cover_excess_m3` is booked into its `clipped_m3` below. The
    # cover's continental share is a share of both.
    stored_material = np.concatenate([p.collect("continental_material_m") for p in plates_in_order])
    prior_material = np.clip(stored_material, 0.0, prior_hc)
    stored_cover = np.concatenate([p.collect("mobile_cover_m") for p in plates_in_order])
    prior_cover = np.clip(stored_cover, 0.0, np.where(has_column, prior_hc, np.inf))
    prior_cover_material = np.clip(
        np.concatenate([p.collect("mobile_cover_continental_m") for p in plates_in_order]),
        0.0,
        np.minimum(prior_cover, prior_material),
    )
    # Capped at the drop to the lowest neighbor so a single step can't erode a node below the
    # valley floor it drains into, and at the Hc headroom. A cap cuts substrate before cover
    # (the cover is on top, so it goes first).
    erosion_cap_m = np.minimum(drop_to_lowest_neighbor_m, removable_m)
    cover_entrained, substrate_detached = _strip_mobile_cover(prior_cover, detach_capacity)
    cover_entrained = np.minimum(cover_entrained, erosion_cap_m)
    erosion_amount = cover_entrained + np.minimum(substrate_detached, erosion_cap_m - cover_entrained) * craton_keep
    cover_m = prior_cover - cover_entrained
    # channel_depth is the terrain's own carved-channel record, so it must track what actually
    # got taken off this point's elevation: river's (and glacier's -- see new_channel_depth)
    # contribution is scaled by applied_scale, what came off over what the substrate laws
    # asked for. The neighbor-drop cap pulls it below 1 -- otherwise a node pinned near its
    # lowest neighbor (a valley floor at grade) would keep "carving" toward
    # MAX_CHANNEL_DEPTH_M while its elevation barely moves, decoupling the two fields. Mobile
    # cover can push it up to MOBILE_COVER_ERODIBILITY_FACTOR: a channel cuts through loose
    # sediment that much faster than through rock.
    applied_scale = np.divide(erosion_amount, raw_erosion_total, out=np.zeros_like(raw_erosion_total), where=raw_erosion_total > 0)
    applied_river = river * applied_scale

    # Split off the sources with genuinely non-water transport (see module comment above
    # GLACIER_TILL_FRACTION/WIND_DEPOSITION_FRACTION): a glacier's own till settles close by
    # (glacier_till, deposited without transport -- same node) while the rest travels with the
    # ice itself (glacier_carried, routed separately below, not through the water pool); wind-
    # eroded material rides the wind rather than water (wind_redeposit_source, carried by
    # _route_wind_deposit below). All fractions are taken from the *applied* (post-neighbor-
    # drop-cap) amount, same applied_scale reasoning as applied_river above, so the split still
    # exactly partitions erosion_amount. Seismic erosion's landslide debris moves by gravity
    # rather than with water (issue #275 phase 4, see MASS_WASTING_*), so it stays out of the
    # water-routed pool and runs out along its own path below.
    applied_weathering = weathering * applied_scale
    applied_glacier = glacier * applied_scale
    applied_seismic = seismic * applied_scale
    glacier_till = applied_glacier * GLACIER_TILL_FRACTION
    glacier_carried = applied_glacier - glacier_till
    wind_redeposit_source = applied_weathering * WIND_DEPOSITION_FRACTION
    weathering_routed = applied_weathering - wind_redeposit_source
    water_routed_amount = (rain * applied_scale) + applied_river + weathering_routed

    # Submarine + coastal erosion: the sea floor's and the shoreline's counterparts to the
    # subaerial sources above (all of which were just zeroed over ocean nodes). Submarine
    # erosion slumps a freshly-uplifted submerged range back down (bottom currents +
    # pressure-driven mass wasting), which is what keeps a range built by two submerged plates
    # colliding growing far slower than a subaerial one; coastal erosion gnaws at the
    # near-sea-level band on both sides of the shoreline (wave attack + frost shattering),
    # including the crest of a mid-ocean range that rises into that band. See the
    # SUBMARINE_EROSION_* / COASTAL_EROSION_* constants and their helper functions. Both draw
    # against whatever drop to the lowest neighbor (and Hc headroom) subaerial erosion didn't
    # already claim, so a single step still can't carve a node below the sea floor / valley it
    # drains into. Their rock sheds seaward onto the surrounding sea floor as marine sediment
    # (_spread_marine_sediment) -- underwater currents and slope failure, not any river's flow
    # graph.
    submarine = submarine_erosion_amount(elevation, slope, is_ocean_node, dt_myr)
    coastal = coastal_erosion_amount(elevation, temperature, dt_myr)
    remaining_drop_m = np.clip(erosion_cap_m - erosion_amount, 0.0, None)
    # ocean_erosion_multiplier scales both the sea-floor slump and the shoreline wave/frost
    # attack together (still capped at the remaining drop to the lowest neighbour). Whatever
    # cover the subaerial sources left goes first, at the cover's higher erodibility.
    sea_cover, sea_substrate = _strip_mobile_cover(cover_m, np.clip((submarine + coastal) * world.ocean_erosion_multiplier, 0.0, None))
    sea_cover = np.minimum(sea_cover, remaining_drop_m)
    sea_side_erosion = sea_cover + np.minimum(sea_substrate * craton_keep, remaining_drop_m - sea_cover)
    cover_entrained = cover_entrained + sea_cover
    cover_m = cover_m - sea_cover

    # Symmetric coastal leveling (see the COASTAL_OPENNESS_* / COASTAL_LEVELING_* / LEVELING_*
    # / INFILL_* / BARRIER_* / PROMINENCE_* constants). `coastal_openness` is this step's
    # wave-exposure proxy; `leveling_datum` is each near-sea-level node's local target
    # elevation. One pass then pushes every band node toward that datum: grind the ones
    # standing above it (`ground_off`, capped here against the erosion already taken this step
    # so land can't be shoved underwater in one step), fill the ones sitting below it
    # (`leveling_fill`). Its fill pool is the ground-off rock + COASTAL_INFILL_MARINE_FRACTION
    # of the submarine/coastal pool + `delta_redirect` (the distributary share of the clumped
    # river-deposition lump); the rest of the sea-side pool still spreads to deep water. All
    # mass-conserving via np.add.at, same as _spread_marine_sediment.
    coastal_openness = _coastal_openness(points, is_ocean_node)
    dist_to_land = world.distance_from_land_approx(points)
    # Each node's height above its own flow-graph neighbourhood mean -- drives the prominence /
    # hollow reweighting that lets the pass actually resolve a checkerboard at node scale
    # (waves plane protrusions, fill hollows; the openness field alone is far too smooth to
    # make a land/ocean call at ~60 km node spacing). Reuses hydrology's k=FLOW_NEIGHBOR_COUNT
    # graph.
    local_relief_m = elevation - elevation[hydro.neighbor_idx].mean(axis=1)
    leveling_datum = leveling_datum_m(world.sea_level_m, coastal_openness, dist_to_land)
    # coastal_leveling_multiplier scales the near-shore planation grind (a prime long-run land
    # drain) -- applied before the one-step "can't be shoved below its datum" safety cap, which
    # then still holds, and feeds through to leveling_source/erosion_amount consistently.
    ground_off = coastal_leveling_grind(elevation, world.sea_level_m, coastal_openness, leveling_datum, dt_myr, local_relief_m)
    ground_off = ground_off * world.coastal_leveling_multiplier * craton_keep
    ground_off = np.minimum(ground_off, np.clip(elevation - erosion_amount - leveling_datum, 0.0, None))
    ground_off = np.minimum(ground_off, np.clip(removable_m - erosion_amount - sea_side_erosion, 0.0, None))

    # Glacial flattening rides the same knob as glacial abrasion -- both are "heavier / more
    # active ice reworks the bed harder". A per-edge send (see _flatten), capped at whatever Hc
    # headroom the removals above left.
    flatten_send = _flatten(
        hydro,
        ice_factor,
        years,
        removable_m=np.clip(removable_m - erosion_amount - sea_side_erosion - ground_off, 0.0, None),
        multiplier=world.glacier_erosion_multiplier,
    )
    # Cratons resist glacial flattening like every other source above: scaling each source
    # node's sends keeps the removal, the transport and the deposits consistent.
    flatten_send = flatten_send * craton_keep[:, None]
    flatten_removed = flatten_send.sum(axis=1)

    # Volume and provenance (issue #275). Every transport below carries *volume*
    # (thickness x the node's own accounting area -- exact cell areas on a quad surface, which
    # is not equal-area) and is converted back to thickness at the receiving node, so moving
    # rock between a large and a small cell no longer creates or destroys material. Alongside
    # each pool travels its continental-derived share (`continental_material_m`, the tracer
    # continental_ledger.py persists): a node gives up continental material first, up to its
    # tracer (it sits on top -- continental sediment draped over an oceanic host is the first
    # thing eroded off it), so `tag` is the continental fraction of everything leaving it this
    # step, and every pathway below moves `volume * tag` along exactly the same transfers.
    #
    # With mobile cover, the continental share of what a node gives up comes from two layers:
    # the cover, well mixed, at its own continental fraction, then the substrate under it, by
    # the same continental-first rule. The leveling grind and glacial flattening take cover
    # first too.
    area = _gather_areas(world, plates_in_order)
    removed_m = erosion_amount + sea_side_erosion + ground_off + flatten_removed
    for removal in (ground_off, flatten_removed):
        taken = np.minimum(cover_m, removal)
        cover_entrained = cover_entrained + taken
        cover_m = cover_m - taken
    cover_fraction = np.divide(prior_cover_material, prior_cover, out=np.zeros(n), where=prior_cover > 0.0)
    substrate_material = prior_material - prior_cover_material
    substrate_removed_m = removed_m - cover_entrained
    continental_removed_m = cover_entrained * cover_fraction + np.minimum(substrate_removed_m, substrate_material)
    tag = np.divide(continental_removed_m, removed_m, out=np.zeros(n), where=removed_m > 0)

    # ocean_deposition_multiplier scales the *settled* marine sediment (beach + marine below),
    # not the pre-spread pool -- scaling the pool would desync the mass-conserving np.add.at
    # spread against the deep-water remainder. At 1.0 this is exact; away from 1.0 it's a
    # deliberate shelf-building / shelf-starving knob. Above 1.0 the extra is ordinary
    # (non-continental) material; below 1.0 the continental share it declines to settle goes
    # to the ledger's declared discarded_marine_sediment_m3 sink rather than vanishing.
    ocean_multiplier = world.ocean_deposition_multiplier
    ocean_tag_keep = min(ocean_multiplier, 1.0)
    discarded_tagged_m3 = 0.0

    # Deposition: wherever a big (water_accum_m > DEPOSITION_MIN_FLOW_M), slow
    # (river_speed < DEPOSITION_SPEED_THRESHOLD) river passes through, DEPOSITION_FRACTION
    # of the material passing through settles right there instead of continuing downstream
    # -- route_downstream still conserves the total exactly either way.
    river_speed = hydrology.compute_river_speed(slope, hydro.flow_accum)
    is_depositing = (river_speed < DEPOSITION_SPEED_THRESHOLD) & (water_accum_m > DEPOSITION_MIN_FLOW_M)
    # river_deposition_multiplier scales the settle-out fraction; clamped below 1.0 so a big
    # multiplier can't make a reach retain more than passes through it (route_downstream still
    # conserves the routed total exactly at any fraction in [0, 1)).
    deposition_fraction = float(np.clip(DEPOSITION_FRACTION * world.river_deposition_multiplier, 0.0, 0.95))
    retain_fraction = np.where(is_depositing, deposition_fraction, 0.0)
    water_routed_vol = water_routed_amount * area
    _, water_routed_deposit = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.flow_target, water_routed_vol, retain_fraction=retain_fraction
    )
    _, water_routed_tagged = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.flow_target, water_routed_vol * tag, retain_fraction=retain_fraction
    )

    # Marine/beach sediment: whatever water_routed_deposit above piled onto a single ocean
    # node (wherever that node's flow path happened to terminate) spreads across nearby
    # shallow coast instead -- see _spread_beach_sediment's own docstring.
    beach_deposit = _spread_beach_sediment(points, elevation, is_ocean_node, np.where(is_ocean_node, water_routed_deposit, 0.0))
    beach_tagged = _spread_beach_sediment(points, elevation, is_ocean_node, np.where(is_ocean_node, water_routed_tagged, 0.0))
    discarded_tagged_m3 += float(beach_tagged.sum()) * (1.0 - ocean_tag_keep)
    # Land-side deposits (floodplain retention, dead-end-basin sinks) spread across a lake/sea's
    # whole flooded extent instead of piling onto its single sink node -- see
    # LAKE_SEDIMENT_UNIFORM_FRACTION/_spread_lake_sediment's own comments. A dry closed basin (no
    # standing water yet) is untouched by this, same old single-point-pile behavior. The beach
    # spread only ever lands on ocean nodes and the lake spread only on land, so their sum is
    # the old per-node `where(is_ocean, beach, lake)` pick without risking dropping either.
    lake_deposit = _spread_lake_sediment(hydro.lake_depth, hydro.neighbor_idx, np.where(is_ocean_node, 0.0, water_routed_deposit))
    lake_tagged = _spread_lake_sediment(hydro.lake_depth, hydro.neighbor_idx, np.where(is_ocean_node, 0.0, water_routed_tagged))
    sediment_vol = beach_deposit * ocean_multiplier + lake_deposit
    sediment_tagged = beach_tagged * ocean_tag_keep + lake_tagged

    wind_vol = wind_redeposit_source * area
    wind_deposit = _route_wind_deposit(points, wind_u_at_nodes, wind_v_at_nodes, wind_vol, world=world)
    wind_tagged = _route_wind_deposit(points, wind_u_at_nodes, wind_v_at_nodes, wind_vol * tag, world=world)

    # Glacial transport: glacier_carried travels along the ice's own real flow path
    # (hydro.ice_flow_target, not water's flow_target -- see GLACIER_TILL_FRACTION's own
    # comment) and settles the moment it reaches a node that's genuinely outside the ice
    # (glacier_depth below hydrology.GLACIER_VISIBLE_DEPTH_M) -- a terminal moraine/outwash
    # deposit at the glacier's melting margin, pushed there by the glacier's own flow rather
    # than left buried under the ice interior. Reuses hydro.glacier_depth (this step's already-
    # flowed, *final* ice depth) to find the margin, not the one-step-lagged prior_glacier_depth
    # ice_factor above is deliberately built from -- this is about where the ice sits *after*
    # this step's own flow, not about damping the erosion-rate formula.
    at_glacier_margin = np.where(hydro.glacier_depth < hydrology.GLACIER_VISIBLE_DEPTH_M, 1.0, 0.0)
    glacier_carried_vol = glacier_carried * area
    _, glacier_transport_deposit = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.ice_flow_target, glacier_carried_vol, retain_fraction=at_glacier_margin
    )
    _, glacier_transport_tagged = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.ice_flow_target, glacier_carried_vol * tag, retain_fraction=at_glacier_margin
    )

    # Mass wasting (issue #275 phase 4): landslide debris runs out downslope by gravity and
    # settles on the foreland, never past a column's Hc cap -- see `_route_mass_wasting`.
    # Debris that falls onto a glacier rides the ice instead: along ice_flow_target, the same
    # path glacial debris takes, out to the ice margin -- and out through a basin's lowest
    # outlet once the ice overtops its rim. Wherever the ice drops it, it runs out by gravity
    # again from there, under the same Hc room, so the ice can't pile it past the cap. Debris
    # landing in a lake spreads across it the same capacity-aware way. What reaches the sea
    # spreads onto the shelf and into the basin like the other marine sediment, under the same
    # ocean_deposition_multiplier.
    mass_wasting_vol = applied_seismic * area
    landslide_room = np.where(has_column, lithosphere.MAX_CRUSTAL_THICKNESS_M - prior_hc + removed_m, np.inf) * area
    headroom_m = np.where(has_column, lithosphere.MAX_CRUSTAL_THICKNESS_M - prior_hc, np.inf)
    landslide_land, landslide_arrival, landslide_ice, landslide_land_tagged, landslide_arrival_tagged, landslide_ice_tagged = _route_mass_wasting(
        elevation,
        is_ocean_node,
        hydro.neighbor_idx,
        slope,
        mass_wasting_vol,
        mass_wasting_vol * tag,
        capacity_vol=landslide_room,
        headroom_m=headroom_m,
        spill_target=hydro.spill_target,
        on_ice=hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M,
    )
    _, ice_dropped = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.ice_flow_target, landslide_ice, retain_fraction=at_glacier_margin
    )
    _, ice_dropped_tagged = hydrology.route_downstream(
        elevation, is_ocean_node, hydro.ice_flow_target, landslide_ice_tagged, retain_fraction=at_glacier_margin
    )
    off_ice_land, off_ice_arrival, _, off_ice_land_tagged, off_ice_arrival_tagged, _ = _route_mass_wasting(
        elevation,
        is_ocean_node,
        hydro.neighbor_idx,
        slope,
        np.zeros(n),
        np.zeros(n),
        capacity_vol=landslide_room - landslide_land,
        headroom_m=headroom_m,
        spill_target=hydro.spill_target,
        inflow_vol=ice_dropped,
        inflow_tagged=ice_dropped_tagged,
    )
    landslide_land = landslide_land + off_ice_land
    landslide_land_tagged = landslide_land_tagged + off_ice_land_tagged
    landslide_arrival = landslide_arrival + off_ice_arrival
    landslide_arrival_tagged = landslide_arrival_tagged + off_ice_arrival_tagged
    landslide_land, landslide_land_tagged = _spread_lake_sediment_capped(
        hydro.lake_depth, hydro.neighbor_idx, landslide_land, landslide_land_tagged, capacity=landslide_room
    )
    landslide_marine_unscaled = _spread_marine_sediment(points, elevation, is_ocean_node, landslide_arrival)
    landslide_marine_tagged = _spread_marine_sediment(points, elevation, is_ocean_node, landslide_arrival_tagged)
    discarded_tagged_m3 += float(landslide_marine_tagged.sum()) * (1.0 - ocean_tag_keep)
    landslide_marine = landslide_marine_unscaled * ocean_multiplier
    landslide_marine_tagged = landslide_marine_tagged * ocean_tag_keep

    # Distributary redirect: route_downstream funnels a slow river's retained load down one
    # discretised channel and drops it on whichever single near-sea-level node that channel
    # runs through. A real delta splits into distributaries and spreads that load across the
    # whole flat fan -- pull DELTA_REDIRECT_FRACTION of it back off those band land nodes and
    # hand it to the leveling fill spread, which scatters it across the band's hollows.
    in_delta_band = (~is_ocean_node) & is_depositing & (np.abs(elevation - world.sea_level_m) <= COASTAL_LEVELING_BAND_M)
    delta_redirect = np.where(in_delta_band, DELTA_REDIRECT_FRACTION * sediment_vol, 0.0)
    delta_redirect_tagged = np.where(in_delta_band, DELTA_REDIRECT_FRACTION * sediment_tagged, 0.0)
    sediment_vol = sediment_vol - delta_redirect
    sediment_tagged = sediment_tagged - delta_redirect_tagged

    leveling_source = (ground_off + COASTAL_INFILL_MARINE_FRACTION * sea_side_erosion) * area + delta_redirect
    leveling_source_tagged = (ground_off + COASTAL_INFILL_MARINE_FRACTION * sea_side_erosion) * area * tag + delta_redirect_tagged
    marine_source = sea_side_erosion * (1.0 - COASTAL_INFILL_MARINE_FRACTION) * area
    marine_unscaled = _spread_marine_sediment(points, elevation, is_ocean_node, marine_source)
    marine_tagged = _spread_marine_sediment(points, elevation, is_ocean_node, marine_source * tag)
    discarded_tagged_m3 += float(marine_tagged.sum()) * (1.0 - ocean_tag_keep)
    marine_deposit = marine_unscaled * ocean_multiplier + landslide_marine
    marine_tagged = marine_tagged * ocean_tag_keep + landslide_marine_tagged
    leveling_fill, leveling_tagged = _spread_coastal_leveling(
        points, elevation, coastal_openness, dist_to_land, world.sea_level_m, leveling_datum, leveling_source, dt_myr, local_relief_m,
        area_m2=area, tagged_amount=leveling_source_tagged,
    )

    # Flattening's per-edge sends land on the neighbour they were sent to.
    flatten_vol = flatten_send * area[:, None]
    flatten_received = np.zeros(n)
    flatten_received_tagged = np.zeros(n)
    np.add.at(flatten_received, hydro.neighbor_idx.ravel(), flatten_vol.ravel())
    np.add.at(flatten_received_tagged, hydro.neighbor_idx.ravel(), (flatten_vol * tag[:, None]).ravel())

    # Everything above is a volume per receiving node; back to thickness there.
    till_vol = glacier_till * area
    plain_deposit_vol = sediment_vol + till_vol + glacier_transport_deposit + wind_deposit + landslide_land
    deposited_vol = plain_deposit_vol + marine_deposit + leveling_fill + flatten_received
    deposited_tagged = (
        sediment_tagged
        + till_vol * tag
        + glacier_transport_tagged
        + wind_tagged
        + landslide_land_tagged
        + marine_tagged
        + leveling_tagged
        + flatten_received_tagged
    )

    # Depositional Hc-cap overflow (issue #288, see OVERFLOW_*). Lake / endorheic-basin
    # siltation counts too: the sediment that settled out of standing water this step
    # (hydrology.step_lakes -> silt_deposited) is folded straight into elevation, so a
    # still-water basin genuinely fills in and stays filled -- a small non-continental source
    # (it isn't drawn from any eroded pool), declared in the budget below. Each node's incoming
    # load is split by pathway (columns in OVERFLOW_PATHWAYS order, then the continental share);
    # whatever exceeds the room its column has left under the cap after this step's own
    # removals is carried on to receivers with room instead of being clipped off.
    is_lake_node = (~is_ocean_node) & (hydro.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
    land_sediment_vol = np.where(is_ocean_node, 0.0, sediment_vol)
    incoming = np.stack(
        [
            np.where(is_lake_node, 0.0, land_sediment_vol),
            np.where(is_lake_node, land_sediment_vol, 0.0),
            hydro.silt_deposited * area,
            till_vol + glacier_transport_deposit + flatten_received,
            wind_deposit,
            landslide_land + landslide_marine,
            leveling_fill,
            np.where(is_ocean_node, sediment_vol, 0.0) + marine_unscaled * ocean_multiplier,
            deposited_tagged,
        ],
        axis=1,
    )
    incoming_vol = incoming[:, :-1].sum(axis=1)
    cap_room_vol = np.where(has_column, np.clip((lithosphere.MAX_CRUSTAL_THICKNESS_M - prior_hc + removed_m) * area, 0.0, None), np.inf)
    excess_vol = np.clip(incoming_vol - cap_room_vol, 0.0, None)
    excess = incoming * np.divide(excess_vol, incoming_vol, out=np.zeros(n), where=excess_vol > 0.0)[:, None]
    # The ice can scour a node only within the drop to its lowest neighbour and the Hc headroom
    # this step's removals left, scaled by craton resistance like every other removal.
    scour_limit_m = (
        np.minimum(
            OVERFLOW_ICE_SCOUR_MAX_M_PER_MYR * dt_myr,
            np.clip(np.minimum(drop_to_lowest_neighbor_m, removable_m) - removed_m, 0.0, None),
        )
        * craton_keep
    )
    carry = _carry_overflow(
        points,
        elevation,
        is_ocean_node,
        hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M,
        hydro.lake_depth,
        hydro.neighbor_idx,
        hydro.spill_target,
        hydro.ice_flow_target,
        area,
        np.where(has_column, np.clip(cap_room_vol - incoming_vol, 0.0, None), np.inf),
        excess,
        scour_limit_m,
        np.clip(substrate_material - substrate_removed_m, 0.0, None),
        world=world,
        scour_cover_m=cover_m,
        scour_cover_material_m=cover_m * cover_fraction,
    )
    settled = incoming - excess + carry.placed
    settled_tagged = settled[:, -1]
    removed_m = removed_m + carry.scour_m
    continental_removed_m = continental_removed_m + carry.scour_tagged_m
    cover_entrained = cover_entrained + carry.scour_cover_m
    cover_m = cover_m - carry.scour_cover_m
    if carry.terminal[-1] > 0.0:
        continental_ledger.record(world, "overloaded_root_delaminated_m3", float(carry.terminal[-1]))

    plain_deposition = plain_deposit_vol / area
    marine_deposit_m = marine_deposit / area
    leveling_fill_m = leveling_fill / area
    flatten_delta = flatten_received / area - flatten_removed
    # Net of the overflow: what full receivers turned away left them, what it was carried to
    # gained it (including the rock the ice scoured on the way, a removal of its own).
    overflow_net_m = (carry.placed[:, :-1] - excess[:, :-1])[:, OVERFLOW_SEDIMENT].sum(axis=1) / area
    total_deposited = plain_deposition + marine_deposit_m + leveling_fill_m + overflow_net_m

    geomorphic_delta = -removed_m + settled[:, :-1].sum(axis=1) / area

    # Erosional isostatic compensation. Every term above moves rock between columns but, on
    # its own, never told isostasy: `elevation` used to absorb the whole change, drifting
    # ever further below isostatic_elevation(Hc, Hm) as coastal + submarine erosion shipped
    # continental crust off to the abyss with no rebound -- which planed every continent flat
    # over a few hundred Myr once orogeny slowed (GitHub issue #120, "Land fraction slowly
    # declines"). Now the full rock-column change books against Hc and `elevation` moves by
    # exactly the resulting Airy response -- so the unloaded crustal root rebounds (only
    # ~1/6 of subaerial erosion, ~1/4 of submarine, survives as a surface drop) and a
    # sediment pile subsides under its own load, the same delta idiom deform() already uses
    # for tectonic Hc/Hm changes. `elevation` stays a faithful readout of the column, so
    # deform()'s mechanism stays exact. v1 PlateWithLines carries no Hc (all-zero) -- those
    # nodes keep the bare 1:1 response.
    rho_c_per_node = np.concatenate(
        [lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in plates_in_order]
    )
    # Clipped to [MIN, MAX]_CRUSTAL_THICKNESS_M (issue #161), the same ceiling `rheology`'s
    # tectonic thickening paths enforce. The removal caps above keep the lower clip from
    # binding on a column that started above the floor, and the overflow carry keeps deposits
    # within the room under the upper one, so either clip only trims round-off now (reported
    # as `hc_clip_residual_m3`, its continental share as numerical_unplaced_m3).
    raw_crustal_thickness = prior_hc + geomorphic_delta
    new_crustal_thickness = np.where(
        has_column,
        np.clip(raw_crustal_thickness, lithosphere.MIN_CRUSTAL_THICKNESS_M, lithosphere.MAX_CRUSTAL_THICKNESS_M),
        prior_hc,
    )
    hc_clip_residual_m = np.where(has_column, np.clip(raw_crustal_thickness - lithosphere.MAX_CRUSTAL_THICKNESS_M, 0.0, None), 0.0)
    isostatic_delta = lithosphere.isostatic_elevation(
        new_crustal_thickness, prior_hm, rho_c_per_node
    ) - lithosphere.isostatic_elevation(prior_hc, prior_hm, rho_c_per_node)
    applied_delta = np.where(has_column, isostatic_delta, geomorphic_delta)
    eroded_elevation = np.clip(elevation + applied_delta, MIN_ELEVATION_M, MAX_ELEVATION_M)

    # Ice loading (issue #275, phase 3) -- see `_apply_ice_load`.
    prior_deflection = collect_all_ice_load_deflection(plates_in_order)
    prior_load = lithosphere.grounded_ice_load_kg_m2(prior_glacier_depth, elevation - prior_deflection, world.sea_level_m)
    new_elevation, new_deflection, ice_load = _apply_ice_load(
        eroded_elevation, prior_deflection, hydro.glacier_depth, new_crustal_thickness, prior_hm, rho_c_per_node, world.sea_level_m
    )
    # Glaciers carve channels too: abrasion and the overflow-laden ice's scour both deepen the
    # trough they cut, which rivers inherit (channel_boost, hydrology's channel-aware routing)
    # once the ice retreats. Like river erosion, only rock actually removed here counts.
    carved_m = applied_river + applied_glacier + carry.scour_m
    new_channel_depth = np.where(is_ocean_node, 0.0, np.clip(prior_channel_depth + carved_m, 0.0, MAX_CHANNEL_DEPTH_M))
    # Width grows with discharge alone (no slope/channel_boost term -- see module constants'
    # own comment for why), same persistent/monotonic/capped shape as depth.
    width_growth = WIDTH_GROWTH_COEFFICIENT * np.power(np.clip(water_accum_m, 0.0, None), WIDTH_FLOW_EXPONENT) * dt_myr
    new_channel_width = np.where(is_ocean_node, 0.0, np.clip(prior_channel_width + width_growth, 0.0, MAX_CHANNEL_WIDTH_M))

    # The tracer moves with its rock and can never exceed the column holding it; whatever the
    # column can't hold (clip round-off, or a no-column v1 node) is declared, not dropped.
    raw_material = prior_material - continental_removed_m + settled_tagged / area
    new_material = np.clip(raw_material, 0.0, new_crustal_thickness)
    unplaced_tagged_m3 = float(np.sum((raw_material - new_material) * area))
    if unplaced_tagged_m3 > 0.0:
        continental_ledger.record(world, "numerical_unplaced_m3", unplaced_tagged_m3)
    if discarded_tagged_m3 > 0.0:
        continental_ledger.record(world, "discarded_marine_sediment_m3", discarded_tagged_m3)

    # Everything that settled this step becomes mobile cover, carrying its continental share.
    # The cover stays inside its column (and its continental share inside both the cover and
    # the tracer) -- removals only reach substrate once the cover is gone, so the clips only
    # trim round-off, reported as `mobile_cover_clip_m3`. Then burial consolidates whatever lies
    # deeper than MOBILE_COVER_CONSOLIDATION_DEPTH_M back into substrate, its continental share
    # with it -- Hc and the tracer are unchanged, only the split moves.
    cover_deposited_m = settled[:, :-1].sum(axis=1) / area
    raw_cover = cover_m + cover_deposited_m
    new_cover = np.clip(raw_cover, 0.0, np.where(has_column, new_crustal_thickness, np.inf))
    new_cover_material = np.clip(
        cover_m * cover_fraction + settled_tagged / area, 0.0, np.minimum(new_cover, new_material)
    )
    consolidated_m = np.clip(new_cover - MOBILE_COVER_CONSOLIDATION_DEPTH_M, 0.0, None) * -np.expm1(
        -dt_myr / MOBILE_COVER_CONSOLIDATION_TIMESCALE_MYR
    )
    consolidated_material_m = consolidated_m * np.divide(new_cover_material, new_cover, out=np.zeros(n), where=new_cover > 0.0)
    settled_cover = new_cover
    new_cover = new_cover - consolidated_m
    new_cover_material = np.clip(new_cover_material - consolidated_material_m, 0.0, new_cover)
    removed_m3 = float(np.sum(removed_m * area))
    # Everything eroded this step that settled on the surface, after the overflow carry --
    # removed_m3 == deposited_m3 + hc_cap_overflow_terminal_m3 less its lake-silt share (with
    # the ocean knob at 1.0).
    deposited_m3 = float(settled[:, :-1][:, OVERFLOW_SEDIMENT].sum())
    lake_silt_m3 = float(np.sum(hydro.silt_deposited * area))
    overflow_hits = excess[:, :-1] > 0.0
    overflow_budget = {
        # Issue #288: depositional Hc-cap overflow -- how many receivers turned load away, how
        # much, how much was carried on to a receiver with room, and how much found none (the
        # overloaded_root_delaminated_m3 reservoir). Per pathway below; placed may exceed
        # attempted for `glacial` by the rock the overflow-laden ice scoured.
        "hc_cap_overflow_nodes": float(np.count_nonzero(overflow_hits.any(axis=1))),
        "hc_cap_overflow_m3": float(excess_vol.sum()),
        "hc_cap_overflow_placed_m3": float(carry.placed[:, :-1].sum()),
        "hc_cap_overflow_terminal_m3": float(carry.terminal[:-1].sum()),
        "hc_cap_overflow_ice_scour_m3": float(np.sum(carry.scour_m * area)),
        "continental_overflow_m3": float(excess[:, -1].sum()),
        "continental_overflow_terminal_m3": float(carry.terminal[-1]),
        "hc_clip_residual_m3": float(np.sum(hc_clip_residual_m * area)),
        **{f"hc_cap_overflow_stuck_{reason}_m3": volume for reason, volume in carry.stuck_m3.items()},
        "hc_cap_overflow_basin_fill_m3": carry.basin_fill_m3,
        "hc_cap_overflow_nearest_search_m3": carry.search_m3,
    }
    for column, pathway in enumerate(OVERFLOW_PATHWAYS):
        overflow_budget[f"hc_cap_overflow_{pathway}_nodes"] = float(np.count_nonzero(overflow_hits[:, column]))
        overflow_budget[f"hc_cap_overflow_{pathway}_m3"] = float(excess[:, column].sum())
        overflow_budget[f"hc_cap_overflow_{pathway}_placed_m3"] = float(carry.placed[:, column].sum())
        overflow_budget[f"hc_cap_overflow_{pathway}_terminal_m3"] = float(carry.terminal[column])
    # Lake silt settles into the mobile cover and is mixed through it like everything else,
    # so the silt record follows it: removals take the prior silt in proportion to the cover
    # they take, this step's settled silt joins, and consolidation takes its share of both.
    # Silt laid down before the cover existed was already substrate (hydrology caps the prior
    # record at the cover when it reads it, which also catches tectonic stripping).
    prior_silt = np.clip(hydro.silt_depth - hydro.silt_deposited, 0.0, None)
    kept_silt = prior_silt * np.divide(cover_m, prior_cover, out=np.zeros(n), where=prior_cover > 0.0)
    settled_silt = settled[:, OVERFLOW_LAKE_SILT] / area
    unconsolidated = np.divide(new_cover, settled_cover, out=np.zeros(n), where=settled_cover > 0.0)
    new_silt_depth = np.minimum((kept_silt + settled_silt) * unconsolidated, new_cover)
    # Keep this step's cached hydrology in step with what the plates now hold.
    hydro.silt_depth = new_silt_depth
    cover_entrained_m3 = float(np.sum(cover_entrained * area))
    mobile_cover_budget = {
        # Issue #297 phase 2: what this step's removals took from the cover and from substrate,
        # what settled into the cover, what burial consolidated, what is left, and the
        # continental share of each. prior + deposited - entrained - consolidated - clip ==
        # remaining.
        "bedrock_detached_m3": removed_m3 - cover_entrained_m3,
        "mobile_cover_prior_m3": float(np.sum(prior_cover * area)),
        "mobile_cover_entrained_m3": cover_entrained_m3,
        "mobile_cover_deposited_m3": float(np.sum(cover_deposited_m * area)),
        "mobile_cover_consolidated_m3": float(np.sum(consolidated_m * area)),
        "mobile_cover_clip_m3": float(np.sum((raw_cover - new_cover - consolidated_m) * area)),
        "mobile_cover_m3": float(np.sum(new_cover * area)),
        "mobile_cover_continental_entrained_m3": float(np.sum(cover_entrained * cover_fraction * area)),
        "mobile_cover_continental_consolidated_m3": float(np.sum(consolidated_material_m * area)),
        "mobile_cover_continental_m3": float(np.sum(new_cover_material * area)),
        "stale_mobile_cover_excess_m3": float(np.sum((stored_cover - prior_cover) * area)),
    }
    mobile_cover.record(world, "deposited_m3", mobile_cover_budget["mobile_cover_deposited_m3"])
    mobile_cover.record(world, "entrained_m3", cover_entrained_m3)
    mobile_cover.record(world, "consolidated_m3", mobile_cover_budget["mobile_cover_consolidated_m3"])
    mobile_cover.record(
        world, "clipped_m3", mobile_cover_budget["mobile_cover_clip_m3"] + mobile_cover_budget["stale_mobile_cover_excess_m3"]
    )
    budget = {
        "removed_m3": removed_m3,
        "deposited_m3": deposited_m3,
        # Net material the ocean_deposition_multiplier knob adds (> 1) or withholds (< 1).
        "ocean_deposition_knob_m3": float(beach_deposit.sum() + marine_unscaled.sum() + landslide_marine_unscaled.sum())
        * (ocean_multiplier - 1.0),
        # Landslide debris (phase 4): mobilized, settled on land (foreland), reaching the sea,
        # handed to the ice.
        "mass_wasting_removed_m3": float(mass_wasting_vol.sum()),
        "mass_wasting_foreland_m3": float(landslide_land.sum()),
        "mass_wasting_marine_m3": float(landslide_arrival.sum()),
        "mass_wasting_on_ice_m3": float(landslide_ice.sum()),
        "lake_silt_m3": lake_silt_m3,
        **overflow_budget,
        "continental_removed_m3": float(np.sum(continental_removed_m * area)),
        "continental_deposited_m3": float(settled_tagged.sum()),
        "continental_discarded_m3": discarded_tagged_m3,
        "continental_unplaced_m3": unplaced_tagged_m3,
        "stale_tracer_excess_m3": float(np.sum((stored_material - prior_material) * area)),
        **mobile_cover_budget,
    }

    # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_* and
    # render_image's "elevReason" view). Group this step's geomorphic contributions and stamp
    # each node with whichever moved it most -- but only where |net change| clears
    # ELEV_CHANGE_MIN_DELTA_M, so a low-relief, low-rainfall node the coastal/erosion passes
    # only brush by sub-metre keeps whatever last genuinely shaped it (tectonics, or NONE --
    # untouched since generation). That gate is the point of the view: it separates land
    # that's actually being planed flat now from land that was simply never built up.
    reason_contrib = np.stack(
        [
            erosion_amount + carry.scour_m,
            np.clip(plain_deposition + overflow_net_m, 0.0, None),
            ground_off + leveling_fill_m,
            sea_side_erosion + marine_deposit_m,
            np.abs(flatten_delta),
            hydro.silt_deposited,
        ],
        axis=-1,
    )
    reason_codes = np.array(
        [
            ELEV_CHANGE_EROSION,
            ELEV_CHANGE_DEPOSITION,
            ELEV_CHANGE_COASTAL_LEVELING,
            ELEV_CHANGE_MARINE,
            ELEV_CHANGE_GLACIAL_FLATTEN,
            ELEV_CHANGE_LAKE_SILT,
        ],
        dtype=float,
    )
    prior_reason = collect_all_elev_change_reason(plates_in_order)
    dominant = reason_codes[np.argmax(reason_contrib, axis=-1)]
    # Provenance tracks which *process* is reshaping the column, so it keys off the raw
    # geomorphic move (rock added/removed), not `new_elevation - elevation` -- isostatic
    # compensation shrinks the surface expression ~5x but doesn't change what's doing the
    # shaping, and the two thresholds below were tuned against the raw pre-compensation move.
    net_abs = np.abs(geomorphic_delta)
    moved = net_abs >= ELEV_CHANGE_MIN_DELTA_M
    # A structural code (deform()/volcanism, re-stamped every step the belt is still active) is
    # sticky against ordinary background wash -- only a large net geomorphic step overrides it.
    prior_structural = (
        ((prior_reason >= ELEV_CHANGE_COLLISION) & (prior_reason <= ELEV_CHANGE_VOLCANO))
        | ((prior_reason >= ELEV_CHANGE_FAULT_NORMAL) & (prior_reason <= ELEV_CHANGE_FAULT_STRIKE_SLIP))
        | (prior_reason == ELEV_CHANGE_LATERAL_MAGMA)
    )
    override_structural = net_abs >= ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR * dt_myr
    overwrite = moved & (~prior_structural | override_structural)
    new_elev_change_reason = np.where(overwrite, dominant, prior_reason)

    # theta (and therefore every other parallel array's shape) is never touched here --
    # writing each changed field straight back via set_fields_on_plate (a vectorized
    # per-plate slice write, no per-node point object) leaves every other persistent field
    # (is_volcano/volcano_active_years_remaining/soil_*/resource deposits, ...) untouched
    # automatically, the same "don't silently wipe fields this site doesn't know about"
    # guarantee line.replace used to provide.
    offset = 0
    for plate in plates_in_order:
        n = plate.node_count()
        plate.set_fields_on_plate(
            elevation=new_elevation[offset : offset + n],
            crustal_thickness_m=new_crustal_thickness[offset : offset + n],
            channel_depth=new_channel_depth[offset : offset + n],
            channel_width=new_channel_width[offset : offset + n],
            lake_depth=hydro.lake_depth[offset : offset + n],
            glacier_depth=hydro.glacier_depth[offset : offset + n],
            silt_depth=new_silt_depth[offset : offset + n],
            elev_change_reason=new_elev_change_reason[offset : offset + n],
            continental_material_m=new_material[offset : offset + n],
            ice_load_deflection_m=new_deflection[offset : offset + n],
            mobile_cover_m=new_cover[offset : offset + n],
            mobile_cover_continental_m=new_cover_material[offset : offset + n],
        )
        offset += n

    return ErosionResult(
        points=points,
        elevation=elevation,
        slope=slope,
        rain=rain,
        river=river,
        weathering=weathering,
        sediment_deposited=total_deposited,
        is_river_depositing=is_depositing,
        net_elevation_change_m=eroded_elevation - elevation,
        temperature_c=temperature,
        precipitation_mm=precipitation_mm,
        budget=budget,
        ice_load_pa=lithosphere.ice_load_pa(ice_load),
        ice_load_change_pa=lithosphere.ice_load_pa(ice_load - prior_load),
        ice_deflection_change_m=new_deflection - prior_deflection,
    )
