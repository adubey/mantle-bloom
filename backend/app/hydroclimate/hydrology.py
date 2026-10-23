"""Rivers, lakes, and glaciers: flow routing, basin/lake detection, ice accumulation/melt/
flow, and downstream flow accumulation over the world's current node cloud.

A fixed grid's regular 8-neighbor adjacency would give flow routing, depression filling, and
downstream accumulation a natural substrate, but mantle-bloom's world is instead an irregular
per-plate node cloud. Steepest-descent flow direction and elevation-ordered downstream
accumulation turn out not to actually need a *grid*, only a *graph*: this module builds one via
a whole-world k-nearest-neighbor query (the same technique reassign.py/erosion.py already use
for their own whole-world passes), then runs both algorithms directly on it. Depression
filling/basin-spill -- deciding where a lake's overflow escapes to once full -- is instead
derived from `lakes.py`'s own catchment+Kruskal depression hierarchy
(`lakes.compute_spill_routing`), not a separately-computed priority-flood pass over this
module's own graph: an earlier version of this module *did* run its own independent
multi-source-Dijkstra priority-flood here (`_compute_basin_spill`, since removed), but it
could drift out of sync with `lakes.py`'s own, physically-correct hierarchy -- see
`lakes.compute_spill_routing`'s own docstring for the real, save-file-confirmed drift this
caused, and `compute_hydrology`'s own docstring for the resulting step order.

Glacier ice is solved by `ice_sheet.py` (issue #353): volume conserved per ice body, shaped
by a plastic ice-sheet profile. This module calls it and owns the result's place in the step
(lake ice in, meltwater out to the rivers, frozen river ice in a step later). Material riding
the ice has its own target (`ice_flow_target`, on HydrologyFields) rather than water's: down
the ice surface on the ice, and down the ground off it, never stopped by freezing, whereas
liquid water's flow_target is deliberately forced shut at a frozen node (rivers freeze over)
-- see `_compute_flow_direction`'s `apply_freeze` parameter.

**Persistence.** A grid-based approach where plates move relative to fixed cells would need a
persistent field like channel_depth to be deliberately advected (semi-Lagrangian, every step)
to keep following the crust -- mantle-bloom's plate-local surface nodes sidestep that
entirely, since they already rotate exactly with their own plate. lake_depth and glacier_depth,
stored as ordinary surface fields right alongside elevation itself (see surface_fields.py),
get that same "just works" persistence for free: no advection scheme needed, since rotating a
plate never touches those arrays at all. flow_target/flow_accum/river_speed are deliberately
*not* persisted -- they're recomputed fresh every step, from that step's real climate, purely
this-step derived quantities, cached on World (see World.hydrology_cache) only so a later
same-turn caller (rendering, stats) doesn't recompute them again.

**Flow routing is computed once, not twice.** Recomputing flow routing separately for
rendering fields and for the water_accum erosion needs would be an accepted redundancy under
a JIT compiler, where the repeated pass is cheap. Flow routing here has no JIT and is
comparatively expensive, so this module computes it once and erosion.py reuses the same
result, rather than duplicating a cost that's real here. This module owns lake and glacier
*state transitions* (this step's final lake_depth/glacier_depth, both returned via
HydrologyFields) -- erosion.py only reads glacier_depth (for its own glacier-erosion/
flattening terms) and channel_depth (which erosion.py itself owns, since it's the module that
actually carves it).

**Lakes are an explicit object tree, not a flat per-node approximation** -- see `lakes.py`
(`lakes.Lake`, `lakes.build_lake_hierarchy`, `lakes.step_lakes`), and that same tree is now
also what drives spill *routing* here, not just water-balance bookkeeping (see above). An
earlier version of this
module tracked lakes purely as a flood-fill over `lake_depth`, seeded from every node that held
water last step and rate-limited per step (a hop cap on growth into new territory, a cap on
per-node depth change) so a large lake's size changed gradually rather than snapping instantly
to whatever a step's basin geometry technically allowed. That bounded *rate* but not *extent* --
a sufficiently large, gently-sloped catchment could keep creeping outward step after step
without ever reaching a real equilibrium (the user-reported bug this rewrite fixes). `lakes.py`
instead computes each basin's actual floor/rim geometry directly (a depression-hierarchy /
priority-flood-plus technique), so a lake's maximum extent is known immediately rather than
discovered incrementally, and still projects down into this same flat `lake_depth` array
afterward so every other consumer here (`is_river` below, `group_rivers`, and
render_image.py/coastline.py/erosion.py outside this module) keeps reading it unchanged. A
river's flow_target can point at a node that's part of a lake (real inflow), but a flooded node
itself is never classified `is_river` -- a river ends at the lake's shore rather than its
classification continuing straight across the water to whatever's on the far side (see
`compute_hydrology`).
Flow *routing* itself also isn't pure steepest descent: `_compute_flow_direction` prefers a
downhill neighbor that already has a real, established channel over the literal lowest one,
so a river stays in its own valley rather than recalculating the mathematically steepest path
every step -- real rivers meander within, not across, their banks.

**A lake that finds a way out erodes its outlet hard.** Once a lake is actively spilling
(`should_spill`, from `_compute_flow_direction`), `compute_hydrology` injects extra "water"
at its sink sized by the lake's own connected surface area (node count), which then rides
`route_downstream` through the whole outflow channel like real precipitation would -- so
erosion.py's existing river-erosion and channel-growth formulas, unchanged, carve that
spillway into a real, deep channel proportional to how big the lake behind it actually is,
not just the trickle of ordinary flow passing through that one point (see
`LAKE_BREACH_EROSION_COEFFICIENT`).
"""

from __future__ import annotations

import heapq
import logging

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from numba import njit
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .. import cratons, geometry, healpix_grid
from . import breaching, lakes
from ..elevation_lines import DEFAULT_NODE_DENSITY, PLANET_RADIUS_KM, line_spacing_rad
from ..plates import (
    Plate,
    cached_node_position_tree,
    collect_all_accounting_areas_m2,
    collect_all_breach_notch_depth,
    collect_all_channel_depth,
    collect_all_channel_width,
    collect_all_elevation,
    collect_all_glacier_depth,
    collect_all_ice_load_deflection,
    collect_all_lake_depth,
    collect_all_silt_depth,
    gather_node_positions,
    query_workers,
)

if TYPE_CHECKING:
    from ..world import World
    from .ice_sheet import IceStep

_logger = logging.getLogger(__name__)

# Same 8-neighbor D8 structure as a regular grid's adjacency -- the graph every algorithm
# below runs on.
FLOW_NEIGHBOR_COUNT = 8

# Top ~7% of land physical discharge counts as "a river" -- used for both rendering and the
# major_river_fraction stat. Raised from an earlier 88.0 (top ~12%): at that lower bar a
# world's full drainage tree, down to its thinnest low-volume headwater tributaries, all
# registered as "a river" -- correct on a topographic map, but visually noisy on this one,
# and it let a tributary earn a place in the River Inspector's ranked list on hardly any flow
# at all. The major stems (which clear this same threshold by a wide margin at every point
# along their course) are unaffected; only the smallest branch tributaries -- exactly what
# should need to prove themselves with real accumulated flow before joining the display -- drop
# out.
RIVER_FLOW_PERCENTILE = 93.0

# How established a downhill neighbor's own channel needs to be before flow direction
# prefers it over the literal steepest candidate -- see _compute_flow_direction. Deliberately
# far below erosion.CHANNEL_BOOST_REFERENCE_M (200m, where channel_boost saturates): this
# only needs to detect "a real channel exists here at all," not "a mature one," so routing
# starts favoring it as soon as there's genuine, non-noise incision.
CHANNEL_PREFERENCE_THRESHOLD_M = 5.0

# Threshold above which a node's lake_depth counts as "a visible lake" -- used both for
# rendering (render_image.py, coastline.py) and for excluding flooded nodes from is_river/
# group_rivers below. The lake water balance itself (`lakes.LakeClimate`, `lakes.balanced_level`)
# lives in lakes.py, alongside the Lake tree that actually computes it -- see this module's own
# docstring and lakes.py's own module docstring for why.
LAKE_MIN_VISIBLE_DEPTH_M = 1.0

# `connected_ocean_mask`'s own "second ocean" threshold: any below-sea-level connected
# component whose node count is at least this fraction of the single largest one is *also*
# classified as ocean, not demoted to a giant lake. Without this, a world with two genuinely
# separate, ocean-scale bodies of water sees the smaller permanently misclassified as an
# endorheic basin -- and worse, right at the strait between them, a few meters of ordinary
# per-step erosion/tectonic drift is enough to open or close that connection, which flips the
# smaller body between "merged into the one ocean component" and "its own giant lake" every
# few steps (confirmed directly against a real save: an 11,151-node, -5856m-floor basin ~186km
# from the 77,748-node main ocean, oscillating in and out of the lake list). Picked from that
# same save's own size gap: the second-largest component there was 14.3% of the largest, while
# every genuine lake was under 0.4% -- 5% sits well inside that gap, generously catching a
# real second ocean while still demoting an actual Caspian-Sea-scale closed basin to a lake.
SECOND_OCEAN_SIZE_FRACTION = 0.05

# ...and an absolute floor alongside that fraction, since a *fraction* of the largest component
# alone can't tell a real second ocean apart from an ordinary small enclosed pit sitting in a
# proportionally tiny world -- this codebase's own unit tests exercise `connected_ocean_mask`
# on toy fixtures as small as ~15-40 total nodes, where an enclosed 3-node pit is already 20%
# of a 15-node "ocean" (see test_connected_ocean_mask_excludes_an_enclosed_interior_pit) despite
# being exactly the kind of small interior depression this function must still demote. A real
# generated world runs tens of thousands of nodes at minimum (node_density in world.py), so a
# floor comfortably above any toy fixture but far below a genuine ocean-scale component costs
# nothing on real worlds while keeping the small-fixture case exact.
SECOND_OCEAN_MIN_NODE_COUNT = 500

# A spilling lake's outlet channel erodes based on the *lake's* own surface area, not just
# the ordinary precip-routed flow passing through that one node. Real overflowing lakes carry
# far more erosive force at their outlet than an equivalent plain river of the same
# instantaneous discharge -- the whole standing body sits behind a narrow breach, not just
# whatever rain happened to fall immediately upstream this step -- and this is what actually
# cuts a lake's spillway down into a real gorge over time rather than a shallow trickle-carved
# groove. Expressed as an extra term added directly into `water_source` (in the same units as
# precipitation) at whichever sink node is spilling this step, sized by node-count of that
# lake's connected component (see `_lake_component_sizes` -- no separate "physical area"
# concept exists anywhere else in this codebase, so node count is the same implicit area
# proxy `flow_accum` itself already relies on). It then rides the ordinary `route_downstream`
# machinery downstream through the whole outflow channel to the sea exactly like real
# precipitation would, so erosion.py's existing river-erosion/channel-width formulas pick it
# up completely unchanged -- both channel_depth and channel_width grow from it automatically,
# "digging a deep channel out" without needing a second, parallel erosion pathway.
#
# This same term is *also* added directly at the lake's own rim node (`rim_breach_source`
# below, keyed by `lakes.compute_spill_routing`'s own `rim_node_idx`) -- not just at the sink
# it's routed from. The sink's own flow_target jumps straight to `spill_target` (the far side
# of the boundary edge, outside the lake), never actually passing through the lake's own rim
# member (`outlet_node_idx`) at all -- so a dam whose *lake-side* rim happens to be the higher
# of the two boundary nodes (confirmed directly: a small synthetic basin reproduces this,
# roughly a coin flip which side is higher) never received one erosive drop from this
# mechanism, and so never wore down -- the lake could spill for the rest of the run without
# its actual cap (`max_depth`, recomputed fresh from `elevation` every step -- see lakes.py's
# own docstring) ever dropping. A real spillway erodes right at the pinch point, whichever
# side of it happens to be taller; injecting the same source at both nodes erodes the true
# bottleneck regardless of which one that is, which is what actually lets a long-spilling lake
# grind its own rim down and, eventually, drain.
LAKE_BREACH_EROSION_COEFFICIENT = 4000.0

# Glaciers use flat per-Myr rates (linear in dt_myr), not an exponential-decay formula like
# the lake evaporation rate above -- so they don't blow up the same way at mantle-bloom's
# larger typical step size -- a bigger step just melts/grows proportionally more, same
# character as every other elevation delta in this codebase (boundary.py's uplift/trench
# rates are also flat per-Myr rates applied directly).
#
# A step spans far more than a year, so seasons aren't stepped through, but their effect on
# melt is imputed. Each node's annual cycle is taken as a sinusoid of its annual-mean
# temperature with a half-amplitude synthesized from latitude, continentality, relief and
# axial tilt (`biomes.seasonal_temp_amplitude`, passed in by erosion.py), and surface melt
# scales with that cycle's mean positive degrees above GLACIER_MELT_THRESHOLD_C
# (`seasonal_positive_degrees`) -- a positive-degree-day melt model. A valley whose annual
# mean is -8C but whose summer peaks at +7C melts every summer; a maritime cell at the same
# -8C with a 4C swing never thaws. With zero amplitude (a tilt-0 world, or a caller that
# passes none) this reduces to melt above an annual mean of 0C.
#
# Melt still competes with accumulation every step, so the ice margin settles where a
# summer's positive degrees just balance that node's snowfall -- for typical precipitation,
# where the warmest-month mean is a degree or two above freezing (about an annual mean of
# -12.5C under a 15C half-amplitude), close to where real glacier equilibrium lines sit.
#
# FREEZE_POINT_C is the real phase-change point: precipitation below it falls as snow/ice
# rather than rain, and standing lake water or flowing river water below it freezes solid
# *this same step* (see compute_hydrology). That stays annual-mean-gated; a node between the
# two (annual mean below 0C, summer above it) freezes its lake/river for the step but the ice
# it forms still melts back every step, so it doesn't build a permanent glacier.
FREEZE_POINT_C = 0.0
GLACIER_MELT_THRESHOLD_C = 0.0
GLACIER_ACCUMULATION_RATE = 0.02
GLACIER_MELT_RATE_M_PER_MYR = 400.0
GLACIER_MELT_REFERENCE_DEGREES_C = 20.0
GLACIER_MELT_MAX_FACTOR = 3.0
GLACIER_VISIBLE_DEPTH_M = 10.0

# A hard ceiling on per-node ice depth, kept as a guard (issue #353). Ice is shaped by
# ice_sheet.py's plastic profile now, whose divide thickness is sqrt(2 c L) for a sheet whose
# margin is L from its divide: 3.3 km at L = 500 km and 4.7 km at L = 1000 km under 100 kPa,
# so only sheets wider than about 1100 km reach it. The excess is shed into the step's
# meltwater (mass-conserving into the river system, not deleted). ~1.5x the thickest ice
# sheet on real Earth.
GLACIER_MAX_DEPTH_M = 5000.0

# During a deep enough ice age -- climate.ice_age_glacial_intensity at or above
# ICE_AGE_SEA_ICE_MIN_INTENSITY -- ice is allowed to pile up on, and persist over, genuinely
# freezing polar *ocean* nodes (surface temperature below SEA_ICE_FORMATION_TEMP_C) rather
# than being discarded at the coastline the way it is the rest of the time (see
# ice_sheet.step_ice's calving of floating ice, gated by an `allow_ocean_ice` mask). This is the "ice
# caps can grow over polar water in an ice age" behaviour: land glaciers calving toward the
# coast, plus frozen precipitation landing on the sea itself, build a floating ice cap over
# the pole. Outside an ice age -- or during its interglacial phase, intensity near 0 -- the
# gate stays shut and the old calve-at-the-coast rule holds exactly (bit-identical). The sea
# ice still melts back through the ordinary melt_factor term, and once the glacial phase
# passes the gate closes and that step discards whatever had accumulated (shed to meltwater,
# mass-conserving into the river/ocean system). SEA_ICE_FORMATION_TEMP_C matches biomes.py's
# own `sst < -1.0` sea-ice cutoff; the intensity floor is half the raised-cosine amplitude,
# i.e. "past the midpoint of a glacial build".
SEA_ICE_FORMATION_TEMP_C = -1.0
ICE_AGE_SEA_ICE_MIN_INTENSITY = 0.5

# River evaporation: a fraction of a river's own downstream flux evaporates at every node it
# passes through this step -- warmer nodes lose a bigger fraction, up to
# RIVER_EVAPORATION_MAX_FRACTION (capped well short of 1.0 so a single very large step can't
# evaporate a river's entire flow in one hop). This is the "evaporation decreases the amount of
# water in rivers" half of the moisture-recycling ask; the *other* half -- that lost water
# actually reappearing as atmospheric humidity -- lives in climate.py's own land-surface
# moisture source, not here: climate.py runs strictly *before* this module each step (it
# produces the precipitation/temperature this module consumes), so this step's real
# evaporation loss can't feed back into this same step's climate. climate.py's own term is a
# same-step stand-in sized from the *persisted* channel_depth instead (see that module's own
# docstring for why that isn't circular). Flat per-Myr rate, the same "linear in dt_myr, no
# exponential blow-up at a large step" character GLACIER_ACCUMULATION_RATE above already has.
RIVER_EVAPORATION_REFERENCE_TEMP_C = 25.0
RIVER_EVAPORATION_RATE_PER_MYR = 0.2
RIVER_EVAPORATION_MAX_FRACTION = 0.6

# River speed -- a stylized, unitless quantity (faster where slope is steeper and discharge is
# larger); discharge is first mapped to the default-density reference scale below. Lives here
# rather than erosion.py so the River Inspector's per-river mouth speed can share it without a
# hydrology.py -> erosion.py reverse import.
RIVER_SPEED_COEFFICIENT = 4.0
RIVER_SPEED_DISCHARGE_EXPONENT = 0.2

# Issue #297 part 1: least-cost depression breaching (see breaching.py). When on, each step's
# depression hierarchy uses notched passage elevations: established channels plus the notches
# this step invents through pits cheap enough to carve. Pits too costly to carve stay closed
# and are filled by lakes.py as before.
BREACH_DEPRESSIONS_ENABLED = True


@dataclass
class HydrologyFields:
    """Everything derived from one whole-world flow-routing pass, all shape (N,) aligned
    with `points`/`elevation` -- one flat array per field over the irregular node cloud, in
    place of a grid. `plates_in_order` is every plate that contributed nodes, in the same
    order those nodes appear in every array here -- see plates.gather_node_positions --
    letting a caller (geology.py) bulk-collect/write-back further per-node fields without
    needing per-node plate/line identity itself. `lake_depth`/`glacier_depth` are already this
    step's *final* values (state-transition owned by this module, see module docstring) -- a
    caller doesn't need to call anything else to get the up-to-date depth."""

    points: np.ndarray
    elevation: np.ndarray
    is_ocean: np.ndarray
    neighbor_idx: np.ndarray  # (N, FLOW_NEIGHBOR_COUNT)
    flow_target: np.ndarray  # (N,) int, index into these same arrays; -1 = sink or ocean
    flow_accum: np.ndarray  # (N,) legacy stacked-mm/yr routing, fallback for caches without discharge_m3_s
    water_deposited: np.ndarray  # (N,) legacy water settled here, same fallback
    filled_elevation: np.ndarray  # (N,) minimal-bottleneck elevation to reach open ocean
    spill_target: np.ndarray  # (N,) one-hop basin-escape neighbor index; -1 = none/ocean
    is_river: np.ndarray  # (N,) bool -- top RIVER_FLOW_PERCENTILE of land discharge_m3_s
    lake_depth: np.ndarray  # (N,) this step's final standing lake depth
    glacier_depth: np.ndarray  # (N,) this step's final accumulated ice depth
    plates_in_order: list[Plate]
    # This step's lake merge/split transitions from lakes.step_lakes, as structured
    # `lakes.LakeEvent`s. erosion.py runs them through `lakes.summarize_lake_events` (which
    # collapses the near-sea-level coastal-pond churn a dithering low-relief shelf produces)
    # and logs the resulting lines via World.log_event -- see compute_hydrology. Defaulted
    # (rather than a required positional field) so every existing HydrologyFields call site
    # that predates lakes.py, including hand-built test fixtures, keeps working unchanged.
    lake_events: list["lakes.LakeEvent"] = field(default_factory=list)
    # Cumulative silt thickness laid down under standing water, end of this step
    # (`prev_silt_depth + silt_deposited`), persisted on plates -- purely an informational
    # record: lake floors are not measured against it, because `silt_deposited` below is folded
    # straight into real `elevation` by erosion.py, as mobile cover. erosion.py caps what it
    # persists at the cell's mobile cover, so stripped silt leaves the record. Also defaulted,
    # see lake_events' own comment for why.
    silt_depth: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # This step's per-node silt *increment* (lakes.step_lakes) -- erosion.py adds it directly
    # into terrain `elevation`, so a lake/endorheic basin genuinely fills in over time. Always
    # >= 0, nonzero only under standing water this step. Defaulted, same reasoning as above.
    silt_deposited: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # This step's fully-resolved depression hierarchy (lakes.step_lakes's own internal
    # `lakes.build_lake_hierarchy` result, built from this step's `elevation` and already
    # reconciled against last step's lake_depth -- see step_lakes's own docstring), exposed
    # here rather than discarded so a caller needing basin floor/outlet/current-water info
    # (main.py's Lake Inspector) can read it straight off these already-resolved `Lake` objects
    # instead of rebuilding its own hierarchy -- which is more than a wasted computation:
    # depression-hierarchy topology is sensitive enough to small elevation shifts (e.g.
    # rebuilding *after* erosion has folded in this step's silt/deposition) that the two
    # rebuilds can genuinely disagree on which basin merged with which -- confirmed directly as
    # a real bug this way: a lake whose reported `outlet_elevation_m` sat *below* its own
    # `floor_elevation_m`. Also defaulted, same backward-compatibility reasoning as
    # `lake_events`/`silt_depth` above.
    lake_forest: list["lakes.Lake"] = field(default_factory=list)
    # The `CellHypsometry` `lake_forest` was resolved against, kept for the same reason: a
    # lake's lowest ground is `hypsometry.dry_level(members)`, below `Lake.floor_elevation`
    # (the lowest cell *mean*) by the floor cell's sub-cell relief, and a lake wet only in that
    # low ground has its surface below `floor_elevation` (issue #344). None until hydrology
    # has run with it. Defaulted, same reasoning as `lake_events` above.
    hypsometry: "lakes.CellHypsometry | None" = None
    # Ice's own downhill flow target -- like flow_target, but never zeroed by is_frozen: glacial
    # ice flows via gravity/internal deformation regardless of whether this particular step is
    # literally below freezing, unlike liquid water/rivers (see _compute_flow_direction's
    # apply_freeze parameter). erosion.py routes glacially-eroded till along this, not
    # flow_target, so it travels with the ice's real flow path out to the glacier's melting
    # margin (a terminal moraine/outwash deposit) instead of being stranded wherever it was
    # scoured. Also defaulted, same backward-compatibility reasoning as lake_events above.
    ice_flow_target: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    # True wherever lake_depth is standing water AND lakes._classify_tier currently reads that
    # body as sea tier (large enough, see lakes.SEA_MIN_FLOODED_AREA_KM2) rather than lake tier
    # -- climate.py's own moisture source and render_image.py's overlay both key off this the
    # same way they already key off lake_depth. Also defaulted, same backward-compatibility
    # reasoning as lake_events above -- a cache loaded from a save written before this field
    # existed has none of it, so every reader treats a shape-0/mismatched array as "no sea
    # anywhere" via the same guard sample_is_ocean already uses for a stale/absent cache.
    is_sea: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    # Prior-step channel_width (erosion.py still owns growing it -- same one-step-stale
    # tolerance as channel_depth elsewhere in this module), aligned with points/elevation.
    # render_image.py's `_draw_rivers` reads this directly so it stays a pure function of
    # HydrologyFields alone, rather than re-deriving it from world.plates the way
    # render_image.py's node-cloud-aligned fields (elevation, channel_depth for hillshade) do --
    # unlike those, HydrologyFields is unit-tested via hand-built fixtures fully decoupled from
    # any real Plate. Also defaulted, same backward-compatibility reasoning as is_sea above --
    # a cache loaded from a save written before this field existed has none of it, so a reader
    # detects the shape-0/mismatched array as "stale/absent" the same way is_sea's own guard does.
    channel_width: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # This step's new breach notches (breaching.breach_depressions), in metres. erosion.py adds
    # them to the persisted breach_notch_depth_m so that the next step's passage elevation
    # keeps them. Empty when breaching is off. Defaulted, for the same backward-compatibility
    # reason as is_sea.
    breach_notch_m: np.ndarray = field(default_factory=lambda: np.zeros(0))
    # This step's glacier solve (issue #353): the ice ledger, per-body budgets, plastic
    # thickness, thermal state, balance flux and sliding speed erosion.py abrades with.
    # `glacier_depth` above already includes it. A plain None default (a class attribute),
    # so a cache pickled before this field existed still reads as "no ice diagnostics".
    ice: "IceStep | None" = None
    # Physical river discharge (issue #355, see "Water units" below): land-node through-flow
    # and ocean-mouth inflow of the m^3/yr runoff-plus-meltwater routing, in m^3/s. Plain None
    # default, like `ice`, so a cache pickled before this existed still loads.
    discharge_m3_s: np.ndarray | None = None
    # Nodes whose `flow_target` is a spill rather than steepest descent: a filled basin's sink
    # routed to its outlet, or a united lake's sinks pooled at its lowest one. Such an edge can
    # span a whole lake, so river drawing follows the spill path instead (`flow_segments`,
    # #520). Plain None default, like `ice`, so a cache pickled before this existed still loads.
    should_spill: np.ndarray | None = None
    # This step's breach-surge sources (`lake_breach_source + rim_breach_source`), in model
    # units. They're kept out of `discharge_m3_s`, but they still reach every consumer of the
    # legacy `flow_accum`/`water_deposited` routing (see "Water units").
    breach_surge: np.ndarray | None = None
    # Routed breach forcing in model units, aligned with the drainage graph. River erosion
    # consumes it separately from physical discharge; other legacy consumers retain the
    # combined routing until their own migration. None keeps older pickled caches readable.
    breach_surge_accum: np.ndarray | None = None

    def __getstate__(self) -> dict:
        # The spill-path drawing cache (`_spill_path_cache`) is rebuilt on first use, so a save
        # never carries paths from older drawing rules.
        state = self.__dict__.copy()
        state.pop("_spill_path_cache", None)
        return state


def _gather_nodes(
    world: "World",
    node_cloud: tuple[np.ndarray, list[Plate]] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[Plate]]:
    """Every node's world position, elevation, prior lake_depth, prior glacier_depth, prior
    channel_depth, prior channel_width, prior silt_depth, and whether it's ocean (the bare
    `elevation <= world.sea_level_m` test -- `compute_hydrology` refines this to a
    connectivity-aware mask once the k-NN graph exists, see `connected_ocean_mask`),
    concatenated, alongside the ordered list of plates that contributed
    them -- same shape as erosion.py's/bathymetry.py's own _gather_nodes. channel_depth/
    channel_width are read-only here (erosion.py still owns growing both) -- flow direction
    just needs to know where an established channel already is (see
    _compute_flow_direction), and channel_width is carried through purely so
    `HydrologyFields.channel_width` (render_image.py's river-color-blend input) doesn't need
    its own separate world.plates read. `node_cloud`, when
    passed (see compute_hydrology), reuses an already-gathered (points, plates_in_order) pair
    instead of re-deriving every node's world position from scratch -- see
    plates.gather_node_positions's own docstring for why."""
    points, plates_in_order = node_cloud if node_cloud is not None else gather_node_positions(world.plates)
    if not plates_in_order:
        empty = np.zeros(0)
        return np.zeros((0, 3)), empty, empty, empty, empty, empty, empty, np.zeros(0, dtype=bool), []
    elevation = collect_all_elevation(plates_in_order)
    prev_lake_depth = collect_all_lake_depth(plates_in_order)
    prev_glacier_depth = collect_all_glacier_depth(plates_in_order)
    prev_channel_depth = collect_all_channel_depth(plates_in_order)
    prev_channel_width = collect_all_channel_width(plates_in_order)
    # The silt record describes lake silt still in the mobile cover (see erosion.py), so it
    # never reads more than the cover holds -- tectonics can strip cover between steps.
    prev_silt_depth = np.minimum(
        collect_all_silt_depth(plates_in_order),
        np.concatenate([p.collect("mobile_cover_m") for p in plates_in_order]),
    )
    is_ocean = elevation <= world.sea_level_m
    return (
        points, elevation, prev_lake_depth, prev_glacier_depth, prev_channel_depth, prev_channel_width,
        prev_silt_depth, is_ocean, plates_in_order,
    )


def _build_neighbor_graph(points: np.ndarray, world: "World | None" = None) -> np.ndarray:
    n = len(points)
    k = min(FLOW_NEIGHBOR_COUNT, n - 1)
    # `world`, when passed (compute_hydrology's own per-step call does), shares this tree via
    # plates.cached_node_position_tree with every other per-step full-node-cloud query
    # (climate.py's own resample, erosion.py's slope/earthquake/wind-deposit queries) instead
    # of rebuilding an equivalent one from scratch -- see that function's own docstring.
    # `world=None` (connected_ocean_mask's own fallback build, and direct unit-test calls)
    # keeps the old balanced_tree=False/compact_nodes=False build -- same build-once/
    # query-once tradeoff as a single-use per-plate tree: built fresh
    # and queried exactly once (one batched k-NN query over every node at once). Either way,
    # an exact k-nearest-neighbor search -- results unchanged.
    tree = (
        cached_node_position_tree(world, points)
        if world is not None
        else cKDTree(points, balanced_tree=False, compact_nodes=False)
    )
    _, neighbor_idx = tree.query(points, k=k + 1, workers=query_workers(n))
    return neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0


def _below_sea_components(below: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Connected-component label per node over the `neighbor_idx` edges joining two
    below-sea-level nodes. Shared by `connected_ocean_mask` and
    `force_earth_strait_connections`, which must agree on the components."""
    n, k = neighbor_idx.shape
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    keep = below[rows] & below[cols]
    graph = coo_matrix((np.ones(keep.sum()), (rows[keep], cols[keep])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def connected_ocean_mask(
    points: np.ndarray,
    elevation: np.ndarray,
    sea_level_m: float,
    neighbor_idx: np.ndarray | None = None,
) -> np.ndarray:
    """Per-node bool, True only where a node is both below sea level *and* part of the world
    ocean -- the single largest connected below-sea-level body, plus any other component within
    `SECOND_OCEAN_SIZE_FRACTION` of its size (see that constant's own comment for why: without
    it, a world with two genuinely separate, ocean-scale bodies of water not only misclassifies
    the smaller as a giant lake, but flickers it in and out of the lake list every few steps as
    the strait between them opens/closes). Found by a connected-components pass
    (`scipy.sparse.csgraph`, same tool `plates.node_components` uses) over the k-NN graph
    restricted to below-sea-level nodes. Every *other* below-sea-level region -- an enclosed
    interior depression that erosion or plate deformation dropped beneath sea level without
    ever connecting it to open water -- comes back False, so hydrology treats it as an ordinary
    endorheic basin (rivers pool there and silt it in, see lakes.py) and biome/climate code
    stops painting it as shallow ocean.

    `neighbor_idx` (this module's own FLOW_NEIGHBOR_COUNT graph) is reused when passed, else a
    fresh graph is built from `points`. Returns the bare `elevation <= sea_level_m` mask
    unchanged when there are too few nodes to build a graph, or when every below-sea-level node
    is already one connected body -- the ordinary single-ocean world, the common case."""
    below = elevation <= sea_level_m
    n = len(elevation)
    if n == 0 or not below.any() or below.all():
        return below
    if neighbor_idx is None:
        if n <= FLOW_NEIGHBOR_COUNT:
            return below
        neighbor_idx = _build_neighbor_graph(points)

    labels = _below_sea_components(below, neighbor_idx)

    below_idx = np.nonzero(below)[0]
    below_labels = labels[below_idx]
    if np.unique(below_labels).size <= 1:
        return below  # every below-sea-level node is one connected ocean already
    counts = np.bincount(below_labels)
    largest_label = counts.argmax()
    largest_count = counts[largest_label]
    # The largest component is always ocean, full stop, regardless of the floor below -- that
    # floor exists only to gate whether some *other* component also counts as ocean, not to
    # second-guess the single biggest below-sea-level body itself (which would otherwise
    # wrongly demote the whole ocean to "no ocean at all" on a small enough world, toy test
    # fixtures included).
    qualifies = (counts >= largest_count * SECOND_OCEAN_SIZE_FRACTION) & (counts >= SECOND_OCEAN_MIN_NODE_COUNT)
    qualifies[largest_label] = True
    ocean_labels = np.nonzero(qualifies)[0]
    return below & np.isin(labels, ocean_labels)


# (lon, lat) degrees on each side of a real strait narrower than this engine can ever resolve
# on its own -- see real_plates.py's own module comment for why Gibraltar specifically is
# narrower than even the "earth" premade world's exact coastline raster, a limit of the
# source cartographic data itself, not a bug in that bake. Only consulted for that one
# premade world (see `force_earth_strait_connections`'s caller) -- a coordinate this specific
# is never meaningful for a hand-drawn or procedurally generated world.
_FORCED_EARTH_STRAITS_LONLAT = [
    ((-6.5, 35.9), (-3.5, 35.9)),  # Gibraltar: just inside the Atlantic, just inside the Mediterranean
]


def force_earth_strait_connections(
    points: np.ndarray, is_ocean: np.ndarray, elevation: np.ndarray, sea_level_m: float, neighbor_idx: np.ndarray
) -> np.ndarray:
    """`is_ocean` (`connected_ocean_mask`'s own result) with every strait in
    `_FORCED_EARTH_STRAITS_LONLAT` counted as connected regardless of whether any node
    actually sits in the gap -- the "earth" premade world's one deliberate exception to
    letting node placement alone decide ocean connectivity (its caller gates this to that one
    premade world). Recomputes the same below-sea-level connected-components labeling
    `connected_ocean_mask` already did (cheap: the k-NN graph is already built) purely to find
    which component each strait's two sides falls in; a no-op for a strait where either side
    isn't actually below sea level this step (ice age glaciation, or plate motion has since
    pushed that coordinate above water), or where they're already connected (a node happens to
    sit in the gap; ordinary `connected_ocean_mask` output already has this right)."""
    n = len(elevation)
    below = elevation <= sea_level_m
    if n == 0 or not np.any(below):
        return is_ocean
    labels = _below_sea_components(below, neighbor_idx)

    tree = cKDTree(points)
    for (lon_a, lat_a), (lon_b, lat_b) in _FORCED_EARTH_STRAITS_LONLAT:
        side_xyz = geometry.latlon_to_xyz(np.radians([lat_a, lat_b]), np.radians([lon_a, lon_b]))
        _, (idx_a, idx_b) = tree.query(side_xyz)
        if not (below[idx_a] and below[idx_b]) or labels[idx_a] == labels[idx_b]:
            continue
        if is_ocean[idx_a] or is_ocean[idx_b]:
            is_ocean = is_ocean | (labels == labels[idx_a]) | (labels == labels[idx_b])
    return is_ocean


# Holds the most recent (hydro, query_xyz, idx, resample_mode) from either
# `sample_is_ocean`/`sample_is_sea` below, so a call to one immediately followed by a call to
# the other *for the same (hydrology cache, query points, resample mode)* -- climate.py's
# `_sample_elevation_and_crust` and render_image.py's `_biome_fields` both do exactly this,
# back-to-back, against the same `world_xyz` -- reuses the first call's tree/index build+query
# instead of repeating an identical one just to read out a different boolean field. Keyed by
# real object identity (`is`, not `id()`, which CPython can recycle once an object is freed)
# held in the tuple itself, so a mismatch can only ever cause an extra (still-correct) rebuild,
# never a stale hit. `resample_mode` (issue #133 phase 2) joins the key alongside identity, not
# just to pick the source below, but because it's plain data -- unlike `hydro`/`query_xyz`,
# comparing it by value rather than `is` is correct and cheap.
_LAST_NEAREST_HYDRO_NODE: tuple[object, object, np.ndarray, str] | None = None

# The cKDTree half of the cache above, held separately (and keyed on `hydro` alone, not also
# `query_xyz`) because a *different* query_xyz this same step -- the climate grid, then later
# the render grid, both resampling the same `hydro.points` -- still shouldn't pay a fresh
# ~130K-node tree build just because the point above's exact-match cache misses on
# `query_xyz`. Not stored as a field on `HydrologyFields` itself: `World.hydrology_cache` is
# part of the pickled save-file object graph (see persistence.py), so a tree living there
# would round-trip through every save/load for no benefit -- a plain module global, like
# `_LAST_NEAREST_HYDRO_NODE` above, never enters that graph at all. Same identity-keyed,
# worst-case-is-an-extra-rebuild safety as that cache.
_LAST_OCEAN_TREE: tuple[object, cKDTree] | None = None

# The "healpix" World.node_cloud_resample_mode analogue of _LAST_OCEAN_TREE above (issue #133
# phase 2) -- same module-global, not-on-World placement (and the same reasoning: a save/load
# round-trip through World.hydrology_cache would carry a HealpixGrid/NodePixelIndex for no
# benefit), split the same way plates.cached_node_healpix_index splits World's own pair: the
# HealpixGrid itself (keyed by nside, a pure function of hydro.points' node *count*) versus the
# per-hydro-cache scatter+fill population, since hydro.points only actually changes identity
# once per step_world (a new HydrologyFields object), the same cadence _LAST_OCEAN_TREE reuses
# at today.
_LAST_OCEAN_HEALPIX_GRID: tuple[int, healpix_grid.HealpixGrid] | None = None
_LAST_OCEAN_HEALPIX_INDEX: tuple[object, healpix_grid.NodePixelIndex] | None = None

# Issue #133 phase 4: phase 2 measured "healpix" mode's is_ocean disagreeing with "kdtree" in a
# dilated/eroded coastal band 16.1% of the time (vs. 0.03% interior) at the default 1:1
# node-count sizing -- real, not a formality, since this is exactly where is_ocean correctness
# matters most. A finer ocean grid was one of the two levers that write-up flagged (the other,
# an is_ocean-aware scatter tie-break, is unimplemented). Measured directly (throwaway spike,
# same seed=0/node_density=4.0/4-step world the phase 2/3 numbers above used): 2x -> 11.4%
# coastal, 4x -> 7.2%, 8x -> 5.4% -- diminishing returns, and 8x's one-time HealpixGrid.build()
# cost (~67s at this node count) makes it impractical. 4x was chosen as the best cost/accuracy
# trade: ~55% relative reduction in coastal mismatch for a per-`step_world` fill cost of
# ~65-90ms (vs. ~24ms at 1x) against a step that otherwise takes several *seconds* -- and the
# grid itself (unlike the fill) is built once per node-count bracket, not every step. Still not
# zero, so this alone doesn't clear the phase-4 default-flip gate; see docs/profiling.md's
# "Phase 4" section for the full numbers and that decision.
_OCEAN_NSIDE_OVERSAMPLE = 4


def _ocean_node_tree(hydro: "HydrologyFields") -> cKDTree:
    global _LAST_OCEAN_TREE
    cached = _LAST_OCEAN_TREE
    if cached is not None and cached[0] is hydro:
        return cached[1]
    tree = cKDTree(hydro.points)
    _LAST_OCEAN_TREE = (hydro, tree)
    return tree


def _ocean_node_healpix_index(hydro: "HydrologyFields") -> healpix_grid.NodePixelIndex:
    """The `"healpix"` counterpart to `_ocean_node_tree` immediately above -- see
    `_LAST_OCEAN_HEALPIX_GRID`/`_LAST_OCEAN_HEALPIX_INDEX`'s own comment for the caching split."""
    global _LAST_OCEAN_HEALPIX_GRID, _LAST_OCEAN_HEALPIX_INDEX
    cached_index = _LAST_OCEAN_HEALPIX_INDEX
    if cached_index is not None and cached_index[0] is hydro:
        return cached_index[1]
    nside = healpix_grid.nside_for_node_count(len(hydro.points), oversample=_OCEAN_NSIDE_OVERSAMPLE)
    cached_grid = _LAST_OCEAN_HEALPIX_GRID
    if cached_grid is not None and cached_grid[0] == nside:
        grid = cached_grid[1]
    else:
        grid = healpix_grid.build(nside)
        _LAST_OCEAN_HEALPIX_GRID = (nside, grid)
    index = healpix_grid.build_node_pixel_index(grid, hydro.points)
    _LAST_OCEAN_HEALPIX_INDEX = (hydro, index)
    return index


def _nearest_hydro_node_idx(hydro: "HydrologyFields", query_xyz: np.ndarray, resample_mode: str = "kdtree") -> np.ndarray:
    """Nearest-node index from last step's `hydro.points` for each `query_xyz` point --
    `resample_mode` (`World.node_cloud_resample_mode`, issue #133 phase 2) picks `_ocean_node_tree`
    (`cKDTree`, default) or `_ocean_node_healpix_index` (`NodePixelIndex`); folded into the
    single-entry cache key below alongside `hydro`/`query_xyz` identity so a mode flip
    mid-process (e.g. a test flipping `World.node_cloud_resample_mode`) can't serve a stale
    cross-mode hit."""
    global _LAST_NEAREST_HYDRO_NODE
    cached = _LAST_NEAREST_HYDRO_NODE
    if cached is not None and cached[0] is hydro and cached[1] is query_xyz and cached[3] == resample_mode:
        return cached[2]
    flat_xyz = np.asarray(query_xyz).reshape(-1, 3)
    source = _ocean_node_healpix_index(hydro) if resample_mode == "healpix" else _ocean_node_tree(hydro)
    _, idx = source.query(flat_xyz, workers=query_workers(len(flat_xyz)))
    _LAST_NEAREST_HYDRO_NODE = (hydro, query_xyz, idx, resample_mode)
    return idx


def sample_is_ocean(world: "World", query_xyz: np.ndarray, fallback_is_ocean: np.ndarray) -> np.ndarray:
    """Connectivity-aware `is_ocean` (see `connected_ocean_mask`) resampled onto arbitrary
    `query_xyz` (any shape, trailing axis 3) from last step's `world.hydrology_cache` --
    nearest-node, the same self-contained cKDTree resample `coastline._lake_mask_on_grid` does
    for lake_depth (a fresh query here, not the caller's own `idx`, so it stays correct even
    when plate topology changed between that cache and this call). `fallback_is_ocean` (a bare
    `elevation <= sea_level_m` mask the caller already has) is returned unchanged before the
    first step, when no hydrology cache exists yet; the result otherwise has
    `fallback_is_ocean`'s shape. Resamples through `_nearest_hydro_node_idx` under
    `world.node_cloud_resample_mode` (issue #133 phase 2) the same way every other node-cloud
    resample in this codebase now does."""
    hydro = getattr(world, "hydrology_cache", None)
    if hydro is None or len(hydro.points) == 0 or len(hydro.is_ocean) != len(hydro.points):
        return fallback_is_ocean
    idx = _nearest_hydro_node_idx(hydro, query_xyz, getattr(world, "node_cloud_resample_mode", "kdtree"))
    return hydro.is_ocean[idx].reshape(fallback_is_ocean.shape)


def sample_is_sea(world: "World", query_xyz: np.ndarray, fallback_is_sea: np.ndarray) -> np.ndarray:
    """Same resample-from-cache contract as `sample_is_ocean` immediately above, for
    `HydrologyFields.is_sea` (lakes.py's sea-tier classification, see `lakes._classify_tier`):
    climate.py's own moisture source and render_image.py's overlay both need this on their own
    grid the same way they already need `is_ocean`/`lake_depth`. `fallback_is_sea` (all-False
    -- nothing is a sea before the world's first hydrology pass) is returned unchanged when no
    cache exists yet, or when a loaded cache predates this field (`HydrologyFields.is_sea`'s own
    comment -- an old save's cached object simply has no such attribute at all)."""
    hydro = getattr(world, "hydrology_cache", None)
    hydro_is_sea = getattr(hydro, "is_sea", None) if hydro is not None else None
    if hydro is None or len(hydro.points) == 0 or hydro_is_sea is None or len(hydro_is_sea) != len(hydro.points):
        return fallback_is_sea
    idx = _nearest_hydro_node_idx(hydro, query_xyz, getattr(world, "node_cloud_resample_mode", "kdtree"))
    return hydro_is_sea[idx].reshape(fallback_is_sea.shape)


def lake_components(is_lake: np.ndarray, neighbor_idx: np.ndarray) -> list[np.ndarray]:
    """Every connected `is_lake` component (as an array of member node indices) via union-find
    over `neighbor_idx`'s edges, treated as **undirected** here even though the k-NN graph
    itself isn't strictly symmetric (a node isn't guaranteed to be among its own neighbor's k
    nearest). Unlike `_flood_fill_lake_extent`'s directed traversal (fine for its own purpose,
    spreading a water *level* outward from a seed), a same-lake grouping needs two lake nodes
    joined by only a one-directional edge to still read as the same lake, not be split apart by
    an accident of which node happens to list the other as one of its own nearest neighbors.
    Shared building block for `_lake_component_sizes` (below) and main.py's Lake Inspector
    (`GET /world/lakes`/`GET /world/lake_at`), which needs the actual per-lake node groups, not
    just their sizes."""
    n = len(is_lake)
    if not np.any(is_lake):
        return []

    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    lake_list = is_lake.tolist()
    neighbor_list = neighbor_idx.tolist()
    for i in range(n):
        if not lake_list[i]:
            continue
        for j in neighbor_list[i]:
            if lake_list[j]:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj

    lake_idx = np.nonzero(is_lake)[0]
    groups: dict[int, list[int]] = {}
    for i in lake_idx.tolist():
        groups.setdefault(find(i), []).append(i)
    return [np.array(members, dtype=np.int64) for members in groups.values()]


def _lake_component_sizes(is_lake: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Per-node size (node count) of the connected `is_lake` component containing that node, 0
    for every non-lake node -- see `lake_components`, which this just tallies."""
    n = len(is_lake)
    sizes = np.zeros(n, dtype=np.int64)
    for members in lake_components(is_lake, neighbor_idx):
        sizes[members] = len(members)
    return sizes


# In-transit river evaporation on real-unit (m^3/yr) fluxes -- the breach gate and the lake
# inflow -- is a property of the river, not of the step. Those paths take the per-node loss
# at this reference step length (`river_transit_loss_fraction`). The model-unit routing of
# `water_source` keeps `river_evaporation_fraction`'s step-scaled rule.
RIVER_TRANSIT_LOSS_REFERENCE_YEARS = 100_000.0


def river_evaporation_fraction(temperature_c: np.ndarray, years: float, excluded: np.ndarray) -> np.ndarray:
    """Share of the water passing through each node lost to in-transit river evaporation this
    step (see RIVER_EVAPORATION_*), zero where `excluded`."""
    years_myr = years / 1_000_000.0
    fraction = np.clip(temperature_c / RIVER_EVAPORATION_REFERENCE_TEMP_C, 0.0, 1.0) * RIVER_EVAPORATION_RATE_PER_MYR * years_myr
    return np.where(excluded, 0.0, np.clip(fraction, 0.0, RIVER_EVAPORATION_MAX_FRACTION))


def river_transit_loss_fraction(temperature_c: np.ndarray, excluded: np.ndarray) -> np.ndarray:
    """Per-node in-transit evaporation for real-unit fluxes, independent of step length (see
    RIVER_TRANSIT_LOSS_REFERENCE_YEARS)."""
    return river_evaporation_fraction(temperature_c, RIVER_TRANSIT_LOSS_REFERENCE_YEARS, excluded)


# Water units (issue #355). Three kinds of water move along `flow_target`:
# - Physical water, in m^3/yr: area-weighted Budyko runoff of liquid precipitation
#   (`breaching.runoff_mm`) plus glacier meltwater (`meltwater_m3_per_yr`). It balances the
#   lakes, and its through-flow is the river discharge, reported in m^3/s
#   (`HydrologyFields.discharge_m3_s`).
# - Glacier meltwater is the #353 ice ledger's per-step volume (`IceLedger.meltwater_m3`),
#   spread over the step. It covers surface and basal melt, overflow, ice shed at the depth
#   guard, and snow that melts on bare ground. Ice and water volumes are booked 1:1, as a
#   freezing lake is. Calved ice goes straight to the ocean, and accumulation and frozen river
#   ice stay ice, so none of them is discharge. The ledger is still in #353's model units, so
#   meltwater is small next to runoff. Rescaling glacier mass balance changes discharge only
#   through the ledger.
# - Frozen river ice is the physical water that settles at a frozen land node, converted to
#   ice at the glacier model's own exchange rate (`WATER_TO_ICE`), one step late.
# - The breach surge (`HydrologyFields.breach_surge`) is a spilling lake's erosive forcing.
#   It's sized by node count (LAKE_BREACH_EROSION_COEFFICIENT) and is in model units. It is
#   never physical discharge.
# `flow_accum`/`water_deposited` are the legacy routing: stacked mm/yr of liquid
# precipitation, melt and breach surge, with step-scaled evaporation. Step 3 moved erosion,
# width, speed, Inspector flow values and the geology river-mouth boost to physical discharge
# (erosion also gets separately-routed `breach_surge_accum`), and lake silt reads the physical
# lake inflow (#363). No live consumer reads the legacy fields now: they remain only as the
# fallback for hand-built fixtures and caches without `discharge_m3_s`.
SECONDS_PER_YEAR = 3.15576e7
# A reference cell area maps physical catchment discharge back to the old stacked-rainfall
# scale at the default node density. Keeping the old scalar calibration on this fixed area
# preserves its default-density meaning while making a fixed catchment independent of the
# actual node density. The square-spacing estimate is intentionally a reference, not a claim
# about each irregular node's true accounting area.
RIVER_REFERENCE_CELL_AREA_M2 = (PLANET_RADIUS_KM * 1000.0 * line_spacing_rad(DEFAULT_NODE_DENSITY)) ** 2


def discharge_m3_s_to_flow_accum_mm_yr(discharge_m3_s: np.ndarray) -> np.ndarray:
    """Map physical discharge to the legacy mm/yr-equivalent scale at default node density."""
    return np.clip(np.asarray(discharge_m3_s, dtype=float), 0.0, None) * SECONDS_PER_YEAR * 1000.0 / RIVER_REFERENCE_CELL_AREA_M2


def flow_accum_mm_yr_to_discharge_m3_s(flow_accum_mm_yr: np.ndarray | float) -> np.ndarray:
    """Convert the legacy accumulated mm/yr scalar to discharge at the reference cell area."""
    return np.clip(np.asarray(flow_accum_mm_yr, dtype=float), 0.0, None) * RIVER_REFERENCE_CELL_AREA_M2 / (SECONDS_PER_YEAR * 1000.0)


def _extend_river_mask_downstream(
    is_river: np.ndarray,
    flow_target: np.ndarray,
    discharge_m3_s: np.ndarray,
    is_ocean: np.ndarray,
    is_lake: np.ndarray,
) -> np.ndarray:
    """Keep a qualifying river connected to its mouth across positive-flow transit losses.

    Physical discharge can decrease downstream through evaporation, so applying the percentile
    independently at every node can punch gaps in an otherwise-real channel. Once a node is
    selected, include its positive-discharge downstream path up to the first ocean, lake, or dry
    node. The shared visited set keeps this linear in the routed graph even when many selected
    headwaters feed the same mainstem.
    """
    selected = np.array(is_river, dtype=bool, copy=True)
    visited: set[int] = set()
    for source in np.flatnonzero(selected):
        node = int(source)
        while True:
            target = int(flow_target[node])
            if (
                target < 0
                or target >= len(selected)
                or is_ocean[target]
                or is_lake[target]
                or discharge_m3_s[target] <= 0.0
                or target in visited
            ):
                break
            visited.add(target)
            selected[target] = True
            node = target
    return selected


def meltwater_m3_per_yr(melt_m: np.ndarray, areas_m2: np.ndarray, years: float) -> np.ndarray:
    """Glacier meltwater as a rate in m^3/yr: `ice_sheet.step_ice`'s per-node melt (metres of
    ice over the cell, summed over the step) times cell area, divided by the step length. Its
    total over all nodes is `IceLedger.meltwater_m3 / years`."""
    return melt_m * areas_m2 / max(years, 1.0)


# Metres of glacier ice per metre of liquid water freezing in place: river water settling at a
# frozen node (issue #355) and a standing lake freezing solid (issue #503). It's
# the rate snowfall builds ice at: GLACIER_ACCUMULATION_RATE turns mm/yr over Myr into metres,
# so a metre of water (1000 mm for 1e-6 Myr-years) is GLACIER_ACCUMULATION_RATE / 1000 m of
# ice, 1/50,000 of real. Frozen river ice therefore keeps the scale of the old
# `GLACIER_ACCUMULATION_RATE * water_deposited` sizing, whose stacked mm/yr is roughly
# catchment volume over cell area, and stays comparable with snowfall in the ledger. Lakes
# froze 1:1 until #503: their water balances in real units (#297), so a freezing lake fed the
# ice about 35x faster than all snowfall, and filling basins within the step grew the ice
# without bound.
# TODO(#353 decision B): when glacier mass balance moves to physical units, convert by mass
# instead of this model-unit factor. Water (1000 kg/m^3) freezing to glacier ice (about
# 917 kg/m^3, denser than the snow and firn it compacts from, and denser still deep in a
# sheet under pressure) gains roughly 9% in volume, and meltwater (`meltwater_m3_per_yr`,
# booked 1:1 today) should use the same density ratio.
WATER_TO_ICE = GLACIER_ACCUMULATION_RATE / 1000.0


def _compute_flow_direction(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    water_surface_offset: np.ndarray,
    filled_elevation: np.ndarray,
    spill_target: np.ndarray,
    prev_channel_depth: np.ndarray,
    is_frozen: np.ndarray,
    apply_freeze: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Flow target per node, among its k nearest neighbors strictly below its own elevation
    (a downhill candidate) -- preferring whichever candidate already has the deepest
    established channel (`prev_channel_depth > CHANNEL_PREFERENCE_THRESHOLD_M`), falling back
    to plain steepest descent only when *no* downhill candidate has a real channel yet. Real
    rivers meander within, and stay inside, their own existing valley rather than
    recalculating the mathematically steepest path every step -- this is that preference
    applied to routing itself, distinct from (and upstream of) erosion.py's own
    channel_boost, which only affects how *fast* a node erodes once water is already flowing
    there. -1 if there's no downhill candidate at all (a sink) -- unless the sink's current
    water surface (elevation + `water_surface_offset`, the previous step's water surface above
    the cell's centre, a one-step-lagged "memory" of the water surface; compute_hydrology
    recovers it from lake_depth through lakes.CellHypsometry) has already reached its
    basin's true spill point, in which case it redirects to spill_target instead of staying
    a dead-end sink forever. Ocean nodes are never routed (they're a destination, not a
    source).

    Also returns `should_spill`, the exact per-node mask of which sinks got redirected this
    way -- `compute_hydrology` needs this specific mask (not something recomputed after the
    fact from `flow_target == spill_target`, which could false-positive on an ordinary
    steepest-descent node that happens to coincide with a neighboring basin's spill_target) to
    know which lakes just found a real outlet, for the breach-erosion term (see its own
    docstring in `compute_hydrology`).

    A frozen node (`is_frozen`, gated at the real freezing point `FREEZE_POINT_C` on the
    annual mean -- unrelated to the seasonal glacier-melt cycle) can't pass
    liquid water downstream at all this step -- rivers freeze over. This override is applied
    *after* the should_spill redirect above, unconditionally, rather than folded into the same
    is_sink/should_spill check that redirect uses: an ordinary, otherwise-unobstructed
    downhill node has `filled_elevation` exactly equal to its own elevation (no real rim ever
    blocks it), so forcing such a node to look like a sink would make it trivially satisfy
    `water_surface >= filled_elevation` and get redirected right back onto its own normal
    downhill neighbor via spill_target -- silently undoing the freeze for most ordinary
    terrain. Applying the override to the final flow_target instead has no such escape hatch;
    `should_spill`'s own returned mask still reflects genuine basin-spill events only, since
    it's computed before this override runs.

    `apply_freeze=False` skips that last override, returning the *un*-frozen-gated target
    instead -- for glacial ice, which flows downhill under its own weight via internal
    deformation regardless of whether this particular step happens to be below freezing,
    unlike liquid water/rivers (which do stop dead at a frozen node). `compute_hydrology` calls
    this twice, once each way, rather than folding an ice-specific branch into the water path:
    reusing the exact same steepest-descent/channel-preference/spill-redirect logic for both
    keeps the two flow directions from silently drifting apart, and callers that only need
    water's own semantics (every existing caller, including this function's own test suite)
    see no change at all, since `apply_freeze` defaults to the prior (always-frozen-gated)
    behavior."""
    n = len(elevation)
    neighbor_elev = elevation[neighbor_idx]  # (n, k)
    own_elev = elevation[:, None]
    lower_mask = neighbor_elev < own_elev

    steepest_masked = np.where(lower_mask, neighbor_elev, np.inf)
    steepest_col = np.argmin(steepest_masked, axis=1)
    rows = np.arange(n)
    has_lower = np.isfinite(steepest_masked[rows, steepest_col])

    neighbor_channel = prev_channel_depth[neighbor_idx]  # (n, k)
    channelized_mask = lower_mask & (neighbor_channel > CHANNEL_PREFERENCE_THRESHOLD_M)
    has_channelized = np.any(channelized_mask, axis=1)
    channel_masked = np.where(channelized_mask, neighbor_channel, -np.inf)
    channel_col = np.argmax(channel_masked, axis=1)

    best_col = np.where(has_channelized, channel_col, steepest_col)
    flow_target = np.where(has_lower, neighbor_idx[rows, best_col], -1).astype(np.int64)

    is_sink = (flow_target < 0) & ~is_ocean
    water_surface = elevation + water_surface_offset
    # A rim at or below the sink's own centre (a breach notch cut to the floor, issue #297)
    # spills whatever water the sink holds. The water surface of a partly wet cell sits below
    # its centre (CellHypsometry), so testing it against a centre-height rim held a breached
    # pit's first shallow water forever -- and with it rivers of up to 13,000 m^3/s (#503).
    should_spill = is_sink & ((water_surface >= filled_elevation) | (filled_elevation <= elevation))
    flow_target = np.where(should_spill, spill_target, flow_target).astype(np.int64)
    if apply_freeze:
        flow_target = np.where(is_frozen & ~is_ocean, -1, flow_target).astype(np.int64)
    flow_target = np.where(is_ocean, -1, flow_target).astype(np.int64)
    return flow_target, should_spill


def _slope_to_flow_target(points: np.ndarray, elevation: np.ndarray, flow_target: np.ndarray) -> np.ndarray:
    """Dimensionless rise/run from each node to its own flow_target (0 for a node with no
    target) -- the slope a river's mouth speed uses (see group_rivers), distinct from
    erosion.py's own SLOPE_NEIGHBOR_COUNT=4 "lowest of k neighbors" slope: it follows
    `flow_target` (this module's own k=FLOW_NEIGHBOR_COUNT graph), not erosion's separate,
    smaller neighborhood."""
    has_target = flow_target >= 0
    safe_target = np.where(has_target, flow_target, 0)
    drop = np.where(has_target, elevation - elevation[safe_target], 0.0)
    run_m = geometry.angular_distance(points, points[safe_target]) * PLANET_RADIUS_KM * 1000.0
    run_m = np.maximum(run_m, 1.0)
    return np.where(has_target, np.clip(drop, 0.0, None) / run_m, 0.0)


def compute_river_speed(slope: np.ndarray, discharge_m3_s: np.ndarray | float) -> np.ndarray:
    """Stylized, unitless speed using physical discharge on the default-density reference scale."""
    flow_scale = discharge_m3_s_to_flow_accum_mm_yr(np.asarray(discharge_m3_s, dtype=float))
    return RIVER_SPEED_COEFFICIENT * np.sqrt(np.clip(slope, 0.0, None)) * np.power(flow_scale, RIVER_SPEED_DISCHARGE_EXPONENT)


def _routing_order(elevation: np.ndarray, is_ocean: np.ndarray, flow_target: np.ndarray) -> tuple[list[int], list[int]]:
    """Land nodes in an order where every node comes before its own flow_target, plus any
    nodes that can't be ordered at all (on, or downstream of, a flow cycle).

    Elevation-descending is such an order whenever every land->land edge runs strictly
    downhill, which is the common case. But spill routing (a filled lake's outlet) and ice
    routing (`ice_flow_target`, judged against elevation + lake + ice) can point a node at a
    higher neighbour, and an elevation sweep then adds that node's flux to a target it already
    visited -- flux that is never passed on or deposited (issue #275: ~0.1% of routed
    sediment and ~9% of glacier-carried sediment silently vanished per step). In that case
    this falls back to Kahn's algorithm, highest available node first, so the order still
    matches the elevation sweep wherever the graph allows."""
    land_indices = np.nonzero(~is_ocean)[0]
    targets = flow_target[land_indices]
    land_edge = targets >= 0
    land_edge[land_edge] = ~is_ocean[targets[land_edge]]
    if not np.any(elevation[targets[land_edge]] >= elevation[land_indices[land_edge]]):
        return land_indices[np.argsort(-elevation[land_indices])].tolist(), []

    n = len(elevation)
    indegree = np.bincount(targets[land_edge], minlength=n)
    heap = [(-float(elevation[i]), int(i)) for i in land_indices[indegree[land_indices] == 0]]
    heapq.heapify(heap)
    flow_target_list = flow_target.tolist()
    is_ocean_list = is_ocean.tolist()
    indegree_list = indegree.tolist()
    elevation_list = elevation.tolist()
    order = []
    while heap:
        _, i = heapq.heappop(heap)
        order.append(i)
        target = flow_target_list[i]
        if target >= 0 and not is_ocean_list[target]:
            indegree_list[target] -= 1
            if indegree_list[target] == 0:
                heapq.heappush(heap, (-elevation_list[target], target))
    if len(order) == len(land_indices):
        return order, []
    placed = np.zeros(n, dtype=bool)
    placed[order] = True
    return order, land_indices[~placed[land_indices]].tolist()


def flow_cycles(flow_target: np.ndarray) -> list[np.ndarray]:
    """Every cycle in the `flow_target` graph, as an array of its nodes. Pointer doubling
    finds one node on each cycle (whatever doesn't end at -1 after `n` hops sits on one); a
    short walk from it lists the rest."""
    n = len(flow_target)
    if n == 0:
        return []
    jump = np.append(np.where(flow_target >= 0, flow_target, n), n)
    hops = 1
    while hops < n:
        jump = jump[jump]
        hops *= 2
    ends = jump[:n]
    cycles: list[np.ndarray] = []
    seen: set[int] = set()
    for start in np.unique(ends[ends < n]).tolist():
        if start in seen:
            continue
        cycle = [start]
        node = int(flow_target[start])
        while node != start:
            cycle.append(node)
            node = int(flow_target[node])
        seen.update(cycle)
        cycles.append(np.asarray(cycle, dtype=np.int64))
    return cycles


# Rounds of the within-step spill/merge fixed point (_spill_within_step). Each round either
# opens a basin's outlet or unites basins one level higher, so real worlds settle in a few;
# this only bounds a pathological case.
SPILL_WITHIN_STEP_MAX_ROUNDS = 32


@dataclass
class SpillBasins:
    """What _spill_within_step needs to know about this step's depression hierarchy."""

    forest: list["lakes.Lake"]
    governing: list["lakes.Lake"]  # lakes.governing_bodies: the rim each leaf's sink obeys
    hypsometry: "lakes.CellHypsometry"
    neighbor_idx: np.ndarray  # the routing graph the hierarchy's catchments were built on
    stored_m3: np.ndarray  # water already standing in each cell at the start of the step
    evaporation_m_per_yr: np.ndarray  # open-water evaporation per cell
    frozen: np.ndarray  # sinks whose water freezes where it settles
    keep: np.ndarray | None = None  # share of passing water each node keeps (1 - transit loss); None keeps all


def _spill_within_step(
    basins: SpillBasins,
    flow_target: np.ndarray,
    should_spill: np.ndarray,
    filled_elevation: np.ndarray,
    spill_target: np.ndarray,
    route,
    years: float,
):
    """Opens, within this step, the outlet of every basin this step's water fills, and unites
    basins whose overflow runs into each other (issue #503).

    `_compute_flow_direction` spills a sink only once *last* step's water reached its rim, and
    `lakes.compute_spill_routing` unites two basins only once last step's water pooled them. At
    100 kyr steps that lag stranded rivers: a river crossing a chain of dry pits advanced one
    pit per step, and two basins that filled to their shared saddle in the same step each spilled
    into the other, so the water circled between them and reached neither the sea nor a lake.

    Each round sweeps the governing basins (`lakes.governing_bodies`) from the highest floor
    down. A basin whose liquid inflow this step, net of the evaporation of a lake at its rim
    (`lakes.lake_evaporation_m3_per_yr`), fills its remaining room opens its outlet, and its
    inflow is passed on to the basin that outlet drains into, so a chain of pits opens in one
    sweep. Then the basins on each flow cycle are united at their lowest common ancestor, which
    is re-tested the same way: it spills at its own rim, or its water pools at its lowest sink.
    A cycle that instead runs through a channel into a neighbouring catchment (routing prefers
    an established channel, while the hierarchy's catchments follow steepest descent) is broken
    by sending that node down its steepest descent, the path the hierarchy drains it by.
    `route(flow_target)` re-routes after any change and returns `(through, settled)`. Repeats to
    a fixed point, at most SPILL_WITHIN_STEP_MAX_ROUNDS rounds; reaching the cap logs a warning.

    Returns `(flow_target, should_spill, filled_elevation, spill_target, through, settled)`."""
    resolver = _SpillResolver(basins, flow_target, should_spill, filled_elevation, spill_target, route, years)
    for _ in range(SPILL_WITHIN_STEP_MAX_ROUNDS):
        opened = resolver.open_filled_basins()
        united = resolver.untangle_cycles()
        if not (opened or united):
            break
    else:
        _logger.warning(
            "basin spills did not settle within %d rounds; some rivers may stop short of an outlet this step",
            SPILL_WITHIN_STEP_MAX_ROUNDS,
        )
    return (
        resolver.flow_target, resolver.should_spill, resolver.filled_elevation, resolver.spill_target,
        resolver.through, resolver.settled,
    )


class _SpillResolver:
    """The state `_spill_within_step` iterates on: the routing, and which basin governs each sink."""

    def __init__(self, basins: SpillBasins, flow_target, should_spill, filled_elevation, spill_target, route, years):
        self.basins = basins
        self.flow_target = flow_target.copy()
        self.should_spill = should_spill.copy()
        self.filled_elevation = filled_elevation.copy()
        self.spill_target = spill_target.copy()
        self.route = route
        self.years = years
        n = len(flow_target)
        self.keep = np.ones(n) if basins.keep is None else basins.keep
        self.climate = lakes.LakeClimate(np.zeros(n), basins.evaporation_m_per_yr)
        self.through, self.settled = route(self.flow_target)

        self.parent: dict[int, "lakes.Lake | None"] = {}
        self.depth: dict[int, int] = {}
        self.leaf_of_node = np.full(n, -1, dtype=np.int64)
        stack: list[tuple["lakes.Lake", "lakes.Lake | None", int]] = [(root, None, 0) for root in basins.forest]
        while stack:
            lake, up, level = stack.pop()
            self.parent[id(lake)] = up
            self.depth[id(lake)] = level
            if lake.children:
                stack.extend((child, lake, level + 1) for child in lake.children)
            else:
                self.leaf_of_node[lake.members] = lake.sink_node_idx
        self.body_of_sink: dict[int, "lakes.Lake"] = {}
        for body in basins.governing:
            for sink in lakes.leaf_sinks(body).tolist():
                self.body_of_sink[sink] = body

        elevation = basins.hypsometry.elevation
        lowest = basins.neighbor_idx[np.arange(n), np.argmin(elevation[basins.neighbor_idx], axis=1)]
        self.steepest = np.where(elevation[lowest] < elevation, lowest, -1)

    def _reroute(self) -> None:
        self.through, self.settled = self.route(self.flow_target)

    def _liquid(self, sinks: np.ndarray) -> np.ndarray:
        """The sinks whose water stays liquid; a frozen sink's water becomes river ice."""
        return sinks[~self.basins.frozen[sinks]]

    def _fills(self, body: "lakes.Lake", inflow_m3_per_yr: float) -> bool:
        if body.max_depth is None or body.outlet_target_idx < 0:
            return False
        hyp = self.basins.hypsometry
        members = body.members
        room = hyp.volume(body.max_depth, members) - float(self.basins.stored_m3[members].sum())
        evaporation = lakes.lake_evaporation_m3_per_yr(body.max_depth, members, hyp, self.climate)
        return (inflow_m3_per_yr - evaporation) * self.years >= room

    def _open_outlet(self, body: "lakes.Lake", sinks: np.ndarray) -> None:
        live = self._liquid(sinks)
        self.flow_target[live] = body.outlet_target_idx
        self.should_spill[live] = True

    def _pass_on(self, start: int, amount: float, pending: np.ndarray) -> None:
        """Adds `amount` (m^3/yr) to the sink where water entering the routing at `start` settles,
        less the transit losses on the way."""
        node, hops = start, 0
        while hops <= len(self.flow_target):
            amount *= self.keep[node]
            target = int(self.flow_target[node])
            if target < 0:
                pending[node] += amount
                return
            node, hops = target, hops + 1

    def open_filled_basins(self) -> bool:
        """Opens the outlet of every governing basin this step's water fills, sweeping from the
        highest floor down so water released upstream counts toward the basins below it."""
        bodies = sorted({id(b): b for b in self.body_of_sink.values()}.values(), key=lambda b: -b.floor_elevation)
        pending = np.zeros(len(self.flow_target))
        opened = False
        for body in bodies:
            live = self._liquid(lakes.leaf_sinks(body))
            if not np.any(self.flow_target[live] < 0):
                continue
            inflow = float(self.settled[live].sum() + pending[live].sum())
            if self._fills(body, inflow):
                self._open_outlet(body, live)
                self._pass_on(body.outlet_target_idx, inflow, pending)
                opened = True
        if opened:
            self._reroute()
        return opened

    def untangle_cycles(self) -> bool:
        """Breaks every flow cycle: a channel crossing into another catchment follows its
        steepest descent again, and basins spilling into each other are united."""
        changed = False
        for cycle in flow_cycles(self.flow_target):
            changed = self._follow_catchment(cycle) or self._unite(cycle) or changed
        if changed:
            self._reroute()
        return changed

    def _follow_catchment(self, cycle: np.ndarray) -> bool:
        # A channel can lead a node downhill into a neighbouring catchment, while the hierarchy
        # drains it down its steepest descent. On a cycle that disagreement closes the loop, so
        # the node follows the hierarchy's own path again.
        elevation = self.basins.hypsometry.elevation
        target = self.flow_target[cycle]
        downhill = cycle[(elevation[target] < elevation[cycle]) & ~self.should_spill[cycle]]
        crossing = downhill[self.leaf_of_node[self.flow_target[downhill]] != self.leaf_of_node[downhill]]
        self.flow_target[crossing] = self.steepest[crossing]
        return len(crossing) > 0

    def _unite(self, cycle: np.ndarray) -> bool:
        sinks_on_cycle = {int(s) for s in self.leaf_of_node[cycle] if s >= 0}
        ancestor = self._lowest_common_ancestor([self.body_of_sink[s] for s in sinks_on_cycle if s in self.body_of_sink])
        if ancestor is None:
            return False
        sinks = lakes.leaf_sinks(ancestor)
        for sink in sinks.tolist():
            self.body_of_sink[sink] = ancestor
        self.filled_elevation[sinks] = ancestor.max_depth if ancestor.max_depth is not None else np.inf
        self.spill_target[sinks] = ancestor.outlet_target_idx
        live = self._liquid(sinks)
        if self._fills(ancestor, float(self.settled[live].sum())):
            self._open_outlet(ancestor, sinks)
        elif len(live):
            # One lake: its water pools at its lowest ground, not at every sink it covers.
            lowest = int(live[np.argmin(self.basins.hypsometry.elevation[live])])
            self.flow_target[live] = lowest
            self.should_spill[live] = True
            self.flow_target[lowest] = -1
            self.should_spill[lowest] = False
        return True

    def _lowest_common_ancestor(self, bodies: list["lakes.Lake"]) -> "lakes.Lake | None":
        nodes = list({id(b): b for b in bodies}.values())
        if not nodes:
            return None
        while len({id(b) for b in nodes}) > 1:
            deepest = max(self.depth[id(b)] for b in nodes)
            lifted = [self.parent[id(b)] if self.depth[id(b)] == deepest else b for b in nodes]
            if any(b is None for b in lifted):
                return None  # different roots: no shared basin to unite into
            nodes = lifted
        return nodes[0]


def route_downstream(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    flow_target: np.ndarray,
    source_amount: np.ndarray,
    retain_fraction: np.ndarray | None = None,
    loss_fraction: np.ndarray | None = None,
    erosion_only_source_amount: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Single forward sweep over land nodes in elevation-descending order, accumulating
    `source_amount` downstream along `flow_target` edges. Correct in one pass because every
    node's target is visited *later* in this same order -- see `_routing_order` for the
    spill/ice edges that point uphill and the topological order used when they exist. `loss_fraction` is an in-transit evaporative loss (see
    hydrology.RIVER_EVAPORATION_* constants): unlike `retain_fraction`, which *deposits* the
    retained share locally (sediment settling out of the flow), a lost share simply vanishes
    -- it isn't added to `deposited` at all, since it left as atmospheric moisture, not as
    material staying in the channel. Applied first, at every node's own turn, before
    retention -- so, like flow_accum itself, the lost fraction compounds along a long flow
    path rather than being taken once at the very end. Returns (through_flux, deposited): each
    node's own accumulated flux passing through it (net of any loss), and how much settled
    there (from `retain_fraction` and/or reaching a sink/ocean). When
    `erosion_only_source_amount` is supplied, it is propagated beside the primary source with
    the same losses but without retention/deposition, and returned as a third through-flux array."""
    n = len(elevation)
    through_flux = np.where(~is_ocean, source_amount, 0.0).astype(np.float64).tolist()
    erosion_only_flux = (
        np.where(~is_ocean, erosion_only_source_amount, 0.0).astype(np.float64).tolist()
        if erosion_only_source_amount is not None
        else None
    )
    deposited = [0.0] * n
    retain = (retain_fraction if retain_fraction is not None else np.zeros(n)).tolist()
    loss = (loss_fraction if loss_fraction is not None else np.zeros(n)).tolist()
    is_ocean_list = is_ocean.tolist()
    flow_target_list = flow_target.tolist()

    order, stuck = _routing_order(elevation, is_ocean, flow_target)
    for i in order:
        evaporated_here = through_flux[i] * loss[i]
        through_flux[i] -= evaporated_here
        if erosion_only_flux is not None:
            erosion_only_flux[i] -= erosion_only_flux[i] * loss[i]

        retained_here = through_flux[i] * retain[i]
        through_flux[i] -= retained_here
        deposited[i] += retained_here

        target = flow_target_list[i]
        if target < 0:
            deposited[i] += through_flux[i]
            continue
        if is_ocean_list[target]:
            deposited[target] += through_flux[i]
        else:
            through_flux[target] += through_flux[i]
            if erosion_only_flux is not None:
                erosion_only_flux[target] += erosion_only_flux[i]
    # A node on a flow cycle has no valid place in the sweep: it still loses/retains its own
    # share like any other node, then keeps the remainder, like a sink, rather than losing it.
    # Same loss -> retain -> dead-end sequence the sweep applies to a sink (target < 0).
    for i in stuck:
        through_flux[i] -= through_flux[i] * loss[i]
        if erosion_only_flux is not None:
            erosion_only_flux[i] -= erosion_only_flux[i] * loss[i]
        retained_here = through_flux[i] * retain[i]
        through_flux[i] -= retained_here
        deposited[i] += retained_here + through_flux[i]

    if erosion_only_flux is not None:
        return np.array(through_flux), np.array(deposited), np.array(erosion_only_flux)
    return np.array(through_flux), np.array(deposited)


def seasonal_positive_degrees(mean_c: np.ndarray, amplitude_c: np.ndarray, threshold_c: float = GLACIER_MELT_THRESHOLD_C) -> np.ndarray:
    """Year-averaged degrees above `threshold_c` of a sinusoidal annual cycle with annual mean
    `mean_c` and half-amplitude `amplitude_c`: mean of max(mean + amplitude * cos(t) -
    threshold, 0) over t. Closed form, with d = mean - threshold and phi = arccos(-d /
    amplitude) the half-width of the above-threshold part of the cycle:
    (d * phi + sqrt(amplitude^2 - d^2)) / pi when |d| < amplitude, d when the whole year is
    above, 0 when it's all below. Zero amplitude reduces to max(d, 0). Smooth and increasing
    in both `mean_c` and `amplitude_c` (for d < amplitude)."""
    d = np.asarray(mean_c, dtype=float) - threshold_c
    amp = np.clip(np.asarray(amplitude_c, dtype=float), 0.0, None)
    d, amp = np.broadcast_arrays(d, amp)
    out = np.clip(d, 0.0, None)
    partial = np.abs(d) < amp
    if np.any(partial):
        dp, ap = d[partial], amp[partial]
        out[partial] = (dp * np.arccos(-dp / ap) + np.sqrt(ap * ap - dp * dp)) / np.pi
    return out


def _ice_surface_flow_target(
    ice_surface: np.ndarray, neighbor_idx: np.ndarray, glacier_depth: np.ndarray, ground_target: np.ndarray
) -> np.ndarray:
    """Where material riding the ice moves next (issue #353 section 6): on visible ice, the
    steepest descent of the ice surface (elevation + ice); elsewhere, or on an ice-surface
    pit, `ground_target` (the ground's own downhill path)."""
    on_ice = glacier_depth >= GLACIER_VISIBLE_DEPTH_M
    neighbor_surface = ice_surface[neighbor_idx]
    lower = np.where(neighbor_surface < ice_surface[:, None], neighbor_surface, np.inf)
    col = np.argmin(lower, axis=1)
    rows = np.arange(len(ice_surface))
    has_lower = np.isfinite(lower[rows, col])
    return np.where(on_ice & has_lower, neighbor_idx[rows, col], ground_target).astype(np.int64)


def compute_hydrology(
    world: "World",
    precipitation_at_nodes: np.ndarray,
    temperature_at_nodes: np.ndarray,
    years: float,
    node_cloud: tuple[np.ndarray, list[Plate]] | None = None,
    seasonal_amplitude_at_nodes: np.ndarray | None = None,
    history_years: float | None = None,
) -> HydrologyFields:
    """Runs the full flow-routing pipeline against the world's current node cloud and this
    step's climate: (lake freeze, if cold enough) -> depression hierarchy + spill routing (see
    `lakes.build_lake_hierarchy`/`lakes.compute_spill_routing`) -> flow direction ->
    glacier mass balance and shape (`ice_sheet.step_ice`) -> precipitation(-minus-frozen)-plus-meltwater-plus-lake-
    breach downstream accumulation -> lake growth/evaporation/merge/split, resolved against the
    *same* depression hierarchy spill routing already used (`lakes.resolve_lakes`) -> river
    classification (excluding flooded lake nodes, so a river ends at a lake's shore
    rather than being treated as continuing through it). A lake that's actively spilling this
    step also adds `LAKE_BREACH_EROSION_COEFFICIENT`-scaled extra "water" into the legacy route
    at its own sink, sized by its connected lake's own node count (see that constant's own
    comment). River erosion receives this forcing from a separate sidecar accumulation;
    physical discharge, channel growth and river classification exclude it. `precipitation_at_nodes`/
    `temperature_at_nodes` are precomputed by the caller (erosion.py, which already needs the
    same climate-grid lookups for its own erosion terms) rather than looked up again here.
    `node_cloud`, likewise, is erosion.py's own already-gathered (points, plates_in_order) pair (see
    plates.gather_node_positions), reused here instead of re-deriving every node's world
    position from scratch a second time this same step.

    A node colder than FREEZE_POINT_C never holds liquid water at all: any precipitation
    there is treated as fully frozen (no partial liquid/frozen split -- matching erosion.py's
    existing "use precipitation is enough" simplification, just gated by the freeze
    threshold), an existing lake sitting there freezes solid into glacier_depth this same
    step, and a river reaching it stops flowing entirely this step (its incoming flux freezes
    in place too, also joining glacier_depth) -- both *before* flow_target is (re)computed, so
    a lake that just froze is correctly treated as a genuine sink again this step (see
    `lakes.step_lakes`'s own docstring for the freeze handling this feeds into), and a frozen
    river node is correctly forced to a sink in `_compute_flow_direction` (see its own
    docstring for why that override can't reuse the ordinary should_spill escape valve).
    Freezing solid for a step doesn't imply building a *permanent* glacier: glacier surface
    melt runs on the node's imputed seasonal cycle (`seasonal_amplitude_at_nodes`, see
    GLACIER_MELT_THRESHOLD_C's comment), so ice formed where only the annual mean is below
    freezing still melts back in summer. `None` means no seasonal swing anywhere.

    `history_years` (default `years`) is the span basins have had to fill and outlets to be
    cut: it sets the breach carving budget and which basins this step's water fills
    (`_spill_within_step`). Everything else still reads `years`. Generation's zero-duration
    bootstrap passes a longer span, since its terrain stands for a drainage history it never
    ran: with no time to fill or carve, every pit ended a river (#218).

    **`should_spill` (below) and a lake's own `is_spilling` (`lakes._resolve`, set inside
    `lakes.resolve_lakes` near the end of this function) are now the same test -- "is this
    lake's water at or above its own rim" -- applied to the same depression-hierarchy geometry,
    just necessarily evaluated at two different points in the step** -- not two
    independently-computed conditions that can drift apart on their own terms (that used to be
    true, back when `should_spill` read a separately-computed Dijkstra basin-spill result --
    see this module's own docstring). `should_spill` has to run early, against *last* step's
    `lake_depth`, because `route_downstream` (which produces this step's `water_deposited`,
    the very thing a lake needs to know its own *this-step* inflow) needs `flow_target` first,
    which needs `should_spill` first. `resolve_lakes`'s own `is_spilling` runs afterward,
    against the freshly-resolved *this-step* level. So a lake right at its rim can still show
    `should_spill=True` while this same step's `is_spilling` says `False` (it hadn't quite
    risen that last bit as of last step's snapshot) or vice versa -- a real, unavoidable
    one-step lag, the same character every other persisted quantity in this codebase already
    has, not a bug to chase further."""
    points, elevation, prev_lake_depth, prev_glacier_depth, prev_channel_depth, prev_channel_width, prev_silt_depth, is_ocean, plates_in_order = _gather_nodes(world, node_cloud=node_cloud)
    n = len(points)
    if n <= FLOW_NEIGHBOR_COUNT:
        empty_i = np.zeros(n, dtype=np.int64)
        empty_f = np.zeros(n)
        return HydrologyFields(
            points, elevation, is_ocean, np.zeros((n, 0), dtype=np.int64), empty_i, empty_f, empty_f, empty_f, empty_i,
            np.zeros(n, dtype=bool), empty_f, empty_f, plates_in_order,
            silt_depth=prev_silt_depth, silt_deposited=np.zeros(n), ice_flow_target=empty_i,
            channel_width=prev_channel_width,
        )

    neighbor_idx = _build_neighbor_graph(points, world=world)
    # `is_ocean` from _gather_nodes is the bare `elevation <= sea_level` test; refine it to
    # "below sea level *and* connected to the world ocean" now that the k-NN graph exists, so an
    # enclosed interior depression that dipped below sea level is routed/pooled/silted as an
    # endorheic basin (see connected_ocean_mask) rather than treated as a vanishing ocean sink.
    is_ocean = connected_ocean_mask(points, elevation, world.sea_level_m, neighbor_idx)
    if world.premade_world_id == "earth":
        is_ocean = force_earth_strait_connections(points, is_ocean, elevation, world.sea_level_m, neighbor_idx)

    # Precipitation falls as snow on a node whose sampled surface temperature is below the
    # real freezing point, OR that already carries a real ice cap (GLACIER_VISIBLE_DEPTH_M
    # of accumulated ice, the same "genuine ice, not a trace" threshold rendering uses): a
    # marginal thaw doesn't turn the snow on kilometres of ice to rain. The glacier's own
    # surface melt (`ice_sheet.step_ice`, temperature-driven) still ablates the cap regardless.
    is_frozen = (temperature_at_nodes < FREEZE_POINT_C) | (prev_glacier_depth >= GLACIER_VISIBLE_DEPTH_M)
    # Liquid water -- a river crossing a node, the water settling at a sink, a standing lake --
    # freezes only where the annual mean is below freezing. The ice clause above doesn't apply
    # to it (issue #503): a river reaching a warm ice node used to stop there and settle as
    # river ice, which kept the node frozen and the river blocked, so a 36,000 m^3/s river
    # rebuilt a single-node, 600 m ice plug every step at +0.4C. The ice's own meltwater
    # refroze in place the same way. Meltwater and rivers now drain past warm ice along the
    # ground. Lakes followed the same rule first (issue #353): such a lake (sub- or proglacial)
    # is refilled every step by runoff and meltwater balanced in real units (#297). A lake that
    # does freeze turns into ice at WATER_TO_ICE, the rate river water does. Computed *before* the
    # depression hierarchy/spill routing below since compute_spill_routing needs
    # lake_depth_adjusted, not just bare terrain -- see that function's own docstring for why.
    water_frozen = temperature_at_nodes < FREEZE_POINT_C
    lake_frozen = water_frozen
    freezing_lake = lake_frozen & (prev_lake_depth > 0.0)
    lake_depth_adjusted = np.where(freezing_lake, 0.0, prev_lake_depth)
    frozen_from_lake = np.where(freezing_lake, WATER_TO_ICE * prev_lake_depth, 0.0)

    # The depression hierarchy is built once, here, from bare terrain alone (no water/inflow
    # input needed -- see build_lake_hierarchy's own signature) so it can drive *both* this
    # step's spill routing (immediately below -- see compute_spill_routing's own docstring for
    # why lakes.py's own hierarchy, not a separate independently-computed basin-spill pass, is
    # now the single source of truth for where a lake's overflow goes) and, later, this same
    # step's lake water-balance resolution (`lakes.resolve_lakes`, once this step's physical water
    # is routed) -- the exact same `forest` object threaded through both, rather than paying for a
    # second, redundant `build_lake_hierarchy` call and risking the two disagreeing.
    #
    # With breaching on, the hierarchy's boundary edges use notched passage elevations, not
    # bare centres. A pit cheap enough to carve through this step then merges toward the ocean
    # at its own floor, so it drains instead of ponding. A costlier pit keeps its centre-height
    # rim and is filled or held as an endorheic lake as before (see breaching.py).
    areas = collect_all_accounting_areas_m2(plates_in_order, line_spacing_rad(world.node_density))
    # One in-transit loss rule for all physical water (issue #355): the breach gate and the
    # routing below both use it. Frozen nodes are dead ends and the ocean isn't a river, so
    # neither loses water in transit. A spilling sink or its rim loses the same share as any
    # other node: the breach surge's own exemption belongs to the legacy routing only.
    transit_loss = river_transit_loss_fraction(temperature_at_nodes, water_frozen | is_ocean)
    fill_years = years if history_years is None else history_years
    breach = None
    interface_pass = None
    if BREACH_DEPRESSIONS_ENABLED:
        craton_strength = cratons.strength(np.concatenate([p.collect("craton_crust_m") for p in plates_in_order]))
        mobile_cover = np.concatenate([p.collect("mobile_cover_m") for p in plates_in_order])
        # The gate sees the same water budget the lakes are balanced with below, as far as it
        # can before routing: liquid precipitation, and the same per-node river losses.
        breach = breaching.breach_depressions(
            elevation, is_ocean, neighbor_idx, prev_channel_depth,
            breaching.bedrock_carve_rate_m_per_myr(craton_strength), fill_years,
            water=breaching.WaterBalance(
                np.where(is_frozen, 0.0, precipitation_at_nodes), temperature_at_nodes, areas,
                loss_fraction=transit_loss,
            ),
            mobile_cover_m=mobile_cover,
            prior_notch_m=collect_all_breach_notch_depth(plates_in_order),
        )
        interface_pass = breaching.interface_pass_elevation(breach.passage_m, neighbor_idx)
    forest = lakes.build_lake_hierarchy(elevation, is_ocean, neighbor_idx, interface_pass_elevation=interface_pass)
    # Partial-cell flooding (issue #297): lake_depth is each cell's mean water depth over its
    # sub-cell relief, so the water surface is recovered through the same hypsometry rather
    # than as elevation + depth.
    hypsometry = lakes.CellHypsometry.from_terrain(elevation, areas, neighbor_idx)
    filled_elevation, spill_target, rim_node_idx = lakes.compute_spill_routing(
        forest, elevation, lake_depth_adjusted, hypsometry
    )
    # Water surface above each wet cell's centre, for _compute_flow_direction's spill test
    # (`elevation + offset >= filled_elevation`). It's the mean depth itself without relief.
    lake_surface_offset = np.zeros(n)
    wet_lake = lake_depth_adjusted > 0.0
    wet_idx = np.flatnonzero(wet_lake)
    lake_surface_offset[wet_idx] = hypsometry.level_from_depth(lake_depth_adjusted[wet_idx], wet_idx) - elevation[wet_idx]

    flow_target, should_spill = _compute_flow_direction(
        elevation, is_ocean, neighbor_idx, lake_surface_offset, filled_elevation, spill_target, prev_channel_depth, water_frozen
    )
    # Off the ice, material riding the ice (glacial till, landslide debris) leaves along the
    # ground's own downhill path, never gated by freezing -- see _compute_flow_direction's
    # apply_freeze docstring. On the ice it follows the ice surface instead (set below, once
    # this step's ice is solved).
    ice_flow_target, _ = _compute_flow_direction(
        elevation, is_ocean, neighbor_idx, lake_surface_offset, filled_elevation, spill_target, prev_channel_depth, is_frozen, apply_freeze=False
    )

    liquid_precip = np.where(is_frozen, 0.0, precipitation_at_nodes)

    # In a deep enough glacial phase, let ice pile up over freezing polar ocean (a floating
    # ice cap) instead of calving at the coast -- see SEA_ICE_FORMATION_TEMP_C. `climate` is
    # imported locally: it imports this module at load time, so a top-level import here would
    # be circular, but by the time this runs climate is fully initialized.
    from . import climate

    # Freezing polar ocean grows an ice cap; an ocean node that *already* holds one keeps it
    # (a bit of thermal-inertia hysteresis -- sea ice doesn't vanish the instant the water
    # edges back above its freezing point) -- but only while the glacial phase itself lasts:
    # once intensity falls back below the gate the cap is shed regardless (it calves, see
    # ice_sheet.step_ice).
    allow_ocean_ice = (
        is_ocean
        & ((temperature_at_nodes < SEA_ICE_FORMATION_TEMP_C) | (prev_glacier_depth >= GLACIER_VISIBLE_DEPTH_M))
        & (climate.ice_age_glacial_intensity(world) >= ICE_AGE_SEA_ICE_MIN_INTENSITY)
    )

    # Issue #353: ice volume is conserved per ice body and shaped by a plastic ice-sheet
    # profile (ice_sheet.py), replacing one-hop bed-slope flow. Imported here for the same
    # reason as climate above: ice_sheet imports this module.
    from . import ice_sheet

    new_glacier_depth, melt, ice = ice_sheet.step_ice(
        points, neighbor_idx, areas, elevation, collect_all_ice_load_deflection(plates_in_order), world.sea_level_m,
        prev_glacier_depth, frozen_from_lake,
        ice_sheet.SurfaceClimate(
            precipitation_at_nodes, temperature_at_nodes,
            seasonal_amplitude_at_nodes if seasonal_amplitude_at_nodes is not None else np.zeros(n), is_frozen,
        ),
        np.concatenate([p.collect("mobile_cover_m") for p in plates_in_order]), is_ocean, allow_ocean_ice, years,
    )
    ice_flow_target = _ice_surface_flow_target(elevation + new_glacier_depth, neighbor_idx, new_glacier_depth, ice_flow_target)

    # A spilling lake's own surface area feeds extra erosive "water" in at its sink, on top of
    # whatever ordinary precip/melt landed there this step -- see LAKE_BREACH_EROSION_
    # COEFFICIENT's own comment for why. Sized from *last* step's lake extent (lake_depth_
    # adjusted), the same one-step-lagged memory should_spill itself is already computed from.
    lake_component_size = _lake_component_sizes(lake_depth_adjusted > LAKE_MIN_VISIBLE_DEPTH_M, neighbor_idx)
    lake_breach_source = np.where(should_spill, LAKE_BREACH_EROSION_COEFFICIENT * lake_component_size, 0.0)

    # ...and, just as importantly, at the governing rim's own lake-side member too
    # (`rim_node_idx`, from `compute_spill_routing`) -- not only at the sink it's routed from.
    # See LAKE_BREACH_EROSION_COEFFICIENT's own comment for why a sink-only surge leaves the
    # true dam untouched whenever the lake's own rim happens to be the higher of the boundary
    # edge's two members. `np.add.at` (not a plain indexed assignment) since a merged, multi-
    # leaf lake's several original sinks can all share one governing rim node, and their
    # surges should sum there rather than overwrite each other. Excludes the (real, if
    # degenerate) case of a single-node lake whose only member is both its own sink *and* its
    # own rim -- there the sink and rim terms would land on the exact same node, doubling the
    # surge for no physical reason rather than covering two genuinely different points.
    rim_breach_source = np.zeros(n)
    spilling_idx = np.nonzero(should_spill)[0]
    has_rim = (rim_node_idx[spilling_idx] >= 0) & (rim_node_idx[spilling_idx] != spilling_idx)
    np.add.at(rim_breach_source, rim_node_idx[spilling_idx[has_rim]], lake_breach_source[spilling_idx[has_rim]])
    is_rim_breach = rim_breach_source > 0.0

    # River evaporation (see RIVER_EVAPORATION_* above) -- excluded at an already-frozen node
    # (its flow_target is already forced closed, see _compute_flow_direction) and, deliberately,
    # at a should_spill node or its own rim: the breach term above models a lake's own
    # concentrated overflow surge, sized purely from its surface area, not an ordinary trickle
    # that should also lose a further fraction to evaporation in the very same step it's cut
    # loose.
    river_evap_fraction = river_evaporation_fraction(
        temperature_at_nodes, years, water_frozen | is_ocean | should_spill | is_rim_breach
    )

    # No live consumer reads this legacy routing any more (issue #355 "Water units"): erosion,
    # geology and river visibility use the physical routing, and lake silt reads the physical
    # lake inflow (#363). It's kept as the fallback for caches without `discharge_m3_s`.
    breach_surge = lake_breach_source + rim_breach_source
    water_source = liquid_precip + melt + breach_surge
    flow_accum, water_deposited, breach_surge_accum = route_downstream(
        elevation,
        is_ocean,
        flow_target,
        water_source,
        loss_fraction=river_evap_fraction,
        erosion_only_source_amount=breach_surge,
    )

    # Physical water (issue #355 "Water units"). Lakes balance in real units (issue #297):
    # area-weighted runoff (the same Budyko curve and open-water evaporation the breach gate
    # uses) plus glacier melt, routed as volume along the same flow paths, with the same
    # transit loss as the gate. Its through-flow is the river discharge; what settles feeds the
    # lakes and frozen river ice below. Basins this water fills open their outlets within the
    # step, which re-routes it (_spill_within_step, issue #503). The legacy routing above keeps
    # the lagged spills: only the breach surge still reads it.
    evaporation_mm = breaching.open_water_evaporation_mm(temperature_at_nodes)
    runoff_m3 = breaching.runoff_mm(liquid_precip, evaporation_mm) / 1000.0 * areas
    physical_source_m3 = runoff_m3 + meltwater_m3_per_yr(melt, areas, years)
    flow_target, should_spill, filled_elevation, spill_target, through_m3, settled_m3 = _spill_within_step(
        SpillBasins(
            forest, lakes.governing_bodies(forest, elevation, lake_depth_adjusted, hypsometry), hypsometry, neighbor_idx,
            lake_depth_adjusted * areas, evaporation_mm / 1000.0, water_frozen, keep=1.0 - transit_loss,
        ),
        flow_target, should_spill, filled_elevation, spill_target,
        lambda targets: route_downstream(elevation, is_ocean, targets, physical_source_m3, loss_fraction=transit_loss),
        fill_years,
    )

    # A river blocked by its own freeze doesn't just vanish -- its water piles up as ice right
    # where it froze (route_downstream already settled it there, since a frozen node's
    # flow_target was forced to -1 above), joining this step's glacier_depth. Added after the
    # ice solve (which needs this routing pass's meltwater, so it can't know what settles)
    # rather than folded into its accumulation -- an explicit one-step lag (issue #353 section
    # 5): this ice joins its body's opening volume next step, and the ledger books it as
    # `frozen_river_ice_m3`.
    # It's the physical water settled over the step, at the glacier model's exchange rate
    # (WATER_TO_ICE, issue #355), so the breach surge never becomes ice. It's then
    # re-clipped to GLACIER_MAX_DEPTH_M so a large river freezing into a cold sink can't
    # create a physically absurd ice column in a single step. What the clip drops is booked
    # (`frozen_river_capped_m3`) and logged as a warning (review of #359): reaching it means a
    # missing sink.
    settled_water_m = np.divide(settled_m3 * years, areas, out=np.zeros(n), where=areas > 0.0)
    frozen_river_ice = WATER_TO_ICE * settled_water_m
    with_river_ice = np.where(
        water_frozen & ~is_ocean,
        np.minimum(new_glacier_depth + frozen_river_ice, GLACIER_MAX_DEPTH_M),
        new_glacier_depth,
    )
    ice.ledger.frozen_river_ice_m3 = float(np.dot(with_river_ice - new_glacier_depth, areas))
    unclipped = np.where(water_frozen & ~is_ocean, new_glacier_depth + frozen_river_ice, new_glacier_depth)
    ice.ledger.frozen_river_capped_m3 = float(np.dot(unclipped - with_river_ice, areas))
    ice_sheet.warn_if_river_ice_capped(ice.ledger)
    new_glacier_depth = with_river_ice

    # is_river isn't known yet -- it needs this step's *final* lake_depth (below) to exclude
    # flooded nodes, so fields is built with a placeholder here and finished after.
    fields = HydrologyFields(
        points, elevation, is_ocean, neighbor_idx, flow_target, flow_accum, water_deposited, filled_elevation, spill_target,
        np.zeros(n, dtype=bool), lake_depth_adjusted, new_glacier_depth, plates_in_order, ice_flow_target=ice_flow_target,
        channel_width=prev_channel_width,
        breach_notch_m=breach.notch_m if breach is not None else np.zeros(0),
        ice=ice,
        should_spill=should_spill,
    )
    # Resolves the *same* forest built early (above, for spill routing) against this step's
    # routed physical water -- not a second lakes.step_lakes call, which would rebuild the
    # hierarchy from scratch a second time and risk it disagreeing with what spill routing just
    # used (see lakes.resolve_lakes's own docstring).
    # A spilling lake's sink passes its water on, so its inflow is the flow through the sink,
    # not what settles there. That through-flow also carries on to any lake downstream
    # without the upstream lake's own evaporation taken off. That overfeeds lake chains a
    # little, but only through lakes already full and spilling. Subtracting it would need
    # lakes resolved during routing (follow-up on #297).
    # Ocean targets hold the land-to-ocean flux in `settled_m3`, since ocean nodes are sinks
    # rather than members of the river routing graph. Expose that mouth discharge on the target
    # too for shelf geology; river networks themselves still report their last land node.
    fields.discharge_m3_s = np.where(is_ocean, settled_m3, through_m3) / SECONDS_PER_YEAR
    fields.breach_surge = breach_surge
    fields.breach_surge_accum = breach_surge_accum
    lake_climate = lakes.LakeClimate(
        inflow_m3_per_yr=settled_m3 + np.where(should_spill, through_m3, 0.0),
        evaporation_m_per_yr=evaporation_mm / 1000.0,
    )
    fields.lake_depth, fields.silt_deposited, fields.is_sea, fields.lake_events = lakes.resolve_lakes(
        forest, elevation, lake_depth_adjusted, years, lake_frozen, lake_climate, hypsometry
    )
    fields.lake_forest = forest
    fields.hypsometry = hypsometry
    # Record only (erosion.py folds `silt_deposited` into real `elevation` and caps the
    # persisted silt_depth at the mobile cover holding it).
    fields.silt_depth = prev_silt_depth + fields.silt_deposited

    # A lake, however wide it's flooded, is never itself "a river" -- excluding it (not just
    # its own former sink node, every flooded node) is what makes a river actually end at a
    # lake's shore instead of its classification jumping straight across to whatever's on
    # the far side of the water.
    land = ~is_ocean
    is_lake = fields.lake_depth > LAKE_MIN_VISIBLE_DEPTH_M
    is_river = np.zeros(n, dtype=bool)
    if np.any(land) and np.any(fields.discharge_m3_s[land] > 0):
        threshold = np.percentile(fields.discharge_m3_s[land], RIVER_FLOW_PERCENTILE)
        is_river = land & (fields.discharge_m3_s > 0) & (fields.discharge_m3_s >= threshold) & ~is_lake

    # Rivers commonly begin right where a glacier's meltwater first emerges from the ice --
    # a fresh headwater has, by definition, no upstream tributaries yet, so it can easily fail
    # to clear the ordinary discharge threshold above on its own even though it's a
    # perfectly real river source. A land node not itself under visible ice but with at least
    # one neighbor that is (GLACIER_VISIBLE_DEPTH_M, the same "real ice, not just a trace"
    # threshold rendering already uses) sits right at a glacier's edge; if it's actually
    # carrying real outflow this step (`discharge_m3_s > 0`, i.e. meltwater genuinely reached here),
    # it's marked is_river directly, giving that stream a visible source at the ice margin
    # instead of only appearing once downstream tributaries happen to push it over the
    # percentile cut.
    has_glacier_neighbor = np.any(new_glacier_depth[neighbor_idx] >= GLACIER_VISIBLE_DEPTH_M, axis=1)
    is_glacial_source = land & ~is_lake & (new_glacier_depth < GLACIER_VISIBLE_DEPTH_M) & has_glacier_neighbor & (fields.discharge_m3_s > 0)
    fields.is_river = _extend_river_mask_downstream(
        is_river | is_glacial_source,
        fields.flow_target,
        fields.discharge_m3_s,
        fields.is_ocean,
        is_lake,
    )
    return fields


@dataclass
class RiverInfo:
    """One connected drainage network within `HydrologyFields.is_river` -- the River
    Inspector's unit of selection. There's no such grouping concept anywhere else in this
    module; `is_river` is otherwise just a flat per-node boolean. `member_idx` indexes back
    into the *same* HydrologyFields arrays (points/elevation/discharge_m3_s/flow_target/...) this
    was built from, so a caller (main.py) resolves geometry/segments straight from those
    rather than this dataclass duplicating them."""

    member_idx: np.ndarray  # (M,) indices into the source HydrologyFields arrays
    mouth_idx: int  # index into the same arrays -- the network's downstream mouth member
    mouth_type: str  # "ocean" | "lake" | "other" (a dry interior sink)
    flow_rate: float  # physical discharge_m3_s at the mouth
    speed: float  # compute_river_speed at the mouth
    num_tributaries: int  # separate headwater branches feeding this network, see group_rivers
    length_m: float | None = None  # summed physical length; None only for legacy fixtures


def _spill_hops(fields: HydrologyFields, sources: np.ndarray) -> np.ndarray:
    """Which of `sources`' flow edges are spill hops: a spilling node routed to a node that
    isn't its k-NN neighbour. A cache pickled before `should_spill` was stored falls back to
    `flow_target == spill_target`, which misses only the united-lake pooling hops."""
    targets = fields.flow_target[sources]
    spilling = fields.should_spill
    if spilling is None or len(spilling) != len(fields.flow_target):
        spilling = fields.flow_target == fields.spill_target
    hop = (targets >= 0) & spilling[sources]
    hop_sources, hop_targets = sources[hop], targets[hop]
    neighbour = (fields.neighbor_idx[hop_sources] == hop_targets[:, None]).any(axis=1)
    neighbour |= (fields.neighbor_idx[hop_targets] == hop_sources[:, None]).any(axis=1)
    hop[hop] = ~neighbour
    return hop


# A backstop on how many nodes one spill-path search may settle before it falls back to the
# straight segment (`_spill_path`); the climb ceiling usually stops a failing search sooner.
SPILL_PATH_MAX_NODES = 50_000


@dataclass
class _SpillPathCache:
    """Each river spill hop's drawn segments, worked out once per HydrologyFields (#520): the
    map, `group_rivers`, the River Inspector and the Lake Inspector all ask for the same hops,
    several times per request. `HydrologyFields.__getstate__` leaves it out of saves, so a save
    never carries paths from older drawing rules."""

    adjacency: tuple[np.ndarray, np.ndarray]
    # source -> (from_idx, to_idx, crosses_lake)
    paths: dict[int, tuple[np.ndarray, np.ndarray, bool]] = field(default_factory=dict)


def _drawn_spill_path(fields: HydrologyFields, adjacency: tuple[np.ndarray, np.ndarray], source: int) -> tuple[np.ndarray, np.ndarray, bool]:
    """One hop's drawn segments and whether it crosses a lake, see `flow_segments`."""
    target = int(fields.flow_target[source])
    path = np.array(_spill_path(fields, adjacency, source, target) or [source, target], dtype=np.int64)
    start, end = path[:-1], path[1:]
    wet = fields.lake_depth[path] > LAKE_MIN_VISIBLE_DEPTH_M
    # Every edge up to the first one under visible water.
    dry = np.logical_and.accumulate(~(wet[:-1] & wet[1:]))
    return start[dry], end[dry], not np.all(dry)


def _spill_path_cache(fields: HydrologyFields, hop_sources: list[int]) -> _SpillPathCache:
    """The fields' cache, with `hop_sources` filled in. The first call also fills in every
    river node's hop, so later callers asking about any river find it there."""
    cache = fields.__dict__.get("_spill_path_cache")
    if cache is None:
        cache = _SpillPathCache(_undirected_neighbours(fields.neighbor_idx))
        river_idx = np.flatnonzero(fields.is_river)
        hop_sources = hop_sources + river_idx[_spill_hops(fields, river_idx)].tolist()
        fields.__dict__["_spill_path_cache"] = cache
    for source in hop_sources:
        if source not in cache.paths:
            cache.paths[source] = _drawn_spill_path(fields, cache.adjacency, source)
    return cache


def _undirected_neighbours(neighbor_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """`(indptr, indices)` of the k-NN graph with every edge in both directions. k-NN is
    directed, and a lake's rim edge is taken from either end (`lakes.build_lake_hierarchy`), so
    the way back over a rim, or from the sea up to its shore, may exist only as a reverse edge."""
    n, k = neighbor_idx.shape
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    graph = coo_matrix((np.ones(2 * n * k, dtype=np.int8), (np.concatenate([rows, cols]), np.concatenate([cols, rows]))), shape=(n, n)).tocsr()
    return graph.indptr, graph.indices


def _spill_path_ceiling(fields: HydrologyFields, source: int) -> float:
    """The highest a spill path from `source` can need to climb: its routed spill level
    (`filled_elevation`) plus the deepest breach notch, since a notched passage
    (`breaching.interface_pass_elevation`) routes the spill below the rim's raw elevation."""
    filled = float(fields.filled_elevation[source]) if len(fields.filled_elevation) == len(fields.elevation) else np.inf
    return filled + breaching.BREACH_MAX_NOTCH_M if np.isfinite(filled) else np.inf


@njit(cache=True)
def _spill_pass_level(indptr, indices, elevation, is_ocean, source, target, ceiling, max_nodes):
    """`_spill_path`'s minimax search, compiled: the lowest pass between `source` and
    `target`, or inf if the search gives up."""
    n = len(elevation)
    best = np.full(n, np.inf)
    best[source] = elevation[source]
    heap = [(elevation[source], np.int64(source))]
    level = np.inf
    found = False
    explored = 0
    while len(heap) > 0:
        pass_m, node = heapq.heappop(heap)
        if node == target:
            level = pass_m
            found = True
            break
        if pass_m > best[node]:
            continue
        # The heap pops passes in rising order, so nothing left can come in under the ceiling.
        explored += 1
        if pass_m > ceiling or explored > max_nodes:
            break
        for j in range(indptr[node], indptr[node + 1]):
            nxt = np.int64(indices[j])
            if is_ocean[nxt] and nxt != target:
                continue
            candidate = max(pass_m, elevation[nxt])
            if candidate < best[nxt]:
                best[nxt] = candidate
                heapq.heappush(heap, (candidate, nxt))
    return level if found else np.inf


@njit(cache=True)
def _fewest_hops_path(indptr, indices, elevation, is_ocean, source, target, level):
    """`_spill_path`'s breadth-first search, compiled: the fewest-hop path from `source` to
    `target` through nodes no higher than `level`."""
    n = len(elevation)
    came_from = np.full(n, -2, dtype=np.int64)
    came_from[source] = -1
    queue = np.empty(n, dtype=np.int64)
    queue[0] = source
    head, tail = 0, 1
    while head < tail and came_from[target] == -2:
        node = queue[head]
        head += 1
        for j in range(indptr[node], indptr[node + 1]):
            nxt = np.int64(indices[j])
            if came_from[nxt] == -2 and (nxt == target or not is_ocean[nxt]) and elevation[nxt] <= level:
                came_from[nxt] = node
                queue[tail] = nxt
                tail += 1
    length = 1
    node = target
    while came_from[node] >= 0:
        node = came_from[node]
        length += 1
    path = np.empty(length, dtype=np.int64)
    node = target
    for i in range(length - 1, -1, -1):
        path[i] = node
        node = came_from[node]
    return path


def _spill_path(fields: HydrologyFields, adjacency: tuple[np.ndarray, np.ndarray], source: int, target: int) -> list[int]:
    """The nodes a spill from `source` to `target` crosses over `adjacency`'s neighbour edges:
    the fewest hops among paths that never rise above the lowest pass between them (found
    first, by a minimax search), so the path runs across the basin floor and out over its rim.
    Ocean nodes other than `target` are excluded.

    Empty if `target` isn't reachable without climbing past `_spill_path_ceiling`, or within
    SPILL_PATH_MAX_NODES: an unreachable target would otherwise search every land node."""
    indptr, indices = adjacency
    elevation = np.asarray(fields.elevation, dtype=np.float64)
    is_ocean = np.asarray(fields.is_ocean, dtype=np.bool_)
    level = _spill_pass_level(
        indptr, indices, elevation, is_ocean, source, target, float(_spill_path_ceiling(fields, source)), SPILL_PATH_MAX_NODES
    )
    if not np.isfinite(level):
        return []
    return _fewest_hops_path(indptr, indices, elevation, is_ocean, source, target, level).tolist()


def flow_segments(fields: HydrologyFields, sources: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The river line segments drawn for the flow edges leaving `sources`, for the map and the
    River Inspector alike. An ordinary edge is one segment to `flow_target`. A spill hop
    (`_spill_hops`) can span a whole lake (#520: a 1,669 km edge drawn across an inland sea),
    so it's drawn along its spill path (`_spill_path`) instead. A path that reaches visible
    lake water is drawn only up to it: the river ends at the lake, and the outflow is its own
    network beyond the rim. Routing and discharge are unchanged; only the drawn geometry is.

    Returns `(from_idx, to_idx, owner, crosses_lake)`: one entry per segment, `owner` being the
    position in `sources` whose edge it draws, and per source, whether its edge crosses visible
    lake water, which `group_rivers` treats as the river's end. A hop with no spill path at all
    keeps its straight segment rather than ending the river."""
    sources = np.asarray(sources, dtype=np.int64)
    targets = fields.flow_target[sources]
    hop = _spill_hops(fields, sources)
    plain = np.flatnonzero((targets >= 0) & ~hop)
    from_parts, to_parts, owner_parts = [sources[plain]], [targets[plain]], [plain]
    crosses_lake = np.zeros(len(sources), dtype=bool)
    hop_positions = np.flatnonzero(hop).tolist()
    if hop_positions:
        cache = _spill_path_cache(fields, sources[hop_positions].tolist())
        for position in hop_positions:
            start, end, crosses_lake[position] = cache.paths[int(sources[position])]
            from_parts.append(start)
            to_parts.append(end)
            owner_parts.append(np.full(len(start), position, dtype=np.int64))
    return np.concatenate(from_parts), np.concatenate(to_parts), np.concatenate(owner_parts), crosses_lake


def group_rivers(fields: HydrologyFields) -> list[RiverInfo]:
    """Groups the flat `is_river` node mask into distinct connected drainage networks, via
    union-find over `flow_target` edges restricted to nodes that are `is_river` on *both*
    ends. Since river selection uses a percentile of physical discharge, joining neighboring
    selected nodes preserves the current flow topology without using the legacy accumulated
    rainfall field.

    A network's **mouth** is its downstream-most member. Its `flow_target` chain determines
    whether it empties into open ocean, a lake (excluded from `is_river` by the `~is_lake`
    term below), or a true dead-end sink. A cycle in a hand-built fixture falls back to the
    member with greatest physical discharge.

    **mouth_type** reads `lake_depth`/`is_ocean` off the chain starting at
    *`flow_target[mouth_idx]`* -- the node(s) the mouth actually empties into -- not off the
    mouth itself: `is_river` above always excludes lake nodes (`~is_lake`), so the mouth of a
    real river can never itself carry standing water, and checking the mouth's own
    (structurally always-dry) lake_depth here made `"lake"` effectively unreachable and
    silently mislabeled every river that actually ends at a lake as `"ocean"` instead (any real
    flow_target at all, lake or open sea, used to read as `"ocean"` outright -- confirmed
    directly: a hand-built fixture modeling a river ending at a lake only ever passed because it
    fabricated `is_river=True` on the flooded node itself, a state real `compute_hydrology`
    output can never produce). The chain walks *past* `flow_target[mouth_idx]` rather than
    stopping there because a wide, gently-sloped basin floor routinely crosses several more dry
    hops of real (currently-above-the-shoreline) catchment before actually reaching standing
    water -- stopping at one hop mislabeled a river ending in an almost-full lake as `"other"`
    (confirmed directly against a real save: a 636-node basin sitting at ~96% of its own rim,
    with only its lowest couple of nodes actually underwater, so the immediate mouth_target read
    as neither ocean nor lake even though the same chain reached real lake water two hops
    later). `"other"` is reserved for what it's meant to mean -- a genuine, still-unresolved dry
    interior sink (the chain runs out at `flow_target == -1` without ever crossing `is_ocean` or
    `LAKE_MIN_VISIBLE_DEPTH_M`, not yet spilling) -- which `lakes.step_lakes` already grows into
    a real lake node by node every step regardless of this classification; once that basin's own
    water crosses `LAKE_MIN_VISIBLE_DEPTH_M` anywhere along the chain this same river correctly
    flips to `"lake"` on its own, no separate fix needed there.

    A spill hop that crosses visible lake water (`flow_segments`) is the one edge between river
    nodes that doesn't join them: the upstream network ends there with `mouth_type` `"lake"`,
    and the outflow is a separate network starting beyond the rim (#520). `length_m` sums the
    drawn segments, so a hop counts its spill path, not the straight line.

    **num_tributaries** (an original definition, with no existing precedent in this codebase
    to follow): counts each member's in-network in-degree (how many *other* members flow
    directly into it); a member with in-degree 0 is a headwater -- a separate source stream
    with nothing upstream of it in this network. A single unbranched channel has exactly one
    headwater (itself) and therefore zero tributaries; each additional headwater is one more
    distinct stream that joins the network somewhere along its course, so
    `num_tributaries = headwater_count - 1`."""
    river_idx = np.nonzero(fields.is_river)[0]
    if len(river_idx) == 0:
        return []

    river_list = river_idx.tolist()
    river_set = set(river_list)
    flow_target_list = fields.flow_target.tolist()

    seg_from, seg_to, seg_owner, crosses_lake = flow_segments(fields, river_idx)
    # A spill across visible lake water ends the river at the lake, as a river meeting a lake's
    # shore does; the outflow is its own network from the outlet on (#520).
    ends_at_lake = set(river_idx[crosses_lake].tolist())

    parent = {i: i for i in river_list}

    def find(x: int) -> int:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:
            parent[x], x = root, parent[x]
        return root

    for i in river_list:
        target = flow_target_list[i]
        if target in river_set and i not in ends_at_lake:
            ra, rb = find(i), find(target)
            if ra != rb:
                parent[ra] = rb

    groups: dict[int, list[int]] = {}
    for i in river_list:
        groups.setdefault(find(i), []).append(i)

    slope_to_target = _slope_to_flow_target(fields.points, fields.elevation, fields.flow_target)
    edge_length_m = np.zeros(len(fields.points), dtype=float)
    np.add.at(
        edge_length_m,
        river_idx[seg_owner],
        geometry.angular_distance(fields.points[seg_from], fields.points[seg_to]) * PLANET_RADIUS_KM * 1000.0,
    )
    discharge_m3_s = fields.discharge_m3_s
    if discharge_m3_s is None or len(discharge_m3_s) != len(fields.flow_accum):
        # Hand-built fixtures and pre-step-1 caches have only the old stacked-mm/yr quantity.
        # Loaded saves with that cache are invalidated by main.load_world before reaching this
        # fallback; it keeps internal callers with legacy fixtures working.
        discharge_m3_s = flow_accum_mm_yr_to_discharge_m3_s(fields.flow_accum)

    rivers: list[RiverInfo] = []
    for members in groups.values():
        member_idx = np.array(members, dtype=np.int64)
        member_set = set(members)
        terminal_members = [member for member in members if flow_target_list[member] not in member_set or member in ends_at_lake]
        mouth_candidates = terminal_members or members
        # A cycle in a hand-built fixture has no terminal; use its greatest-discharge member.
        mouth_idx = max(mouth_candidates, key=lambda member: discharge_m3_s[member])

        # See this function's own docstring for why this reads the *target*'s state, not the
        # (structurally always-dry) mouth's own -- and why it walks the chain past that first
        # hop rather than stopping there. A wide, gently-sloped basin floor routinely crosses
        # several more dry hops (real, currently-above-the-shoreline catchment, not a bug -- see
        # lakes.py's own "members include dry higher ground" docstring) before actually reaching
        # standing water, so checking only mouth_target's own lake_depth mislabeled a river
        # ending in an almost-full lake as "other" (confirmed directly against a real save: a
        # 636-node basin sitting at ~96% of its own rim, with only its lowest couple of nodes
        # actually underwater -- mouth_target was still dry basin floor, but the same chain
        # reached real lake water two hops later). `visited` guards this since it walks live
        # per-step state rather than trusting the usual downhill-monotonic invariant to
        # terminate.
        mouth_type = "lake" if mouth_idx in ends_at_lake else "other"
        node = -1 if mouth_idx in ends_at_lake else int(fields.flow_target[mouth_idx])
        visited: set[int] = set()
        while node >= 0 and node not in visited:
            if fields.is_ocean[node]:
                mouth_type = "ocean"
                break
            if fields.lake_depth[node] > LAKE_MIN_VISIBLE_DEPTH_M:
                mouth_type = "lake"
                break
            visited.add(node)
            node = int(fields.flow_target[node])

        in_degree = {m: 0 for m in members}
        for m in members:
            t = flow_target_list[m]
            if t in member_set and m not in ends_at_lake:
                in_degree[t] += 1
        headwater_count = sum(1 for m in members if in_degree[m] == 0)

        length_m = float(np.sum(edge_length_m[member_idx]))

        rivers.append(
            RiverInfo(
                member_idx=member_idx,
                mouth_idx=mouth_idx,
                mouth_type=mouth_type,
                flow_rate=float(discharge_m3_s[mouth_idx]),
                speed=float(compute_river_speed(slope_to_target[mouth_idx], discharge_m3_s[mouth_idx])),
                num_tributaries=max(0, headwater_count - 1),
                length_m=length_m,
            )
        )
    return rivers


def river_at(fields: HydrologyFields, rivers: list[RiverInfo], query_xyz: np.ndarray) -> int | None:
    """Nearest-river-node lookup across every network's own member points -- the River
    Inspector's click hit-test, mirroring plates.nearest_plate_id. Returns an index into
    `rivers` (not a persistent id: rivers are regrouped fresh from the current
    hydrology_cache on every request, see main.py), or None if there are no rivers this
    step."""
    if not rivers:
        return None
    all_idx = np.concatenate([r.member_idx for r in rivers])
    owner = np.concatenate([np.full(len(r.member_idx), i, dtype=np.int64) for i, r in enumerate(rivers)])
    tree = cKDTree(fields.points[all_idx])
    _, nearest = tree.query(query_xyz.reshape(1, 3), k=1)
    return int(owner[int(nearest[0])])
