"""Rivers, lakes, and glaciers: flow routing, basin/lake detection, ice accumulation/melt/
flow, and downstream flow accumulation over the world's current node cloud.

A fixed grid's regular 8-neighbor adjacency would give flow routing, depression filling, and
downstream accumulation a natural substrate, but mantle-bloom's world is instead an irregular
per-plate node cloud. All three core flow-routing algorithms -- steepest-descent flow
direction, priority-flood basin-spill, and elevation-ordered downstream accumulation -- turn
out not to actually need a *grid*, only a *graph*: this module builds one via a whole-world
k-nearest-neighbor query (the same
technique reassign.py/erosion.py already use for their own whole-world passes), then runs
the same algorithms directly on it. Glacier flow (accumulated ice moving one hop downhill
per step) reuses the same `flow_target` graph water uses.

**Persistence.** A grid-based approach where plates move relative to fixed cells would need a
persistent field like channel_depth to be deliberately advected (semi-Lagrangian, every step)
to keep following the crust -- mantle-bloom's elevation-line nodes sidestep that entirely,
since they already rotate exactly with their own plate. lake_depth and glacier_depth, stored
as ordinary parallel arrays on ElevationLine right alongside elevation itself (see plates.py),
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
(`lakes.Lake`, `lakes.build_lake_hierarchy`, `lakes.step_lakes`). An earlier version of this
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
from collections import deque
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, lakes
from .plates import PLANET_RADIUS_KM, ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

# Same 8-neighbor D8 structure as a regular grid's adjacency -- the graph every algorithm
# below runs on.
FLOW_NEIGHBOR_COUNT = 8

# Top decile of land flow_accum counts as "a river" -- used for both rendering and the
# major_river_fraction stat.
RIVER_FLOW_PERCENTILE = 90.0

# How established a downhill neighbor's own channel needs to be before flow direction
# prefers it over the literal steepest candidate -- see _compute_flow_direction. Deliberately
# far below erosion.CHANNEL_BOOST_REFERENCE_M (200m, where channel_boost saturates): this
# only needs to detect "a real channel exists here at all," not "a mature one," so routing
# starts favoring it as soon as there's genuine, non-noise incision.
CHANNEL_PREFERENCE_THRESHOLD_M = 5.0

# Threshold above which a node's lake_depth counts as "a visible lake" -- used both for
# rendering (render_image.py, coastline.py) and for excluding flooded nodes from is_river/
# group_rivers below. Lake growth/evaporation tuning itself (LAKE_FILL_RATE and friends) now
# lives in lakes.py, alongside the Lake tree that actually computes it -- see this module's own
# docstring and lakes.py's own module docstring for why.
LAKE_MIN_VISIBLE_DEPTH_M = 1.0

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
LAKE_BREACH_EROSION_COEFFICIENT = 4000.0

# Glaciers use flat per-Myr rates (linear in dt_myr), not an exponential-decay formula like
# the lake evaporation rate above -- so they don't blow up the same way at mantle-bloom's
# larger typical step size -- a bigger step just melts/grows proportionally more, same
# character as every other elevation delta in this codebase (boundary.py's uplift/trench
# rates are also flat per-Myr rates applied directly). No seasons modeled here either
# (consistent with mantle-bloom's climate.py generally) -- GLACIER_ACCUMULATION_TEMP_C
# being well below 0C is what keeps this meaning *permanent* accumulation, not a place with
# ordinary seasonal snow that would fully melt over a real year.
GLACIER_ACCUMULATION_TEMP_C = -10.0
GLACIER_ACCUMULATION_RATE = 0.02
GLACIER_MELT_RATE_M_PER_MYR = 400.0
GLACIER_MELT_REFERENCE_DEGREES_C = 20.0
GLACIER_MELT_MAX_FACTOR = 3.0
GLACIER_FLOW_RATE_PER_MYR = 0.15
GLACIER_MAX_FLOW_FRACTION = 0.5
GLACIER_VISIBLE_DEPTH_M = 10.0

# River speed -- a stylized, unitless quantity (faster where slope is steeper and where more
# water has accumulated); meaningful only relative to other nodes in the same world, not a
# real physical speed. Lives here rather than erosion.py (which only *consumes* it for its
# deposition threshold) so the River Inspector's per-river mouth speed (see group_rivers)
# can share it without a hydrology.py -> erosion.py reverse import.
RIVER_SPEED_COEFFICIENT = 4.0
RIVER_SPEED_DISCHARGE_EXPONENT = 0.2


@dataclass
class HydrologyFields:
    """Everything derived from one whole-world flow-routing pass, all shape (N,) aligned
    with `points`/`elevation` -- one flat array per field over the irregular node cloud, in
    place of a grid. `line_refs` is (plate, line_index, start, end) per line, letting a
    caller slice any of these flat arrays back onto a specific line's own node range.
    `lake_depth`/`glacier_depth` are already this step's *final* values (state-transition
    owned by this module, see module docstring) -- a caller doesn't need to call anything
    else to get the up-to-date depth."""

    points: np.ndarray
    elevation: np.ndarray
    is_ocean: np.ndarray
    neighbor_idx: np.ndarray  # (N, FLOW_NEIGHBOR_COUNT)
    flow_target: np.ndarray  # (N,) int, index into these same arrays; -1 = sink or ocean
    flow_accum: np.ndarray  # (N,) precipitation-weighted downstream water accumulation
    water_deposited: np.ndarray  # (N,) water settled here (nonzero only at a true sink)
    filled_elevation: np.ndarray  # (N,) minimal-bottleneck elevation to reach open ocean
    spill_target: np.ndarray  # (N,) one-hop basin-escape neighbor index; -1 = none/ocean
    is_river: np.ndarray  # (N,) bool -- top RIVER_FLOW_PERCENTILE of land flow_accum
    lake_depth: np.ndarray  # (N,) this step's final standing lake depth
    glacier_depth: np.ndarray  # (N,) this step's final accumulated ice depth
    line_refs: list[tuple[Plate, int, int, int]]
    # Human-readable lake merge/split transition messages from lakes.step_lakes, for the
    # caller (erosion.py) to log via World.log_event -- see compute_hydrology. Defaulted
    # (rather than a required positional field) so every existing HydrologyFields call site
    # that predates lakes.py, including hand-built test fixtures, keeps working unchanged.
    lake_events: list[str] = field(default_factory=list)
    # This step's final silt depth (lakes.step_lakes's own state transition, same "already the
    # final value" character as lake_depth/glacier_depth above) -- also defaulted, see
    # lake_events' own comment for why.
    silt_depth: np.ndarray = field(default_factory=lambda: np.zeros(0))


def _gather_nodes(
    world: "World",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position, elevation, prior lake_depth, prior glacier_depth, prior
    channel_depth, prior silt_depth, and whether it's ocean (elevation <= world.sea_level_m, the
    live-adjustable sea-level convention used everywhere else in this codebase -- see
    World.sea_level_m), concatenated, alongside (plate, line_index, start, end) references
    -- same shape as erosion.py's/bathymetry.py's own _gather_nodes.
    channel_depth is read-only here (erosion.py still owns growing it) -- flow direction
    just needs to know where an established channel already is, see
    _compute_flow_direction."""
    points_list, elev_list, lake_list, glacier_list, channel_list, silt_list = [], [], [], [], [], []
    line_refs: list[tuple[Plate, int, int, int]] = []
    offset = 0
    for plate in world.plates:
        for line_index, line in enumerate(plate.lines):
            n = len(line.theta)
            if n == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            elev_list.append(line.elevation)
            lake_list.append(line.lake_depth)
            glacier_list.append(line.glacier_depth)
            channel_list.append(line.channel_depth)
            silt_list.append(line.silt_depth)
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        empty = np.zeros(0)
        return np.zeros((0, 3)), empty, empty, empty, empty, empty, np.zeros(0, dtype=bool), []
    points = np.concatenate(points_list, axis=0)
    elevation = np.concatenate(elev_list, axis=0)
    prev_lake_depth = np.concatenate(lake_list, axis=0)
    prev_glacier_depth = np.concatenate(glacier_list, axis=0)
    prev_channel_depth = np.concatenate(channel_list, axis=0)
    prev_silt_depth = np.concatenate(silt_list, axis=0)
    is_ocean = elevation <= world.sea_level_m
    return points, elevation, prev_lake_depth, prev_glacier_depth, prev_channel_depth, prev_silt_depth, is_ocean, line_refs


def _build_neighbor_graph(points: np.ndarray) -> np.ndarray:
    n = len(points)
    k = min(FLOW_NEIGHBOR_COUNT, n - 1)
    tree = cKDTree(points)
    _, neighbor_idx = tree.query(points, k=k + 1)
    return neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0


def _compute_basin_spill(elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Priority-flood depression filling: a multi-source Dijkstra seeded from every ocean
    node, relaxing by *max* (the highest point a path is forced to cross) rather than by
    sum -- a minimax path cost, not a shortest path. The same priority-flood algorithm
    generalizes cleanly from grid 8-neighbor adjacency to the k-NN graph's edges. Returns
    (filled_elevation, spill_target): the minimal-bottleneck elevation to
    reach open ocean from each node, and a one-hop escape neighbor toward it (-1 for ocean
    or a node with no path to any ocean at all)."""
    n = len(elevation)
    elevation_list = elevation.tolist()  # plain-list access is much faster than per-element
    is_ocean_list = is_ocean.tolist()  # numpy indexing inside this loop's hot path
    neighbor_list = neighbor_idx.tolist()

    cost = [np.inf] * n
    spill_target = [-1] * n
    visited = [False] * n

    heap: list[tuple[float, int]] = []
    for i in range(n):
        if is_ocean_list[i]:
            cost[i] = elevation_list[i]
            heap.append((cost[i], i))
    heapq.heapify(heap)

    while heap:
        d, i = heapq.heappop(heap)
        if visited[i]:
            continue
        visited[i] = True
        for j in neighbor_list[i]:
            if visited[j] or is_ocean_list[j]:
                continue
            neighbor_elev = elevation_list[j]
            candidate = neighbor_elev if neighbor_elev > d else d
            if candidate < cost[j]:
                cost[j] = candidate
                if candidate > d:
                    spill_target[j] = i
                else:
                    inherited = spill_target[i]
                    spill_target[j] = inherited if inherited >= 0 else i
                heapq.heappush(heap, (candidate, j))

    return np.array(cost), np.array(spill_target, dtype=np.int64)


def _lake_component_sizes(is_lake: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Per-node size (node count) of the connected `is_lake` component containing that node,
    0 for every non-lake node -- a same-lake grouping via union-find over `neighbor_idx`'s
    edges, treated as **undirected** here even though the k-NN graph itself isn't strictly
    symmetric (a node isn't guaranteed to be among its own neighbor's k nearest). Unlike
    `_flood_fill_lake_extent`'s directed traversal (fine for its own purpose, spreading a
    water *level* outward from a seed), a size count needs two lake nodes joined by only a
    one-directional edge to still read as the same lake, not be split apart by an accident of
    which node happens to list the other as one of its own nearest neighbors."""
    n = len(is_lake)
    if not np.any(is_lake):
        return np.zeros(n, dtype=np.int64)

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
    roots = np.array([find(i) for i in lake_idx.tolist()])
    _, inverse, counts = np.unique(roots, return_inverse=True, return_counts=True)

    sizes = np.zeros(n, dtype=np.int64)
    sizes[lake_idx] = counts[inverse]
    return sizes


def _compute_flow_direction(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    prev_lake_depth: np.ndarray,
    filled_elevation: np.ndarray,
    spill_target: np.ndarray,
    prev_channel_depth: np.ndarray,
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
    water surface (elevation + the *previous* step's lake_depth, a one-step-lagged "memory"
    of the water surface) has already reached its
    basin's true spill point, in which case it redirects to spill_target instead of staying
    a dead-end sink forever. Ocean nodes are never routed (they're a destination, not a
    source).

    Also returns `should_spill`, the exact per-node mask of which sinks got redirected this
    way -- `compute_hydrology` needs this specific mask (not something recomputed after the
    fact from `flow_target == spill_target`, which could false-positive on an ordinary
    steepest-descent node that happens to coincide with a neighboring basin's spill_target) to
    know which lakes just found a real outlet, for the breach-erosion term (see its own
    docstring in `compute_hydrology`)."""
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
    water_surface = elevation + prev_lake_depth
    should_spill = is_sink & (water_surface >= filled_elevation)
    flow_target = np.where(should_spill, spill_target, flow_target).astype(np.int64)
    flow_target = np.where(is_ocean, -1, flow_target).astype(np.int64)
    return flow_target, should_spill


def _slope_to_flow_target(points: np.ndarray, elevation: np.ndarray, flow_target: np.ndarray) -> np.ndarray:
    """Dimensionless rise/run from each node to its own flow_target (0 for a node with no
    target) -- the slope glacier flow moves along, distinct from erosion.py's own
    SLOPE_NEIGHBOR_COUNT=4 "lowest of k neighbors" slope: glacier flow specifically follows
    `flow_target` (this module's own k=FLOW_NEIGHBOR_COUNT graph), not erosion's separate,
    smaller neighborhood."""
    has_target = flow_target >= 0
    safe_target = np.where(has_target, flow_target, 0)
    drop = np.where(has_target, elevation - elevation[safe_target], 0.0)
    run_m = geometry.angular_distance(points, points[safe_target]) * PLANET_RADIUS_KM * 1000.0
    run_m = np.maximum(run_m, 1.0)
    return np.where(has_target, np.clip(drop, 0.0, None) / run_m, 0.0)


def compute_river_speed(slope: np.ndarray, flow_accum: np.ndarray) -> np.ndarray:
    """Stylized, unitless river speed -- faster where slope is steeper and where more water
    has accumulated; meaningful only relative to other nodes in the same world, not a real
    physical speed."""
    return RIVER_SPEED_COEFFICIENT * np.sqrt(np.clip(slope, 0.0, None)) * np.power(np.clip(flow_accum, 0.0, None), RIVER_SPEED_DISCHARGE_EXPONENT)


def route_downstream(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    flow_target: np.ndarray,
    source_amount: np.ndarray,
    retain_fraction: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Single forward sweep over land nodes in elevation-descending order, accumulating
    `source_amount` downstream along `flow_target` edges. Correct in one pass because every
    node's target is guaranteed strictly lower in elevation, so it's always visited *later*
    in this same order. No in-transit river evaporation term (a loss_fraction consumed along
    the way) is modeled here -- not asked for, and this module already drops
    temperature-driven effects on hydrology, see LAKE_EVAPORATION_* above. Returns
    (through_flux, deposited): each node's own
    accumulated flux passing through it, and how much settled there (from `retain_fraction`
    and/or reaching a sink/ocean)."""
    n = len(elevation)
    through_flux = np.where(~is_ocean, source_amount, 0.0).astype(np.float64).tolist()
    deposited = [0.0] * n
    retain = (retain_fraction if retain_fraction is not None else np.zeros(n)).tolist()
    is_ocean_list = is_ocean.tolist()
    flow_target_list = flow_target.tolist()

    land_indices = np.nonzero(~is_ocean)[0]
    order = land_indices[np.argsort(-elevation[land_indices])].tolist()

    for i in order:
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

    return np.array(through_flux), np.array(deposited)


def _update_glaciers(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    flow_target: np.ndarray,
    slope_to_target: np.ndarray,
    prev_glacier_depth: np.ndarray,
    frozen_precip: np.ndarray,
    frozen_from_lake: np.ndarray,
    is_accumulating: np.ndarray,
    temperature: np.ndarray,
    years: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Grows/melts/flows glacier ice, generalized from grid 8-neighbor D8 routing to this
    module's flow_target graph. Returns (new_glacier_depth, melt): melt feeds back into this
    step's water source (see compute_hydrology) as real meltwater, not a separate accounting
    bucket -- glacial melt becomes real river discharge rather than a parallel, disconnected
    accounting path."""
    years_myr = years / 1_000_000.0
    accumulation = np.where(is_accumulating, GLACIER_ACCUMULATION_RATE * frozen_precip * years_myr, 0.0)
    depth_before_flow = prev_glacier_depth + frozen_from_lake + accumulation

    melt_factor = np.clip(
        (temperature - GLACIER_ACCUMULATION_TEMP_C) / GLACIER_MELT_REFERENCE_DEGREES_C, 0.0, GLACIER_MELT_MAX_FACTOR
    )
    melt = np.minimum(GLACIER_MELT_RATE_M_PER_MYR * melt_factor * years_myr, depth_before_flow)
    depth_after_melt = depth_before_flow - melt

    flow_fraction = np.clip(GLACIER_FLOW_RATE_PER_MYR * slope_to_target * years_myr, 0.0, GLACIER_MAX_FLOW_FRACTION)
    has_target = flow_target >= 0
    outflow = np.where(has_target, depth_after_melt * flow_fraction, 0.0)
    remaining = depth_after_melt - outflow

    inflow = np.zeros_like(elevation)
    np.add.at(inflow, flow_target[has_target], outflow[has_target])

    # Ice reaching (or accumulating on) an ocean node is discarded, not piled up -- real sea
    # ice is a different, thinner, seasonal phenomenon this model doesn't represent; this
    # same guard reads as calving where a glacier reaches a coast.
    new_depth = np.where(is_ocean, 0.0, np.clip(remaining + inflow, 0.0, None))
    return new_depth, melt


def compute_hydrology(world: "World", precipitation_at_nodes: np.ndarray, temperature_at_nodes: np.ndarray, years: float) -> HydrologyFields:
    """Runs the full flow-routing pipeline against the world's current node cloud and this
    step's climate: basin-spill -> (lake freeze, if cold enough) -> flow direction ->
    glacier accumulation/melt/flow -> precipitation(-minus-frozen)-plus-meltwater-plus-lake-
    breach downstream accumulation -> lake growth/evaporation/merge/split (see `lakes.py`) ->
    river classification (excluding flooded lake nodes, so a river ends at a lake's shore
    rather than being treated as continuing through it). A lake that's actively spilling this
    step also adds `LAKE_BREACH_EROSION_COEFFICIENT`-scaled extra "water" into the routed total
    at its own sink, sized by its connected lake's own node count (see that constant's own
    comment) -- so `flow_accum`, and everything erosion.py derives from it (river erosion,
    channel_depth/channel_width growth), reflects a spilling lake's real erosive force, not
    just the ordinary precip passing through that one node. `precipitation_at_nodes`/
    `temperature_at_nodes` are precomputed by the caller (erosion.py, which already needs the
    same climate-grid lookups for its own erosion terms) rather than looked up again here.

    A node colder than GLACIER_ACCUMULATION_TEMP_C never holds liquid water at all: any
    precipitation there is treated as fully frozen (no partial liquid/frozen split --
    matching erosion.py's existing "use precipitation is enough" simplification, just gated
    by the same accumulation threshold),
    and an existing lake sitting there freezes solid into glacier_depth this same step,
    *before* flow_target is (re)computed -- so a lake that just froze is correctly treated as
    a genuine sink again this step (see `lakes.step_lakes`'s own docstring for the freeze
    handling this feeds into)."""
    points, elevation, prev_lake_depth, prev_glacier_depth, prev_channel_depth, prev_silt_depth, is_ocean, line_refs = _gather_nodes(world)
    n = len(points)
    if n <= FLOW_NEIGHBOR_COUNT:
        empty_i = np.zeros(n, dtype=np.int64)
        empty_f = np.zeros(n)
        return HydrologyFields(
            points, elevation, is_ocean, np.zeros((n, 0), dtype=np.int64), empty_i, empty_f, empty_f, empty_f, empty_i,
            np.zeros(n, dtype=bool), empty_f, empty_f, line_refs, silt_depth=prev_silt_depth,
        )

    neighbor_idx = _build_neighbor_graph(points)
    filled_elevation, spill_target = _compute_basin_spill(elevation, is_ocean, neighbor_idx)

    is_accumulating = temperature_at_nodes < GLACIER_ACCUMULATION_TEMP_C
    freezing_lake = is_accumulating & (prev_lake_depth > 0.0)
    lake_depth_adjusted = np.where(freezing_lake, 0.0, prev_lake_depth)
    frozen_from_lake = np.where(freezing_lake, prev_lake_depth, 0.0)

    flow_target, should_spill = _compute_flow_direction(
        elevation, is_ocean, neighbor_idx, lake_depth_adjusted, filled_elevation, spill_target, prev_channel_depth
    )

    frozen_precip = np.where(is_accumulating, precipitation_at_nodes, 0.0)
    liquid_precip = np.where(is_accumulating, 0.0, precipitation_at_nodes)

    slope_to_target = _slope_to_flow_target(points, elevation, flow_target)
    new_glacier_depth, melt = _update_glaciers(
        elevation, is_ocean, flow_target, slope_to_target, prev_glacier_depth, frozen_precip, frozen_from_lake, is_accumulating, temperature_at_nodes, years
    )

    # A spilling lake's own surface area feeds extra erosive "water" in at its sink, on top of
    # whatever ordinary precip/melt landed there this step -- see LAKE_BREACH_EROSION_
    # COEFFICIENT's own comment for why. Sized from *last* step's lake extent (lake_depth_
    # adjusted), the same one-step-lagged memory should_spill itself is already computed from.
    lake_component_size = _lake_component_sizes(lake_depth_adjusted > LAKE_MIN_VISIBLE_DEPTH_M, neighbor_idx)
    lake_breach_source = np.where(should_spill, LAKE_BREACH_EROSION_COEFFICIENT * lake_component_size, 0.0)

    water_source = liquid_precip + melt + lake_breach_source
    flow_accum, water_deposited = route_downstream(elevation, is_ocean, flow_target, water_source)

    # is_river isn't known yet -- it needs this step's *final* lake_depth (below) to exclude
    # flooded nodes, so fields is built with a placeholder here and finished after.
    fields = HydrologyFields(
        points, elevation, is_ocean, neighbor_idx, flow_target, flow_accum, water_deposited, filled_elevation, spill_target,
        np.zeros(n, dtype=bool), lake_depth_adjusted, new_glacier_depth, line_refs,
    )
    fields.lake_depth, fields.silt_depth, _lake_forest, fields.lake_events = lakes.step_lakes(
        elevation, is_ocean, neighbor_idx, lake_depth_adjusted, prev_silt_depth, water_deposited, years, is_accumulating
    )

    # A lake, however wide it's flooded, is never itself "a river" -- excluding it (not just
    # its own former sink node, every flooded node) is what makes a river actually end at a
    # lake's shore instead of its classification jumping straight across to whatever's on
    # the far side of the water.
    land = ~is_ocean
    is_lake = fields.lake_depth > LAKE_MIN_VISIBLE_DEPTH_M
    is_river = np.zeros(n, dtype=bool)
    if np.any(land) and np.any(flow_accum[land] > 0):
        threshold = np.percentile(flow_accum[land], RIVER_FLOW_PERCENTILE)
        is_river = land & (flow_accum > 0) & (flow_accum >= threshold) & ~is_lake
    fields.is_river = is_river
    return fields


@dataclass
class RiverInfo:
    """One connected drainage network within `HydrologyFields.is_river` -- the River
    Inspector's unit of selection. There's no such grouping concept anywhere else in this
    module; `is_river` is otherwise just a flat per-node boolean. `member_idx` indexes back
    into the *same* HydrologyFields arrays (points/elevation/flow_accum/flow_target/...) this
    was built from, so a caller (main.py) resolves geometry/segments straight from those
    rather than this dataclass duplicating them."""

    member_idx: np.ndarray  # (M,) indices into the source HydrologyFields arrays
    mouth_idx: int  # index into the same arrays -- the network's max-flow_accum member
    mouth_type: str  # "ocean" | "lake" | "other" (a dry interior sink)
    flow_rate: float  # flow_accum at the mouth
    speed: float  # compute_river_speed at the mouth
    num_tributaries: int  # separate headwater branches feeding this network, see group_rivers


def group_rivers(fields: HydrologyFields) -> list[RiverInfo]:
    """Groups the flat `is_river` node mask into distinct connected drainage networks, via
    union-find over `flow_target` edges restricted to nodes that are `is_river` on *both*
    ends. This is exact (not a heuristic) because flow_accum is monotonically non-decreasing
    downhill: once a node clears
    RIVER_FLOW_PERCENTILE, every node further downstream in its own chain does too, so two
    river nodes joined by a flow_target edge always belong to the same real drainage network,
    and two nodes in different connected components never do.

    A network's **mouth** is always its own max-flow_accum member -- flow only accumulates
    downhill, so nothing in the group can out-flow it -- which also guarantees the mouth's
    own `flow_target` (if any) points *outside* the group: a same-or-higher-flow_accum land
    node downstream would itself have cleared the threshold and been unioned in already, so
    the only things left outside are open ocean or a true dead-end sink.

    **mouth_type** checks `lake_depth` first (a river can end at a still-spilling/draining
    lake, which reads more usefully as `"lake"` than `"ocean"` even though it also has a
    `flow_target`), then whether that `flow_target` lands on the ocean, else `"other"` (a dry
    interior sink -- `flow_target == -1` with no lake standing there).

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
        if target in river_set:
            ra, rb = find(i), find(target)
            if ra != rb:
                parent[ra] = rb

    groups: dict[int, list[int]] = {}
    for i in river_list:
        groups.setdefault(find(i), []).append(i)

    slope_to_target = _slope_to_flow_target(fields.points, fields.elevation, fields.flow_target)

    rivers: list[RiverInfo] = []
    for members in groups.values():
        member_idx = np.array(members, dtype=np.int64)
        mouth_idx = int(member_idx[np.argmax(fields.flow_accum[member_idx])])

        if fields.lake_depth[mouth_idx] > LAKE_MIN_VISIBLE_DEPTH_M:
            mouth_type = "lake"
        elif fields.flow_target[mouth_idx] >= 0:
            mouth_type = "ocean"
        else:
            mouth_type = "other"

        member_set = set(members)
        in_degree = {m: 0 for m in members}
        for m in members:
            t = flow_target_list[m]
            if t in member_set:
                in_degree[t] += 1
        headwater_count = sum(1 for m in members if in_degree[m] == 0)

        rivers.append(
            RiverInfo(
                member_idx=member_idx,
                mouth_idx=mouth_idx,
                mouth_type=mouth_type,
                flow_rate=float(fields.flow_accum[mouth_idx]),
                speed=float(compute_river_speed(slope_to_target[mouth_idx], fields.flow_accum[mouth_idx])),
                num_tributaries=max(0, headwater_count - 1),
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
