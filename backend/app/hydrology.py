"""Rivers and lakes: flow routing, basin/lake detection, and downstream flow accumulation
over the world's current node cloud.

Ported from plate-sim's hydrology.py, adapted from its fixed grid (whose regular 8-neighbor
adjacency gives flow routing, depression filling, and downstream accumulation a natural
substrate) to mantle-bloom's irregular per-plate node cloud. All three of plate-sim's core
algorithms -- steepest-descent flow direction, priority-flood basin-spill, and
elevation-ordered downstream accumulation -- turn out not to actually need a *grid*, only a
*graph*: this module builds one via a whole-world k-nearest-neighbor query (the same
technique reassign.py/erosion.py already use for their own whole-world passes), then runs
the same three algorithms directly on it.

**Persistence.** Unlike plate-sim -- whose plates move relative to a fixed grid, so a
persistent field like channel_depth needs deliberate semi-Lagrangian advection every step to
keep following the crust -- mantle-bloom's elevation-line nodes already rotate exactly with
their own plate. channel_depth and lake_depth, stored as ordinary parallel arrays on
ElevationLine right alongside elevation itself (see plates.py), get that same "just works"
persistence for free: no advection scheme needed, since rotating a plate never touches those
arrays at all. flow_target/flow_accum/river_speed are deliberately *not* persisted, mirroring
plate-sim (where they're recomputed fresh every step too, from that step's real climate) --
they're purely this-step derived quantities, cached on World (see World.hydrology_cache) only
so a later same-turn caller (rendering, stats) doesn't recompute them again.

**Deliberate deviation from plate-sim**: plate-sim computes flow routing *twice* per step
(once in hydrology.py for the real river_flow/rendering fields, again inside erosion.py for
the water_accum erosion itself needs) -- an accepted redundancy there, cheap under numba JIT.
Flow routing here has no JIT and is comparatively expensive, so this module computes it once
and erosion.py reuses the same result, rather than duplicating a cost that's real here.

**Lakes are single nodes**, same simplification plate-sim itself documents: a lake is its own
sink node plus a `lake_depth`, not a flood-filled multi-node shoreline region.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from .plates import ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

# Matches plate-sim's 8-neighbor D8 structure -- the graph every algorithm below runs on.
FLOW_NEIGHBOR_COUNT = 8

# Top decile of land flow_accum counts as "a river" -- same threshold plate-sim's own
# _major_river_mask uses for both rendering and the major_river_fraction stat.
RIVER_FLOW_PERCENTILE = 90.0

# Lake growth/evaporation -- meaning ported from plate-sim, but NOT its evaporation rate
# values verbatim: plate-sim takes much finer internal substeps than mantle-bloom's typical
# 1-10 Myr step, so its own per-Myr exponential-decay rate, applied over a full mantle-bloom
# step, evaporated a lake almost entirely in a single step (confirmed directly -- a lake at
# 15m with no new inflow dropped to exactly 0 after one 1 Myr step). Rescaled down to be
# consistent with this codebase's other already-tuned relaxation rates at its own actual
# step granularity (boundary.DIVERGENT_RELAX_RATE_PER_MYR = 0.5,
# bathymetry.BATHYMETRY_RELAX_RATE_PER_MYR = 0.3). MIN/MAX evaporation temperature scaling
# is dropped entirely: mantle-bloom's erosion.py already decided "use precipitation is
# enough," no temperature-gated liquid/snow split, so lake evaporation follows the same
# simplification rather than reading temperature at all.
LAKE_MIN_VISIBLE_DEPTH_M = 1.0
LAKE_FILL_RATE = 0.12
LAKE_EVAPORATION_RATE_PER_MYR = 0.3
LAKE_EVAPORATION_BASELINE_M_PER_MYR = 0.5


@dataclass
class HydrologyFields:
    """Everything derived from one whole-world flow-routing pass, all shape (N,) aligned
    with `points`/`elevation` -- the irregular-node-cloud analogue of plate-sim's
    HydrologyField grid. `line_refs` is (plate, line_index, start, end) per line, letting a
    caller slice any of these flat arrays back onto a specific line's own node range."""

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
    line_refs: list[tuple[Plate, int, int, int]]


def _gather_nodes(world: "World") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position, elevation, prior lake_depth, and whether it's ocean
    (elevation <= 0, the sea-level convention used everywhere else in this codebase),
    concatenated, alongside (plate, line_index, start, end) references -- same shape as
    erosion.py's/bathymetry.py's own _gather_nodes."""
    points_list, elev_list, lake_list = [], [], []
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
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        empty = np.zeros(0)
        return np.zeros((0, 3)), empty, empty, np.zeros(0, dtype=bool), []
    points = np.concatenate(points_list, axis=0)
    elevation = np.concatenate(elev_list, axis=0)
    prev_lake_depth = np.concatenate(lake_list, axis=0)
    is_ocean = elevation <= 0.0
    return points, elevation, prev_lake_depth, is_ocean, line_refs


def _build_neighbor_graph(points: np.ndarray) -> np.ndarray:
    n = len(points)
    k = min(FLOW_NEIGHBOR_COUNT, n - 1)
    tree = cKDTree(points)
    _, neighbor_idx = tree.query(points, k=k + 1)
    return neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0


def _compute_basin_spill(elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Priority-flood depression filling: a multi-source Dijkstra seeded from every ocean
    node, relaxing by *max* (the highest point a path is forced to cross) rather than by
    sum -- a minimax path cost, not a shortest path. Direct port of plate-sim's own
    compute_basin_spill, generalized from grid 8-neighbor adjacency to the k-NN graph's
    edges. Returns (filled_elevation, spill_target): the minimal-bottleneck elevation to
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


def _compute_flow_direction(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    prev_lake_depth: np.ndarray,
    filled_elevation: np.ndarray,
    spill_target: np.ndarray,
) -> np.ndarray:
    """Steepest-descent flow target per node: the lowest of its k nearest neighbors,
    strictly below its own elevation, or -1 if none (a sink) -- unless the sink's current
    water surface (elevation + the *previous* step's lake_depth, the same one-step-lagged
    "memory" plate-sim's own compute_flow_direction relies on) has already reached its
    basin's true spill point, in which case it redirects to spill_target instead of staying
    a dead-end sink forever. Ocean nodes are never routed (they're a destination, not a
    source)."""
    n = len(elevation)
    neighbor_elev = elevation[neighbor_idx]  # (n, k)
    own_elev = elevation[:, None]
    lower_mask = neighbor_elev < own_elev
    masked = np.where(lower_mask, neighbor_elev, np.inf)
    best_col = np.argmin(masked, axis=1)
    rows = np.arange(n)
    has_lower = np.isfinite(masked[rows, best_col])
    flow_target = np.where(has_lower, neighbor_idx[rows, best_col], -1).astype(np.int64)

    is_sink = (flow_target < 0) & ~is_ocean
    water_surface = elevation + prev_lake_depth
    should_spill = is_sink & (water_surface >= filled_elevation)
    flow_target = np.where(should_spill, spill_target, flow_target).astype(np.int64)
    return np.where(is_ocean, -1, flow_target).astype(np.int64)


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
    in this same order -- direct port of plate-sim's own route_downstream (its loss_fraction
    parameter, used there only for in-transit river evaporation, isn't ported -- not asked
    for and this module already drops temperature-driven effects on hydrology, see
    LAKE_EVAPORATION_* above). Returns (through_flux, deposited): each node's own
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


def compute_hydrology(world: "World", precipitation_at_nodes: np.ndarray) -> HydrologyFields:
    """Runs the full flow-routing pipeline against the world's current node cloud and this
    step's climate: basin-spill -> flow direction -> precipitation-weighted downstream
    accumulation -> river classification. `precipitation_at_nodes` is precomputed by the
    caller (erosion.py, which already needs the same climate-grid lookup for its own rain
    erosion term) rather than looked up again here."""
    points, elevation, prev_lake_depth, is_ocean, line_refs = _gather_nodes(world)
    n = len(points)
    if n <= FLOW_NEIGHBOR_COUNT:
        empty_i = np.zeros(n, dtype=np.int64)
        empty_f = np.zeros(n)
        return HydrologyFields(
            points, elevation, is_ocean, np.zeros((n, 0), dtype=np.int64), empty_i, empty_f, empty_f, empty_f, empty_i, np.zeros(n, dtype=bool), line_refs
        )

    neighbor_idx = _build_neighbor_graph(points)
    filled_elevation, spill_target = _compute_basin_spill(elevation, is_ocean, neighbor_idx)
    flow_target = _compute_flow_direction(elevation, is_ocean, neighbor_idx, prev_lake_depth, filled_elevation, spill_target)
    flow_accum, water_deposited = route_downstream(elevation, is_ocean, flow_target, precipitation_at_nodes)

    land = ~is_ocean
    is_river = np.zeros(n, dtype=bool)
    if np.any(land) and np.any(flow_accum[land] > 0):
        threshold = np.percentile(flow_accum[land], RIVER_FLOW_PERCENTILE)
        is_river = land & (flow_accum > 0) & (flow_accum >= threshold)

    return HydrologyFields(
        points, elevation, is_ocean, neighbor_idx, flow_target, flow_accum, water_deposited, filled_elevation, spill_target, is_river, line_refs
    )


def update_lakes(fields: HydrologyFields, prev_lake_depth: np.ndarray, water_deposited: np.ndarray, years: float) -> np.ndarray:
    """New lake_depth per node -- grows at sink nodes fed by route_downstream's own
    terminal deposition there (`water_deposited`, the water-routing pass's `deposited`
    output -- nonzero only at a true sink, since a through-flowing node's water all
    continues downstream), evaporates elsewhere, capped at the basin's true spill depth
    (`filled_elevation - elevation`) and pinned there once reached rather than evaporating
    back down -- see this module's docstring and boundary.py's own divergent-relaxation
    style for why a persistent field needs this kind of steady-state pinning, and
    plate-sim's own docs for the oscillation bug this specifically avoids (a full lake that
    evaporates even slightly un-redirects its flow_target back to a sink, refills, and
    oscillates forever instead of settling as a real, continuously-draining lake)."""
    years_myr = years / 1_000_000.0
    is_sink = (fields.flow_target < 0) & ~fields.is_ocean
    retention = np.exp(-LAKE_EVAPORATION_RATE_PER_MYR * years_myr)
    inflow_gain = np.where(is_sink, LAKE_FILL_RATE * water_deposited * years_myr, 0.0)
    baseline_loss = LAKE_EVAPORATION_BASELINE_M_PER_MYR * years_myr
    cap = np.clip(fields.filled_elevation - fields.elevation, 0.0, None)

    grown = np.clip(prev_lake_depth * retention + inflow_gain - baseline_loss, 0.0, cap)
    not_sink_depth = np.zeros_like(prev_lake_depth)
    result = np.where(is_sink, grown, not_sink_depth)
    return np.where(fields.is_ocean, 0.0, result)


def write_lake_depth(fields: HydrologyFields, new_lake_depth: np.ndarray) -> None:
    """Writes new_lake_depth (aligned with fields.points/elevation) back onto each line's
    own lake_depth array, matching fields.line_refs' (plate, line_index, start, end)
    slices."""
    for plate, line_index, start, end in fields.line_refs:
        line = plate.lines[line_index]
        plate.lines[line_index] = ElevationLine(
            phi=line.phi,
            theta=line.theta,
            elevation=line.elevation,
            channel_depth=line.channel_depth,
            lake_depth=new_lake_depth[start:end],
        )
