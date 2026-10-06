"""Least-cost depression breaching over the hydrology graph (issue #297, part 1).

A coarse terrain cell's single elevation is its centre sample, not the lowest ground it
contains. A shallow pit between coarse samples is therefore often an artefact. In the real
landscape a river would cut through the low barrier long before a lake filled to the
neighbouring cell's centre height. This module decides, once per step and for every node at
once, which pits get cut through and which stay closed:

1. **Passage elevation.** Water crosses a cell at its *passage* elevation: the centre
   elevation, lowered by its established channel (`channel_depth`) and by any breach notch
   cut on an earlier step (`breach_notch_depth_m`). The notch floor never goes below the
   cell's lowest neighbour's centre, because a corridor through a cell only has to connect
   that cell's neighbours.
2. **Least climb to the ocean.** One multi-source Dijkstra from every connected-ocean node,
   over the hydrology k-NN graph, gives every node the least total height its water would have
   to climb to reach the ocean. A step from `i` onto `j` costs the time to cut `j` down by
   that climb, `max(0, passage[j] - passage[i])`, expressed as metres of reference-strength
   rock. The cut goes through `j`'s loose mobile cover first, at the weak-rock rate, then
   through bedrock at `j`'s own rate. A node on an ordinary downhill path costs zero.
3. **Breach or fill.** A pit (a land node with no strictly lower neighbour) is breached only
   if two things hold. First, its cost fits this step's carving budget,
   `BREACH_REFERENCE_CARVE_M_PER_MYR * dt`. Second, every node on its least-cost path can be
   cut down to the water level arriving from the pit within that node's own carving limit
   for the step: its cover at the weak-rock rate, then bedrock at its own rate for the time
   left (`_carve_limit_m`). The second check matters because the cost counts only climbs, while
   draining also means cutting nodes past the rim that still stand above the pit floor, and
   those can be in stronger rock. A breached pit's path is notched until it never rises, so
   the pit always drains. Each reach of the river carves at the same time, so each node's
   limit applies to its own cut, not to the sum.
4. **Only where water overflows.** A river carves an outlet only if the basin overflows. A
   pit is breached only if its catchment's runoff exceeds what a lake filled to its rim would
   evaporate. Otherwise the lake would level off below the rim, so the pit stays a real
   endorheic basin, like the Great Basin, the Caspian or the Dead Sea. A dry pit always fails
   this test. Runoff is precipitation minus evapotranspiration, from Fu's (1981) form of the
   Budyko curve. Open-water evaporation rises with temperature.

Every pit not breached stays closed. It is left to `lakes.py`, which fills it, spills it, or
keeps it as an endorheic lake, depending on its water balance.

This is Lindsay's (2016) breach-first depression treatment with a physical cost cap. The cap
uses the carving rate of a river: 10-100 m per 100 kyr, depending on rock strength.

The notches are returned per node. `hydrology.compute_hydrology` passes the notched passages to
`lakes.build_lake_hierarchy` as interface passes. `erosion.py` keeps them in their own
persisted field, `breach_notch_depth_m`, which records sub-cell relief, not volume removed
from the cell mean. Only the passage elevation reads that field. `channel_depth`, which drives
erosion's channel boost, channel-preferring flow, river evaporation and rendering, keeps
recording only rock that was actually carved. Notches are written only for pits that fully
drain, so a breach that needs more than one step's carving never accumulates on its own. It
opens only once the lake fills and spills, and the existing breach-erosion term cuts its rim.
Later phases of #297 will carry the notch volume in the sediment ledger.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

# River incision rates into bedrock: roughly 10-100 m per 100 kyr, from strong to weak rock.
# The reference rate is ordinary continental rock. Cratons are strong, and loose mobile cover
# (sediment, lake silt, regolith: erosion.py's `mobile_cover_m`) is weak. A step of `dt` Myr
# can cut rate * dt metres, so longer steps cut deeper notches.
BREACH_REFERENCE_CARVE_M_PER_MYR = 300.0
BREACH_STRONG_CARVE_M_PER_MYR = 100.0
BREACH_WEAK_CARVE_M_PER_MYR = 1000.0
# Same ceiling as erosion.MAX_CHANNEL_DEPTH_M, which caps every notch persisted into
# channel_depth. It's repeated here because erosion.py imports hydrology.py, which imports this.
BREACH_MAX_NOTCH_M = 2000.0
# Water-balance gate (step 4 of the module docstring). Open-water (lake) evaporation is
# roughly linear in mean annual temperature: about 0.3 m/yr at 0 C, 1.2 m/yr at 15 C, and
# 1.8 m/yr at 25 C. Those are typical of the Caspian, the Great Basin and the Dead Sea.
OPEN_WATER_EVAPORATION_MM_AT_0C = 300.0
OPEN_WATER_EVAPORATION_MM_PER_C = 60.0
OPEN_WATER_EVAPORATION_MIN_MM = 50.0
# Fu's Budyko-curve shape parameter. 2.6 is the commonly fitted global mean.
BUDYKO_FU_OMEGA = 2.6
# How much more water a pit's catchment must deliver than its rim-level lake would evaporate
# before the pit is breached. Raise it to keep more basins closed.
ENDORHEIC_DEMAND_FACTOR = 1.0
# Added to every graph step so that Dijkstra keeps zero-climb edges, and so that ties go to the
# path with fewer hops. It also makes each node's cost strictly larger than its next hop's,
# which `breach_depressions` relies on to visit paths in order.
_HOP_EPSILON_M = 1e-6
# A pit needs breaching only if its water climbs at least this much. Below it, the cost is just
# the accumulated `_HOP_EPSILON_M`, as on a pit that an earlier step's notch already drains.
_MIN_BREACH_CLIMB_M = 0.01


@dataclass
class BreachResult:
    """Per-node output of `breach_depressions`. All arrays have shape (N,)."""

    # Least total climb to the ocean, in metres of reference-strength rock. 0 at the ocean,
    # inf where no path exists.
    cost_m: np.ndarray
    # The neighbour that the least-cost path visits next, or -1 at the ocean or where no
    # path exists.
    next_hop: np.ndarray
    # The passage elevation before this step's invented notches (centre minus channel notch).
    channel_passage_m: np.ndarray
    # Notch invented this step, in metres (>= 0).
    notch_m: np.ndarray
    # `channel_passage_m - notch_m`: the passage elevation the lake hierarchy should use.
    passage_m: np.ndarray
    # Pit nodes that were breached this step.
    breached_pits: np.ndarray
    # Pit nodes whose cost was over budget, which stay closed and are left to lakes.py.
    closed_pits: np.ndarray
    # Pit nodes cheap enough to breach whose water balance keeps them closed (endorheic).
    # These stay closed too. Empty when no climate is given.
    endorheic_pits: np.ndarray = None

    def __post_init__(self) -> None:
        if self.endorheic_pits is None:
            self.endorheic_pits = np.zeros(0, dtype=np.int64)


@dataclass
class WaterBalance:
    """Per-node climate inputs for the breach gate, all with shape (N,).

    `precipitation_mm` should be liquid precipitation, as the lake balance uses.
    `loss_fraction` is the in-transit river evaporation each node takes from water passing
    through it (`hydrology.river_transit_loss_fraction`, independent of step length). The gate routes runoff to each pit
    along steepest descent, losing that fraction at every node on the way, the same rule
    `hydrology.route_downstream` applies when lakes are balanced. Glacier melt isn't included:
    it's computed after breaching, from the ice routing that breaching feeds."""

    precipitation_mm: np.ndarray
    temperature_c: np.ndarray
    area_m2: np.ndarray
    loss_fraction: np.ndarray | None = None


def open_water_evaporation_mm(temperature_c: np.ndarray) -> np.ndarray:
    """Annual evaporation from a lake surface (mm/yr), linear in mean temperature."""
    rate = OPEN_WATER_EVAPORATION_MM_AT_0C + OPEN_WATER_EVAPORATION_MM_PER_C * np.asarray(temperature_c, dtype=float)
    return np.maximum(rate, OPEN_WATER_EVAPORATION_MIN_MM)


def runoff_mm(precipitation_mm: np.ndarray, potential_evaporation_mm: np.ndarray) -> np.ndarray:
    """Annual runoff (mm/yr) from Fu's Budyko curve. Evapotranspiration is
    `P * (1 + phi - (1 + phi**w) ** (1/w))` with aridity `phi = PET / P`, so runoff tends to
    `P` in a cold, wet climate and to nothing in a hot, dry one."""
    p = np.clip(np.asarray(precipitation_mm, dtype=float), 0.0, None)
    pet = np.asarray(potential_evaporation_mm, dtype=float)
    phi = np.divide(pet, p, out=np.full_like(p, np.inf), where=p > 0.0)
    w = BUDYKO_FU_OMEGA
    with np.errstate(over="ignore", invalid="ignore"):
        et_ratio = 1.0 + phi - np.power(1.0 + np.power(phi, w), 1.0 / w)
    et_ratio = np.where(np.isfinite(et_ratio), np.clip(et_ratio, 0.0, 1.0), 1.0)
    return p * (1.0 - et_ratio)


def _routed_to_sinks(
    elevation: np.ndarray, neighbor_idx: np.ndarray, is_ocean: np.ndarray, keep: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """`(sink, carried)` for each node. `sink` is the node its steepest descent ends at: a
    pit, or an ocean node. `carried` is the share of water starting at the node that reaches
    the sink, after every node on the way, the sink included, keeps only `keep` of it. Uses
    vectorised pointer doubling, so it takes about log2(longest path) passes."""
    n = len(elevation)
    neighbor_elev = elevation[neighbor_idx]
    col = np.argmin(neighbor_elev, axis=1)
    lowest = neighbor_idx[np.arange(n), col]
    receiver = np.where((neighbor_elev[np.arange(n), col] < elevation) & ~is_ocean, lowest, np.arange(n))
    is_root = receiver == np.arange(n)
    # `share[i]` is the product of `keep` from i up to, not including, `receiver[i]`.
    share = np.where(is_root, 1.0, keep)
    while not np.array_equal(receiver[receiver], receiver):
        share = share * share[receiver]
        receiver = receiver[receiver]
    return receiver, share * keep[receiver]


def bedrock_carve_rate_m_per_myr(craton_strength: np.ndarray) -> np.ndarray:
    """Each node's bedrock river-incision rate, interpolated in log space from the reference
    rate (ordinary rock) toward the strong-rock rate by craton strength (0..1,
    `cratons.strength`). Loose mobile cover on top carves at BREACH_WEAK_CARVE_M_PER_MYR; see
    `_cut_time_myr`."""
    strong = np.clip(np.asarray(craton_strength, dtype=float), 0.0, 1.0)
    return np.exp(
        np.log(BREACH_REFERENCE_CARVE_M_PER_MYR)
        + strong * np.log(BREACH_STRONG_CARVE_M_PER_MYR / BREACH_REFERENCE_CARVE_M_PER_MYR)
    )


def _cut_time_myr(depth_m: np.ndarray, cover_m: np.ndarray, bedrock_rate: np.ndarray) -> np.ndarray:
    """Time to cut `depth_m` down through a cell: its loose cover first, at the weak-rock
    rate, then bedrock at `bedrock_rate`. Only the share of the cut that's actually in cover
    gets the weak rate."""
    depth_m = np.clip(depth_m, 0.0, None)
    in_cover = np.minimum(depth_m, cover_m)
    return in_cover / BREACH_WEAK_CARVE_M_PER_MYR + (depth_m - in_cover) / bedrock_rate


def _carve_limit_m(cover_m: np.ndarray, bedrock_rate: np.ndarray, dt_myr: float) -> np.ndarray:
    """How deep a river can cut each cell in `dt_myr`: through its cover at the weak-rock rate,
    then into bedrock at its own rate for the time left. Capped at BREACH_MAX_NOTCH_M."""
    cover_time = cover_m / BREACH_WEAK_CARVE_M_PER_MYR
    limit = np.where(
        cover_time >= dt_myr,
        BREACH_WEAK_CARVE_M_PER_MYR * dt_myr,
        cover_m + (dt_myr - cover_time) * bedrock_rate,
    )
    return np.minimum(limit, BREACH_MAX_NOTCH_M)


def channel_passage_elevation(elevation: np.ndarray, channel_depth: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Centre elevation, lowered by `channel_depth`: the cell's established channel plus any
    persisted breach notch. The result is never below the lowest neighbour's centre and never
    above the cell's own centre."""
    elevation = np.asarray(elevation, dtype=float)
    notched = elevation - np.clip(np.asarray(channel_depth, dtype=float), 0.0, None)
    if neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        return elevation.copy()
    floor = elevation[neighbor_idx].min(axis=1)
    return np.minimum(elevation, np.maximum(notched, floor))


def least_climb_to_ocean(
    passage: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    bedrock_rate: np.ndarray,
    cover_m: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """`(cost_m, next_hop)` for every node, from one multi-source Dijkstra rooted at every
    ocean node. Water moves from `i` to each `j` in `neighbor_idx[i]`. That step costs the time
    to cut `j` down by `max(0, passage[j] - passage[i])` (`_cut_time_myr`), times the reference
    rate, so it's in metres of reference-strength rock. `_HOP_EPSILON_M` is added to every step. Dijkstra runs
    over the transposed graph, so a distance *from* the ocean there is a distance *to* the
    ocean here, and its predecessor is the next hop."""
    n = len(passage)
    cost = np.full(n, np.inf)
    next_hop = np.full(n, -1, dtype=np.int64)
    sources = np.flatnonzero(is_ocean)
    if n == 0 or len(sources) == 0 or neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        cost[sources] = 0.0
        return cost, next_hop
    k = neighbor_idx.shape[1]
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.reshape(-1)
    cover_m = np.zeros(n) if cover_m is None else np.clip(np.asarray(cover_m, dtype=float), 0.0, None)
    bedrock_rate = np.asarray(bedrock_rate, dtype=float)
    climb = np.clip(passage[cols] - passage[rows], 0.0, None)
    weight = _cut_time_myr(climb, cover_m[cols], bedrock_rate[cols]) * BREACH_REFERENCE_CARVE_M_PER_MYR + _HOP_EPSILON_M
    # Water never leaves the ocean, so no edge starts there. This also stops a path from
    # running back out of the ocean.
    keep = ~is_ocean[rows]
    rows, cols, weight = rows[keep], cols[keep], weight[keep]
    # A padded neighbour row can repeat an edge, and a sparse matrix sums duplicates, so keep
    # one copy of each.
    _, first = np.unique(rows * n + cols, return_index=True)
    transposed = csr_matrix((weight[first], (cols[first], rows[first])), shape=(n, n))
    dist, pred, _ = dijkstra(transposed, directed=True, indices=sources, min_only=True, return_predecessors=True)
    cost[:] = dist
    cost[sources] = 0.0
    next_hop[:] = np.where(pred >= 0, pred, -1)
    next_hop[sources] = -1
    return cost, next_hop


def _pits(elevation: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Land nodes with no strictly lower neighbour. This is the same sink test as
    `hydrology._compute_flow_direction` uses, and every lake catchment drains to one."""
    has_lower = np.any(elevation[neighbor_idx] < elevation[:, None], axis=1)
    return np.flatnonzero(~has_lower & ~is_ocean)


def breach_depressions(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    channel_depth: np.ndarray,
    carve_rate: np.ndarray,
    years: float,
    water: WaterBalance | None = None,
    prior_notch_m: np.ndarray | None = None,
    mobile_cover_m: np.ndarray | None = None,
) -> BreachResult:
    """Breach every pit whose least climb to the ocean fits this step's carving budget. See
    the module docstring.

    The result doesn't depend on node order. Each node's notch is the deepest cut that any
    breached path through it needs. That cut is set by the lowest water level arriving at the
    node from upstream. Nodes are visited from the highest cost down, which is a valid
    upstream-to-downstream order because each node's cost is strictly larger than its next
    hop's (see `_HOP_EPSILON_M`). One pass visits each node on the breached paths once. It is
    not a separate path search per pit.

    `carve_rate` is each node's bedrock incision rate (`bedrock_carve_rate_m_per_myr`), and
    `mobile_cover_m` the loose cover on top of it, which carves at the weak-rock rate.

    `prior_notch_m` is the breach notch persisted from earlier steps. It lowers the passage
    along with `channel_depth`. `notch_m` in the result is only this step's new cut.

    The drain check (step 3 of the module docstring) walks each candidate pit's least-cost
    path once along `next_hop`. That costs the path's length per candidate, not a search.

    With `water`, pits that pass the cost and drain tests also have to pass the water-balance
    test (step 4 of the module docstring). The rim level is the highest passage on the pit's
    least-cost path, found in the same walk as the drain check, in metres, unlike the
    rock-weighted cost. The rim-level lake is every catchment member below it. That path
    minimises total climb, not the highest point, so its top can sit above the basin's real
    spill point but never below it. The estimate therefore errs toward keeping a basin closed."""
    elevation = np.asarray(elevation, dtype=float)
    n = len(elevation)
    dt_myr = years / 1_000_000.0
    depth = np.asarray(channel_depth, dtype=float)
    if prior_notch_m is not None:
        depth = depth + np.clip(np.asarray(prior_notch_m, dtype=float), 0.0, None)
    channel_passage = channel_passage_elevation(elevation, depth, neighbor_idx)
    bedrock_rate = np.asarray(carve_rate, dtype=float)
    cover = np.zeros(n) if mobile_cover_m is None else np.clip(np.asarray(mobile_cover_m, dtype=float), 0.0, None)
    cost, next_hop = least_climb_to_ocean(channel_passage, is_ocean, neighbor_idx, bedrock_rate, cover)
    notch = np.zeros(n)
    empty = np.zeros(0, dtype=np.int64)
    if n == 0 or neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        return BreachResult(cost, next_hop, channel_passage, notch, channel_passage.copy(), empty, empty)

    budget = min(BREACH_REFERENCE_CARVE_M_PER_MYR * dt_myr, BREACH_MAX_NOTCH_M)
    pits = _pits(elevation, is_ocean, neighbor_idx)
    pit_cost = cost[pits]
    breachable = (pit_cost > _MIN_BREACH_CLIMB_M) & (pit_cost <= budget)
    closed_pits = pits[np.isfinite(pit_cost) & (pit_cost > budget)]
    endorheic_pits = np.zeros(0, dtype=np.int64)
    max_cut = _carve_limit_m(cover, bedrock_rate, dt_myr)
    if np.any(breachable):
        candidates = np.flatnonzero(breachable)
        drains, rim_level = _walk_paths(pits[candidates], channel_passage, next_hop, max_cut)
        closed_pits = np.sort(np.concatenate([closed_pits, pits[candidates[~drains]]]))
        breachable[candidates[~drains]] = False
        if water is not None and np.any(drains):
            draining = candidates[drains]
            overflows = _overflows(elevation, is_ocean, neighbor_idx, pits[draining], rim_level[drains], water)
            endorheic_pits = pits[draining[~overflows]]
            breachable[draining[~overflows]] = False
    breached_pits = pits[breachable]

    if len(breached_pits):
        # Mark every node on a breached path with vectorised pointer-chasing. Each round
        # advances one hop on every path that isn't done yet.
        on_path = np.zeros(n, dtype=bool)
        frontier = breached_pits
        while len(frontier):
            frontier = frontier[~on_path[frontier]]
            on_path[frontier] = True
            frontier = next_hop[frontier]
            frontier = np.unique(frontier[frontier >= 0])
        path_nodes = np.flatnonzero(on_path & ~is_ocean)
        order = path_nodes[np.argsort(-cost[path_nodes], kind="stable")]

        # The lowest water level arriving at each node from a breached pit upstream. The
        # node's cut takes its passage down to that level. `_paths_drain` already checked
        # that every breached pit's own cut fits each node's carving limit, and the combined
        # cut is the largest of those, so it fits too.
        arriving = np.full(n, np.inf)
        arriving[breached_pits] = channel_passage[breached_pits]
        arriving_list = arriving.tolist()
        passage_list = channel_passage.tolist()
        next_list = next_hop.tolist()
        notch_list = notch.tolist()
        for i in order.tolist():
            level = arriving_list[i]
            cut = max(passage_list[i] - level, 0.0)
            notch_list[i] = cut
            out_level = passage_list[i] - cut
            nxt = next_list[i]
            if nxt >= 0 and out_level < arriving_list[nxt]:
                arriving_list[nxt] = out_level
        notch = np.asarray(notch_list)

    return BreachResult(cost, next_hop, channel_passage, notch, channel_passage - notch, breached_pits, closed_pits, endorheic_pits)


def _walk_paths(
    pits: np.ndarray, passage: np.ndarray, next_hop: np.ndarray, max_cut: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """`(drains, rim_level)` for each pit, from one walk along its least-cost path.
    `drains`: whether the path can be cut until it never rises, with every node cut no deeper
    than its own `max_cut`. The water level starts at the pit's passage and only falls along
    the path. `rim_level`: the highest passage on the path, in metres."""
    passage_list = passage.tolist()
    next_list = next_hop.tolist()
    max_cut_list = max_cut.tolist()
    drains = []
    rims = []
    for pit in pits.tolist():
        level = passage_list[pit]
        rim = level
        node = next_list[pit]
        ok = True
        while node >= 0:
            here = passage_list[node]
            rim = max(rim, here)
            if here - level > max_cut_list[node] + 1e-9:
                ok = False
            level = min(level, here)
            node = next_list[node]
        drains.append(ok)
        rims.append(rim)
    return np.asarray(drains, dtype=bool), np.asarray(rims, dtype=float)


def _overflows(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    pits: np.ndarray,
    rim_level: np.ndarray,
    water: WaterBalance,
) -> np.ndarray:
    """Whether each pit's catchment runoff exceeds the evaporation of a lake filled to
    `rim_level` (metres), scaled by ENDORHEIC_DEMAND_FACTOR. Both sides are in m^3/yr."""
    n = len(elevation)
    loss = np.zeros(n) if water.loss_fraction is None else np.asarray(water.loss_fraction, dtype=float)
    sink, carried = _routed_to_sinks(elevation, neighbor_idx, is_ocean, 1.0 - loss)
    area = np.asarray(water.area_m2, dtype=float)
    evaporation_mm = open_water_evaporation_mm(water.temperature_c)
    runoff_m3 = runoff_mm(water.precipitation_mm, evaporation_mm) / 1000.0 * area
    inflow = np.bincount(sink, weights=runoff_m3 * carried, minlength=n)
    level = np.full(n, -np.inf)
    level[pits] = rim_level
    flooded = elevation < level[sink]
    demand = np.bincount(sink[flooded], weights=(evaporation_mm / 1000.0 * area)[flooded], minlength=n)
    return inflow[pits] > ENDORHEIC_DEMAND_FACTOR * demand[pits]


def interface_pass_elevation(passage: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Pass elevations aligned with `neighbor_idx`, in the form
    `lakes.build_lake_hierarchy(interface_pass_elevation=...)` takes. Water crossing from `i`
    to `j` must clear the higher of the two passages."""
    return np.maximum(passage[:, None], passage[neighbor_idx])
