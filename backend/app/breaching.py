"""Least-cost depression breaching over the hydrology graph (issue #297, part 1).

A coarse terrain cell's single elevation is its centre sample, not the lowest ground it
contains. A shallow pit between coarse samples is therefore often an artefact. In the real
landscape a river would cut through the low barrier long before a lake filled to the
neighbouring cell's centre height. This module decides, once per step and for every node at
once, which pits get cut through and which stay closed:

1. **Passage elevation.** Water crosses a cell at its *passage* elevation: the centre elevation,
   lowered by an established channel (`channel_depth`) that already notches the cell. The
   notch floor never goes below the cell's lowest neighbour's centre, because a corridor
   through a cell only has to connect that cell's neighbours.
2. **Least climb to the ocean.** One multi-source Dijkstra from every connected-ocean node,
   over the hydrology k-NN graph, gives every node the least total height its water would have
   to climb to reach the ocean. A step from `i` onto `j` costs `max(0, passage[j] -
   passage[i])`, divided by `j`'s relative erodibility, so the cost is in metres of
   reference-strength rock. A node on an ordinary downhill path costs zero.
3. **Breach or fill.** A pit (a land node with no strictly lower neighbour) whose cost fits this
   step's carving budget, `BREACH_REFERENCE_CARVE_M_PER_MYR * dt`, is breached. Every node on
   its least-cost path is notched down to the water level arriving from upstream, so the path
   becomes non-increasing and the pit drains. The cuts can add up to more than the total
   climb, because a node past the rim can still stand above the pit floor. But no single cut
   is larger than the total climb. Each reach of the river carves its own cut at the same
   time, so the budget limits each cut, not their sum. A pit that costs more stays closed. It is left
   to `lakes.py`, which fills it, spills it, or keeps it as an endorheic lake, depending on its
   water balance.

This is Lindsay's (2016) breach-first depression treatment with a physical cost cap. The cap
uses the carving rate of a river: 10-100 m per 100 kyr, depending on rock strength.

The notches are returned per node. `hydrology.compute_hydrology` passes the notched passages to
`lakes.build_lake_hierarchy` as interface passes. `erosion.py` adds the invented notches to
`channel_depth`, which records them as sub-cell relief, not volume removed from the cell mean.
A breach that needs more than one step's budget therefore never accumulates on its own. It
opens only once the lake fills and spills, and the existing breach-erosion term cuts its rim.
Later phases of #297 will carry the notch volume in the sediment ledger.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

# River incision rates into bedrock: roughly 10-100 m per 100 kyr, from strong to weak rock.
# The reference rate is ordinary continental rock. Cratons are strong, and lake silt is loose
# and weak. A step of `dt` Myr can cut rate * dt metres, so longer steps cut deeper notches.
BREACH_REFERENCE_CARVE_M_PER_MYR = 300.0
BREACH_STRONG_CARVE_M_PER_MYR = 100.0
BREACH_WEAK_CARVE_M_PER_MYR = 1000.0
# Silt this thick or more makes a cell fully weak.
BREACH_WEAK_SILT_REFERENCE_M = 10.0
# Same ceiling as erosion.MAX_CHANNEL_DEPTH_M, which caps every notch persisted into
# channel_depth. It's repeated here because erosion.py imports hydrology.py, which imports this.
BREACH_MAX_NOTCH_M = 2000.0
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


def carve_rate_m_per_myr(craton_strength: np.ndarray, silt_depth_m: np.ndarray) -> np.ndarray:
    """Each node's river-incision rate, interpolated in log space. Ordinary rock carves at the
    reference rate. Craton strength (0..1, `cratons.strength`) moves the rate toward the
    strong-rock rate. Silt cover moves the remaining range toward the weak-rock rate."""
    strong = np.clip(np.asarray(craton_strength, dtype=float), 0.0, 1.0)
    weak = np.clip(np.asarray(silt_depth_m, dtype=float) / BREACH_WEAK_SILT_REFERENCE_M, 0.0, 1.0) * (1.0 - strong)
    log_rate = (
        np.log(BREACH_REFERENCE_CARVE_M_PER_MYR)
        + strong * np.log(BREACH_STRONG_CARVE_M_PER_MYR / BREACH_REFERENCE_CARVE_M_PER_MYR)
        + weak * np.log(BREACH_WEAK_CARVE_M_PER_MYR / BREACH_REFERENCE_CARVE_M_PER_MYR)
    )
    return np.exp(log_rate)


def channel_passage_elevation(elevation: np.ndarray, channel_depth: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Centre elevation, lowered by the cell's established channel. The result is never below
    the lowest neighbour's centre and never above the cell's own centre."""
    elevation = np.asarray(elevation, dtype=float)
    notched = elevation - np.clip(np.asarray(channel_depth, dtype=float), 0.0, None)
    if neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        return elevation.copy()
    floor = elevation[neighbor_idx].min(axis=1)
    return np.minimum(elevation, np.maximum(notched, floor))


def least_climb_to_ocean(
    passage: np.ndarray, is_ocean: np.ndarray, neighbor_idx: np.ndarray, resistance: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """`(cost_m, next_hop)` for every node, from one multi-source Dijkstra rooted at every
    ocean node. Water moves from `i` to each `j` in `neighbor_idx[i]`. That step costs
    `max(0, passage[j] - passage[i]) * resistance[j]`, plus `_HOP_EPSILON_M`. Dijkstra runs
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
    climb = np.clip(passage[cols] - passage[rows], 0.0, None) * resistance[cols]
    weight = climb + _HOP_EPSILON_M
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
) -> BreachResult:
    """Breach every pit whose least climb to the ocean fits this step's carving budget. See
    the module docstring.

    The result doesn't depend on node order. Each node's notch is the deepest cut that any
    breached path through it needs. That cut is set by the lowest water level arriving at the
    node from upstream. Nodes are visited from the highest cost down, which is a valid
    upstream-to-downstream order because each node's cost is strictly larger than its next
    hop's (see `_HOP_EPSILON_M`). One pass visits each node on the breached paths once. It is
    not a separate path search per pit."""
    elevation = np.asarray(elevation, dtype=float)
    n = len(elevation)
    dt_myr = years / 1_000_000.0
    channel_passage = channel_passage_elevation(elevation, channel_depth, neighbor_idx)
    resistance = BREACH_REFERENCE_CARVE_M_PER_MYR / np.asarray(carve_rate, dtype=float)
    cost, next_hop = least_climb_to_ocean(channel_passage, is_ocean, neighbor_idx, resistance)
    notch = np.zeros(n)
    empty = np.zeros(0, dtype=np.int64)
    if n == 0 or neighbor_idx.ndim != 2 or neighbor_idx.shape[1] == 0:
        return BreachResult(cost, next_hop, channel_passage, notch, channel_passage.copy(), empty, empty)

    budget = min(BREACH_REFERENCE_CARVE_M_PER_MYR * dt_myr, BREACH_MAX_NOTCH_M)
    pits = _pits(elevation, is_ocean, neighbor_idx)
    pit_cost = cost[pits]
    breachable = (pit_cost > _MIN_BREACH_CLIMB_M) & (pit_cost <= budget)
    breached_pits = pits[breachable]
    closed_pits = pits[np.isfinite(pit_cost) & (pit_cost > budget)]

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
        # node's cut takes its passage down to that level, but no further than its own rock
        # can be carved this step.
        arriving = np.full(n, np.inf)
        arriving[breached_pits] = channel_passage[breached_pits]
        max_cut = np.minimum(np.asarray(carve_rate, dtype=float) * dt_myr, BREACH_MAX_NOTCH_M)
        arriving_list = arriving.tolist()
        passage_list = channel_passage.tolist()
        max_cut_list = max_cut.tolist()
        next_list = next_hop.tolist()
        notch_list = notch.tolist()
        for i in order.tolist():
            level = arriving_list[i]
            cut = min(max(passage_list[i] - level, 0.0), max_cut_list[i])
            notch_list[i] = cut
            out_level = passage_list[i] - cut
            nxt = next_list[i]
            if nxt >= 0 and out_level < arriving_list[nxt]:
                arriving_list[nxt] = out_level
        notch = np.asarray(notch_list)

    return BreachResult(cost, next_hop, channel_passage, notch, channel_passage - notch, breached_pits, closed_pits)


def interface_pass_elevation(passage: np.ndarray, neighbor_idx: np.ndarray) -> np.ndarray:
    """Pass elevations aligned with `neighbor_idx`, in the form
    `lakes.build_lake_hierarchy(interface_pass_elevation=...)` takes. Water crossing from `i`
    to `j` must clear the higher of the two passages."""
    return np.maximum(passage[:, None], passage[neighbor_idx])
