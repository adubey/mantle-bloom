"""Plates laid over a cratons-first world's continents (RFC 0009, #462).

The second stage of a cratons-first layout, after `continent_layout` has placed the
continents, cores and belts on the whole-sphere `layout_lattice.Lattice`:

1. **Sites and the starting partition.** `num_plates` sites, Poisson-disc sampled outside
   the cores, take every point by multi-source Dijkstra over the 8-connected lattice (path
   length through cores, not around them). Each core then goes whole to the plate holding
   most of its area. Sites are redrawn if a plate is disconnected or under the node floor,
   and of `PARTITION_CANDIDATES` valid draws the one with the lowest J is kept.
2. **Boundary refinement.** Boundary points move between adjacent plates while that lowers
   `J`, the length of plate-crossing links weighted by the mean `PLATE_COST_*` of their
   endpoints (ocean and belt cheap, margin and core dear), plus `PLATE_COST_LAND_LINK` where
   both ends are land. Each plate stays within
   `PLATE_AREA_BAND` of its starting area and above the node floor, cores never move, and a
   local simple-point check keeps every plate connected. Sweeps run in node order until one
   lowers nothing, or `REFINE_MAX_SWEEPS`.

The RFC's J doesn't see land, so on its own it lets boundaries split continents' interiors:
boundary length with land on both sides was 21-35% over 40 seeds, against the RFC's 25%
(analysis/issue462). Two changes address it (#462 review). `PLATE_COST_LAND_LINK` adds a cost
to every link with land at both ends, from the land mask the generator passes in, and
`PARTITION_CANDIDATES` valid starting partitions are drawn and the one with the lowest J kept.
The land cost does nothing through refinement alone, which only reaches a local minimum:
choosing between whole partitions is what lets it move a boundary off a continent.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from .layout_lattice import Lattice, poisson_disc

# Boundary-link costs: a plate boundary through ocean or along a belt is cheap; through a
# continental margin or along a core's edge it costs PLATE_COST_MARGIN. Ocean and belt are
# fixed at 1; the margin cost is tuned to the land-edge measure (RFC 0009, Knobs).
PLATE_COST_OCEAN = 1.0
PLATE_COST_BELT = 1.0
PLATE_COST_MARGIN = 3.0
PLATE_COST_CORE = PLATE_COST_MARGIN
# Extra cost, per unit link length, of a boundary link with land at both ends. 10 against 5,
# both with 8 candidates, met the land, land-edge and mixed-plate targets together on 33 and
# 30 of 40 seeds (10 of 40 without either change); a cheaper Dijkstra path over continents,
# tried with it, made 29 (#462).
PLATE_COST_LAND_LINK = 10.0
# Valid starting partitions drawn; the lowest-J one is kept. 16 did no better than 8.
PARTITION_CANDIDATES = 8
# Refinement keeps each plate within this share of its starting (Dijkstra) area.
PLATE_AREA_BAND = (0.75, 1.25)
REFINE_MAX_SWEEPS = 20
# Plate sites are at least this fraction of the mean site spacing, sqrt(4 pi / n), apart.
PLATE_SITE_SPACING_FRACTION = 0.6
PLATE_SITE_MAX_TRIES = 5000
MAX_PARTITION_ATTEMPTS = 50


@dataclass
class PlatePartition:
    """`plate` labels each lattice point (0..num_plates-1); `initial_plate` is the starting
    partition before refinement; `sites` are the plate sites. `refine_costs` is J after each
    refinement sweep, first entry the starting partition; `attempts` counts the site draws."""

    plate: np.ndarray
    initial_plate: np.ndarray
    sites: np.ndarray
    refine_costs: list[float]
    attempts: int


def node_floor(node_density: float) -> int:
    from .merge_split import SPLIT_MIN_NODES

    return int(np.ceil(SPLIT_MIN_NODES * node_density))


def point_costs(continent: np.ndarray, core: np.ndarray, belt: np.ndarray) -> np.ndarray:
    """`PLATE_COST_*` per lattice point."""
    cost = np.full(len(continent), PLATE_COST_OCEAN)
    cost[continent >= 0] = PLATE_COST_MARGIN
    cost[belt] = PLATE_COST_BELT
    cost[core >= 0] = PLATE_COST_CORE
    return cost


def _link_costs(lat: Lattice, cost: np.ndarray, land: np.ndarray) -> np.ndarray:
    """Per 4-connected link, its length times the mean endpoint cost, plus
    `PLATE_COST_LAND_LINK` per unit length when both ends are land."""
    rows, cols = lat.adj_rows, lat.adj_neighbours
    return lat.adj_length_km * (0.5 * (cost[rows] + cost[cols]) + PLATE_COST_LAND_LINK * (land[rows] & land[cols]))


def boundary_cost(lat: Lattice, labels: np.ndarray, cost: np.ndarray, land: np.ndarray | None = None) -> float:
    """J: the summed cost of every plate-crossing link (see `_link_costs`)."""
    land = np.zeros(lat.size, dtype=bool) if land is None else land
    crossing = labels[lat.adj_rows] != labels[lat.adj_neighbours]
    # Each link appears once from each end.
    return 0.5 * float(_link_costs(lat, cost, land)[crossing].sum())


def label_components(lat: Lattice, labels: np.ndarray) -> int:
    """Connected components of the 4-connected lattice with links cut between labels."""
    rows = lat.adj_rows
    same = labels[rows] == labels[lat.adj_neighbours]
    graph = coo_matrix((np.ones(int(same.sum())), (rows[same], lat.adj_neighbours[same])), shape=(lat.size, lat.size))
    return int(connected_components(graph, directed=False)[0])


def _initial_partition(lat: Lattice, site_index: np.ndarray, core: np.ndarray) -> np.ndarray:
    """Nearest site by path length over the 8-connected lattice, then each core whole to the plate holding
    most of its area (lowest plate id on a tie)."""
    _, _, sources = dijkstra(lat.graph(ring=True), directed=False, indices=site_index, min_only=True, return_predecessors=True)
    plate_of_site = np.full(lat.size, -1, dtype=np.int64)
    plate_of_site[site_index] = np.arange(len(site_index))
    labels = plate_of_site[sources]
    for core_id in np.unique(core[core >= 0]):
        members = core == core_id
        shares = np.bincount(labels[members], weights=lat.area_km2[members], minlength=len(site_index))
        labels[members] = int(np.argmax(shares))
    return labels


def _plate_sites(rng: np.random.Generator, lat: Lattice, num_plates: int, core: np.ndarray) -> np.ndarray:
    """Lattice indices of `num_plates` Poisson-disc sites outside every core."""
    open_points = np.flatnonzero(core < 0)
    sep = PLATE_SITE_SPACING_FRACTION * np.sqrt(4.0 * np.pi / num_plates)
    while True:
        sites = poisson_disc(rng, num_plates, sep, lat.xyz[open_points], lat.area_km2[open_points], PLATE_SITE_MAX_TRIES)
        if len(sites) == num_plates:
            return open_points[cKDTree(lat.xyz[open_points]).query(sites)[1]]
        sep *= 0.9


@njit(cache=True)
def _is_simple(p, label, labels, adj_offsets, adj_neighbours, ring_offsets, ring_neighbours):  # pragma: no cover - numba
    """Whether `p`'s `label` edge neighbours stay connected through its 1-ring once `p`
    leaves, moving only along edge links between same-label ring points."""
    ring = ring_neighbours[ring_offsets[p] : ring_offsets[p + 1]]
    members = np.empty(len(ring), dtype=np.int64)
    m = 0
    for q in ring:
        if labels[q] == label:
            members[m] = q
            m += 1
    if m == 0:
        return True
    seen = np.zeros(m, dtype=np.bool_)
    stack = np.empty(m, dtype=np.int64)
    # Start the search from an edge neighbour of p with the label.
    start = -1
    for idx in range(m):
        q = members[idx]
        for e in range(adj_offsets[p], adj_offsets[p + 1]):
            if adj_neighbours[e] == q:
                start = idx
                break
        if start >= 0:
            break
    if start < 0:
        return True
    seen[start] = True
    stack[0] = start
    top = 1
    while top > 0:
        top -= 1
        q = members[stack[top]]
        for e in range(adj_offsets[q], adj_offsets[q + 1]):
            r = adj_neighbours[e]
            for idx in range(m):
                if not seen[idx] and members[idx] == r:
                    seen[idx] = True
                    stack[top] = idx
                    top += 1
    for e in range(adj_offsets[p], adj_offsets[p + 1]):
        q = adj_neighbours[e]
        if labels[q] != label:
            continue
        for idx in range(m):
            if members[idx] == q and not seen[idx]:
                return False
    return True


@njit(cache=True)
def _refine(
    labels, frozen, link_cost, area, adj_offsets, adj_neighbours, ring_offsets, ring_neighbours, area_lo, area_hi, node_floor, max_sweeps, j_start
):  # pragma: no cover - numba
    """Greedy boundary sweeps lowering J; returns J after each sweep (see the module docstring)."""
    n_plates = len(area_lo)
    plate_area = np.zeros(n_plates)
    plate_nodes = np.zeros(n_plates, dtype=np.int64)
    for p in range(len(labels)):
        plate_area[labels[p]] += area[p]
        plate_nodes[labels[p]] += 1
    history = [j_start]
    j = j_start
    candidates = np.empty(8, dtype=np.int64)
    released = np.empty(8)
    for _ in range(max_sweeps):
        moved = False
        for p in range(len(labels)):
            if frozen[p]:
                continue
            a = labels[p]
            # Moving p from `a` to a neighbouring plate `b` makes p's links to `a` start
            # crossing (`exposed`) and its links to `b` stop crossing (`released[b]`).
            n_cand = 0
            exposed = 0.0
            for e in range(adj_offsets[p], adj_offsets[p + 1]):
                q = adj_neighbours[e]
                if labels[q] == a:
                    exposed += link_cost[e]
                    continue
                found = -1
                for k in range(n_cand):
                    if candidates[k] == labels[q]:
                        found = k
                if found < 0:
                    candidates[n_cand] = labels[q]
                    released[n_cand] = 0.0
                    found = n_cand
                    n_cand += 1
                released[found] += link_cost[e]
            if n_cand == 0:
                continue
            if plate_area[a] - area[p] < area_lo[a] or plate_nodes[a] - 1 < node_floor:
                continue
            # The move that lowers J most (negative change), lowest plate id on a tie.
            best = -1
            best_cost_change = -1e-9
            for k in range(n_cand):
                b = candidates[k]
                cost_change = exposed - released[k]
                if plate_area[b] + area[p] > area_hi[b]:
                    continue
                if cost_change < best_cost_change or (best >= 0 and cost_change == best_cost_change and b < candidates[best]):
                    best = k
                    best_cost_change = cost_change
            if best < 0:
                continue
            if not _is_simple(p, a, labels, adj_offsets, adj_neighbours, ring_offsets, ring_neighbours):
                continue
            b = candidates[best]
            labels[p] = b
            plate_area[a] -= area[p]
            plate_area[b] += area[p]
            plate_nodes[a] -= 1
            plate_nodes[b] += 1
            j += best_cost_change
            moved = True
        history.append(j)
        if not moved:
            break
    return history


def refine_boundaries(
    lat: Lattice, labels: np.ndarray, cost: np.ndarray, frozen: np.ndarray, floor: int, land: np.ndarray | None = None
) -> tuple[np.ndarray, list[float]]:
    """Step 2: `labels` after the boundary sweeps, and J before and after each sweep."""
    land = np.zeros(lat.size, dtype=bool) if land is None else land
    labels = labels.astype(np.int64).copy()
    n_plates = int(labels.max()) + 1
    start_area = np.bincount(labels, weights=lat.area_km2, minlength=n_plates)
    history = _refine(
        labels,
        frozen.astype(np.bool_),
        _link_costs(lat, cost.astype(np.float64), land),
        lat.area_km2.astype(np.float64),
        lat.adj_offsets,
        lat.adj_neighbours,
        lat.ring_offsets,
        lat.ring_neighbours,
        PLATE_AREA_BAND[0] * start_area,
        PLATE_AREA_BAND[1] * start_area,
        floor,
        REFINE_MAX_SWEEPS,
        boundary_cost(lat, labels, cost, land),
    )
    return labels, list(history)


def partition_plates(
    rng: np.random.Generator,
    lat: Lattice,
    num_plates: int,
    continent: np.ndarray,
    core: np.ndarray,
    belt: np.ndarray,
    node_density: float = 1.0,
    land: np.ndarray | None = None,
) -> PlatePartition:
    """Plates over the continents `continent`/`core`/`belt` (per lattice point, as
    `continent_layout` builds them), deterministic in `rng`. `land` is the generator's land
    mask, which `PLATE_COST_LAND_LINK` reads; without it no link counts as land."""
    land = np.zeros(lat.size, dtype=bool) if land is None else land
    floor = node_floor(node_density)
    cost = point_costs(continent, core, belt)
    candidates: list[tuple[float, int, np.ndarray, np.ndarray]] = []
    for attempt in range(1, MAX_PARTITION_ATTEMPTS + 1):
        site_index = _plate_sites(rng, lat, num_plates, core)
        labels = _initial_partition(lat, site_index, core)
        if np.bincount(labels, minlength=num_plates).min() >= floor and label_components(lat, labels) == num_plates:
            candidates.append((boundary_cost(lat, labels, cost, land), attempt, labels, site_index))
            if len(candidates) == PARTITION_CANDIDATES:
                break
    if not candidates:
        raise RuntimeError(f"no valid {num_plates}-plate partition in {MAX_PARTITION_ATTEMPTS} attempts")
    _, attempts, initial, site_index = min(candidates, key=lambda c: (c[0], c[1]))
    labels, history = refine_boundaries(lat, initial, cost, core >= 0, floor, land)
    return PlatePartition(plate=labels, initial_plate=initial, sites=lat.xyz[site_index], refine_costs=history, attempts=attempts)
