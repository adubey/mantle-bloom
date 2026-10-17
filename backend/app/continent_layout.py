"""Cratons-first world layout: cores, then continents, then plates (RFC 0009, #462).

Random worlds used to be built plates first: a Voronoi tiling, each plate wholly continental
or oceanic, and a craton in the middle of every continental plate. Every continental edge was
then a plate boundary, and the tiling alone set each craton's place, size and shape. On
Earth the order runs the other way: cratons are the oldest parts of the continents, the
continents were assembled around them along orogenic belts, and today's plates are younger
and cut across that assembly. This module lays a world out in that order, before any plate
exists, on a whole-sphere lattice at the world's node spacing (`Lattice`):

1. **Continents.** `CONTINENT_COUNT_RANGE` centres, Poisson-disc sampled
   `CONTINENT_SPACING_KM` apart, split the continental area (`continental_fraction` of the
   sphere) between them by a Dirichlet draw (`CONTINENT_SIZE_CONCENTRATION`).
2. **Cores.** Each continent's core nuclei chain out from its centre, at least
   `CRATON_CORE_SPACING_KM` from every other nucleus. A core is the points within a
   noise-perturbed `CRATON_CORE_RADIUS_KM` of its nucleus, overlapping claims going to the
   nearest nucleus, less a band either side of the line halfway between two nuclei, so cores
   are at least `CRATON_BELT_WIDTH_KM` apart. Nuclei are added until the cores can cover
   `CRATON_SHARE_OF_CONTINENT` of the continent; the radius is then trimmed to hit that
   share. A continent smaller than one core disk gets none.
3. **Belts and margins.** Each continent grows outward from its cores (from its centre if it
   has none) by a noise-perturbed distance until it reaches its area. The closure of its cores
   (dilated, then eroded, by `BELT_CLOSURE_KM`) is the core-and-belt interior: everything in
   it that isn't core is belt, so the whole gap between neighbouring cores is belt and no
   unmarked gap is left to form craton later. Continental points outside it are margin.
4. **Plates.** `num_plates` sites, Poisson-disc sampled outside the cores, take every point
   by multi-source Dijkstra over the 8-connected lattice (path length through cores, not
   around them). Each core then goes whole to the plate holding most of its area. Sites are
   resampled if a plate is disconnected or under the node floor.
5. **Boundary refinement.** Boundary points move between adjacent plates while that lowers
   `J`, the length of plate-crossing links weighted by the mean `PLATE_COST_*` of their
   endpoints (ocean and belt cheap, margin and core dear). Each plate stays within
   `PLATE_AREA_BAND` of its starting area and above the node floor, cores never move, and a
   local simple-point check keeps every plate connected. Sweeps run in node order until one
   lowers nothing, or `REFINE_MAX_SWEEPS`.

`lithosphere_plate.generate_plates` builds mixed-crust plates from the result, behind
`CRATONS_FIRST_DEFAULT` until RFC 0009's gate G3 decides whether it becomes the default.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from .elevation_lines import PLANET_RADIUS_KM
from .noise import SphereNoise
from .sparse_quad_patch import PlateWithSparseQuadPatch

# Off until RFC 0009's gate G3 picks the layout: main's random worlds stay plates-first.
CRATONS_FIRST_DEFAULT = False

# --- Continents --------------------------------------------------------------------------

# Design choices (RFC 0009 knobs), tuned to the layout measures over 20 seeds.
CONTINENT_COUNT_RANGE = (2, 6)
CONTINENT_SPACING_KM = 4000.0
# Dirichlet concentration: 2 gives clearly unequal continents without slivers.
CONTINENT_SIZE_CONCENTRATION = 2.0
# Continental area share of the sphere is clamped here so there is always some of each.
CONTINENTAL_FRACTION_BOUNDS = (0.05, 0.9)
# Amplitude of the noise on the distance a continent grows outward from its cores, so
# coastlines and margins aren't parallel to the core outlines.
CONTINENT_OUTLINE_NOISE_KM = 400.0

# --- Cores and belts ---------------------------------------------------------------------

# Fixed to today's mean starting craton share (37.6-50.7% of genuine continental area on
# seeds 7, 11 and 2024, #403) so craton calibrations carry over. Capped so every continent
# keeps some non-cratonic crust.
CRATON_SHARE_OF_CONTINENT = 0.4
CRATON_SHARE_MAX = 0.7
# Minimum nucleus separation -- not a count of geological cratons.
CRATON_CORE_SPACING_KM = 1500.0
# About 3.1 M km^2 per full disk.
CRATON_CORE_RADIUS_KM = 1000.0
# Orogens welded onto shields are a few hundred km wide (Cawood et al., 2009).
CRATON_BELT_WIDTH_KM = 300.0
# Noise amplitude on a point's distance to its nucleus, for irregular core outlines.
CRATON_CORE_OUTLINE_NOISE_KM = 250.0
# A new nucleus is placed this far (in multiples of the core spacing) from an existing
# nucleus of the same continent, so a continent's cores stay a connected assembly.
NUCLEUS_STEP_RANGE = (1.0, 1.25)
NUCLEUS_MAX_TRIES = 200
# Closing radius for the core-and-belt interior: bridges gaps between cores up to twice
# this wide, more than the widest gap the nucleus step and a trimmed core radius leave. At
# 600 km belts were 9.7-16% of continental area over 20 seeds, under the RFC's 10-30%; at
# 800 km, 10.8-18% (#462).
BELT_CLOSURE_KM = 800.0
# The sinusoid noise's lowest frequency (cycles per radian): wavelengths of ~5,000 km down
# to ~1,250 km over its three octaves, the scale of a core's or continent's outline.
OUTLINE_NOISE_BASE_FREQ = 8.0

# --- Plates ------------------------------------------------------------------------------

# Boundary-link costs: a plate boundary through ocean or along a belt is cheap; through a
# continental margin or along a core's edge it costs PLATE_COST_MARGIN. Ocean and belt are
# fixed at 1; the margin cost is tuned to the land-edge measure (RFC 0009, Knobs).
PLATE_COST_OCEAN = 1.0
PLATE_COST_BELT = 1.0
PLATE_COST_MARGIN = 3.0
PLATE_COST_CORE = PLATE_COST_MARGIN
# Refinement keeps each plate within this share of its starting (Dijkstra) area.
PLATE_AREA_BAND = (0.75, 1.25)
REFINE_MAX_SWEEPS = 20
# Plate sites are at least this fraction of the mean site spacing, sqrt(4 pi / n), apart.
PLATE_SITE_SPACING_FRACTION = 0.6
PLATE_SITE_MAX_TRIES = 5000
MAX_PARTITION_ATTEMPTS = 50


def _arc_km(chord: np.ndarray) -> np.ndarray:
    return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)) * PLANET_RADIUS_KM


def _cap_area_km2(radius_km: float) -> float:
    return 2.0 * np.pi * PLANET_RADIUS_KM**2 * (1.0 - np.cos(radius_km / PLANET_RADIUS_KM))


@dataclass(frozen=True)
class Lattice:
    """The whole sphere's cells at one node spacing, in the identity frame.

    `adj_offsets`/`adj_neighbours` are the 4-connected (edge-sharing) CSR adjacency and
    `adj_length_km` each link's centre-to-centre length; `ring_offsets`/`ring_neighbours`
    add the diagonal neighbours (the 8-connected 1-ring), with `ring_length_km`."""

    xyz: np.ndarray
    area_km2: np.ndarray
    adj_offsets: np.ndarray
    adj_neighbours: np.ndarray
    adj_length_km: np.ndarray
    ring_offsets: np.ndarray
    ring_neighbours: np.ndarray
    ring_length_km: np.ndarray

    @property
    def size(self) -> int:
        return len(self.xyz)

    @functools.cached_property
    def tree(self) -> cKDTree:
        return cKDTree(self.xyz)

    def graph(self, ring: bool = False) -> csr_matrix:
        offsets, neighbours, lengths = (
            (self.ring_offsets, self.ring_neighbours, self.ring_length_km) if ring else (self.adj_offsets, self.adj_neighbours, self.adj_length_km)
        )
        return csr_matrix((lengths, neighbours, offsets), shape=(self.size, self.size))


@functools.lru_cache(maxsize=4)
def lattice(spacing_rad: float) -> Lattice:
    """The `Lattice` at `spacing_rad`, cached: one full-sphere cell patch per node density."""
    patch = PlateWithSparseQuadPatch.from_lattice(0, np.eye(3), "oceanic", spacing_rad, lambda pts: np.ones(len(pts), dtype=bool))
    xyz = patch.surface_nodes().world_xyz.copy()
    n = len(xyz)
    adjacency = patch.adjacency()
    offsets, neighbours = np.asarray(adjacency.offsets), np.asarray(adjacency.neighbours)
    rows = np.repeat(np.arange(n), np.diff(offsets))
    a4 = coo_matrix((np.ones(len(rows)), (rows, neighbours)), shape=(n, n)).tocsr()
    # A diagonal neighbour is two edge hops away by two different paths; a straight two-hop
    # neighbour is reached by only one.
    two_hop = (a4 @ a4).tocoo()
    diagonal = (two_hop.data >= 2) & (two_hop.row != two_hop.col)
    d_rows, d_cols = two_hop.row[diagonal], two_hop.col[diagonal]
    keep = ~np.asarray(a4[d_rows, d_cols]).ravel().astype(bool)
    ring = coo_matrix(
        (np.ones(len(rows) + int(keep.sum())), (np.concatenate([rows, d_rows[keep]]), np.concatenate([neighbours, d_cols[keep]]))), shape=(n, n)
    ).tocsr()
    ring.sort_indices()
    ring_rows = np.repeat(np.arange(n), np.diff(ring.indptr))

    def lengths(r: np.ndarray, c: np.ndarray) -> np.ndarray:
        return _arc_km(np.linalg.norm(xyz[r] - xyz[c], axis=1))

    arrays = dict(
        xyz=xyz,
        area_km2=patch.node_areas_m2() / 1e6,
        adj_offsets=offsets.astype(np.int64),
        adj_neighbours=neighbours.astype(np.int64),
        adj_length_km=lengths(rows, neighbours),
        ring_offsets=ring.indptr.astype(np.int64),
        ring_neighbours=ring.indices.astype(np.int64),
        ring_length_km=lengths(ring_rows, ring.indices),
    )
    for values in arrays.values():
        values.flags.writeable = False
    return Lattice(**arrays)


@dataclass
class ContinentLayout:
    """Per lattice point: `continent` (-1 for ocean), `core` (core id, -1 outside cores),
    `belt` and `margin` masks, `plate` (0..num_plates-1) and `initial_plate`, the partition
    before boundary refinement. `sites` are the plate sites,
    `nuclei` the core nuclei (`nucleus_continent` names each one's continent), and
    `continent_area_km2` each continent's target area. `refine_costs` is J after each
    refinement sweep, first entry the starting partition."""

    lattice: Lattice
    continent: np.ndarray
    core: np.ndarray
    belt: np.ndarray
    margin: np.ndarray
    plate: np.ndarray
    initial_plate: np.ndarray
    sites: np.ndarray
    nuclei: np.ndarray
    nucleus_continent: np.ndarray
    continent_area_km2: np.ndarray
    refine_costs: list[float]
    partition_attempts: int

    @property
    def is_continental(self) -> np.ndarray:
        return self.continent >= 0

    @property
    def is_core(self) -> np.ndarray:
        return self.core >= 0

    @property
    def num_plates(self) -> int:
        return len(self.sites)

    def nearest(self, points_xyz: np.ndarray) -> np.ndarray:
        """Index of the lattice point nearest each of `points_xyz`."""
        return self.lattice.tree.query(points_xyz, workers=-1)[1]

    def plate_mask(self, plate_id: int, points_xyz: np.ndarray) -> np.ndarray:
        """Whether each of `points_xyz` falls in plate `plate_id`: its nearest lattice point
        is labelled with it. Points outside the plate's bounding cap around its site (padded
        by a cell) are rejected before the nearest-point query."""
        owned = np.zeros(len(points_xyz), dtype=bool)
        near = np.flatnonzero(points_xyz @ self.sites[plate_id] >= self._cap_cos[plate_id])
        if len(near):
            owned[near] = self.plate[self.nearest(points_xyz[near])] == plate_id
        return owned

    @functools.cached_property
    def _cap_cos(self) -> np.ndarray:
        pad = 2.0 * float(np.sqrt(self.lattice.area_km2.max())) / PLANET_RADIUS_KM
        angle = np.arccos(np.clip(np.einsum("ij,ij->i", self.lattice.xyz, self.sites[self.plate]), -1.0, 1.0))
        reach = np.zeros(self.num_plates)
        np.maximum.at(reach, self.plate, angle)
        return np.cos(np.minimum(reach + pad, np.pi))

    def diagnostics(self) -> dict:
        area = self.lattice.area_km2
        continental = float(area[self.is_continental].sum())
        per_continent = []
        disk = _cap_area_km2(CRATON_CORE_RADIUS_KM)
        for k, target in enumerate(self.continent_area_km2):
            in_k = self.continent == k
            per_continent.append(
                {
                    "area_km2": float(area[in_k].sum()),
                    "target_km2": float(target),
                    "cores": int(len(np.unique(self.core[in_k & self.is_core]))),
                    "room_for_two_cores": bool(CRATON_SHARE_OF_CONTINENT * target > disk),
                }
            )
        return {
            "continental_area_fraction": continental / float(area.sum()),
            "core_share_of_continent": float(area[self.is_core].sum()) / max(continental, 1.0),
            "belt_share_of_continent": float(area[self.belt].sum()) / max(continental, 1.0),
            "margin_share_of_continent": float(area[self.margin].sum()) / max(continental, 1.0),
            "continents": per_continent,
            "refine_sweeps": len(self.refine_costs) - 1,
            "partition_attempts": self.partition_attempts,
        }


# --- Sampling helpers --------------------------------------------------------------------


def _poisson_disc(
    rng: np.random.Generator, count: int, min_sep_rad: float, candidates: np.ndarray, weights: np.ndarray, max_tries: int
) -> np.ndarray:
    """Up to `count` of `candidates` (unit vectors), drawn by `weights` and kept when at
    least `min_sep_rad` from every point already kept -- dart throwing."""
    p = weights / weights.sum()
    kept: list[np.ndarray] = []
    min_dot = np.cos(min_sep_rad)
    for _ in range(max_tries):
        if len(kept) == count:
            break
        candidate = candidates[rng.choice(len(candidates), p=p)]
        if not kept or float(np.max(np.asarray(kept) @ candidate)) < min_dot:
            kept.append(candidate)
    return np.asarray(kept).reshape(-1, 3)


def _step_from(point: np.ndarray, distance_rad: float, bearing: float) -> np.ndarray:
    """The point `distance_rad` from `point` along `bearing` (any fixed tangent basis)."""
    helper = np.array([0.0, 0.0, 1.0]) if abs(point[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(point, helper)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(point, e1)
    direction = np.cos(bearing) * e1 + np.sin(bearing) * e2
    return np.cos(distance_rad) * point + np.sin(distance_rad) * direction


def _fill_to_area(order_values: np.ndarray, area: np.ndarray, target: float) -> np.ndarray:
    """Boolean mask of the points with the smallest `order_values` whose area first reaches
    `target` (all of them when they can't)."""
    order = np.argsort(order_values, kind="stable")
    cumulative = np.cumsum(area[order])
    count = int(np.searchsorted(cumulative, target)) + 1
    mask = np.zeros(len(order_values), dtype=bool)
    mask[order[: min(count, len(order))]] = True
    return mask


# --- 1-3: continents, cores, belts -------------------------------------------------------


class _NucleusDistances:
    """Per lattice point, the distance in km to its nearest and second-nearest nucleus and
    which nucleus is nearest, updated one nucleus at a time."""

    def __init__(self, lat: Lattice, core_noise_km: np.ndarray) -> None:
        self.lat = lat
        self.noise = core_noise_km
        self.d1 = np.full(lat.size, np.inf)
        self.d2 = np.full(lat.size, np.inf)
        self.nearest = np.zeros(lat.size, dtype=np.int64)
        self.count = 0

    def add(self, nucleus: np.ndarray) -> None:
        d = np.arccos(np.clip(self.lat.xyz @ nucleus, -1.0, 1.0)) * PLANET_RADIUS_KM
        closer = d < self.d1
        self.d2 = np.where(closer, self.d1, np.minimum(self.d2, d))
        self.d1 = np.where(closer, d, self.d1)
        self.nearest = np.where(closer, self.count, self.nearest)
        self.count += 1

    def eligible(self, owners: np.ndarray, k: int) -> np.ndarray:
        """Points continent `k`'s cores can claim: nearest to one of its nuclei, within the
        perturbed core radius, and clear of the belt band. Half the difference of the two
        nearest distances is at most the distance to their bisector, so cores stay at least
        the belt width apart."""
        clears_belt = (self.d2 - self.d1) / 2.0 >= CRATON_BELT_WIDTH_KM / 2.0
        return (owners[self.nearest] == k) & clears_belt & (self.perturbed <= CRATON_CORE_RADIUS_KM)

    @property
    def perturbed(self) -> np.ndarray:
        return self.d1 + self.noise


def _place_nuclei(rng: np.random.Generator, field: _NucleusDistances, centres: np.ndarray, core_targets: np.ndarray):
    """Nuclei for every continent with a core target, chained out from its centre until its
    cores could cover the target. A target larger than one core disk always gets at least two
    nuclei, even where the noise would let one core stretch to cover it."""
    nuclei: list[np.ndarray] = []
    owners: list[int] = []
    spacing_rad = CRATON_CORE_SPACING_KM / PLANET_RADIUS_KM
    min_dot = np.cos(spacing_rad)
    disk = _cap_area_km2(CRATON_CORE_RADIUS_KM)

    def accept(candidate: np.ndarray, k: int) -> bool:
        nuclei.append(candidate)
        owners.append(k)
        field.add(candidate)
        if core_targets[k] > disk and owners.count(k) < 2:
            return False
        return float(field.lat.area_km2[field.eligible(np.asarray(owners), k)].sum()) >= core_targets[k]

    for k, target in enumerate(core_targets):
        if target <= 0.0:
            continue
        if accept(centres[k], k):
            continue
        for _ in range(NUCLEUS_MAX_TRIES):
            mine = [i for i, o in enumerate(owners) if o == k]
            parent = nuclei[mine[int(rng.integers(len(mine)))]]
            step = rng.uniform(*NUCLEUS_STEP_RANGE) * spacing_rad
            candidate = _step_from(parent, step, rng.uniform(0.0, 2.0 * np.pi))
            if float(np.max(np.asarray(nuclei) @ candidate)) >= min_dot:
                continue
            if int(np.argmax(centres @ candidate)) != k:
                continue
            if accept(candidate, k):
                break
    return np.asarray(nuclei).reshape(-1, 3), np.asarray(owners, dtype=np.int64)


def _cores(field: _NucleusDistances, owners: np.ndarray, core_targets: np.ndarray) -> np.ndarray:
    """Core id (nucleus index) per lattice point, -1 outside every core: each continent's
    eligible points nearest their nucleus, up to its core target."""
    core = np.full(field.lat.size, -1, dtype=np.int64)
    perturbed = field.perturbed
    for k in np.unique(owners):
        eligible = np.flatnonzero(field.eligible(owners, k))
        chosen = eligible[_fill_to_area(perturbed[eligible], field.lat.area_km2[eligible], core_targets[k])]
        core[chosen] = field.nearest[chosen]
    return core


def _continents(lat: Lattice, centres: np.ndarray, areas: np.ndarray, core: np.ndarray, owners: np.ndarray, outline_noise_km: np.ndarray):
    """(continent id per point, -1 for ocean; raw distance in km from each point to its own
    continent's cores, or centre when it has none)."""
    n_c = len(centres)
    distance = np.empty((lat.size, n_c))
    for k in range(n_c):
        in_core = (core >= 0) & (owners[np.maximum(core, 0)] == k)
        if np.any(in_core):
            distance[:, k] = _arc_km(cKDTree(lat.xyz[in_core]).query(lat.xyz, workers=-1)[0])
        else:
            distance[:, k] = np.arccos(np.clip(lat.xyz @ centres[k], -1.0, 1.0)) * PLANET_RADIUS_KM
    perturbed = distance + outline_noise_km[:, None]
    owner = np.argmin(perturbed, axis=1)
    continent = np.full(lat.size, -1, dtype=np.int64)
    for k in range(n_c):
        in_core = (core >= 0) & (owners[np.maximum(core, 0)] == k)
        continent[in_core] = k
        rest = np.flatnonzero((owner == k) & ~in_core)
        remaining = areas[k] - float(lat.area_km2[in_core].sum())
        if remaining > 0.0 and len(rest):
            continent[rest[_fill_to_area(perturbed[rest, k], lat.area_km2[rest], remaining)]] = k
    return continent, distance[np.arange(lat.size), np.maximum(continent, 0)]


def _belts(lat: Lattice, continent: np.ndarray, core: np.ndarray, distance_to_cores_km: np.ndarray) -> np.ndarray:
    """The closure of each continent's cores, less the cores themselves."""
    belt = np.zeros(lat.size, dtype=bool)
    has_core = core >= 0
    for k in np.unique(continent[has_core]):
        # The distance is to each point's own continent's cores, so only k's points dilate.
        dilated = (continent == k) & (distance_to_cores_km <= BELT_CLOSURE_KM)
        outside = ~dilated
        if not np.any(outside):
            continue
        candidates = np.flatnonzero(dilated & ~has_core)
        if len(candidates) == 0:
            continue
        gap = _arc_km(cKDTree(lat.xyz[outside]).query(lat.xyz[candidates], workers=-1)[0])
        belt[candidates[gap > BELT_CLOSURE_KM]] = True
    return belt


# --- 4-5: plates -------------------------------------------------------------------------


def _node_floor(node_density: float) -> int:
    from .merge_split import SPLIT_MIN_NODES

    return int(np.ceil(SPLIT_MIN_NODES * node_density))


def _label_components(lat: Lattice, labels: np.ndarray) -> int:
    """Connected components of the 4-connected lattice with links cut between labels."""
    rows = np.repeat(np.arange(lat.size), np.diff(lat.adj_offsets))
    same = labels[rows] == labels[lat.adj_neighbours]
    graph = coo_matrix((np.ones(int(same.sum())), (rows[same], lat.adj_neighbours[same])), shape=(lat.size, lat.size))
    return int(connected_components(graph, directed=False)[0])


def _initial_partition(lat: Lattice, site_index: np.ndarray, core: np.ndarray) -> np.ndarray:
    """Nearest site by path length over the 8-connected lattice, then each core whole to the
    plate holding most of its area (lowest plate id on a tie)."""
    _, _, sources = dijkstra(lat.graph(ring=True), directed=False, indices=site_index, min_only=True, return_predecessors=True)
    plate_of_site = np.full(lat.size, -1, dtype=np.int64)
    plate_of_site[site_index] = np.arange(len(site_index))
    labels = plate_of_site[sources]
    for core_id in np.unique(core[core >= 0]):
        members = core == core_id
        shares = np.bincount(labels[members], weights=lat.area_km2[members], minlength=len(site_index))
        labels[members] = int(np.argmax(shares))
    return labels


def point_costs(continent: np.ndarray, core: np.ndarray, belt: np.ndarray) -> np.ndarray:
    """`PLATE_COST_*` per lattice point."""
    cost = np.full(len(continent), PLATE_COST_OCEAN)
    cost[continent >= 0] = PLATE_COST_MARGIN
    cost[belt] = PLATE_COST_BELT
    cost[core >= 0] = PLATE_COST_CORE
    return cost


def boundary_cost(lat: Lattice, labels: np.ndarray, cost: np.ndarray) -> float:
    """J: the length of every plate-crossing link times the mean cost of its endpoints."""
    rows = np.repeat(np.arange(lat.size), np.diff(lat.adj_offsets))
    crossing = labels[rows] != labels[lat.adj_neighbours]
    weights = lat.adj_length_km * 0.5 * (cost[rows] + cost[lat.adj_neighbours])
    # Each link appears once from each end.
    return 0.5 * float(weights[crossing].sum())


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
    labels, frozen, cost, area, adj_offsets, adj_neighbours, adj_length, ring_offsets, ring_neighbours, area_lo, area_hi, node_floor, max_sweeps, j_start
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
    gains = np.empty(8)
    for _ in range(max_sweeps):
        moved = False
        for p in range(len(labels)):
            if frozen[p]:
                continue
            a = labels[p]
            # Candidate destination labels among p's edge neighbours, with each move's change
            # in J: links to `a` start crossing, links to `b` stop.
            n_cand = 0
            to_a = 0.0
            for e in range(adj_offsets[p], adj_offsets[p + 1]):
                q = adj_neighbours[e]
                c = adj_length[e] * 0.5 * (cost[p] + cost[q])
                if labels[q] == a:
                    to_a += c
                    continue
                found = -1
                for k in range(n_cand):
                    if candidates[k] == labels[q]:
                        found = k
                if found < 0:
                    candidates[n_cand] = labels[q]
                    gains[n_cand] = 0.0
                    found = n_cand
                    n_cand += 1
                gains[found] += c
            if n_cand == 0:
                continue
            if plate_area[a] - area[p] < area_lo[a] or plate_nodes[a] - 1 < node_floor:
                continue
            # Best move first (largest drop, then lowest plate id).
            best = -1
            best_delta = -1e-9
            for k in range(n_cand):
                b = candidates[k]
                delta = to_a - gains[k]
                if plate_area[b] + area[p] > area_hi[b]:
                    continue
                if delta < best_delta or (best >= 0 and delta == best_delta and b < candidates[best]):
                    best = k
                    best_delta = delta
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
            j += best_delta
            moved = True
        history.append(j)
        if not moved:
            break
    return history


def refine_boundaries(lat: Lattice, labels: np.ndarray, cost: np.ndarray, frozen: np.ndarray, node_floor: int) -> tuple[np.ndarray, list[float]]:
    """Step 5: `labels` after the boundary sweeps, and J before and after each sweep."""
    labels = labels.astype(np.int64).copy()
    n_plates = int(labels.max()) + 1
    start_area = np.bincount(labels, weights=lat.area_km2, minlength=n_plates)
    history = _refine(
        labels,
        frozen.astype(np.bool_),
        cost.astype(np.float64),
        lat.area_km2.astype(np.float64),
        lat.adj_offsets,
        lat.adj_neighbours,
        lat.adj_length_km,
        lat.ring_offsets,
        lat.ring_neighbours,
        PLATE_AREA_BAND[0] * start_area,
        PLATE_AREA_BAND[1] * start_area,
        node_floor,
        REFINE_MAX_SWEEPS,
        boundary_cost(lat, labels, cost),
    )
    return labels, list(history)


def _plate_sites(rng: np.random.Generator, lat: Lattice, num_plates: int, core: np.ndarray) -> np.ndarray:
    """Lattice indices of `num_plates` Poisson-disc sites outside every core."""
    open_points = np.flatnonzero(core < 0)
    sep = PLATE_SITE_SPACING_FRACTION * np.sqrt(4.0 * np.pi / num_plates)
    while True:
        sites = _poisson_disc(rng, num_plates, sep, lat.xyz[open_points], lat.area_km2[open_points], PLATE_SITE_MAX_TRIES)
        if len(sites) == num_plates:
            return open_points[cKDTree(lat.xyz[open_points]).query(sites)[1]]
        sep *= 0.9


def _partition(rng: np.random.Generator, lat: Lattice, num_plates: int, core: np.ndarray, node_floor: int) -> tuple[np.ndarray, np.ndarray, int]:
    """Step 4: (labels, site lattice indices, attempts used)."""
    for attempt in range(1, MAX_PARTITION_ATTEMPTS + 1):
        site_index = _plate_sites(rng, lat, num_plates, core)
        labels = _initial_partition(lat, site_index, core)
        counts = np.bincount(labels, minlength=num_plates)
        if counts.min() >= node_floor and _label_components(lat, labels) == num_plates:
            return labels, site_index, attempt
    raise RuntimeError(f"no valid {num_plates}-plate partition in {MAX_PARTITION_ATTEMPTS} attempts")


def build_layout(rng: np.random.Generator, spacing_rad: float, num_plates: int, continental_fraction: float, node_density: float = 1.0) -> ContinentLayout:
    """The cratons-first layout (see the module docstring), deterministic in `rng`."""
    lat = lattice(spacing_rad)
    continental_fraction = float(np.clip(continental_fraction, *CONTINENTAL_FRACTION_BOUNDS))
    share = min(CRATON_SHARE_OF_CONTINENT, CRATON_SHARE_MAX)

    core_noise_km = CRATON_CORE_OUTLINE_NOISE_KM * SphereNoise(rng, base_freq=OUTLINE_NOISE_BASE_FREQ).sample(lat.xyz)
    outline_noise_km = CONTINENT_OUTLINE_NOISE_KM * SphereNoise(rng, base_freq=OUTLINE_NOISE_BASE_FREQ).sample(lat.xyz)

    n_c = int(rng.integers(CONTINENT_COUNT_RANGE[0], CONTINENT_COUNT_RANGE[1] + 1))
    centres = _poisson_disc(rng, n_c, CONTINENT_SPACING_KM / PLANET_RADIUS_KM, lat.xyz, lat.area_km2, PLATE_SITE_MAX_TRIES)
    total_area = continental_fraction * float(lat.area_km2.sum())
    areas = total_area * rng.dirichlet(np.full(len(centres), CONTINENT_SIZE_CONCENTRATION))
    disk = _cap_area_km2(CRATON_CORE_RADIUS_KM)
    core_targets = np.where(areas >= disk, share * areas, 0.0)

    field = _NucleusDistances(lat, core_noise_km)
    nuclei, owners = _place_nuclei(rng, field, centres, core_targets)
    core = _cores(field, owners, core_targets)
    continent, distance_to_cores = _continents(lat, centres, areas, core, owners, outline_noise_km)
    belt = _belts(lat, continent, core, distance_to_cores)
    margin = (continent >= 0) & (core < 0) & ~belt

    floor = _node_floor(node_density)
    initial, site_index, attempts = _partition(rng, lat, num_plates, core, floor)
    labels, history = refine_boundaries(lat, initial, point_costs(continent, core, belt), core >= 0, floor)
    return ContinentLayout(
        lattice=lat,
        continent=continent,
        core=core,
        belt=belt,
        margin=margin,
        plate=labels,
        initial_plate=initial,
        sites=lat.xyz[site_index],
        nuclei=nuclei,
        nucleus_continent=owners,
        continent_area_km2=areas,
        refine_costs=history,
        partition_attempts=attempts,
    )
