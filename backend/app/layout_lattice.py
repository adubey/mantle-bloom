"""The whole-sphere lattice a cratons-first world is laid out on (RFC 0009, #462).

`continent_layout` (continents, cores, belts) and `plate_partition` (the plates over them)
both work on one cell lattice covering the sphere at the world's node spacing, built once
per spacing and cached. This module holds it and the sampling helpers both stages share.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.spatial import cKDTree

from .elevation_lines import PLANET_RADIUS_KM
from .sparse_quad_patch import PlateWithSparseQuadPatch


def arc_km(chord: np.ndarray) -> np.ndarray:
    """Chord length between unit vectors -> great-circle distance in km."""
    return 2.0 * np.arcsin(np.clip(chord / 2.0, 0.0, 1.0)) * PLANET_RADIUS_KM


def cap_area_km2(radius_km: float) -> float:
    """Area of a spherical cap of geodesic radius `radius_km`."""
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

    @functools.cached_property
    def adj_rows(self) -> np.ndarray:
        """The first endpoint of each 4-connected link, aligned with `adj_neighbours`."""
        return np.repeat(np.arange(self.size), np.diff(self.adj_offsets))

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
        return arc_km(np.linalg.norm(xyz[r] - xyz[c], axis=1))

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


def poisson_disc(rng: np.random.Generator, count: int, min_sep_rad: float, candidates: np.ndarray, weights: np.ndarray, max_tries: int) -> np.ndarray:
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
