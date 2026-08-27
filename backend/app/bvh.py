"""A binary tree of bounding spherical caps over a plate's own node cloud (spec section
4.3), for accelerating plate-vs-plate nearest-node/overlap queries.

Scope note: v1's `PlateWithLines.contains_batch` (the row-lookup + winding-number fallback,
see plates.py) is already O(log rows) and independently validated bit-exact against the
winding-number test across 100k+ real production points -- rewriting *that* specific query
onto a BVH would be replacing a working, already-fast, already-correct path with equivalent
performance for real risk of regression, not a genuine improvement. What this module targets
instead is the query v1 has no equivalent accelerator for at all: a real **tree-vs-tree**
traversal between two plates' own node clouds (`torque.py`'s boundary-force gathering,
`deform()`'s own neighbour-distance classification) that can prune a neighbour's *entire*
distant subtree at once, rather than -- as v1's `cKDTree`-over-concatenated-neighbour-points
approach does -- visiting a single flat index built from every neighbour's points with no
per-neighbour subtree structure to skip early. `contains_batch`/point-in-polygon containment
is untouched, inherited as-is from `PlateWithLines` (see lithosphere_plate.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import geometry

LEAF_SIZE = 16


@dataclass
class BVHNode:
    center: np.ndarray  # (3,) unit vector, mean direction of member points
    radius: float  # max angular distance (rad) from `center` to any member point
    indices: np.ndarray | None  # leaf: original indices into the points array this tree was built from
    left: "BVHNode | None" = None
    right: "BVHNode | None" = None

    @property
    def is_leaf(self) -> bool:
        return self.indices is not None


def _bounding_cap(points_xyz: np.ndarray) -> tuple[np.ndarray, float]:
    center = geometry.normalize(points_xyz.mean(axis=0, keepdims=True))[0]
    radius = float(geometry.angular_distance(points_xyz, center[None, :]).max())
    return center, radius


def _build(points_xyz: np.ndarray, indices: np.ndarray, leaf_size: int) -> BVHNode:
    center, radius = _bounding_cap(points_xyz)
    if len(indices) <= leaf_size:
        return BVHNode(center=center, radius=radius, indices=indices)

    # Median split along whichever xyz axis has the greatest spread among member points --
    # the standard BVH build heuristic (cheap, and self-balancing regardless of how the
    # points happen to be distributed across the plate's own footprint).
    spread = points_xyz.max(axis=0) - points_xyz.min(axis=0)
    axis = int(np.argmax(spread))
    order = np.argsort(points_xyz[:, axis])
    mid = len(order) // 2
    left_sel, right_sel = order[:mid], order[mid:]

    left = _build(points_xyz[left_sel], indices[left_sel], leaf_size)
    right = _build(points_xyz[right_sel], indices[right_sel], leaf_size)
    return BVHNode(center=center, radius=radius, indices=None, left=left, right=right)


def build(points_xyz: np.ndarray, leaf_size: int = LEAF_SIZE) -> BVHNode | None:
    if len(points_xyz) == 0:
        return None
    return _build(points_xyz, np.arange(len(points_xyz)), leaf_size)


def _cap_min_distance(node: BVHNode, query_xyz: np.ndarray) -> float:
    """A lower bound on the true angular distance from `query_xyz` to any point this node's
    subtree contains -- the pruning criterion: a candidate branch can be skipped outright
    once this bound already exceeds the best distance found so far."""
    center_dist = float(geometry.angular_distance(query_xyz[None, :], node.center[None, :])[0])
    return max(0.0, center_dist - node.radius)


def _query_nearest_one(node: BVHNode, points_xyz: np.ndarray, query_xyz: np.ndarray, best_idx: int, best_dist: float) -> tuple[int, float]:
    if _cap_min_distance(node, query_xyz) > best_dist:
        return best_idx, best_dist
    if node.is_leaf:
        dists = geometry.angular_distance(points_xyz[node.indices], query_xyz[None, :])
        local_best = int(np.argmin(dists))
        d = float(dists[local_best])
        if d < best_dist:
            return int(node.indices[local_best]), d
        return best_idx, best_dist

    # Descend into whichever child's cap is closer first -- makes it far more likely the
    # second child gets pruned outright by the `best_dist` bound already tightened by the
    # first, rather than visiting both subtrees' full depth regardless of order.
    left_bound = _cap_min_distance(node.left, query_xyz)
    right_bound = _cap_min_distance(node.right, query_xyz)
    first, second = (node.left, node.right) if left_bound <= right_bound else (node.right, node.left)
    best_idx, best_dist = _query_nearest_one(first, points_xyz, query_xyz, best_idx, best_dist)
    best_idx, best_dist = _query_nearest_one(second, points_xyz, query_xyz, best_idx, best_dist)
    return best_idx, best_dist


def query_nearest(tree: BVHNode | None, tree_points_xyz: np.ndarray, query_points_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For every point in `query_points_xyz`, the index (into `tree_points_xyz`) and angular
    distance (rad) of its nearest point in `tree` -- `tree` built via `build(tree_points_xyz)`.
    `(-1, inf)` per query point if `tree` is empty."""
    n = len(query_points_xyz)
    if tree is None or n == 0:
        return np.full(n, -1, dtype=int), np.full(n, np.inf)
    idx = np.empty(n, dtype=int)
    dist = np.empty(n)
    for i in range(n):
        idx[i], dist[i] = _query_nearest_one(tree, tree_points_xyz, query_points_xyz[i], -1, np.inf)
    return idx, dist


def _cap_cap_min_distance(a: BVHNode, b: BVHNode) -> float:
    center_dist = float(geometry.angular_distance(a.center[None, :], b.center[None, :])[0])
    return max(0.0, center_dist - a.radius - b.radius)


def query_nearest_cross(
    tree_a: BVHNode | None, points_a: np.ndarray, tree_b: BVHNode | None, points_b: np.ndarray, reach_rad: float
) -> list[tuple[int, int, float]]:
    """Every (index into `points_a`, index into `points_b`, angular distance) pair whose two
    points sit within `reach_rad` of each other -- the genuine tree-vs-tree traversal this
    module adds over v1's per-plate flat-KDTree approach: a whole subtree pair gets pruned
    together the moment their two bounding caps' closest possible approach already exceeds
    `reach_rad`, without either subtree's individual points ever being visited."""
    results: list[tuple[int, int, float]] = []
    if tree_a is None or tree_b is None:
        return results

    def recurse(node_a: BVHNode, node_b: BVHNode) -> None:
        if _cap_cap_min_distance(node_a, node_b) > reach_rad:
            return
        if node_a.is_leaf and node_b.is_leaf:
            pa = points_a[node_a.indices]
            pb = points_b[node_b.indices]
            d = geometry.angular_distance(pa[:, None, :], pb[None, :, :])
            hits = np.argwhere(d <= reach_rad)
            for i, j in hits:
                results.append((int(node_a.indices[i]), int(node_b.indices[j]), float(d[i, j])))
            return
        if node_a.is_leaf:
            recurse(node_a, node_b.left)
            recurse(node_a, node_b.right)
        elif node_b.is_leaf:
            recurse(node_a.left, node_b)
            recurse(node_a.right, node_b)
        else:
            recurse(node_a.left, node_b.left)
            recurse(node_a.left, node_b.right)
            recurse(node_a.right, node_b.left)
            recurse(node_a.right, node_b.right)

    recurse(tree_a, tree_b)
    return results
