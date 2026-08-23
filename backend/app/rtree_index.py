"""A minimal, self-contained R-tree spatial index over 2D points.

Bulk-loaded via STR (Sort-Tile-Recursive) packing rather than built by repeated insertion:
`PlateWithRTree` (see plates.py) rebuilds its whole index from scratch whenever its node set
changes, so there's no need for the insert/delete/rebalance machinery a general-purpose,
incrementally-updated R-tree needs -- STR packing is both simpler and produces a
better-fitted tree (less overlap between sibling boxes) than repeated insertion would anyway.
No external dependency is pulled in for this (e.g. the `rtree` package, a wrapper around the
libspatialindex C library) since STR packing plus box/nearest-neighbor queries is only a few
dozen lines of numpy.

Immutable once built (`RTree.build`) -- construct a fresh index whenever the underlying
point set changes, rather than mutating one in place.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Max entries per node before an STR-packed leaf/internal node splits into siblings -- a
# typical R-tree fanout: high enough to keep tree depth shallow (a plate's node count is at
# most a few tens of thousands, so this keeps depth to 3-4 levels), low enough that a range
# query's box-intersection test doesn't waste time on large, poorly-fitting boxes.
MAX_ENTRIES_PER_NODE = 16


@dataclass
class _Node:
    min_xy: np.ndarray
    max_xy: np.ndarray
    children: list["_Node"] | None  # None for a leaf
    point_indices: np.ndarray | None  # only set for a leaf, into the tree's own `points`


def _leaf_bbox(points: np.ndarray, indices: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    pts = points[indices]
    return pts.min(axis=0), pts.max(axis=0)


def _str_partition(items: np.ndarray, key_points: np.ndarray, m: int) -> list[np.ndarray]:
    """Partition `items` (parallel to `key_points`, one 2D point per item) into groups of at
    most `m` via Sort-Tile-Recursive: sort by x into ceil(sqrt(ceil(n/m))) vertical slices,
    then sort each slice by y and cut into groups of m. Standard STR bulk-loading, applied
    both to raw points (building leaves) and to child-node bbox centers (building each next
    level up) -- same partitioning logic either way, just a different "point" per item."""
    n = len(items)
    num_leaves = max(1, -(-n // m))  # ceil(n / m)
    num_slices = max(1, int(np.ceil(np.sqrt(num_leaves))))
    slice_size = max(1, -(-n // num_slices))  # ceil(n / num_slices)

    order_x = np.argsort(key_points[:, 0], kind="stable")
    groups = []
    for s in range(0, n, slice_size):
        slice_order = order_x[s : s + slice_size]
        slice_sorted_y = slice_order[np.argsort(key_points[slice_order, 1], kind="stable")]
        for g in range(0, len(slice_sorted_y), m):
            groups.append(items[slice_sorted_y[g : g + m]])
    return groups


def _build_leaves(points: np.ndarray, indices: np.ndarray) -> list[_Node]:
    groups = _str_partition(indices, points[indices], MAX_ENTRIES_PER_NODE)
    leaves = []
    for g in groups:
        min_xy, max_xy = _leaf_bbox(points, g)
        leaves.append(_Node(min_xy=min_xy, max_xy=max_xy, children=None, point_indices=g))
    return leaves


def _build_parent_level(nodes: list[_Node]) -> _Node:
    if len(nodes) == 1:
        return nodes[0]
    node_indices = np.arange(len(nodes))
    centers = np.array([(n.min_xy + n.max_xy) / 2.0 for n in nodes])
    groups = _str_partition(node_indices, centers, MAX_ENTRIES_PER_NODE)
    parents = []
    for g in groups:
        children = [nodes[i] for i in g]
        min_xy = np.min([c.min_xy for c in children], axis=0)
        max_xy = np.max([c.max_xy for c in children], axis=0)
        parents.append(_Node(min_xy=min_xy, max_xy=max_xy, children=children, point_indices=None))
    return _build_parent_level(parents)


def _box_intersects(node_min: np.ndarray, node_max: np.ndarray, q_min: np.ndarray, q_max: np.ndarray) -> bool:
    return bool(np.all(node_max >= q_min) and np.all(node_min <= q_max))


def _query_box(node: _Node, points: np.ndarray, q_min: np.ndarray, q_max: np.ndarray, found: list[np.ndarray]) -> None:
    if not _box_intersects(node.min_xy, node.max_xy, q_min, q_max):
        return
    if node.children is None:
        idx = node.point_indices
        pts = points[idx]
        # The leaf's own bbox already overlaps the query box (the check above), but that
        # doesn't mean every point inside the leaf does -- filter precisely rather than
        # returning the whole leaf.
        inside = np.all((pts >= q_min) & (pts <= q_max), axis=1)
        if np.any(inside):
            found.append(idx[inside])
        return
    for child in node.children:
        _query_box(child, points, q_min, q_max, found)


def _min_dist_to_box(box_min: np.ndarray, box_max: np.ndarray, q: np.ndarray) -> float:
    d = np.maximum(np.maximum(box_min - q, q - box_max), 0.0)
    return float(np.hypot(d[0], d[1]))


def _nearest_one(node: _Node, points: np.ndarray, q: np.ndarray, exclude_index: int | None, best: list) -> None:
    """`best`: mutable [best_dist, best_idx], updated in place as a better candidate is
    found. Prunes any subtree whose own bounding-box distance to `q` already exceeds the
    current best -- the same lower-bound-pruning principle real R-tree nearest-neighbor
    search uses, just depth-first (closest child first) rather than a global best-first
    priority queue, which is simpler and just as effective at this index's scale (a single
    plate's own nodes, at most a few tens of thousands)."""
    if _min_dist_to_box(node.min_xy, node.max_xy, q) > best[0]:
        return
    if node.children is None:
        idx = node.point_indices
        dist = np.hypot(*(points[idx] - q).T)
        if exclude_index is not None:
            dist = np.where(idx == exclude_index, np.inf, dist)
        local_best = int(np.argmin(dist))
        if dist[local_best] < best[0]:
            best[0] = float(dist[local_best])
            best[1] = int(idx[local_best])
        return
    # Visit the closer child first -- maximizes how often the pruning check above skips a
    # farther sibling subtree entirely once a good candidate has already been found nearby.
    for child in sorted(node.children, key=lambda c: _min_dist_to_box(c.min_xy, c.max_xy, q)):
        _nearest_one(child, points, q, exclude_index, best)


class RTree:
    """A static, bulk-loaded R-tree over 2D points -- indexed by point (a degenerate box),
    supporting box-range and nearest-neighbor queries. See the module docstring for why this
    is bulk-loaded (STR packing) and immutable rather than incrementally updated."""

    def __init__(self, points: np.ndarray, root: _Node | None) -> None:
        self._points = points
        self._root = root

    @property
    def points(self) -> np.ndarray:
        return self._points

    @classmethod
    def build(cls, points: np.ndarray) -> "RTree":
        """`points`: (n, 2) array. Empty input gives an index with no root -- every query on
        it just returns nothing, rather than needing a special case at every call site."""
        points = np.asarray(points, dtype=float)
        if len(points) == 0:
            return cls(points, None)
        indices = np.arange(len(points))
        leaves = _build_leaves(points, indices)
        return cls(points, _build_parent_level(leaves))

    def query_box(self, min_xy: np.ndarray, max_xy: np.ndarray) -> np.ndarray:
        """Indices of every point inside [min_xy, max_xy] (inclusive)."""
        if self._root is None:
            return np.zeros(0, dtype=int)
        found: list[np.ndarray] = []
        _query_box(self._root, self._points, np.asarray(min_xy, dtype=float), np.asarray(max_xy, dtype=float), found)
        if not found:
            return np.zeros(0, dtype=int)
        return np.concatenate(found)

    def count_in_box(self, min_xy: np.ndarray, max_xy: np.ndarray) -> int:
        return len(self.query_box(min_xy, max_xy))

    def nearest_one(self, query_xy: np.ndarray, exclude_index: int | None = None) -> tuple[int, float] | None:
        """The single nearest point to `query_xy` (optionally excluding one point index, for
        "nearest to point i other than i itself"). `None` for an empty index."""
        if self._root is None:
            return None
        best = [np.inf, -1]
        _nearest_one(self._root, self._points, np.asarray(query_xy, dtype=float), exclude_index, best)
        if best[1] < 0:
            return None
        return int(best[1]), float(best[0])
