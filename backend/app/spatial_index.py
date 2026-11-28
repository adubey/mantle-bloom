"""Shared 3D point index with `cKDTree`'s query semantics, backed by Rust when available (#1042).

`SpatialIndex(points)` answers the `scipy.spatial.cKDTree` queries the simulation and renderer
use -- `query` (k nearest, optionally bounded), `query_ball_point` and `query_pairs` -- with the
same results, and adds `query_ball_point_csr`, which returns one flat CSR `(indptr, indices)` pair
instead of a Python list per query point. With `mantle_bloom_rust` installed it is the native
`KdIndex` (`backend/rust/src/spatial/`), whose batch queries run across every core without the
GIL; without it, or with `MANTLE_BLOOM_DISABLE_RUST=1`, it wraps a `cKDTree`.

What both backends guarantee, so either can stand in for the other or for a bare `cKDTree`:

* Distances are `cKDTree`'s, bit for bit: the square root of `dx² + dy² + dz²` summed in that
  order with the installed scipy's own rounding, which `_scipy_fuses` reads once at import.
  Every threshold test therefore agrees too: `query` keeps points strictly closer than
  `distance_upper_bound` (misses are `inf` with index `n`), and radius queries keep points at
  most `r` away.
* `query` rows are sorted by distance; `query_ball_point` rows are ascending, single-point
  queries included (scipy leaves those in traversal order); `query_pairs(output_type="ndarray")`
  is in lexicographic order (scipy's is traversal order).
* Exact ties go to the lowest index, on both backends: `query` rows are ordered by
  `(squared distance, index)`, and when tied points straddle the `k`-th slot the lowest-index
  ones are kept. Ordering is decided on the squared distance, as `cKDTree` decides it, so two
  points whose distances differ only before the square root are not a tie. A bare `cKDTree` keeps whichever tied point its traversal meets first, which depends on
  its leaf size and build options, so its rows can differ from these exactly where two distinct
  points are equidistant from the query. That is not rare: lattice and grid symmetry makes
  neighbours bit-identically equidistant somewhere in every step (see
  `bin/debug/compare_spatial_backends.py` and `analysis/issue1042/README.md`).

Only Euclidean (`p=2`), exact (`eps=0`) queries on `(n, 3)` points are supported. `workers` is
accepted for drop-in compatibility; the Rust backend picks its own parallelism from the batch
size, which never changes a result.
"""

from __future__ import annotations

import itertools
import logging
import math
import os

import numpy as np
from scipy.spatial import cKDTree

_logger = logging.getLogger(__name__)


def _import_kdindex():
    if os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1":
        return None
    try:
        import mantle_bloom_rust
    except ImportError:
        return None
    kdindex = getattr(mantle_bloom_rust, "KdIndex", None)
    if kdindex is None:
        _logger.warning(
            "mantle_bloom_rust has no KdIndex (stale build?); spatial queries fall back to scipy. "
            "Rebuild with `backend/.venv/bin/python bin/build_rust.py` or `./bin/setup.sh`."
        )
    return kdindex


def _squared_distances(queries: np.ndarray, points: np.ndarray, fused: bool) -> np.ndarray:
    """`cKDTree`'s squared distance between paired rows, rounded as the installed scipy rounds
    it: `(dx*dx + dy*dy) + dz*dz`, or with the last two products fused into multiply-adds (see
    `kdtree::Metric` in the Rust crate). numpy has no fused multiply-add, so the fused form runs
    `math.fma` per row; the fallback calls this only for rows whose rounded distances collide."""
    diff = queries - points
    if not fused:
        return diff[:, 0] * diff[:, 0] + diff[:, 1] * diff[:, 1] + diff[:, 2] * diff[:, 2]
    fma = math.fma
    return np.array([fma(dz, dz, fma(dy, dy, dx * dx)) for dx, dy, dz in diff.tolist()], dtype=np.float64)


def _scipy_fuses() -> bool | None:
    """Whether the installed scipy fuses its squared-distance sum, read by comparing `cKDTree`
    distances on a fixed random cloud against both roundings. About one query in nine rounds
    differently under the two, so 512 queries always tell them apart; `None` if neither matches."""
    rng = np.random.default_rng(1042)
    points = rng.normal(size=(64, 3))
    queries = rng.normal(size=(512, 3))
    expected, idx = cKDTree(points).query(queries)
    for fused in (True, False):
        if np.array_equal(np.sqrt(_squared_distances(queries, points[idx], fused)), expected):
            return fused
    _logger.warning("scipy's cKDTree distances match neither rounding; spatial-index tie order may differ by backend.")
    return None


def _native_matches_scipy(kdindex, fused: bool) -> bool:
    rng = np.random.default_rng(1042)
    points = rng.normal(size=(64, 3))
    queries = rng.normal(size=(512, 3))
    try:
        got, _ = kdindex(points, fused).query(queries, 1, np.inf)
    except Exception as exc:  # noqa: BLE001 -- a stale build with another signature must not break import
        _logger.warning("mantle_bloom_rust KdIndex failed its self-check (%s; stale build?); spatial queries fall back to scipy.", exc)
        return False
    if np.array_equal(got, cKDTree(points).query(queries)[0]):
        return True
    _logger.warning("mantle_bloom_rust KdIndex distances differ from scipy's; spatial queries fall back to scipy.")
    return False


FUSED: bool | None = _scipy_fuses()
_KdIndex = _import_kdindex()
if _KdIndex is not None and (FUSED is None or not _native_matches_scipy(_KdIndex, FUSED)):
    _KdIndex = None
RUST_AVAILABLE: bool = _KdIndex is not None


def _check_supported(p: float, eps: float) -> None:
    if p != 2 or eps != 0:
        raise ValueError("SpatialIndex supports only exact Euclidean queries (p=2, eps=0)")


def _as_queries(x) -> tuple[np.ndarray, tuple[int, ...], bool]:
    """`x` as contiguous `(m, 3)` float64 rows, plus its leading shape and whether it was one point."""
    arr = np.asarray(x, dtype=np.float64)
    if arr.ndim == 0 or arr.shape[-1] != 3:
        raise ValueError(f"query points must have a last dimension of 3, got shape {arr.shape}")
    return np.ascontiguousarray(arr.reshape(-1, 3)), arr.shape[:-1], arr.ndim == 1


def _radii(r, lead: tuple[int, ...]) -> np.ndarray:
    radii = np.asarray(r, dtype=np.float64)
    if radii.size == 1:  # a scalar, or one radius for every query point as cKDTree allows
        return radii.reshape(1)
    return np.ascontiguousarray(np.broadcast_to(radii, lead).reshape(-1))


def _exact_order_rows(index: "SpatialIndex", rows: np.ndarray, dist: np.ndarray, idx: np.ndarray):
    """Re-order candidate rows by `(squared distance, index)`, misses last."""
    found = np.isfinite(dist)
    squared = np.full(dist.shape, np.inf)
    r, c = np.nonzero(found)
    squared[r, c] = _squared_distances(rows[r], index.data[idx[r, c]], bool(FUSED))
    order = np.lexsort((idx, squared), axis=1)
    return np.take_along_axis(dist, order, axis=1), np.take_along_axis(idx, order, axis=1)


def _tie_reaches_end(dist: np.ndarray, k: int) -> np.ndarray:
    """Rows whose last returned neighbour is a finite tie with the `k`-th: more tied points may
    lie beyond what was asked for."""
    return np.isfinite(dist[:, k - 1]) & (dist[:, -1] == dist[:, k - 1])


def _rows_to_lists(indptr: np.ndarray, indices: np.ndarray, lead: tuple[int, ...]):
    flat = indices.tolist()
    bounds = indptr.tolist()
    rows = [flat[a:b] for a, b in zip(bounds[:-1], bounds[1:])]
    out = np.empty(len(rows), dtype=object)
    out[:] = rows
    return out.reshape(lead)


class SpatialIndex:
    """Immutable index over `(n, 3)` points; see the module docstring for its semantics."""

    __slots__ = ("data", "n", "m", "_native", "_tree")

    def __init__(self, points) -> None:
        # An owned copy: `cKDTree` keeps a view of its caller's array, so a caller editing that
        # array afterwards would change `data` under a tree built from the old coordinates.
        data = np.array(points, dtype=np.float64, order="C", copy=True)
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"SpatialIndex needs (n, 3) points, got shape {data.shape}")
        data.flags.writeable = False
        self.data = data
        self.n = int(data.shape[0])
        self.m = 3
        if _KdIndex is not None:
            self._native = _KdIndex(data, FUSED)
            self._tree = None
        else:
            self._native = None
            self._tree = cKDTree(data)

    @property
    def backend(self) -> str:
        return "rust" if self._native is not None else "scipy"

    def __reduce__(self):
        return (SpatialIndex, (np.array(self.data),))

    def __copy__(self) -> "SpatialIndex":
        return self

    def __deepcopy__(self, memo: dict) -> "SpatialIndex":
        # Immutable, so a deep copy may share it -- as copies of a cached world do.
        return self

    def query(self, x, k=1, eps=0, p=2, distance_upper_bound=np.inf, workers=1):
        """`cKDTree.query`: distances and indices of the `k` nearest points (`k` an integer, or a
        sequence of 1-based ranks as in scipy), exact ties to the lowest index."""
        _check_supported(p, eps)
        if np.ndim(k) == 0:
            if float(k) != int(k):
                raise ValueError("k must be a positive integer or a sequence of positive ranks")
            k = int(k)  # cKDTree accepts integral floats too
        ranks = None if isinstance(k, int) else [int(r) for r in k]
        kmax = int(k) if ranks is None else max(ranks, default=0)
        if kmax < 1 or (ranks is not None and min(ranks) < 1):
            raise ValueError("k must be a positive integer or a sequence of positive ranks")
        rows, lead, single = _as_queries(x)
        if self._native is not None:
            dist, idx = self._native.query(rows, kmax, float(distance_upper_bound))
        else:
            dist, idx = self._scipy_query(rows, kmax, distance_upper_bound, workers)
        if ranks is not None:
            cols = np.asarray(ranks) - 1
            shape = lead + (len(ranks),)
            return dist.reshape(-1, kmax)[:, cols].reshape(shape), idx.reshape(-1, kmax)[:, cols].reshape(shape)
        if kmax == 1:
            if single:
                return float(dist[0]), int(idx[0])
            return dist.reshape(lead), idx.reshape(lead)
        return dist.reshape(lead + (kmax,)), idx.reshape(lead + (kmax,))

    def _scipy_query(self, rows: np.ndarray, k: int, bound: float, workers: int) -> tuple[np.ndarray, np.ndarray]:
        """The Rust backend's `query` answer from `cKDTree`: ask for one more neighbour than
        needed, and more again for rows where an equal distance still reaches the last slot
        asked for, so every point tied with the `k`-th is among the candidates. `cKDTree` already
        orders each row by squared distance; rows where two returned distances are equal (an
        exact tie, or distinct squared distances that round to the same root) are re-ordered by
        `(squared distance, index)` before keeping `k`. Distances are `cKDTree`'s own."""
        m = len(rows)
        if m == 0:
            return np.zeros(0), np.zeros(0, dtype=np.int64)
        extra = min(k + 1, self.n) if self.n else 1
        dist, idx = (a.reshape(m, -1) for a in self._tree.query(rows, k=extra, distance_upper_bound=bound, workers=workers))
        open_rows = np.flatnonzero(_tie_reaches_end(dist, k)) if extra > k else np.zeros(0, dtype=np.intp)
        while len(open_rows) and extra < self.n:
            extra = min(2 * extra, self.n)
            more_d, more_i = (a.reshape(len(open_rows), -1) for a in self._tree.query(rows[open_rows], k=extra, distance_upper_bound=bound, workers=workers))
            dist = np.pad(dist, ((0, 0), (0, extra - dist.shape[1])), constant_values=np.inf)
            idx = np.pad(idx, ((0, 0), (0, extra - idx.shape[1])), constant_values=self.n)
            dist[open_rows], idx[open_rows] = more_d, more_i
            open_rows = open_rows[_tie_reaches_end(more_d, k)]
        collide = np.flatnonzero(np.any((dist[:, 1:] == dist[:, :-1]) & np.isfinite(dist[:, 1:]), axis=1))
        if len(collide):
            dist[collide], idx[collide] = self._exact_order(rows[collide], dist[collide], idx[collide])
        dist, idx = dist[:, :k], idx[:, :k]
        if dist.shape[1] < k:  # fewer points than k: pad as cKDTree does
            dist = np.pad(dist, ((0, 0), (0, k - dist.shape[1])), constant_values=np.inf)
            idx = np.pad(idx, ((0, 0), (0, k - idx.shape[1])), constant_values=self.n)
        return dist.reshape(-1), idx.reshape(-1).astype(np.int64, copy=False)

    def _exact_order(self, rows, dist, idx):
        return _exact_order_rows(self, rows, dist, idx)

    def query_ball_point_csr(self, x, r) -> tuple[np.ndarray, np.ndarray]:
        """Every point within `r` (scalar, or one per query point) of each row of `x`, as int64
        CSR arrays: row `i`'s ascending indices are `indices[indptr[i]:indptr[i + 1]]`."""
        rows, lead, _single = _as_queries(x)
        radii = _radii(r, lead)
        if self._native is not None:
            return self._native.query_ball_point(rows, radii)
        lists = self._tree.query_ball_point(rows, radii if radii.size > 1 else float(radii[0]), return_sorted=True)
        counts = np.fromiter(map(len, lists), dtype=np.int64, count=len(lists))
        indptr = np.zeros(len(lists) + 1, dtype=np.int64)
        np.cumsum(counts, out=indptr[1:])
        indices = np.fromiter(itertools.chain.from_iterable(lists), dtype=np.int64, count=int(indptr[-1]))
        return indptr, indices

    def query_ball_point(self, x, r, p=2.0, eps=0, workers=1, return_length=False):
        """`cKDTree.query_ball_point` (without `return_sorted`: rows always come back ascending,
        see the module docstring)."""
        _check_supported(p, eps)
        if self._native is None:
            return self._tree.query_ball_point(x, r, workers=workers, return_sorted=True, return_length=return_length)
        rows, lead, single = _as_queries(x)
        radii = _radii(r, lead)
        if return_length:
            counts = self._native.query_ball_point_count(rows, radii)
            return counts[0] if single else counts.reshape(lead)
        indptr, indices = self._native.query_ball_point(rows, radii)
        if single:
            return indices.tolist()
        return _rows_to_lists(indptr, indices, lead)

    def query_pairs(self, r, p=2.0, eps=0, output_type="set"):
        """`cKDTree.query_pairs`; the `"ndarray"` form is in lexicographic order."""
        _check_supported(p, eps)
        if output_type not in ("set", "ndarray"):
            raise ValueError("output_type must be 'set' or 'ndarray'")
        if self._native is not None:
            pairs = self._native.query_pairs(float(r)).reshape(-1, 2)
        elif output_type == "set":
            return self._tree.query_pairs(r)
        else:
            pairs = self._tree.query_pairs(r, output_type="ndarray")
            pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))] if len(pairs) else pairs.reshape(0, 2)
        if output_type == "set":
            return set(map(tuple, pairs.tolist()))
        return pairs
