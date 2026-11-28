"""Shared 3D point index with `cKDTree`'s query semantics, backed by Rust when available (#1042).

`SpatialIndex(points)` answers the `scipy.spatial.cKDTree` queries the simulation and renderer
use -- `query` (k nearest, optionally bounded), `query_ball_point` and `query_pairs` -- with the
same results, and adds `query_ball_point_csr`, which returns one flat CSR `(indptr, indices)` pair
instead of a Python list per query point. With `mantle_bloom_rust` installed it is the native
`KdIndex` (`backend/rust/src/spatial/`), whose batch queries run across every core without the
GIL; without it, or with `MANTLE_BLOOM_DISABLE_RUST=1`, it wraps a `cKDTree`.

What both backends guarantee, so either can stand in for the other or for a bare `cKDTree`:

* Distances are `cKDTree`'s, bit for bit: the square root of `dx² + dy² + dz²` summed in that
  order with the installed scipy's own rounding, which `_probe_fused` reads once at import.
  Every threshold test therefore agrees too: `query` keeps points strictly closer than
  `distance_upper_bound` (misses are `inf` with index `n`), and radius queries keep points at
  most `r` away.
* `query` rows are sorted by distance; `query_ball_point` rows are ascending, single-point
  queries included (scipy leaves those in traversal order); `query_pairs(output_type="ndarray")`
  is in lexicographic order (scipy's is traversal order).
* Exact distance ties go to the lowest index in the Rust backend. `cKDTree` keeps whichever tied
  point its traversal meets first, which depends on its leaf size and build options. The two
  agree on every query without an exact tie; on the node clouds the simulation queries, exact
  ties need coincident points.

Only Euclidean (`p=2`), exact (`eps=0`) queries on `(n, 3)` points are supported. `workers` is
accepted for drop-in compatibility; the Rust backend picks its own parallelism from the batch
size, which never changes a result.
"""

from __future__ import annotations

import itertools
import logging
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


def _probe_fused(kdindex) -> bool | None:
    """Whether the installed scipy fuses its squared-distance sum (see `kdtree::Metric` in the
    Rust crate), read by comparing `cKDTree` distances against both native roundings on a fixed
    random cloud. About one query in nine rounds differently under the two, so 512 queries
    always tell them apart. `None` (fall back to scipy) if neither rounding reproduces scipy."""
    rng = np.random.default_rng(1042)
    points = rng.normal(size=(64, 3))
    queries = rng.normal(size=(512, 3))
    expected, _ = cKDTree(points).query(queries)
    for fused in (True, False):
        got, _ = kdindex(points, fused).query(queries, 1, np.inf)
        if np.array_equal(got, expected):
            return fused
    _logger.warning("mantle_bloom_rust KdIndex distances differ from scipy's; spatial queries fall back to scipy.")
    return None


_KdIndex = _import_kdindex()
FUSED: bool | None = _probe_fused(_KdIndex) if _KdIndex is not None else None
if FUSED is None:
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
    if radii.ndim == 0:
        return radii.reshape(1)
    return np.ascontiguousarray(np.broadcast_to(radii, lead).reshape(-1))


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
        data = np.ascontiguousarray(points, dtype=np.float64)
        if data.ndim != 2 or data.shape[1] != 3:
            raise ValueError(f"SpatialIndex needs (n, 3) points, got shape {data.shape}")
        data = data.view()
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
        """`cKDTree.query` for an integer `k`: distances and indices of the `k` nearest points."""
        _check_supported(p, eps)
        if not isinstance(k, (int, np.integer)) or k < 1:
            raise ValueError("k must be a positive integer")
        if self._native is None:
            return self._tree.query(x, k=k, distance_upper_bound=distance_upper_bound, workers=workers)
        rows, lead, single = _as_queries(x)
        dist, idx = self._native.query(rows, int(k), float(distance_upper_bound))
        if k == 1:
            if single:
                return float(dist[0]), int(idx[0])
            return dist.reshape(lead), idx.reshape(lead)
        return dist.reshape(lead + (int(k),)), idx.reshape(lead + (int(k),))

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
