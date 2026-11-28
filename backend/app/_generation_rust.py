"""Optional Rust world-generation kernels (#1042)."""

import logging
import os

import numpy as np

try:
    if os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1":
        raise ImportError
    import mantle_bloom_rust as _rust
except ImportError:
    _rust = None

_REQUIRED_RUST_EXPORTS = ("priority_flood_rs", "nearest_distance_rs", "label_component_count_rs")

# Brute-force nearest in Rust beats a `cKDTree` build plus query only while the
# source set is small (measured crossover near 500 sources at lattice scale,
# #1042); larger source sets keep the tree.
NEAREST_RUST_MAX_SOURCES = 512


def _warn_missing_exports(module):
    if module is not None:
        missing = [name for name in _REQUIRED_RUST_EXPORTS if not hasattr(module, name)]
        if missing:
            logging.getLogger(__name__).warning(
                "mantle_bloom_rust is missing generation exports %s (stale build?); affected kernels "
                "fall back to Python. Rebuild with `backend/.venv/bin/python bin/build_rust.py`.",
                missing,
            )


_warn_missing_exports(_rust)
RUST_AVAILABLE = _rust is not None and hasattr(_rust, "priority_flood_rs")


def has_kernel(name):
    """Return whether an individual native kernel is available."""
    return RUST_AVAILABLE and hasattr(_rust, name)


def use_nearest(sources) -> bool:
    """Whether the Rust brute-force nearest should replace a `cKDTree` query."""
    return has_kernel("nearest_distance_rs") and len(sources) <= NEAREST_RUST_MAX_SOURCES


def _f64(a):
    return np.ascontiguousarray(a, dtype=np.float64)


def _i64(a):
    return np.ascontiguousarray(a, dtype=np.int64)


def _xyz(a):
    return np.ascontiguousarray(a, dtype=np.float64).reshape(-1, 3)


def priority_flood(indptr, indices, rise, seeds_level, elevation):
    """Lowest drainage level per node; `rise` may be None for a level flood."""
    return _rust.priority_flood_rs(
        _i64(indptr), _i64(indices), None if rise is None else _f64(rise), _f64(seeds_level), _f64(elevation)
    )


def nearest_distance(sources, points):
    """(distance km, source index) of the nearest source per point."""
    dist, index = _rust.nearest_distance_rs(_xyz(sources), _xyz(points))
    return dist, index


def label_component_count(offsets, neighbours, labels):
    """Connected components of the lattice graph with links cut between labels."""
    return int(_rust.label_component_count_rs(_i64(offsets), _i64(neighbours), _i64(labels)))
