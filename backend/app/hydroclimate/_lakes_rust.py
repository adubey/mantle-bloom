"""Thin Python bridge to `mantle_bloom_rust` lake and component kernels (#1046)."""

from __future__ import annotations

import logging
import os

import numpy as np

_logger = logging.getLogger(__name__)
_DISABLE_RUST = os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1"

_REQUIRED_RUST_EXPORTS = (
    "lake_components_rs",
    "catchment_roots_rs",
    "lake_hierarchy_events_rs",
)


def _resolve_rust_module(mod):
    if mod is None:
        return None
    missing = [name for name in _REQUIRED_RUST_EXPORTS if not hasattr(mod, name)]
    if missing:
        _logger.warning(
            "mantle_bloom_rust is missing lake exports %s (stale build?); falling back to Python. "
            "Rebuild with `backend/.venv/bin/python bin/build_rust.py` or `./bin/setup.sh`.",
            missing,
        )
        return None
    return mod


if _DISABLE_RUST:
    _imported_rust = None
else:
    try:
        import mantle_bloom_rust as _imported_rust
    except ImportError:
        _imported_rust = None

_rust = _resolve_rust_module(_imported_rust)
RUST_AVAILABLE: bool = _rust is not None


def _f64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.float64)


def _i64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.int64)


def _bool(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.bool_)


def lake_components(is_lake: np.ndarray, neighbor_idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Connected lake-node groups as CSR `(indptr, members)` via Rust."""
    return _rust.lake_components_rs(_bool(is_lake), _i64(neighbor_idx))


def catchment_roots(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    ocean_catchment: int,
) -> np.ndarray:
    """Pure steepest-descent catchment roots via Rust."""
    return _rust.catchment_roots_rs(_f64(elevation), _bool(is_ocean), _i64(neighbor_idx), int(ocean_catchment))


def lake_hierarchy_events(
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    neighbor_idx: np.ndarray,
    interface_pass_elevation: np.ndarray | None = None,
    ocean_catchment: int = -2,
):
    """Catchment roots and ordered Kruskal merge events via Rust."""
    pass_arr = None if interface_pass_elevation is None else _f64(interface_pass_elevation)
    return _rust.lake_hierarchy_events_rs(
        _f64(elevation),
        _bool(is_ocean),
        _i64(neighbor_idx),
        pass_arr,
        int(ocean_catchment),
    )
