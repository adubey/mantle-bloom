"""Optional Rust flow-routing kernels (#1041)."""

import os

import numpy as np

try:
    if os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1":
        raise ImportError
    import mantle_bloom_rust as _rust
except ImportError:
    _rust = None

RUST_AVAILABLE = _rust is not None and hasattr(_rust, "flow_direction_rs")


def flow_direction(elevation, is_ocean, neighbor_idx, offset, filled, spill, channel, threshold):
    """Copy inputs to contiguous arrays before invoking the native kernel."""
    def f64(a):
        return np.ascontiguousarray(a, dtype=np.float64)

    def i64(a):
        return np.ascontiguousarray(a, dtype=np.int64)
    return _rust.flow_direction_rs(
        f64(elevation), np.ascontiguousarray(is_ocean, dtype=bool), i64(neighbor_idx),
        f64(offset), f64(filled), i64(spill), f64(channel), float(threshold),
    )
