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


def route_downstream(order, stuck, target, ocean, source, retain, loss, erosion):
    """Accumulate flux using the Python routing order."""
    def i64(a):
        return np.ascontiguousarray(a, dtype=np.int64)

    def f64(a):
        return np.ascontiguousarray(a, dtype=np.float64)

    return _rust.route_downstream_rs(
        i64(order), i64(stuck), i64(target), np.ascontiguousarray(ocean, dtype=bool),
        f64(source), f64(retain), f64(loss), None if erosion is None else f64(erosion),
    )


def routing_order(elevation, ocean, target):
    """Order an uphill flow graph using the native priority queue."""
    return _rust.routing_order_rs(
        np.ascontiguousarray(elevation, dtype=np.float64),
        np.ascontiguousarray(ocean, dtype=bool), np.ascontiguousarray(target, dtype=np.int64),
    )
