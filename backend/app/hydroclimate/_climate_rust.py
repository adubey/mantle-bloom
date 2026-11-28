"""Optional native diagnostic climate kernels (#1046)."""

import os

import numpy as np

try:
    if os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1":
        raise ImportError
    import mantle_bloom_rust as _rust
except ImportError:
    _rust = None


def available():
    """Allow old extension builds to fall back until rebuilt."""
    return has_kernel("humidity_zonal_sweep_rs")


def humidity_zonal_sweep(ocean, elevation, ceiling, factor, source, direction, retention, *constants):
    """Pass contiguous arrays to the native sweep."""
    return _rust.humidity_zonal_sweep_rs(
        np.ascontiguousarray(ocean, dtype=bool),
        *(np.ascontiguousarray(a, dtype=np.float64) for a in
          (elevation, ceiling, factor, source, direction, retention)),
        *constants,
    )


def has_kernel(name):
    """Check each export independently for compatibility with older builds."""
    return _rust is not None and hasattr(_rust, name)


def smooth_along_coast(u, v, ocean, lat, ambient_u, ambient_v, *constants):
    """Run the coastal Jacobi iterations on contiguous arrays."""
    return _rust.smooth_along_coast_rs(
        np.ascontiguousarray(u, dtype=np.float64), np.ascontiguousarray(v, dtype=np.float64),
        np.ascontiguousarray(ocean, dtype=bool),
        *(np.ascontiguousarray(a, dtype=np.float64) for a in (lat, ambient_u, ambient_v)),
        *constants,
    )


def supports_grid(ocean, *fields):
    """Keep NumPy semantics for nonstandard dtypes and shapes."""
    return (
        ocean.ndim == 2 and ocean.dtype == np.bool_ and ocean.size > 0
        and all(a.shape == ocean.shape and a.dtype == np.float64 for a in fields)
    )
