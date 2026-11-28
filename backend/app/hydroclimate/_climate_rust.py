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


def can_sample_offset(field, dir_u, dir_v, distance, lat):
    """Keep NumPy's casting behavior for unusual or non-finite inputs."""
    return (
        has_kernel("offset_indices_rs") and field.ndim == 2 and field.size > 0
        and dir_u.shape == dir_v.shape == field.shape
        and dir_u.dtype == dir_v.dtype == lat.dtype == np.float64
        and lat.shape == (field.shape[0],) and np.isfinite(distance) and abs(distance) <= 3600
        and np.all(np.isfinite(lat)) and np.all(np.abs(lat) <= 3600)
        and np.all(np.isfinite(dir_u)) and np.all(np.abs(dir_u) <= 1000)
        and np.all(np.isfinite(dir_v)) and np.all(np.abs(dir_v) <= 1000)
    )


def sample_at_offset(field, dir_u, dir_v, distance, lat):
    """Gather nearest cells after computing offset indices in Rust."""
    height, width = field.shape
    lon = -180.0 + (np.arange(width) + 0.5) * (360.0 / width)
    cos_lat = np.clip(np.cos(np.radians(lat)), 0.15, 1.0)
    rows, cols = _rust.offset_indices_rs(
        np.ascontiguousarray(dir_u), np.ascontiguousarray(dir_v),
        np.ascontiguousarray(lat), np.ascontiguousarray(cos_lat), lon, float(distance),
    )
    return field[rows, cols]


def supports_float_grid(*fields):
    """Accept matching float64 grids and leave other NumPy dtypes on the fallback."""
    return (fields[0].ndim == 2 and fields[0].size > 0
            and all(a.shape == fields[0].shape and a.dtype == np.float64 for a in fields))


def moisture_divergence(humidity, u, v, cos_lat, radius):
    """Compute metric-correct centered moisture divergence."""
    return _rust.moisture_divergence_rs(
        *(np.ascontiguousarray(a) for a in (humidity, u, v, cos_lat)), float(radius),
    )
