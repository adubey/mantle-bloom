"""Optional short-trace distance kernel for fault deformation (#1046)."""

import os

import numpy as np

try:
    if os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1":
        raise ImportError
    import mantle_bloom_rust as _rust
except ImportError:
    _rust = None


def nearest_trace(points, trace):
    """Return distances and first nearest trace indices, with a NumPy fallback."""
    if (
        _rust is not None and hasattr(_rust, "nearest_trace_rs")
        and points.dtype == trace.dtype == np.float64
    ):
        return _rust.nearest_trace_rs(np.ascontiguousarray(points), np.ascontiguousarray(trace))
    distances = np.linalg.norm(points[:, None, :] - trace[None, :, :], axis=-1)
    indices = np.argmin(distances, axis=1)
    return distances[np.arange(len(points)), indices], indices


def shear_backtrace(points, trace, dip, reach, slip, sense):
    """Return moved point indices and their normalized upstream positions."""
    if (_rust is not None and hasattr(_rust, "shear_backtrace_rs")
            and points.dtype == trace.dtype == dip.dtype == np.float64
            and np.isfinite(points).all() and np.isfinite(trace).all()
            and np.isfinite(dip).all() and np.isfinite(reach) and reach > 0 and np.isfinite(slip)):
        return _rust.shear_backtrace_rs(
            np.ascontiguousarray(points), np.ascontiguousarray(trace),
            np.ascontiguousarray(dip), reach, slip, float(sense),
        )
    d, nearest = nearest_trace(points, trace)
    taper = np.clip(1.0 - d / reach, 0.0, 1.0)
    lo = np.clip(nearest - 1, 0, len(trace) - 1)
    hi = np.clip(nearest + 1, 0, len(trace) - 1)
    tangent = trace[hi] - trace[lo]
    tangent = tangent - np.sum(tangent * points, axis=-1, keepdims=True) * points
    norm = np.linalg.norm(tangent, axis=-1, keepdims=True)
    valid = norm[:, 0] > 1e-12
    tangent = np.divide(tangent, norm, out=np.zeros_like(tangent), where=norm > 1e-12)
    side = np.sign((points - trace[len(trace) // 2]) @ dip)
    side[side == 0] = 1.0
    shift = taper * side * sense * slip
    shift[~valid] = 0.0
    moved = np.flatnonzero(np.abs(shift) > 1e-12)
    upstream = points[moved] - shift[moved, None] * tangent[moved]
    norm = np.linalg.norm(upstream, axis=-1, keepdims=True)
    return moved, upstream / np.where(norm < 1e-12, 1.0, norm)


def influence(distances, reach, floor):
    """Map nearest distances to fault-deformation weights."""
    if _rust is not None and hasattr(_rust, "influence_rs") and distances.dtype == np.float64 and reach != 0:
        return _rust.influence_rs(np.ascontiguousarray(distances), reach, floor)
    return np.clip(1.0 - distances / reach, floor, 1.0)


def foreign_nearest(points, outline, owners, indptr, indices, owner):
    """Find the first closest foreign point in each sorted neighborhood."""
    if _rust is not None and hasattr(_rust, "foreign_nearest_rs"):
        return _rust.foreign_nearest_rs(
            np.ascontiguousarray(points, dtype=np.float64), np.ascontiguousarray(outline, dtype=np.float64),
            np.ascontiguousarray(owners, dtype=np.int64), np.ascontiguousarray(indptr, dtype=np.int64),
            np.ascontiguousarray(indices, dtype=np.int64), int(owner),
        )
    distances = np.full(len(outline), np.inf)
    best = np.full(len(outline), -1, dtype=np.int64)
    for i in range(len(outline)):
        candidates = indices[indptr[i]:indptr[i + 1]]
        foreign = candidates[owners[candidates] != owner]
        if len(foreign):
            d = np.linalg.norm(points[foreign] - outline[i], axis=1)
            best[i] = foreign[np.argmin(d)]
            distances[i] = np.min(d)
    return distances, best


def relief_contrib(points, distances, mid, dip, params):
    """Calculate requested elevation changes before conservative redistribution."""
    if (_rust is not None and hasattr(_rust, "relief_contrib_rs")
            and all(a.dtype == np.float64 and np.isfinite(a).all() for a in (points, distances, mid, dip))
            and np.isfinite(params).all() and params[0] > 0):
        return _rust.relief_contrib_rs(
            np.ascontiguousarray(points), np.ascontiguousarray(distances), np.ascontiguousarray(mid),
            np.ascontiguousarray(dip), np.asarray(params, dtype=np.float64),
        )
    reach, slip, years, rate, scale, positive, negative, bend = params
    mag = np.clip(1.0 - distances / reach, 0.0, 1.0) * slip * years * rate * scale
    side = (points - mid) @ dip > 0.0
    return np.where(side, positive, negative) * bend * mag


def relief_transfer(hc, material, areas, kinds, requested, minimum, maximum):
    """Redistribute crust and its tracer within crust types, conserving volume."""
    if (_rust is not None and hasattr(_rust, "relief_transfer_rs")
            and kinds.dtype == np.int8
            and all(a.dtype == np.float64 and np.isfinite(a).all() for a in (hc, material, areas, requested))
            and np.isfinite(minimum) and np.isfinite(maximum) and 0 < minimum <= maximum
            and np.all(areas >= 0)):
        return _rust.relief_transfer_rs(
            *(np.ascontiguousarray(a) for a in (hc, material, areas, kinds, requested)), minimum, maximum,
        )
    hc, material = hc.copy(), material.copy()
    for kind in np.unique(kinds):
        group = kinds == kind
        donor = group & (requested < 0) & (hc > minimum)
        receiver = group & (requested > 0) & (hc < maximum)
        debit = np.minimum(-requested[donor], hc[donor] - minimum)
        credit = np.minimum(requested[receiver], maximum - hc[receiver])
        debit_volume = float(debit @ areas[donor])
        credit_volume = float(credit @ areas[receiver])
        amount = min(debit_volume, credit_volume)
        if amount <= 0:
            continue
        taken = debit * (amount / debit_volume)
        given = credit * (amount / credit_volume)
        moved = taken / hc[donor] * material[donor]
        material[donor] -= moved
        material[receiver] += float(moved @ areas[donor]) * given / amount
        hc[donor] -= taken
        hc[receiver] += given
    return hc, material


def boundary_runs(regime, owners):
    """Identify consecutive classified stretches sharing the same neighbor."""
    if _rust is not None and hasattr(_rust, "boundary_runs_rs"):
        return _rust.boundary_runs_rs(np.ascontiguousarray(regime, dtype=np.int8),
                                      np.ascontiguousarray(owners, dtype=np.int64))
    starts, ends = [], []
    i = 0
    while i < len(regime):
        if regime[i] == 0:
            i += 1
            continue
        end = i + 1
        while end < len(regime) and regime[end] == regime[i] and owners[end] == owners[i]:
            end += 1
        starts.append(i)
        ends.append(end)
        i = end
    return np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)


def unloading_stress(load, sense, ratio, sine, cosine, friction):
    """Resolve load changes to Coulomb stress, preserving scalar and grid fallback."""
    if (_rust is not None and hasattr(_rust, "unloading_stress_rs")
            and load.ndim == 1 and load.dtype == np.float64):
        return _rust.unloading_stress_rs(np.ascontiguousarray(load), sense, ratio, sine, cosine, friction)
    shear = sense * (1.0 - ratio) * load * sine * cosine
    normal = load * (cosine * cosine + ratio * sine * sine)
    return shear - friction * normal


def quake_magnitude(length, slip, overlap, jitter, constants):
    """Calculate magnitude while retaining NumPy's seeded random draws in Python."""
    if (_rust is not None and hasattr(_rust, "quake_magnitude_rs")
            and np.isfinite((length, slip, jitter, *constants)).all() and constants[4] <= constants[5]):
        return _rust.quake_magnitude_rs(length, slip, overlap, jitter, np.asarray(constants))
    base, lc, sc, bonus, lower, upper = constants
    value = (base + lc * np.log10(max(length, 1.0)) + sc * np.log10(max(slip, 1.0))
             + (bonus if overlap else 0.0) + jitter)
    return float(np.clip(value, lower, upper))


def boundary_regime(distances, closing, reach, threshold):
    """Classify nearby boundary points by signed closing rate."""
    if (_rust is not None and hasattr(_rust, "boundary_regime_rs")
            and distances.dtype == closing.dtype == np.float64
            and np.isfinite(reach) and np.isfinite(threshold) and threshold >= 0):
        return _rust.boundary_regime_rs(np.ascontiguousarray(distances), np.ascontiguousarray(closing), reach, threshold)
    near = distances <= reach
    regime = np.zeros(len(distances), dtype=np.int8)
    regime[near & (closing > threshold)] = 1
    regime[near & (closing < -threshold)] = 2
    regime[near & (np.abs(closing) <= threshold)] = 3
    return regime


def moment_magnitude(length, slip, width, shear, minimum, maximum):
    """Calculate unloading-rupture moment magnitude."""
    if (_rust is not None and hasattr(_rust, "moment_magnitude_rs")
            and np.isfinite((length, slip, width, shear, minimum, maximum)).all() and minimum <= maximum):
        return _rust.moment_magnitude_rs(length, slip, width, shear, minimum, maximum)
    moment = max(shear * length * min(width, length) * slip, 1.0)
    return float(np.clip((2.0 / 3.0) * (np.log10(moment) - 9.1), minimum, maximum))


def resample_boundary(points, positions):
    """Interpolate a boundary polyline and project onto the unit sphere."""
    if (_rust is not None and hasattr(_rust, "resample_boundary_rs")
            and points.dtype == positions.dtype == np.float64 and np.isfinite(points).all()):
        return _rust.resample_boundary_rs(np.ascontiguousarray(points), np.ascontiguousarray(positions))
    lo = np.floor(positions).astype(int)
    hi = np.minimum(lo + 1, len(points) - 1)
    frac = (positions - lo)[:, None]
    out = points[lo] * (1.0 - frac) + points[hi] * frac
    norm = np.linalg.norm(out, axis=-1, keepdims=True)
    return out / np.where(norm < 1e-12, 1.0, norm)


def boundary_strands(master, centroid, count, spacing):
    """Construct tangent-plane inward vectors and normalized offset traces."""
    if (_rust is not None and hasattr(_rust, "boundary_strands_rs")
            and master.dtype == centroid.dtype == np.float64 and count >= 0):
        return _rust.boundary_strands_rs(np.ascontiguousarray(master), np.ascontiguousarray(centroid), count, spacing)
    offset = centroid[None, :] - master
    inward = offset - np.sum(offset * master, axis=-1, keepdims=True) * master
    norm = np.linalg.norm(inward, axis=-1, keepdims=True)
    inward = inward / np.where(norm < 1e-12, 1.0, norm)
    strands = []
    for k in range(count):
        strand = master + (k + 1) * spacing * inward
        norm = np.linalg.norm(strand, axis=-1, keepdims=True)
        strands.append(strand / np.where(norm < 1e-12, 1.0, norm))
    return inward, np.asarray(strands).reshape(count, len(master), 3)
