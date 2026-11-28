//! Short-trace distance calculations for fault deformation (#1046).

use numpy::{PyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Match NumPy's first-minimum tie breaking for each point against a fault trace.
#[pyfunction]
pub fn nearest_trace_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    trace: PyReadonlyArray2<'py, f64>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i64>>)> {
    if points.shape()[1] != 3 || trace.shape()[1] != 3 || trace.shape()[0] == 0 {
        return Err(PyValueError::new_err(
            "points and nonempty trace must have shape (n, 3)",
        ));
    }
    let points = points.as_slice()?;
    let trace = trace.as_slice()?;
    let mut distances = Vec::with_capacity(points.len() / 3);
    let mut indices = Vec::with_capacity(points.len() / 3);
    for point in points.chunks_exact(3) {
        let mut best = f64::INFINITY;
        let mut index = 0;
        for (i, node) in trace.chunks_exact(3).enumerate() {
            let x = point[0] - node[0];
            let y = point[1] - node[1];
            let z = point[2] - node[2];
            let d = (x * x + y * y + z * z).sqrt();
            if d < best || (d.is_nan() && !best.is_nan()) {
                best = d;
                index = i;
            }
        }
        distances.push(best);
        indices.push(index as i64);
    }
    Ok((
        PyArray1::from_vec(py, distances),
        PyArray1::from_vec(py, indices),
    ))
}

/// Compute semi-Lagrangian backtraces in the same affected-node order as NumPy.
#[pyfunction]
pub fn shear_backtrace_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    trace: PyReadonlyArray2<'py, f64>,
    dip: numpy::PyReadonlyArray1<'py, f64>,
    reach: f64,
    slip: f64,
    sense: f64,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, numpy::PyArray2<f64>>)> {
    use numpy::IntoPyArray;
    if points.shape()[1] != 3 || trace.shape()[1] != 3 || trace.shape()[0] == 0 || dip.len() != 3 {
        return Err(PyValueError::new_err("invalid shear geometry"));
    }
    let pts = points.as_slice()?;
    let line = trace.as_slice()?;
    let dip = dip.as_slice()?;
    let n = line.len() / 3;
    let mid = &line[3 * (n / 2)..][..3];
    let mut indices = Vec::new();
    let mut upstream = Vec::new();
    for (i, p) in pts.chunks_exact(3).enumerate() {
        let mut best = f64::INFINITY;
        let mut nearest = 0;
        for (j, q) in line.chunks_exact(3).enumerate() {
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            if d < best {
                best = d;
                nearest = j;
            }
        }
        let lo = nearest.saturating_sub(1);
        let hi = (nearest + 1).min(n - 1);
        let mut tangent = [0.0; 3];
        for k in 0..3 {
            tangent[k] = line[3 * hi + k] - line[3 * lo + k];
        }
        let dot = tangent[0] * p[0] + tangent[1] * p[1] + tangent[2] * p[2];
        for k in 0..3 {
            tangent[k] -= dot * p[k];
        }
        let norm =
            (tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]).sqrt();
        if norm <= 1e-12 {
            continue;
        }
        let side = (p[0] - mid[0]) * dip[0] + (p[1] - mid[1]) * dip[1] + (p[2] - mid[2]) * dip[2];
        let shift = (1.0 - best / reach).clamp(0.0, 1.0)
            * if side < 0.0 { -1.0 } else { 1.0 }
            * sense
            * slip;
        if shift.abs() <= 1e-12 {
            continue;
        }
        let mut u = [0.0; 3];
        for k in 0..3 {
            u[k] = p[k] - shift * (tangent[k] / norm);
        }
        let length = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let length = if length < 1e-12 { 1.0 } else { length };
        indices.push(i as i64);
        upstream.extend(u.iter().map(|v| v / length));
    }
    let array = numpy::ndarray::Array2::from_shape_vec((indices.len(), 3), upstream).unwrap();
    Ok((PyArray1::from_vec(py, indices), array.into_pyarray(py)))
}

/// Turn nearest-trace chord distances into deformation weights.
#[pyfunction]
pub fn influence_rs<'py>(
    py: Python<'py>,
    distances: numpy::PyReadonlyArray1<'py, f64>,
    reach: f64,
    floor: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let values = distances
        .as_slice()?
        .iter()
        .map(|&d| {
            let x = 1.0 - d / reach;
            if x.is_nan() || floor.is_nan() {
                f64::NAN
            } else {
                x.max(floor).min(1.0)
            }
        })
        .collect();
    Ok(PyArray1::from_vec(py, values))
}

/// Select the first nearest foreign-owned node in each sorted CSR neighborhood.
#[pyfunction]
pub fn foreign_nearest_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    outline: PyReadonlyArray2<'py, f64>,
    owners: numpy::PyReadonlyArray1<'py, i64>,
    indptr: numpy::PyReadonlyArray1<'py, i64>,
    indices: numpy::PyReadonlyArray1<'py, i64>,
    owner: i64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i64>>)> {
    let n = points.shape()[0];
    let m = outline.shape()[0];
    if points.shape()[1] != 3
        || outline.shape()[1] != 3
        || owners.len() != n
        || indptr.len() != m + 1
    {
        return Err(PyValueError::new_err("invalid foreign-neighbor shapes"));
    }
    let pts = points.as_slice()?;
    let out = outline.as_slice()?;
    let owners = owners.as_slice()?;
    let offsets = indptr.as_slice()?;
    let ids = indices.as_slice()?;
    if offsets.first() != Some(&0)
        || offsets.last() != Some(&(ids.len() as i64))
        || offsets.windows(2).any(|w| w[0] < 0 || w[0] > w[1])
        || ids.iter().any(|&id| id < 0 || id as usize >= n)
    {
        return Err(PyValueError::new_err("invalid foreign-neighbor CSR"));
    }
    let mut distances = vec![f64::INFINITY; m];
    let mut best_ids = vec![-1_i64; m];
    for row in 0..m {
        let p = &out[3 * row..][..3];
        for &id in &ids[offsets[row] as usize..offsets[row + 1] as usize] {
            let j = id as usize;
            if owners[j] == owner {
                continue;
            }
            let q = &pts[3 * j..][..3];
            let d = ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt();
            if d < distances[row] || best_ids[row] == -1 {
                distances[row] = d;
                best_ids[row] = id;
            }
        }
    }
    Ok((
        PyArray1::from_vec(py, distances),
        PyArray1::from_vec(py, best_ids),
    ))
}

/// Compute one trace's requested elevation change for each affected point.
#[pyfunction]
pub fn relief_contrib_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    distances: numpy::PyReadonlyArray1<'py, f64>,
    mid: numpy::PyReadonlyArray1<'py, f64>,
    dip: numpy::PyReadonlyArray1<'py, f64>,
    params: numpy::PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    if points.shape()[1] != 3
        || points.shape()[0] != distances.len()
        || mid.len() != 3
        || dip.len() != 3
        || params.len() != 8
    {
        return Err(PyValueError::new_err("invalid fault relief arrays"));
    }
    let pts = points.as_slice()?;
    let ds = distances.as_slice()?;
    let mid = mid.as_slice()?;
    let dip = dip.as_slice()?;
    let p = params.as_slice()?;
    let mut result = Vec::with_capacity(ds.len());
    for (point, &d) in pts.chunks_exact(3).zip(ds) {
        let mag = (1.0 - d / p[0]).clamp(0.0, 1.0) * p[1] * p[2] * p[3] * p[4];
        let side = (point[0] - mid[0]) * dip[0]
            + (point[1] - mid[1]) * dip[1]
            + (point[2] - mid[2]) * dip[2]
            > 0.0;
        result.push(if side { p[5] } else { p[6] } * p[7] * mag);
    }
    Ok(PyArray1::from_vec(py, result))
}

/// Balance a single trace's transfer within each crust type, in affected-node order.
#[pyfunction]
pub fn relief_transfer_rs<'py>(
    py: Python<'py>,
    hc: numpy::PyReadonlyArray1<'py, f64>,
    material: numpy::PyReadonlyArray1<'py, f64>,
    areas: numpy::PyReadonlyArray1<'py, f64>,
    kinds: numpy::PyReadonlyArray1<'py, i8>,
    requested: numpy::PyReadonlyArray1<'py, f64>,
    minimum: f64,
    maximum: f64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let n = hc.len();
    if material.len() != n || areas.len() != n || kinds.len() != n || requested.len() != n {
        return Err(PyValueError::new_err(
            "relief transfer arrays must have equal lengths",
        ));
    }
    let mut hc = hc.as_slice()?.to_vec();
    let mut material = material.as_slice()?.to_vec();
    let areas = areas.as_slice()?;
    let kinds = kinds.as_slice()?;
    let requested = requested.as_slice()?;
    let mut unique = kinds.to_vec();
    unique.sort_unstable();
    unique.dedup();
    for kind in unique {
        let mut donors = Vec::new();
        let mut receivers = Vec::new();
        let mut debit_volume = 0.0;
        let mut credit_volume = 0.0;
        for i in 0..n {
            if kinds[i] != kind {
                continue;
            }
            if requested[i] < 0.0 && hc[i] > minimum {
                let debit = (-requested[i]).min(hc[i] - minimum);
                donors.push((i, debit));
                debit_volume += debit * areas[i];
            } else if requested[i] > 0.0 && hc[i] < maximum {
                let credit = requested[i].min(maximum - hc[i]);
                receivers.push((i, credit));
                credit_volume += credit * areas[i];
            }
        }
        let amount = debit_volume.min(credit_volume);
        if amount <= 0.0 {
            continue;
        }
        let mut moved_volume = 0.0;
        for (i, debit) in donors {
            let taken = debit * (amount / debit_volume);
            let moved = taken / hc[i] * material[i];
            moved_volume += moved * areas[i];
            material[i] -= moved;
            hc[i] -= taken;
        }
        for (i, credit) in receivers {
            let given = credit * (amount / credit_volume);
            material[i] += moved_volume * given / amount;
            hc[i] += given;
        }
    }
    Ok((PyArray1::from_vec(py, hc), PyArray1::from_vec(py, material)))
}

/// Consecutive nonzero boundary classifications, split when the neighbor changes.
#[pyfunction]
pub fn boundary_runs_rs<'py>(
    py: Python<'py>,
    regime: numpy::PyReadonlyArray1<'py, i8>,
    owners: numpy::PyReadonlyArray1<'py, i64>,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<i64>>)> {
    let regime = regime.as_slice()?;
    let owners = owners.as_slice()?;
    if regime.len() != owners.len() {
        return Err(PyValueError::new_err("boundary arrays differ in length"));
    }
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    let mut i = 0;
    while i < regime.len() {
        if regime[i] == 0 {
            i += 1;
            continue;
        }
        let mut end = i + 1;
        while end < regime.len() && regime[end] == regime[i] && owners[end] == owners[i] {
            end += 1;
        }
        starts.push(i as i64);
        ends.push(end as i64);
        i = end;
    }
    Ok((PyArray1::from_vec(py, starts), PyArray1::from_vec(py, ends)))
}

/// Resolve an unloading stress change onto fault planes in the caller's arithmetic order.
#[pyfunction]
pub fn unloading_stress_rs<'py>(
    py: Python<'py>,
    loads: numpy::PyReadonlyArray1<'py, f64>,
    sense: f64,
    ratio: f64,
    sine: f64,
    cosine: f64,
    friction: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let out = loads
        .as_slice()?
        .iter()
        .map(|&load| {
            let shear = sense * (1.0 - ratio) * load * sine * cosine;
            let normal = load * (cosine * cosine + ratio * sine * sine);
            shear - friction * normal
        })
        .collect();
    Ok(PyArray1::from_vec(py, out))
}

/// Characteristic magnitude; random jitter is drawn by NumPy in its existing stream.
#[pyfunction]
pub fn quake_magnitude_rs(
    length: f64,
    slip: f64,
    overlap: bool,
    jitter: f64,
    constants: numpy::PyReadonlyArray1<'_, f64>,
) -> PyResult<f64> {
    let c = constants.as_slice()?;
    if c.len() != 6 {
        return Err(PyValueError::new_err("expected six magnitude constants"));
    }
    let value = c[0]
        + c[1] * length.max(1.0).log10()
        + c[2] * slip.max(1.0).log10()
        + if overlap { c[3] } else { 0.0 }
        + jitter;
    Ok(value.clamp(c[4], c[5]))
}

/// Classify near-boundary closing rates using the transform threshold.
#[pyfunction]
pub fn boundary_regime_rs<'py>(
    py: Python<'py>,
    distances: numpy::PyReadonlyArray1<'py, f64>,
    closing: numpy::PyReadonlyArray1<'py, f64>,
    reach: f64,
    threshold: f64,
) -> PyResult<Bound<'py, PyArray1<i8>>> {
    if distances.len() != closing.len() {
        return Err(PyValueError::new_err("boundary arrays differ in length"));
    }
    let result = distances
        .as_slice()?
        .iter()
        .zip(closing.as_slice()?)
        .map(|(&d, &v)| {
            if d > reach || d.is_nan() {
                0
            } else if v > threshold {
                1
            } else if v < -threshold {
                2
            } else if v.abs() <= threshold {
                3
            } else {
                0
            }
        })
        .collect();
    Ok(PyArray1::from_vec(py, result))
}

/// Moment magnitude with bounded seismogenic width.
#[pyfunction]
pub fn moment_magnitude_rs(
    length: f64,
    slip: f64,
    width: f64,
    shear: f64,
    minimum: f64,
    maximum: f64,
) -> f64 {
    let moment = (shear * length * width.min(length) * slip).max(1.0);
    ((2.0 / 3.0) * (moment.log10() - 9.1)).clamp(minimum, maximum)
}

/// Interpolate and normalize a boundary polyline at caller-supplied fractional indices.
#[pyfunction]
pub fn resample_boundary_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    positions: numpy::PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
    use numpy::IntoPyArray;
    let n = points.shape()[0];
    if points.shape()[1] != 3 || n == 0 {
        return Err(PyValueError::new_err("invalid boundary polyline"));
    }
    let pts = points.as_slice()?;
    let positions = positions.as_slice()?;
    let mut output = Vec::with_capacity(3 * positions.len());
    for &position in positions {
        if !position.is_finite() || position < 0.0 || position > (n - 1) as f64 {
            return Err(PyValueError::new_err("invalid interpolation position"));
        }
        let lo = position.floor() as usize;
        let hi = (lo + 1).min(n - 1);
        let fraction = position - lo as f64;
        let mut v = [0.0; 3];
        for k in 0..3 {
            v[k] = pts[3 * lo + k] * (1.0 - fraction) + pts[3 * hi + k] * fraction;
        }
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let norm = if norm < 1e-12 { 1.0 } else { norm };
        output.extend(v.iter().map(|x| x / norm));
    }
    Ok(
        numpy::ndarray::Array2::from_shape_vec((positions.len(), 3), output)
            .unwrap()
            .into_pyarray(py),
    )
}

/// Project inward directions and generate the offset strands of a boundary segment.
#[pyfunction]
pub fn boundary_strands_rs<'py>(
    py: Python<'py>,
    master: PyReadonlyArray2<'py, f64>,
    centroid: numpy::PyReadonlyArray1<'py, f64>,
    count: usize,
    spacing: f64,
) -> PyResult<(
    Bound<'py, numpy::PyArray2<f64>>,
    Bound<'py, numpy::PyArray3<f64>>,
)> {
    use numpy::IntoPyArray;
    if master.shape()[1] != 3 || centroid.len() != 3 {
        return Err(PyValueError::new_err("invalid boundary geometry"));
    }
    let n = master.shape()[0];
    let master = master.as_slice()?;
    let centroid = centroid.as_slice()?;
    let mut inward = Vec::with_capacity(n * 3);
    for p in master.chunks_exact(3) {
        let t = [centroid[0] - p[0], centroid[1] - p[1], centroid[2] - p[2]];
        let dot = t[0] * p[0] + t[1] * p[1] + t[2] * p[2];
        let v = [t[0] - dot * p[0], t[1] - dot * p[1], t[2] - dot * p[2]];
        let norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let norm = if norm < 1e-12 { 1.0 } else { norm };
        inward.extend(v.iter().map(|x| x / norm));
    }
    let mut strands = Vec::with_capacity(count * n * 3);
    for k in 0..count {
        let offset = (k + 1) as f64 * spacing;
        for i in 0..n {
            let p = &master[3 * i..][..3];
            let v = &inward[3 * i..][..3];
            let q = [
                p[0] + offset * v[0],
                p[1] + offset * v[1],
                p[2] + offset * v[2],
            ];
            let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
            let norm = if norm < 1e-12 { 1.0 } else { norm };
            strands.extend(q.iter().map(|x| x / norm));
        }
    }
    Ok((
        numpy::ndarray::Array2::from_shape_vec((n, 3), inward)
            .unwrap()
            .into_pyarray(py),
        numpy::ndarray::Array3::from_shape_vec((count, n, 3), strands)
            .unwrap()
            .into_pyarray(py),
    ))
}
