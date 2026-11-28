//! Numerical kernels for tectonic transport (#1042).

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::spatial::kdtree::{Fused, KdTree, Metric, Plain};

fn weighted_pairs<M: Metric>(
    origins: &[f64],
    points: &[f64],
    hc: &[f64],
    radius: f64,
    planet_km: f64,
    range_km: f64,
    reference_hc: f64,
) -> (Vec<i64>, Vec<i64>, Vec<f64>) {
    let cloud: Vec<[f64; 3]> = points.chunks_exact(3).map(|p| [p[0], p[1], p[2]]).collect();
    let tree = KdTree::build(&cloud);
    let (mut parcels, mut destinations, mut weights) = (Vec::new(), Vec::new(), Vec::new());
    let mut nearby = Vec::new();
    for (p, xyz) in origins.chunks_exact(3).enumerate() {
        nearby.clear();
        tree.ball::<M>([xyz[0], xyz[1], xyz[2]], radius * radius, &mut nearby);
        for &j in &nearby {
            let j = j as usize;
            let dx = xyz[0] - cloud[j][0];
            let dy = xyz[1] - cloud[j][1];
            let dz = xyz[2] - cloud[j][2];
            // NumPy's norm uses the unfused sum, independently of the tree metric.
            let distance = ((dx * dx + dy * dy) + dz * dz).sqrt() * planet_km;
            let thinness = (reference_hc - hc[j]).max(0.0);
            let weight = thinness * (1.0 - distance / range_km).clamp(0.0, 1.0);
            if weight > 0.0 {
                parcels.push(p as i64);
                destinations.push(j as i64);
                weights.push(weight);
            }
        }
    }
    (parcels, destinations, weights)
}

#[pyfunction]
pub fn magma_destination_pairs_rs<'py>(
    py: Python<'py>,
    origins: PyReadonlyArray2<'py, f64>,
    points: PyReadonlyArray2<'py, f64>,
    hc: PyReadonlyArray1<'py, f64>,
    radius: f64,
    planet_km: f64,
    range_km: f64,
    reference_hc: f64,
    fused: bool,
) -> PyResult<(
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
)> {
    if origins.shape()[1] != 3
        || points.shape()[1] != 3
        || points.shape()[0] != hc.len()
        || points.shape()[0] >= u32::MAX as usize
        || !radius.is_finite()
        || radius < 0.0
        || !range_km.is_finite()
        || range_km <= 0.0
    {
        return Err(PyValueError::new_err("invalid magma destination inputs"));
    }
    let (origins, points, hc) = (origins.as_slice()?, points.as_slice()?, hc.as_slice()?);
    let (p, d, w) = py.detach(|| {
        if fused {
            weighted_pairs::<Fused>(
                origins,
                points,
                hc,
                radius,
                planet_km,
                range_km,
                reference_hc,
            )
        } else {
            weighted_pairs::<Plain>(
                origins,
                points,
                hc,
                radius,
                planet_km,
                range_km,
                reference_hc,
            )
        }
    });
    Ok((
        PyArray1::from_vec(py, p),
        PyArray1::from_vec(py, d),
        PyArray1::from_vec(py, w),
    ))
}

/// Pick the first neighbour with the minimum distance, as NumPy's strict `<` update does.
#[pyfunction]
pub fn boundary_nearest_rs<'py>(
    py: Python<'py>,
    distances: Vec<PyReadonlyArray1<'py, f64>>,
    indices: Vec<PyReadonlyArray1<'py, i64>>,
) -> PyResult<(
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
)> {
    let n = distances.first().map_or(0, |a| a.len());
    if distances.len() != indices.len()
        || distances.iter().any(|a| a.len() != n)
        || indices.iter().any(|a| a.len() != n)
    {
        return Err(PyValueError::new_err("distance and index shapes differ"));
    }
    let ds = distances
        .iter()
        .map(|a| a.as_slice())
        .collect::<Result<Vec<_>, _>>()?;
    let ids = indices
        .iter()
        .map(|a| a.as_slice())
        .collect::<Result<Vec<_>, _>>()?;
    let (best, owner, index) = py.detach(|| {
        let (mut best, mut owner, mut index) = (vec![f64::INFINITY; n], vec![-1; n], vec![-1; n]);
        for neighbour in 0..ds.len() {
            for i in 0..n {
                if ds[neighbour][i] < best[i] {
                    best[i] = ds[neighbour][i];
                    owner[i] = neighbour as i64;
                    index[i] = ids[neighbour][i];
                }
            }
        }
        (best, owner, index)
    });
    Ok((
        PyArray1::from_vec(py, best),
        PyArray1::from_vec(py, owner),
        PyArray1::from_vec(py, index),
    ))
}

fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

fn front_labels<M: Metric>(points: &[f64], radius: f64) -> Vec<i64> {
    let cloud: Vec<[f64; 3]> = points.chunks_exact(3).map(|p| [p[0], p[1], p[2]]).collect();
    let tree = KdTree::build(&cloud);
    let mut parent: Vec<usize> = (0..cloud.len()).collect();
    let mut nearby = Vec::new();
    for (i, &point) in cloud.iter().enumerate() {
        nearby.clear();
        tree.ball::<M>(point, radius * radius, &mut nearby);
        for &j in &nearby {
            if (j as usize) > i {
                let a = root(&mut parent, i);
                let b = root(&mut parent, j as usize);
                if a != b {
                    parent[b] = a;
                }
            }
        }
    }
    let mut root_label = vec![-1; cloud.len()];
    let mut labels = vec![-1; cloud.len()];
    let mut next = 0;
    for (i, label) in labels.iter_mut().enumerate() {
        let r = root(&mut parent, i);
        if root_label[r] < 0 {
            root_label[r] = next;
            next += 1;
        }
        *label = root_label[r];
    }
    labels
}

#[pyfunction]
pub fn collision_front_labels_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    radius: f64,
    fused: bool,
) -> PyResult<Bound<'py, PyArray1<i64>>> {
    if points.shape()[1] != 3
        || points.shape()[0] >= u32::MAX as usize
        || !radius.is_finite()
        || radius < 0.0
    {
        return Err(PyValueError::new_err("invalid collision front inputs"));
    }
    let points = points.as_slice()?;
    let labels = py.detach(|| {
        if fused {
            front_labels::<Fused>(points, radius)
        } else {
            front_labels::<Plain>(points, radius)
        }
    });
    Ok(PyArray1::from_vec(py, labels))
}

/// Test adaptive cubed-sphere cells; boundary-adjacent rows are returned for Python to resolve.
#[pyfunction]
pub fn gap_contains_rs<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    keys: PyReadonlyArray1<'py, i64>,
    cells_per_edge: usize,
    max_level: usize,
) -> PyResult<(Bound<'py, PyArray1<bool>>, Bound<'py, PyArray1<bool>>)> {
    if points.shape()[1] != 3
        || cells_per_edge == 0
        || max_level >= 24
        || cells_per_edge
            .checked_shl(max_level as u32)
            .is_none_or(|n| n >= (1 << 24))
    {
        return Err(PyValueError::new_err("invalid gap cell inputs"));
    }
    let (points, keys) = (points.as_slice()?, keys.as_slice()?);
    let (covered, uncertain) = py.detach(|| {
        let mut covered = vec![false; points.len() / 3];
        let mut uncertain = vec![false; points.len() / 3];
        for (row, p) in points.chunks_exact(3).enumerate() {
            if !p.iter().all(|v| v.is_finite()) {
                uncertain[row] = true;
                continue;
            }
            let axis = if p[0].abs() >= p[1].abs() && p[0].abs() >= p[2].abs() {
                0
            } else if p[1].abs() >= p[2].abs() {
                1
            } else {
                2
            };
            let face = match (axis, p[axis] > 0.0) {
                (0, true) => 0,
                (1, true) => 1,
                (0, false) => 2,
                (1, false) => 3,
                (2, true) => 4,
                _ => 5,
            };
            let (u, v) = match face {
                0 => (p[1], p[2]),
                1 => (-p[0], p[2]),
                2 => (-p[1], p[2]),
                3 => (p[0], p[2]),
                4 => (p[1], -p[0]),
                _ => (p[1], p[0]),
            };
            let depth = p[axis].abs();
            if depth == 0.0 || !depth.is_finite() {
                uncertain[row] = true;
                continue;
            }
            let (alpha, beta) = ((u / depth).atan(), (v / depth).atan());
            for level in 0..=max_level {
                let n = cells_per_edge << level;
                let step = std::f64::consts::FRAC_PI_2 / n as f64;
                let a = (alpha + std::f64::consts::FRAC_PI_4) / step;
                let b = (beta + std::f64::consts::FRAC_PI_4) / step;
                // At fine levels, scaling the angular calculation by n magnifies
                // atan and division rounding; defer every near-edge point to NumPy (#1042).
                let tolerance = 1e-10_f64.max(32.0 * n as f64 * f64::EPSILON);
                if (a - a.round()).abs() < tolerance || (b - b.round()).abs() < tolerance {
                    uncertain[row] = true;
                }
                let i = (a.floor() as i64).clamp(0, n as i64 - 1);
                let j = (b.floor() as i64).clamp(0, n as i64 - 1);
                let key = ((face as i64) << 56) | ((level as i64) << 48) | (j << 24) | i;
                if keys.binary_search(&key).is_ok() {
                    covered[row] = true;
                }
            }
        }
        (covered, uncertain)
    });
    Ok((
        PyArray1::from_vec(py, covered),
        PyArray1::from_vec(py, uncertain),
    ))
}

/// Allocate parcel volumes with one shared cap per destination, preserving pair order.
#[pyfunction]
pub fn magma_allocate_rs<'py>(
    py: Python<'py>,
    parcel_idx: PyReadonlyArray1<'py, i64>,
    dest_idx: PyReadonlyArray1<'py, i64>,
    weights: PyReadonlyArray1<'py, f64>,
    volumes: PyReadonlyArray1<'py, f64>,
    areas: PyReadonlyArray1<'py, f64>,
    hc: PyReadonlyArray1<'py, f64>,
    cap: f64,
    max_hc: f64,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<f64>>)> {
    let (p, d, w, v, area, hc) = (
        parcel_idx.as_slice()?,
        dest_idx.as_slice()?,
        weights.as_slice()?,
        volumes.as_slice()?,
        areas.as_slice()?,
        hc.as_slice()?,
    );
    if p.len() != d.len()
        || p.len() != w.len()
        || area.len() != hc.len()
        || p.iter().any(|&i| i < 0 || i as usize >= v.len())
        || d.iter().any(|&i| i < 0 || i as usize >= area.len())
    {
        return Err(PyValueError::new_err("invalid magma allocation inputs"));
    }
    let (placed, realized) = py.detach(|| {
        let mut sums = vec![0.0; v.len()];
        for i in 0..p.len() {
            sums[p[i] as usize] += w[i];
        }
        let mut requested = vec![0.0; p.len()];
        let mut totals = vec![0.0; area.len()];
        for i in 0..p.len() {
            let (pi, di) = (p[i] as usize, d[i] as usize);
            requested[i] = v[pi] * (w[i] / sums[pi]) / area[di];
            totals[di] += requested[i];
        }
        let mut realized = vec![0.0; area.len()];
        let mut ratios = vec![0.0; area.len()];
        for i in 0..area.len() {
            realized[i] = totals[i].min(cap).min((max_hc - hc[i]).max(0.0));
            if totals[i] > 0.0 {
                ratios[i] = realized[i] / totals[i];
            }
        }
        let mut placed = vec![0.0; v.len()];
        for i in 0..p.len() {
            let di = d[i] as usize;
            placed[p[i] as usize] += (requested[i] * ratios[di]) * area[di];
        }
        (placed, realized)
    });
    Ok((
        PyArray1::from_vec(py, placed),
        PyArray1::from_vec(py, realized),
    ))
}
