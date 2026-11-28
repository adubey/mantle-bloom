//! World-generation kernels for `continent_layout`, `plate_partition` and
//! `synthetic_relief` (#1042). Every kernel is deterministic and consumes no
//! randomness: all `rng` draws stay in Python, so a seeded world is unchanged.

use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

use crate::spatial::kdtree::KdTree;

pub const PLANET_RADIUS_KM: f64 = 6371.0;

/// Total-order wrapper for heap entries. Pop order need not match Python's
/// `(level, node)` tuples (and `total_cmp` separates `-0.0` from `0.0`, where
/// Python falls through to the node index): with non-negative rises the flood
/// is a monotone Dijkstra whose fixpoint is independent of pop order.
#[derive(Clone, Copy, PartialEq)]
struct OrdF64(f64);

impl Eq for OrdF64 {}
impl PartialOrd for OrdF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for OrdF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn chord_to_arc_km(chord: f64) -> f64 {
    2.0 * (chord / 2.0).clamp(0.0, 1.0).asin() * PLANET_RADIUS_KM
}

/// Priority flood: the lowest level each node can drain at, rising by `rise`
/// per edge (`synthetic_relief._flood`). `indptr`/`indices` are CSR reverse
/// adjacency; a `None` rise means every edge rises by zero.
#[pyfunction]
pub fn priority_flood_rs<'py>(
    py: Python<'py>,
    indptr: PyReadonlyArray1<'py, i64>,
    indices: PyReadonlyArray1<'py, i64>,
    rise: Option<PyReadonlyArray1<'py, f64>>,
    seeds_level: PyReadonlyArray1<'py, f64>,
    elevation: PyReadonlyArray1<'py, f64>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let n = seeds_level.len();
    if elevation.len() != n || indptr.len() != n + 1 {
        return Err(PyValueError::new_err("incompatible flood shapes"));
    }
    let indptr = indptr.as_slice()?.to_vec();
    if indptr[0] != 0
        || indptr.windows(2).any(|w| w[1] < w[0])
        || indptr[n] as usize != indices.len()
    {
        return Err(PyValueError::new_err("indptr is not a valid CSR pointer array"));
    }
    let indices = indices.as_slice()?.to_vec();
    let rise: Vec<f64> = match &rise {
        Some(r) => {
            if r.len() != indices.len() {
                return Err(PyValueError::new_err("rise must match indices"));
            }
            r.as_slice()?.to_vec()
        }
        None => vec![0.0; indices.len()],
    };
    let seeds = seeds_level.as_slice()?.to_vec();
    let z = elevation.as_slice()?.to_vec();
    // Python's `max(z, reach)` keeps a NaN elevation; Rust's `f64::max` drops
    // it. Reject NaN instead of diverging (cf. `route_downstream_rs`).
    if seeds.iter().chain(z.iter()).any(|v| v.is_nan()) {
        return Err(PyValueError::new_err("flood levels must not be NaN"));
    }
    if indices.iter().any(|&j| j < 0 || j as usize >= n) {
        return Err(PyValueError::new_err("flood index out of range"));
    }
    let level = py.detach(move || {
        let mut level = seeds;
        let mut heap: BinaryHeap<Reverse<(OrdF64, usize)>> = BinaryHeap::new();
        for (i, &lv) in level.iter().enumerate() {
            if lv.is_finite() {
                heap.push(Reverse((OrdF64(lv), i)));
            }
        }
        while let Some(Reverse((OrdF64(here), i))) = heap.pop() {
            if here > level[i] {
                continue;
            }
            for e in indptr[i] as usize..indptr[i + 1] as usize {
                let j = indices[e] as usize;
                let reach = z[j].max(here + rise[e]);
                if reach < level[j] {
                    level[j] = reach;
                    heap.push(Reverse((OrdF64(reach), j)));
                }
            }
        }
        level
    });
    Ok(PyArray1::from_vec(py, level))
}

/// Great-circle distance (km) and index of the nearest of `sources` for every
/// point (`synthetic_relief._distance_km`, `continent_layout` queries).
/// Sources and points are unit vectors; ties keep the lowest source index.
#[pyfunction]
pub fn nearest_distance_rs<'py>(
    py: Python<'py>,
    sources: PyReadonlyArray2<'py, f64>,
    points: PyReadonlyArray2<'py, f64>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i64>>)> {
    let s_shape = sources.shape();
    let p_shape = points.shape();
    if s_shape[1] != 3 || p_shape[1] != 3 {
        return Err(PyValueError::new_err("sources and points must be (m, 3) and (n, 3)"));
    }
    let (m, n) = (s_shape[0], p_shape[0]);
    let sources = sources.as_slice()?.to_vec();
    let points = points.as_slice()?.to_vec();
    if sources.iter().chain(points.iter()).any(|v| v.is_nan()) {
        return Err(PyValueError::new_err("nearest inputs must not be NaN"));
    }
    let (dist, index) = py.detach(move || {
        if m == 0 {
            return (vec![f64::INFINITY; n], vec![0_i64; n]);
        }
        let (dist, index): (Vec<f64>, Vec<i64>) = (0..n)
            .into_par_iter()
            .map(|i| {
                let (px, pz) = (points[i * 3], points[i * 3 + 2]);
                let py_ = points[i * 3 + 1];
                let mut best_chord2 = f64::INFINITY;
                let mut best = 0_usize;
                for s in 0..m {
                    let dx = px - sources[s * 3];
                    let dy = py_ - sources[s * 3 + 1];
                    let dz = pz - sources[s * 3 + 2];
                    let chord2 = dx * dx + dy * dy + dz * dz;
                    if chord2 < best_chord2 {
                        best_chord2 = chord2;
                        best = s;
                    }
                }
                (chord_to_arc_km(best_chord2.sqrt()), best as i64)
            })
            .unzip();
        (dist, index)
    });
    Ok((PyArray1::from_vec(py, dist), PyArray1::from_vec(py, index)))
}

/// Seed assignment (`plate_partition`'s Voronoi step): the owner label of the
/// nearest seed for every lattice cell, from a balanced k-d tree over the
/// seeds (`spatial::kdtree::KdTree`, the shared tree). Exact ties resolve
/// to the lowest seed index; `cKDTree` promises no tie order.
#[pyfunction]
pub fn nearest_seed_labels_rs<'py>(
    py: Python<'py>,
    xyz: PyReadonlyArray2<'py, f64>,
    seeds: PyReadonlyArray2<'py, f64>,
    seed_owner: PyReadonlyArray1<'py, i64>,
) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let x_shape = xyz.shape();
    let s_shape = seeds.shape();
    if x_shape[1] != 3 || s_shape[1] != 3 || seed_owner.len() != s_shape[0] {
        return Err(PyValueError::new_err("incompatible seed-label shapes"));
    }
    let (n, m) = (x_shape[0], s_shape[0]);
    if m == 0 {
        return Err(PyValueError::new_err("at least one seed is required"));
    }
    let xyz = xyz.as_slice()?.to_vec();
    let seeds = seeds.as_slice()?.to_vec();
    if xyz.iter().chain(seeds.iter()).any(|v| v.is_nan()) {
        return Err(PyValueError::new_err("seed-label inputs must not be NaN"));
    }
    let owners = seed_owner.as_slice()?.to_vec();
    let labels = py.detach(move || {
        let seed_pts: Vec<[f64; 3]> = (0..m).map(|s| [seeds[s * 3], seeds[s * 3 + 1], seeds[s * 3 + 2]]).collect();
        let tree = KdTree::build(&seed_pts);
        (0..n)
            .into_par_iter()
            .map(|i| owners[tree.query_nearest([xyz[i * 3], xyz[i * 3 + 1], xyz[i * 3 + 2]])])
            .collect::<Vec<i64>>()
    });
    Ok(PyArray1::from_vec(py, labels))
}

/// Connected components of the lattice graph with links cut between labels
/// (`plate_partition.label_components`): the component count only.
#[pyfunction]
pub fn label_component_count_rs(
    offsets: PyReadonlyArray1<'_, i64>,
    neighbours: PyReadonlyArray1<'_, i64>,
    labels: PyReadonlyArray1<'_, i64>,
) -> PyResult<i64> {
    let n = labels.len();
    if offsets.len() != n + 1 {
        return Err(PyValueError::new_err("incompatible component shapes"));
    }
    let offsets = offsets.as_slice()?;
    if offsets[0] != 0
        || offsets.windows(2).any(|w| w[1] < w[0])
        || offsets[n] as usize != neighbours.len()
    {
        return Err(PyValueError::new_err("offsets is not a valid CSR pointer array"));
    }
    let neighbours = neighbours.as_slice()?;
    let labels = labels.as_slice()?;
    if neighbours.iter().any(|&j| j < 0 || j as usize >= n) {
        return Err(PyValueError::new_err("neighbour index out of range"));
    }
    let mut seen = vec![false; n];
    let mut stack: Vec<usize> = Vec::new();
    let mut count = 0_i64;
    for start in 0..n {
        if seen[start] {
            continue;
        }
        count += 1;
        seen[start] = true;
        stack.push(start);
        while let Some(u) = stack.pop() {
            for e in offsets[u] as usize..offsets[u + 1] as usize {
                let v = neighbours[e] as usize;
                if !seen[v] && labels[v] == labels[u] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_matches_chord_conversion() {
        // Antipodal unit vectors: chord 2 gives half the circumference.
        assert!((chord_to_arc_km(2.0) - std::f64::consts::PI * PLANET_RADIUS_KM).abs() < 1e-9);
        assert_eq!(chord_to_arc_km(0.0), 0.0);
    }

    #[test]
    fn ord_f64_orders_like_python_floats() {
        assert!(OrdF64(1.0) < OrdF64(2.0));
        assert!(OrdF64(f64::INFINITY) > OrdF64(1e308));
        assert_eq!(OrdF64(0.0).cmp(&OrdF64(-0.0)), Ordering::Greater);
    }
}
