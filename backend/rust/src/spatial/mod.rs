//! Shared 3D point index for simulation and rendering (#1042).
//!
//! `KdIndex` wraps the static k-d tree in `kdtree.rs` as a Python class. Its batch queries
//! return flat arrays (a CSR `(indptr, indices)` pair for radius queries) instead of
//! `cKDTree.query_ball_point`'s Python lists, and run across the rayon pool without the GIL once
//! a batch is large enough to pay for it. `backend/app/spatial_index.py` is its Python face.

pub mod kdtree;

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

use kdtree::{Fused, KdTree, Metric, Plain};

/// Batches smaller than this run on the calling thread with the GIL held: below it, handing
/// work to the rayon pool and releasing the GIL cost more than they save. Most step-time calls
/// are either a few points along a fault trace or tens of thousands of nodes.
const PARALLEL_MIN_QUERIES: usize = 1024;
/// Query points per rayon task.
const CHUNK: usize = 256;

#[inline(always)]
fn at(xs: &[f64], i: usize) -> [f64; 3] {
    [xs[3 * i], xs[3 * i + 1], xs[3 * i + 2]]
}

fn rows<'a>(arr: &'a PyReadonlyArray2<'_, f64>, what: &str) -> PyResult<&'a [f64]> {
    if arr.shape()[1] != 3 {
        return Err(PyValueError::new_err(format!("{what} must have shape (n, 3)")));
    }
    Ok(arr.as_slice()?)
}

/// Run `f` over `[0, m)` in `CHUNK`-sized pieces, in parallel when `m` is large, and return the
/// per-chunk results in order.
fn chunked<T: Send>(m: usize, f: impl Fn(usize, usize) -> T + Sync) -> Vec<T> {
    let n_chunks = m.div_ceil(CHUNK);
    let piece = |c: usize| f(c * CHUNK, ((c + 1) * CHUNK).min(m));
    if m >= PARALLEL_MIN_QUERIES {
        (0..n_chunks).into_par_iter().map(piece).collect()
    } else {
        (0..n_chunks).map(piece).collect()
    }
}

fn query_impl<M: Metric>(tree: &KdTree, xs: &[f64], k: usize, ub2: f64) -> (Vec<f64>, Vec<i64>) {
    let m = xs.len() / 3;
    let missing = tree.len() as i64;
    let mut dist = vec![f64::INFINITY; m * k];
    let mut idx = vec![missing; m * k];
    if k == 0 || m == 0 {
        return (dist, idx);
    }
    let fill = |first: usize, d: &mut [f64], i: &mut [i64]| {
        if k == 1 {
            for (row, (dr, ir)) in d.iter_mut().zip(i.iter_mut()).enumerate() {
                if let Some((d2, id)) = tree.nearest::<M>(at(xs, first + row), ub2) {
                    *dr = d2.sqrt();
                    *ir = id as i64;
                }
            }
            return;
        }
        let mut best = Vec::with_capacity(k + 1);
        for (row, (dr, ir)) in d.chunks_mut(k).zip(i.chunks_mut(k)).enumerate() {
            tree.knn::<M>(at(xs, first + row), k, ub2, &mut best);
            for (slot, &(d2, id)) in best.iter().enumerate() {
                dr[slot] = d2.sqrt();
                ir[slot] = id as i64;
            }
        }
    };
    if m >= PARALLEL_MIN_QUERIES {
        dist.par_chunks_mut(k * CHUNK)
            .zip(idx.par_chunks_mut(k * CHUNK))
            .enumerate()
            .for_each(|(c, (d, i))| fill(c * CHUNK, d, i));
    } else {
        fill(0, &mut dist, &mut idx);
    }
    (dist, idx)
}

#[inline(always)]
fn radius_sq(radii: &[f64], row: usize) -> f64 {
    let r = if radii.len() == 1 { radii[0] } else { radii[row] };
    r * r
}

fn ball_impl<M: Metric>(tree: &KdTree, xs: &[f64], radii: &[f64]) -> (Vec<i64>, Vec<i64>) {
    let m = xs.len() / 3;
    let parts = chunked(m, |lo, hi| {
        let mut counts = Vec::with_capacity(hi - lo);
        let mut ids = Vec::new();
        for row in lo..hi {
            let before = ids.len();
            tree.ball::<M>(at(xs, row), radius_sq(radii, row), &mut ids);
            counts.push(ids.len() - before);
        }
        (counts, ids)
    });
    let total: usize = parts.iter().map(|p| p.1.len()).sum();
    let mut indptr = Vec::with_capacity(m + 1);
    let mut indices = Vec::with_capacity(total);
    indptr.push(0_i64);
    for (counts, ids) in parts {
        for c in counts {
            let last = *indptr.last().unwrap();
            indptr.push(last + c as i64);
        }
        indices.extend(ids.into_iter().map(i64::from));
    }
    (indptr, indices)
}

fn count_impl<M: Metric>(tree: &KdTree, xs: &[f64], radii: &[f64]) -> Vec<i64> {
    let m = xs.len() / 3;
    chunked(m, |lo, hi| {
        (lo..hi)
            .map(|row| tree.ball_count::<M>(at(xs, row), radius_sq(radii, row)) as i64)
            .collect::<Vec<_>>()
    })
    .concat()
}

fn pairs_impl<M: Metric>(tree: &KdTree, r2: f64) -> Vec<i64> {
    let pts = tree.points_by_id();
    let flat: Vec<f64> = pts.iter().flatten().copied().collect();
    chunked(pts.len(), |lo, hi| {
        let mut out = Vec::new();
        let mut near = Vec::new();
        for i in lo..hi {
            near.clear();
            tree.ball::<M>(at(&flat, i), r2, &mut near);
            for &j in near.iter().filter(|&&j| j as usize > i) {
                out.push(i as i64);
                out.push(j as i64);
            }
        }
        out
    })
    .concat()
}

/// Immutable 3D point index with `scipy.spatial.cKDTree`'s query semantics (#1042).
///
/// `fused` picks the squared-distance rounding (see `kdtree::Metric`); the Python wrapper
/// probes the installed scipy for it. Results are flat: `query` returns `(m * k)` distances and
/// indices (misses are `inf` and `n`), `query_ball_point` returns CSR `(indptr, indices)` with
/// each row ascending, and `query_pairs` returns `i < j` pairs in lexicographic order.
#[pyclass(frozen, module = "mantle_bloom_rust")]
pub struct KdIndex {
    tree: KdTree,
    fused: bool,
}

impl KdIndex {
    /// Run `f` without the GIL when the batch is big enough to share across threads.
    fn run<T: Send>(py: Python<'_>, m: usize, f: impl FnOnce() -> T + Send) -> T {
        if m >= PARALLEL_MIN_QUERIES {
            py.detach(f)
        } else {
            f()
        }
    }
}

#[pymethods]
impl KdIndex {
    #[new]
    fn new(py: Python<'_>, points: PyReadonlyArray2<'_, f64>, fused: bool) -> PyResult<Self> {
        let xs = rows(&points, "points")?;
        let n = xs.len() / 3;
        if n >= u32::MAX as usize {
            return Err(PyValueError::new_err("point cloud exceeds u32 index range"));
        }
        let tree = Self::run(py, n, || {
            let pts: Vec<[f64; 3]> = (0..n).map(|i| at(xs, i)).collect();
            KdTree::build(&pts)
        });
        Ok(Self { tree, fused })
    }

    #[getter]
    fn n(&self) -> usize {
        self.tree.len()
    }

    #[getter]
    fn fused(&self) -> bool {
        self.fused
    }

    /// The `k` nearest points to each row of `x` closer than `distance_upper_bound`.
    fn query<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray2<'py, f64>,
        k: usize,
        distance_upper_bound: f64,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i64>>)> {
        let xs = rows(&x, "x")?;
        let ub2 = distance_upper_bound * distance_upper_bound;
        let (d, i) = Self::run(py, xs.len() / 3, || {
            if self.fused {
                query_impl::<Fused>(&self.tree, xs, k, ub2)
            } else {
                query_impl::<Plain>(&self.tree, xs, k, ub2)
            }
        });
        Ok((PyArray1::from_vec(py, d), PyArray1::from_vec(py, i)))
    }

    /// Every point within `r` (one radius, or one per row) of each row of `x`, as CSR.
    fn query_ball_point<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray2<'py, f64>,
        r: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<i64>>)> {
        let xs = rows(&x, "x")?;
        let radii = r.as_slice()?;
        let m = xs.len() / 3;
        if radii.len() != 1 && radii.len() != m {
            return Err(PyValueError::new_err("r must hold one radius or one per query point"));
        }
        let (indptr, indices) = Self::run(py, m, || {
            if self.fused {
                ball_impl::<Fused>(&self.tree, xs, radii)
            } else {
                ball_impl::<Plain>(&self.tree, xs, radii)
            }
        });
        Ok((PyArray1::from_vec(py, indptr), PyArray1::from_vec(py, indices)))
    }

    /// How many points lie within `r` of each row of `x`.
    fn query_ball_point_count<'py>(
        &self,
        py: Python<'py>,
        x: PyReadonlyArray2<'py, f64>,
        r: PyReadonlyArray1<'py, f64>,
    ) -> PyResult<Bound<'py, PyArray1<i64>>> {
        let xs = rows(&x, "x")?;
        let radii = r.as_slice()?;
        let m = xs.len() / 3;
        if radii.len() != 1 && radii.len() != m {
            return Err(PyValueError::new_err("r must hold one radius or one per query point"));
        }
        let counts = Self::run(py, m, || {
            if self.fused {
                count_impl::<Fused>(&self.tree, xs, radii)
            } else {
                count_impl::<Plain>(&self.tree, xs, radii)
            }
        });
        Ok(PyArray1::from_vec(py, counts))
    }

    /// Every `i < j` pair within `r` of each other, flattened as `[i0, j0, i1, j1, ...]`.
    fn query_pairs<'py>(&self, py: Python<'py>, r: f64) -> Bound<'py, PyArray1<i64>> {
        let r2 = r * r;
        let flat = Self::run(py, self.tree.len(), || {
            if self.fused {
                pairs_impl::<Fused>(&self.tree, r2)
            } else {
                pairs_impl::<Plain>(&self.tree, r2)
            }
        });
        PyArray1::from_vec(py, flat)
    }
}
