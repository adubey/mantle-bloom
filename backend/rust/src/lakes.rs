use std::cmp::Ordering;

use numpy::{Element, PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[inline]
fn find_halving(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

#[inline]
fn find_compress(parent: &mut [usize], mut x: usize) -> usize {
    let mut root = x;
    while parent[root] != root {
        root = parent[root];
    }
    while parent[x] != root {
        let next = parent[x];
        parent[x] = root;
        x = next;
    }
    root
}

#[inline]
fn copy_1d<T: Copy + Element>(arr: &PyReadonlyArray1<'_, T>) -> PyResult<Vec<T>> {
    if arr.is_c_contiguous() {
        Ok(arr.as_slice()?.to_vec())
    } else {
        Ok(arr.as_array().iter().copied().collect())
    }
}

#[inline]
fn copy_2d_row_major<T: Copy + Element>(arr: &PyReadonlyArray2<'_, T>) -> PyResult<Vec<T>> {
    if arr.is_c_contiguous() {
        Ok(arr.as_slice()?.to_vec())
    } else {
        Ok(arr.as_array().iter().copied().collect())
    }
}

#[inline]
fn max_endpoint_elev(ea: f64, eb: f64) -> f64 {
    if ea.is_nan() || eb.is_nan() {
        f64::NAN
    } else if ea > eb {
        ea
    } else if eb > ea {
        eb
    } else if ea.is_sign_positive() {
        ea
    } else {
        eb
    }
}

#[inline]
fn cmp_f64_numpy_stable(a: f64, b: f64) -> Ordering {
    match (a.is_nan(), b.is_nan()) {
        (true, true) => Ordering::Equal,
        (true, false) => Ordering::Greater,
        (false, true) => Ordering::Less,
        (false, false) => a.partial_cmp(&b).unwrap_or(Ordering::Equal),
    }
}

pub fn lake_components_core(
    is_lake: &[bool],
    neighbor_idx: &[i64],
    k: usize,
) -> Result<(Vec<i64>, Vec<i64>), String> {
    let n = is_lake.len();
    if neighbor_idx.len() != n * k {
        return Err(format!(
            "neighbor_idx length {} does not match n*k ({n}*{k})",
            neighbor_idx.len()
        ));
    }
    if !is_lake.iter().any(|&b| b) {
        return Ok((vec![0], Vec::new()));
    }

    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        if !is_lake[i] {
            continue;
        }
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb in row {
            if nb < 0 || (nb as usize) >= n {
                return Err(format!("neighbor_idx entry {nb} out of range [0, {n})"));
            }
            let j = nb as usize;
            if is_lake[j] {
                let ri = find_halving(&mut parent, i);
                let rj = find_halving(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }

    let mut group_id = vec![u32::MAX; n];
    let mut group_sizes: Vec<i64> = Vec::new();
    let mut node_groups: Vec<(usize, usize)> = Vec::new();

    for i in 0..n {
        if !is_lake[i] {
            continue;
        }
        let r = find_halving(&mut parent, i);
        let gid = if group_id[r] == u32::MAX {
            let id = group_sizes.len() as u32;
            group_id[r] = id;
            group_sizes.push(0);
            id as usize
        } else {
            group_id[r] as usize
        };
        group_sizes[gid] += 1;
        node_groups.push((i, gid));
    }

    let mut indptr = Vec::with_capacity(group_sizes.len() + 1);
    let mut acc = 0_i64;
    indptr.push(0);
    for &sz in &group_sizes {
        acc += sz;
        indptr.push(acc);
    }

    let mut cursor: Vec<usize> = indptr[..group_sizes.len()]
        .iter()
        .map(|&v| v as usize)
        .collect();
    let mut members = vec![0_i64; acc as usize];
    for (node, gid) in node_groups {
        members[cursor[gid]] = node as i64;
        cursor[gid] += 1;
    }

    Ok((indptr, members))
}

pub fn catchment_roots_core(
    elevation: &[f64],
    is_ocean: &[bool],
    neighbor_idx: &[i64],
    k: usize,
    ocean_catchment: i64,
) -> Result<Vec<i64>, String> {
    let n = elevation.len();
    if is_ocean.len() != n {
        return Err(format!(
            "is_ocean length {} does not match elevation length {n}",
            is_ocean.len()
        ));
    }
    if neighbor_idx.len() != n * k {
        return Err(format!(
            "neighbor_idx length {} does not match n*k ({n}*{k})",
            neighbor_idx.len()
        ));
    }
    if n == 0 {
        return Ok(Vec::new());
    }
    if k == 0 {
        return Err("neighbor_idx must have at least one column when n > 0".to_string());
    }

    let mut flow_target = vec![-1_i64; n];
    for i in 0..n {
        let own_z = elevation[i];
        let mut best_val = f64::INFINITY;
        let mut best_nb: i64 = -1;
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb_raw in row {
            if nb_raw < 0 || (nb_raw as usize) >= n {
                return Err(format!("neighbor_idx entry {nb_raw} out of range [0, {n})"));
            }
            let nb = nb_raw as usize;
            let nz = elevation[nb];
            if nz < own_z && nz < best_val {
                best_val = nz;
                best_nb = nb_raw;
            }
        }
        if best_val.is_finite() {
            flow_target[i] = best_nb;
        }
    }

    let mut catchment_root = vec![-1_i64; n];
    let mut path: Vec<usize> = Vec::new();
    for start in 0..n {
        if catchment_root[start] != -1 {
            continue;
        }
        path.clear();
        let mut x = start;
        loop {
            if catchment_root[x] != -1 {
                break;
            }
            path.push(x);
            if is_ocean[x] {
                catchment_root[x] = ocean_catchment;
                break;
            }
            let target = flow_target[x];
            if target < 0 {
                catchment_root[x] = x as i64;
                break;
            }
            x = target as usize;
        }
        let result = catchment_root[x];
        for &node in &path {
            catchment_root[node] = result;
        }
    }

    Ok(catchment_root)
}

#[derive(Copy, Clone)]
struct BoundaryEdge {
    a: u32,
    b: u32,
    w: f64,
}

pub struct HierarchyEvents {
    pub catchment_root: Vec<i64>,
    pub ev_kind: Vec<i8>,
    pub ev_a: Vec<i64>,
    pub ev_b: Vec<i64>,
    pub ev_w: Vec<f64>,
    pub ev_ra: Vec<i64>,
    pub ev_rb: Vec<i64>,
}

pub fn lake_hierarchy_events_core(
    elevation: &[f64],
    is_ocean: &[bool],
    neighbor_idx: &[i64],
    k: usize,
    interface_pass: Option<&[f64]>,
    ocean_catchment: i64,
) -> Result<HierarchyEvents, String> {
    let n = elevation.len();
    if n > (u32::MAX as usize) {
        return Err("node count exceeds u32 index range".to_string());
    }
    if let Some(pass) = interface_pass {
        if pass.len() != n * k {
            return Err(format!(
                "interface_pass length {} does not match n*k ({n}*{k})",
                pass.len()
            ));
        }
    }
    if n == 0 || k == 0 {
        return Ok(HierarchyEvents {
            catchment_root: Vec::new(),
            ev_kind: Vec::new(),
            ev_a: Vec::new(),
            ev_b: Vec::new(),
            ev_w: Vec::new(),
            ev_ra: Vec::new(),
            ev_rb: Vec::new(),
        });
    }

    let catchment_root = catchment_roots_core(elevation, is_ocean, neighbor_idx, k, ocean_catchment)?;

    let ocean_idx = n;
    let mut parent: Vec<usize> = (0..=n).collect();
    for node in 0..n {
        let root = catchment_root[node];
        if root != ocean_catchment {
            let r = root as usize;
            if node != r {
                let rn = find_compress(&mut parent, node);
                let rr = find_compress(&mut parent, r);
                parent[rn] = rr;
            }
        }
    }
    for node in 0..n {
        if catchment_root[node] == ocean_catchment {
            let rn = find_compress(&mut parent, node);
            let ro = find_compress(&mut parent, ocean_idx);
            parent[rn] = ro;
        }
    }
    let ocean_root = find_compress(&mut parent, ocean_idx);

    let mut edges: Vec<BoundaryEdge> = Vec::new();
    for i in 0..n {
        let cr_i = catchment_root[i];
        let row_offset = i * k;
        for c in 0..k {
            let idx = row_offset + c;
            let j = neighbor_idx[idx] as usize;
            if cr_i != catchment_root[j] {
                let w = match interface_pass {
                    Some(pass) => pass[idx],
                    None => max_endpoint_elev(elevation[i], elevation[j]),
                };
                edges.push(BoundaryEdge {
                    a: i as u32,
                    b: j as u32,
                    w,
                });
            }
        }
    }

    edges.sort_by(|x, y| cmp_f64_numpy_stable(x.w, y.w));

    let mut ev_kind: Vec<i8> = Vec::new();
    let mut ev_a: Vec<i64> = Vec::new();
    let mut ev_b: Vec<i64> = Vec::new();
    let mut ev_w: Vec<f64> = Vec::new();
    let mut ev_ra: Vec<i64> = Vec::new();
    let mut ev_rb: Vec<i64> = Vec::new();

    for edge in edges {
        let a = edge.a as usize;
        let b = edge.b as usize;
        let w = edge.w;
        let ra = find_compress(&mut parent, a);
        let rb = find_compress(&mut parent, b);
        if ra == rb {
            continue;
        }
        let a_ocean = ra == ocean_root;
        let b_ocean = rb == ocean_root;
        if a_ocean && b_ocean {
            continue;
        }

        if a_ocean || b_ocean {
            let land_root = if a_ocean { rb } else { ra };
            let land_node = if a_ocean { b } else { a };
            let ocean_node = if a_ocean { a } else { b };
            ev_kind.push(0);
            ev_a.push(land_node as i64);
            ev_b.push(ocean_node as i64);
            ev_w.push(w);
            ev_ra.push(land_root as i64);
            ev_rb.push(ocean_root as i64);
            parent[land_root] = ocean_root;
            continue;
        }

        ev_kind.push(1);
        ev_a.push(a as i64);
        ev_b.push(b as i64);
        ev_w.push(w);
        ev_ra.push(ra as i64);
        ev_rb.push(rb as i64);
        parent[ra] = rb;
    }

    Ok(HierarchyEvents {
        catchment_root,
        ev_kind,
        ev_a,
        ev_b,
        ev_w,
        ev_ra,
        ev_rb,
    })
}

#[pyfunction]
pub fn lake_components_rs<'py>(
    py: Python<'py>,
    is_lake: PyReadonlyArray1<'py, bool>,
    neighbor_idx: PyReadonlyArray2<'py, i64>,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<i64>>)> {
    let nb_shape = neighbor_idx.shape();
    let n = is_lake.len();
    if nb_shape[0] != n {
        return Err(PyValueError::new_err(format!(
            "neighbor_idx first dimension {} does not match is_lake length {n}",
            nb_shape[0]
        )));
    }
    let k = nb_shape[1];
    let lake_vec = copy_1d(&is_lake)?;
    if !lake_vec.iter().any(|&b| b) {
        return Ok((
            PyArray1::from_vec(py, vec![0_i64]),
            PyArray1::from_vec(py, Vec::<i64>::new()),
        ));
    }
    let nb_vec = copy_2d_row_major(&neighbor_idx)?;
    let (indptr, members) = py
        .detach(|| lake_components_core(&lake_vec, &nb_vec, k))
        .map_err(PyValueError::new_err)?;
    Ok((
        PyArray1::from_vec(py, indptr),
        PyArray1::from_vec(py, members),
    ))
}

#[pyfunction]
pub fn catchment_roots_rs<'py>(
    py: Python<'py>,
    elevation: PyReadonlyArray1<'py, f64>,
    is_ocean: PyReadonlyArray1<'py, bool>,
    neighbor_idx: PyReadonlyArray2<'py, i64>,
    ocean_catchment: i64,
) -> PyResult<Bound<'py, PyArray1<i64>>> {
    let n = elevation.len();
    if is_ocean.len() != n {
        return Err(PyValueError::new_err(format!(
            "is_ocean length {} does not match elevation length {n}",
            is_ocean.len()
        )));
    }
    let nb_shape = neighbor_idx.shape();
    if nb_shape[0] != n {
        return Err(PyValueError::new_err(format!(
            "neighbor_idx first dimension {} does not match elevation length {n}",
            nb_shape[0]
        )));
    }
    let k = nb_shape[1];
    let elev_vec = copy_1d(&elevation)?;
    let ocean_vec = copy_1d(&is_ocean)?;
    let nb_vec = copy_2d_row_major(&neighbor_idx)?;
    let roots = py
        .detach(|| catchment_roots_core(&elev_vec, &ocean_vec, &nb_vec, k, ocean_catchment))
        .map_err(PyValueError::new_err)?;
    Ok(PyArray1::from_vec(py, roots))
}

#[pyfunction]
#[pyo3(signature = (elevation, is_ocean, neighbor_idx, interface_pass=None, ocean_catchment=-2))]
pub fn lake_hierarchy_events_rs<'py>(
    py: Python<'py>,
    elevation: PyReadonlyArray1<'py, f64>,
    is_ocean: PyReadonlyArray1<'py, bool>,
    neighbor_idx: PyReadonlyArray2<'py, i64>,
    interface_pass: Option<PyReadonlyArray2<'py, f64>>,
    ocean_catchment: i64,
) -> PyResult<(
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i8>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<f64>>,
    Bound<'py, PyArray1<i64>>,
    Bound<'py, PyArray1<i64>>,
)> {
    let n = elevation.len();
    if is_ocean.len() != n {
        return Err(PyValueError::new_err(format!(
            "is_ocean length {} does not match elevation length {n}",
            is_ocean.len()
        )));
    }
    let nb_shape = neighbor_idx.shape();
    if nb_shape[0] != n {
        return Err(PyValueError::new_err(format!(
            "neighbor_idx first dimension {} does not match elevation length {n}",
            nb_shape[0]
        )));
    }
    let k = nb_shape[1];
    let pass_vec = match &interface_pass {
        Some(p) => {
            let ps = p.shape();
            if ps != nb_shape {
                return Err(PyValueError::new_err(
                    "interface_pass shape must match neighbor_idx",
                ));
            }
            Some(copy_2d_row_major(p)?)
        }
        None => None,
    };
    let elev_vec = copy_1d(&elevation)?;
    let ocean_vec = copy_1d(&is_ocean)?;
    let nb_vec = copy_2d_row_major(&neighbor_idx)?;

    let ev = py
        .detach(|| {
            lake_hierarchy_events_core(
                &elev_vec,
                &ocean_vec,
                &nb_vec,
                k,
                pass_vec.as_deref(),
                ocean_catchment,
            )
        })
        .map_err(PyValueError::new_err)?;

    Ok((
        PyArray1::from_vec(py, ev.catchment_root),
        PyArray1::from_vec(py, ev.ev_kind),
        PyArray1::from_vec(py, ev.ev_a),
        PyArray1::from_vec(py, ev.ev_b),
        PyArray1::from_vec(py, ev.ev_w),
        PyArray1::from_vec(py, ev.ev_ra),
        PyArray1::from_vec(py, ev.ev_rb),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lake_components_undirected_and_ordered() {
        // 0->1 (one-directional), 2 isolated, 3->4 (4->3), 5 non-lake
        let is_lake = vec![true, true, true, true, true, false];
        let nb = vec![
            1, 5, // 0 -> 1, 5
            5, 5, // 1 -> 5 (not back to 0)
            5, 5, // 2 -> 5
            4, 5, // 3 -> 4
            3, 5, // 4 -> 3
            0, 1, // 5 -> 0, 1
        ];
        let (indptr, members) = lake_components_core(&is_lake, &nb, 2).unwrap();
        assert_eq!(indptr, vec![0, 2, 3, 5]);
        assert_eq!(members, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_lake_components_all_false_and_empty() {
        let (ip0, m0) = lake_components_core(&[], &[], 4).unwrap();
        assert_eq!(ip0, vec![0]);
        assert!(m0.is_empty());

        let (ip1, m1) = lake_components_core(&[false, false], &[1, 0], 1).unwrap();
        assert_eq!(ip1, vec![0]);
        assert!(m1.is_empty());
    }

    #[test]
    fn test_catchment_roots_ties_and_ocean() {
        // Node 0 (30) -> min, Node 1 (50) -> ties between 0 (30) and 2 (30): picks col 0 (node 0).
        // Node 2 (30) -> 3 (-10 ocean) -> -2.
        let elev = vec![30.0, 50.0, 30.0, -10.0];
        let is_ocean = vec![false, false, false, true];
        let nb = vec![1, 1, 0, 2, 1, 3, 2, 2];
        let roots = catchment_roots_core(&elev, &is_ocean, &nb, 2, -2).unwrap();
        assert_eq!(roots, vec![0, 0, -2, -2]);
    }

    #[test]
    fn test_lake_hierarchy_events_merge_and_ocean() {
        // Two depressions: {0, 1} (floor 2.0) and {2, 3} (floor 4.0), saddle 1<->2 at 12.0,
        // rim 3<->4 at 30.0, node 4 drains to ocean 5.
        let elev = vec![2.0, 12.0, 9.0, 4.0, 30.0, -10.0];
        let is_ocean = vec![false, false, false, false, false, true];
        let nb = vec![1, 1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 4];
        let ev = lake_hierarchy_events_core(&elev, &is_ocean, &nb, 2, None, -2).unwrap();
        assert_eq!(ev.catchment_root, vec![0, 0, 3, 3, -2, -2]);
        assert_eq!(ev.ev_kind, vec![1, 0]);
        assert_eq!(ev.ev_a, vec![1, 3]);
        assert_eq!(ev.ev_b, vec![2, 4]);
        assert_eq!(ev.ev_w, vec![12.0, 30.0]);
        assert_eq!(ev.ev_ra, vec![0, 3]);
        assert_eq!(ev.ev_rb, vec![3, 6]);
    }

    #[test]
    fn test_signed_zero_endpoint_max_prefers_positive_zero() {
        let elev = vec![-0.0, 0.0];
        let is_ocean = vec![false, true];
        let nb = vec![1, 0];
        let ev = lake_hierarchy_events_core(&elev, &is_ocean, &nb, 1, None, -2).unwrap();
        assert_eq!(ev.ev_w.len(), 1);
        assert_eq!(ev.ev_w[0].to_bits(), 0.0_f64.to_bits());
    }
}
