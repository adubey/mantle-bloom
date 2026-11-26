use std::collections::VecDeque;

use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;

pub fn build_flood_adjacency_from_neighbor_idx(neighbor_idx: &[i64], n: usize, k: usize) -> (Vec<u32>, Vec<u32>) {
    let mut degree = vec![0_u32; n + 1];
    for i in 0..n {
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb in row {
            if nb >= 0 {
                let u = nb as usize;
                degree[i + 1] += 1;
                degree[u + 1] += 1;
            }
        }
    }
    let mut indptr = vec![0_u32; n + 1];
    let mut acc = 0_u32;
    for i in 0..=n {
        acc += degree[i];
        indptr[i] = acc;
    }
    let mut fill = indptr[..n].to_vec();
    let mut adjacency = vec![0_u32; indptr[n] as usize];
    for i in 0..n {
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb in row {
            if nb >= 0 {
                let u = nb as usize;
                adjacency[fill[i] as usize] = u as u32;
                fill[i] += 1;
                adjacency[fill[u] as usize] = i as u32;
                fill[u] += 1;
            }
        }
    }
    (indptr, adjacency)
}

#[pyclass]
pub struct MarchGraph {
    pub n: usize,
    pub indptr: Vec<u32>,
    pub indices: Vec<u32>,
    pub lengths: Vec<f64>,
    pub step_m: f64,
    pub flood_indptr: Option<Vec<u32>>,
    pub flood_edges: Option<Vec<u32>>,
    pub neighbor_idx: Option<Vec<i32>>,
    pub neighbor_k: usize,
}

#[pymethods]
impl MarchGraph {
    #[new]
    #[pyo3(signature = (indptr, indices, data, step_m, neighbor_idx=None))]
    pub fn new(
        indptr: PyReadonlyArray1<'_, i64>,
        indices: PyReadonlyArray1<'_, i64>,
        data: PyReadonlyArray1<'_, f64>,
        step_m: f64,
        neighbor_idx: Option<PyReadonlyArray2<'_, i64>>,
    ) -> PyResult<Self> {
        let indptr_slice = indptr.as_slice()?;
        let indices_slice = indices.as_slice()?;
        let data_slice = data.as_slice()?;
        let n = indptr_slice.len().saturating_sub(1);
        if n > (u32::MAX as usize) || indices_slice.len() > (u32::MAX as usize) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "graph exceeds u32 index range",
            ));
        }

        let indptr_vec: Vec<u32> = indptr_slice.iter().map(|&v| v as u32).collect();
        let indices_vec: Vec<u32> = indices_slice.iter().map(|&v| v as u32).collect();
        let lengths_vec: Vec<f64> = data_slice.to_vec();

        let (flood_indptr, flood_edges, neighbor_vec, neighbor_k) = if let Some(nb_arr) = neighbor_idx {
            let shape = nb_arr.shape();
            let nb_n = shape[0];
            let k = shape[1];
            let nb_slice = nb_arr.as_slice()?;
            let (fi, fe) = build_flood_adjacency_from_neighbor_idx(nb_slice, nb_n, k);
            let nb_i32: Vec<i32> = nb_slice.iter().map(|&v| v as i32).collect();
            (Some(fi), Some(fe), Some(nb_i32), k)
        } else {
            (None, None, None, 0)
        };

        Ok(Self {
            n,
            indptr: indptr_vec,
            indices: indices_vec,
            lengths: lengths_vec,
            step_m,
            flood_indptr,
            flood_edges,
            neighbor_idx: neighbor_vec,
            neighbor_k,
        })
    }
}

/// Connected components over CSR edges restricted to `members`, numbered in ascending order of
/// each component's lowest node index (`-1` for non-members). Returns `(comp, num_bodies)`.
pub fn connected_components_masked(
    n: usize,
    indptr: &[u32],
    indices: &[u32],
    members: &[bool],
) -> (Vec<i32>, usize) {
    let mut comp = vec![-1_i32; n];
    let mut queue = VecDeque::new();
    let mut next_id = 0_i32;
    for start in 0..n {
        if !members[start] || comp[start] >= 0 {
            continue;
        }
        comp[start] = next_id;
        queue.push_back(start as u32);
        while let Some(u) = queue.pop_front() {
            let uu = u as usize;
            let start_e = indptr[uu] as usize;
            let end_e = indptr[uu + 1] as usize;
            for &v in &indices[start_e..end_e] {
                let vu = v as usize;
                if members[vu] && comp[vu] < 0 {
                    comp[vu] = next_id;
                    queue.push_back(v);
                }
            }
        }
        next_id += 1;
    }
    (comp, next_id as usize)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_connected_components_lowest_node_order() {
        let indptr = vec![0, 1, 2, 3, 4];
        let indices = vec![2, 3, 0, 1];
        let members = vec![false, true, true, true];
        let (comp, nb) = connected_components_masked(4, &indptr, &indices, &members);
        assert_eq!(nb, 2);
        assert_eq!(comp, vec![-1, 0, 1, 0]);
    }
}
