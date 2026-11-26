use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;

use crate::climate::{clip_min, max_nan, min_nan, RHO_ICE, RHO_WATER};
use crate::graph::{build_flood_adjacency_from_neighbor_idx, MarchGraph};

const SOLVE_TOLERANCE_M: f64 = 1.0;
const SOLVE_MAX_ITERS: usize = 80;
const ENCLOSED_OCEAN_KEEP_OVERLAP_FRACTION: f64 = 0.5;

#[derive(Copy, Clone, Debug)]
struct FloodEntry {
    level: f64,
    node: u32,
}

impl PartialEq for FloodEntry {
    fn eq(&self, other: &Self) -> bool {
        self.level == other.level && self.node == other.node
    }
}

impl Eq for FloodEntry {}

impl Ord for FloodEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .level
            .total_cmp(&self.level)
            .then_with(|| other.node.cmp(&self.node))
    }
}

impl PartialOrd for FloodEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Priority-flood sea-level thresholds per node from `seed` over CSR undirected `adjacency`,
/// optionally pruned at `max_level` (nodes needing a higher threshold stay `INFINITY`).
pub fn seeded_flood_levels_bounded(
    indptr: &[u32],
    adjacency: &[u32],
    elevations: &[f64],
    seed: &[bool],
    max_level: Option<f64>,
) -> Vec<f64> {
    let n = elevations.len();
    let mut level = vec![f64::INFINITY; n];
    let mut heap = BinaryHeap::new();
    let mut plain_queue = VecDeque::new();

    for i in 0..n {
        if seed[i] && !elevations[i].is_nan() {
            let z = elevations[i];
            if max_level.map_or(true, |lim| z <= lim) {
                level[i] = z;
            }
        }
    }

    for i in 0..n {
        if !seed[i] || level[i].is_infinite() {
            continue;
        }
        let reached = level[i];
        let start = indptr[i] as usize;
        let end = indptr[i + 1] as usize;
        for &nb in &adjacency[start..end] {
            let nbu = nb as usize;
            if seed[nbu] {
                continue;
            }
            let z = elevations[nbu];
            if z.is_nan() {
                continue;
            }
            let cand = if reached > z { reached } else { z };
            if max_level.map_or(true, |lim| cand <= lim) && cand < level[nbu] {
                level[nbu] = cand;
                heap.push(FloodEntry {
                    level: cand,
                    node: nb,
                });
            }
        }
    }

    while let Some(FloodEntry {
        level: reached,
        node: start_node,
    }) = heap.pop()
    {
        if reached > level[start_node as usize] {
            continue;
        }
        plain_queue.push_back(start_node);
        while let Some(iu) = plain_queue.pop_front() {
            let i = iu as usize;
            let start = indptr[i] as usize;
            let end = indptr[i + 1] as usize;
            for &nb in &adjacency[start..end] {
                let nbu = nb as usize;
                let z = elevations[nbu];
                if z.is_nan() {
                    continue;
                }
                let cand = if reached > z { reached } else { z };
                if max_level.map_or(true, |lim| cand <= lim) && cand < level[nbu] {
                    level[nbu] = cand;
                    if z <= reached {
                        plain_queue.push_back(nb);
                    } else {
                        heap.push(FloodEntry {
                            level: cand,
                            node: nb,
                        });
                    }
                }
            }
        }
    }

    level
}

pub fn seeded_flood_levels_core(
    indptr: &[u32],
    adjacency: &[u32],
    elevations: &[f64],
    seed: &[bool],
) -> Vec<f64> {
    seeded_flood_levels_bounded(indptr, adjacency, elevations, seed, None)
}

#[inline]
fn volume_at_level(levels: &[f64], floor: &[f64], areas: &[f64], h: f64) -> f64 {
    let mut total = 0.0_f64;
    for i in 0..levels.len() {
        if levels[i] <= h {
            let d = h - floor[i];
            if d > 0.0 {
                total += areas[i] * d;
            }
        }
    }
    total
}

#[inline]
fn volume_at_bracketed(
    v_blo_always: f64,
    area_always: f64,
    blo: f64,
    coastal: &[(f64, f64, f64)],
    h: f64,
) -> f64 {
    let mut total = v_blo_always + area_always * (h - blo);
    for &(lvl, flr, area) in coastal {
        if lvl <= h {
            let d = h - flr;
            if d > 0.0 {
                total += area * d;
            }
        }
    }
    total
}

#[allow(clippy::too_many_arguments)]
pub fn solve_sea_level_connected_core(
    indptr: &[u32],
    adjacency: &[u32],
    elevations: &[f64],
    areas: &[f64],
    ice_water_m: Option<&[f64]>,
    water_volume_m3: f64,
    seed: &[bool],
    bracket: Option<(f64, f64)>,
) -> f64 {
    let n = elevations.len();
    if n == 0 {
        return 0.0;
    }
    let mut min_z = f64::INFINITY;
    let mut max_z = f64::NEG_INFINITY;
    let mut total_area = 0.0_f64;
    for i in 0..n {
        let z = elevations[i];
        min_z = min_nan(min_z, z);
        max_z = max_nan(max_z, z);
        total_area += areas[i];
    }
    if water_volume_m3 <= 0.0 || total_area <= 0.0 {
        return min_z;
    }

    let floor: Vec<f64> = match ice_water_m {
        Some(iw) => elevations
            .iter()
            .zip(iw.iter())
            .map(|(&z, &w)| z + clip_min(w, 0.0))
            .collect(),
        None => elevations.to_vec(),
    };
    let tolerance_m3 = SOLVE_TOLERANCE_M * total_area / (n as f64);

    if let Some((blo, bhi)) = bracket {
        if blo < bhi {
            let bounded = seeded_flood_levels_bounded(indptr, adjacency, elevations, seed, Some(bhi));
            let mut v_blo_always = 0.0_f64;
            let mut area_always = 0.0_f64;
            let mut coastal: Vec<(f64, f64, f64)> = Vec::new();
            for i in 0..n {
                let lvl = bounded[i];
                if lvl <= bhi {
                    let flr = floor[i];
                    let a = areas[i];
                    if lvl <= blo && flr <= blo {
                        v_blo_always += a * (blo - flr);
                        area_always += a;
                    } else {
                        coastal.push((lvl, flr, a));
                    }
                }
            }
            let v_blo = volume_at_bracketed(v_blo_always, area_always, blo, &coastal, blo);
            let v_bhi = volume_at_bracketed(v_blo_always, area_always, blo, &coastal, bhi);
            if v_blo <= water_volume_m3 && water_volume_m3 <= v_bhi {
                let mut lo = blo;
                let mut hi = bhi;
                for _ in 0..SOLVE_MAX_ITERS {
                    let mid = 0.5 * (lo + hi);
                    let v = volume_at_bracketed(v_blo_always, area_always, blo, &coastal, mid);
                    if (v - water_volume_m3).abs() < tolerance_m3 {
                        return mid;
                    }
                    if v < water_volume_m3 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                return 0.5 * (lo + hi);
            }
        }
    }

    let levels = seeded_flood_levels_core(indptr, adjacency, elevations, seed);
    let mut lo = min_z;
    let mut hi = max_z + water_volume_m3 / total_area + 1.0;

    for _ in 0..SOLVE_MAX_ITERS {
        let mid = 0.5 * (lo + hi);
        let v = volume_at_level(&levels, &floor, areas, mid);
        if (v - water_volume_m3).abs() < tolerance_m3 {
            return mid;
        }
        if v < water_volume_m3 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[inline]
fn uf_find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let p = parent[x as usize];
        parent[x as usize] = parent[p as usize];
        x = parent[x as usize];
    }
    x
}

pub fn below_sea_components_core(below: &[bool], neighbor_idx: &[i32], n: usize, k: usize) -> Vec<i32> {
    let mut parent: Vec<u32> = (0..n as u32).collect();
    for i in 0..n {
        if !below[i] {
            continue;
        }
        let mut ri = uf_find(&mut parent, i as u32);
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb in row {
            if nb < 0 {
                continue;
            }
            let u = nb as usize;
            if !below[u] {
                continue;
            }
            let rj = uf_find(&mut parent, u as u32);
            if ri != rj {
                parent[ri as usize] = rj;
                ri = rj;
            }
        }
    }
    let mut labels = vec![0_i32; n];
    let mut root_label = vec![-1_i32; n];
    let mut count = 0_i32;
    for i in 0..n {
        let r = uf_find(&mut parent, i as u32) as usize;
        if root_label[r] < 0 {
            root_label[r] = count;
            count += 1;
        }
        labels[i] = root_label[r];
    }
    labels
}

#[allow(clippy::too_many_arguments)]
pub fn connected_ocean_mask_core(
    elevation: &[f64],
    sea_level_m: f64,
    neighbor_idx: &[i32],
    k: usize,
    areas_m2: Option<&[f64]>,
    prior_ocean: Option<&[bool]>,
    enclosed_min_area_km2: f64,
    enclosed_keep_area_km2: f64,
    planet_radius_km: f64,
) -> Vec<bool> {
    let n = elevation.len();
    let mut below = vec![false; n];
    let mut any_below = false;
    let mut all_below = true;
    for i in 0..n {
        let b = elevation[i] <= sea_level_m;
        below[i] = b;
        if b {
            any_below = true;
        } else {
            all_below = false;
        }
    }
    if n == 0 || !any_below || all_below {
        return below;
    }

    let mut parent: Vec<u32> = (0..n as u32).collect();
    for i in 0..n {
        if !below[i] {
            continue;
        }
        let mut ri = uf_find(&mut parent, i as u32);
        let row = &neighbor_idx[i * k..(i + 1) * k];
        for &nb in row {
            if nb < 0 {
                continue;
            }
            let u = nb as usize;
            if !below[u] {
                continue;
            }
            let rj = uf_find(&mut parent, u as u32);
            if ri != rj {
                parent[ri as usize] = rj;
                ri = rj;
            }
        }
    }

    // Quick check for single below-sea-level component before allocating component maps.
    let mut first_root = u32::MAX;
    let mut single_comp = true;
    for i in 0..n {
        if !below[i] {
            continue;
        }
        let r = uf_find(&mut parent, i as u32);
        if first_root == u32::MAX {
            first_root = r;
        } else if r != first_root {
            single_comp = false;
            break;
        }
    }
    if single_comp {
        return below;
    }

    let mut root_to_comp = vec![-1_i32; n];
    let mut node_comp = vec![-1_i32; n];
    let mut n_comps = 0_usize;
    let default_area = 4.0 * std::f64::consts::PI * (planet_radius_km * 1000.0).powi(2) / (n as f64);

    let mut comp_area: Vec<f64> = Vec::new();
    let mut comp_prior_area: Vec<f64> = Vec::new();

    for i in 0..n {
        if !below[i] {
            continue;
        }
        let r = uf_find(&mut parent, i as u32) as usize;
        let c = if root_to_comp[r] < 0 {
            let idx = n_comps as i32;
            root_to_comp[r] = idx;
            comp_area.push(0.0);
            comp_prior_area.push(0.0);
            n_comps += 1;
            idx
        } else {
            root_to_comp[r]
        };
        node_comp[i] = c;
        let cu = c as usize;
        let a = areas_m2.map_or(default_area, |arr| arr[i]);
        comp_area[cu] += a;
        if let Some(po) = prior_ocean {
            if po[i] {
                comp_prior_area[cu] += a;
            }
        }
    }

    let mut largest_idx = 0usize;
    let mut largest_area = comp_area[0];
    for c in 1..n_comps {
        if comp_area[c] > largest_area {
            largest_area = comp_area[c];
            largest_idx = c;
        }
    }

    let mut qualifies = vec![false; n_comps];
    for c in 0..n_comps {
        let was_ocean = prior_ocean.is_some()
            && comp_prior_area[c] >= ENCLOSED_OCEAN_KEEP_OVERLAP_FRACTION * comp_area[c];
        let threshold_m2 = if was_ocean {
            enclosed_keep_area_km2
        } else {
            enclosed_min_area_km2
        } * 1e6;
        if comp_area[c] >= threshold_m2 {
            qualifies[c] = true;
        }
    }
    qualifies[largest_idx] = true;

    let mut result = vec![false; n];
    for i in 0..n {
        if below[i] {
            let c = node_comp[i] as usize;
            if qualifies[c] {
                result[i] = true;
            }
        }
    }
    result
}

#[pyfunction]
#[pyo3(signature = (elevations, seed, indptr=None, adjacency=None, neighbor_idx=None, rust_graph=None))]
pub fn seeded_flood_levels<'py>(
    py: Python<'py>,
    elevations: PyReadonlyArray1<'py, f64>,
    seed: PyReadonlyArray1<'py, bool>,
    indptr: Option<PyReadonlyArray1<'py, i64>>,
    adjacency: Option<PyReadonlyArray1<'py, i64>>,
    neighbor_idx: Option<PyReadonlyArray2<'py, i64>>,
    rust_graph: Option<&MarchGraph>,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let elev = elevations.as_slice()?;
    let seed_slice = seed.as_slice()?;

    if let Some(rg) = rust_graph {
        if let (Some(fi), Some(fe)) = (&rg.flood_indptr, &rg.flood_edges) {
            let levels = py.detach(|| seeded_flood_levels_core(fi, fe, elev, seed_slice));
            return Ok(PyArray1::from_vec(py, levels));
        }
    }

    if let (Some(ip), Some(adj)) = (indptr, adjacency) {
        let ip_slice = ip.as_slice()?;
        let adj_slice = adj.as_slice()?;
        let levels = py.detach(|| {
            let ip_u: Vec<u32> = ip_slice.iter().map(|&v| v as u32).collect();
            let adj_u: Vec<u32> = adj_slice.iter().map(|&v| v as u32).collect();
            seeded_flood_levels_core(&ip_u, &adj_u, elev, seed_slice)
        });
        return Ok(PyArray1::from_vec(py, levels));
    }

    if let Some(nb) = neighbor_idx {
        let shape = nb.shape();
        let n = shape[0];
        let k = shape[1];
        let nb_slice = nb.as_slice()?;
        let levels = py.detach(|| {
            let (fi, fe) = build_flood_adjacency_from_neighbor_idx(nb_slice, n, k);
            seeded_flood_levels_core(&fi, &fe, elev, seed_slice)
        });
        return Ok(PyArray1::from_vec(py, levels));
    }

    Err(pyo3::exceptions::PyValueError::new_err(
        "seeded_flood_levels requires rust_graph, (indptr, adjacency), or neighbor_idx",
    ))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (elevations, areas_m2, water_volume_m3, seed, ice_water_m=None, bracket=None, indptr=None, adjacency=None, neighbor_idx=None, rust_graph=None))]
pub fn solve_sea_level_connected<'py>(
    py: Python<'py>,
    elevations: PyReadonlyArray1<'py, f64>,
    areas_m2: PyReadonlyArray1<'py, f64>,
    water_volume_m3: f64,
    seed: PyReadonlyArray1<'py, bool>,
    ice_water_m: Option<PyReadonlyArray1<'py, f64>>,
    bracket: Option<(f64, f64)>,
    indptr: Option<PyReadonlyArray1<'py, i64>>,
    adjacency: Option<PyReadonlyArray1<'py, i64>>,
    neighbor_idx: Option<PyReadonlyArray2<'py, i64>>,
    rust_graph: Option<&MarchGraph>,
) -> PyResult<f64> {
    let elev = elevations.as_slice()?;
    let areas = areas_m2.as_slice()?;
    let seed_slice = seed.as_slice()?;
    let iw_slice = match &ice_water_m {
        Some(arr) => Some(arr.as_slice()?),
        None => None,
    };

    if let Some(rg) = rust_graph {
        if let (Some(fi), Some(fe)) = (&rg.flood_indptr, &rg.flood_edges) {
            return Ok(py.detach(|| {
                solve_sea_level_connected_core(fi, fe, elev, areas, iw_slice, water_volume_m3, seed_slice, bracket)
            }));
        }
    }

    if let (Some(ip), Some(adj)) = (indptr, adjacency) {
        let ip_slice = ip.as_slice()?;
        let adj_slice = adj.as_slice()?;
        return Ok(py.detach(|| {
            let ip_u: Vec<u32> = ip_slice.iter().map(|&v| v as u32).collect();
            let adj_u: Vec<u32> = adj_slice.iter().map(|&v| v as u32).collect();
            solve_sea_level_connected_core(
                &ip_u,
                &adj_u,
                elev,
                areas,
                iw_slice,
                water_volume_m3,
                seed_slice,
                bracket,
            )
        }));
    }

    if let Some(nb) = neighbor_idx {
        let shape = nb.shape();
        let n = shape[0];
        let k = shape[1];
        let nb_slice = nb.as_slice()?;
        return Ok(py.detach(|| {
            let (fi, fe) = build_flood_adjacency_from_neighbor_idx(nb_slice, n, k);
            solve_sea_level_connected_core(
                &fi,
                &fe,
                elev,
                areas,
                iw_slice,
                water_volume_m3,
                seed_slice,
                bracket,
            )
        }));
    }

    Err(pyo3::exceptions::PyValueError::new_err(
        "solve_sea_level_connected requires rust_graph, (indptr, adjacency), or neighbor_idx",
    ))
}

#[pyfunction]
#[pyo3(signature = (below, neighbor_idx=None, rust_graph=None))]
pub fn below_sea_components<'py>(
    py: Python<'py>,
    below: PyReadonlyArray1<'py, bool>,
    neighbor_idx: Option<PyReadonlyArray2<'py, i64>>,
    rust_graph: Option<&MarchGraph>,
) -> PyResult<Bound<'py, PyArray1<i32>>> {
    let below_slice = below.as_slice()?;
    let labels = if let Some(rg) = rust_graph {
        if let Some(nb) = &rg.neighbor_idx {
            let n = below_slice.len();
            let k = rg.neighbor_k;
            py.detach(|| below_sea_components_core(below_slice, nb, n, k))
        } else {
            return Err(pyo3::exceptions::PyValueError::new_err("rust_graph missing neighbor_idx"));
        }
    } else if let Some(nb) = neighbor_idx {
        let shape = nb.shape();
        let nb_slice = nb.as_slice()?;
        py.detach(|| {
            let nb_i32: Vec<i32> = nb_slice.iter().map(|&v| v as i32).collect();
            below_sea_components_core(below_slice, &nb_i32, shape[0], shape[1])
        })
    } else {
        return Err(pyo3::exceptions::PyValueError::new_err("neighbor_idx or rust_graph required"));
    };
    Ok(PyArray1::from_vec(py, labels))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (elevation, sea_level_m, enclosed_min_area_km2, enclosed_keep_area_km2, planet_radius_km, neighbor_idx=None, areas_m2=None, prior_ocean=None, rust_graph=None))]
pub fn connected_ocean_mask<'py>(
    py: Python<'py>,
    elevation: PyReadonlyArray1<'py, f64>,
    sea_level_m: f64,
    enclosed_min_area_km2: f64,
    enclosed_keep_area_km2: f64,
    planet_radius_km: f64,
    neighbor_idx: Option<PyReadonlyArray2<'py, i64>>,
    areas_m2: Option<PyReadonlyArray1<'py, f64>>,
    prior_ocean: Option<PyReadonlyArray1<'py, bool>>,
    rust_graph: Option<&MarchGraph>,
) -> PyResult<Bound<'py, PyArray1<bool>>> {
    let elev_slice = elevation.as_slice()?;
    let areas_slice = match &areas_m2 {
        Some(a) => Some(a.as_slice()?),
        None => None,
    };
    let prior_slice = match &prior_ocean {
        Some(p) => Some(p.as_slice()?),
        None => None,
    };

    let mask = if let Some(rg) = rust_graph {
        if let Some(nb) = &rg.neighbor_idx {
            let k = rg.neighbor_k;
            py.detach(|| {
                connected_ocean_mask_core(
                    elev_slice,
                    sea_level_m,
                    nb,
                    k,
                    areas_slice,
                    prior_slice,
                    enclosed_min_area_km2,
                    enclosed_keep_area_km2,
                    planet_radius_km,
                )
            })
        } else {
            return Err(pyo3::exceptions::PyValueError::new_err("rust_graph missing neighbor_idx"));
        }
    } else if let Some(nb) = neighbor_idx {
        let shape = nb.shape();
        let nb_slice = nb.as_slice()?;
        py.detach(|| {
            let nb_i32: Vec<i32> = nb_slice.iter().map(|&v| v as i32).collect();
            connected_ocean_mask_core(
                elev_slice,
                sea_level_m,
                &nb_i32,
                shape[1],
                areas_slice,
                prior_slice,
                enclosed_min_area_km2,
                enclosed_keep_area_km2,
                planet_radius_km,
            )
        })
    } else {
        return Err(pyo3::exceptions::PyValueError::new_err("neighbor_idx or rust_graph required"));
    };

    let _ = (RHO_ICE, RHO_WATER);
    Ok(PyArray1::from_vec(py, mask))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_flood(indptr: &[u32], adj: &[u32], elev: &[f64], seed: &[bool]) -> Vec<f64> {
        let n = elev.len();
        let mut level = vec![f64::INFINITY; n];
        let mut heap = BinaryHeap::new();
        for i in 0..n {
            if seed[i] && !elev[i].is_nan() {
                level[i] = elev[i];
                heap.push(FloodEntry {
                    level: level[i],
                    node: i as u32,
                });
            }
        }
        while let Some(FloodEntry { level: reached, node: iu }) = heap.pop() {
            let i = iu as usize;
            if reached > level[i] {
                continue;
            }
            for &nb in &adj[indptr[i] as usize..indptr[i + 1] as usize] {
                let nbu = nb as usize;
                let z = elev[nbu];
                if z.is_nan() {
                    continue;
                }
                let cand = if reached > z { reached } else { z };
                if cand < level[nbu] {
                    level[nbu] = cand;
                    heap.push(FloodEntry { level: cand, node: nb });
                }
            }
        }
        level
    }

    #[test]
    fn test_flood_matches_naive_with_ties_and_nan() {
        let nb = vec![1, -1, 0, 2, 1, 3, 2, 4, 3, -1];
        let (ip, adj) = build_flood_adjacency_from_neighbor_idx(&nb, 5, 2);
        let elev = vec![-100.0, -50.0, -50.0, f64::NAN, -200.0];
        let seed = vec![true, false, false, false, false];
        let fast = seeded_flood_levels_core(&ip, &adj, &elev, &seed);
        let slow = naive_flood(&ip, &adj, &elev, &seed);
        assert_eq!(fast.len(), slow.len());
        for (a, b) in fast.iter().zip(slow.iter()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    #[test]
    fn test_bounded_warm_bracket_matches_unlimited_and_fallback() {
        let nb = vec![1, -1, 0, 2, 1, 3, 2, 4, 3, 5, 4, -1];
        let (ip, adj) = build_flood_adjacency_from_neighbor_idx(&nb, 6, 2);
        let elev = vec![-200.0, -150.0, -80.0, -10.0, 200.0, 500.0];
        let areas = vec![1e10; 6];
        let seed = vec![true, false, false, false, false, false];
        let ice_w = vec![0.0, 0.0, 20.0, 5.0, 0.0, 0.0];
        let budget = 2.0e12;

        let unbracketed = solve_sea_level_connected_core(&ip, &adj, &elev, &areas, Some(&ice_w), budget, &seed, None);
        let warm_hit = solve_sea_level_connected_core(
            &ip, &adj, &elev, &areas, Some(&ice_w), budget, &seed, Some((unbracketed - 25.0, unbracketed + 25.0)),
        );
        let stale_miss = solve_sea_level_connected_core(
            &ip, &adj, &elev, &areas, Some(&ice_w), budget, &seed, Some((-300.0, -250.0)),
        );
        assert!((warm_hit - unbracketed).abs() < 1e-10);
        assert_eq!(stale_miss.to_bits(), unbracketed.to_bits());
    }
}
