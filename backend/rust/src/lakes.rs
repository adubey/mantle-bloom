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
    if n == 0 {
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
    if k == 0 {
        return Err("neighbor_idx must have at least one column when n > 0".to_string());
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

#[inline]
fn max_zero_numpy(v: f64) -> f64 {
    if v.is_nan() {
        f64::NAN
    } else if v > 0.0 {
        v
    } else {
        0.0
    }
}

#[allow(clippy::too_many_arguments)]
pub fn integrate_body_volume_core(
    available: f64,
    members: &[i64],
    cap_level: f64,
    dry_level: f64,
    dt_step: f64,
    elevation: &[f64],
    relief_m: &[f64],
    area_m2: &[f64],
    evaporation: &[f64],
    sublimation: Option<&[f64]>,
    shallow_evaporation_boost: f64,
    max_iterations: usize,
    newton_iterations: usize,
    storage_tolerance: f64,
) -> Result<(f64, f64, f64, f64), String> {
    let n = elevation.len();
    if relief_m.len() != n || area_m2.len() != n || evaporation.len() != n {
        return Err("terrain/evaporation arrays must have equal length".to_string());
    }
    if let Some(sub) = sublimation {
        if sub.len() != n {
            return Err("sublimation length must match elevation length".to_string());
        }
    }
    let m = members.len();
    let mut z = Vec::with_capacity(m);
    let mut r = Vec::with_capacity(m);
    let mut a = Vec::with_capacity(m);
    let mut ev = Vec::with_capacity(m);
    let mut su = Vec::with_capacity(m);
    let mut area_sum = 0.0_f64;

    for &idx_raw in members {
        if idx_raw < 0 || (idx_raw as usize) >= n {
            return Err(format!("member index {idx_raw} out of range [0, {n})"));
        }
        let idx = idx_raw as usize;
        z.push(elevation[idx]);
        r.push(relief_m[idx]);
        let area_i = area_m2[idx];
        a.push(area_i);
        area_sum += area_i;
        ev.push(evaporation[idx]);
        su.push(match sublimation {
            Some(sub) => max_zero_numpy(sub[idx]),
            None => 0.0,
        });
    }

    let budget = |level: f64| -> (f64, f64, f64, f64) {
        let mut volume = 0.0_f64;
        let mut evap_rate_sum = 0.0_f64;
        let mut wet_area = 0.0_f64;
        let mut evap_slope = 0.0_f64;

        for i in 0..m {
            let zi = z[i];
            let ri = r[i];
            let ai = a[i];
            let evi = ev[i];
            let sui = su[i];
            let full = level >= zi + ri;
            let depth = if full {
                max_zero_numpy(level - zi)
            } else if ri > 0.0 {
                let above_low = max_zero_numpy(level - (zi - ri));
                (above_low * above_low) / (4.0 * ri)
            } else {
                0.0
            };
            let wet = if full {
                if level > zi - ri { 1.0 } else { 0.0 }
            } else if ri > 0.0 {
                ((level - (zi - ri)) / (2.0 * ri)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            volume += depth * ai;
            evap_rate_sum += (ai * wet) * (evi * (1.0 + shallow_evaporation_boost * (1.0 - wet)) + sui);
            wet_area += ai * wet;
            if ri > 0.0 && wet > 0.0 && wet < 1.0 {
                let d_wet = 1.0 / (2.0 * ri);
                evap_slope += ai
                    * (evi * (1.0 + shallow_evaporation_boost - 2.0 * shallow_evaporation_boost * wet) + sui)
                    * d_wet;
            }
        }
        (
            volume,
            evap_rate_sum * dt_step,
            wet_area + evap_slope * dt_step,
            wet_area,
        )
    };

    let (cap_vol_raw, cap_evaporation, _, _) = budget(cap_level);
    let cap_volume = if cap_level.is_finite() { cap_vol_raw } else { f64::INFINITY };
    let surplus = available - cap_volume - cap_evaporation;
    if surplus > 0.0 {
        return Ok((cap_volume, cap_evaporation, surplus, cap_level));
    }

    let mut low = dry_level;
    let mut high = cap_level;
    let denom = if area_sum > 1.0 { area_sum } else { 1.0 };
    let mut level = dry_level + available / denom;
    if level < low {
        level = low;
    }
    if level > high {
        level = high;
    }
    let (mut volume, mut evaporated, mut slope, mut wet_area) = budget(level);

    let max_abs = dry_level
        .abs()
        .max(if cap_level.is_finite() { cap_level.abs() } else { 0.0 })
        .max(1.0);
    let level_floor = 8.0 * f64::EPSILON * max_abs;

    for iteration in 0..max_iterations {
        let residual = volume + evaporated - available;
        let level_error = if slope > 0.0 { residual.abs() / slope } else { f64::INFINITY };
        let vol_scale = if volume > 1.0 { volume } else { 1.0 };
        if level_error * wet_area <= storage_tolerance * vol_scale
            || level_error <= level_floor
            || high - low <= level_floor
        {
            break;
        }
        if residual > 0.0 {
            high = level;
        } else {
            low = level;
        }
        let candidate = if slope > 0.0 && iteration < newton_iterations {
            level - residual / slope
        } else {
            f64::NAN
        };
        level = if candidate > low && candidate < high {
            candidate
        } else {
            0.5 * (low + high)
        };
        (volume, evaporated, slope, wet_area) = budget(level);
    }

    let loss = (available - volume).max(0.0).min(available);
    Ok((volume, loss, 0.0, level))
}

#[pyfunction]
#[pyo3(signature = (
    available,
    members,
    cap_level,
    dry_level,
    dt_step,
    elevation,
    relief_m,
    area_m2,
    evaporation,
    sublimation=None,
    shallow_evaporation_boost=0.3,
    max_iterations=80,
    newton_iterations=12,
    storage_tolerance=1e-13,
))]
#[allow(clippy::too_many_arguments)]
pub fn integrate_body_volume_rs<'py>(
    available: f64,
    members: PyReadonlyArray1<'py, i64>,
    cap_level: f64,
    dry_level: f64,
    dt_step: f64,
    elevation: PyReadonlyArray1<'py, f64>,
    relief_m: PyReadonlyArray1<'py, f64>,
    area_m2: PyReadonlyArray1<'py, f64>,
    evaporation: PyReadonlyArray1<'py, f64>,
    sublimation: Option<PyReadonlyArray1<'py, f64>>,
    shallow_evaporation_boost: f64,
    max_iterations: usize,
    newton_iterations: usize,
    storage_tolerance: f64,
) -> PyResult<(f64, f64, f64, f64)> {
    let mem_vec = copy_1d(&members)?;
    let elev_slice = elevation.as_slice()?;
    let relief_slice = relief_m.as_slice()?;
    let area_slice = area_m2.as_slice()?;
    let evap_slice = evaporation.as_slice()?;
    let sub_slice = match &sublimation {
        Some(s) => Some(s.as_slice()?),
        None => None,
    };
    integrate_body_volume_core(
        available,
        &mem_vec,
        cap_level,
        dry_level,
        dt_step,
        elev_slice,
        relief_slice,
        area_slice,
        evap_slice,
        sub_slice,
        shallow_evaporation_boost,
        max_iterations,
        newton_iterations,
        storage_tolerance,
    )
    .map_err(PyValueError::new_err)
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
