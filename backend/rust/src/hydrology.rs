use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Choose downhill channels and redirect filled sinks (#1041).
#[pyfunction]
pub fn flow_direction_rs<'py>(
    py: Python<'py>,
    elevation: PyReadonlyArray1<'py, f64>,
    ocean: PyReadonlyArray1<'py, bool>,
    neighbors: PyReadonlyArray2<'py, i64>,
    offset: PyReadonlyArray1<'py, f64>,
    filled: PyReadonlyArray1<'py, f64>,
    spill: PyReadonlyArray1<'py, i64>,
    channel: PyReadonlyArray1<'py, f64>,
    threshold: f64,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<bool>>)> {
    let n = elevation.len();
    let shape = neighbors.shape();
    if ocean.len() != n
        || offset.len() != n
        || filled.len() != n
        || spill.len() != n
        || channel.len() != n
        || shape[0] != n
        || shape[1] == 0
    {
        return Err(PyValueError::new_err("incompatible flow direction shapes"));
    }
    let k = shape[1];
    let elev = elevation.as_slice()?;
    let ocean = ocean.as_slice()?;
    let nb = neighbors.as_slice()?;
    let offset = offset.as_slice()?;
    let filled = filled.as_slice()?;
    let spill = spill.as_slice()?;
    let channel = channel.as_slice()?;
    if nb.iter().any(|&j| j < -(n as i64) || j >= n as i64) {
        return Err(PyValueError::new_err("neighbor index out of range"));
    }
    let (targets, spills) = py.detach(|| {
        let mut targets = vec![-1_i64; n];
        let mut spills = vec![false; n];
        for i in 0..n {
            if ocean[i] {
                continue;
            }
            let mut lowest = f64::INFINITY;
            let mut best_col = 0;
            let mut deepest = f64::NEG_INFINITY;
            let mut channel_col = 0;
            for col in 0..k {
                let raw = nb[i * k + col];
                let j = if raw < 0 {
                    (raw + n as i64) as usize
                } else {
                    raw as usize
                };
                if elev[j] < elev[i] {
                    if elev[j] < lowest {
                        lowest = elev[j];
                        best_col = col;
                    }
                    if channel[j] > threshold && channel[j] > deepest {
                        deepest = channel[j];
                        channel_col = col;
                    }
                }
            }
            if lowest.is_finite() {
                let col = if deepest > f64::NEG_INFINITY {
                    channel_col
                } else {
                    best_col
                };
                targets[i] = nb[i * k + col];
            }
            spills[i] = targets[i] < 0
                && !ocean[i]
                && (elev[i] + offset[i] >= filled[i] || filled[i] <= elev[i]);
            if spills[i] {
                targets[i] = spill[i];
            }
            if ocean[i] {
                targets[i] = -1;
            }
        }
        (targets, spills)
    });
    Ok((
        PyArray1::from_vec(py, targets),
        PyArray1::from_vec(py, spills),
    ))
}

/// Accumulate routed material in the exact order chosen by Python (#1046).
#[pyfunction]
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn route_downstream_rs<'py>(
    py: Python<'py>,
    order: PyReadonlyArray1<'py, i64>,
    stuck: PyReadonlyArray1<'py, i64>,
    target: PyReadonlyArray1<'py, i64>,
    ocean: PyReadonlyArray1<'py, bool>,
    source: PyReadonlyArray1<'py, f64>,
    retain: PyReadonlyArray1<'py, f64>,
    loss: PyReadonlyArray1<'py, f64>,
    erosion: Option<PyReadonlyArray1<'py, f64>>,
) -> PyResult<Vec<Bound<'py, PyArray1<f64>>>> {
    let n = source.len();
    if target.len() != n
        || ocean.len() != n
        || retain.len() != n
        || loss.len() != n
        || erosion.as_ref().is_some_and(|a| a.len() != n)
    {
        return Err(PyValueError::new_err("incompatible routing shapes"));
    }
    let order = order.as_slice()?;
    let stuck = stuck.as_slice()?;
    let target = target.as_slice()?;
    let ocean = ocean.as_slice()?;
    let source = source.as_slice()?;
    let retain = retain.as_slice()?;
    let loss = loss.as_slice()?;
    let erosion = erosion.as_ref().map(|a| a.as_slice()).transpose()?;
    if order.iter().chain(stuck).any(|&i| i < 0 || i as usize >= n)
        || target.iter().any(|&i| i >= n as i64)
    {
        return Err(PyValueError::new_err("routing index out of range"));
    }
    let arrays = py.detach(|| {
        let mut through: Vec<f64> = (0..n)
            .map(|i| if ocean[i] { 0.0 } else { source[i] })
            .collect();
        let mut erode: Vec<f64> = (0..n)
            .map(|i| {
                if ocean[i] {
                    0.0
                } else {
                    erosion.map_or(0.0, |e| e[i])
                }
            })
            .collect();
        let mut deposited = vec![0.0; n];
        let mut retained = vec![0.0; n];
        let mut sink = vec![0.0; n];
        let mut export = vec![0.0; n];
        let mut cycle = vec![0.0; n];
        let mut evaporated = vec![0.0; n];
        for (&raw, is_stuck) in order
            .iter()
            .map(|i| (i, false))
            .chain(stuck.iter().map(|i| (i, true)))
        {
            let i = raw as usize;
            let evaporated_here = through[i] * loss[i];
            through[i] -= evaporated_here;
            evaporated[i] += evaporated_here;
            if erosion.is_some() {
                erode[i] -= erode[i] * loss[i];
            }
            let retained_here = through[i] * retain[i];
            through[i] -= retained_here;
            deposited[i] += retained_here;
            retained[i] += retained_here;
            if is_stuck {
                deposited[i] += through[i];
                cycle[i] += through[i];
                if erosion.is_some() {
                    erode[i] = 0.0;
                }
            } else if target[i] < 0 {
                deposited[i] += through[i];
                sink[i] += through[i];
            } else {
                let j = target[i] as usize;
                if ocean[j] {
                    deposited[j] += through[i];
                    export[j] += through[i];
                } else {
                    through[j] += through[i];
                    if erosion.is_some() {
                        erode[j] += erode[i];
                    }
                }
            }
        }
        vec![
            through, deposited, retained, sink, export, cycle, evaporated, erode,
        ]
    });
    Ok(arrays
        .into_iter()
        .map(|a| PyArray1::from_vec(py, a))
        .collect())
}

/// Topologically order uphill routing edges, preferring highest ready nodes (#1046).
#[pyfunction]
pub fn routing_order_rs<'py>(
    py: Python<'py>,
    elevation: PyReadonlyArray1<'py, f64>,
    ocean: PyReadonlyArray1<'py, bool>,
    target: PyReadonlyArray1<'py, i64>,
) -> PyResult<(Vec<usize>, Vec<usize>)> {
    use std::cmp::Ordering;
    use std::collections::BinaryHeap;
    #[derive(PartialEq)]
    struct Ready(f64, usize);
    impl Eq for Ready {}
    impl PartialOrd for Ready {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    impl Ord for Ready {
        fn cmp(&self, other: &Self) -> Ordering {
            self.0
                .partial_cmp(&other.0)
                .unwrap_or(Ordering::Equal)
                .then_with(|| other.1.cmp(&self.1))
        }
    }
    let n = elevation.len();
    if ocean.len() != n || target.len() != n {
        return Err(PyValueError::new_err("incompatible routing shapes"));
    }
    let elevation = elevation.as_slice()?;
    let ocean = ocean.as_slice()?;
    let target = target.as_slice()?;
    if target.iter().any(|&i| i >= n as i64) {
        return Err(PyValueError::new_err("routing index out of range"));
    }
    Ok(py.detach(|| {
        let mut indegree = vec![0_usize; n];
        for i in 0..n {
            if !ocean[i] && target[i] >= 0 && !ocean[target[i] as usize] {
                indegree[target[i] as usize] += 1;
            }
        }
        let mut heap = BinaryHeap::new();
        for i in 0..n {
            if !ocean[i] && indegree[i] == 0 {
                heap.push(Ready(elevation[i], i));
            }
        }
        let mut order = Vec::new();
        while let Some(Ready(_, i)) = heap.pop() {
            order.push(i);
            let j = target[i];
            if j >= 0 && !ocean[j as usize] {
                let j = j as usize;
                indegree[j] -= 1;
                if indegree[j] == 0 {
                    heap.push(Ready(elevation[j], j));
                }
            }
        }
        let stuck = (0..n).filter(|&i| !ocean[i] && indegree[i] > 0).collect();
        (order, stuck)
    }))
}

/// Sweep linear river reservoirs in upstream-to-downstream order (#1046).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn sweep_river_reaches_rs<'py>(
    py: Python<'py>,
    order: PyReadonlyArray1<'py, i64>,
    target: PyReadonlyArray1<'py, i64>,
    ocean: PyReadonlyArray1<'py, bool>,
    stores: PyReadonlyArray1<'py, bool>,
    tau: PyReadonlyArray1<'py, f64>,
    dt: f64,
    storage: PyReadonlyArray1<'py, f64>,
    mut runoff: numpy::PyReadwriteArray1<'py, f64>,
    mut melt: numpy::PyReadwriteArray1<'py, f64>,
    mut prior: numpy::PyReadwriteArray1<'py, f64>,
    runoff_loss: PyReadonlyArray1<'py, f64>,
    melt_loss: PyReadonlyArray1<'py, f64>,
) -> PyResult<Vec<Bound<'py, PyArray1<f64>>>> {
    let n = target.len();
    if ocean.len() != n
        || stores.len() != n
        || tau.len() != n
        || storage.len() != n
        || runoff.len() != n
        || melt.len() != n
        || prior.len() != n
        || runoff_loss.len() != n
        || melt_loss.len() != n
    {
        return Err(PyValueError::new_err("incompatible river reach shapes"));
    }
    let order = order.as_slice()?;
    let target = target.as_slice()?;
    if order.iter().any(|&i| i < 0 || i as usize >= n) || target.iter().any(|&i| i >= n as i64) {
        return Err(PyValueError::new_err("river reach index out of range"));
    }
    let ocean = ocean.as_slice()?;
    let stores = stores.as_slice()?;
    let tau = tau.as_slice()?;
    let storage = storage.as_slice()?;
    let runoff_loss = runoff_loss.as_slice()?;
    let melt_loss = melt_loss.as_slice()?;
    let runoff = runoff.as_slice_mut()?;
    let melt = melt.as_slice_mut()?;
    let prior = prior.as_slice_mut()?;
    let output = py.detach(|| {
        let mut out = vec![vec![0.0; n]; 10];
        for &raw in order {
            let i = raw as usize;
            let total = runoff[i] + melt[i] + prior[i];
            out[8][i] = total;
            let mut pass = 1.0;
            let mut released = storage[i];
            if stores[i] && dt > 0.0 {
                let ratio = dt / tau[i];
                let drained = -(-ratio).exp_m1();
                released = storage[i] * drained;
                let held = total * drained / ratio;
                out[0][i] = storage[i] - released + held;
                pass = 1.0 - drained / ratio;
            } else if stores[i] {
                released = 0.0;
                pass = 0.0;
                out[0][i] = storage[i];
            }
            let out_runoff = runoff[i] * pass;
            let out_melt = melt[i] * pass;
            let out_prior = prior[i] * pass + released;
            out[1][i] = out_runoff + out_melt + out_prior;
            out[9][i] = storage[i] + total - out[0][i] - out[1][i];
            let kept_runoff = out_runoff * (1.0 - runoff_loss[i]);
            let kept_melt = out_melt * (1.0 - melt_loss[i]);
            let kept_prior = out_prior * (1.0 - runoff_loss[i]);
            out[2][i] = kept_runoff + kept_melt + kept_prior;
            out[7][i] = out[1][i] - out[2][i];
            let j = target[i];
            if j < 0 || ocean[j as usize] {
                let j = if j < 0 { i } else { j as usize };
                out[3][j] += kept_runoff;
                out[4][j] += kept_melt;
                out[5][j] += kept_prior;
                if target[i] >= 0 {
                    out[6][j] += out[2][i];
                }
            } else {
                let j = j as usize;
                runoff[j] += kept_runoff;
                melt[j] += kept_melt;
                prior[j] += kept_prior;
            }
        }
        out
    });
    Ok(output
        .into_iter()
        .map(|a| PyArray1::from_vec(py, a))
        .collect())
}
