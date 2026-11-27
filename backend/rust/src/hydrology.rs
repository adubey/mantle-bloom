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
