use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadwriteArray1};
use pyo3::prelude::*;

pub const GRAVITY_M_S2: f64 = 9.81;
pub const RHO_WATER: f64 = 1000.0;
pub const RHO_ICE: f64 = 917.0;
pub const RHO_LITHOSPHERE_MANTLE: f64 = 3300.0;
pub const RHO_ASTHENOSPHERE: f64 = 3250.0;
pub const ISOSTATIC_REFERENCE_OFFSET_M: f64 = -4184.615384615388;
pub const MIN_ELEVATION_M: f64 = -11000.0;
pub const MAX_ELEVATION_M: f64 = 9000.0;
pub const LAPSE_RATE_C_PER_KM: f64 = 6.5;

pub const ICE_YIELD_STRESS_HARD_PA: f64 = 100_000.0;
pub const ICE_YIELD_STRESS_SOFT_PA: f64 = 30_000.0;
pub const ICE_SOFT_BED_COVER_M: f64 = 20.0;

pub const ICE_GEOTHERMAL_FLUX_W_M2: f64 = 0.06;
pub const ICE_THERMAL_CONDUCTIVITY_W_M_K: f64 = 2.1;
pub const ICE_PRESSURE_MELT_C_PER_M: f64 = 8.7e-4;

pub const ICE_LATENT_HEAT_J_KG: f64 = 334_000.0;
pub const SECONDS_PER_YEAR: f64 = 3.15576e7;
pub const W_M2_TO_M_ICE_PER_YR: f64 = SECONDS_PER_YEAR / (RHO_ICE * ICE_LATENT_HEAT_J_KG);

pub const DDF_SNOW_MM_WATER_PER_DAY_K: f64 = 4.0;
pub const DDF_ICE_MM_WATER_PER_DAY_K: f64 = 8.0;
pub const RAIN_SNOW_THRESHOLD_C: f64 = 1.0;
pub const ICE_SNOW_RAMP_C: f64 = 1.0;

pub const GLACIER_MAX_DEPTH_M: f64 = 5000.0;
pub const GLACIER_VISIBLE_DEPTH_M: f64 = 10.0;
pub const FREEZE_POINT_C: f64 = 0.0;

pub const GLACIAL_ABRASION_COEFFICIENT: f64 = 1.0e-4;
pub const GLACIAL_ABRASION_EXPONENT: f64 = 0.65;

pub const ICE_BODY_MIN_NODES: usize = 3;
pub const ICE_SLIDING_FRACTION: f64 = 0.7;
pub const ICE_VELOCITY_MIN_THICKNESS_M: f64 = 10.0;

pub const ICE_ROOT_TOLERANCE_M: f64 = 100.0;
pub const ICE_ROOT_MAX_ITERATIONS: usize = 40;
pub const ICE_FLOTATION_PASSES: usize = 6;
pub const FLOTATION_RATIO: f64 = RHO_WATER / RHO_ICE;

#[inline]
pub fn clip_min(x: f64, lo: f64) -> f64 {
    if x.is_nan() {
        f64::NAN
    } else if x < lo {
        lo
    } else {
        x
    }
}

#[inline]
pub fn clip_range(x: f64, lo: f64, hi: f64) -> f64 {
    if x.is_nan() {
        f64::NAN
    } else if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

#[inline]
pub fn min_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a < b {
        a
    } else {
        b
    }
}

#[inline]
pub fn max_nan(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else if a > b {
        a
    } else {
        b
    }
}

#[derive(Clone, Copy, Default)]
pub struct NeumaierSum {
    sum: f64,
    c: f64,
}

impl NeumaierSum {
    #[inline]
    pub fn add(&mut self, x: f64) {
        let t = self.sum + x;
        if self.sum.abs() >= x.abs() {
            self.c += (self.sum - t) + x;
        } else {
            self.c += (x - t) + self.sum;
        }
        self.sum = t;
    }

    #[inline]
    pub fn finish(self) -> f64 {
        self.sum + self.c
    }
}

/// Neumaier compensated dot product for high-precision summation matching NumPy/BLAS `np.dot`.
pub fn compensated_dot(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut acc = NeumaierSum::default();
    for i in 0..a.len() {
        acc.add(a[i] * b[i]);
    }
    acc.finish()
}

/// Neumaier compensated dot product over an index slice.
pub fn compensated_dot_indexed(indices: &[usize], a: &[f64], b: &[f64]) -> f64 {
    let mut acc = NeumaierSum::default();
    for &u in indices {
        acc.add(a[u] * b[u]);
    }
    acc.finish()
}

#[inline]
pub fn air_temperature_at_elevation_scalar(
    air_temperature_c: f64,
    base_elevation_m: f64,
    target_elevation_m: f64,
) -> f64 {
    air_temperature_c - LAPSE_RATE_C_PER_KM * (target_elevation_m - base_elevation_m) / 1000.0
}

#[inline]
pub fn snow_fraction_scalar(temperature_c: f64, amplitude_c: f64, threshold_c: f64) -> f64 {
    let eff = max_nan(amplitude_c, ICE_SNOW_RAMP_C / 2.0);
    let ratio = clip_range((threshold_c - temperature_c) / max_nan(eff, 1e-6), -1.0, 1.0);
    let seasonal = 0.5 + ratio.asin() / std::f64::consts::PI;
    let ramp = clip_range(0.5 + (threshold_c - temperature_c) / ICE_SNOW_RAMP_C, 0.0, 1.0);
    if amplitude_c > 1e-3 {
        seasonal
    } else {
        ramp
    }
}

#[inline]
pub fn seasonal_positive_degrees_scalar(mean_c: f64, amplitude_c: f64, threshold_c: f64) -> f64 {
    let d = mean_c - threshold_c;
    let amp = clip_min(amplitude_c, 0.0);
    if d.abs() < amp {
        (d * (-d / amp).acos() + (amp * amp - d * d).sqrt()) / std::f64::consts::PI
    } else {
        clip_min(d, 0.0)
    }
}

#[inline]
pub fn basal_heat_scalar(thickness_m: f64, surface_c: f64, frictional_w_m2: f64) -> f64 {
    if thickness_m > 0.0 {
        let t_pmp = -ICE_PRESSURE_MELT_C_PER_M * thickness_m;
        let conducted =
            ICE_THERMAL_CONDUCTIVITY_W_M_K * clip_min(t_pmp - surface_c, 0.0) / thickness_m.max(1e-9);
        ICE_GEOTHERMAL_FLUX_W_M2 + frictional_w_m2 - conducted
    } else {
        f64::NEG_INFINITY
    }
}

#[inline]
pub fn bare_rates_scalar(
    temperature_c: f64,
    seasonal_amplitude_c: f64,
    precipitation_mm: f64,
    bed_frozen_snow_mm: Option<f64>,
) -> (f64, f64) {
    let to_ice = RHO_WATER / RHO_ICE;
    let mut snow_mm =
        snow_fraction_scalar(temperature_c, seasonal_amplitude_c, RAIN_SNOW_THRESHOLD_C)
            * precipitation_mm;
    if let Some(bfs) = bed_frozen_snow_mm {
        snow_mm += bfs;
    }
    let snow_we = snow_mm * 1.0e-3;
    let pdd = 365.25 * seasonal_positive_degrees_scalar(temperature_c, seasonal_amplitude_c, 0.0);
    let ddf_snow_m = DDF_SNOW_MM_WATER_PER_DAY_K * 1.0e-3;
    let snow_capacity = ddf_snow_m * pdd;
    let snow_melt = min_nan(snow_we, snow_capacity);
    (snow_we * to_ice * 1.0e6, snow_melt * to_ice * 1.0e6)
}

#[inline]
pub fn node_rates_scalar(
    thickness_m: f64,
    temperature_c: f64,
    seasonal_amplitude_c: f64,
    precipitation_mm: f64,
    frictional_heat_w_m2: f64,
    bed_frozen_snow_mm: Option<f64>,
) -> (f64, f64, f64) {
    let thickness = clip_min(thickness_m, 0.0);
    let surface_c = air_temperature_at_elevation_scalar(temperature_c, 0.0, thickness);
    let to_ice = RHO_WATER / RHO_ICE;
    let mut snow_mm =
        snow_fraction_scalar(surface_c, seasonal_amplitude_c, RAIN_SNOW_THRESHOLD_C) * precipitation_mm;
    if let Some(bfs) = bed_frozen_snow_mm {
        snow_mm += bfs;
    }
    let snow_we = snow_mm * 1.0e-3;
    let pdd = 365.25 * seasonal_positive_degrees_scalar(surface_c, seasonal_amplitude_c, 0.0);
    let ddf_snow_m = DDF_SNOW_MM_WATER_PER_DAY_K * 1.0e-3;
    let snow_capacity = ddf_snow_m * pdd;
    let snow_melt = min_nan(snow_we, snow_capacity);
    let pdd_left = if DDF_SNOW_MM_WATER_PER_DAY_K > 0.0 {
        clip_min(pdd - snow_melt / ddf_snow_m, 0.0)
    } else {
        pdd
    };
    let ice_melt = if thickness > 0.0 {
        DDF_ICE_MM_WATER_PER_DAY_K * 1.0e-3 * pdd_left
    } else {
        0.0
    };
    let basal = W_M2_TO_M_ICE_PER_YR
        * clip_min(basal_heat_scalar(thickness, surface_c, frictional_heat_w_m2), 0.0);
    (
        snow_we * to_ice * 1.0e6,
        (snow_melt + ice_melt) * to_ice * 1.0e6,
        basal * 1.0e6,
    )
}

#[inline]
pub fn warm_based_scalar(thickness_m: f64, surface_temp_c: f64, frictional_w_m2: f64) -> bool {
    basal_heat_scalar(thickness_m, surface_temp_c, frictional_w_m2) >= 0.0
}

#[inline]
pub fn yield_stress_scalar(mobile_cover_m: f64, is_warm: bool) -> f64 {
    let softness = if is_warm {
        clip_range(mobile_cover_m / ICE_SOFT_BED_COVER_M, 0.0, 1.0)
    } else {
        0.0
    };
    ICE_YIELD_STRESS_HARD_PA + (ICE_YIELD_STRESS_SOFT_PA - ICE_YIELD_STRESS_HARD_PA) * softness
}

#[inline]
pub fn floating_water_depth_scalar(bed_elevation_m: f64, sea_level_m: f64, is_ocean: bool) -> f64 {
    if is_ocean {
        clip_min(sea_level_m - bed_elevation_m, 0.0)
    } else {
        0.0
    }
}

#[inline]
pub fn is_grounded_scalar(depth_m: f64, unloaded_bed_m: f64, sea_level_m: f64, is_ocean: bool) -> bool {
    let water_depth = floating_water_depth_scalar(unloaded_bed_m, sea_level_m, is_ocean);
    RHO_ICE * clip_min(depth_m, 0.0) - RHO_WATER * water_depth > 0.0
}

#[inline]
pub fn isostatic_elevation_scalar(hc_m: f64, hm_m: f64, rho_c: f64, surface_load_kg_m2: f64) -> f64 {
    let bracket = hc_m * (1.0 - rho_c / RHO_ASTHENOSPHERE)
        + hm_m * (1.0 - RHO_LITHOSPHERE_MANTLE / RHO_ASTHENOSPHERE)
        - surface_load_kg_m2 / RHO_ASTHENOSPHERE;
    let shifted_bracket = bracket + ISOSTATIC_REFERENCE_OFFSET_M;
    let z = if shifted_bracket <= 0.0 {
        shifted_bracket * (RHO_ASTHENOSPHERE / (RHO_ASTHENOSPHERE - RHO_WATER))
    } else {
        shifted_bracket
    };
    clip_range(z, MIN_ELEVATION_M, MAX_ELEVATION_M)
}

#[inline]
pub fn ice_load_deflection_scalar(hc_m: f64, hm_m: f64, rho_c: f64, load_kg_m2: f64) -> f64 {
    if hc_m > 0.0 {
        isostatic_elevation_scalar(hc_m, hm_m, rho_c, load_kg_m2)
            - isostatic_elevation_scalar(hc_m, hm_m, rho_c, 0.0)
    } else {
        -load_kg_m2 / RHO_ASTHENOSPHERE
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn march_substep_prep<'py>(
    py: Python<'py>,
    working_elevation: PyReadonlyArray1<'py, f64>,
    deflection_m: PyReadonlyArray1<'py, f64>,
    ref_depth: PyReadonlyArray1<'py, f64>,
    sea_level: f64,
    ocean_mask: PyReadonlyArray1<'py, bool>,
    crust_h_m: PyReadonlyArray1<'py, f64>,
    mantle_h_m: PyReadonlyArray1<'py, f64>,
    rho_crust: PyReadonlyArray1<'py, f64>,
    temperature_at_nodes: PyReadonlyArray1<'py, f64>,
) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<bool>>)> {
    let elev = working_elevation.as_slice()?;
    let defl = deflection_m.as_slice()?;
    let depth = ref_depth.as_slice()?;
    let ocean = ocean_mask.as_slice()?;
    let hc = crust_h_m.as_slice()?;
    let hm = mantle_h_m.as_slice()?;
    let rho = rho_crust.as_slice()?;
    let temp0 = temperature_at_nodes.as_slice()?;
    let n = elev.len();

    let (temps, frozen) = py.detach(|| {
        let mut temps = Vec::with_capacity(n);
        let mut frozen = Vec::with_capacity(n);
        for i in 0..n {
            let unloaded_bed = elev[i] - defl[i];
            let wd = floating_water_depth_scalar(unloaded_bed, sea_level, ocean[i]);
            let load = clip_min(RHO_ICE * clip_min(depth[i], 0.0) - RHO_WATER * wd, 0.0);
            let target = ice_load_deflection_scalar(hc[i], hm[i], rho[i], load);
            let t = air_temperature_at_elevation_scalar(temp0[i], defl[i], target);
            temps.push(t);
            frozen.push((t < FREEZE_POINT_C) || (depth[i] >= GLACIER_VISIBLE_DEPTH_M));
        }
        (temps, frozen)
    });

    Ok((PyArray1::from_vec(py, temps), PyArray1::from_vec(py, frozen)))
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn apply_march_abrasion<'py>(
    py: Python<'py>,
    sliding_speed_m_per_yr: PyReadonlyArray1<'py, f64>,
    dt: f64,
    multiplier: f64,
    ocean_mask: PyReadonlyArray1<'py, bool>,
    mut cover_remaining: PyReadwriteArray1<'py, f64>,
    mut substrate_headroom: PyReadwriteArray1<'py, f64>,
    craton_keep: PyReadonlyArray1<'py, f64>,
    mut raw_abrasion_depth: PyReadwriteArray1<'py, f64>,
    mut abrasion_depth: PyReadwriteArray1<'py, f64>,
    mut working_elevation: PyReadwriteArray1<'py, f64>,
) -> PyResult<f64> {
    let sliding = sliding_speed_m_per_yr.as_slice()?;
    let ocean = ocean_mask.as_slice()?;
    let craton = craton_keep.as_slice()?;
    let cover = cover_remaining.as_slice_mut()?;
    let headroom = substrate_headroom.as_slice_mut()?;
    let raw_out = raw_abrasion_depth.as_slice_mut()?;
    let abr_out = abrasion_depth.as_slice_mut()?;
    let elev = working_elevation.as_slice_mut()?;
    let n = sliding.len();

    let mult = multiplier.max(0.0);
    let max_cut = py.detach(|| {
        let mut max_cut = 0.0_f64;
        for i in 0..n {
            let raw = if ocean[i] {
                0.0
            } else {
                let s = clip_min(sliding[i], 0.0);
                GLACIAL_ABRASION_COEFFICIENT * mult * s.powf(GLACIAL_ABRASION_EXPONENT) * dt
            };
            raw_out[i] += raw;
            let cover_cut = min_nan(raw, cover[i]);
            let substrate_demand = clip_min(raw - cover_cut, 0.0);
            let substrate_take = min_nan(substrate_demand, headroom[i]);
            let substrate_cut = substrate_take * craton[i];
            let interval_abr = cover_cut + substrate_cut;
            cover[i] -= cover_cut;
            headroom[i] -= substrate_cut;
            abr_out[i] += interval_abr;
            elev[i] -= interval_abr;
            if interval_abr > max_cut {
                max_cut = interval_abr;
            }
        }
        max_cut
    });

    Ok(max_cut)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rates_known_values() {
        let (acc0, sm0, bm0) = node_rates_scalar(0.0, -10.0, 5.0, 500.0, 0.0, None);
        assert!((acc0 - 500.0e-3 * (RHO_WATER / RHO_ICE) * 1.0e6).abs() < 1e-6);
        assert_eq!(sm0, 0.0);
        assert_eq!(bm0, 0.0);

        let (acc1, sm1, bm1) = node_rates_scalar(1500.0, -10.0, 5.0, 500.0, 0.05, Some(10.0));
        assert!(acc1 > acc0);
        assert_eq!(sm1, 0.0);
        assert!(bm1 > 0.0);
    }
}

/// Two periodic laps of the diagnostic zonal moisture recurrence (#1046).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn humidity_zonal_sweep_rs<'py>(
    py: Python<'py>,
    ocean: PyReadonlyArray2<'py, bool>,
    elevation: PyReadonlyArray2<'py, f64>,
    ceiling: PyReadonlyArray2<'py, f64>,
    factor: PyReadonlyArray2<'py, f64>,
    source: PyReadonlyArray2<'py, f64>,
    direction: PyReadonlyArray1<'py, f64>,
    retention: PyReadonlyArray1<'py, f64>,
    min_wind: f64,
    max_humidity: f64,
    lift_scale: f64,
    shadow_factor: f64,
) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>)> {
    let ocean = ocean.as_array();
    let elevation = elevation.as_array();
    let ceiling = ceiling.as_array();
    let factor = factor.as_array();
    let source = source.as_array();
    let direction = direction.as_array();
    let retention = retention.as_array();
    let (height, width) = ocean.dim();
    if height == 0
        || width == 0
        || elevation.dim() != (height, width)
        || ceiling.dim() != (height, width)
        || factor.dim() != (height, width)
        || source.dim() != (height, width)
        || direction.len() != height
        || retention.len() != height
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "climate sweep shape mismatch",
        ));
    }
    let (humidity, orographic) = py.detach(|| {
        let mut humidity = numpy::ndarray::Array2::<f64>::zeros((height, width));
        let mut orographic = numpy::ndarray::Array2::<f64>::zeros((height, width));
        for row in 0..height {
            let mut moisture = 0.0;
            let mut previous = 0.0;
            for lap in 0..2 {
                for i in 0..width {
                    let col = if direction[row] == 1.0 {
                        i
                    } else {
                        width - 1 - i
                    };
                    let elev = elevation[[row, col]];
                    let wind_factor = min_wind + (1.0 - min_wind) * factor[[row, col]];
                    let land = min_nan(
                        moisture * (retention[row] * wind_factor) + source[[row, col]],
                        max_humidity,
                    );
                    let after_source = if ocean[[row, col]] {
                        ceiling[[row, col]]
                    } else {
                        land
                    };
                    let gain = elev - previous;
                    let retained = if gain > 0.0 {
                        let saturation = 1.0 - (-gain / lift_scale).exp();
                        1.0 - saturation * (1.0 - shadow_factor)
                    } else {
                        1.0
                    };
                    let next = if ocean[[row, col]] {
                        after_source
                    } else {
                        after_source * retained
                    };
                    if lap == 1 {
                        humidity[[row, col]] = next;
                        orographic[[row, col]] = if ocean[[row, col]] {
                            0.0
                        } else {
                            after_source - next
                        };
                    }
                    moisture = next;
                    previous = elev;
                }
            }
        }
        (humidity, orographic)
    });
    Ok((
        PyArray2::from_owned_array(py, humidity),
        PyArray2::from_owned_array(py, orographic),
    ))
}

/// Ocean-only Jacobi smoothing and coastal redirection (#1046).
#[pyfunction]
#[allow(clippy::too_many_arguments)]
pub fn smooth_along_coast_rs<'py>(
    py: Python<'py>,
    u: PyReadonlyArray2<'py, f64>,
    v: PyReadonlyArray2<'py, f64>,
    ocean: PyReadonlyArray2<'py, bool>,
    lat: PyReadonlyArray1<'py, f64>,
    ambient_u: PyReadonlyArray2<'py, f64>,
    ambient_v: PyReadonlyArray2<'py, f64>,
    iterations: usize,
    weight: f64,
    restore: f64,
    speedup: f64,
    eps: f64,
) -> PyResult<(Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>)> {
    let ocean = ocean.as_array();
    let lat = lat.as_array();
    let ambient_u = ambient_u.as_array();
    let ambient_v = ambient_v.as_array();
    let u_in = u.as_array();
    let v_in = v.as_array();
    let (h, w) = ocean.dim();
    if u_in.dim() != (h, w)
        || v_in.dim() != (h, w)
        || ambient_u.dim() != (h, w)
        || ambient_v.dim() != (h, w)
        || lat.len() != h
        || h == 0
        || w == 0
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coast smoothing shape mismatch",
        ));
    }
    let (u, v) = py.detach(|| {
        let mut u = u_in.to_owned();
        let mut v = v_in.to_owned();
        let mut next_u = numpy::ndarray::Array2::<f64>::zeros((h, w));
        let mut next_v = numpy::ndarray::Array2::<f64>::zeros((h, w));
        for _ in 0..iterations {
            for r in 0..h {
                let north = (r + h - 1) % h;
                let south = (r + 1) % h;
                for c in 0..w {
                    let east = (c + 1) % w;
                    let west = (c + w - 1) % w;
                    let neighbors = [(r, east), (r, west), (south, c), (north, c)];
                    let mut count = 0.0;
                    let mut total_u = 0.0;
                    let mut total_v = 0.0;
                    for &(nr, nc) in &neighbors {
                        let mask = if ocean[[nr, nc]] { 1.0 } else { 0.0 };
                        count += mask;
                        total_u += u[[nr, nc]] * mask;
                        total_v += v[[nr, nc]] * mask;
                    }
                    let avg_u = if count > 0.0 { total_u / count } else { 0.0 };
                    let avg_v = if count > 0.0 { total_v / count } else { 0.0 };
                    let mut uu = u[[r, c]] + weight * (avg_u - u[[r, c]]);
                    let mut vv = v[[r, c]] + weight * (avg_v - v[[r, c]]);
                    uu += restore * (ambient_u[[r, c]] - uu);
                    vv += restore * (ambient_v[[r, c]] - vv);
                    let nx = (!ocean[[r, west]]) as i32 - (!ocean[[r, east]]) as i32;
                    let ny = (!ocean[[south, c]]) as i32 - (!ocean[[north, c]]) as i32;
                    let mag = (nx as f64).hypot(ny as f64);
                    if ocean[[r, c]] && mag > eps {
                        let nx = nx as f64 / mag;
                        let ny = ny as f64 / mag;
                        let hemisphere = if lat[r] >= 0.0 { 1.0 } else { -1.0 };
                        let into = max_nan(-(uu * nx + vv * ny), 0.0);
                        let redirected = into * speedup;
                        next_u[[r, c]] = uu + into * nx + redirected * (-hemisphere * ny);
                        next_v[[r, c]] = vv + into * ny + redirected * (hemisphere * nx);
                    } else {
                        next_u[[r, c]] = uu;
                        next_v[[r, c]] = vv;
                    }
                }
            }
            std::mem::swap(&mut u, &mut next_u);
            std::mem::swap(&mut v, &mut next_v);
        }
        (u, v)
    });
    Ok((
        PyArray2::from_owned_array(py, u),
        PyArray2::from_owned_array(py, v),
    ))
}

/// Compute NumPy-compatible nearest offset indices (#1046).
#[pyfunction]
pub fn offset_indices_rs<'py>(
    py: Python<'py>,
    u: PyReadonlyArray2<'py, f64>,
    v: PyReadonlyArray2<'py, f64>,
    lat: PyReadonlyArray1<'py, f64>,
    cos_lat: PyReadonlyArray1<'py, f64>,
    lon: PyReadonlyArray1<'py, f64>,
    distance: f64,
) -> PyResult<Option<(Bound<'py, PyArray2<i64>>, Bound<'py, PyArray2<i64>>)>> {
    let u = u.as_array();
    let v = v.as_array();
    let lat = lat.as_array();
    let cos_lat = cos_lat.as_array();
    let lon = lon.as_array();
    let (h, w) = u.dim();
    if h == 0
        || w == 0
        || v.dim() != (h, w)
        || lat.len() != h
        || cos_lat.len() != h
        || lon.len() != w
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "offset shape mismatch",
        ));
    }
    let indices = py.detach(|| {
        let mut rows = numpy::ndarray::Array2::<i64>::zeros((h, w));
        let mut cols = numpy::ndarray::Array2::<i64>::zeros((h, w));
        for r in 0..h {
            for c in 0..w {
                let du = u[[r, c]];
                let dv = v[[r, c]];
                if !du.is_finite() || !dv.is_finite() || du.abs() > 1000.0 || dv.abs() > 1000.0 {
                    return None;
                }
                let dst_lat = lat[r] + dv * distance;
                let dst_lon = lon[c] + du * distance / cos_lat[r];
                let row = ((90.0 - dst_lat) / (180.0 / h as f64) - 0.5).round_ties_even() as i64;
                let col = ((dst_lon + 180.0) / (360.0 / w as f64) - 0.5).round_ties_even() as i64;
                rows[[r, c]] = row.clamp(0, h as i64 - 1);
                cols[[r, c]] = col.rem_euclid(w as i64);
            }
        }
        Some((rows, cols))
    });
    Ok(indices.map(|(rows, cols)| {
        (
            PyArray2::from_owned_array(py, rows),
            PyArray2::from_owned_array(py, cols),
        )
    }))
}

/// Spherical centered moisture divergence before smoothing (#1046).
#[pyfunction]
pub fn moisture_divergence_rs<'py>(
    py: Python<'py>,
    humidity: PyReadonlyArray2<'py, f64>,
    u: PyReadonlyArray2<'py, f64>,
    v: PyReadonlyArray2<'py, f64>,
    cos_lat: PyReadonlyArray1<'py, f64>,
    radius: f64,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let q = humidity.as_array();
    let u = u.as_array();
    let v = v.as_array();
    let cos = cos_lat.as_array();
    let (h, w) = q.dim();
    if h == 0 || w == 0 || u.dim() != (h, w) || v.dim() != (h, w) || cos.len() != h {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "moisture divergence shape mismatch",
        ));
    }
    let dlon = 2.0 * std::f64::consts::PI / w as f64;
    let dlat = std::f64::consts::PI / h as f64;
    let result = py.detach(|| {
        let mut result = numpy::ndarray::Array2::<f64>::zeros((h, w));
        for r in 0..h {
            let n = (r + h - 1) % h;
            let s = (r + 1) % h;
            for c in 0..w {
                let e = (c + 1) % w;
                let west = (c + w - 1) % w;
                let dqu = (q[[r, e]] * u[[r, e]] - q[[r, west]] * u[[r, west]]) / 2.0;
                let dqv =
                    ((q[[n, c]] * v[[n, c]]) * cos[n] - (q[[s, c]] * v[[s, c]]) * cos[s]) / 2.0;
                result[[r, c]] = (dqu / dlon + dqv / dlat) / (cos[r] * radius);
            }
        }
        result
    });
    Ok(PyArray2::from_owned_array(py, result))
}

/// Periodic Jacobi smoothing, preserving NumPy's neighbor addition order (#1046).
#[pyfunction]
pub fn smooth_field_rs<'py>(
    py: Python<'py>,
    field: PyReadonlyArray2<'py, f64>,
    iterations: usize,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let (h, w) = field.as_array().dim();
    if h == 0 || w == 0 {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "smoothing shape mismatch",
        ));
    }
    let mut values = field.as_slice()?.to_vec();
    let values = py.detach(|| {
        let mut scratch = vec![0.0; h * w];
        for _ in 0..iterations {
            for r in 0..h {
                let north = if r == 0 { h - 1 } else { r - 1 };
                let south = if r + 1 == h { 0 } else { r + 1 };
                for c in 0..w {
                    let east = if c + 1 == w { 0 } else { c + 1 };
                    let west = if c == 0 { w - 1 } else { c - 1 };
                    scratch[r * w + c] = 0.25
                        * (values[r * w + east]
                            + values[r * w + west]
                            + values[south * w + c]
                            + values[north * w + c]);
                }
            }
            std::mem::swap(&mut values, &mut scratch);
        }
        values
    });
    let out = numpy::ndarray::Array2::from_shape_vec((h, w), values)
        .expect("smoothing output has the validated grid shape");
    Ok(PyArray2::from_owned_array(py, out))
}
