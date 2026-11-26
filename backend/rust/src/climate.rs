use numpy::{PyArray1, PyReadonlyArray1, PyReadwriteArray1};
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
