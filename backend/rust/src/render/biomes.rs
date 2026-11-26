use std::cmp::Ordering;

use rayon::prelude::*;

pub const PLANET_RADIUS_KM: f64 = 6371.0;
pub const CONTINENTALITY_SATURATION_RAD: f64 = 12.0 * std::f64::consts::PI / 180.0;
pub const SHELF_RANGE_RAD: f64 = 2.0 * std::f64::consts::PI / 180.0;

pub const N_BIOMES: usize = 41;
pub const ICE_CAP: u8 = 30;
pub const POLAR_SEA_ICE: u8 = 40;
pub const COMBINED_LAKE_ID_CODE: u8 = 42;
pub const COMBINED_SEA_ID_CODE: u8 = 44;

pub const BIOME_SHADE_AMPLITUDE: f64 = 0.25;
pub const BIOME_SHADE_MIN: f64 = 1.0 - BIOME_SHADE_AMPLITUDE;
pub const BIOME_SHADE_MAX: f64 = 1.0 + BIOME_SHADE_AMPLITUDE;

pub const BIOME_VOTE_MIN_NEIGHBOUR_FRACTION: f64 = 0.65;
pub const BIOME_VOTE_ITERATIONS: usize = 2;
pub const BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS: usize = 6;

const TROPICAL_TEMP_C: f64 = 18.0;
const SEASON_LAT_AMPLITUDE_C: f64 = 15.5;
const SEASON_LAT_EXPONENT: f64 = 1.3;
const SEASON_CONTINENTAL_GAIN: f64 = 0.9;
const SEASON_CONTINENTAL_OFFSET_C: f64 = 8.0;
const SEASON_TILT_REFERENCE_DEG: f64 = 23.5;
const SEASON_MAX_AMPLITUDE_C: f64 = 33.0;
const SEASON_RELIEF_GAIN_C_PER_KM: f64 = 2.5;
const SEASON_RELIEF_MAX_C: f64 = 8.0;

const SOFT_TEMP_BLEND_C: f64 = 1.4;
const SOFT_SHARE_BLEND: f64 = 0.04;
const SOFT_MONTHS_BLEND: f64 = 0.5;
const SOFT_LAT_BLEND_DEG: f64 = 1.5;
const SOFT_CONTINENTALITY_BLEND: f64 = 0.05;
const SOFT_PRECIP_MONTH_BLEND_MM: f64 = 12.0;
const SOFT_PRECIP_FRACTION: f64 = 0.08;
const SOFT_PRECIP_FLOOR_MM: f64 = 20.0;
const SOFT_SHELF_DIST_BLEND_RAD: f64 = 0.4 * std::f64::consts::PI / 180.0;
pub const SOFT_SEA_ICE_FRACTION_BLEND: f64 = 0.05;
pub const MIN_DOMINANT_SHARE: f64 = 0.5 + 1e-4;

const PELAGIC_SST_EDGES: [f64; 4] = [2.0, 12.0, 19.0, 24.0];
const EQUATORIAL_LAT_DEG: f64 = 6.0;

pub const BIOME_COLORS: [[u8; 3]; N_BIOMES] = [
    // Koppen (31)
    [26, 51, 21],    // 0: Af
    [34, 66, 29],    // 1: Am
    [74, 92, 44],    // 2: Aw
    [92, 103, 52],   // 3: As
    [176, 141, 99],  // 4: BWh
    [150, 136, 104], // 5: BWk
    [156, 134, 78],  // 6: BSh
    [124, 128, 82],  // 7: BSk
    [150, 138, 66],  // 8: Csa
    [128, 136, 76],  // 9: Csb
    [112, 126, 86],  // 10: Csc
    [94, 116, 54],   // 11: Cwa
    [100, 122, 66],  // 12: Cwb
    [106, 124, 80],  // 13: Cwc
    [58, 96, 44],    // 14: Cfa
    [72, 104, 56],   // 15: Cfb
    [92, 112, 76],   // 16: Cfc
    [132, 126, 80],  // 17: Dsa
    [124, 124, 84],  // 18: Dsb
    [120, 122, 92],  // 19: Dsc
    [128, 126, 104], // 20: Dsd
    [96, 110, 56],   // 21: Dwa
    [104, 116, 68],  // 22: Dwb
    [116, 118, 84],  // 23: Dwc
    [128, 122, 100], // 24: Dwd
    [84, 104, 52],   // 25: Dfa
    [100, 114, 66],  // 26: Dfb
    [122, 112, 80],  // 27: Dfc
    [140, 128, 100], // 28: Dfd
    [150, 152, 132], // 29: ET
    [234, 238, 240], // 30: EF
    // Pelagic (10)
    [7, 17, 43],     // 31: Tropical Open Ocean
    [4, 12, 36],     // 32: Subtropical Gyre
    [14, 40, 78],    // 33: Equatorial Divergence
    [34, 82, 122],   // 34: Tropical Coastal Waters
    [10, 26, 58],    // 35: Temperate Open Ocean
    [30, 70, 106],   // 36: Temperate Shelf
    [18, 42, 72],    // 37: Cold-Temperate Open Ocean
    [38, 74, 104],   // 38: Cold-Temperate Shelf
    [32, 62, 94],    // 39: Polar Ocean
    [226, 233, 237], // 40: Polar Sea Ice
];

/// Pure-Rust `erf(x)` accurate to < 1e-15 (Sun Microsystems / fdlibm `s_erf.c`).
#[inline]
pub fn erf_f64(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let ax = x.abs();
    let sign = if x < 0.0 { -1.0 } else { 1.0 };

    if ax < 0.84375 {
        if ax < 1.4901161193847656e-8 {
            return x * (2.0 / std::f64::consts::PI.sqrt());
        }
        let z = x * x;
        const PP: [f64; 5] = [
            1.28379167095512558561e-01,
            -3.25042107247001499370e-01,
            -2.84817495755985104766e-02,
            -5.77027029648944159157e-03,
            -2.37630166566501626084e-05,
        ];
        const QQ: [f64; 6] = [
            1.0,
            3.97917223959155352819e-01,
            6.50222499887672944485e-02,
            5.08130628187576562776e-03,
            1.32494738004321644526e-04,
            -3.96022827877536812320e-06,
        ];
        let r = PP[0] + z * (PP[1] + z * (PP[2] + z * (PP[3] + z * PP[4])));
        let s = QQ[0] + z * (QQ[1] + z * (QQ[2] + z * (QQ[3] + z * (QQ[4] + z * QQ[5]))));
        return x + x * (r / s);
    }

    if ax < 1.25 {
        let s2 = ax - 1.0;
        const PA: [f64; 7] = [
            -2.36211856075265944077e-03,
            4.14856118683748331666e-01,
            -3.72207876035701323847e-01,
            3.18346619901161753674e-01,
            -1.10894694282396677476e-01,
            3.54783043256182359371e-02,
            -2.16637559486879084300e-03,
        ];
        const QA: [f64; 7] = [
            1.0,
            1.06420880400844228286e-01,
            5.40397917702171048935e-01,
            7.18286544141962662868e-02,
            1.26171219808761642112e-01,
            1.36370839120290507362e-02,
            1.19844998467991074170e-02,
        ];
        let p = PA[0] + s2 * (PA[1] + s2 * (PA[2] + s2 * (PA[3] + s2 * (PA[4] + s2 * (PA[5] + s2 * PA[6])))));
        let q = QA[0] + s2 * (QA[1] + s2 * (QA[2] + s2 * (QA[3] + s2 * (QA[4] + s2 * (QA[5] + s2 * QA[6])))));
        return sign * (0.845062911510467529297 + p / q);
    }

    if ax >= 6.0 {
        return sign;
    }

    let z = 1.0 / (ax * ax);
    let (r, s) = if ax < 2.857142857142857 {
        const RA: [f64; 8] = [
            -9.86494403484714822705e-03,
            -6.93858572707181764372e-01,
            -1.05586262253232909814e+01,
            -6.23753324503260060396e+01,
            -1.62396669462573470355e+02,
            -1.84605092906711035994e+02,
            -8.12874355063065934246e+01,
            -9.81432934416914548592e+00,
        ];
        const SA: [f64; 9] = [
            1.0,
            1.96512716674392571292e+01,
            1.37657754143519042600e+02,
            4.34565877475229228821e+02,
            6.45387271733267880336e+02,
            4.29008140027567833386e+02,
            1.08635005541779435134e+02,
            6.57024977031928170135e+00,
            -6.04244152148580987438e-02,
        ];
        let r = RA[0] + z * (RA[1] + z * (RA[2] + z * (RA[3] + z * (RA[4] + z * (RA[5] + z * (RA[6] + z * RA[7]))))));
        let s = SA[0] + z * (SA[1] + z * (SA[2] + z * (SA[3] + z * (SA[4] + z * (SA[5] + z * (SA[6] + z * (SA[7] + z * SA[8])))))));
        (r, s)
    } else {
        const RB: [f64; 7] = [
            -9.86494292470009928597e-03,
            -7.99283237680523006574e-01,
            -1.77579549177547519889e+01,
            -1.60636384855821916062e+02,
            -6.37566443368389627722e+02,
            -1.02509513161107724954e+03,
            -4.83519191608651397019e+02,
        ];
        const SB: [f64; 8] = [
            1.0,
            3.03380607434824582924e+01,
            3.25792512996573918826e+02,
            1.53672958608443695994e+03,
            3.19985821950859553908e+03,
            2.55305040643316442583e+03,
            4.74528541206955367215e+02,
            -2.24409524465858183362e+01,
        ];
        let r = RB[0] + z * (RB[1] + z * (RB[2] + z * (RB[3] + z * (RB[4] + z * (RB[5] + z * RB[6])))));
        let s = SB[0] + z * (SB[1] + z * (SB[2] + z * (SB[3] + z * (SB[4] + z * (SB[5] + z * (SB[6] + z * SB[7]))))));
        (r, s)
    };
    let erfc = (-ax * ax - 0.5625).exp() * (r / s).exp() / ax;
    sign * (1.0 - erfc)
}

#[inline]
fn edt_1d(f: &[f64], d: &mut [f64], v: &mut [usize], z: &mut [f64]) {
    const INF: f64 = 1e20;
    let n = f.len();
    if n == 0 {
        return;
    }
    let mut k: isize = -1;
    for q in 0..n {
        if f[q] >= INF {
            continue;
        }
        if k < 0 {
            k = 0;
            v[0] = q;
            z[0] = f64::NEG_INFINITY;
            z[1] = f64::INFINITY;
            continue;
        }
        loop {
            let vk = v[k as usize];
            let s = ((f[q] + (q * q) as f64) - (f[vk] + (vk * vk) as f64))
                / (2.0 * (q as f64 - vk as f64));
            if s <= z[k as usize] {
                k -= 1;
                if k < 0 {
                    k = 0;
                    v[0] = q;
                    z[0] = f64::NEG_INFINITY;
                    z[1] = f64::INFINITY;
                    break;
                }
            } else {
                k += 1;
                v[k as usize] = q;
                z[k as usize] = s;
                z[(k + 1) as usize] = f64::INFINITY;
                break;
            }
        }
    }
    if k < 0 {
        d.fill(INF);
        return;
    }
    let max_k = k as usize;
    let mut idx = 0usize;
    for q in 0..n {
        while idx < max_k && z[idx + 1] < q as f64 {
            idx += 1;
        }
        let dq = q as f64 - v[idx] as f64;
        d[q] = dq * dq + f[v[idx]];
    }
}

/// Exact 2D Euclidean distance transform matching `scipy.ndimage.distance_transform_edt(mask)`.
///
/// Distance is 0 where `mask[i]` is `false`, and Euclidean distance (in cells) to the nearest
/// `false` cell where `mask[i]` is `true`. Returns all 0.0 when no `false` cell exists.
pub fn edt_2d_bool(mask: &[bool], h: usize, w: usize) -> Vec<f64> {
    const INF: f64 = 1e20;
    let n = h * w;
    if n == 0 {
        return Vec::new();
    }
    let any_false = mask.par_iter().any(|&m| !m);
    let any_true = mask.par_iter().any(|&m| m);
    if !any_false || !any_true {
        return vec![0.0; n];
    }

    let mut col_dist = vec![0.0f64; n];
    // Transpose-friendly parallel column pass: compute squared distance vertically per column.
    let col_results: Vec<Vec<f64>> = (0..w)
        .into_par_iter()
        .map(|c| {
            let mut f = vec![0.0f64; h];
            let mut d = vec![0.0f64; h];
            let mut v = vec![0usize; h];
            let mut z = vec![0.0f64; h + 1];
            for r in 0..h {
                f[r] = if mask[r * w + c] { INF } else { 0.0 };
            }
            edt_1d(&f, &mut d, &mut v, &mut z);
            d
        })
        .collect();

    for (c, col) in col_results.iter().enumerate() {
        for r in 0..h {
            col_dist[r * w + c] = col[r];
        }
    }

    let mut out = vec![0.0f64; n];
    out.par_chunks_mut(w)
        .zip(col_dist.par_chunks(w))
        .for_each(|(out_row, in_row)| {
            let mut v = vec![0usize; w];
            let mut z = vec![0.0f64; w + 1];
            edt_1d(in_row, out_row, &mut v, &mut z);
            for val in out_row.iter_mut() {
                *val = val.sqrt();
            }
        });
    out
}

pub fn grid_continentality_dist(is_ocean: &[bool], h: usize, w: usize) -> Vec<f64> {
    let land: Vec<bool> = is_ocean.iter().map(|&o| !o).collect();
    if !land.iter().any(|&l| l) {
        return vec![0.0; h * w];
    }
    let eq_cell_km = 2.0 * std::f64::consts::PI * PLANET_RADIUS_KM / w as f64;
    let scale = (eq_cell_km / PLANET_RADIUS_KM) / CONTINENTALITY_SATURATION_RAD;
    let mut dist = edt_2d_bool(&land, h, w);
    dist.par_iter_mut()
        .for_each(|d| *d = (*d * scale).clamp(0.0, 1.0));
    dist
}

pub fn grid_dist_to_land_rad(is_ocean: &[bool], h: usize, w: usize) -> Vec<f64> {
    let eq_cell_km = 2.0 * std::f64::consts::PI * PLANET_RADIUS_KM / w as f64;
    let scale = eq_cell_km / PLANET_RADIUS_KM;
    let mut dist = edt_2d_bool(is_ocean, h, w);
    dist.par_iter_mut().for_each(|d| *d *= scale);
    dist
}

#[inline]
pub fn seasonal_temp_amplitude(
    lat_deg: f64,
    continentality: f64,
    axial_tilt_deg: f64,
    relief_m: Option<f64>,
) -> f64 {
    let cont = continentality.clamp(0.0, 1.0);
    let tilt_factor = (axial_tilt_deg / SEASON_TILT_REFERENCE_DEG).clamp(0.0, 1.7);
    let lat_term = SEASON_LAT_AMPLITUDE_C * (lat_deg.abs() / 90.0).powf(SEASON_LAT_EXPONENT);
    let mut amp =
        lat_term * (1.0 + SEASON_CONTINENTAL_GAIN * cont) + SEASON_CONTINENTAL_OFFSET_C * cont;
    if let Some(rm) = relief_m {
        let relief_km = rm.max(0.0) / 1000.0;
        amp += (SEASON_RELIEF_GAIN_C_PER_KM * relief_km).min(SEASON_RELIEF_MAX_C);
    }
    tilt_factor * amp.min(SEASON_MAX_AMPLITUDE_C)
}

#[inline]
fn precip_season(lat_deg: f64, season_strength: f64) -> (f64, f64) {
    let lat = lat_deg.abs();
    let equator_ramp = (lat / 8.0).clamp(0.0, 1.0);
    let monsoon_wet = (-((lat - 14.0) / 13.0).powi(2)).exp() * equator_ramp;
    let subtropical_dry_summer = (-((lat - 35.0) / 7.0).powi(2)).exp();
    let summer_share =
        (0.5 + season_strength * (0.36 * monsoon_wet - 0.24 * subtropical_dry_summer))
            .clamp(0.12, 0.9);

    let monsoon = (-((lat - 14.0) / 15.0).powi(2)).exp() * equator_ramp;
    let subtropical = (-((lat - 35.0) / 10.0).powi(2)).exp();
    let conc = (0.10 + season_strength * (0.66 * monsoon + 0.28 * subtropical)).clamp(0.05, 0.9);
    (summer_share, conc)
}

#[inline]
fn months_above(t_annual: f64, amplitude: f64, threshold_c: f64) -> f64 {
    let amp = amplitude.max(1e-6);
    let x = ((threshold_c - t_annual) / amp).clamp(-1.0, 1.0);
    12.0 * x.acos() / std::f64::consts::PI
}

struct SeasonalExtremes {
    t_cold: f64,
    t_warm: f64,
    amplitude: f64,
    driest_month: f64,
    wettest_summer_month: f64,
    driest_summer_month: f64,
    wettest_winter_month: f64,
    summer_share: f64,
}

#[inline]
fn seasonal_extremes(
    t_annual: f64,
    precip_mm: f64,
    lat_deg: f64,
    continentality: f64,
    axial_tilt_deg: f64,
) -> SeasonalExtremes {
    let amplitude = seasonal_temp_amplitude(lat_deg, continentality, axial_tilt_deg, None);
    let t_cold = t_annual - amplitude;
    let t_warm = t_annual + amplitude;
    let season_strength = (axial_tilt_deg / SEASON_TILT_REFERENCE_DEG).clamp(0.0, 1.3);
    let (summer_share, conc) = precip_season(lat_deg, season_strength);
    let monthly_mean = precip_mm / 12.0;
    let summer_month = precip_mm * summer_share / 6.0;
    let winter_month = precip_mm * (1.0 - summer_share) / 6.0;
    let wet_is_summer = summer_share >= 0.5;
    let driest_other = (if wet_is_summer {
        winter_month
    } else {
        summer_month
    }) * (1.0 - conc);
    let driest_month = driest_other.min(monthly_mean * (1.0 - conc));
    let wettest_summer_month =
        summer_month * (1.0 + 4.0 * conc * if wet_is_summer { 1.0 } else { 0.0 });
    let driest_summer_month = summer_month * (1.0 - conc);
    let wettest_winter_month =
        winter_month * (1.0 + 4.0 * conc * if wet_is_summer { 0.0 } else { 1.0 });
    SeasonalExtremes {
        t_cold,
        t_warm,
        amplitude,
        driest_month,
        wettest_summer_month,
        driest_summer_month,
        wettest_winter_month,
        summer_share,
    }
}

#[inline]
pub fn classify_koppen_cell(
    t_annual: f64,
    precip_mm: f64,
    lat_deg: f64,
    continentality: f64,
    axial_tilt_deg: f64,
) -> u8 {
    let ex = seasonal_extremes(t_annual, precip_mm, lat_deg, continentality, axial_tilt_deg);
    let pth_offset = if ex.summer_share >= 0.7 {
        280.0
    } else if ex.summer_share <= 0.3 {
        0.0
    } else {
        140.0
    };
    let pth = 20.0 * t_annual + pth_offset;
    let is_b = precip_mm < pth;
    let is_bw = is_b && (precip_mm < 0.5 * pth);
    let is_hot = t_annual >= TROPICAL_TEMP_C;

    if is_bw && is_hot {
        return 4; // BWh
    }
    if is_bw {
        return 5; // BWk
    }
    if is_b && is_hot {
        return 6; // BSh
    }
    if is_b {
        return 7; // BSk
    }

    let is_e = ex.t_warm < 10.0;
    let is_ef = is_e && (ex.t_warm < 0.0);
    if is_ef {
        return 30; // EF
    }
    if is_e {
        return 29; // ET
    }

    let is_a = ex.t_cold >= 18.0;
    if is_a {
        if ex.driest_month >= 60.0 {
            return 0; // Af
        }
        if ex.driest_month >= 100.0 - precip_mm / 25.0 {
            return 1; // Am
        }
        if ex.summer_share >= 0.5 {
            return 2; // Aw
        }
        return 3; // As
    }

    let abs_lat = lat_deg.abs();
    let mediterranean_band = abs_lat >= 27.0 && abs_lat <= 46.0 && continentality < 0.4;
    let summer_dry = mediterranean_band
        && precip_mm < 900.0
        && ex.summer_share < 0.45
        && ex.driest_summer_month < 40.0
        && ex.driest_summer_month < ex.wettest_winter_month / 3.0;
    let winter_dry =
        ex.driest_month < ex.wettest_summer_month / 10.0 && ex.summer_share >= 0.62;

    let warm_a = ex.t_warm >= 22.0;
    let months_10 = months_above(t_annual, ex.amplitude, 10.0);
    let warm_b = !warm_a && months_10 >= 4.0;
    let very_cold_d = ex.t_cold < -38.0;

    let is_d = ex.t_cold < 0.0;
    let third = if warm_a {
        0u8
    } else if warm_b {
        1u8
    } else {
        2u8
    };

    if !is_d {
        if summer_dry {
            return 8 + third; // Csa, Csb, Csc
        }
        if winter_dry {
            return 11 + third; // Cwa, Cwb, Cwc
        }
        return 14 + third; // Cfa, Cfb, Cfc
    }

    if summer_dry {
        return if third < 2 {
            17 + third
        } else if very_cold_d {
            20
        } else {
            19
        };
    }
    if winter_dry {
        return if third < 2 {
            21 + third
        } else if very_cold_d {
            24
        } else {
            23
        };
    }
    if third < 2 {
        25 + third
    } else if very_cold_d {
        28
    } else {
        27
    }
}

#[inline]
pub fn classify_pelagic_cell(
    sst: f64,
    lat_deg: f64,
    dist_to_land_rad: f64,
    has_sea_ice: bool,
) -> u8 {
    if has_sea_ice {
        return POLAR_SEA_ICE;
    }
    let abs_lat = lat_deg.abs();
    let is_shelf = dist_to_land_rad <= SHELF_RANGE_RAD;
    let [e0, e1, e2, e3] = PELAGIC_SST_EDGES;
    if sst < e0 {
        return 39; // Polar Ocean
    }
    if sst < e1 {
        return if is_shelf { 38 } else { 37 };
    }
    if sst < e2 {
        return if is_shelf { 36 } else { 35 };
    }
    if is_shelf {
        return 34; // Tropical Coastal Waters
    }
    if abs_lat <= EQUATORIAL_LAT_DEG {
        return 33; // Equatorial Divergence
    }
    if sst < e3 {
        return 32; // Subtropical Gyre
    }
    31 // Tropical Open Ocean
}

#[inline]
fn soft_lt(x: f64, threshold: f64, scale: f64) -> f64 {
    let s = scale.max(1e-9);
    0.5 * (1.0 + erf_f64((threshold - x) / (s * std::f64::consts::SQRT_2)))
}

#[inline]
fn soft_ge(x: f64, threshold: f64, scale: f64) -> f64 {
    1.0 - soft_lt(x, threshold, scale)
}

#[inline]
fn soft_precip_scale(threshold: f64) -> f64 {
    (SOFT_PRECIP_FRACTION * threshold.abs()).max(SOFT_PRECIP_FLOOR_MM)
}

#[inline]
pub fn soft_sea_ice_weight_scalar(fraction: f64, threshold: f64) -> f64 {
    soft_ge(fraction, threshold, SOFT_SEA_ICE_FRACTION_BLEND)
}

const KOPPEN_CASCADE_CODES: [usize; 31] = [
    4, 5, 6, 7, 30, 29, 0, 1, 2, 3, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
    24, 25, 26, 27, 28,
];

#[inline]
fn soft_koppen_weights(
    t_annual: f64,
    precip_mm: f64,
    lat_deg: f64,
    continentality: f64,
    axial_tilt_deg: f64,
) -> [f64; N_BIOMES] {
    let ex = seasonal_extremes(t_annual, precip_mm, lat_deg, continentality, axial_tilt_deg);
    let tc = SOFT_TEMP_BLEND_C;

    let w_wet = soft_ge(ex.summer_share, 0.7, SOFT_SHARE_BLEND);
    let w_dry = soft_ge(0.3, ex.summer_share, SOFT_SHARE_BLEND);
    let w_even = (1.0 - w_wet - w_dry).clamp(0.0, 1.0);
    let pth_offset = w_wet * 280.0 + w_even * 140.0;
    let pth = 20.0 * t_annual + pth_offset;
    let is_b = soft_lt(precip_mm, pth, soft_precip_scale(pth));
    let is_bw = is_b * soft_lt(precip_mm, 0.5 * pth, soft_precip_scale(0.5 * pth));
    let is_hot = soft_ge(t_annual, TROPICAL_TEMP_C, tc);

    let is_e = soft_lt(ex.t_warm, 10.0, tc);
    let is_ef = is_e * soft_lt(ex.t_warm, 0.0, tc);

    let is_a = soft_ge(ex.t_cold, 18.0, tc) * (1.0 - is_b);
    let af = is_a * soft_ge(ex.driest_month, 60.0, SOFT_PRECIP_MONTH_BLEND_MM);
    let am = is_a
        * (1.0 - af)
        * soft_ge(
            ex.driest_month,
            100.0 - precip_mm / 25.0,
            SOFT_PRECIP_MONTH_BLEND_MM,
        );
    let a_dry = is_a * (1.0 - af) * (1.0 - am);
    let aw = a_dry * soft_ge(ex.summer_share, 0.5, SOFT_SHARE_BLEND);
    let a_s = a_dry * (1.0 - aw);

    let abs_lat = lat_deg.abs();
    let lat_band = soft_ge(abs_lat, 27.0, SOFT_LAT_BLEND_DEG)
        * soft_lt(abs_lat, 46.0, SOFT_LAT_BLEND_DEG);
    let med_band = lat_band * soft_lt(continentality, 0.4, SOFT_CONTINENTALITY_BLEND);
    let summer_dry = med_band
        * soft_lt(precip_mm, 900.0, soft_precip_scale(900.0))
        * soft_lt(ex.summer_share, 0.45, SOFT_SHARE_BLEND)
        * soft_lt(ex.driest_summer_month, 40.0, SOFT_PRECIP_MONTH_BLEND_MM)
        * soft_lt(
            ex.driest_summer_month,
            ex.wettest_winter_month / 3.0,
            SOFT_PRECIP_MONTH_BLEND_MM,
        );
    let winter_dry = soft_lt(
        ex.driest_month,
        ex.wettest_summer_month / 10.0,
        SOFT_PRECIP_MONTH_BLEND_MM,
    ) * soft_ge(ex.summer_share, 0.62, SOFT_SHARE_BLEND);

    let warm_a = soft_ge(ex.t_warm, 22.0, tc);
    let months_10 = months_above(t_annual, ex.amplitude, 10.0);
    let warm_b = (1.0 - warm_a) * soft_ge(months_10, 4.0, SOFT_MONTHS_BLEND);
    let very_cold_d = soft_lt(ex.t_cold, -38.0, tc);

    let is_cd = (1.0 - is_a) * (1.0 - is_b) * (1.0 - is_e) * soft_ge(ex.t_warm, 10.0, tc);
    let is_d = is_cd * soft_lt(ex.t_cold, 0.0, tc);
    let is_c = is_cd * (1.0 - is_d);

    let not_warm_a = 1.0 - warm_a;
    let not_warm_b = 1.0 - warm_b;
    let not_very_cold = 1.0 - very_cold_d;
    let not_summer_dry = 1.0 - summer_dry;
    let not_winter_dry = 1.0 - winter_dry;

    let c_s = is_c * summer_dry;
    let c_w = is_c * not_summer_dry * winter_dry;
    let c_f = is_c * not_summer_dry * not_winter_dry;
    let d_s = is_d * summer_dry;
    let d_w = is_d * not_summer_dry * winter_dry;
    let d_f = is_d * not_summer_dry * not_winter_dry;

    let conds: [f64; 31] = [
        is_bw * is_hot,
        is_bw,
        is_b * is_hot,
        is_b,
        is_ef,
        is_e,
        af,
        am,
        aw,
        a_s,
        c_s * warm_a,
        c_s * not_warm_a * warm_b,
        c_s * not_warm_a * not_warm_b,
        c_w * warm_a,
        c_w * not_warm_a * warm_b,
        c_w * not_warm_a * not_warm_b,
        c_f * warm_a,
        c_f * not_warm_a * warm_b,
        c_f * not_warm_a * not_warm_b,
        d_s * warm_a,
        d_s * not_warm_a * warm_b,
        d_s * not_warm_a * not_warm_b * not_very_cold,
        d_s * not_warm_a * not_warm_b * very_cold_d,
        d_w * warm_a,
        d_w * not_warm_a * warm_b,
        d_w * not_warm_a * not_warm_b * not_very_cold,
        d_w * not_warm_a * not_warm_b * very_cold_d,
        d_f * warm_a,
        d_f * not_warm_a * warm_b,
        d_f * not_warm_a * not_warm_b * not_very_cold,
        d_f * not_warm_a * not_warm_b * very_cold_d,
    ];

    let mut weights = [0.0f64; N_BIOMES];
    let mut remaining = 1.0f64;
    for (cond, &choice) in conds.iter().zip(KOPPEN_CASCADE_CODES.iter()) {
        let c = cond.clamp(0.0, 1.0);
        weights[choice] += c * remaining;
        remaining *= 1.0 - c;
    }
    weights[16] += remaining; // Cfc default
    weights
}

const PELAGIC_CASCADE_CODES: [usize; 10] = [40, 39, 38, 37, 36, 35, 34, 33, 32, 31];

#[inline]
fn soft_pelagic_weights(
    sst: f64,
    lat_deg: f64,
    dist_to_land_rad: f64,
    sea_ice_weight: f64,
) -> [f64; N_BIOMES] {
    let abs_lat = lat_deg.abs();
    let is_shelf = soft_lt(dist_to_land_rad, SHELF_RANGE_RAD, SOFT_SHELF_DIST_BLEND_RAD);
    let [e0, e1, e2, e3] = PELAGIC_SST_EDGES;
    let tc = SOFT_TEMP_BLEND_C;
    let polar = soft_lt(sst, e0, tc);
    let cold_temp = soft_ge(sst, e0, tc) * soft_lt(sst, e1, tc);
    let temp = soft_ge(sst, e1, tc) * soft_lt(sst, e2, tc);
    let subtropical = soft_ge(sst, e2, tc) * soft_lt(sst, e3, tc);
    let tropical = soft_ge(sst, e2, tc);

    let conds: [f64; 10] = [
        sea_ice_weight,
        polar,
        cold_temp * is_shelf,
        cold_temp,
        temp * is_shelf,
        temp,
        tropical * is_shelf,
        tropical * soft_lt(abs_lat, EQUATORIAL_LAT_DEG, SOFT_LAT_BLEND_DEG),
        subtropical,
        tropical,
    ];

    let mut weights = [0.0f64; N_BIOMES];
    let mut remaining = 1.0f64;
    for (cond, &choice) in conds.iter().zip(PELAGIC_CASCADE_CODES.iter()) {
        let c = cond.clamp(0.0, 1.0);
        weights[choice] += c * remaining;
        remaining *= 1.0 - c;
    }
    weights[35] += remaining; // Temperate Open Ocean default
    weights
}

#[inline]
pub fn soft_top2_cell(
    t_c: f64,
    precip_mm: f64,
    is_ocean: bool,
    lat_deg: f64,
    axial_tilt_deg: f64,
    continentality: f64,
    dist_to_land_rad: f64,
    sea_ice_weight: f64,
    winner_id: u8,
    overridden: bool,
) -> (f64, u8) {
    let mut weights = if is_ocean {
        soft_pelagic_weights(t_c, lat_deg, dist_to_land_rad, sea_ice_weight)
    } else {
        soft_koppen_weights(t_c, precip_mm, lat_deg, continentality, axial_tilt_deg)
    };
    let w_idx = winner_id as usize;
    let winner_w = weights[w_idx];
    weights[w_idx] = -1.0;

    let mut runner_up = 0usize;
    let mut runner_up_w = f64::NEG_INFINITY;
    for (idx, &w) in weights.iter().enumerate() {
        if w > runner_up_w {
            runner_up_w = w;
            runner_up = idx;
        }
    }
    let share = if overridden {
        1.0
    } else {
        let denom = (winner_w + runner_up_w).max(1e-9);
        (winner_w / denom).clamp(MIN_DOMINANT_SHARE, 1.0)
    };
    (share, runner_up as u8)
}

pub fn neighbor_vote_2pass(biome_ids: &mut [u8], is_ocean: &[bool], h: usize, w: usize) -> Vec<bool> {
    let n = h * w;
    if n == 0 {
        return Vec::new();
    }
    let offsets: [(isize, isize); 8] = [
        (-1, -1),
        (-1, 0),
        (-1, 1),
        (0, -1),
        (0, 1),
        (1, -1),
        (1, 0),
        (1, 1),
    ];

    // Compute eligibility once before the 2 iterations, matching biomes.py.
    let eligible: Vec<bool> = (0..h)
        .into_par_iter()
        .flat_map_iter(|r| {
            let mut row_elig = Vec::with_capacity(w);
            for c in 0..w {
                let idx = r * w + c;
                if is_ocean[idx] {
                    row_elig.push(false);
                    continue;
                }
                let self_id = biome_ids[idx];
                let mut disagree = 0usize;
                for &(dy, dx) in &offsets {
                    let nr = (r as isize + dy).rem_euclid(h as isize) as usize;
                    let nc = (c as isize + dx).rem_euclid(w as isize) as usize;
                    let nidx = nr * w + nc;
                    if !is_ocean[nidx] && biome_ids[nidx] != self_id {
                        disagree += 1;
                    }
                }
                row_elig.push(disagree >= BIOME_RAGGED_MIN_DISAGREEING_NEIGHBOURS);
            }
            row_elig
        })
        .collect();

    let mut overridden = vec![false; n];
    let mut next_ids = biome_ids.to_vec();

    for _ in 0..BIOME_VOTE_ITERATIONS {
        let prev_ids: &[u8] = biome_ids;
        let row_updates: Vec<(Vec<u8>, Vec<bool>)> = (0..h)
            .into_par_iter()
            .map(|r| {
                let mut row_ids = Vec::with_capacity(w);
                let mut row_chg = Vec::with_capacity(w);
                for c in 0..w {
                    let idx = r * w + c;
                    let cur_id = prev_ids[idx];
                    if !eligible[idx] || is_ocean[idx] {
                        row_ids.push(cur_id);
                        row_chg.push(false);
                        continue;
                    }
                    let mut nbs = [0u8; 8];
                    let mut valid = 0usize;
                    for &(dy, dx) in &offsets {
                        let nr = (r as isize + dy).rem_euclid(h as isize) as usize;
                        let nc = (c as isize + dx).rem_euclid(w as isize) as usize;
                        let nidx = nr * w + nc;
                        if !is_ocean[nidx] {
                            nbs[valid] = prev_ids[nidx];
                            valid += 1;
                        }
                    }
                    if valid == 0 {
                        row_ids.push(cur_id);
                        row_chg.push(false);
                        continue;
                    }
                    nbs[..valid].sort_unstable();
                    let mut best_count = 0usize;
                    let mut best_biome = cur_id;
                    let mut run = 0usize;
                    let mut prev_nb = u8::MAX;
                    for &nb in &nbs[..valid] {
                        if nb == prev_nb {
                            run += 1;
                        } else {
                            run = 1;
                            prev_nb = nb;
                        }
                        if run > best_count {
                            best_count = run;
                            best_biome = nb;
                        }
                    }
                    let threshold = (BIOME_VOTE_MIN_NEIGHBOUR_FRACTION * valid as f64).ceil() as usize;
                    if best_count >= threshold && best_biome != cur_id {
                        row_ids.push(best_biome);
                        row_chg.push(true);
                    } else {
                        row_ids.push(cur_id);
                        row_chg.push(false);
                    }
                }
                (row_ids, row_chg)
            })
            .collect();

        for (r, (row_ids, row_chg)) in row_updates.into_iter().enumerate() {
            let base = r * w;
            next_ids[base..base + w].copy_from_slice(&row_ids);
            for c in 0..w {
                if row_chg[c] {
                    overridden[base + c] = true;
                }
            }
        }
        biome_ids.copy_from_slice(&next_ids);
    }

    overridden
}

pub fn smooth_biome_field_blend_core(
    temperature_c: &[f64],
    precipitation_mm: &[f64],
    is_ocean: &[bool],
    lat_deg: &[f64],
    h: usize,
    w: usize,
    axial_tilt_deg: f64,
    has_sea_ice: Option<&[bool]>,
    sea_ice_weight: Option<&[f64]>,
) -> (Vec<u8>, Vec<f64>, Vec<u8>) {
    let n = h * w;
    let continentality = grid_continentality_dist(is_ocean, h, w);
    let dist_to_land = grid_dist_to_land_rad(is_ocean, h, w);

    let mut biome_ids: Vec<u8> = (0..n)
        .into_par_iter()
        .map(|i| {
            let lat = lat_deg[i / w];
            if is_ocean[i] {
                let ice = has_sea_ice.map_or(false, |s| s[i]);
                classify_pelagic_cell(temperature_c[i], lat, dist_to_land[i], ice)
            } else {
                classify_koppen_cell(
                    temperature_c[i],
                    precipitation_mm[i],
                    lat,
                    continentality[i],
                    axial_tilt_deg,
                )
            }
        })
        .collect();

    let overridden = neighbor_vote_2pass(&mut biome_ids, is_ocean, h, w);

    let mut winner_share = vec![0.0f64; n];
    let mut runner_up_ids = vec![0u8; n];
    winner_share
        .par_iter_mut()
        .zip(runner_up_ids.par_iter_mut())
        .enumerate()
        .for_each(|(i, (sh, ru))| {
            let lat = lat_deg[i / w];
            let ice_w = if let Some(sw) = sea_ice_weight {
                sw[i]
            } else if let Some(hs) = has_sea_ice {
                if hs[i] { 1.0 } else { 0.0 }
            } else {
                0.0
            };
            let (s, r) = soft_top2_cell(
                temperature_c[i],
                precipitation_mm[i],
                is_ocean[i],
                lat,
                axial_tilt_deg,
                continentality[i],
                dist_to_land[i],
                ice_w,
                biome_ids[i],
                overridden[i],
            );
            *sh = s;
            *ru = r;
        });

    (biome_ids, winner_share, runner_up_ids)
}

/// Rank-based shaded-relief factor per cell matching `scipy.stats.rankdata(..., method="average")`.
pub fn biome_relative_shade_factor(biome_ids: &[u8], elevation_m: &[f64]) -> Vec<f64> {
    let n = biome_ids.len();
    let mut buckets: Vec<Vec<(f64, usize)>> = vec![Vec::new(); N_BIOMES];
    for (i, (&bid, &elev)) in biome_ids.iter().zip(elevation_m.iter()).enumerate() {
        let b = bid as usize;
        if b < N_BIOMES {
            buckets[b].push((elev, i));
        }
    }

    let bucket_factors: Vec<Vec<(usize, f64)>> = buckets
        .into_par_iter()
        .map(|mut items| {
            let m = items.len();
            if m <= 1 {
                return items.into_iter().map(|(_, idx)| (idx, 1.0)).collect();
            }
            items.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
            let mut out = Vec::with_capacity(m);
            let denom = (m - 1) as f64;
            let span = BIOME_SHADE_MAX - BIOME_SHADE_MIN;
            let mut start = 0usize;
            while start < m {
                let val = items[start].0;
                let mut end = start + 1;
                while end < m && items[end].0 == val {
                    end += 1;
                }
                // Average 1-based rank over [start+1 .. end] is 1 + 0.5 * (start + end - 1),
                // so (avg_rank - 1) / (m - 1) = 0.5 * (start + end - 1) / (m - 1).
                let frac = 0.5 * ((start + end - 1) as f64) / denom;
                let factor = BIOME_SHADE_MIN + span * frac;
                for &(_, idx) in &items[start..end] {
                    out.push((idx, factor));
                }
                start = end;
            }
            out
        })
        .collect();

    let mut factors = vec![1.0f64; n];
    for entries in bucket_factors {
        for (idx, f) in entries {
            factors[idx] = f;
        }
    }
    factors
}

// ---------------------------------------------------------------------------
// Sea-ice & climate grid helpers matching sea_ice.py / climate.py / hydrology.py
// ---------------------------------------------------------------------------

pub const SEAWATER_FREEZING_C: f64 = -1.8;
pub const OCEAN_HEAT_FLUX_W_M2: f64 = 2.0;
pub const POLAR_SEA_ICE_FRACTION: f64 = 0.5;
pub const PERENNIAL_COVER_LOW: f64 = 0.65;
pub const PERENNIAL_COVER_HIGH: f64 = 1.25;
pub const SUNLIGHT: f64 = 1.0;
pub const INSOLATION_FLOOR: f64 = 0.03;
pub const AXIAL_TILT_DECLINATION_SAMPLES: usize = 24;
pub const LAND_TEMP_MIN_C: f64 = -28.0;
pub const LAND_TEMP_RANGE_C: f64 = 55.0;
pub const WATER_TEMP_MIN_C: f64 = -2.0;
pub const WATER_TEMP_RANGE_C: f64 = 30.0;

const RHO_ICE: f64 = 917.0;
const RHO_WATER: f64 = 1000.0;
const ICE_THERMAL_CONDUCTIVITY_W_M_K: f64 = 2.1;
const ICE_LATENT_HEAT_J_KG: f64 = 334_000.0;
const SECONDS_PER_YEAR: f64 = 3.15576e7;
const W_M2_TO_M_ICE_PER_YR: f64 = SECONDS_PER_YEAR / (RHO_ICE * ICE_LATENT_HEAT_J_KG);
const DDF_ICE_MM_WATER_PER_DAY_K: f64 = 8.0;

#[inline]
pub fn compute_heating_index_scalar(
    lat_deg: f64,
    axial_tilt_deg: f64,
    solar_multiplier: f64,
) -> f64 {
    let sunlight = SUNLIGHT * solar_multiplier;
    let lat_r = lat_deg.to_radians();
    if axial_tilt_deg <= 1e-6 {
        return lat_r.cos().clamp(INSOLATION_FLOOR, 1.0) * sunlight;
    }
    let mut sum = 0.0f64;
    let n = AXIAL_TILT_DECLINATION_SAMPLES;
    for k in 0..n {
        let dec_deg = -axial_tilt_deg + (2.0 * axial_tilt_deg) * (k as f64) / ((n - 1) as f64);
        let cz = (lat_r - dec_deg.to_radians()).cos().max(0.0);
        sum += cz;
    }
    let mean = sum / (n as f64);
    mean.clamp(INSOLATION_FLOOR, 1.0) * sunlight
}

#[inline]
pub fn seasonal_positive_degrees_scalar(mean_c: f64, amplitude_c: f64, threshold_c: f64) -> f64 {
    let d = mean_c - threshold_c;
    let amp = amplitude_c.max(0.0);
    if d.abs() < amp {
        let phi = (-d / amp).acos();
        (d * phi + (amp * amp - d * d).sqrt()) / std::f64::consts::PI
    } else {
        d.max(0.0)
    }
}

#[inline]
pub fn seasonal_fraction_above_scalar(mean_c: f64, amplitude_c: f64, threshold_c: f64) -> f64 {
    let d = mean_c - threshold_c;
    let amp = amplitude_c.max(0.0);
    let frac = if amp > 0.0 && d.abs() < amp {
        (-d / amp).clamp(-1.0, 1.0).acos() / std::f64::consts::PI
    } else if d > 0.0 {
        1.0
    } else if d == 0.0 && amp == 0.0 {
        0.5
    } else {
        0.0
    };
    frac.clamp(0.0, 1.0)
}

#[derive(Clone, Copy)]
pub struct SeaIceCell {
    pub temperature_c: f64,
    pub perennial: bool,
    pub fraction: f64,
    pub perennial_cover: f64,
}

#[inline]
pub fn sea_ice_annual_balance_scalar(
    ocean_temperature_c: f64,
    lat_deg: f64,
    axial_tilt_deg: f64,
    solar_multiplier: f64,
    is_ocean: bool,
) -> SeaIceCell {
    let heating = compute_heating_index_scalar(lat_deg, axial_tilt_deg, solar_multiplier);
    let ocean_baseline_c = WATER_TEMP_MIN_C + WATER_TEMP_RANGE_C * heating;
    let sea_level_air_c = LAND_TEMP_MIN_C + LAND_TEMP_RANGE_C * heating;
    let temperature_c = ocean_temperature_c + (sea_level_air_c - ocean_baseline_c);
    let amplitude_c = seasonal_temp_amplitude(lat_deg, 0.0, axial_tilt_deg, None);

    let freezing_days =
        365.25 * seasonal_positive_degrees_scalar(-SEAWATER_FREEZING_C - (-temperature_c), amplitude_c, 0.0);
    //Wait: seasonal_freezing_degrees(mean, amp, thresh) = seasonal_positive_degrees(thresh - mean, amp, 0.0)
    // Check: above expression simplified to thresh - mean! Let's write it explicitly to be crystal clear:
    let _ = freezing_days;
    let freezing_days =
        365.25 * seasonal_positive_degrees_scalar(SEAWATER_FREEZING_C - temperature_c, amplitude_c, 0.0);
    let thawing_days = 365.25 * seasonal_positive_degrees_scalar(temperature_c, amplitude_c, 0.0);
    let growth_m = (2.0 * ICE_THERMAL_CONDUCTIVITY_W_M_K * freezing_days * 86_400.0
        / (RHO_ICE * ICE_LATENT_HEAT_J_KG))
        .sqrt();
    let surface_loss_m =
        DDF_ICE_MM_WATER_PER_DAY_K * 1e-3 * thawing_days * RHO_WATER / RHO_ICE;
    let bottom_loss_m = OCEAN_HEAT_FLUX_W_M2 * W_M2_TO_M_ICE_PER_YR;
    let total_loss_m = surface_loss_m + bottom_loss_m;

    let freezes = is_ocean && freezing_days > 0.0;
    let perennial = freezes && (growth_m > total_loss_m);
    let seasonal = freezes && !perennial;
    let winter_fraction =
        seasonal_fraction_above_scalar(-temperature_c, amplitude_c, -SEAWATER_FREEZING_C);
    let fraction = if perennial {
        1.0
    } else if seasonal {
        winter_fraction
    } else {
        0.0
    };
    let ratio = growth_m / total_loss_m;
    let ramp = ((ratio - PERENNIAL_COVER_LOW) / (PERENNIAL_COVER_HIGH - PERENNIAL_COVER_LOW))
        .clamp(0.0, 1.0);
    let perennial_cover = if freezes {
        ramp * ramp * (3.0 - 2.0 * ramp)
    } else {
        0.0
    };
    SeaIceCell {
        temperature_c,
        perennial,
        fraction,
        perennial_cover,
    }
}

#[inline]
pub fn bilinear_sample_climate(
    field: &[f64],
    src_h: usize,
    src_w: usize,
    lat_deg: f64,
    lon_deg: f64,
) -> f64 {
    let row_step = 180.0 / src_h as f64;
    let row_f = ((90.0 - lat_deg) / row_step - 0.5).clamp(0.0, (src_h - 1) as f64);
    let row0 = row_f.floor() as usize;
    let row1 = (row0 + 1).min(src_h - 1);
    let row_t = row_f - row0 as f64;

    let col_step = 360.0 / src_w as f64;
    let lon_norm = (lon_deg + 180.0).rem_euclid(360.0) - 180.0;
    let col_f = (lon_norm + 180.0) / col_step - 0.5;
    let col0 = (col_f.floor() as isize).rem_euclid(src_w as isize) as usize;
    let col1 = (col0 + 1) % src_w;
    let col_t = col_f - col_f.floor();

    let f00 = field[row0 * src_w + col0];
    let f01 = field[row0 * src_w + col1];
    let f10 = field[row1 * src_w + col0];
    let f11 = field[row1 * src_w + col1];
    let top = f00 * (1.0 - col_t) + f01 * col_t;
    let bot = f10 * (1.0 - col_t) + f11 * col_t;
    top * (1.0 - row_t) + bot * row_t
}

#[inline]
pub fn grid_cell_indices_scalar(lat_deg: f64, lon_deg: f64, h: usize, w: usize) -> (usize, usize) {
    let row = (((90.0 - lat_deg) / (180.0 / h as f64)).floor() as isize)
        .clamp(0, (h - 1) as isize) as usize;
    let col = (((lon_deg + 180.0) / (360.0 / w as f64)).floor() as isize)
        .rem_euclid(w as isize) as usize;
    (row, col)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erf_known_values() {
        assert!((erf_f64(0.0)).abs() < 1e-15);
        assert!((erf_f64(1.0) - 0.8427007929497149).abs() < 1e-14);
        assert!((erf_f64(-1.0) + 0.8427007929497149).abs() < 1e-14);
        assert!((erf_f64(2.5) - 0.999593047982555).abs() < 1e-14);
    }

    #[test]
    fn edt_matches_3x4_grid() {
        // False at (0, 0), True everywhere else.
        let mut mask = vec![true; 12];
        mask[0] = false;
        let d = edt_2d_bool(&mask, 3, 4);
        assert_eq!(d[0], 0.0);
        assert!((d[1] - 1.0).abs() < 1e-12);
        assert!((d[4 + 1] - 2.0f64.sqrt()).abs() < 1e-12);
        assert!((d[2 * 4 + 3] - (4.0f64 + 9.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn neighbor_vote_outvotes_lone_speckle() {
        let mut ids = vec![14u8; 25];
        ids[12] = 4; // lone speckle in center of 5x5
        let ocean = vec![false; 25];
        let ov = neighbor_vote_2pass(&mut ids, &ocean, 5, 5);
        assert_eq!(ids[12], 14);
        assert!(ov[12]);
    }
}
