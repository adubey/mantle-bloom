use std::f64::consts::PI;
use std::sync::{Arc, Mutex, OnceLock};

use rayon::prelude::*;

use crate::render::raster::CELL_OVERLAP_FACTOR;

pub const BEHRMANN_STANDARD_PARALLEL: f64 = 30.0 * PI / 180.0;
const ECKERT4_NEWTON_MAX_ITERS: usize = 12;
const ECKERT4_NEWTON_TOL: f64 = 1e-14;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Projection {
    Behrmann,
    Eckert4,
    Globe,
}

impl Projection {
    pub fn parse(name: &str) -> Result<Self, String> {
        match name {
            "behrmann" => Ok(Self::Behrmann),
            "eckert4" => Ok(Self::Eckert4),
            "globe" => Ok(Self::Globe),
            _ => Err(format!("Unknown projection {name:?}")),
        }
    }
}

#[inline]
pub fn xyz_to_latlon(p: [f64; 3]) -> (f64, f64) {
    let lat = p[2].clamp(-1.0, 1.0).asin();
    let lon = p[1].atan2(p[0]);
    (lat, lon)
}

#[inline]
pub fn latlon_to_xyz(lat: f64, lon: f64) -> [f64; 3] {
    let cl = lat.cos();
    [cl * lon.cos(), cl * lon.sin(), lat.sin()]
}

#[inline]
pub fn rotate_point(p: [f64; 3], rot: &[f64; 9]) -> [f64; 3] {
    [
        p[0] * rot[0] + p[1] * rot[1] + p[2] * rot[2],
        p[0] * rot[3] + p[1] * rot[4] + p[2] * rot[5],
        p[0] * rot[6] + p[1] * rot[7] + p[2] * rot[8],
    ]
}

#[inline]
pub fn eckert4_theta(lat: f64) -> f64 {
    let rhs_scale = 2.0 + PI / 2.0;
    if (lat.abs() - PI / 2.0).abs() <= 1e-8 + 1e-5 * (PI / 2.0) {
        return lat.signum() * (PI / 2.0);
    }
    let target = rhs_scale * lat.sin();
    let mut theta = lat * 0.5;
    for _ in 0..ECKERT4_NEWTON_MAX_ITERS {
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        let f = theta + sin_t * cos_t + 2.0 * sin_t - target;
        let fprime = 2.0 * cos_t * (cos_t + 1.0);
        let safe_fp = if fprime.abs() < 1e-9 { 1e-9 } else { fprime };
        let step = f / safe_fp;
        theta -= step;
        if step.abs() < ECKERT4_NEWTON_TOL {
            break;
        }
    }
    theta
}

#[inline]
pub fn project_latlon(proj: Projection, lat: f64, lon: f64) -> [f64; 2] {
    match proj {
        Projection::Behrmann => {
            let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
            [lon * cos_phi0, lat.sin() / cos_phi0]
        }
        Projection::Eckert4 => {
            let c = 2.0 / (PI * (4.0 + PI)).sqrt();
            let theta = eckert4_theta(lat);
            [c * lon * (1.0 + theta.cos()), c * PI * theta.sin()]
        }
        Projection::Globe => [lat.cos() * lon.sin(), lat.sin()],
    }
}

#[inline]
pub fn project_point(proj: Projection, p: [f64; 3]) -> [f64; 2] {
    let (lat, lon) = xyz_to_latlon(p);
    project_latlon(proj, lat, lon)
}

#[inline]
pub fn project_offset(proj: Projection, p: [f64; 3], reference_lon: f64) -> [f64; 2] {
    let (lat, lon) = xyz_to_latlon(p);
    let wrapped = reference_lon + ((lon - reference_lon + PI).rem_euclid(2.0 * PI) - PI);
    project_latlon(proj, lat, wrapped)
}

pub fn sphere_extent(proj: Projection) -> [f64; 4] {
    static EXTENTS: OnceLock<[[f64; 4]; 3]> = OnceLock::new();
    let table = EXTENTS.get_or_init(|| {
        let compute = |pr: Projection| {
            let mut min_x = f64::INFINITY;
            let mut min_y = f64::INFINITY;
            let mut max_x = f64::NEG_INFINITY;
            let mut max_y = f64::NEG_INFINITY;
            for r in 0..181 {
                let lat = (-90.0 + r as f64).to_radians();
                for c in 0..361 {
                    let lon = (-180.0 + c as f64).to_radians();
                    let [x, y] = project_latlon(pr, lat, lon);
                    if x < min_x {
                        min_x = x;
                    }
                    if x > max_x {
                        max_x = x;
                    }
                    if y < min_y {
                        min_y = y;
                    }
                    if y > max_y {
                        max_y = y;
                    }
                }
            }
            [min_x, min_y, max_x, max_y]
        };
        [
            compute(Projection::Behrmann),
            compute(Projection::Eckert4),
            compute(Projection::Globe),
        ]
    });
    match proj {
        Projection::Behrmann => table[0],
        Projection::Eckert4 => table[1],
        Projection::Globe => table[2],
    }
}

#[inline]
pub fn fit_view(extent: [f64; 4], width: usize, height: usize, padding_px: f64) -> (f64, f64, f64) {
    let [min_x, min_y, max_x, max_y] = extent;
    let data_w = (max_x - min_x).max(1e-9);
    let data_h = (max_y - min_y).max(1e-9);
    let w = width as f64;
    let h = height as f64;
    let scale = ((w - 2.0 * padding_px) / data_w).min((h - 2.0 * padding_px) / data_h);
    let offset_x = w * 0.5 - scale * (min_x + max_x) * 0.5;
    let offset_y = h * 0.5 + scale * (min_y + max_y) * 0.5;
    (scale, offset_x, offset_y)
}

#[inline]
pub fn to_pixel(scale: f64, offset_x: f64, offset_y: f64, xy: [f64; 2]) -> [f64; 2] {
    [scale * xy[0] + offset_x, -scale * xy[1] + offset_y]
}

#[derive(Clone, PartialEq, Eq)]
struct GridProjKey {
    grid_h: usize,
    grid_w: usize,
    proj_id: u8,
    rot_bits: [u64; 9],
    width: usize,
    height: usize,
    padding_micros: i64,
}

pub struct ProjectedClimateGrid {
    pub x0: Vec<i32>,
    pub x1: Vec<i32>,
    pub y0: Vec<i32>,
    pub y1: Vec<i32>,
    pub scale: f64,
    pub offset_x: f64,
    pub offset_y: f64,
}

static PROJ_GRID_CACHE: OnceLock<Mutex<Vec<(GridProjKey, Arc<ProjectedClimateGrid>)>>> =
    OnceLock::new();

pub fn projected_climate_grid(
    grid_h: usize,
    grid_w: usize,
    proj: Projection,
    rot: &[f64; 9],
    width: usize,
    height: usize,
    padding_px: f64,
) -> Arc<ProjectedClimateGrid> {
    let mut rot_bits = [0u64; 9];
    for (i, &v) in rot.iter().enumerate() {
        // Normalize -0.0 to +0.0 just in case.
        let clean = if v == 0.0 { 0.0 } else { v };
        rot_bits[i] = clean.to_bits();
    }
    let proj_id = match proj {
        Projection::Behrmann => 0,
        Projection::Eckert4 => 1,
        Projection::Globe => 2,
    };
    let key = GridProjKey {
        grid_h,
        grid_w,
        proj_id,
        rot_bits,
        width,
        height,
        padding_micros: (padding_px * 1_000_000.0).round() as i64,
    };
    let cache = PROJ_GRID_CACHE.get_or_init(|| Mutex::new(Vec::new()));
    {
        let guard = cache.lock().unwrap();
        for (k, v) in guard.iter() {
            if *k == key {
                return Arc::clone(v);
            }
        }
    }

    let built = Arc::new(build_projected_climate_grid(
        grid_h,
        grid_w,
        proj,
        rot,
        width,
        height,
        padding_px,
    ));
    let mut guard = cache.lock().unwrap();
    for (k, v) in guard.iter() {
        if *k == key {
            return Arc::clone(v);
        }
    }
    if guard.len() >= 8 {
        guard.clear();
    }
    guard.push((key, Arc::clone(&built)));
    built
}

fn build_projected_climate_grid(
    grid_h: usize,
    grid_w: usize,
    proj: Projection,
    rot: &[f64; 9],
    width: usize,
    height: usize,
    padding_px: f64,
) -> ProjectedClimateGrid {
    let (scale, offset_x, offset_y) = fit_view(sphere_extent(proj), width, height, padding_px);
    let n = grid_h * grid_w;
    let half_dlat_deg = 180.0 / grid_h as f64 / 2.0;
    let half_dlon_deg = 360.0 / grid_w as f64 / 2.0;
    let iw = width as i32;
    let ih = height as i32;

    let mut x0 = vec![0i32; n];
    let mut x1 = vec![0i32; n];
    let mut y0 = vec![0i32; n];
    let mut y1 = vec![0i32; n];

    let rows: Vec<(Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>)> = (0..grid_h)
        .into_par_iter()
        .map(|r| {
            let lat_deg = 90.0 - (r as f64 + 0.5) * (180.0 / grid_h as f64);
            let lat_rad = lat_deg.to_radians();
            // For Eckert4, center_theta only depends on rotated latitude, which varies with rotation.
            let mut rx0 = Vec::with_capacity(grid_w);
            let mut rx1 = Vec::with_capacity(grid_w);
            let mut ry0 = Vec::with_capacity(grid_w);
            let mut ry1 = Vec::with_capacity(grid_w);
            for c in 0..grid_w {
                let lon_deg = -180.0 + (c as f64 + 0.5) * (360.0 / grid_w as f64);
                let lon_rad = lon_deg.to_radians();
                let world_p = latlon_to_xyz(lat_rad, lon_rad);
                let rot_p = rotate_point(world_p, rot);
                let (rot_lat, center_lon) = xyz_to_latlon(rot_p);
                let xy = project_latlon(proj, rot_lat, center_lon);
                let px = to_pixel(scale, offset_x, offset_y, xy);

                let mut max_dx = 0.0f64;
                let mut max_dy = 0.0f64;
                for dlat_sign in [1.0f64, -1.0f64] {
                    let clat = (lat_deg + dlat_sign * half_dlat_deg).to_radians();
                    for dlon_sign in [1.0f64, -1.0f64] {
                        let clon = (lon_deg + dlon_sign * half_dlon_deg).to_radians();
                        let cp = rotate_point(latlon_to_xyz(clat, clon), rot);
                        let cxy = project_offset(proj, cp, center_lon);
                        let dx = (cxy[0] - xy[0]).abs();
                        let dy = (cxy[1] - xy[1]).abs();
                        if dx > max_dx {
                            max_dx = dx;
                        }
                        if dy > max_dy {
                            max_dy = dy;
                        }
                    }
                }
                let hw = max_dx * scale * CELL_OVERLAP_FACTOR;
                let hh = max_dy * scale * CELL_OVERLAP_FACTOR;
                rx0.push(((px[0] - hw).round_ties_even() as i32).clamp(0, iw));
                rx1.push(((px[0] + hw).round_ties_even() as i32).clamp(0, iw));
                ry0.push(((px[1] - hh).round_ties_even() as i32).clamp(0, ih));
                ry1.push(((px[1] + hh).round_ties_even() as i32).clamp(0, ih));
            }
            (rx0, rx1, ry0, ry1)
        })
        .collect();

    for (r, (rx0, rx1, ry0, ry1)) in rows.into_iter().enumerate() {
        let base = r * grid_w;
        x0[base..base + grid_w].copy_from_slice(&rx0);
        x1[base..base + grid_w].copy_from_slice(&rx1);
        y0[base..base + grid_w].copy_from_slice(&ry0);
        y1[base..base + grid_w].copy_from_slice(&ry1);
    }

    ProjectedClimateGrid {
        x0,
        x1,
        y0,
        y1,
        scale,
        offset_x,
        offset_y,
    }
}

/// Build the elevation view's local lattice (`iter_local_lattice(eye, spacing_rad)`),
/// preserving exact row and cell iteration order for deterministic last-write-wins fill.
pub fn build_elevation_lattice(spacing_rad: f64) -> Vec<(f64, f64, [f64; 3])> {
    let max_abs_phi = PI * 0.5 - spacing_rad * 0.5;
    let mut phi_values = Vec::new();
    let mut k = 0usize;
    loop {
        let phi = -max_abs_phi + (k as f64) * spacing_rad;
        if phi >= max_abs_phi {
            break;
        }
        phi_values.push(phi);
        k += 1;
    }

    let mut out = Vec::new();
    for phi in phi_values {
        let dtheta = spacing_rad / phi.cos().max(1e-3);
        let n_theta = ((2.0 * PI / dtheta).round_ties_even() as usize).max(1);
        let step = 2.0 * PI / (n_theta as f64);
        for j in 0..n_theta {
            let theta = -PI + (j as f64) * step;
            let pt = latlon_to_xyz(phi, theta);
            out.push((phi, theta, pt));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eckert4_and_behrmann_roundtrip_sanity() {
        let [bx, by] = project_latlon(Projection::Behrmann, 0.5, -1.0);
        assert!((bx + (30.0f64.to_radians()).cos()).abs() < 1e-12);
        assert!((by - 0.5f64.sin() / (30.0f64.to_radians()).cos()).abs() < 1e-12);
        let [ex, ey] = project_latlon(Projection::Eckert4, PI / 2.0, 0.0);
        assert!(ex.abs() < 1e-12);
        assert!(ey > 1.3);
    }
}
