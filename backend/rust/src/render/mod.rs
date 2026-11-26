pub mod biomes;
pub mod kdtree;
pub mod projection;
pub mod raster;

use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};
use std::sync::OnceLock;

use numpy::{
    PyArray1, PyArray2, PyArray3, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2,
    PyReadwriteArrayDyn, PyUntypedArrayMethods,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use rayon::prelude::*;

use crate::render::biomes::{
    bilinear_sample_climate, biome_relative_shade_factor, biome_relative_shade_factor_i64,
    grid_cell_indices_scalar, grid_continentality_dist, sea_ice_annual_balance_scalar,
    seasonal_fraction_above_scalar, seasonal_temp_amplitude, smooth_biome_field_blend_core,
    soft_sea_ice_weight_scalar, BIOME_COLORS, COMBINED_LAKE_ID_CODE, COMBINED_SEA_ID_CODE,
    POLAR_SEA_ICE_FRACTION,
};
use crate::render::kdtree::KdTree;
use crate::render::projection::{
    build_elevation_lattice, fit_view, latlon_to_xyz, project_offset, project_point,
    projected_climate_grid, rotate_point, to_pixel, xyz_to_latlon, BEHRMANN_STANDARD_PARALLEL,
    Projection,
};
use crate::render::raster::{
    crop_even_rgb24, draw_rivers_rgb, elevation_colors_float, elevation_colors_u8, encode_png_rgb,
    encode_png_rgba, fill_rects_generic, fill_rects_rgb, fill_rects_rgb_and_alpha,
    gaussian_blur_rgb, globe_source_size, resample_globe_rgb, resample_globe_rgba, RiverSegments,
    BACKGROUND_RGB, CELL_BLUR_RADIUS_PX, CELL_OVERLAP_FACTOR, GLACIER_COLOR_RGB, LAKE_COLOR_RGB,
    MOUNTAIN_OVERLAY_RGB, OCEAN_PELAGIC_RELIEF_BLEND, PADDING_PX, PLAINS_PLATEAU_OVERLAY_RGB,
    REFERENCE_WIDTH_PX, RELIEF_BLEND_MAX, RELIEF_ELEVATION_RANGE_M, SEA_COLOR_RGB,
    TERRAIN_MOUNTAIN, TERRAIN_OVERLAY_ALPHA, TERRAIN_PLAINS_PLATEAU,
};

const GLACIER_VISIBLE_DEPTH_M: f64 = 10.0;
const LAKE_MIN_VISIBLE_DEPTH_M: f64 = 1.0;
const PLANET_RADIUS_M: f64 = 6371.0 * 1000.0;
const CHANNEL_INCISION_MIN_DEPTH_M: f64 = 80.0;
const CHANNEL_INCISION_MIN_WIDTH_M: f64 = 250.0;
const CHANNEL_INCISION_WIDTH_RANGE_M: f64 = 400.0;
const HILLSHADE_NEIGHBOR_COUNT: usize = 8;
const HILLSHADE_SUN_AZIMUTH_DEG: f64 = 315.0;
const HILLSHADE_SUN_ALTITUDE_DEG: f64 = 45.0;
const HILLSHADE_VERTICAL_EXAGGERATION: f64 = 30.0;
const HILLSHADE_STRENGTH: f64 = 0.85;
const HILLSHADE_MIN_MULT: f64 = 0.45;
const HILLSHADE_MAX_MULT: f64 = 1.5;

static RENDER_POOL: OnceLock<(usize, rayon::ThreadPool)> = OnceLock::new();

fn render_pool() -> &'static (usize, rayon::ThreadPool) {
    RENDER_POOL.get_or_init(|| {
        let threads = std::env::var("MANTLE_BLOOM_RENDER_THREADS")
            .ok()
            .and_then(|s| s.parse::<usize>().ok())
            .map(|n| n.clamp(1, 16))
            .unwrap_or_else(|| {
                let avail = std::thread::available_parallelism()
                    .map(|v| v.get())
                    .unwrap_or(2);
                if avail <= 4 { 1 } else { 2 }
            });
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .thread_name(|i| format!("mb-render-{i}"))
            .build()
            .expect("failed to build dedicated render ThreadPool");
        (threads, pool)
    })
}

#[pyfunction]
pub fn render_pool_threads() -> usize {
    render_pool().0
}

#[pyfunction]
pub fn grid_edt(py: Python<'_>, mask: PyReadonlyArray2<'_, bool>) -> PyResult<Py<PyArray2<f64>>> {
    let view = mask.as_array();
    let (h, w) = (view.nrows(), view.ncols());
    let flat: Vec<bool> = view.iter().copied().collect();
    let dist = py.detach(|| render_pool().1.install(|| biomes::edt_2d_bool(&flat, h, w)));
    let arr = PyArray2::from_vec2(
        py,
        &dist.chunks(w).map(|c| c.to_vec()).collect::<Vec<_>>(),
    )?;
    Ok(arr.unbind())
}

#[pyfunction]
#[pyo3(signature = (
    temperature_c,
    precipitation_mm,
    is_ocean,
    lat_deg,
    axial_tilt_deg,
    has_sea_ice=None,
    sea_ice_weight=None,
))]
pub fn smooth_biome_field_blend_rs(
    py: Python<'_>,
    temperature_c: PyReadonlyArray2<'_, f64>,
    precipitation_mm: PyReadonlyArray2<'_, f64>,
    is_ocean: PyReadonlyArray2<'_, bool>,
    lat_deg: PyReadonlyArray1<'_, f64>,
    axial_tilt_deg: f64,
    has_sea_ice: Option<PyReadonlyArray2<'_, bool>>,
    sea_ice_weight: Option<PyReadonlyArray2<'_, f64>>,
) -> PyResult<(Py<PyArray2<i64>>, Py<PyArray2<f64>>, Py<PyArray2<i64>>)> {
    let t_arr = temperature_c.as_array();
    let (h, w) = (t_arr.nrows(), t_arr.ncols());
    let t_vec: Vec<f64> = t_arr.iter().copied().collect();
    let p_vec: Vec<f64> = precipitation_mm.as_array().iter().copied().collect();
    let o_vec: Vec<bool> = is_ocean.as_array().iter().copied().collect();
    let lat_vec: Vec<f64> = lat_deg.as_array().iter().copied().collect();
    let ice_vec: Option<Vec<bool>> = has_sea_ice.map(|a| a.as_array().iter().copied().collect());
    let sw_vec: Option<Vec<f64>> = sea_ice_weight.map(|a| a.as_array().iter().copied().collect());

    if lat_vec.len() != h || p_vec.len() != h * w || o_vec.len() != h * w {
        return Err(PyValueError::new_err("shape mismatch in smooth_biome_field_blend_rs"));
    }
    if let Some(ref v) = ice_vec {
        if v.len() != h * w {
            return Err(PyValueError::new_err(
                "has_sea_ice shape mismatch in smooth_biome_field_blend_rs",
            ));
        }
    }
    if let Some(ref v) = sw_vec {
        if v.len() != h * w {
            return Err(PyValueError::new_err(
                "sea_ice_weight shape mismatch in smooth_biome_field_blend_rs",
            ));
        }
    }

    let (ids, shares, runner_ups) = py.detach(|| {
        render_pool().1.install(|| {
            smooth_biome_field_blend_core(
                &t_vec,
                &p_vec,
                &o_vec,
                &lat_vec,
                h,
                w,
                axial_tilt_deg,
                ice_vec.as_deref(),
                sw_vec.as_deref(),
            )
        })
    });

    let py_ids = PyArray2::from_vec2(
        py,
        &ids.chunks(w)
            .map(|c| c.iter().map(|&v| v as i64).collect())
            .collect::<Vec<Vec<i64>>>(),
    )?;
    let py_shares = PyArray2::from_vec2(
        py,
        &shares.chunks(w).map(|c| c.to_vec()).collect::<Vec<_>>(),
    )?;
    let py_ru = PyArray2::from_vec2(
        py,
        &runner_ups
            .chunks(w)
            .map(|c| c.iter().map(|&v| v as i64).collect())
            .collect::<Vec<Vec<i64>>>(),
    )?;
    Ok((py_ids.unbind(), py_shares.unbind(), py_ru.unbind()))
}

#[pyfunction]
pub fn biome_relative_shade_factor_rs(
    py: Python<'_>,
    biome_ids: PyReadonlyArray1<'_, i64>,
    elevation_m: PyReadonlyArray1<'_, f64>,
) -> PyResult<Py<PyArray1<f64>>> {
    let ids = biome_ids.as_slice()?.to_vec();
    let elev = elevation_m.as_slice()?.to_vec();
    if ids.len() != elev.len() {
        return Err(PyValueError::new_err("biome_ids and elevation_m length mismatch"));
    }
    let out = py.detach(|| {
        render_pool()
            .1
            .install(|| biome_relative_shade_factor_i64(&ids, &elev))
    });
    Ok(PyArray1::from_vec(py, out).unbind())
}

#[pyfunction]
pub fn fill_rects_rs(
    py: Python<'_>,
    mut pixels: PyReadwriteArrayDyn<'_, u8>,
    x0: PyReadonlyArray1<'_, i64>,
    x1: PyReadonlyArray1<'_, i64>,
    y0: PyReadonlyArray1<'_, i64>,
    y1: PyReadonlyArray1<'_, i64>,
    colors: PyReadonlyArray2<'_, u8>,
) -> PyResult<()> {
    let shape = pixels.shape().to_vec();
    if shape.len() != 3 {
        return Err(PyValueError::new_err("pixels must be 3D (H, W, C)"));
    }
    let width = shape[1];
    let channels = shape[2];
    let x0s = x0.as_slice()?;
    let x1s = x1.as_slice()?;
    let y0s = y0.as_slice()?;
    let y1s = y1.as_slice()?;
    let cols = colors.as_slice()?;
    if cols.len() != x0s.len() * channels {
        return Err(PyValueError::new_err("colors shape mismatch"));
    }
    let pix_slice = pixels.as_slice_mut()?;
    // Keep serial pointer write inside detach without aliasing.
    let x0_v = x0s.to_vec();
    let x1_v = x1s.to_vec();
    let y0_v = y0s.to_vec();
    let y1_v = y1s.to_vec();
    let cols_v = cols.to_vec();
    py.detach(|| {
        fill_rects_generic(
            pix_slice,
            width,
            channels,
            &x0_v,
            &x1_v,
            &y0_v,
            &y1_v,
            &cols_v,
        );
    });
    Ok(())
}

#[derive(Clone, Copy)]
struct HeapItem {
    pass_m: f64,
    node: usize,
}

impl PartialEq for HeapItem {
    fn eq(&self, other: &Self) -> bool {
        self.pass_m == other.pass_m && self.node == other.node
    }
}
impl Eq for HeapItem {}
impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .pass_m
            .partial_cmp(&self.pass_m)
            .unwrap_or(Ordering::Equal)
            .then_with(|| other.node.cmp(&self.node))
    }
}
impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[pyfunction]
pub fn spill_path_rs(
    py: Python<'_>,
    indptr: PyReadonlyArray1<'_, i64>,
    indices: PyReadonlyArray1<'_, i64>,
    elevation: PyReadonlyArray1<'_, f64>,
    is_ocean: PyReadonlyArray1<'_, bool>,
    source: usize,
    target: usize,
    ceiling: f64,
    max_nodes: usize,
) -> PyResult<Vec<i64>> {
    let indptr_s = indptr.as_slice()?;
    let indices_s = indices.as_slice()?;
    let elev_s = elevation.as_slice()?;
    let ocean_s = is_ocean.as_slice()?;

    let path = py.detach(|| {
        let n = elev_s.len();
        let mut best = vec![f64::INFINITY; n];
        best[source] = elev_s[source];
        let mut heap = BinaryHeap::new();
        heap.push(HeapItem {
            pass_m: elev_s[source],
            node: source,
        });
        let mut level = f64::INFINITY;
        let mut found = false;
        let mut explored = 0usize;
        while let Some(HeapItem { pass_m, node }) = heap.pop() {
            if node == target {
                level = pass_m;
                found = true;
                break;
            }
            if pass_m > best[node] {
                continue;
            }
            explored += 1;
            if pass_m > ceiling || explored > max_nodes {
                break;
            }
            for j in indptr_s[node]..indptr_s[node + 1] {
                let nxt = indices_s[j as usize] as usize;
                if ocean_s[nxt] && nxt != target {
                    continue;
                }
                let cand = pass_m.max(elev_s[nxt]);
                if cand < best[nxt] {
                    best[nxt] = cand;
                    heap.push(HeapItem {
                        pass_m: cand,
                        node: nxt,
                    });
                }
            }
        }
        if !found || !level.is_finite() {
            return Vec::new();
        }
        let mut came_from = vec![-2isize; n];
        came_from[source] = -1;
        let mut queue = VecDeque::new();
        queue.push_back(source);
        while let Some(node) = queue.pop_front() {
            if came_from[target] != -2 {
                break;
            }
            for j in indptr_s[node]..indptr_s[node + 1] {
                let nxt = indices_s[j as usize] as usize;
                if came_from[nxt] == -2
                    && (nxt == target || !ocean_s[nxt])
                    && elev_s[nxt] <= level
                {
                    came_from[nxt] = node as isize;
                    queue.push_back(nxt);
                }
            }
        }
        let mut out = Vec::new();
        let mut cur = target as isize;
        while cur >= 0 {
            out.push(cur as i64);
            cur = came_from[cur as usize];
        }
        out.reverse();
        out
    });
    Ok(path)
}

fn to_pts3(arr: &PyReadonlyArray2<'_, f64>) -> PyResult<Vec<[f64; 3]>> {
    let v = arr.as_array();
    if v.ncols() != 3 {
        return Err(PyValueError::new_err("expected (N, 3) point array"));
    }
    Ok(v.rows().into_iter().map(|r| [r[0], r[1], r[2]]).collect())
}

fn to_rot9(arr: &PyReadonlyArray2<'_, f64>) -> PyResult<[f64; 9]> {
    let v = arr.as_array();
    if v.shape() != [3, 3] {
        return Err(PyValueError::new_err("view_rotation must be (3, 3)"));
    }
    Ok([
        v[[0, 0]],
        v[[0, 1]],
        v[[0, 2]],
        v[[1, 0]],
        v[[1, 1]],
        v[[1, 2]],
        v[[2, 0]],
        v[[2, 1]],
        v[[2, 2]],
    ])
}

fn compute_hillshade_from_tree(
    points: &[[f64; 3]],
    elevation: &[f64],
    glacier_depth: &[f64],
    channel_depth: &[f64],
    channel_width: &[f64],
    tree: &KdTree,
) -> Vec<f64> {
    let n = points.len();
    if n <= HILLSHADE_NEIGHBOR_COUNT {
        return vec![1.0; n];
    }
    let az = HILLSHADE_SUN_AZIMUTH_DEG.to_radians();
    let alt = HILLSHADE_SUN_ALTITUDE_DEG.to_radians();
    let sun = [alt.cos() * az.sin(), alt.cos() * az.cos(), alt.sin()];
    let flat_illum = alt.sin();

    (0..n)
        .into_par_iter()
        .map(|i| {
            let p = points[i];
            let knn = tree.query_knn(p, HILLSHADE_NEIGHBOR_COUNT + 1);
            let nbs: Vec<usize> = knn.into_iter().filter(|&idx| idx != i).take(HILLSHADE_NEIGHBOR_COUNT).collect();
            if nbs.len() < HILLSHADE_NEIGHBOR_COUNT {
                return 1.0;
            }
            let up = [0.0, 0.0, 1.0];
            let dot_up = p[2];
            let mut north = [up[0] - dot_up * p[0], up[1] - dot_up * p[1], up[2] - dot_up * p[2]];
            let mut nnorm = (north[0] * north[0] + north[1] * north[1] + north[2] * north[2]).sqrt();
            if nnorm < 1e-6 {
                let dot_x = p[0];
                north = [1.0 - dot_x * p[0], -dot_x * p[1], -dot_x * p[2]];
                nnorm = (north[0] * north[0] + north[1] * north[1] + north[2] * north[2]).sqrt();
            }
            north = [north[0] / nnorm, north[1] / nnorm, north[2] / nnorm];
            let mut east = [
                north[1] * p[2] - north[2] * p[1],
                north[2] * p[0] - north[0] * p[2],
                north[0] * p[1] - north[1] * p[0],
            ];
            let enorm = (east[0] * east[0] + east[1] * east[1] + east[2] * east[2]).sqrt();
            east = [east[0] / enorm, east[1] / enorm, east[2] / enorm];

            let surf_i = if glacier_depth[i] > 0.0 {
                elevation[i] + glacier_depth[i]
            } else {
                let d_ex = (channel_depth[i] - CHANNEL_INCISION_MIN_DEPTH_M).max(0.0);
                let w_vis = ((channel_width[i] - CHANNEL_INCISION_MIN_WIDTH_M)
                    / CHANNEL_INCISION_WIDTH_RANGE_M)
                    .clamp(0.0, 1.0);
                elevation[i] - d_ex * w_vis
            };

            let (mut sxx, mut syy, mut sxy, mut sxz, mut syz) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for &nb in &nbs {
                let np = points[nb];
                let off = [np[0] - p[0], np[1] - p[1], np[2] - p[2]];
                let de = (off[0] * east[0] + off[1] * east[1] + off[2] * east[2]) * PLANET_RADIUS_M;
                let dn = (off[0] * north[0] + off[1] * north[1] + off[2] * north[2]) * PLANET_RADIUS_M;
                let surf_nb = if glacier_depth[nb] > 0.0 {
                    elevation[nb] + glacier_depth[nb]
                } else {
                    let d_ex = (channel_depth[nb] - CHANNEL_INCISION_MIN_DEPTH_M).max(0.0);
                    let w_vis = ((channel_width[nb] - CHANNEL_INCISION_MIN_WIDTH_M)
                        / CHANNEL_INCISION_WIDTH_RANGE_M)
                        .clamp(0.0, 1.0);
                    elevation[nb] - d_ex * w_vis
                };
                let dz = surf_nb - surf_i;
                sxx += de * de;
                syy += dn * dn;
                sxy += de * dn;
                sxz += de * dz;
                syz += dn * dz;
            }
            let det = sxx * syy - sxy * sxy;
            let (gx, gy) = if det.abs() < 1e-6 {
                (0.0, 0.0)
            } else {
                (
                    ((syy * sxz - sxy * syz) / det) * HILLSHADE_VERTICAL_EXAGGERATION,
                    ((sxx * syz - sxy * sxz) / det) * HILLSHADE_VERTICAL_EXAGGERATION,
                )
            };
            let nlen = (gx * gx + gy * gy + 1.0).sqrt();
            let normal = [-gx / nlen, -gy / nlen, 1.0 / nlen];
            let illum = (normal[0] * sun[0] + normal[1] * sun[1] + normal[2] * sun[2]).max(0.0);
            let ratio = illum / flat_illum;
            (1.0 + HILLSHADE_STRENGTH * (ratio - 1.0)).clamp(HILLSHADE_MIN_MULT, HILLSHADE_MAX_MULT)
        })
        .collect()
}

fn package_outputs(
    py: Python<'_>,
    rgb: &[u8],
    alpha: Option<&[u8]>,
    width: usize,
    height: usize,
    want_png: bool,
    want_rgb: bool,
) -> PyResult<(Option<Py<PyBytes>>, Option<Py<PyArray3<u8>>>)> {
    let png_bytes = if want_png {
        let data = py
            .detach(|| {
                render_pool().1.install(|| match alpha {
                    Some(a) => encode_png_rgba(rgb, a, width, height),
                    None => encode_png_rgb(rgb, width, height),
                })
            })
            .map_err(PyValueError::new_err)?;
        Some(PyBytes::new(py, &data).unbind())
    } else {
        None
    };

    let rgb_arr = if want_rgb {
        let (cropped, enc_w, enc_h) = crop_even_rgb24(rgb, width, height);
        let arr = PyArray3::<u8>::zeros(py, [enc_h, enc_w, 3], false);
        unsafe {
            arr.as_slice_mut()?.copy_from_slice(&cropped);
        }
        Some(arr.unbind())
    } else {
        None
    };

    Ok((png_bytes, rgb_arr))
}

#[pyfunction]
#[pyo3(signature = (
    view,
    projection,
    width,
    height,
    view_rotation,
    grid_h,
    grid_w,
    sea_level_m,
    axial_tilt_deg,
    solar_multiplier,
    smooth_sea_ice,
    soft_sea_ice_edge,
    points,
    elevation,
    lake_depth,
    lake_ice,
    glacier_depth,
    hillshade,
    channel_depth,
    channel_width,
    hydro_points_same,
    hydro_points,
    hydro_is_ocean,
    hydro_is_sea,
    hydro_perennial_lake_ice,
    clim_air_temp,
    clim_ocean_temp,
    clim_precip,
    clim_lat_deg,
    clim_elev_m,
    clim_is_ocean,
    clim_h,
    clim_w,
    river_from_xyz,
    river_to_xyz,
    river_width_frac,
    river_color_alpha,
    river_canyon_fill,
    want_png,
    want_rgb,
    stroke_scale_override=1.0,
))]
#[allow(clippy::too_many_arguments)]
pub fn render_biome_or_combined(
    py: Python<'_>,
    view: &str,
    projection: &str,
    width: usize,
    height: usize,
    view_rotation: PyReadonlyArray2<'_, f64>,
    grid_h: usize,
    grid_w: usize,
    sea_level_m: f64,
    axial_tilt_deg: f64,
    solar_multiplier: f64,
    smooth_sea_ice: bool,
    soft_sea_ice_edge: bool,
    points: PyReadonlyArray2<'_, f64>,
    elevation: PyReadonlyArray1<'_, f64>,
    lake_depth: PyReadonlyArray1<'_, f64>,
    lake_ice: PyReadonlyArray1<'_, f64>,
    glacier_depth: PyReadonlyArray1<'_, f64>,
    hillshade: Option<PyReadonlyArray1<'_, f64>>,
    channel_depth: Option<PyReadonlyArray1<'_, f64>>,
    channel_width: Option<PyReadonlyArray1<'_, f64>>,
    hydro_points_same: bool,
    hydro_points: Option<PyReadonlyArray2<'_, f64>>,
    hydro_is_ocean: Option<PyReadonlyArray1<'_, bool>>,
    hydro_is_sea: Option<PyReadonlyArray1<'_, bool>>,
    hydro_perennial_lake_ice: Option<PyReadonlyArray1<'_, bool>>,
    clim_air_temp: PyReadonlyArray1<'_, f64>,
    clim_ocean_temp: PyReadonlyArray1<'_, f64>,
    clim_precip: PyReadonlyArray1<'_, f64>,
    clim_lat_deg: PyReadonlyArray1<'_, f64>,
    clim_elev_m: Option<PyReadonlyArray1<'_, f64>>,
    clim_is_ocean: Option<PyReadonlyArray1<'_, bool>>,
    clim_h: usize,
    clim_w: usize,
    river_from_xyz: PyReadonlyArray2<'_, f64>,
    river_to_xyz: PyReadonlyArray2<'_, f64>,
    river_width_frac: PyReadonlyArray1<'_, f64>,
    river_color_alpha: PyReadonlyArray1<'_, f64>,
    river_canyon_fill: PyReadonlyArray1<'_, i32>,
    want_png: bool,
    want_rgb: bool,
    stroke_scale_override: f64,
) -> PyResult<(Option<Py<PyBytes>>, Option<Py<PyArray3<u8>>>)> {
    let is_combined = match view {
        "combined" => true,
        "biome" => false,
        _ => return Err(PyValueError::new_err(format!("Unsupported view {view:?}"))),
    };
    let req_proj = Projection::parse(projection).map_err(PyValueError::new_err)?;
    let rot = to_rot9(&view_rotation)?;
    let pts_v = to_pts3(&points)?;
    let elev_v = elevation.as_slice()?.to_vec();
    let lake_v = lake_depth.as_slice()?.to_vec();
    let lake_ice_v = lake_ice.as_slice()?.to_vec();
    let glac_v = glacier_depth.as_slice()?.to_vec();
    let hill_v = hillshade.map(|a| a.as_slice().unwrap().to_vec());
    let cdepth_v = channel_depth.map(|a| a.as_slice().unwrap().to_vec());
    let cwidth_v = channel_width.map(|a| a.as_slice().unwrap().to_vec());
    let hpts_v = hydro_points.map(|a| to_pts3(&a)).transpose()?;
    let hocean_v = hydro_is_ocean.map(|a| a.as_slice().unwrap().to_vec());
    let hsea_v = hydro_is_sea.map(|a| a.as_slice().unwrap().to_vec());
    let hlake_ice_v = hydro_perennial_lake_ice.map(|a| a.as_slice().unwrap().to_vec());

    let c_air = clim_air_temp.as_slice()?.to_vec();
    let c_sst = clim_ocean_temp.as_slice()?.to_vec();
    let c_pre = clim_precip.as_slice()?.to_vec();
    let c_lat = clim_lat_deg.as_slice()?.to_vec();
    let c_elev = clim_elev_m.map(|a| a.as_slice().unwrap().to_vec());
    let c_ocean = clim_is_ocean.map(|a| a.as_slice().unwrap().to_vec());

    let r_from = to_pts3(&river_from_xyz)?;
    let r_to = to_pts3(&river_to_xyz)?;
    let r_wf = river_width_frac.as_slice()?.to_vec();
    let r_ca = river_color_alpha.as_slice()?.to_vec();
    let r_cf = river_canyon_fill.as_slice()?.to_vec();

    let (out_rgb, out_alpha) = py.detach(|| {
        render_pool().1.install(|| {
            let (render_proj, map_w, map_h, stroke_mult) = if req_proj == Projection::Globe {
                let (mw, mh) = globe_source_size(width, height);
                let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
                (
                    Projection::Behrmann,
                    mw,
                    mh,
                    width as f64 / (mw as f64 * cos_phi0),
                )
            } else {
                (req_proj, width, height, stroke_scale_override)
            };

            let n_cells = grid_h * grid_w;
            let tree = KdTree::build(&pts_v);
            let hydro_tree = if !hydro_points_same {
                hpts_v.as_ref().map(|hp| KdTree::build(hp))
            } else {
                None
            };
            let resolved_hill = if is_combined {
                match &hill_v {
                    Some(h) => h.clone(),
                    None => compute_hillshade_from_tree(
                        &pts_v,
                        &elev_v,
                        &glac_v,
                        cdepth_v.as_deref().unwrap_or(&[]),
                        cwidth_v.as_deref().unwrap_or(&[]),
                        &tree,
                    ),
                }
            } else {
                Vec::new()
            };

            let any_lake_cap = lake_ice_v.iter().any(|&v| v > 0.0);
            let clim_cont = if any_lake_cap && hlake_ice_v.is_none() {
                c_ocean
                    .as_ref()
                    .map(|co| grid_continentality_dist(co, clim_h, clim_w))
            } else {
                None
            };

            // Precompute coarse sea-ice grid if smooth_sea_ice == false.
            let coarse_sea_ice: Option<Vec<(f64, bool, f64, f64)>> = if !smooth_sea_ice {
                let mut g = Vec::with_capacity(clim_h * clim_w);
                for r in 0..clim_h {
                    let clat = c_lat[r];
                    for c in 0..clim_w {
                        let res = sea_ice_annual_balance_scalar(
                            c_sst[r * clim_w + c],
                            clat,
                            axial_tilt_deg,
                            solar_multiplier,
                            true,
                        );
                        g.push((
                            res.temperature_c,
                            res.perennial,
                            res.fraction,
                            res.perennial_cover,
                        ));
                    }
                }
                Some(g)
            } else {
                None
            };

            let lat_deg_rows: Vec<f64> = (0..grid_h)
                .map(|r| 90.0 - (r as f64 + 0.5) * (180.0 / grid_h as f64))
                .collect();

            struct SampledCell {
                bed_m: f64,
                is_ocean: bool,
                is_sea: bool,
                lake_depth_m: f64,
                glacier_depth_m: f64,
                hillshade: f64,
                display_temp_c: f64,
                precip_mm: f64,
                lake_frozen: bool,
                sea_perennial: bool,
                sea_fraction: f64,
                sea_perennial_cover: f64,
                glacier_on: bool,
                floating_ice: bool,
                climate_ocean: bool,
                has_sea_ice: bool,
                sea_ice_weight: f64,
            }

            let sampled_rows: Vec<Vec<SampledCell>> = (0..grid_h)
                .into_par_iter()
                .map(|r| {
                    let lat_d = lat_deg_rows[r];
                    let lat_r = lat_d.to_radians();
                    let sin_lat = lat_r.sin();
                    let cos_lat = lat_r.cos();
                    let mut row = Vec::with_capacity(grid_w);
                    for c in 0..grid_w {
                        let lon_d = -180.0 + (c as f64 + 0.5) * (360.0 / grid_w as f64);
                        let lon_r = lon_d.to_radians();
                        let xyz = [cos_lat * lon_r.cos(), cos_lat * lon_r.sin(), sin_lat];
                        let idx = tree.query_nearest(xyz);
                        let bed_m = elev_v[idx];
                        let l_depth = lake_v[idx];
                        let l_cap = lake_ice_v[idx];
                        let g_depth = glac_v[idx];
                        let hs = if is_combined { resolved_hill[idx] } else { 1.0 };

                        let hidx = if hydro_points_same {
                            idx
                        } else if let Some(ht) = hydro_tree.as_ref() {
                            ht.query_nearest(xyz)
                        } else {
                            idx
                        };
                        let is_ocean = hocean_v
                            .as_ref()
                            .map_or(bed_m <= sea_level_m, |ho| ho[hidx]);
                        let is_sea = hsea_v.as_ref().map_or(false, |hs_arr| hs_arr[hidx]);

                        let air_t = bilinear_sample_climate(&c_air, clim_h, clim_w, lat_d, lon_d);
                        let ocean_t = bilinear_sample_climate(&c_sst, clim_h, clim_w, lat_d, lon_d);
                        let precip_mm =
                            bilinear_sample_climate(&c_pre, clim_h, clim_w, lat_d, lon_d);

                        let (crow, ccol) = grid_cell_indices_scalar(lat_d, lon_d, clim_h, clim_w);
                        let cflat = crow * clim_w + ccol;

                        let lake_frozen_raw = if l_cap <= 0.0 {
                            false
                        } else if let Some(pl) = hlake_ice_v.as_ref() {
                            pl[hidx]
                        } else {
                            let water_elev = bed_m + l_depth + l_cap;
                            let base_elev = c_elev.as_ref().map_or(0.0, |ce| ce[cflat]);
                            let temp = c_air[cflat] - 6.5 * (water_elev - base_elev) / 1000.0;
                            let cont = clim_cont.as_ref().map_or(0.0, |cc| cc[cflat]);
                            let amp = seasonal_temp_amplitude(
                                c_lat[crow],
                                cont,
                                axial_tilt_deg,
                                Some(bed_m - sea_level_m),
                            );
                            seasonal_fraction_above_scalar(temp, amp, 0.0) == 0.0
                        };
                        let lake_frozen = lake_frozen_raw && !is_ocean;

                        let ocean_water = is_ocean && g_depth == 0.0;
                        let sea_cell = if smooth_sea_ice {
                            sea_ice_annual_balance_scalar(
                                ocean_t,
                                lat_d,
                                axial_tilt_deg,
                                solar_multiplier,
                                ocean_water,
                            )
                        } else {
                            let cg = coarse_sea_ice.as_ref().unwrap()[cflat];
                            biomes::SeaIceCell {
                                temperature_c: cg.0,
                                perennial: ocean_water && cg.1,
                                fraction: if ocean_water { cg.2 } else { 0.0 },
                                perennial_cover: if ocean_water { cg.3 } else { 0.0 },
                            }
                        };

                        let glacier_air = air_t - 6.5 * g_depth / 1000.0;
                        let open_surf = if is_ocean && g_depth == 0.0 {
                            ocean_t
                        } else {
                            glacier_air
                        };
                        let display_temp_c = if sea_cell.perennial {
                            sea_cell.temperature_c
                        } else {
                            open_surf
                        };
                        let glacier_on = g_depth >= GLACIER_VISIBLE_DEPTH_M;
                        let floating_ice = (lake_frozen || sea_cell.perennial) && !glacier_on;
                        let has_sea_ice = sea_cell.fraction >= POLAR_SEA_ICE_FRACTION;
                        let climate_ocean = is_ocean && g_depth == 0.0;
                        let sea_ice_weight = if soft_sea_ice_edge {
                            soft_sea_ice_weight_scalar(sea_cell.fraction, POLAR_SEA_ICE_FRACTION)
                        } else if has_sea_ice {
                            1.0
                        } else {
                            0.0
                        };

                        row.push(SampledCell {
                            bed_m,
                            is_ocean,
                            is_sea,
                            lake_depth_m: l_depth,
                            glacier_depth_m: g_depth,
                            hillshade: hs,
                            display_temp_c,
                            precip_mm,
                            lake_frozen,
                            sea_perennial: sea_cell.perennial,
                            sea_fraction: sea_cell.fraction,
                            sea_perennial_cover: sea_cell.perennial_cover,
                            glacier_on,
                            floating_ice,
                            climate_ocean,
                            has_sea_ice,
                            sea_ice_weight,
                        });
                    }
                    row
                })
                .collect();

            let mut disp_t = Vec::with_capacity(n_cells);
            let mut precip_all = Vec::with_capacity(n_cells);
            let mut clim_ocean_all = Vec::with_capacity(n_cells);
            let mut has_ice_all = Vec::with_capacity(n_cells);
            let mut ice_w_all = Vec::with_capacity(n_cells);
            for row in &sampled_rows {
                for cell in row {
                    disp_t.push(cell.display_temp_c);
                    precip_all.push(cell.precip_mm);
                    clim_ocean_all.push(cell.climate_ocean);
                    has_ice_all.push(cell.has_sea_ice);
                    ice_w_all.push(cell.sea_ice_weight);
                }
            }

            let (winner_ids, mut winner_share, runner_up_ids) = smooth_biome_field_blend_core(
                &disp_t,
                &precip_all,
                &clim_ocean_all,
                &lat_deg_rows,
                grid_h,
                grid_w,
                axial_tilt_deg,
                Some(&has_ice_all),
                if soft_sea_ice_edge {
                    Some(&ice_w_all)
                } else {
                    None
                },
            );

            let mut display_ids = vec![0u8; n_cells];
            let mut flat_cells = Vec::with_capacity(n_cells);
            for row in sampled_rows {
                flat_cells.extend(row);
            }

            for i in 0..n_cells {
                let c = &flat_cells[i];
                let mut id = if c.sea_fraction >= POLAR_SEA_ICE_FRACTION {
                    biomes::POLAR_SEA_ICE
                } else {
                    winner_ids[i]
                };
                if c.glacier_on || c.lake_frozen {
                    id = biomes::ICE_CAP;
                }
                display_ids[i] = id;
                let flat_ice = if soft_sea_ice_edge {
                    c.glacier_on || c.lake_frozen
                } else {
                    c.glacier_on || c.lake_frozen || c.has_sea_ice
                };
                if flat_ice {
                    winner_share[i] = 1.0;
                }
            }

            let (colors, cell_alpha): (Vec<[u8; 3]>, Vec<u8>) = if !is_combined {
                let cols: Vec<[u8; 3]> = (0..n_cells)
                    .into_par_iter()
                    .map(|i| {
                        let w_c = BIOME_COLORS[display_ids[i] as usize];
                        let r_c = BIOME_COLORS[runner_up_ids[i] as usize];
                        let sh = winner_share[i];
                        let inv = 1.0 - sh;
                        [
                            (w_c[0] as f64 * sh + r_c[0] as f64 * inv)
                                .round_ties_even()
                                .clamp(0.0, 255.0) as u8,
                            (w_c[1] as f64 * sh + r_c[1] as f64 * inv)
                                .round_ties_even()
                                .clamp(0.0, 255.0) as u8,
                            (w_c[2] as f64 * sh + r_c[2] as f64 * inv)
                                .round_ties_even()
                                .clamp(0.0, 255.0) as u8,
                        ]
                    })
                    .collect();
                let alp: Vec<u8> = display_ids.iter().map(|&id| 255u8 - (id + 1)).collect();
                (cols, alp)
            } else {
                let surface_elev: Vec<f64> = flat_cells
                    .iter()
                    .map(|c| c.bed_m + c.glacier_depth_m)
                    .collect();
                let shade = biome_relative_shade_factor(&display_ids, &surface_elev);
                let mut cols = Vec::with_capacity(n_cells);
                let mut alp = Vec::with_capacity(n_cells);
                let results: Vec<([u8; 3], u8)> = (0..n_cells)
                    .into_par_iter()
                    .map(|i| {
                        let c = &flat_cells[i];
                        let w_c = BIOME_COLORS[display_ids[i] as usize];
                        let r_c = BIOME_COLORS[runner_up_ids[i] as usize];
                        let sh = winner_share[i];
                        let inv = 1.0 - sh;
                        let b_rgb = [
                            w_c[0] as f64 * sh + r_c[0] as f64 * inv,
                            w_c[1] as f64 * sh + r_c[1] as f64 * inv,
                            w_c[2] as f64 * sh + r_c[2] as f64 * inv,
                        ];
                        let elev_m = surface_elev[i];
                        let t_rgb = elevation_colors_float(elev_m, sea_level_m);
                        let mut col_f = if c.is_ocean {
                            [
                                b_rgb[0] * (1.0 - OCEAN_PELAGIC_RELIEF_BLEND)
                                    + t_rgb[0] * OCEAN_PELAGIC_RELIEF_BLEND,
                                b_rgb[1] * (1.0 - OCEAN_PELAGIC_RELIEF_BLEND)
                                    + t_rgb[1] * OCEAN_PELAGIC_RELIEF_BLEND,
                                b_rgb[2] * (1.0 - OCEAN_PELAGIC_RELIEF_BLEND)
                                    + t_rgb[2] * OCEAN_PELAGIC_RELIEF_BLEND,
                            ]
                        } else {
                            let sf = shade[i];
                            let shaded = [
                                (b_rgb[0] * sf).clamp(0.0, 255.0),
                                (b_rgb[1] * sf).clamp(0.0, 255.0),
                                (b_rgb[2] * sf).clamp(0.0, 255.0),
                            ];
                            let rt = ((elev_m - sea_level_m) / RELIEF_ELEVATION_RANGE_M)
                                .clamp(0.0, 1.0);
                            let blend = rt * RELIEF_BLEND_MAX;
                            [
                                ((shaded[0] * (1.0 - blend) + t_rgb[0] * blend) * c.hillshade)
                                    .clamp(0.0, 255.0),
                                ((shaded[1] * (1.0 - blend) + t_rgb[1] * blend) * c.hillshade)
                                    .clamp(0.0, 255.0),
                                ((shaded[2] * (1.0 - blend) + t_rgb[2] * blend) * c.hillshade)
                                    .clamp(0.0, 255.0),
                            ]
                        };
                        let is_lake = !c.is_ocean && c.lake_depth_m > LAKE_MIN_VISIBLE_DEPTH_M;
                        let is_sea_cell = is_lake && c.is_sea;
                        if is_lake {
                            col_f = if is_sea_cell {
                                SEA_COLOR_RGB
                            } else {
                                LAKE_COLOR_RGB
                            };
                        }
                        let sea_weight = if c.glacier_on || c.lake_frozen {
                            0.0
                        } else if soft_sea_ice_edge {
                            c.sea_perennial_cover
                        } else if c.sea_perennial {
                            1.0
                        } else {
                            0.0
                        };
                        if sea_weight > 0.0 {
                            col_f = [
                                col_f[0] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[0] * sea_weight,
                                col_f[1] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[1] * sea_weight,
                                col_f[2] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[2] * sea_weight,
                            ];
                        }
                        if c.glacier_on || c.lake_frozen {
                            let s = if c.glacier_on { c.hillshade } else { 1.0 };
                            col_f = [
                                (GLACIER_COLOR_RGB[0] * s).clamp(0.0, 255.0),
                                (GLACIER_COLOR_RGB[1] * s).clamp(0.0, 255.0),
                                (GLACIER_COLOR_RGB[2] * s).clamp(0.0, 255.0),
                            ];
                        }

                        let mut id_code = display_ids[i] + 1;
                        if is_lake {
                            id_code = COMBINED_LAKE_ID_CODE;
                        }
                        if is_sea_cell {
                            id_code = COMBINED_SEA_ID_CODE;
                        }
                        if c.glacier_on || c.floating_ice {
                            id_code = display_ids[i] + 1;
                        }
                        (
                            [
                                col_f[0].round_ties_even().clamp(0.0, 255.0) as u8,
                                col_f[1].round_ties_even().clamp(0.0, 255.0) as u8,
                                col_f[2].round_ties_even().clamp(0.0, 255.0) as u8,
                            ],
                            255u8 - id_code,
                        )
                    })
                    .collect();
                for (c, a) in results {
                    cols.push(c);
                    alp.push(a);
                }
                (cols, alp)
            };

            let padding_px = PADDING_PX * (map_w as f64) / REFERENCE_WIDTH_PX;
            let pgrid = projected_climate_grid(
                grid_h, grid_w, render_proj, &rot, map_w, map_h, padding_px,
            );
            let mut rgb = vec![0u8; map_w * map_h * 3];
            for px in rgb.chunks_exact_mut(3) {
                px.copy_from_slice(&BACKGROUND_RGB);
            }
            let mut alpha = vec![255u8; map_w * map_h];
            fill_rects_rgb_and_alpha(
                &mut rgb,
                &mut alpha,
                map_w,
                &pgrid.x0,
                &pgrid.x1,
                &pgrid.y0,
                &pgrid.y1,
                &colors,
                &cell_alpha,
            );

            let blur_radius = CELL_BLUR_RADIUS_PX * (map_w as f64) / REFERENCE_WIDTH_PX;
            gaussian_blur_rgb(&mut rgb, map_w, map_h, blur_radius);

            let pixel_scale = (map_w as f64 / REFERENCE_WIDTH_PX) * stroke_mult;
            let rivers = RiverSegments {
                from_xyz: &r_from,
                to_xyz: &r_to,
                width_frac: &r_wf,
                color_alpha: &r_ca,
                canyon_fill: &r_cf,
            };
            draw_rivers_rgb(
                &mut rgb,
                map_w,
                map_h,
                render_proj,
                &rot,
                pgrid.scale,
                pgrid.offset_x,
                pgrid.offset_y,
                pixel_scale,
                &rivers,
            );

            if req_proj == Projection::Globe {
                resample_globe_rgba(&rgb, &alpha, map_w, map_h, width, height)
            } else {
                (rgb, alpha)
            }
        })
    });

    package_outputs(py, &out_rgb, Some(&out_alpha), width, height, want_png, want_rgb)
}

#[pyfunction]
#[pyo3(signature = (
    projection,
    width,
    height,
    view_rotation,
    grid_spacing_rad,
    sea_level_m,
    axial_tilt_deg,
    solar_multiplier,
    smooth_sea_ice,
    soft_sea_ice_edge,
    show_mountains,
    show_plains_plateaus,
    tectonics_bbox,
    points,
    elevation,
    lake_depth,
    lake_ice,
    glacier_depth,
    hillshade,
    channel_depth,
    channel_width,
    terrain_relief,
    hydro_points_same,
    hydro_points,
    hydro_is_ocean,
    hydro_is_sea,
    hydro_perennial_lake_ice,
    clim_air_temp,
    clim_ocean_temp,
    clim_lat_deg,
    clim_elev_m,
    clim_is_ocean,
    clim_h,
    clim_w,
    river_from_xyz,
    river_to_xyz,
    river_width_frac,
    river_color_alpha,
    river_canyon_fill,
    want_png,
    want_rgb,
    stroke_scale_override=1.0,
))]
#[allow(clippy::too_many_arguments)]
pub fn render_elevation(
    py: Python<'_>,
    projection: &str,
    width: usize,
    height: usize,
    view_rotation: PyReadonlyArray2<'_, f64>,
    grid_spacing_rad: f64,
    sea_level_m: f64,
    axial_tilt_deg: f64,
    solar_multiplier: f64,
    smooth_sea_ice: bool,
    soft_sea_ice_edge: bool,
    show_mountains: bool,
    show_plains_plateaus: bool,
    tectonics_bbox: Option<[f64; 4]>,
    points: PyReadonlyArray2<'_, f64>,
    elevation: PyReadonlyArray1<'_, f64>,
    lake_depth: PyReadonlyArray1<'_, f64>,
    lake_ice: PyReadonlyArray1<'_, f64>,
    glacier_depth: PyReadonlyArray1<'_, f64>,
    hillshade: Option<PyReadonlyArray1<'_, f64>>,
    channel_depth: Option<PyReadonlyArray1<'_, f64>>,
    channel_width: Option<PyReadonlyArray1<'_, f64>>,
    terrain_relief: Option<PyReadonlyArray1<'_, u8>>,
    hydro_points_same: bool,
    hydro_points: Option<PyReadonlyArray2<'_, f64>>,
    hydro_is_ocean: Option<PyReadonlyArray1<'_, bool>>,
    hydro_is_sea: Option<PyReadonlyArray1<'_, bool>>,
    hydro_perennial_lake_ice: Option<PyReadonlyArray1<'_, bool>>,
    clim_air_temp: PyReadonlyArray1<'_, f64>,
    clim_ocean_temp: PyReadonlyArray1<'_, f64>,
    clim_lat_deg: PyReadonlyArray1<'_, f64>,
    clim_elev_m: Option<PyReadonlyArray1<'_, f64>>,
    clim_is_ocean: Option<PyReadonlyArray1<'_, bool>>,
    clim_h: usize,
    clim_w: usize,
    river_from_xyz: PyReadonlyArray2<'_, f64>,
    river_to_xyz: PyReadonlyArray2<'_, f64>,
    river_width_frac: PyReadonlyArray1<'_, f64>,
    river_color_alpha: PyReadonlyArray1<'_, f64>,
    river_canyon_fill: PyReadonlyArray1<'_, i32>,
    want_png: bool,
    want_rgb: bool,
    stroke_scale_override: f64,
) -> PyResult<(Option<Py<PyBytes>>, Option<Py<PyArray3<u8>>>)> {
    let req_proj = Projection::parse(projection).map_err(PyValueError::new_err)?;
    let rot = to_rot9(&view_rotation)?;
    let pts_v = to_pts3(&points)?;
    let elev_v = elevation.as_slice()?.to_vec();
    let lake_v = lake_depth.as_slice()?.to_vec();
    let lake_ice_v = lake_ice.as_slice()?.to_vec();
    let glac_v = glacier_depth.as_slice()?.to_vec();
    let hill_v = hillshade.map(|a| a.as_slice().unwrap().to_vec());
    let cdepth_v = channel_depth.map(|a| a.as_slice().unwrap().to_vec());
    let cwidth_v = channel_width.map(|a| a.as_slice().unwrap().to_vec());
    let trelief_v = terrain_relief.map(|a| a.as_slice().unwrap().to_vec());

    let hpts_v = hydro_points.map(|a| to_pts3(&a)).transpose()?;
    let hocean_v = hydro_is_ocean.map(|a| a.as_slice().unwrap().to_vec());
    let hsea_v = hydro_is_sea.map(|a| a.as_slice().unwrap().to_vec());
    let hlake_ice_v = hydro_perennial_lake_ice.map(|a| a.as_slice().unwrap().to_vec());

    let c_air = clim_air_temp.as_slice()?.to_vec();
    let c_sst = clim_ocean_temp.as_slice()?.to_vec();
    let c_lat = clim_lat_deg.as_slice()?.to_vec();
    let c_elev = clim_elev_m.map(|a| a.as_slice().unwrap().to_vec());
    let c_ocean = clim_is_ocean.map(|a| a.as_slice().unwrap().to_vec());

    let r_from = to_pts3(&river_from_xyz)?;
    let r_to = to_pts3(&river_to_xyz)?;
    let r_wf = river_width_frac.as_slice()?.to_vec();
    let r_ca = river_color_alpha.as_slice()?.to_vec();
    let r_cf = river_canyon_fill.as_slice()?.to_vec();

    let out_rgb = py.detach(|| {
        render_pool().1.install(|| {
            let (render_proj, map_w, map_h, stroke_mult) = if req_proj == Projection::Globe {
                let (mw, mh) = globe_source_size(width, height);
                let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
                (
                    Projection::Behrmann,
                    mw,
                    mh,
                    width as f64 / (mw as f64 * cos_phi0),
                )
            } else {
                (req_proj, width, height, stroke_scale_override)
            };

            let tree = KdTree::build(&pts_v);
            let hydro_tree = if !hydro_points_same {
                hpts_v.as_ref().map(|hp| KdTree::build(hp))
            } else {
                None
            };
            let resolved_hill = match &hill_v {
                Some(h) => h.clone(),
                None => compute_hillshade_from_tree(
                    &pts_v,
                    &elev_v,
                    &glac_v,
                    cdepth_v.as_deref().unwrap_or(&[]),
                    cwidth_v.as_deref().unwrap_or(&[]),
                    &tree,
                ),
            };
            let any_lake_cap = lake_ice_v.iter().any(|&v| v > 0.0);
            let clim_cont = if any_lake_cap && hlake_ice_v.is_none() {
                c_ocean
                    .as_ref()
                    .map(|co| grid_continentality_dist(co, clim_h, clim_w))
            } else {
                None
            };
            let coarse_sea_ice: Option<Vec<(f64, bool, f64, f64)>> = if !smooth_sea_ice {
                let mut g = Vec::with_capacity(clim_h * clim_w);
                for r in 0..clim_h {
                    let clat = c_lat[r];
                    for c in 0..clim_w {
                        let res = sea_ice_annual_balance_scalar(
                            c_sst[r * clim_w + c],
                            clat,
                            axial_tilt_deg,
                            solar_multiplier,
                            true,
                        );
                        g.push((
                            res.temperature_c,
                            res.perennial,
                            res.fraction,
                            res.perennial_cover,
                        ));
                    }
                }
                Some(g)
            } else {
                None
            };

            let lattice = build_elevation_lattice(grid_spacing_rad);
            let half_dphi = grid_spacing_rad * 0.5;

            struct ElevCell {
                xy: [f64; 2],
                hw: f64,
                hh: f64,
                color: [u8; 3],
            }

            let cells: Vec<ElevCell> = lattice
                .par_iter()
                .map(|&(phi, theta, world_pt)| {
                    let idx = tree.query_nearest(world_pt);
                    let rot_pt = rotate_point(world_pt, &rot);
                    let xy = project_point(render_proj, rot_pt);
                    let (_, center_lon) = xyz_to_latlon(rot_pt);
                    let half_dtheta = grid_spacing_rad / phi.cos().max(1e-3) * 0.5;

                    let mut max_dx = 0.0f64;
                    let mut max_dy = 0.0f64;
                    for dphi_s in [1.0f64, -1.0f64] {
                        let cphi = phi + dphi_s * half_dphi;
                        for dtheta_s in [1.0f64, -1.0f64] {
                            let ctheta = theta + dtheta_s * half_dtheta;
                            let clocal = latlon_to_xyz(cphi, ctheta);
                            let crot = rotate_point(clocal, &rot);
                            let cxy = project_offset(render_proj, crot, center_lon);
                            max_dx = max_dx.max((cxy[0] - xy[0]).abs());
                            max_dy = max_dy.max((cxy[1] - xy[1]).abs());
                        }
                    }

                    let bed = elev_v[idx];
                    let g_depth = glac_v[idx];
                    let l_depth = lake_v[idx];
                    let l_cap = lake_ice_v[idx];
                    let hs = resolved_hill[idx];
                    let elev = bed + g_depth;

                    let hidx = if hydro_points_same {
                        idx
                    } else if let Some(ht) = hydro_tree.as_ref() {
                        ht.query_nearest(world_pt)
                    } else {
                        idx
                    };
                    let is_ocean = hocean_v
                        .as_ref()
                        .map_or(bed <= sea_level_m, |ho| ho[hidx]);
                    let is_sea = hsea_v.as_ref().map_or(false, |hs_arr| hs_arr[hidx]);

                    let (lat_r, lon_r) = xyz_to_latlon(world_pt);
                    let lat_d = lat_r.to_degrees();
                    let lon_d = lon_r.to_degrees();
                    let (crow, ccol) = grid_cell_indices_scalar(lat_d, lon_d, clim_h, clim_w);
                    let cflat = crow * clim_w + ccol;

                    let lake_frozen_raw = if l_cap <= 0.0 {
                        false
                    } else if let Some(pl) = hlake_ice_v.as_ref() {
                        pl[hidx]
                    } else {
                        let water_elev = bed + l_depth + l_cap;
                        let base_elev = c_elev.as_ref().map_or(0.0, |ce| ce[cflat]);
                        let temp = c_air[cflat] - 6.5 * (water_elev - base_elev) / 1000.0;
                        let cont = clim_cont.as_ref().map_or(0.0, |cc| cc[cflat]);
                        let amp = seasonal_temp_amplitude(
                            c_lat[crow],
                            cont,
                            axial_tilt_deg,
                            Some(bed - sea_level_m),
                        );
                        seasonal_fraction_above_scalar(temp, amp, 0.0) == 0.0
                    };
                    let lake_frozen = lake_frozen_raw && !is_ocean;

                    let ocean_water = is_ocean && g_depth == 0.0;
                    let sea_cell = if smooth_sea_ice {
                        let sst = bilinear_sample_climate(&c_sst, clim_h, clim_w, lat_d, lon_d);
                        sea_ice_annual_balance_scalar(
                            sst,
                            lat_d,
                            axial_tilt_deg,
                            solar_multiplier,
                            ocean_water,
                        )
                    } else {
                        let cg = coarse_sea_ice.as_ref().unwrap()[cflat];
                        biomes::SeaIceCell {
                            temperature_c: cg.0,
                            perennial: ocean_water && cg.1,
                            fraction: if ocean_water { cg.2 } else { 0.0 },
                            perennial_cover: if ocean_water { cg.3 } else { 0.0 },
                        }
                    };

                    let base_u8 = elevation_colors_u8(elev, sea_level_m);
                    let mut c_u8 = if hs != 1.0 {
                        [
                            (base_u8[0] as f64 * hs).round_ties_even().clamp(0.0, 255.0) as u8,
                            (base_u8[1] as f64 * hs).round_ties_even().clamp(0.0, 255.0) as u8,
                            (base_u8[2] as f64 * hs).round_ties_even().clamp(0.0, 255.0) as u8,
                        ]
                    } else {
                        base_u8
                    };

                    if let Some(tr) = trelief_v.as_ref() {
                        let code = tr[idx];
                        if show_mountains && code == TERRAIN_MOUNTAIN {
                            c_u8 = [
                                (c_u8[0] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + MOUNTAIN_OVERLAY_RGB[0] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                                (c_u8[1] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + MOUNTAIN_OVERLAY_RGB[1] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                                (c_u8[2] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + MOUNTAIN_OVERLAY_RGB[2] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                            ];
                        }
                        if show_plains_plateaus && code == TERRAIN_PLAINS_PLATEAU {
                            c_u8 = [
                                (c_u8[0] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + PLAINS_PLATEAU_OVERLAY_RGB[0] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                                (c_u8[1] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + PLAINS_PLATEAU_OVERLAY_RGB[1] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                                (c_u8[2] as f64 * (1.0 - TERRAIN_OVERLAY_ALPHA)
                                    + PLAINS_PLATEAU_OVERLAY_RGB[2] * TERRAIN_OVERLAY_ALPHA)
                                    .round_ties_even()
                                    .clamp(0.0, 255.0) as u8,
                            ];
                        }
                    }

                    if l_depth > LAKE_MIN_VISIBLE_DEPTH_M {
                        let lk = if is_sea { SEA_COLOR_RGB } else { LAKE_COLOR_RGB };
                        c_u8 = [lk[0] as u8, lk[1] as u8, lk[2] as u8];
                    }

                    let glacier_on = g_depth >= GLACIER_VISIBLE_DEPTH_M;
                    let mut col_f = [c_u8[0] as f64, c_u8[1] as f64, c_u8[2] as f64];
                    let sea_weight = if glacier_on || lake_frozen {
                        0.0
                    } else if soft_sea_ice_edge {
                        sea_cell.perennial_cover
                    } else if sea_cell.perennial {
                        1.0
                    } else {
                        0.0
                    };
                    if sea_weight > 0.0 {
                        col_f = [
                            col_f[0] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[0] * sea_weight,
                            col_f[1] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[1] * sea_weight,
                            col_f[2] * (1.0 - sea_weight) + GLACIER_COLOR_RGB[2] * sea_weight,
                        ];
                    }
                    if glacier_on || lake_frozen {
                        let glacier_tone =
                            0.93 + 0.07 * ((elev - sea_level_m) / 4500.0).clamp(0.0, 1.0);
                        let shade = if glacier_on { hs * glacier_tone } else { 1.0 };
                        col_f = [
                            (GLACIER_COLOR_RGB[0] * shade).clamp(0.0, 255.0),
                            (GLACIER_COLOR_RGB[1] * shade).clamp(0.0, 255.0),
                            (GLACIER_COLOR_RGB[2] * shade).clamp(0.0, 255.0),
                        ];
                    }
                    let final_u8 = [
                        col_f[0].round_ties_even().clamp(0.0, 255.0) as u8,
                        col_f[1].round_ties_even().clamp(0.0, 255.0) as u8,
                        col_f[2].round_ties_even().clamp(0.0, 255.0) as u8,
                    ];

                    ElevCell {
                        xy,
                        hw: max_dx,
                        hh: max_dy,
                        color: final_u8,
                    }
                })
                .collect();

            let mut min_x = f64::INFINITY;
            let mut min_y = f64::INFINITY;
            let mut max_x = f64::NEG_INFINITY;
            let mut max_y = f64::NEG_INFINITY;
            for c in &cells {
                min_x = min_x.min(c.xy[0]);
                min_y = min_y.min(c.xy[1]);
                max_x = max_x.max(c.xy[0]);
                max_y = max_y.max(c.xy[1]);
            }
            if let Some([tb0, tb1, tb2, tb3]) = tectonics_bbox {
                min_x = min_x.min(tb0);
                min_y = min_y.min(tb1);
                max_x = max_x.max(tb2);
                max_y = max_y.max(tb3);
            }
            let padding_px = PADDING_PX * (map_w as f64) / REFERENCE_WIDTH_PX;
            let (scale, offset_x, offset_y) =
                fit_view([min_x, min_y, max_x, max_y], map_w, map_h, padding_px);

            let iw = map_w as i32;
            let ih = map_h as i32;
            let mut x0 = Vec::with_capacity(cells.len());
            let mut x1 = Vec::with_capacity(cells.len());
            let mut y0 = Vec::with_capacity(cells.len());
            let mut y1 = Vec::with_capacity(cells.len());
            let mut colors = Vec::with_capacity(cells.len());
            let mult = scale * CELL_OVERLAP_FACTOR;
            for c in &cells {
                let px = to_pixel(scale, offset_x, offset_y, c.xy);
                let hw_px = c.hw * mult;
                let hh_px = c.hh * mult;
                x0.push(((px[0] - hw_px).round_ties_even() as i32).clamp(0, iw));
                x1.push(((px[0] + hw_px).round_ties_even() as i32).clamp(0, iw));
                y0.push(((px[1] - hh_px).round_ties_even() as i32).clamp(0, ih));
                y1.push(((px[1] + hh_px).round_ties_even() as i32).clamp(0, ih));
                colors.push(c.color);
            }

            let mut rgb = vec![0u8; map_w * map_h * 3];
            for px in rgb.chunks_exact_mut(3) {
                px.copy_from_slice(&BACKGROUND_RGB);
            }
            fill_rects_rgb(&mut rgb, map_w, &x0, &x1, &y0, &y1, &colors);

            let blur_radius = CELL_BLUR_RADIUS_PX * (map_w as f64) / REFERENCE_WIDTH_PX;
            gaussian_blur_rgb(&mut rgb, map_w, map_h, blur_radius);

            let pixel_scale = (map_w as f64 / REFERENCE_WIDTH_PX) * stroke_mult;
            let rivers = RiverSegments {
                from_xyz: &r_from,
                to_xyz: &r_to,
                width_frac: &r_wf,
                color_alpha: &r_ca,
                canyon_fill: &r_cf,
            };
            draw_rivers_rgb(
                &mut rgb,
                map_w,
                map_h,
                render_proj,
                &rot,
                scale,
                offset_x,
                offset_y,
                pixel_scale,
                &rivers,
            );

            if req_proj == Projection::Globe {
                resample_globe_rgb(&rgb, map_w, map_h, width, height)
            } else {
                rgb
            }
        })
    });

    package_outputs(py, &out_rgb, None, width, height, want_png, want_rgb)
}
