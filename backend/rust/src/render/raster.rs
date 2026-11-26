use std::io::Cursor;

use rayon::prelude::*;

use crate::render::projection::{
    fit_view, project_offset, project_point, rotate_point, sphere_extent, to_pixel, xyz_to_latlon,
    BEHRMANN_STANDARD_PARALLEL, Projection,
};

pub const REFERENCE_WIDTH_PX: f64 = 1100.0;
pub const PADDING_PX: f64 = 20.0;
pub const CELL_OVERLAP_FACTOR: f64 = 1.15;
pub const CELL_BLUR_RADIUS_PX: f64 = 1.0;

pub const BACKGROUND_RGB: [u8; 3] = [11, 16, 32];
pub const LAKE_COLOR_RGB: [f64; 3] = [58.0, 92.0, 122.0];
pub const SEA_COLOR_RGB: [f64; 3] = [42.0, 128.0, 140.0];
pub const GLACIER_COLOR_RGB: [f64; 3] = [221.0, 240.0, 245.0];

pub const MOUNTAIN_OVERLAY_RGB: [f64; 3] = [196.0, 64.0, 56.0];
pub const PLAINS_PLATEAU_OVERLAY_RGB: [f64; 3] = [214.0, 188.0, 74.0];
pub const TERRAIN_OVERLAY_ALPHA: f64 = 0.35;
pub const TERRAIN_MOUNTAIN: u8 = 1;
pub const TERRAIN_PLAINS_PLATEAU: u8 = 2;

pub const OCEAN_PELAGIC_RELIEF_BLEND: f64 = 0.5;
pub const RELIEF_BLEND_MAX: f64 = 0.55;
pub const RELIEF_ELEVATION_RANGE_M: f64 = 6000.0;

pub const RIVER_COLOR_RGB: [f64; 3] = [77.0, 216.0, 230.0];
pub const RIVER_LINE_WIDTH_PX: f64 = 1.0;
pub const RIVER_LINE_WIDTH_MAX_MULT: f64 = 2.0;
pub const RIVER_BLUR_RADIUS_PX: f64 = 0.6;
pub const RIVER_CORE_WIDTH_FRACTION: f64 = 0.5;
pub const RIVER_HALO_FILL_FRACTION: f64 = 0.55;
pub const RIVER_CANYON_SHADE: f64 = 0.15;

pub const ELEVATION_STOP_E: [f64; 11] = [
    -11000.0, -4000.0, -1500.0, -200.0, 0.0, 40.0, 200.0, 1200.0, 3000.0, 6000.0, 9000.0,
];
pub const ELEVATION_STOP_RGB: [[f64; 3]; 11] = [
    [10.0, 10.0, 40.0],
    [15.0, 40.0, 110.0],
    [40.0, 110.0, 190.0],
    [90.0, 160.0, 205.0],
    [150.0, 195.0, 222.0],
    [205.0, 205.0, 155.0],
    [90.0, 150.0, 60.0],
    [170.0, 160.0, 90.0],
    [120.0, 90.0, 60.0],
    [195.0, 188.0, 178.0],
    [222.0, 217.0, 210.0],
];

#[inline]
pub fn elevation_colors_float(elev_m: f64, sea_level_m: f64) -> [f64; 3] {
    let x = elev_m - sea_level_m;
    if x <= ELEVATION_STOP_E[0] {
        return ELEVATION_STOP_RGB[0];
    }
    let last = ELEVATION_STOP_E.len() - 1;
    if x >= ELEVATION_STOP_E[last] {
        return ELEVATION_STOP_RGB[last];
    }
    let mut i = 0usize;
    while i + 1 < ELEVATION_STOP_E.len() && x > ELEVATION_STOP_E[i + 1] {
        i += 1;
    }
    let t = (x - ELEVATION_STOP_E[i]) / (ELEVATION_STOP_E[i + 1] - ELEVATION_STOP_E[i]);
    let a = ELEVATION_STOP_RGB[i];
    let b = ELEVATION_STOP_RGB[i + 1];
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}

#[inline]
pub fn elevation_colors_u8(elev_m: f64, sea_level_m: f64) -> [u8; 3] {
    let f = elevation_colors_float(elev_m, sea_level_m);
    [
        f[0].round_ties_even().clamp(0.0, 255.0) as u8,
        f[1].round_ties_even().clamp(0.0, 255.0) as u8,
        f[2].round_ties_even().clamp(0.0, 255.0) as u8,
    ]
}

/// Serial last-write-wins rectangle fill matching `_fill_rects_kernel` in render_image.py.
pub fn fill_rects_rgb_and_alpha(
    rgb: &mut [u8],
    alpha: &mut [u8],
    width: usize,
    x0: &[i32],
    x1: &[i32],
    y0: &[i32],
    y1: &[i32],
    colors: &[[u8; 3]],
    cell_alpha: &[u8],
) {
    let n = x0.len();
    for i in 0..n {
        let sx0 = x0[i];
        let sx1 = x1[i];
        let sy0 = y0[i];
        let sy1 = y1[i];
        if sx1 > sx0 && sy1 > sy0 {
            let c = colors[i];
            let a = cell_alpha[i];
            for y in sy0..sy1 {
                let row = y as usize * width;
                for x in sx0..sx1 {
                    let idx = row + x as usize;
                    let rgb_idx = idx * 3;
                    rgb[rgb_idx] = c[0];
                    rgb[rgb_idx + 1] = c[1];
                    rgb[rgb_idx + 2] = c[2];
                    alpha[idx] = a;
                }
            }
        }
    }
}

pub fn fill_rects_rgb(
    rgb: &mut [u8],
    width: usize,
    x0: &[i32],
    x1: &[i32],
    y0: &[i32],
    y1: &[i32],
    colors: &[[u8; 3]],
) {
    let n = x0.len();
    for i in 0..n {
        let sx0 = x0[i];
        let sx1 = x1[i];
        let sy0 = y0[i];
        let sy1 = y1[i];
        if sx1 > sx0 && sy1 > sy0 {
            let c = colors[i];
            for y in sy0..sy1 {
                let row = y as usize * width;
                for x in sx0..sx1 {
                    let rgb_idx = (row + x as usize) * 3;
                    rgb[rgb_idx] = c[0];
                    rgb[rgb_idx + 1] = c[1];
                    rgb[rgb_idx + 2] = c[2];
                }
            }
        }
    }
}

pub fn fill_rects_generic(
    pixels: &mut [u8],
    width: usize,
    channels: usize,
    x0: &[i64],
    x1: &[i64],
    y0: &[i64],
    y1: &[i64],
    colors: &[u8],
) {
    let n = x0.len();
    for i in 0..n {
        let sx0 = x0[i];
        let sx1 = x1[i];
        let sy0 = y0[i];
        let sy1 = y1[i];
        if sx1 > sx0 && sy1 > sy0 {
            let c_base = i * channels;
            for y in sy0..sy1 {
                let row = y as usize * width;
                for x in sx0..sx1 {
                    let p_base = (row + x as usize) * channels;
                    pixels[p_base..p_base + channels]
                        .copy_from_slice(&colors[c_base..c_base + channels]);
                }
            }
        }
    }
}

/// Gwosdek et al. (2011) 3-pass extended box filter matching Pillow's `ImageFilter.GaussianBlur`.
fn gwosdek_box_params(sigma: f64) -> (isize, f64) {
    if sigma <= 0.0 {
        return (0, 0.0);
    }
    let sigma2 = (sigma * sigma) / 3.0;
    for l in 0..64 {
        let lf = l as f64;
        let s_l = lf * (lf + 1.0) * (2.0 * lf + 1.0) / 3.0;
        let denom = 2.0 * ((lf + 1.0) * (lf + 1.0) - sigma2);
        if denom <= 0.0 {
            continue;
        }
        let alpha = (sigma2 * (2.0 * lf + 1.0) - s_l) / denom;
        if (0.0..=1.0).contains(&alpha) {
            return (l, alpha);
        }
    }
    let l = ((12.0 * sigma2 + 1.0).sqrt() - 1.0) * 0.5;
    let li = l.floor().max(0.0) as isize;
    (li, (l - li as f64).clamp(0.0, 1.0))
}

pub fn gaussian_blur_rgb(rgb: &mut [u8], width: usize, height: usize, radius: f64) {
    if radius <= 0.0 || width == 0 || height == 0 {
        return;
    }
    let (l, alpha) = gwosdek_box_params(radius);
    if l == 0 && alpha < 1e-9 {
        return;
    }
    let norm = 1.0 / ((2 * l + 1) as f64 + 2.0 * alpha);
    let n = width * height * 3;
    let mut buf_a: Vec<f64> = rgb.par_iter().map(|&v| v as f64).collect();
    let mut buf_b = vec![0.0f64; n];

    let max_x = width as isize - 1;
    let max_y = height as isize - 1;

    for _ in 0..3 {
        // Horizontal pass: buf_a -> buf_b
        buf_b
            .par_chunks_mut(width * 3)
            .zip(buf_a.par_chunks(width * 3))
            .for_each(|(out_row, in_row)| {
                for x in 0..width {
                    let xi = x as isize;
                    let mut s0 = 0.0;
                    let mut s1 = 0.0;
                    let mut s2 = 0.0;
                    for dx in -l..=l {
                        let sx = (xi + dx).clamp(0, max_x) as usize * 3;
                        s0 += in_row[sx];
                        s1 += in_row[sx + 1];
                        s2 += in_row[sx + 2];
                    }
                    let sx_l = (xi - l - 1).clamp(0, max_x) as usize * 3;
                    let sx_r = (xi + l + 1).clamp(0, max_x) as usize * 3;
                    s0 += alpha * (in_row[sx_l] + in_row[sx_r]);
                    s1 += alpha * (in_row[sx_l + 1] + in_row[sx_r + 1]);
                    s2 += alpha * (in_row[sx_l + 2] + in_row[sx_r + 2]);
                    let ox = x * 3;
                    out_row[ox] = s0 * norm;
                    out_row[ox + 1] = s1 * norm;
                    out_row[ox + 2] = s2 * norm;
                }
            });

        // Vertical pass: buf_b -> buf_a (chunked by rows in parallel)
        let row_stride = width * 3;
        buf_a
            .par_chunks_mut(row_stride)
            .enumerate()
            .for_each(|(y, out_row)| {
                let yi = y as isize;
                let yl = (yi - l - 1).clamp(0, max_y) as usize * row_stride;
                let yr = (yi + l + 1).clamp(0, max_y) as usize * row_stride;
                for i in 0..row_stride {
                    out_row[i] = alpha * (buf_b[yl + i] + buf_b[yr + i]);
                }
                for dy in -l..=l {
                    let sy = (yi + dy).clamp(0, max_y) as usize * row_stride;
                    for i in 0..row_stride {
                        out_row[i] += buf_b[sy + i];
                    }
                }
                for val in out_row.iter_mut() {
                    *val *= norm;
                }
            });
    }

    rgb.par_iter_mut()
        .zip(buf_a.par_iter())
        .for_each(|(dst, &v)| {
            *dst = v.round_ties_even().clamp(0.0, 255.0) as u8;
        });
}

pub fn gaussian_blur_gray(gray: &[u8], width: usize, height: usize, radius: f64) -> Vec<u8> {
    if width == 0 || height == 0 {
        return Vec::new();
    }
    if radius <= 0.0 {
        return gray.to_vec();
    }
    let (l, alpha) = gwosdek_box_params(radius);
    if l == 0 && alpha < 1e-9 {
        return gray.to_vec();
    }
    let norm = 1.0 / ((2 * l + 1) as f64 + 2.0 * alpha);
    let n = width * height;
    let mut buf_a: Vec<f64> = gray.par_iter().map(|&v| v as f64).collect();
    let mut buf_b = vec![0.0f64; n];

    let max_x = width as isize - 1;
    let max_y = height as isize - 1;

    for _ in 0..3 {
        buf_b
            .par_chunks_mut(width)
            .zip(buf_a.par_chunks(width))
            .for_each(|(out_row, in_row)| {
                for x in 0..width {
                    let xi = x as isize;
                    let mut s = 0.0;
                    for dx in -l..=l {
                        s += in_row[(xi + dx).clamp(0, max_x) as usize];
                    }
                    let sl = (xi - l - 1).clamp(0, max_x) as usize;
                    let sr = (xi + l + 1).clamp(0, max_x) as usize;
                    s += alpha * (in_row[sl] + in_row[sr]);
                    out_row[x] = s * norm;
                }
            });

        buf_a
            .par_chunks_mut(width)
            .enumerate()
            .for_each(|(y, out_row)| {
                let yi = y as isize;
                let yl = (yi - l - 1).clamp(0, max_y) as usize * width;
                let yr = (yi + l + 1).clamp(0, max_y) as usize * width;
                for x in 0..width {
                    out_row[x] = alpha * (buf_b[yl + x] + buf_b[yr + x]);
                }
                for dy in -l..=l {
                    let sy = (yi + dy).clamp(0, max_y) as usize * width;
                    for x in 0..width {
                        out_row[x] += buf_b[sy + x];
                    }
                }
                for val in out_row.iter_mut() {
                    *val *= norm;
                }
            });
    }

    buf_a
        .par_iter()
        .map(|&v| v.round_ties_even().clamp(0.0, 255.0) as u8)
        .collect()
}

#[inline]
fn set_gray_px(buf: &mut [u8], width: usize, height: usize, x: i32, y: i32, val: u8) {
    if x >= 0 && y >= 0 && (x as usize) < width && (y as usize) < height {
        buf[y as usize * width + x as usize] = val;
    }
}

/// Rasterize a line matching Pillow `ImageDraw.line` for 1px and thick strokes (`width >= 1`).
fn draw_gray_line(
    buf: &mut [u8],
    width: usize,
    height: usize,
    x1f: f64,
    y1f: f64,
    x2f: f64,
    y2f: f64,
    val: u8,
    line_w: i32,
) {
    if line_w <= 1 {
        let mut x0 = x1f as i32;
        let mut y0 = y1f as i32;
        let x1 = x2f as i32;
        let y1 = y2f as i32;
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -1 };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -1 };
        let mut err = dx + dy;
        loop {
            set_gray_px(buf, width, height, x0, y0, val);
            if x0 == x1 && y0 == y1 {
                break;
            }
            let e2 = 2 * err;
            if e2 >= dy {
                err += dy;
                x0 += sx;
            }
            if e2 <= dx {
                err += dx;
                y0 += sy;
            }
        }
        return;
    }

    // Thick segment matching Pillow's thick-line footprint: filled capsule/quad + circular caps.
    let dx = x2f - x1f;
    let dy = y2f - y1f;
    let len = (dx * dx + dy * dy).sqrt();
    let r = (line_w as f64) * 0.5;
    let r2 = r * r + 1e-6;
    let min_x = (x1f.min(x2f) - r - 1.0).floor() as i32;
    let max_x = (x1f.max(x2f) + r + 1.0).ceil() as i32;
    let min_y = (y1f.min(y2f) - r - 1.0).floor() as i32;
    let max_y = (y1f.max(y2f) + r + 1.0).ceil() as i32;

    let min_xc = min_x.max(0);
    let max_xc = max_x.min(width as i32 - 1);
    let min_yc = min_y.max(0);
    let max_yc = max_y.min(height as i32 - 1);

    for py in min_yc..=max_yc {
        let cy = py as f64 + 0.5;
        for px in min_xc..=max_xc {
            let cx = px as f64 + 0.5;
            let d2 = if len < 1e-9 {
                (cx - x1f).powi(2) + (cy - y1f).powi(2)
            } else {
                let t = (((cx - x1f) * dx + (cy - y1f) * dy) / (len * len)).clamp(0.0, 1.0);
                let proj_x = x1f + t * dx;
                let proj_y = y1f + t * dy;
                (cx - proj_x).powi(2) + (cy - proj_y).powi(2)
            };
            if d2 <= r2 {
                buf[py as usize * width + px as usize] = val;
            }
        }
    }
}

pub struct RiverSegments<'a> {
    pub from_xyz: &'a [[f64; 3]],
    pub to_xyz: &'a [[f64; 3]],
    pub width_frac: &'a [f64],
    pub color_alpha: &'a [f64],
    pub canyon_fill: &'a [i32],
}

pub fn draw_rivers_rgb(
    rgb: &mut [u8],
    width: usize,
    height: usize,
    proj: Projection,
    rot: &[f64; 9],
    scale: f64,
    offset_x: f64,
    offset_y: f64,
    pixel_scale: f64,
    rivers: &RiverSegments<'_>,
) {
    let m = rivers.from_xyz.len();
    if m == 0 {
        return;
    }
    let mut mask = vec![0u8; width * height];
    let mut canyon_mask = vec![0u8; width * height];

    for i in 0..m {
        let from_p = rotate_point(rivers.from_xyz[i], rot);
        let to_p = rotate_point(rivers.to_xyz[i], rot);
        let (_, from_lon) = xyz_to_latlon(from_p);
        let from_xy = project_point(proj, from_p);
        let to_xy = project_offset(proj, to_p, from_lon);
        let from_px = to_pixel(scale, offset_x, offset_y, from_xy);
        let to_px = to_pixel(scale, offset_x, offset_y, to_xy);

        let canyon = rivers.canyon_fill[i];
        if canyon > 0 {
            draw_gray_line(
                &mut canyon_mask,
                width,
                height,
                from_px[0],
                from_px[1],
                to_px[0],
                to_px[1],
                canyon.clamp(0, 255) as u8,
                1,
            );
        }

        let f = (255.0 * rivers.color_alpha[i]).round_ties_even().clamp(0.0, 255.0) as i32;
        if f == 0 {
            continue;
        }
        let w = (pixel_scale
            * RIVER_LINE_WIDTH_PX
            * (1.0 + (RIVER_LINE_WIDTH_MAX_MULT - 1.0) * rivers.width_frac[i]))
            .round_ties_even()
            .max(1.0) as i32;
        if w > 1 {
            let halo_f = ((f as f64) * RIVER_HALO_FILL_FRACTION)
                .round_ties_even()
                .clamp(0.0, 255.0) as u8;
            let core_w = ((w as f64) * RIVER_CORE_WIDTH_FRACTION)
                .round_ties_even()
                .max(1.0) as i32;
            draw_gray_line(
                &mut mask,
                width,
                height,
                from_px[0],
                from_px[1],
                to_px[0],
                to_px[1],
                halo_f,
                w,
            );
            draw_gray_line(
                &mut mask,
                width,
                height,
                from_px[0],
                from_px[1],
                to_px[0],
                to_px[1],
                f as u8,
                core_w,
            );
        } else {
            draw_gray_line(
                &mut mask,
                width,
                height,
                from_px[0],
                from_px[1],
                to_px[0],
                to_px[1],
                f as u8,
                w,
            );
        }
    }

    let blur_r = RIVER_BLUR_RADIUS_PX * pixel_scale;
    let blurred_canyon = gaussian_blur_gray(&canyon_mask, width, height, blur_r);
    let blurred_alpha = gaussian_blur_gray(&mask, width, height, blur_r);

    rgb.par_chunks_mut(width * 3)
        .zip(blurred_canyon.par_chunks(width))
        .zip(blurred_alpha.par_chunks(width))
        .for_each(|((rgb_row, can_row), alp_row)| {
            for x in 0..width {
                let c = can_row[x] as f64 / 255.0;
                let a = alp_row[x] as f64 / 255.0;
                if c == 0.0 && a == 0.0 {
                    continue;
                }
                let canyon_mult = 1.0 - RIVER_CANYON_SHADE * c;
                let inv_a = 1.0 - a;
                let ox = x * 3;
                rgb_row[ox] = ((rgb_row[ox] as f64 * canyon_mult) * inv_a
                    + RIVER_COLOR_RGB[0] * a)
                    .round_ties_even()
                    .clamp(0.0, 255.0) as u8;
                rgb_row[ox + 1] = ((rgb_row[ox + 1] as f64 * canyon_mult) * inv_a
                    + RIVER_COLOR_RGB[1] * a)
                    .round_ties_even()
                    .clamp(0.0, 255.0) as u8;
                rgb_row[ox + 2] = ((rgb_row[ox + 2] as f64 * canyon_mult) * inv_a
                    + RIVER_COLOR_RGB[2] * a)
                    .round_ties_even()
                    .clamp(0.0, 255.0) as u8;
            }
        });
}

pub fn globe_transform(width: usize, height: usize) -> (f64, f64, f64) {
    let padding_px = PADDING_PX * (width as f64) / REFERENCE_WIDTH_PX;
    fit_view(sphere_extent(Projection::Globe), width, height, padding_px)
}

pub fn globe_source_size(width: usize, height: usize) -> (usize, usize) {
    let radius_px = globe_transform(width, height).0.max(1.0);
    let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
    let mut map_width =
        ((2.0 * std::f64::consts::PI * radius_px) / (1.0 - 2.0 * PADDING_PX / REFERENCE_WIDTH_PX))
            .ceil() as usize;
    map_width = map_width.min(2 * width.max(height));
    let map_padding_px = PADDING_PX * (map_width as f64) / REFERENCE_WIDTH_PX;
    let map_height = (((map_width as f64 - 2.0 * map_padding_px)
        / (std::f64::consts::PI * cos_phi0 * cos_phi0))
        + 2.0 * map_padding_px)
        .ceil() as usize;
    (map_width, map_height)
}

pub fn resample_globe_rgba(
    src_rgb: &[u8],
    src_alpha: &[u8],
    map_width: usize,
    map_height: usize,
    width: usize,
    height: usize,
) -> (Vec<u8>, Vec<u8>) {
    let (scale, offset_x, offset_y) = globe_transform(width, height);
    let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
    let map_padding_px = PADDING_PX * (map_width as f64) / REFERENCE_WIDTH_PX;
    let (map_scale, map_offset_x, map_offset_y) = fit_view(
        sphere_extent(Projection::Behrmann),
        map_width,
        map_height,
        map_padding_px,
    );

    let bg_rgb = [src_rgb[0], src_rgb[1], src_rgb[2]];
    let bg_a = src_alpha[0];
    let mut out_rgb = vec![0u8; width * height * 3];
    let mut out_alpha = vec![bg_a; width * height];

    out_rgb
        .par_chunks_mut(width * 3)
        .zip(out_alpha.par_chunks_mut(width))
        .enumerate()
        .for_each(|(py, (rgb_row, alp_row))| {
            let y = (offset_y - (py as f64 + 0.5)) / scale;
            let y2 = y * y;
            let my = ((map_offset_y - map_scale * y / cos_phi0).floor() as isize)
                .clamp(0, map_height as isize - 1) as usize;
            let src_row = my * map_width;
            for px in 0..width {
                let x = (px as f64 + 0.5 - offset_x) / scale;
                let r2 = x * x + y2;
                let ox = px * 3;
                if r2 <= 1.0 {
                    let depth = (1.0 - r2).clamp(0.0, 1.0).sqrt();
                    let lon = x.atan2(depth);
                    let mx = ((map_scale * lon * cos_phi0 + map_offset_x).floor() as isize)
                        .clamp(0, map_width as isize - 1) as usize;
                    let sidx = src_row + mx;
                    let srgb = sidx * 3;
                    rgb_row[ox] = src_rgb[srgb];
                    rgb_row[ox + 1] = src_rgb[srgb + 1];
                    rgb_row[ox + 2] = src_rgb[srgb + 2];
                    alp_row[px] = src_alpha[sidx];
                } else {
                    rgb_row[ox] = bg_rgb[0];
                    rgb_row[ox + 1] = bg_rgb[1];
                    rgb_row[ox + 2] = bg_rgb[2];
                }
            }
        });
    (out_rgb, out_alpha)
}

pub fn resample_globe_rgb(
    src_rgb: &[u8],
    map_width: usize,
    map_height: usize,
    width: usize,
    height: usize,
) -> Vec<u8> {
    let (scale, offset_x, offset_y) = globe_transform(width, height);
    let cos_phi0 = BEHRMANN_STANDARD_PARALLEL.cos();
    let map_padding_px = PADDING_PX * (map_width as f64) / REFERENCE_WIDTH_PX;
    let (map_scale, map_offset_x, map_offset_y) = fit_view(
        sphere_extent(Projection::Behrmann),
        map_width,
        map_height,
        map_padding_px,
    );

    let bg_rgb = [src_rgb[0], src_rgb[1], src_rgb[2]];
    let mut out_rgb = vec![0u8; width * height * 3];

    out_rgb
        .par_chunks_mut(width * 3)
        .enumerate()
        .for_each(|(py, rgb_row)| {
            let y = (offset_y - (py as f64 + 0.5)) / scale;
            let y2 = y * y;
            let my = ((map_offset_y - map_scale * y / cos_phi0).floor() as isize)
                .clamp(0, map_height as isize - 1) as usize;
            let src_row = my * map_width;
            for px in 0..width {
                let x = (px as f64 + 0.5 - offset_x) / scale;
                let r2 = x * x + y2;
                let ox = px * 3;
                if r2 <= 1.0 {
                    let depth = (1.0 - r2).clamp(0.0, 1.0).sqrt();
                    let lon = x.atan2(depth);
                    let mx = ((map_scale * lon * cos_phi0 + map_offset_x).floor() as isize)
                        .clamp(0, map_width as isize - 1) as usize;
                    let srgb = (src_row + mx) * 3;
                    rgb_row[ox] = src_rgb[srgb];
                    rgb_row[ox + 1] = src_rgb[srgb + 1];
                    rgb_row[ox + 2] = src_rgb[srgb + 2];
                } else {
                    rgb_row[ox] = bg_rgb[0];
                    rgb_row[ox + 1] = bg_rgb[1];
                    rgb_row[ox + 2] = bg_rgb[2];
                }
            }
        });
    out_rgb
}

pub fn encode_png_rgb(rgb: &[u8], width: usize, height: usize) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(width * height);
    {
        let mut enc = png::Encoder::new(Cursor::new(&mut out), width as u32, height as u32);
        enc.set_color(png::ColorType::Rgb);
        enc.set_depth(png::BitDepth::Eight);
        enc.set_compression(png::Compression::Fast);
        enc.set_filter(png::FilterType::Sub);
        let mut writer = enc.write_header().map_err(|e| e.to_string())?;
        writer.write_image_data(rgb).map_err(|e| e.to_string())?;
    }
    Ok(out)
}

pub fn encode_png_rgba(
    rgb: &[u8],
    alpha: &[u8],
    width: usize,
    height: usize,
) -> Result<Vec<u8>, String> {
    let n = width * height;
    let mut rgba = vec![0u8; n * 4];
    rgba.par_chunks_mut(width * 4)
        .zip(rgb.par_chunks(width * 3))
        .zip(alpha.par_chunks(width))
        .for_each(|((dst, r_row), a_row)| {
            for x in 0..width {
                let d = x * 4;
                let s = x * 3;
                dst[d] = r_row[s];
                dst[d + 1] = r_row[s + 1];
                dst[d + 2] = r_row[s + 2];
                dst[d + 3] = a_row[x];
            }
        });

    let mut out = Vec::with_capacity(n);
    {
        let mut enc = png::Encoder::new(Cursor::new(&mut out), width as u32, height as u32);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        enc.set_compression(png::Compression::Fast);
        enc.set_filter(png::FilterType::Sub);
        let mut writer = enc.write_header().map_err(|e| e.to_string())?;
        writer.write_image_data(&rgba).map_err(|e| e.to_string())?;
    }
    Ok(out)
}

pub fn crop_even_rgb24(rgb: &[u8], width: usize, height: usize) -> (Vec<u8>, usize, usize) {
    let enc_w = width - (width % 2);
    let enc_h = height - (height % 2);
    if enc_w == width && enc_h == height {
        return (rgb.to_vec(), enc_w, enc_h);
    }
    let mut out = vec![0u8; enc_h * enc_w * 3];
    for y in 0..enc_h {
        let src = y * width * 3;
        let dst = y * enc_w * 3;
        out[dst..dst + enc_w * 3].copy_from_slice(&rgb[src..src + enc_w * 3]);
    }
    (out, enc_w, enc_h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blur_constant_image_is_invariant() {
        let mut rgb = vec![100u8; 16 * 16 * 3];
        gaussian_blur_rgb(&mut rgb, 16, 16, 1.0);
        assert!(rgb.iter().all(|&v| v == 100));
    }

    #[test]
    fn png_roundtrip_header_is_valid() {
        let rgb = vec![10u8; 8 * 6 * 3];
        let png = encode_png_rgb(&rgb, 8, 6).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
    }
}
