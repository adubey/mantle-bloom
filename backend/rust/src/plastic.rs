use std::collections::HashMap;

use rayon::prelude::*;

use crate::climate::{
    clip_min, compensated_dot, floating_water_depth_scalar, is_grounded_scalar, node_rates_scalar,
    NeumaierSum, FLOTATION_RATIO, GLACIER_MAX_DEPTH_M, ICE_FLOTATION_PASSES,
    ICE_ROOT_MAX_ITERATIONS, ICE_ROOT_TOLERANCE_M,
};
use crate::distance::IndexedHeap;

pub struct LocalDomain {
    pub nodes: Vec<usize>,
    pub indptr: Vec<u32>,
    pub nbr: Vec<i32>,
    pub length: Vec<f64>,
    pub nbr_phi: Vec<f64>,
    pub nbr_bed: Vec<f64>,
    pub nbr_cap: Vec<f64>,
    pub nbr_ocean: Vec<bool>,
    pub phi: Vec<f64>,
    pub bed: Vec<f64>,
    pub c: Vec<f64>,
    pub area: Vec<f64>,
    pub water_depth: Vec<f64>,
    pub ocean: Vec<bool>,
    pub precip: Vec<f64>,
    pub temp_c: Vec<f64>,
    pub amp_c: Vec<f64>,
    pub frictional: Option<Vec<f64>>,
    pub bed_frozen_snow: Option<Vec<f64>>,
    pub bare_acc: Vec<f64>,
    pub bare_sm: Vec<f64>,
    pub has_wall: bool,
    pub phi_min: f64,
    pub phi_max: f64,
}

pub struct SweepScratch {
    pub surface: Vec<f64>,
    pub done: Vec<bool>,
    pub heap: IndexedHeap,
    pub phi: Vec<f64>,
    pub active: Vec<bool>,
}

impl SweepScratch {
    pub fn new(n: usize) -> Self {
        Self {
            surface: vec![f64::NAN; n],
            done: vec![false; n],
            heap: IndexedHeap::new(n),
            phi: vec![0.0; n],
            active: vec![false; n],
        }
    }
}

#[derive(Clone)]
pub struct IceShape {
    pub d: f64,
    pub active: Vec<bool>,
    pub thickness: Vec<f64>,
    pub depth: Vec<f64>,
    pub cover: Vec<f64>,
    pub volume: f64,
    pub rate_acc: Vec<f64>,
    pub rate_sm: Vec<f64>,
    pub rate_bm: Vec<f64>,
}

impl IceShape {
    fn scaled(&self, factor: f64) -> Self {
        Self {
            d: self.d,
            active: self.active.clone(),
            thickness: self.thickness.iter().map(|&v| v * factor).collect(),
            depth: self.depth.iter().map(|&v| v * factor).collect(),
            cover: self.cover.clone(),
            volume: self.volume * factor,
            rate_acc: self.rate_acc.clone(),
            rate_sm: self.rate_sm.clone(),
            rate_bm: self.rate_bm.clone(),
        }
    }

    pub fn budget_rates(&self, dom: &LocalDomain) -> (f64, f64, f64) {
        let mut acc_cov = NeumaierSum::default();
        let mut acc_bare = NeumaierSum::default();
        let mut sm_cov = NeumaierSum::default();
        let mut sm_bare = NeumaierSum::default();
        let mut bm_cov = NeumaierSum::default();
        for i in 0..dom.nodes.len() {
            let cov = self.cover[i];
            let area = dom.area[i];
            let covered = cov * area;
            acc_cov.add(covered * self.rate_acc[i]);
            sm_cov.add(covered * self.rate_sm[i]);
            bm_cov.add(covered * self.rate_bm[i]);
            let acc0 = dom.bare_acc[i];
            let sm0 = dom.bare_sm[i];
            if acc0 - sm0 > 0.0 {
                let bare = (1.0 - cov.min(1.0)) * area;
                acc_bare.add(bare * acc0);
                sm_bare.add(bare * sm0);
            }
        }
        (
            acc_cov.finish() + acc_bare.finish(),
            sm_cov.finish() + sm_bare.finish(),
            bm_cov.finish(),
        )
    }
}

#[inline]
fn upwind(s_up: f64, b_up: f64, b_i: f64, c_i: f64, step_m: f64) -> f64 {
    let x = (2.0 * s_up - b_i - b_up).max(0.0);
    let s = 0.5 * ((b_i + b_up) + (x * x + 8.0 * c_i * step_m).sqrt());
    s.max(b_i)
}

#[inline]
#[allow(clippy::too_many_arguments)]
fn margin_point(
    i: usize,
    e: usize,
    j: i32,
    phi: &[f64],
    nbr_phi: &[f64],
    nbr_bed: &[f64],
    nbr_cap: &[f64],
    nbr_ocean: &[bool],
    ocean: &[bool],
    bed: &[f64],
    sea_level: f64,
    d: f64,
) -> (f64, f64, f64) {
    let mut t = if j >= 0 && phi[j as usize] <= d {
        0.5
    } else if nbr_phi[e] > phi[i] {
        (d - phi[i]) / (nbr_phi[e] - phi[i])
    } else {
        1.0
    };
    if t < 0.0 {
        t = 0.0;
    }
    if t > nbr_cap[e] {
        t = nbr_cap[e];
    }
    let b_m = bed[i] + t * (nbr_bed[e] - bed[i]);
    let wet = nbr_ocean[e] || ocean[i];
    let h_m = if wet {
        FLOTATION_RATIO * (sea_level - b_m).max(0.0)
    } else {
        0.0
    };
    (t, b_m, h_m)
}

#[inline]
fn segment_mean(h0: f64, h1: f64) -> f64 {
    let a = h0 * h0;
    let b = h1 * h1;
    if (b - a).abs() <= 1e-9 * (a + b) || a + b == 0.0 {
        0.5 * (h0 + h1)
    } else {
        (2.0 / 3.0) * (b * b.sqrt() - a * a.sqrt()) / (b - a)
    }
}

/// Fills `scratch.surface` in-place and returns `(depth, cover)`.
fn plastic_sweep(
    dom: &LocalDomain,
    sea_level: f64,
    d: f64,
    max_thickness: f64,
    scratch: &mut SweepScratch,
) -> (Vec<f64>, Vec<f64>) {
    let n = dom.nodes.len();
    let phi = &scratch.phi;
    let active = &scratch.active;
    scratch.surface.fill(f64::INFINITY);
    scratch.done.fill(false);
    scratch.heap.reset();

    for i in 0..n {
        if !active[i] {
            continue;
        }
        let mut best = f64::INFINITY;
        let start = dom.indptr[i] as usize;
        let end = dom.indptr[i + 1] as usize;
        for e in start..end {
            let j = dom.nbr[e];
            if j >= 0 && active[j as usize] {
                continue;
            }
            let (t, b_m, h_m) = margin_point(
                i,
                e,
                j,
                phi,
                &dom.nbr_phi,
                &dom.nbr_bed,
                &dom.nbr_cap,
                &dom.nbr_ocean,
                &dom.ocean,
                &dom.bed,
                sea_level,
                d,
            );
            let cand = upwind(b_m + h_m, b_m, dom.bed[i], dom.c[i], t * dom.length[e]);
            if cand < best {
                best = cand;
            }
        }
        if best < f64::INFINITY {
            scratch.surface[i] = best;
            scratch.heap.push_or_decrease(i as u32, &scratch.surface);
        }
    }

    while let Some(i_u32) = scratch.heap.pop(&scratch.surface) {
        let i = i_u32 as usize;
        scratch.done[i] = true;
        let s = scratch.surface[i];
        let start = dom.indptr[i] as usize;
        let end = dom.indptr[i + 1] as usize;
        for e in start..end {
            let j_raw = dom.nbr[e];
            if j_raw < 0 {
                continue;
            }
            let j = j_raw as usize;
            if !active[j] || scratch.done[j] {
                continue;
            }
            let candidate = upwind(s, dom.bed[i], dom.bed[j], dom.c[j], dom.length[e]);
            if candidate < scratch.surface[j] {
                scratch.surface[j] = candidate;
                scratch.heap.push_or_decrease(j as u32, &scratch.surface);
            }
        }
    }

    for i in 0..n {
        if !active[i] {
            scratch.surface[i] = f64::NAN;
        } else if !scratch.done[i] {
            scratch.surface[i] = dom.bed[i];
        } else {
            scratch.surface[i] = scratch.surface[i].min(dom.bed[i] + max_thickness);
        }
    }

    let mut depth = vec![0.0_f64; n];
    let mut cover = vec![0.0_f64; n];
    for i in 0..n {
        if !active[i] {
            continue;
        }
        let h_i = (scratch.surface[i] - dom.bed[i]).max(0.0);
        let start = dom.indptr[i] as usize;
        let end = dom.indptr[i + 1] as usize;
        let deg_usize = end - start;
        if deg_usize == 0 {
            depth[i] += h_i;
            cover[i] += 1.0;
            continue;
        }
        let degree = deg_usize as f64;
        for e in start..end {
            let j = dom.nbr[e];
            if j >= 0 && active[j as usize] {
                let ju = j as usize;
                let h_j = (scratch.surface[ju] - dom.bed[ju]).max(0.0);
                depth[i] += segment_mean(h_i, (0.5 * (h_i * h_i + h_j * h_j)).sqrt()) / degree;
                cover[i] += 1.0 / degree;
                continue;
            }
            let (t, _, h_m) = margin_point(
                i,
                e,
                j,
                phi,
                &dom.nbr_phi,
                &dom.nbr_bed,
                &dom.nbr_cap,
                &dom.nbr_ocean,
                &dom.ocean,
                &dom.bed,
                sea_level,
                d,
            );
            if t <= 0.5 {
                depth[i] += 2.0 * t * segment_mean(h_i, h_m) / degree;
                cover[i] += 2.0 * t / degree;
                continue;
            }
            let ju = j as usize;
            let h_half = (h_i * h_i + (0.5 / t) * (h_m * h_m - h_i * h_i))
                .max(0.0)
                .sqrt();
            depth[i] += segment_mean(h_i, h_half) / degree;
            cover[i] += 1.0 / degree;
            let degree_j = (dom.indptr[ju + 1] - dom.indptr[ju]) as f64;
            depth[ju] += 2.0 * (t - 0.5) * segment_mean(h_half, h_m) / degree_j;
            cover[ju] += 2.0 * (t - 0.5) / degree_j;
        }
    }

    (depth, cover)
}

fn build_shape(dom: &LocalDomain, d: f64, sea_level_m: f64, scratch: &mut SweepScratch) -> IceShape {
    let n = dom.nodes.len();
    scratch.phi.copy_from_slice(&dom.phi);
    for i in 0..n {
        scratch.active[i] = dom.phi[i] <= d;
    }

    let mut depth = Vec::new();
    let mut cover = Vec::new();

    for attempt in 0..=ICE_FLOTATION_PASSES {
        let (mut dep, mut cov) = plastic_sweep(dom, sea_level_m, d, GLACIER_MAX_DEPTH_M, scratch);
        let mut any_floating = false;
        for i in 0..n {
            if dep[i] > 0.0
                && dom.water_depth[i] > 0.0
                && !is_grounded_scalar(dep[i], dom.bed[i], sea_level_m, dom.ocean[i])
            {
                any_floating = true;
                scratch.active[i] = false;
                if attempt < ICE_FLOTATION_PASSES {
                    scratch.phi[i] = f64::NEG_INFINITY;
                } else {
                    dep[i] = 0.0;
                    cov[i] = 0.0;
                }
            }
        }
        depth = dep;
        cover = cov;
        if !any_floating {
            break;
        }
    }

    let active = scratch.active.clone();
    let mut thickness = vec![0.0_f64; n];
    let mut rate_acc = Vec::with_capacity(n);
    let mut rate_sm = Vec::with_capacity(n);
    let mut rate_bm = Vec::with_capacity(n);

    for i in 0..n {
        let holds = active[i] || (cover[i] > 0.0 && scratch.phi[i] > f64::NEG_INFINITY);
        if !holds {
            depth[i] = 0.0;
            cover[i] = 0.0;
        }
        let cov_i = cover[i];
        let mean_t = if cov_i > 0.0 { depth[i] / cov_i } else { 0.0 };
        thickness[i] = if active[i] {
            clip_min(scratch.surface[i] - dom.bed[i], 0.0)
        } else {
            mean_t
        };
        if mean_t == 0.0 {
            rate_acc.push(dom.bare_acc[i]);
            rate_sm.push(dom.bare_sm[i]);
            rate_bm.push(0.0);
        } else {
            let fric = dom.frictional.as_ref().map_or(0.0, |v| v[i]);
            let bfs = dom.bed_frozen_snow.as_ref().map(|v| v[i]);
            let (acc, sm, bm) = node_rates_scalar(
                mean_t,
                dom.temp_c[i],
                dom.amp_c[i],
                dom.precip[i],
                fric,
                bfs,
            );
            rate_acc.push(acc);
            rate_sm.push(sm);
            rate_bm.push(bm);
        }
    }
    let volume = compensated_dot(&depth, &dom.area);

    IceShape {
        d,
        active,
        thickness,
        depth,
        cover,
        volume,
        rate_acc,
        rate_sm,
        rate_bm,
    }
}

struct ShapeCache<'a> {
    dom: &'a LocalDomain,
    sea_level_m: f64,
    scratch: SweepScratch,
    map: HashMap<u64, usize>,
    shapes: Vec<IceShape>,
}

impl<'a> ShapeCache<'a> {
    fn new(dom: &'a LocalDomain, sea_level_m: f64) -> Self {
        Self {
            dom,
            sea_level_m,
            scratch: SweepScratch::new(dom.nodes.len()),
            map: HashMap::with_capacity(24),
            shapes: Vec::with_capacity(24),
        }
    }

    fn get(&mut self, d: f64) -> &IceShape {
        let key = d.to_bits();
        if let Some(&idx) = self.map.get(&key) {
            return &self.shapes[idx];
        }
        let s = build_shape(self.dom, d, self.sea_level_m, &mut self.scratch);
        let idx = self.shapes.len();
        self.shapes.push(s);
        self.map.insert(key, idx);
        &self.shapes[idx]
    }
}

struct VolumeShortfall<'a, 'b> {
    cache: &'a mut ShapeCache<'b>,
    volume_m3: f64,
    short_up_to: f64,
    holds_from: f64,
}

impl<'a, 'b> VolumeShortfall<'a, 'b> {
    fn is_short(&mut self, d: f64) -> bool {
        if d <= self.short_up_to {
            return true;
        }
        if d >= self.holds_from {
            return false;
        }
        let short = self.cache.get(d).volume < self.volume_m3;
        if short {
            self.short_up_to = d;
        } else {
            self.holds_from = d;
        }
        short
    }
}

fn bisection_cell(target: f64, d_min: f64, d_max: f64, tolerance: f64) -> (f64, f64) {
    let mut lo = d_min;
    let mut hi = d_max;
    while hi - lo > tolerance {
        let mid = 0.5 * (lo + hi);
        if mid < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (lo, hi)
}

fn probe_near(short: &mut VolumeShortfall<'_, '_>, hint: f64, d_min: f64, d_max: f64, tolerance: f64) {
    let clamped_hint = hint.max(d_min).min(d_max);
    let (lo, hi) = bisection_cell(clamped_hint, d_min, d_max, tolerance);
    let mut reach = (hi - lo).max(tolerance);
    if short.is_short(hi) {
        let mut probe = hi;
        while probe < d_max && short.is_short(probe) {
            probe = bisection_cell((probe + reach).min(d_max), d_min, d_max, tolerance).1;
            reach *= 2.0;
        }
        return;
    }
    let mut probe = lo;
    while probe > d_min && !short.is_short(probe) {
        let target = probe - reach;
        reach *= 2.0;
        if target <= d_min {
            return;
        }
        probe = bisection_cell(target, d_min, d_max, tolerance).1;
    }
}

fn volume_extent(
    cache: &mut ShapeCache<'_>,
    volume_m3: f64,
    d_min: f64,
    d_max: f64,
    tolerance: f64,
    hint: Option<f64>,
) -> f64 {
    if volume_m3 <= 0.0 {
        return d_min;
    }
    let mut short = VolumeShortfall {
        cache,
        volume_m3,
        short_up_to: d_min,
        holds_from: f64::INFINITY,
    };
    if let Some(h) = hint {
        if h.is_finite() {
            probe_near(&mut short, h, d_min, d_max, tolerance);
        }
    }
    if short.is_short(d_max) {
        return d_max;
    }
    let mut lo = d_min;
    let mut hi = d_max;
    while hi - lo > tolerance {
        let mid = 0.5 * (lo + hi);
        if short.is_short(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    hi
}

fn bracketed_root<F: FnMut(f64) -> f64>(
    mut f: F,
    mut lo: f64,
    mut hi: f64,
    mut f_lo: f64,
    mut f_hi: f64,
    tolerance: f64,
) -> (f64, f64) {
    let mut side = 0_i8;
    for _ in 0..ICE_ROOT_MAX_ITERATIONS {
        if hi - lo <= tolerance {
            break;
        }
        let mut x = if f_hi != f_lo {
            hi - f_hi * (hi - lo) / (f_hi - f_lo)
        } else {
            0.5 * (lo + hi)
        };
        if !(lo < x && x < hi) {
            x = 0.5 * (lo + hi);
        }
        let fx = f(x);
        if fx <= 0.0 {
            lo = x;
            f_lo = fx;
            if side == -1 {
                f_hi *= 0.5;
            }
            side = -1;
        } else {
            hi = x;
            f_hi = fx;
            if side == 1 {
                f_lo *= 0.5;
            }
            side = 1;
        }
    }
    (lo, hi)
}

pub fn solve_plastic_body(
    dom: &LocalDomain,
    opening_m3: f64,
    dt_myr: f64,
    sea_level_m: f64,
    step_m: f64,
    extent_hint: Option<f64>,
) -> (IceShape, f64, f64, f64, f64, bool) {
    let mut cache = ShapeCache::new(dom, sea_level_m);
    let d_min = dom.phi_min - 1.0;
    let d_max = (dom.phi_max - if dom.has_wall { 0.0 } else { 1.0 }).max(d_min);
    let d0 = volume_extent(&mut cache, opening_m3, d_min, d_max, step_m, extent_hint);

    let eval_f = |cache: &mut ShapeCache<'_>, d: f64| -> f64 {
        let s = cache.get(d);
        let (acc, sm, bm) = s.budget_rates(dom);
        s.volume - opening_m3 - dt_myr * (acc - sm - bm)
    };

    let f0 = eval_f(&mut cache, d0);
    let mut saturated = false;
    let mut lo;
    let mut hi;
    let mut f_lo;
    let mut f_hi;

    if f0 <= 0.0 {
        lo = d0;
        f_lo = f0;
        loop {
            hi = (lo + step_m).min(d_max);
            f_hi = eval_f(&mut cache, hi);
            if f_hi > 0.0 {
                break;
            }
            lo = hi;
            f_lo = f_hi;
            if hi >= d_max {
                saturated = true;
                break;
            }
        }
    } else {
        hi = d0;
        f_hi = f0;
        loop {
            lo = (hi - step_m).max(d_min);
            f_lo = eval_f(&mut cache, lo);
            if f_lo <= 0.0 {
                break;
            }
            hi = lo;
            f_hi = f_lo;
        }
    }

    let mut final_shape = if saturated {
        cache.get(d_max).clone()
    } else {
        (lo, hi) = bracketed_root(|x| eval_f(&mut cache, x), lo, hi, f_lo, f_hi, ICE_ROOT_TOLERANCE_M);
        let flo_final = eval_f(&mut cache, lo);
        let fhi_final = eval_f(&mut cache, hi);
        let lo_vol = cache.get(lo).volume;
        if lo_vol <= 0.0 || fhi_final.abs() < flo_final.abs() {
            cache.get(hi).clone()
        } else {
            cache.get(lo).clone()
        }
    };

    let (r_acc, r_sm, r_bm) = final_shape.budget_rates(dom);
    let acc = dt_myr * r_acc;
    let mut sm = dt_myr * r_sm;
    let mut bm = dt_myr * r_bm;
    let mut target = opening_m3 + acc - sm - bm;
    let mut excess = 0.0_f64;
    if saturated {
        excess = (target - final_shape.volume).max(0.0);
        target -= excess;
    }
    if target < 0.0 {
        let shortfall = -target;
        let sm_share = sm / (sm + bm).max(1e-30);
        sm -= shortfall * sm_share;
        bm -= shortfall * (1.0 - sm_share);
        target = 0.0;
    }
    if final_shape.volume > 0.0 {
        let factor = target / final_shape.volume;
        final_shape = final_shape.scaled(factor);
    }

    (final_shape, acc, sm, bm, excess, saturated)
}

#[allow(clippy::too_many_arguments)]
pub fn solve_uniform_body(
    nodes: &[usize],
    opening_m3: f64,
    area: &[f64],
    precip: &[f64],
    temp_c: &[f64],
    amp_c: &[f64],
    frictional: Option<&[f64]>,
    bed_frozen_snow: Option<&[f64]>,
    dt_myr: f64,
) -> (f64, f64, f64, f64, f64) {
    let mut total_area_acc = NeumaierSum::default();
    for &u in nodes {
        total_area_acc.add(area[u]);
    }
    let total_area = total_area_acc.finish();

    let budget = |h: f64| -> (f64, f64, f64) {
        let mut acc_sum = NeumaierSum::default();
        let mut sm_sum = NeumaierSum::default();
        let mut bm_sum = NeumaierSum::default();
        for &u in nodes {
            let fric = frictional.map_or(0.0, |v| v[u]);
            let bfs = bed_frozen_snow.map(|v| v[u]);
            let (a, s, b) = node_rates_scalar(h, temp_c[u], amp_c[u], precip[u], fric, bfs);
            acc_sum.add(area[u] * a);
            sm_sum.add(area[u] * s);
            bm_sum.add(area[u] * b);
        }
        (
            dt_myr * acc_sum.finish(),
            dt_myr * sm_sum.finish(),
            dt_myr * bm_sum.finish(),
        )
    };

    let g = |h: f64| -> f64 {
        let (acc, sm, bm) = budget(h);
        h * total_area - opening_m3 - (acc - sm - bm)
    };

    let g0 = g(0.0);
    if g0 >= 0.0 {
        let (acc, _, _) = budget(0.0);
        return (0.0, acc, opening_m3 + acc, 0.0, 0.0);
    }
    let h_max = GLACIER_MAX_DEPTH_M;
    let g_max = g(h_max);
    if g_max < 0.0 {
        let (acc, sm, bm) = budget(h_max);
        let excess = opening_m3 + acc - sm - bm - h_max * total_area;
        return (h_max, acc, sm, bm, excess);
    }
    let (lo, hi) = bracketed_root(g, 0.0, h_max, g0, g_max, 1e-3);
    let (mut acc, mut sm, mut bm) = budget(hi);
    let mut target = opening_m3 + acc - sm - bm;
    if target < 0.0 {
        (acc, sm, bm) = budget(lo);
        target = (opening_m3 + acc - sm - bm).max(0.0);
    }
    (target / total_area, acc, sm, bm, 0.0)
}

pub fn balance_flux(
    order: &[usize],
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    surface: &[f64],
    is_ice: &[bool],
    source: &[f64],
) -> Vec<f64> {
    let mut flux = vec![0.0_f64; surface.len()];
    for &i in order {
        let q = (flux[i] + source[i]).max(0.0);
        flux[i] = q;
        let start = indptr[i] as usize;
        let end = indptr[i + 1] as usize;
        let mut total = 0.0_f64;
        for e in start..end {
            let j = indices[e] as usize;
            if is_ice[j] && surface[j] < surface[i] {
                total += (surface[i] - surface[j]) / lengths[e];
            }
        }
        if total <= 0.0 || q <= 0.0 {
            continue;
        }
        for e in start..end {
            let j = indices[e] as usize;
            if is_ice[j] && surface[j] < surface[i] {
                flux[j] += q * ((surface[i] - surface[j]) / lengths[e]) / total;
            }
        }
    }
    flux
}

#[allow(clippy::too_many_arguments)]
pub fn build_local_domain(
    dom_nodes: &[usize],
    global_to_local: &[i32],
    indptr: &[u32],
    indices: &[u32],
    lengths: &[f64],
    phi: &[f64],
    bed: &[f64],
    c: &[f64],
    area: &[f64],
    sea_level_m: f64,
    is_ocean: &[bool],
    precip: &[f64],
    temp_c: &[f64],
    amp_c: &[f64],
    frictional: Option<&[f64]>,
    bed_frozen_snow: Option<&[f64]>,
    bare_acc_global: &[f64],
    bare_sm_global: &[f64],
) -> LocalDomain {
    let m = dom_nodes.len();
    let mut sub_indptr = Vec::with_capacity(m + 1);
    sub_indptr.push(0_u32);

    let mut total_edges = 0usize;
    for &u in dom_nodes {
        total_edges += (indptr[u + 1] - indptr[u]) as usize;
        sub_indptr.push(total_edges as u32);
    }

    // The gathers below are per-node independent and dominate the march's domain
    // build on large bodies (#1041): collect each node's edges and fields in
    // parallel, then concatenate in node order, so the result is identical to
    // the serial build. Small domains stay serial to avoid rayon overhead.
    let edge_part = |&u: &usize| {
        let phi_u = phi[u];
        let start = indptr[u] as usize;
        let end = indptr[u + 1] as usize;
        let deg = end - start;
        let mut p_nbr = Vec::with_capacity(deg);
        let mut p_len = Vec::with_capacity(deg);
        let mut p_phi = Vec::with_capacity(deg);
        let mut p_bed = Vec::with_capacity(deg);
        let mut p_cap = Vec::with_capacity(deg);
        let mut p_ocean = Vec::with_capacity(deg);
        let mut wall = false;
        for e in start..end {
            let v = indices[e] as usize;
            let lv = global_to_local[v];
            let len_e = lengths[e];
            p_nbr.push(lv);
            p_len.push(len_e);
            if lv < 0 {
                wall = true;
                p_phi.push(phi_u + len_e);
                p_cap.push(0.5);
            } else {
                p_phi.push(phi[v]);
                p_cap.push(1.0);
            }
            p_bed.push(bed[v]);
            p_ocean.push(is_ocean[v]);
        }
        (p_nbr, p_len, p_phi, p_bed, p_cap, p_ocean, wall)
    };
    let node_part = |&u: &usize| {
        (
            phi[u],
            bed[u],
            c[u],
            area[u],
            floating_water_depth_scalar(bed[u], sea_level_m, is_ocean[u]),
            is_ocean[u],
            precip[u],
            temp_c[u],
            amp_c[u],
            frictional.map(|f| f[u]),
            bed_frozen_snow.map(|b| b[u]),
            bare_acc_global[u],
            bare_sm_global[u],
        )
    };
    // Chunked parallel gather: one task per block of nodes (each block builds
    // its parts serially), so task and allocation overhead stay amortized.
    let (edge_parts, node_parts) = if m >= 4096 {
        let block = 2048;
        let blocks: Vec<(usize, usize)> = (0..m)
            .step_by(block)
            .map(|s| (s, (s + block).min(m)))
            .collect();
        rayon::join(
            || {
                blocks
                    .par_iter()
                    .map(|&(s, e)| {
                        let cap = (sub_indptr[e] - sub_indptr[s]) as usize;
                        let mut b_nbr = Vec::with_capacity(cap);
                        let mut b_len = Vec::with_capacity(cap);
                        let mut b_phi = Vec::with_capacity(cap);
                        let mut b_bed = Vec::with_capacity(cap);
                        let mut b_cap = Vec::with_capacity(cap);
                        let mut b_ocean = Vec::with_capacity(cap);
                        let mut wall = false;
                        for &u in &dom_nodes[s..e] {
                            let phi_u = phi[u];
                            let start = indptr[u] as usize;
                            let end = indptr[u + 1] as usize;
                            for e in start..end {
                                let v = indices[e] as usize;
                                let lv = global_to_local[v];
                                let len_e = lengths[e];
                                b_nbr.push(lv);
                                b_len.push(len_e);
                                if lv < 0 {
                                    wall = true;
                                    b_phi.push(phi_u + len_e);
                                    b_cap.push(0.5);
                                } else {
                                    b_phi.push(phi[v]);
                                    b_cap.push(1.0);
                                }
                                b_bed.push(bed[v]);
                                b_ocean.push(is_ocean[v]);
                            }
                        }
                        (b_nbr, b_len, b_phi, b_bed, b_cap, b_ocean, wall)
                    })
                    .collect::<Vec<_>>()
            },
            || {
                blocks
                    .par_iter()
                    .map(|&(s, e)| dom_nodes[s..e].iter().map(node_part).collect::<Vec<_>>())
                    .collect::<Vec<_>>()
                    .into_iter()
                    .flatten()
                    .collect::<Vec<_>>()
            },
        )
    } else {
        (
            dom_nodes.iter().map(edge_part).collect::<Vec<_>>(),
            dom_nodes.iter().map(node_part).collect::<Vec<_>>(),
        )
    };

    let mut nbr = Vec::with_capacity(total_edges);
    let mut length = Vec::with_capacity(total_edges);
    let mut nbr_phi = Vec::with_capacity(total_edges);
    let mut nbr_bed = Vec::with_capacity(total_edges);
    let mut nbr_cap = Vec::with_capacity(total_edges);
    let mut nbr_ocean = Vec::with_capacity(total_edges);
    let mut has_wall = false;
    for (p_nbr, p_len, p_phi, p_bed, p_cap, p_ocean, wall) in edge_parts {
        nbr.extend(p_nbr);
        length.extend(p_len);
        nbr_phi.extend(p_phi);
        nbr_bed.extend(p_bed);
        nbr_cap.extend(p_cap);
        nbr_ocean.extend(p_ocean);
        has_wall |= wall;
    }

    let mut l_phi = Vec::with_capacity(m);
    let mut l_bed = Vec::with_capacity(m);
    let mut l_c = Vec::with_capacity(m);
    let mut l_area = Vec::with_capacity(m);
    let mut l_wd = Vec::with_capacity(m);
    let mut l_ocean = Vec::with_capacity(m);
    let mut l_precip = Vec::with_capacity(m);
    let mut l_temp = Vec::with_capacity(m);
    let mut l_amp = Vec::with_capacity(m);
    let mut l_fric = frictional.map(|_| Vec::with_capacity(m));
    let mut l_bfs = bed_frozen_snow.map(|_| Vec::with_capacity(m));
    let mut bare_acc = Vec::with_capacity(m);
    let mut bare_sm = Vec::with_capacity(m);

    let mut phi_min = f64::INFINITY;
    let mut phi_max = f64::NEG_INFINITY;

    for (pu, bd, cc, ar, wd, oc, pr, te, am, fr, bf, ba, bs) in node_parts {
        if pu < phi_min {
            phi_min = pu;
        }
        if pu > phi_max {
            phi_max = pu;
        }
        l_phi.push(pu);
        l_bed.push(bd);
        l_c.push(cc);
        l_area.push(ar);
        l_wd.push(wd);
        l_ocean.push(oc);
        l_precip.push(pr);
        l_temp.push(te);
        l_amp.push(am);
        if let (Some(ref mut vf), Some(f)) = (&mut l_fric, fr) {
            vf.push(f);
        }
        if let (Some(ref mut vb), Some(b)) = (&mut l_bfs, bf) {
            vb.push(b);
        }
        bare_acc.push(ba);
        bare_sm.push(bs);
    }

    LocalDomain {
        nodes: dom_nodes.to_vec(),
        indptr: sub_indptr,
        nbr,
        length,
        nbr_phi,
        nbr_bed,
        nbr_cap,
        nbr_ocean,
        phi: l_phi,
        bed: l_bed,
        c: l_c,
        area: l_area,
        water_depth: l_wd,
        ocean: l_ocean,
        precip: l_precip,
        temp_c: l_temp,
        amp_c: l_amp,
        frictional: l_fric,
        bed_frozen_snow: l_bfs,
        bare_acc,
        bare_sm,
        has_wall,
        phi_min,
        phi_max,
    }
}
