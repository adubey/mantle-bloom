use numpy::{PyArray1, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rayon::prelude::*;

use crate::climate::{
    air_temperature_at_elevation_scalar, bare_rates_scalar, clip_min, compensated_dot,
    compensated_dot_indexed, is_grounded_scalar, node_rates_scalar, snow_fraction_scalar,
    warm_based_scalar, yield_stress_scalar, NeumaierSum, FLOTATION_RATIO, GLACIER_MAX_DEPTH_M,
    GRAVITY_M_S2, ICE_BODY_MIN_NODES, ICE_SLIDING_FRACTION, ICE_VELOCITY_MIN_THICKNESS_M, RHO_ICE,
};
use crate::distance::signed_distance;
use crate::graph::{connected_components_masked, MarchGraph};
use crate::plastic::{
    balance_flux, build_local_domain, solve_plastic_body, solve_uniform_body, IceShape, LocalDomain,
};

struct BodyOut {
    nodes: Vec<i64>,
    domain_nodes: Vec<i64>,
    opening_m3: f64,
    accumulation_m3: f64,
    surface_melt_m3: f64,
    basal_melt_m3: f64,
    calved_m3: f64,
    overflow_m3: f64,
    closing_m3: f64,
    capped_m3: f64,
    margin_m: f64,
    saturated: bool,
    dom_thickness: Vec<f64>,
    dom_covered: Vec<f64>,
    dom_depth: Vec<f64>,
    dom_surface: Vec<f64>,
    dom_melt_m3: Vec<f64>,
    dom_ice_release_m3: Vec<f64>,
    dom_surface_melt_m3: Vec<f64>,
    dom_bare_snowmelt_m3: Vec<f64>,
    dom_retained_snow: Vec<f64>,
    excess_calved: Vec<(usize, f64)>,
    excess_melt: Vec<(usize, f64)>,
    extra_bare_snow_melted: f64,
}

fn body_margin_hint(margin_hint_m: Option<&[f64]>, body_nodes: &[usize]) -> Option<f64> {
    let hints = margin_hint_m?;
    let mut known: Vec<f64> = body_nodes
        .iter()
        .map(|&u| hints[u])
        .filter(|v| v.is_finite())
        .collect();
    if known.is_empty() {
        return None;
    }
    known.sort_by(|a, b| a.total_cmp(b));
    let mid = known.len() / 2;
    if known.len() % 2 == 1 {
        Some(known[mid])
    } else {
        Some(0.5 * (known[mid - 1] + known[mid]))
    }
}

fn touches_ocean(dom: &LocalDomain, active: &[bool]) -> bool {
    for (i, &act) in active.iter().enumerate() {
        if !act {
            continue;
        }
        let start = dom.indptr[i] as usize;
        let end = dom.indptr[i + 1] as usize;
        for e in start..end {
            if dom.nbr_ocean[e] {
                return true;
            }
        }
    }
    false
}

fn shed_excess_items(
    excess: f64,
    margin_nodes: &[usize],
    marine: bool,
    area: &[f64],
) -> (Vec<(usize, f64)>, Vec<(usize, f64)>, f64, f64) {
    if excess <= 0.0 {
        return (Vec::new(), Vec::new(), 0.0, 0.0);
    }
    let mut total_area_acc = NeumaierSum::default();
    for &u in margin_nodes {
        total_area_acc.add(area[u]);
    }
    let total_area = total_area_acc.finish();
    let mut items = Vec::with_capacity(margin_nodes.len());
    for &u in margin_nodes {
        items.push((u, excess * (area[u] / total_area)));
    }
    if marine {
        (items, Vec::new(), excess, 0.0)
    } else {
        (Vec::new(), items, 0.0, excess)
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_one_uniform_body(
    body_nodes: &[usize],
    dom_nodes: &[usize],
    opening: f64,
    area: &[f64],
    bed: &[f64],
    is_ocean: &[bool],
    precip: &[f64],
    temp_c: &[f64],
    amp_c: &[f64],
    frictional: Option<&[f64]>,
    bed_frozen_snow: Option<&[f64]>,
    bare_acc: &[f64],
    dt_myr: f64,
    is_member: &[bool],
) -> BodyOut {
    let (h, acc, sm, bm, excess) = solve_uniform_body(
        body_nodes,
        opening,
        area,
        precip,
        temp_c,
        amp_c,
        frictional,
        bed_frozen_snow,
        dt_myr,
    );
    let mut body_area_acc = NeumaierSum::default();
    for &u in body_nodes {
        body_area_acc.add(area[u]);
    }
    let body_area_sum = body_area_acc.finish();
    let marine = body_nodes.iter().any(|&u| is_ocean[u]);
    let (excess_calved, excess_melt, calved_m3, overflow_m3) =
        shed_excess_items(excess, body_nodes, marine, area);

    let m = dom_nodes.len();
    let mut dom_thickness = vec![0.0; m];
    let mut dom_covered = vec![0.0; m];
    let mut dom_depth = vec![0.0; m];
    let mut dom_surface = vec![f64::NAN; m];
    let mut dom_melt_m3 = vec![0.0; m];
    let mut dom_ice_release_m3 = vec![0.0; m];
    let mut dom_surface_melt_m3 = vec![0.0; m];
    let mut dom_bare_snowmelt_m3 = vec![0.0; m];
    let mut dom_retained_snow = vec![0.0; m];
    let mut extra_bare_snow_melted = NeumaierSum::default();

    for (idx, &u) in dom_nodes.iter().enumerate() {
        if is_member[u] {
            dom_thickness[idx] = h;
            dom_covered[idx] = if h > 0.0 { 1.0 } else { 0.0 };
            dom_depth[idx] = h;
            dom_retained_snow[idx] = 1.0;
            dom_surface[idx] = if h > 0.0 { bed[u] + h } else { f64::NAN };
            dom_ice_release_m3[idx] = (sm + bm) * area[u] / body_area_sum;
            dom_surface_melt_m3[idx] = sm * area[u] / body_area_sum;
        } else {
            let w = if is_ocean[u] {
                0.0
            } else {
                bare_acc[u] * area[u] * dt_myr
            };
            dom_bare_snowmelt_m3[idx] = w;
            extra_bare_snow_melted.add(w);
        }
        dom_melt_m3[idx] = dom_ice_release_m3[idx] + dom_bare_snowmelt_m3[idx];
    }

    BodyOut {
        nodes: body_nodes.iter().map(|&u| u as i64).collect(),
        domain_nodes: dom_nodes.iter().map(|&u| u as i64).collect(),
        opening_m3: opening,
        accumulation_m3: acc,
        surface_melt_m3: sm,
        basal_melt_m3: bm,
        calved_m3,
        overflow_m3,
        closing_m3: 0.0,
        capped_m3: 0.0,
        margin_m: 0.0,
        saturated: excess > 0.0,
        dom_thickness,
        dom_covered,
        dom_depth,
        dom_surface,
        dom_melt_m3,
        dom_ice_release_m3,
        dom_surface_melt_m3,
        dom_bare_snowmelt_m3,
        dom_retained_snow,
        excess_calved,
        excess_melt,
        extra_bare_snow_melted: extra_bare_snow_melted.finish(),
    }
}

#[allow(clippy::too_many_arguments)]
fn solve_one_plastic_body(
    body_nodes: &[usize],
    dom: &LocalDomain,
    opening: f64,
    area: &[f64],
    dt_myr: f64,
    sea_level_m: f64,
    step_m: f64,
    hint: Option<f64>,
) -> BodyOut {
    let (shape, acc, sm, bm, excess, saturated): (IceShape, f64, f64, f64, f64, bool) =
        solve_plastic_body(dom, opening, dt_myr, sea_level_m, step_m, hint);
    let m = dom.nodes.len();
    let dom_thickness = shape.thickness;
    let mut dom_covered = Vec::with_capacity(m);
    let mut dom_surface = Vec::with_capacity(m);
    let mut where_surface = Vec::with_capacity(m);
    let mut where_basal = Vec::with_capacity(m);
    let mut where_surface_acc = NeumaierSum::default();
    let mut where_basal_acc = NeumaierSum::default();

    for i in 0..m {
        let cov = shape.cover[i].min(1.0);
        dom_covered.push(cov);
        dom_surface.push(if shape.active[i] {
            dom.bed[i] + dom_thickness[i]
        } else {
            f64::NAN
        });
        let feeds = (dom.bare_acc[i] - dom.bare_sm[i]) > 0.0;
        let ws = shape.cover[i] * dom.area[i] * shape.rate_sm[i]
            + (1.0 - cov) * dom.area[i] * if feeds { dom.bare_sm[i] } else { 0.0 };
        let wb = shape.cover[i] * dom.area[i] * shape.rate_bm[i];
        where_surface.push(ws);
        where_basal.push(wb);
        where_surface_acc.add(ws);
        where_basal_acc.add(wb);
    }
    let where_surface_sum = where_surface_acc.finish();
    let where_basal_sum = where_basal_acc.finish();

    let mut dom_melt_m3 = vec![0.0_f64; m];
    let mut dom_ice_release_m3 = vec![0.0_f64; m];
    let mut dom_surface_melt_m3 = vec![0.0_f64; m];
    let mut dom_bare_snowmelt_m3 = vec![0.0_f64; m];
    if where_surface_sum > 0.0 {
        for i in 0..m {
            let surface_melt = sm * where_surface[i] / where_surface_sum;
            dom_surface_melt_m3[i] += surface_melt;
            dom_ice_release_m3[i] += surface_melt;
        }
    } else if sm + bm > 0.0 {
        let mut body_area_acc = NeumaierSum::default();
        for &u in body_nodes {
            body_area_acc.add(area[u]);
        }
        let body_area_sum = body_area_acc.finish();
        for (i, &u) in dom.nodes.iter().enumerate() {
            if body_nodes.binary_search(&u).is_ok() {
                let surface_melt = sm * area[u] / body_area_sum;
                dom_surface_melt_m3[i] += surface_melt;
                dom_ice_release_m3[i] += surface_melt;
            }
        }
    }
    if where_basal_sum > 0.0 {
        for i in 0..m {
            dom_ice_release_m3[i] += bm * where_basal[i] / where_basal_sum;
        }
    } else if bm > 0.0 {
        let mut body_area_acc = NeumaierSum::default();
        for &u in body_nodes {
            body_area_acc.add(area[u]);
        }
        let body_area_sum = body_area_acc.finish();
        for (i, &u) in dom.nodes.iter().enumerate() {
            if body_nodes.binary_search(&u).is_ok() {
                dom_ice_release_m3[i] += bm * area[u] / body_area_sum;
            }
        }
    }

    let mut dom_retained_snow = Vec::with_capacity(m);
    let mut extra_bare_snow_melted = NeumaierSum::default();
    for i in 0..m {
        let cov = dom_covered[i];
        let feeds = (dom.bare_acc[i] - dom.bare_sm[i]) > 0.0;
        let bare_snow = if dom.ocean[i] {
            0.0
        } else {
            (1.0 - cov) * dom.area[i] * if feeds { 0.0 } else { dom.bare_acc[i] } * dt_myr
        };
        dom_bare_snowmelt_m3[i] += bare_snow;
        extra_bare_snow_melted.add(bare_snow);
        dom_retained_snow.push(if feeds { 1.0 } else { cov });
        dom_melt_m3[i] = dom_ice_release_m3[i] + dom_bare_snowmelt_m3[i];
    }

    let mut frontier_nodes = Vec::new();
    let mut active_nodes = Vec::new();
    let mut marine = touches_ocean(dom, &shape.active);
    for i in 0..m {
        if shape.active[i] {
            let u = dom.nodes[i];
            active_nodes.push(u);
            if dom_covered[i] < 1.0 {
                frontier_nodes.push(u);
            }
            if dom.water_depth[i] > 0.0 {
                marine = true;
            }
        }
    }
    let margin = if !frontier_nodes.is_empty() {
        &frontier_nodes[..]
    } else if !active_nodes.is_empty() {
        &active_nodes[..]
    } else {
        body_nodes
    };
    let (excess_calved, excess_melt, calved_m3, overflow_m3) =
        shed_excess_items(excess, margin, marine, area);

    BodyOut {
        nodes: body_nodes.iter().map(|&u| u as i64).collect(),
        domain_nodes: dom.nodes.iter().map(|&u| u as i64).collect(),
        opening_m3: opening,
        accumulation_m3: acc,
        surface_melt_m3: sm,
        basal_melt_m3: bm,
        calved_m3,
        overflow_m3,
        closing_m3: 0.0,
        capped_m3: 0.0,
        margin_m: shape.d,
        saturated,
        dom_thickness,
        dom_covered,
        dom_depth: shape.depth,
        dom_surface,
        dom_melt_m3,
        dom_ice_release_m3,
        dom_surface_melt_m3,
        dom_bare_snowmelt_m3,
        dom_retained_snow,
        excess_calved,
        excess_melt,
        extra_bare_snow_melted: extra_bare_snow_melted.finish(),
    }
}

#[pyfunction]
#[allow(clippy::too_many_arguments)]
#[pyo3(signature = (
    rust_graph,
    area_m2,
    elevation,
    prior_deflection_m,
    sea_level_m,
    prev_depth,
    precipitation_mm,
    temperature_c,
    seasonal_amplitude_c,
    mobile_cover_m,
    is_ocean,
    years,
    frictional_heat_w_m2=None,
    dry_snow_mm=None,
    bed_frozen_snow_mm=None,
    margin_hint_m=None,
    dijkstra_limit_m=None,
))]
pub fn step_ice<'py>(
    py: Python<'py>,
    rust_graph: &MarchGraph,
    area_m2: PyReadonlyArray1<'py, f64>,
    elevation: PyReadonlyArray1<'py, f64>,
    prior_deflection_m: PyReadonlyArray1<'py, f64>,
    sea_level_m: f64,
    prev_depth: PyReadonlyArray1<'py, f64>,
    precipitation_mm: PyReadonlyArray1<'py, f64>,
    temperature_c: PyReadonlyArray1<'py, f64>,
    seasonal_amplitude_c: PyReadonlyArray1<'py, f64>,
    mobile_cover_m: PyReadonlyArray1<'py, f64>,
    is_ocean: PyReadonlyArray1<'py, bool>,
    years: f64,
    frictional_heat_w_m2: Option<PyReadonlyArray1<'py, f64>>,
    dry_snow_mm: Option<PyReadonlyArray1<'py, f64>>,
    bed_frozen_snow_mm: Option<PyReadonlyArray1<'py, f64>>,
    margin_hint_m: Option<PyReadonlyArray1<'py, f64>>,
    dijkstra_limit_m: Option<f64>,
) -> PyResult<Py<PyAny>> {
    let area = area_m2.as_slice()?;
    let elev = elevation.as_slice()?;
    let defl = prior_deflection_m.as_slice()?;
    let prev = prev_depth.as_slice()?;
    let precip = precipitation_mm.as_slice()?;
    let temp = temperature_c.as_slice()?;
    let amp = seasonal_amplitude_c.as_slice()?;
    let mobile = mobile_cover_m.as_slice()?;
    let ocean = is_ocean.as_slice()?;
    let fric = match &frictional_heat_w_m2 {
        Some(a) => Some(a.as_slice()?),
        None => None,
    };
    let dry_snow = match &dry_snow_mm {
        Some(a) => Some(a.as_slice()?),
        None => None,
    };
    let bfs = match &bed_frozen_snow_mm {
        Some(a) => Some(a.as_slice()?),
        None => None,
    };
    let hints = match &margin_hint_m {
        Some(a) => Some(a.as_slice()?),
        None => None,
    };

    let n = area.len();
    let dt_myr = years / 1_000_000.0;

    // Opt-in per-section timing for march profiling (#1041); off by default and
    // read once per call so a probe run can toggle it between steps.
    let probe = std::env::var_os("MANTLE_BLOOM_ICE_PROFILE").is_some_and(|v| v == "1");
    let result = py.detach(|| {
        let t0 = std::time::Instant::now();
        let mut bed = vec![0.0_f64; n];
        let mut opening_depth = vec![0.0_f64; n];
        let mut h0 = vec![0.0_f64; n];
        let mut eligible = vec![false; n];
        let mut bare_acc = vec![0.0_f64; n];
        let mut bare_sm = vec![0.0_f64; n];
        let mut ref_balance = vec![0.0_f64; n];
        let mut positive = vec![false; n];
        let mut members = vec![false; n];
        let mut c_coef = vec![0.0_f64; n];
        let mut calved = vec![0.0_f64; n];

        let setup_node = |i: usize| {
            let b = elev[i] - defl[i];
            let od = clip_min(prev[i], 0.0);
            let wd = if ocean[i] {
                clip_min(sea_level_m - b, 0.0)
            } else {
                0.0
            };
            let elig = FLOTATION_RATIO * wd <= GLACIER_MAX_DEPTH_M;
            let stray = od > 0.0 && (!elig || !is_grounded_scalar(od, b, sea_level_m, ocean[i]));
            let (calv, h_val) = if stray {
                (od * area[i], 0.0)
            } else {
                (0.0, od)
            };
            let fric_i = fric.map_or(0.0, |v| v[i]);
            let bfs_i = bfs.map(|v| v[i]);
            let (ba, bs) = bare_rates_scalar(temp[i], amp[i], precip[i], bfs_i);
            let rb_net = if h_val == 0.0 {
                ba - bs
            } else {
                let (ra, rs, rb) =
                    node_rates_scalar(h_val, temp[i], amp[i], precip[i], fric_i, bfs_i);
                ra - rs - rb
            };
            let pos = elig && (rb_net > 0.0) && (!ocean[i] || h_val > 0.0);
            let mem = elig && (h_val > 0.0 || pos);
            let surf_t = air_temperature_at_elevation_scalar(temp[i], 0.0, h_val);
            let warm0 = warm_based_scalar(h_val, surf_t, fric_i);
            let c_val = yield_stress_scalar(mobile[i], warm0) / (RHO_ICE * GRAVITY_M_S2);
            (b, od, h_val, elig, calv, ba, bs, rb_net, pos, mem, c_val)
        };

        let any_members = if n >= 4096 {
            let chunk = 4096;
            (
                (
                    bed.par_chunks_mut(chunk),
                    opening_depth.par_chunks_mut(chunk),
                    h0.par_chunks_mut(chunk),
                    eligible.par_chunks_mut(chunk),
                    calved.par_chunks_mut(chunk),
                    bare_acc.par_chunks_mut(chunk),
                ),
                (
                    bare_sm.par_chunks_mut(chunk),
                    ref_balance.par_chunks_mut(chunk),
                    positive.par_chunks_mut(chunk),
                    members.par_chunks_mut(chunk),
                    c_coef.par_chunks_mut(chunk),
                ),
            )
                .into_par_iter()
                .enumerate()
                .map(
                    |(
                        chunk_idx,
                        ((b_c, od_c, h0_c, el_c, cv_c, ba_c), (bs_c, rb_c, pos_c, mem_c, cc_c)),
                    )| {
                        let base = chunk_idx * chunk;
                        let mut local_any = false;
                        for j in 0..b_c.len() {
                            let (b, od, hv, el, cv, ba, bs, rbn, pos, mem, cc) =
                                setup_node(base + j);
                            b_c[j] = b;
                            od_c[j] = od;
                            h0_c[j] = hv;
                            el_c[j] = el;
                            cv_c[j] = cv;
                            ba_c[j] = ba;
                            bs_c[j] = bs;
                            rb_c[j] = rbn;
                            pos_c[j] = pos;
                            mem_c[j] = mem;
                            cc_c[j] = cc;
                            local_any |= mem;
                        }
                        local_any
                    },
                )
                .reduce(|| false, |a, b| a || b)
        } else {
            let mut local_any = false;
            for i in 0..n {
                let (b, od, hv, el, cv, ba, bs, rbn, pos, mem, cc) = setup_node(i);
                bed[i] = b;
                opening_depth[i] = od;
                h0[i] = hv;
                eligible[i] = el;
                calved[i] = cv;
                bare_acc[i] = ba;
                bare_sm[i] = bs;
                ref_balance[i] = rbn;
                positive[i] = pos;
                members[i] = mem;
                c_coef[i] = cc;
                local_any |= mem;
            }
            local_any
        };

        let t_setup = t0.elapsed();
        let opening_m3 = compensated_dot(&opening_depth, area);

        let mut depth = vec![0.0_f64; n];
        let mut melt_m3 = vec![0.0_f64; n];
        let mut ice_release_m3 = vec![0.0_f64; n];
        let mut surface_melt_by_node_m3 = vec![0.0_f64; n];
        let mut bare_snowmelt_m3 = vec![0.0_f64; n];
        let mut thickness = vec![0.0_f64; n];
        let mut covered = vec![0.0_f64; n];
        let mut surface = vec![f64::NAN; n];
        let mut retained_snow = vec![0.0_f64; n];
        let mut in_domain = vec![false; n];
        let mut body_outs: Vec<BodyOut> = Vec::new();
        let mut snow_melted_acc = NeumaierSum::default();

        if any_members {
            let (comp, num_bodies) =
                connected_components_masked(n, &rust_graph.indptr, &rust_graph.indices, &members);

            let mut body_nodes_all: Vec<Vec<usize>> = vec![Vec::new(); num_bodies];
            for i in 0..n {
                let c = comp[i];
                if c >= 0 {
                    body_nodes_all[c as usize].push(i);
                }
            }

            let mut seeds = positive.clone();

            for b_nodes in &body_nodes_all {
                let mut any_pos = false;
                let mut all_pos = true;
                for &u in b_nodes {
                    if positive[u] {
                        any_pos = true;
                    } else {
                        all_pos = false;
                    }
                }
                // Only an all-positive body reads whether its nodes touch a
                // non-positive neighbour, so that scan runs lazily here over
                // the body's own nodes instead of the whole cloud (PR #1062
                // review); same per-node edge sum, same result.
                let any_touch = all_pos
                    && b_nodes.iter().any(|&u| {
                        let start = rust_graph.indptr[u] as usize;
                        let end = rust_graph.indptr[u + 1] as usize;
                        let mut sum_out = 0.0_f64;
                        for e in start..end {
                            if !positive[rust_graph.indices[e] as usize] {
                                sum_out += rust_graph.lengths[e];
                            }
                        }
                        sum_out > 0.0
                    });
                if !any_pos {
                    let mut best_u = b_nodes[0];
                    let mut best_val = ref_balance[best_u];
                    for &u in &b_nodes[1..] {
                        if ref_balance[u] > best_val {
                            best_val = ref_balance[u];
                            best_u = u;
                        }
                    }
                    seeds[best_u] = true;
                } else if all_pos && !any_touch {
                    let mut best_u = b_nodes[0];
                    for &u in &b_nodes[1..] {
                        let cmp_rb = ref_balance[u].total_cmp(&ref_balance[best_u]);
                        if cmp_rb == std::cmp::Ordering::Less
                            || (cmp_rb == std::cmp::Ordering::Equal
                                && bed[u].total_cmp(&bed[best_u]) == std::cmp::Ordering::Less)
                        {
                            best_u = u;
                        }
                    }
                    seeds[best_u] = false;
                }
            }

            let t_comp = t0.elapsed();
            let (phi, sources) = signed_distance(
                n,
                &rust_graph.indptr,
                &rust_graph.indices,
                &rust_graph.lengths,
                &eligible,
                &seeds,
                &members,
                dijkstra_limit_m,
            );

            let t_dist = t0.elapsed();
            let mut dom_nodes_all: Vec<Vec<usize>> = vec![Vec::new(); num_bodies];
            for i in 0..n {
                let lbl = if members[i] {
                    comp[i]
                } else if eligible[i] && sources[i] >= 0 {
                    comp[sources[i] as usize]
                } else {
                    -1
                };
                if lbl >= 0 {
                    in_domain[i] = true;
                    dom_nodes_all[lbl as usize].push(i);
                }
            }

            let mut global_to_local = vec![-1_i32; n];
            let mut plastic_specs: Vec<(usize, LocalDomain, f64, Option<f64>)> = Vec::new();
            let mut uniform_indices: Vec<usize> = Vec::new();
            let mut total_dom_nodes = 0usize;

            for b in 0..num_bodies {
                let d_nodes = &dom_nodes_all[b];
                let b_nodes = &body_nodes_all[b];
                let opening = compensated_dot_indexed(b_nodes, &h0, area);
                if d_nodes.len() < ICE_BODY_MIN_NODES {
                    uniform_indices.push(b);
                } else {
                    for (local_i, &u) in d_nodes.iter().enumerate() {
                        global_to_local[u] = local_i as i32;
                    }
                    let dom = build_local_domain(
                        d_nodes,
                        &global_to_local,
                        &rust_graph.indptr,
                        &rust_graph.indices,
                        &rust_graph.lengths,
                        &phi,
                        &bed,
                        &c_coef,
                        area,
                        sea_level_m,
                        ocean,
                        precip,
                        temp,
                        amp,
                        fric,
                        bfs,
                        &bare_acc,
                        &bare_sm,
                    );
                    for &u in d_nodes {
                        global_to_local[u] = -1;
                    }
                    total_dom_nodes += d_nodes.len();
                    let hint = body_margin_hint(hints, b_nodes);
                    plastic_specs.push((b, dom, opening, hint));
                }
            }

            let t_dom = t0.elapsed();
            let mut ordered_outs: Vec<Option<BodyOut>> = (0..num_bodies).map(|_| None).collect();

            for &b in &uniform_indices {
                let b_nodes = &body_nodes_all[b];
                let d_nodes = &dom_nodes_all[b];
                let opening = compensated_dot_indexed(b_nodes, &h0, area);
                let out_b = solve_one_uniform_body(
                    b_nodes, d_nodes, opening, area, &bed, ocean, precip, temp, amp, fric, bfs,
                    &bare_acc, dt_myr, &members,
                );
                ordered_outs[b] = Some(out_b);
            }

            let plastic_results: Vec<(usize, BodyOut)> =
                if plastic_specs.len() >= 2 && total_dom_nodes >= 256 {
                    plastic_specs
                        .into_par_iter()
                        .map(|(b, dom, opening, hint)| {
                            let out_b = solve_one_plastic_body(
                                &body_nodes_all[b],
                                &dom,
                                opening,
                                area,
                                dt_myr,
                                sea_level_m,
                                rust_graph.step_m,
                                hint,
                            );
                            (b, out_b)
                        })
                        .collect()
                } else {
                    plastic_specs
                        .into_iter()
                        .map(|(b, dom, opening, hint)| {
                            let out_b = solve_one_plastic_body(
                                &body_nodes_all[b],
                                &dom,
                                opening,
                                area,
                                dt_myr,
                                sea_level_m,
                                rust_graph.step_m,
                                hint,
                            );
                            (b, out_b)
                        })
                        .collect()
                };

            for (b, out_b) in plastic_results {
                ordered_outs[b] = Some(out_b);
            }

            for opt_b in ordered_outs {
                let mut out_b = opt_b.unwrap();
                for (i, &u_i64) in out_b.domain_nodes.iter().enumerate() {
                    let u = u_i64 as usize;
                    thickness[u] = out_b.dom_thickness[i];
                    covered[u] = out_b.dom_covered[i];
                    depth[u] = out_b.dom_depth[i];
                    surface[u] = out_b.dom_surface[i];
                    melt_m3[u] += out_b.dom_melt_m3[i];
                    ice_release_m3[u] += out_b.dom_ice_release_m3[i];
                    surface_melt_by_node_m3[u] += out_b.dom_surface_melt_m3[i];
                    bare_snowmelt_m3[u] += out_b.dom_bare_snowmelt_m3[i];
                    retained_snow[u] = out_b.dom_retained_snow[i];
                }
                for (u, v) in &out_b.excess_calved {
                    calved[*u] += *v;
                }
                for (u, v) in &out_b.excess_melt {
                    melt_m3[*u] += *v;
                    ice_release_m3[*u] += *v;
                }
                snow_melted_acc.add(out_b.extra_bare_snow_melted);

                let mut capped_here = 0.0_f64;
                let mut calved_guard = 0.0_f64;
                let mut closing_acc = NeumaierSum::default();
                for &u_i64 in &out_b.domain_nodes {
                    let u = u_i64 as usize;
                    if depth[u] > GLACIER_MAX_DEPTH_M {
                        let over = depth[u] - GLACIER_MAX_DEPTH_M;
                        let kept = GLACIER_MAX_DEPTH_M;
                        if depth[u] > 0.0 {
                            thickness[u] *= kept / depth[u];
                        }
                        depth[u] = kept;
                        melt_m3[u] += over * area[u];
                        ice_release_m3[u] += over * area[u];
                        capped_here += over * area[u];
                    }
                    if depth[u] > 0.0
                        && !is_grounded_scalar(depth[u], bed[u], sea_level_m, ocean[u])
                    {
                        let dv = depth[u] * area[u];
                        calved[u] += dv;
                        calved_guard += dv;
                        depth[u] = 0.0;
                        thickness[u] = 0.0;
                        covered[u] = 0.0;
                        surface[u] = f64::NAN;
                    }
                    closing_acc.add(depth[u] * area[u]);
                }
                out_b.capped_m3 += capped_here;
                out_b.calved_m3 += calved_guard;
                out_b.closing_m3 = closing_acc.finish();
                body_outs.push(out_b);
            }
            if probe {
                eprintln!(
                    "MANTLE_BLOOM_ICE_PROFILE setup_us={} comp_us={} dist_us={} dom_us={} solve_scatter_us={} nbodies={} dom_nodes={}",
                    t_setup.as_micros(),
                    (t_comp - t_setup).as_micros(),
                    (t_dist - t_comp).as_micros(),
                    (t_dom - t_dist).as_micros(),
                    (t0.elapsed() - t_dom).as_micros(),
                    num_bodies,
                    total_dom_nodes
                );
            }
        }

        for i in 0..n {
            if !in_domain[i] && !ocean[i] {
                let water = bare_acc[i] * area[i] * dt_myr;
                melt_m3[i] += water;
                bare_snowmelt_m3[i] += water;
                snow_melted_acc.add(water);
            }
        }
        if probe {
            eprintln!(
                "MANTLE_BLOOM_ICE_PROFILE total_us={} bodies_done",
                t0.elapsed().as_micros()
            );
        }
        let snow_melted_m3 = snow_melted_acc.finish();

        let mut accumulation_m3 = 0.0_f64;
        let mut surface_melt_m3 = 0.0_f64;
        let mut basal_melt_m3 = 0.0_f64;
        let mut overflow_m3 = 0.0_f64;
        let mut capped_m3 = 0.0_f64;
        for b in &body_outs {
            accumulation_m3 += b.accumulation_m3;
            surface_melt_m3 += b.surface_melt_m3;
            basal_melt_m3 += b.basal_melt_m3;
            overflow_m3 += b.overflow_m3;
            capped_m3 += b.capped_m3;
        }
        let mut calved_acc = NeumaierSum::default();
        for &cv in &calved {
            calved_acc.add(cv);
        }
        let calved_m3 = calved_acc.finish();
        let closing_m3: f64 = compensated_dot(&depth, area);
        let cap_thresh = GLACIER_MAX_DEPTH_M * (1.0 - 1e-9);
        let mut capped_nodes = 0_usize;
        for i in 0..n {
            if depth[i] > 0.0 && thickness[i] >= cap_thresh {
                capped_nodes += 1;
            }
        }

        let mut is_warm = Vec::with_capacity(n);
        let mut is_ice = Vec::with_capacity(n);
        let mut ice_idx = Vec::new();
        for i in 0..n {
            let surf_t = air_temperature_at_elevation_scalar(temp[i], 0.0, thickness[i]);
            let fric_i = fric.map_or(0.0, |v| v[i]);
            is_warm.push(warm_based_scalar(thickness[i], surf_t, fric_i));
            let ice_here = surface[i].is_finite() && thickness[i] > 0.0;
            is_ice.push(ice_here);
            if ice_here {
                ice_idx.push(i);
            }
        }

        let mut flux = vec![0.0_f64; n];
        let mut sliding = vec![0.0_f64; n];
        if !ice_idx.is_empty() && years > 0.0 {
            let mut source = vec![0.0_f64; n];
            let mut masked_surface = vec![f64::NEG_INFINITY; n];
            for &i in &ice_idx {
                let fric_i = fric.map_or(0.0, |v| v[i]);
                let bfs_i = bfs.map(|v| v[i]);
                let (acc, sm, bm) =
                    node_rates_scalar(thickness[i], temp[i], amp[i], precip[i], fric_i, bfs_i);
                let bal_m_yr = covered[i] * (acc - sm - bm) / 1.0e6;
                source[i] = (bal_m_yr - (depth[i] - opening_depth[i]) / years) * area[i];
                masked_surface[i] = surface[i];
            }
            ice_idx.sort_by(|&a, &b| surface[b].total_cmp(&surface[a]).then_with(|| a.cmp(&b)));
            flux = balance_flux(
                &ice_idx,
                &rust_graph.indptr,
                &rust_graph.indices,
                &rust_graph.lengths,
                &masked_surface,
                &is_ice,
                &source,
            );
            for &i in &ice_idx {
                let width = clip_min(area[i], 1.0).sqrt();
                let vel = flux[i] / (width * thickness[i].max(ICE_VELOCITY_MIN_THICKNESS_M));
                if is_warm[i] {
                    sliding[i] = ICE_SLIDING_FRACTION * vel;
                }
            }
        }

        let mut melt_m = Vec::with_capacity(n);
        let mut ice_release_m = Vec::with_capacity(n);
        let mut surface_melt_by_node_m = Vec::with_capacity(n);
        let mut bare_snowmelt_m = Vec::with_capacity(n);
        let mut land_acc = NeumaierSum::default();
        let mut ocean_retained_acc = NeumaierSum::default();
        let mut retained_acc = NeumaierSum::default();
        let mut bare_melted_acc = NeumaierSum::default();
        for i in 0..n {
            melt_m.push(melt_m3[i] / if area[i] > 0.0 { area[i] } else { 1.0 });
            ice_release_m.push(ice_release_m3[i] / if area[i] > 0.0 { area[i] } else { 1.0 });
            surface_melt_by_node_m
                .push(surface_melt_by_node_m3[i] / if area[i] > 0.0 { area[i] } else { 1.0 });
            bare_snowmelt_m.push(bare_snowmelt_m3[i] / if area[i] > 0.0 { area[i] } else { 1.0 });
            let mut smm = match dry_snow {
                Some(ds) => ds[i],
                None => snow_fraction_scalar(temp[i], amp[i], 1.0) * precip[i],
            };
            if let Some(b_arr) = bfs {
                smm += b_arr[i];
            }
            let f_m3 = clip_min(smm, 0.0) / 1000.0 * area[i] * years;
            let r_m3 = f_m3 * retained_snow[i];
            retained_acc.add(r_m3);
            if ocean[i] {
                ocean_retained_acc.add(r_m3);
            } else {
                land_acc.add(f_m3);
                bare_melted_acc.add(f_m3 - r_m3);
            }
        }
        let land_m3 = land_acc.finish();
        let ocean_retained_m3 = ocean_retained_acc.finish();
        let retained_m3 = retained_acc.finish();
        let bare_melted_m3 = bare_melted_acc.finish();
        if probe {
            eprintln!(
                "MANTLE_BLOOM_ICE_PROFILE final_us={}",
                t0.elapsed().as_micros()
            );
        }

        (
            depth,
            melt_m,
            ice_release_m,
            surface_melt_by_node_m,
            bare_snowmelt_m,
            thickness,
            covered,
            is_warm,
            flux,
            sliding,
            calved,
            (
                opening_m3,
                accumulation_m3,
                surface_melt_m3,
                basal_melt_m3,
                calved_m3,
                overflow_m3,
                closing_m3,
                snow_melted_m3,
                capped_m3,
                capped_nodes,
            ),
            (land_m3, ocean_retained_m3, retained_m3, bare_melted_m3),
            body_outs,
        )
    });

    let (
        depth,
        melt_m,
        ice_release_m,
        surface_melt_by_node_m,
        bare_snowmelt_m,
        thickness,
        covered,
        is_warm,
        flux,
        sliding,
        calved,
        ledger_tuple,
        frozen_tuple,
        body_outs,
    ) = result;

    let mut py_bodies = Vec::with_capacity(body_outs.len());
    for b in body_outs {
        let tuple = PyTuple::new(
            py,
            &[
                PyArray1::from_vec(py, b.nodes).into_any(),
                PyArray1::from_vec(py, b.domain_nodes).into_any(),
                b.opening_m3.into_pyobject(py)?.into_any(),
                b.accumulation_m3.into_pyobject(py)?.into_any(),
                b.surface_melt_m3.into_pyobject(py)?.into_any(),
                b.basal_melt_m3.into_pyobject(py)?.into_any(),
                b.calved_m3.into_pyobject(py)?.into_any(),
                b.overflow_m3.into_pyobject(py)?.into_any(),
                b.closing_m3.into_pyobject(py)?.into_any(),
                b.capped_m3.into_pyobject(py)?.into_any(),
                b.margin_m.into_pyobject(py)?.into_any(),
                b.saturated.into_pyobject(py)?.to_owned().into_any(),
            ],
        )?;
        py_bodies.push(tuple.into_any());
    }

    let out = PyTuple::new(
        py,
        &[
            PyArray1::from_vec(py, depth).into_any(),
            PyArray1::from_vec(py, melt_m).into_any(),
            PyArray1::from_vec(py, ice_release_m).into_any(),
            PyArray1::from_vec(py, surface_melt_by_node_m).into_any(),
            PyArray1::from_vec(py, bare_snowmelt_m).into_any(),
            PyArray1::from_vec(py, thickness).into_any(),
            PyArray1::from_vec(py, covered).into_any(),
            PyArray1::from_vec(py, is_warm).into_any(),
            PyArray1::from_vec(py, flux).into_any(),
            PyArray1::from_vec(py, sliding).into_any(),
            PyArray1::from_vec(py, calved).into_any(),
            ledger_tuple.into_pyobject(py)?.into_any(),
            frozen_tuple.into_pyobject(py)?.into_any(),
            py_bodies.into_pyobject(py)?.into_any(),
        ],
    )?;

    Ok(out.into_any().unbind())
}
