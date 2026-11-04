#!/usr/bin/env python3
"""Measure RFC 0018 graph-rift candidates at current-source split opportunities.

The runner uses normal production world generation and stepping, while monkeypatching only
the split decision to collect diagnostics before calling the original production function.
It never substitutes the diagnostic cut or otherwise changes split outcomes. Run controls
and the opening-residual rule are documented in analysis/issue797/README.md.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_rfc0018_g0.py \
        --out-json /private/tmp/rfc0018-g0-results.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse import diags
from scipy.sparse.csgraph import connected_components
from scipy.cluster.vq import kmeans2
from scipy.sparse.linalg import eigsh
from scipy.spatial import cKDTree

from app import cratons, faults, lithosphere, mantle, merge_split, world as world_mod
from app.elevation_lines import line_spacing_rad
from app.sparse_quad_patch import lattice_points, unpack_cell_keys

MAX_CANDIDATES = 12
MIN_COMPONENT_SHARE = 1.0
FAULT_REACH_KM = 250.0


def _shared_edges(plate):
    """Return unique node pairs and their exact spherical interface length in metres.

    `_probe_neighbour_indices` exposes both half-edges, which preserves length across 2:1
    adaptive cells. Equal-level pairs appear twice and their segment lengths are summed.
    """
    keys = plate.cell_keys.view(np.int64)
    face, level, ci, cj = unpack_cell_keys(keys)
    probes = plate._probe_neighbour_indices()
    segments = {}
    for node in range(len(keys)):
        scale = 1 << int(level[node])
        for side in range(4):
            for half, neighbor in enumerate(probes[node, side]):
                if neighbor < 0 or node >= neighbor:
                    continue
                t0, t1 = half * 0.5, (half + 1) * 0.5
                if side == 0:
                    a0, b0, a1, b1 = ci[node] + t0, cj[node], ci[node] + t1, cj[node]
                elif side == 1:
                    a0, b0, a1, b1 = ci[node] + 1, cj[node] + t0, ci[node] + 1, cj[node] + t1
                elif side == 2:
                    a0, b0, a1, b1 = ci[node] + 1 - t0, cj[node] + 1, ci[node] + 1 - t1, cj[node] + 1
                else:
                    a0, b0, a1, b1 = ci[node], cj[node] + 1 - t0, ci[node], cj[node] + 1 - t1
                p = lattice_points(np.array([face[node], face[node]]),
                                   np.array([a0, a1]) / scale, np.array([b0, b1]) / scale,
                                   plate.cells_per_edge)
                angle = math.acos(float(np.clip(np.dot(p[0], p[1]), -1.0, 1.0)))
                pair = (min(node, int(neighbor)), max(node, int(neighbor)))
                segments[pair] = segments.get(pair, 0.0) + angle * lithosphere.PLANET_RADIUS_M
    pairs = sorted(segments)
    i = np.fromiter((p[0] for p in pairs), dtype=np.int64)
    j = np.fromiter((p[1] for p in pairs), dtype=np.int64)
    lengths = np.fromiter((segments[p] for p in pairs), dtype=float)
    return i, j, lengths


def _fault_mask(world, plate, points: np.ndarray) -> np.ndarray:
    # Faults remain stored as scars after becoming inactive; include both active and historic
    # traces. The saved trace is world-space when available, otherwise reconstruct its frame.
    traces = []
    for fault in world.faults:
        if fault.plate_id != plate.plate_id:
            continue
        trace = fault.world_polyline
        if trace is None or len(trace) == 0:
            trace = faults.fault_world_points(fault, plate)
        if trace is not None and len(trace):
            traces.append(np.asarray(trace, dtype=float))
    if not traces:
        return np.zeros(len(points), dtype=bool)
    tree = cKDTree(np.concatenate(traces))
    chord = 2.0 * np.sin(0.5 * FAULT_REACH_KM / lithosphere.PLANET_RADIUS_KM)
    distance, _ = tree.query(points, distance_upper_bound=chord)
    return np.isfinite(distance)


def _node_resistance(world, plate, points: np.ndarray, flow: np.ndarray,
                     cluster_omega: np.ndarray, candidate_side: np.ndarray | None = None):
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    crust_density = lithosphere.node_crust_density(plate.collect("crust_type_code"), plate.crust_type)
    areal_mass = hc * crust_density + hm * lithosphere.RHO_LITHOSPHERE_MANTLE
    median_mass = max(float(np.median(areal_mass)), 1.0)

    # Thicker columns cost more to break. This is normalized so it remains comparable across
    # plates and resolutions; exact area enters the daughter masses below.
    resistance = 0.6 + 0.4 * np.clip(areal_mass / median_mass, 0.0, 4.0)
    craton = plate.collect("craton_crust_m") > 0.0
    resistance *= np.where(craton, 2.5, 1.0)

    # Residual to the closer of the two fitted mantle-flow regimes: coherent nodes resist,
    # while nodes that depart from either rigid-motion fit are easier to break.
    if candidate_side is None:
        predicted = [np.cross(omega, points) for omega in cluster_omega]
        residual = np.minimum(np.linalg.norm(flow - predicted[0], axis=1),
                              np.linalg.norm(flow - predicted[1], axis=1))
    else:
        omega_a = mantle.fit_euler_pole(points[candidate_side], flow[candidate_side])
        omega_b = mantle.fit_euler_pole(points[~candidate_side], flow[~candidate_side])
        predicted = np.empty_like(flow)
        predicted[candidate_side] = np.cross(omega_a, points[candidate_side])
        predicted[~candidate_side] = np.cross(omega_b, points[~candidate_side])
        residual = np.linalg.norm(flow - predicted, axis=1)
    scale = max(float(np.median(residual)), 1.0e-15)
    resistance /= 1.0 + np.clip(residual / scale, 0.0, 4.0)

    # Inherited belts/sutures and accumulated weakening signals lower edge resistance.
    belt = cratons.belt_mask(plate, world.elapsed_years)
    # The current source tree stores craton/belt history but has no per-node suture field or
    # suture-mask API. Preserve the fixed RFC multiplier while recording it as inactive; do
    # not guess that every belt is a suture.
    suture = np.zeros(len(points), dtype=bool)
    volcano = plate.collect("is_volcano").astype(bool)
    plume = plate.collect("plume_heat_myr") > 0.0
    fault = _fault_mask(world, plate, points)
    thin = areal_mass < median_mass
    weakening = np.ones(len(points), dtype=float)
    for mask, reduction in ((belt, 0.35), (suture, 0.45), (volcano, 0.35),
                            (plume, 0.35), (fault, 0.40), (thin, 0.25)):
        weakening[mask] *= 1.0 - reduction
    resistance *= np.clip(weakening, 0.12, 1.0)
    return np.maximum(resistance, 0.02), areal_mass, {
        "craton_nodes": int(craton.sum()), "belt_nodes": int(belt.sum()),
        "suture_nodes": int(suture.sum()), "suture_signal_available": False,
        "volcano_nodes": int(volcano.sum()),
        "plume_heated_nodes": int(plume.sum()), "historic_or_active_fault_nodes": int(fault.sum()),
        "thinned_nodes_below_plate_median_mass": int(thin.sum()),
    }


def _fiedler_vector(n, i, j, edge_cost):
    """Approximate a low-resistance balanced partition with the graph's second mode."""
    conductance = 1.0 / np.maximum(edge_cost, 1e-10)
    rows = np.concatenate([i, j])
    cols = np.concatenate([j, i])
    values = np.concatenate([conductance, conductance])
    adjacency = coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()
    degree = np.asarray(adjacency.sum(axis=1)).ravel()
    inverse_root = 1.0 / np.sqrt(np.maximum(degree, 1e-15))
    normalized = diags(inverse_root) @ (diags(degree) - adjacency) @ diags(inverse_root)
    indices = np.arange(1, n + 1, dtype=float)
    v0 = np.cos(indices * 0.6180339887498948) + np.sin(indices * 0.1732050807568877)
    _eigenvalues, eigenvectors = eigsh(normalized, k=2, which="SM", tol=1e-5, v0=v0)
    return eigenvectors[:, 1]


def _stable_u32(*parts) -> float:
    payload = "|".join(str(part) for part in parts).encode("utf-8")
    value = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return (value + 0.5) / 2**64


def _largest_component_share(n, i, j, side):
    shares = []
    for selected in (side, ~side):
        local = selected[i] & selected[j]
        indices = np.flatnonzero(selected)
        remap = np.full(n, -1, dtype=np.int64)
        remap[indices] = np.arange(len(indices))
        graph = coo_matrix((np.ones(local.sum()), (remap[i[local]], remap[j[local]])),
                           shape=(len(indices), len(indices)))
        _, labels = connected_components(graph, directed=False)
        shares.append(float(np.bincount(labels).max() / len(indices)) if len(indices) else 0.0)
    return shares


def _trace_stats(points, i, j, lengths, crossing, boundary_nodes, spacing_rad):
    """Approximate the cut's dual trace by clustering crossed-edge midpoints on the sphere."""
    cut_i, cut_j = i[crossing], j[crossing]
    midpoint = points[cut_i] + points[cut_j]
    midpoint /= np.maximum(np.linalg.norm(midpoint, axis=1)[:, None], 1e-15)
    if not len(midpoint):
        return 0, 0.0, 0, np.zeros(0, dtype=bool)
    # Neighboring cell interfaces meet within about one nominal cell spacing. The tolerance
    # accommodates the adaptive quad mesh while keeping distant boundary fragments separate.
    join_chord = 2.0 * np.sin(0.75 * spacing_rad)
    pairs = cKDTree(midpoint).query_pairs(join_chord, output_type="ndarray")
    graph = coo_matrix((np.ones(2 * len(pairs)),
                        (np.concatenate([pairs[:, 0], pairs[:, 1]]),
                         np.concatenate([pairs[:, 1], pairs[:, 0]]))),
                       shape=(len(midpoint), len(midpoint)))
    _count, labels = connected_components(graph, directed=False)
    component_lengths = np.bincount(labels, weights=lengths[crossing])
    dominant = int(np.argmax(component_lengths))
    dominant_mask = labels == dominant
    main_edges = np.unique(np.concatenate([cut_i[dominant_mask], cut_j[dominant_mask]]))
    contact_nodes = main_edges[boundary_nodes[main_edges]]
    if len(contact_nodes):
        contact_mask = np.zeros(len(points), dtype=bool)
        contact_mask[contact_nodes] = True
        contact_edges = contact_mask[i] & contact_mask[j]
        remap = np.full(len(points), -1, dtype=np.int64)
        remap[contact_nodes] = np.arange(len(contact_nodes))
        contact_graph = coo_matrix(
            (np.ones(2 * int(contact_edges.sum())),
             (np.concatenate([remap[i[contact_edges]], remap[j[contact_edges]]]),
              np.concatenate([remap[j[contact_edges]], remap[i[contact_edges]]]))),
            shape=(len(contact_nodes), len(contact_nodes)))
        contact_count = int(connected_components(contact_graph, directed=False)[0])
    else:
        contact_count = 0
    share = float(component_lengths[dominant] / max(lengths[crossing].sum(), 1e-15))
    return int(_count), share, contact_count, dominant_mask


def _opening(points, flow, i, j, edge_length, side, mass, parent_omega):
    crossing = side[i] != side[j]
    ii, jj = i[crossing], j[crossing]
    lengths = edge_length[crossing]
    if not len(ii):
        return None
    # Fit each daughter as a rigid Euler rotation to the local mantle flow.
    omega_a = mantle.fit_euler_pole(points[side], flow[side])
    omega_b = mantle.fit_euler_pole(points[~side], flow[~side])
    # Orient each crossed edge from side A to side B; the tangent projection measures opening.
    a, b = np.where(side[ii], ii, jj), np.where(side[ii], jj, ii)
    p, q = points[a], points[b]
    tangent = q - np.sum(q * p, axis=1)[:, None] * p
    tangent /= np.maximum(np.linalg.norm(tangent, axis=1)[:, None], 1e-15)
    va = np.cross(omega_a, p) * lithosphere.PLANET_RADIUS_M
    vb = np.cross(omega_b, q) * lithosphere.PLANET_RADIUS_M
    speeds = np.sum((vb - va) * tangent, axis=1)
    residual_a = flow[a] - np.cross(omega_a, p)
    residual_b = flow[b] - np.cross(omega_b, q)
    residual_difference = lithosphere.PLANET_RADIUS_M * np.sum((residual_b - residual_a) * tangent, axis=1)
    length_total = float(lengths.sum())
    net_opening_m_per_yr = float(np.dot(speeds, lengths) / length_total)
    mean_positive_m_per_yr = float(np.dot(np.maximum(speeds, 0.0), lengths) / length_total)
    open_fraction = float(lengths[speeds > 0.0].sum() / length_total)
    sigma_m_per_yr = float(np.sqrt(np.dot(np.square(residual_difference), lengths) / length_total))
    mass_a, mass_b = float(mass[side].sum()), float(mass[~side].sum())
    reduced_mass = mass_a * mass_b / max(mass_a + mass_b, 1.0)
    parent_v = np.cross(parent_omega[None, :], points) * lithosphere.PLANET_RADIUS_M
    daughter_v = np.empty_like(parent_v)
    daughter_v[side] = np.cross(omega_a[None, :], points[side]) * lithosphere.PLANET_RADIUS_M
    daughter_v[~side] = np.cross(omega_b[None, :], points[~side]) * lithosphere.PLANET_RADIUS_M
    momentum = float(np.dot(mass, np.linalg.norm(daughter_v - parent_v, axis=1)))
    small_share = min(mass_a, mass_b) / max(mass_a + mass_b, 1.0)
    balance = 4.0 * small_share * (1.0 - small_share)
    return {
        "omega_a_rad_per_yr": omega_a.tolist(), "omega_b_rad_per_yr": omega_b.tolist(),
        "mass_a_kg": mass_a, "mass_b_kg": mass_b,
        "smaller_mass_fraction": min(mass_a, mass_b) / max(mass_a + mass_b, 1.0),
        "reduced_mass_kg": reduced_mass,
        "cut_length_km": float(length_total * lithosphere.PLANET_RADIUS_M / 1000.0),
        "outward_opening_m_per_yr": mean_positive_m_per_yr,
        "outward_opening_cm_per_yr": net_opening_m_per_yr * 100.0,
        "net_opening_cm_per_yr": net_opening_m_per_yr * 100.0,
        "opening_uncertainty_cm_per_yr": 1.96 * sigma_m_per_yr * 100.0,
        "residual_envelope_rms_cm_per_yr": sigma_m_per_yr * 100.0,
        "opening_above_uncertainty": net_opening_m_per_yr > 1.96 * sigma_m_per_yr,
        "opening_edge_fraction": open_fraction,
        "potential_momentum_kg_m_per_yr": momentum,
        "reduced_mass_opening_kg_m_per_yr": reduced_mass * net_opening_m_per_yr,
        "balance_factor": balance,
        "momentum_scores_by_mu": {str(mu): momentum * balance**mu for mu in (0.0, 0.5, 1.0, 2.0)},
        "median_signed_opening_cm_per_yr": float(np.median(speeds) * 100.0),
        "crossing_i": ii.tolist(), "crossing_j": jj.tolist(),
    }


def choose_cut(world, plate, max_candidates=MAX_CANDIDATES, seed_tag=0):
    points, _ = plate.all_points_and_elevation()
    flow = mantle.flow_at(points, world.mantle_centers)
    centres, labels = kmeans2(flow, k=2, minit="++", seed=plate.plate_id + seed_tag)
    if len(np.unique(labels)) < 2:
        raise ValueError("mantle-flow clustering produced only one regime")
    cluster_omega = [mantle.fit_euler_pole(points[labels == k], flow[labels == k]) for k in (0, 1)]
    i, j, lengths = _shared_edges(plate)
    if not len(i):
        raise ValueError("plate adjacency graph has no edges")
    node_resistance, areal_mass, signals = _node_resistance(world, plate, points, flow, cluster_omega)
    spacing_rad = line_spacing_rad(world.node_density)
    median_edge_length = max(float(np.median(lengths)), 1.0)
    edge_cost = lengths / median_edge_length * 0.5 * (node_resistance[i] + node_resistance[j])
    unlike = labels[i] != labels[j]
    cluster_edges = np.flatnonzero(unlike)
    if not len(cluster_edges):
        raise ValueError("no adjacency edges cross the two flow clusters")
    areas = plate.node_areas_m2()
    mass = areal_mass * areas
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    _, radius = merge_split.geometry.bounding_sphere(points)
    relax = max(min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    pole_threshold = merge_split.SPLIT_MIN_POLE_SEPARATION * (1.0 - relax)
    fiedler = _fiedler_vector(len(points), i, j, edge_cost)
    # Favor weak edges that also bridge distinct parts of the low-conductance partition.
    fiedler_gap = np.abs(fiedler[i] - fiedler[j])
    eligible = np.flatnonzero(unlike & (fiedler_gap > 1e-12))
    if not len(eligible):
        raise ValueError("no weak adjacency edge crosses the spectral partition")
    weights = fiedler_gap[eligible] / np.maximum(edge_cost[eligible], 1e-10)
    weights /= weights.max()
    count = min(max_candidates, len(eligible))
    node_ids = plate.cell_keys.view(np.uint64)
    priorities = np.array([
        -math.log(_stable_u32(world.seed, plate.plate_id, round(world.elapsed_years),
                              int(node_ids[i[e]]), int(node_ids[j[e]]))) / max(float(w), 1e-15)
        for e, w in zip(eligible, weights)
    ])
    selected = eligible[np.argsort(priorities, kind="stable")[:count]]
    variants = []
    boundary_nodes = plate.exposed_sides().any(axis=1)
    for edge_idx in selected:
        threshold = 0.5 * (fiedler[i[edge_idx]] + fiedler[j[edge_idx]])
        side = fiedler > threshold
        # The RFC's motion-fit resistance depends on the daughter assignment. Iterate the
        # candidate's fit, resistance, spectral cut and sampled-edge threshold to a fixed
        # assignment (or stop after eight bounded updates if it cycles).
        fit_iterations = 0
        for fit_iterations in range(1, 9):
            candidate_resistance, _, _ = _node_resistance(
                world, plate, points, flow, np.asarray(cluster_omega), side)
            candidate_cost = lengths / median_edge_length * 0.5 * (
                candidate_resistance[i] + candidate_resistance[j])
            candidate_fiedler = _fiedler_vector(len(points), i, j, candidate_cost)
            next_side = candidate_fiedler > 0.5 * (
                candidate_fiedler[i[edge_idx]] + candidate_fiedler[j[edge_idx]])
            if np.array_equal(next_side, side):
                break
            side = next_side
        candidate_resistance, _, _ = _node_resistance(
            world, plate, points, flow, np.asarray(cluster_omega), side)
        candidate_cost = lengths / median_edge_length * 0.5 * (
            candidate_resistance[i] + candidate_resistance[j])
        if not side.any() or side.all():
            continue
        components = _largest_component_share(len(points), i, j, side)
        opening = _opening(points, flow, i, j, lengths, side, mass, plate.omega)
        if opening is None:
            continue
        cut = side[i] != side[j]
        cut_nodes = np.unique(np.concatenate([i[cut], j[cut]]))
        boundary_touches = int(boundary_nodes[cut_nodes].sum())
        trace_count, trace_share, main_boundary_touches, main_trace = _trace_stats(
            points, i, j, lengths, cut, boundary_nodes, spacing_rad)
        omega_a = mantle.fit_euler_pole(points[side], flow[side])
        omega_b = mantle.fit_euler_pole(points[~side], flow[~side])
        daughter_size_ok = int(side.sum()) >= min_nodes and int((~side).sum()) >= min_nodes
        pole_separation = float(np.linalg.norm(omega_a - omega_b))
        pole_separation_ok = pole_separation >= pole_threshold
        connected = min(components) >= MIN_COMPONENT_SHARE
        trace_ok = trace_share >= 0.90 and main_boundary_touches >= 2
        usable = (connected and trace_ok and daughter_size_ok and pole_separation_ok
                  and opening["opening_above_uncertainty"])
        variants.append({
            "seed_edge": [int(node_ids[i[edge_idx]]), int(node_ids[j[edge_idx]])],
            "candidate_fit_iterations": fit_iterations,
            "cut_cost": float(candidate_cost[cut].sum()),
            "min_component_share": min(components),
            "components": components,
            "admissible_connected": connected,
            "boundary_touch_nodes": boundary_touches,
            "reaches_plate_edge_twice": boundary_touches >= 2,
            "cut_trace_components": trace_count,
            "largest_trace_length_share": trace_share,
            "largest_trace_boundary_touch_nodes": main_boundary_touches,
            "single_edge_to_edge_trace": trace_share >= 0.90 and main_boundary_touches >= 2,
            "daughter_size_ok": daughter_size_ok,
            "daughter_node_counts": [int(side.sum()), int((~side).sum())],
            "pole_separation_rad_per_yr": pole_separation,
            "pole_separation_threshold_rad_per_yr": pole_threshold,
            "pole_separation_ok": pole_separation_ok,
            "usable_candidate": usable,
            "opening": opening,
            "side_mask": side,
            "cut_edges_i": i[cut], "cut_edges_j": j[cut],
        })
    usable = [v for v in variants if v["usable_candidate"]]
    pool = usable if usable else variants
    selected_cut = max(pool, key=lambda v: v["opening"]["momentum_scores_by_mu"]["1.0"])
    selected_cut["accepted"] = bool(selected_cut["usable_candidate"])
    if selected_cut["accepted"]:
        selected_cut["rejected_reason"] = None
    elif not selected_cut["admissible_connected"]:
        selected_cut["rejected_reason"] = "daughter side is disconnected"
    elif not selected_cut["single_edge_to_edge_trace"]:
        selected_cut["rejected_reason"] = "cut lacks a dominant edge-to-edge trace"
    elif not selected_cut["daughter_size_ok"]:
        selected_cut["rejected_reason"] = "daughter fails minimum-size check"
    elif not selected_cut["pole_separation_ok"]:
        selected_cut["rejected_reason"] = "daughter pole separation is below the existing threshold"
    else:
        selected_cut["rejected_reason"] = "net opening does not exceed fit-residual uncertainty"
    return points, selected_cut, variants, signals


PRIMARY_SEEDS = [1629738962, 1097710202, 591453549, 969845080, 763144949, 870996356,
                 494568917, 1854897345, 954941394, 600370317, 1643146570, 146046339,
                 1397817354, 1630166320, 604997962, 704718994, 1097624328, 487506346,
                 106533960, 1439649077]
REPLACEMENTS = [1227913777, 303885793, 2045331463, 157959777, 1208093067]


def _is_opportunity(world, plate):
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    if plate.node_count() < min_nodes or plate.age_years < merge_split.SPLIT_MIN_AGE_YEARS:
        return False
    points, _ = plate.all_points_and_elevation()
    flow = mantle.flow_at(points, world.mantle_centers)
    _, radius = merge_split.geometry.bounding_sphere(points)
    relax = max(min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    threshold = merge_split.SPLIT_RMS_RESIDUAL_THRESHOLD * (1.0 - relax)
    return merge_split._fit_residual_rms(points, flow, plate.omega) >= threshold


def _jsonable_candidates(variants):
    result = []
    for candidate in variants:
        opening = {k: v for k, v in candidate["opening"].items()
                   if k not in ("crossing_i", "crossing_j")}
        result.append({k: v for k, v in candidate.items()
                       if k not in ("side_mask", "cut_edges_i", "cut_edges_j", "opening")}
                      | {"opening": opening})
    return result


def _save_report(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(report, separators=(",", ":")) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out-json", type=Path, default=Path("/private/tmp/rfc0018-g0-results.json"))
    parser.add_argument("--horizon-myr", type=int, default=2)
    parser.add_argument("--step-myr", type=int, default=1)
    args = parser.parse_args()
    report = {"protocol": "RFC0018-G0-v1", "source_commit": "63290b053de15b72280fe28fe0086ff37affdba2",
              "world_defaults": True, "horizon_myr": args.horizon_myr, "step_myr": args.step_myr,
              "primary_seeds": PRIMARY_SEEDS, "replacement_seeds": REPLACEMENTS,
              "worlds": [], "baseline_modified": False}
    original = merge_split.maybe_split_plate
    for seed in PRIMARY_SEEDS + REPLACEMENTS:
        if len([w for w in report["worlds"] if w.get("constructible")]) >= 20:
            break
        world_started = time.perf_counter()
        row = {"seed": seed, "constructible": False, "opportunities": []}
        try:
            world = world_mod.generate_world(seed)
            if world is None or not world.plates:
                raise RuntimeError("generation returned no world or no plates")
            row["constructible"] = True
        except Exception as exc:
            row["generation_error"] = f"{type(exc).__name__}: {exc}"
            row["runtime_s"] = time.perf_counter() - world_started
            report["worlds"].append(row)
            _save_report(args.out_json, report)
            continue

        original_for_world = original
        diagnosed_step = None
        step_opportunities = {}
        def instrumented(current_world, plate, step_years=0.0):
            nonlocal diagnosed_step, step_opportunities
            step_key = round(current_world.elapsed_years)
            if diagnosed_step != step_key:
                # Capture every eligible plate from the same pre-split state. Production only
                # permits one successful rift per step; probing inside its plate loop alone
                # would omit later opportunities whenever an earlier plate split.
                diagnosed_step = step_key
                step_opportunities = {}
                for probe_plate in list(current_world.plates):
                    if not _is_opportunity(current_world, probe_plate):
                        continue
                    started = time.perf_counter()
                    opportunity = {"seed": int(current_world.seed), "plate_id": int(probe_plate.plate_id),
                                   "elapsed_myr": current_world.elapsed_years / 1e6,
                                   "node_count": probe_plate.node_count(),
                                   "production_evaluated": False}
                    try:
                        _, chosen, candidates, signals = choose_cut(current_world, probe_plate)
                        opportunity.update({"status": "measured", "edge_signals": signals,
                                            "selected": {k: v for k, v in chosen.items()
                                                         if k not in ("side_mask", "cut_edges_i", "cut_edges_j", "opening")}
                                                         | {"opening": {k: v for k, v in chosen["opening"].items()
                                                                        if k not in ("crossing_i", "crossing_j")}},
                                            "candidates": _jsonable_candidates(candidates)})
                    except Exception as exc:
                        opportunity.update({"status": "diagnostic_error",
                                            "error": f"{type(exc).__name__}: {exc}"})
                    opportunity["diagnostic_runtime_s"] = time.perf_counter() - started
                    step_opportunities[probe_plate.plate_id] = opportunity
                    row["opportunities"].append(opportunity)

            prior_age = plate.age_years
            prior_event_count = len(current_world.events)
            actual = original_for_world(current_world, plate, step_years)
            opportunity = step_opportunities.get(plate.plate_id)
            if opportunity is not None:
                new_events = [message for _, message in current_world.events[prior_event_count:]]
                failed_rift = any("rift failed" in message.lower() or "began rifting" in message.lower()
                                  for message in new_events)
                opportunity["production_evaluated"] = True
                opportunity["production_split"] = actual is not None
                opportunity["production_result"] = (
                    "split" if actual is not None else
                    "failed_rift_stochastic" if failed_rift or plate.age_years < prior_age else
                    "rejected_by_existing_gate"
                )
                opportunity["production_events"] = new_events
            return actual

        merge_split.maybe_split_plate = instrumented
        stepping_started = time.perf_counter()
        step_years = args.step_myr * 1_000_000
        try:
            while world.elapsed_years < args.horizon_myr * 1_000_000:
                world_mod.step_world(world, step_years)
        finally:
            merge_split.maybe_split_plate = original
        for opportunity in row["opportunities"]:
            if not opportunity["production_evaluated"]:
                opportunity["production_result"] = "not_evaluated_after_earlier_split"
        row["elapsed_myr"] = world.elapsed_years / 1e6
        row["step_runtime_s"] = time.perf_counter() - stepping_started
        row["runtime_s"] = time.perf_counter() - world_started
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        row["max_rss_mb_process"] = rss / (1024.0 * 1024.0) if sys.platform == "darwin" else rss / 1024.0
        report["worlds"].append(row)
        report["opportunity_count"] = sum(len(w["opportunities"]) for w in report["worlds"])
        _save_report(args.out_json, report)
        print(json.dumps({"seed": seed, "constructible": True,
                          "opportunities": len(row["opportunities"]), "runtime_s": row["runtime_s"]}), flush=True)
    report["opportunity_count"] = sum(len(w["opportunities"]) for w in report["worlds"])
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    report["peak_rss_mb_process"] = rss / (1024.0 * 1024.0) if sys.platform == "darwin" else rss / 1024.0
    _save_report(args.out_json, report)
    print(json.dumps({"worlds": len(report["worlds"]), "constructible": sum(w["constructible"] for w in report["worlds"]),
                      "opportunities": report["opportunity_count"], "out_json": str(args.out_json)}, indent=2))


if __name__ == "__main__":
    main()
