#!/usr/bin/env python3
"""Couple a production slab/midline load snapshot to local damage and a crack-gated split.

This is a diagnostic experiment, not production rifting: it freezes one generated plate's
geometry and measured force pattern, ramps only the midline-source force, integrates the
Gate 3a local scalar damage law at each node, and tests whether failed cells form a connected
graph separator before constructing two material-retaining daughter masks.

Usage from the repository root:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_crack_coupled.py \
        --out analysis/issue699/crack-coupled-edge-links-seed2024-density4.json
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from app import cratons, elevation_lines, geometry, lithosphere, mantle, thin_shell, torque, world as world_mod
import probe_gate3b_component_shell as shell_audit
import probe_gate3b_midplate_current as current_audit
import probe_gate3b_slab_pull as slab_screen
import replay_gate2_capacity as capacity_audit

SEED = 2024
DENSITY = 4.0
STARTUP_YEARS = 2_000_000.0
STEP_YEARS = 2_000_000.0
MAX_HISTORY_YEARS = 1_000_000_000.0
SHELL_CAP = 512
EFFECTIVE_MODULUS_PA = 70.0e9
POISSON_RATIO = 0.25
# Gate 3a local-feedback diagnostic values; they remain scenario assumptions, not calibrated
# production constants.
EPSILON_ONSET = 0.25
EPSILON_FAIL = 1.25
DAMAGE_VISCOSITY_PA_S = 9.0e21
RESIDUAL_STRENGTH_FRACTION = 1.0 / 3.0
# Each tuple is elapsed Myr, effective slab-sink-length multiplier, midline-current multiplier.
# The native 250 km sink is a production constitutive length; the larger values are explicit
# load brackets, not dynamic slab geometry.
LOAD_RAMP = (
    (0.0, 1.0, 1.0),
    (10.0, 2.0, 1.0),
    (20.0, 4.0, 1.0),
    (30.0, 8.0, 1.0),
    (40.0, 8.0, 2.0),
    (50.0, 8.0, 4.0),
)
MIN_DAUGHTER_AREA_FRACTION = 0.10
MIN_INTACT_AREA_FRACTION = 0.80


def _capacity_fields(continental_mask):
    raw = capacity_audit._physical_capacities()
    names = ("p10", "median", "p90")
    continental = np.asarray(raw["continental"], dtype=float)
    oceanic = np.asarray(raw["oceanic_mature"], dtype=float)
    fields = {
        name: np.where(continental_mask, continental[index], oceanic[index])
        for index, name in enumerate(names)
    }
    references = {
        name: {
            "continental_n_per_m": float(continental[index]),
            "oceanic_mature_n_per_m": float(oceanic[index]),
        }
        for index, name in enumerate(names)
    }
    return fields, references


def _advance_damage(epsilon, stress_pa, yield_pa, duration_s):
    """Exact scalar update for the Gate 3a two-stage linear-softening damage law."""
    old = np.asarray(epsilon, dtype=float)
    stress = np.asarray(stress_pa, dtype=float)
    strength = np.asarray(yield_pa, dtype=float)
    new = old.copy()
    onset, failure = EPSILON_ONSET, EPSILON_FAIL
    span = failure - onset
    residual = RESIDUAL_STRENGTH_FRACTION

    below = old < onset
    positive = below & (stress > strength)
    time_to_onset = np.full(old.shape, np.inf)
    time_to_onset[positive] = (onset - old[positive]) * DAMAGE_VISCOSITY_PA_S / (
        stress[positive] - strength[positive]
    )
    partial = positive & (time_to_onset >= duration_s)
    new[partial] = old[partial] + (
        (stress[partial] - strength[partial]) * duration_s / DAMAGE_VISCOSITY_PA_S
    )
    through = positive & (time_to_onset < duration_s)
    new[through] = onset
    remaining = np.zeros(old.shape, dtype=float)
    remaining[through] = duration_s - time_to_onset[through]

    active = (~below) & (old < failure)
    remaining[active] = duration_s
    active |= through
    if np.any(active):
        beta = strength * (1.0 - residual) / span
        x = np.maximum(new - onset, 0.0)
        # Integrate dx/dt = (stress - strength + beta*x) / eta exactly while
        # the softened material remains above its current yield strength.
        driving = active & (remaining > 0.0) & (beta > 0.0)
        offset = np.zeros_like(old)
        offset[driving] = (stress[driving] - strength[driving]) / beta[driving]
        advancing = driving & (stress > strength - beta * x)
        exponent = np.zeros_like(old)
        exponent[advancing] = np.minimum(
            beta[advancing] * remaining[advancing] / DAMAGE_VISCOSITY_PA_S, 700.0
        )
        denom = x + offset
        valid = advancing & (denom > 0.0) & (span + offset > 0.0)
        x_new = x.copy()
        x_new[valid] = denom[valid] * np.exp(exponent[valid]) - offset[valid]
        new[advancing] = onset + np.maximum(x_new[advancing], 0.0)
        new = np.clip(new, 0.0, failure)
    return new


def _opening_capacity_demand(plate, source_force):
    """Map the diagnostic shell's principal tensile resultant to every plate cell."""
    all_indices = np.arange(plate.node_count(), dtype=np.int64)
    points, _ = plate.all_points_and_elevation()
    row = shell_audit._plate_shell(
        plate, all_indices, points, source_force, SHELL_CAP, include_demand=True
    )
    demand = np.asarray(row["component_demand_n_per_m"], dtype=float)
    direction = np.asarray(row["component_principal_extension_xyz"], dtype=float)
    if demand.shape != (plate.node_count(),):
        raise RuntimeError("shell demand did not map back to all production cells")
    if direction.shape != (plate.node_count(), 3):
        raise RuntimeError("shell principal-extension direction did not map to every cell")
    return row, demand, geometry.normalize(direction)


def _edge_cut_candidate(plate, cracked, extension_direction, stress_pa, yield_pa, areas):
    """A split candidate requires supported failed cell links to disconnect the plate graph."""
    adjacency = plate.adjacency()
    n = plate.node_count()
    rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(adjacency.offsets))
    cols = adjacency.neighbours.astype(np.int64, copy=False)
    pairs = np.unique(np.sort(np.column_stack((rows, cols)), axis=1), axis=0)
    points, _ = plate.all_points_and_elevation()
    edge_direction = points[pairs[:, 1]] - (
        np.sum(points[pairs[:, 1]] * points[pairs[:, 0]], axis=1)[:, None] * points[pairs[:, 0]]
    )
    edge_direction = geometry.normalize(edge_direction)
    opening_a = np.abs(np.sum(edge_direction * extension_direction[pairs[:, 0]], axis=1))
    opening_b = np.abs(np.sum(edge_direction * extension_direction[pairs[:, 1]], axis=1))
    # Use the edge-local opening direction and residual-strength support at both cells.
    supported = (opening_a > 0.5) & (opening_b > 0.5)
    supported &= (
        stress_pa[pairs[:, 0]] > RESIDUAL_STRENGTH_FRACTION * yield_pa[pairs[:, 0]]
    ) & (
        stress_pa[pairs[:, 1]] > RESIDUAL_STRENGTH_FRACTION * yield_pa[pairs[:, 1]]
    )
    failed_pairs = supported & cracked[pairs[:, 0]] & cracked[pairs[:, 1]]
    surviving = pairs[~failed_pairs]
    graph_rows = np.concatenate((surviving[:, 0], surviving[:, 1]))
    graph_cols = np.concatenate((surviving[:, 1], surviving[:, 0]))
    graph = csr_matrix((np.ones(len(graph_rows)), (graph_rows, graph_cols)), shape=(n, n))
    if not len(surviving):
        return None
    count, labels = connected_components(graph, directed=False)
    component_area = np.bincount(labels, weights=areas, minlength=count)
    order = np.argsort(component_area)[::-1]
    total_area = float(areas.sum())
    if count != 2:
        return None
    first, second = int(order[0]), int(order[1])
    area_a, area_b = float(component_area[first]), float(component_area[second])
    retained_fraction = (area_a + area_b) / max(total_area, 1.0)
    if min(area_a, area_b) < MIN_DAUGHTER_AREA_FRACTION * total_area:
        return None
    if retained_fraction < MIN_INTACT_AREA_FRACTION:
        return None
    side_a = labels == first
    side_b = labels == second
    daughters = plate._plates_from_node_masks([side_a, side_b], [plate.plate_id, plate.plate_id + 1])
    if len(daughters) != 2 or sum(d.node_count() for d in daughters) != n:
        return None
    daughter_areas = [float(d.node_areas_m2().sum()) for d in daughters]
    if min(daughter_areas) < MIN_DAUGHTER_AREA_FRACTION * total_area:
        return None
    return {
        "intact_connected_components": int(count),
        "failed_adjacency_link_count": int(failed_pairs.sum()),
        "daughter_area_km2": [value / 1.0e6 for value in daughter_areas],
        "daughter_area_fraction": [value / total_area for value in daughter_areas],
        "retained_cell_area_fraction": retained_fraction,
        "cracked_cell_count": int(cracked.sum()),
        "node_count_conserved": sum(d.node_count() for d in daughters) == n,
    }


def _connectivity_summary(plate, cracked, extension_direction, stress_pa, yield_pa, areas):
    adjacency = plate.adjacency()
    n = plate.node_count()
    rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(adjacency.offsets))
    cols = adjacency.neighbours.astype(np.int64, copy=False)
    pairs = np.unique(np.sort(np.column_stack((rows, cols)), axis=1), axis=0)
    points, _ = plate.all_points_and_elevation()
    edge_direction = points[pairs[:, 1]] - (
        np.sum(points[pairs[:, 1]] * points[pairs[:, 0]], axis=1)[:, None] * points[pairs[:, 0]]
    )
    edge_direction = geometry.normalize(edge_direction)
    opening_a = np.abs(np.sum(edge_direction * extension_direction[pairs[:, 0]], axis=1))
    opening_b = np.abs(np.sum(edge_direction * extension_direction[pairs[:, 1]], axis=1))
    supported = (opening_a > 0.5) & (opening_b > 0.5)
    supported &= (
        stress_pa[pairs[:, 0]] > RESIDUAL_STRENGTH_FRACTION * yield_pa[pairs[:, 0]]
    ) & (
        stress_pa[pairs[:, 1]] > RESIDUAL_STRENGTH_FRACTION * yield_pa[pairs[:, 1]]
    )
    failed_pairs = supported & cracked[pairs[:, 0]] & cracked[pairs[:, 1]]
    surviving = pairs[~failed_pairs]
    graph_rows = np.concatenate((surviving[:, 0], surviving[:, 1]))
    graph_cols = np.concatenate((surviving[:, 1], surviving[:, 0]))
    graph = csr_matrix((np.ones(len(graph_rows)), (graph_rows, graph_cols)), shape=(n, n))
    count, labels = connected_components(graph, directed=False)
    component_areas = np.bincount(labels, weights=areas, minlength=count)
    sorted_areas = np.sort(component_areas)[::-1]
    crack_nodes = np.flatnonzero(cracked)
    if len(crack_nodes):
        crack_edge_mask = cracked[rows] & cracked[cols]
        crack_graph = csr_matrix((np.ones(int(crack_edge_mask.sum())),
                                  (rows[crack_edge_mask], cols[crack_edge_mask])), shape=(n, n))
        crack_count, _ = connected_components(crack_graph[crack_nodes][:, crack_nodes], directed=False)
    else:
        crack_count = 0
    return {
        "failed_adjacency_link_count": int(failed_pairs.sum()),
        "failed_cell_count": int(cracked.sum()),
        "failed_cell_area_fraction": float(areas[cracked].sum() / max(float(areas.sum()), 1.0)),
        "failed_cell_component_count": int(crack_count),
        "plate_components_after_failed_links": int(count),
        "largest_daughter_area_fraction": float(sorted_areas[0] / max(float(areas.sum()), 1.0)),
        "second_daughter_area_fraction": float(sorted_areas[1] / max(float(areas.sum()), 1.0)) if count > 1 else 0.0,
        "two_largest_area_fraction": float(sorted_areas[:2].sum() / max(float(areas.sum()), 1.0)),
    }


def _connected_overcapacity(plate, demand, domain, areas, capacity):
    """Measure the largest connected over-capacity area on the eligible continental domain."""
    over = domain & (demand >= capacity)
    indices = np.flatnonzero(over)
    domain_area = float(areas[domain].sum())
    if not len(indices) or domain_area <= 0.0:
        return {
            "overcapacity_area_fraction": 0.0,
            "largest_connected_overcapacity_fraction": 0.0,
            "overcapacity_component_count": 0,
        }
    adjacency = plate.adjacency()
    n = plate.node_count()
    rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(adjacency.offsets))
    cols = adjacency.neighbours.astype(np.int64, copy=False)
    keep = over[rows] & over[cols]
    graph = csr_matrix((np.ones(int(keep.sum())), (rows[keep], cols[keep])), shape=(n, n))
    count, labels = connected_components(graph[indices][:, indices], directed=False)
    component_areas = np.bincount(labels, weights=areas[indices], minlength=count)
    return {
        "overcapacity_area_fraction": float(areas[indices].sum() / domain_area),
        "largest_connected_overcapacity_fraction": float(component_areas.max() / domain_area),
        "overcapacity_component_count": int(count),
    }


def _source_and_target(world):
    plate, candidates = slab_screen._select(world)
    if plate is None or plate.plate_id != 0:
        raise RuntimeError("seed-2024 target plate 0 did not reproduce")
    candidate = slab_screen._bilateral_candidate(world, plate)
    if not candidate:
        raise RuntimeError("target plate no longer meets bilateral slab-pull precondition")
    native_centers = list(world.mantle_centers)
    native_p95 = float(mantle.rad_per_yr_to_cm_per_yr(np.percentile(
        np.linalg.norm(mantle.flow_at(candidate["points"], native_centers), axis=1), 95.0
    )))
    source_geometry = current_audit._source_geometry(
        plate, shell_audit.line_spacing_rad(DENSITY)
    )
    centers, scale = current_audit._scaled_source_centers(
        source_geometry["source_centers_unit_strength"], candidate["points"], native_p95
    )
    world.mantle_centers = native_centers + centers
    field, inputs, subducting = slab_screen._field_and_inputs(world, plate)
    source_velocity_m_s = (
        mantle.flow_at(geometry.normalize(field.points_xyz), centers)
        * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    )
    midline_force = (
        torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
        * field.areas_m2[:, None] * source_velocity_m_s
    )
    # Keep the complete production force field (including the bilateral slab terms) and
    # scale only the new midline-current traction in the later sweep arms.
    return plate, candidates, candidate, native_p95, centers, scale, field, midline_force, inputs, subducting


def run(out_path: Path):
    started = time.perf_counter()
    world = world_mod.generate_world(SEED, node_density=DENSITY, cratons_first=True, synthetic_relief=False)
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    world_mod.step_world(world, STARTUP_YEARS)
    (plate, candidates, candidate, native_p95, centers, source_scale, field,
     midline_force, inputs, subducting) = _source_and_target(world)
    source_base = np.asarray(field.total_n, dtype=float)
    shell_by_factor = {}
    demand_by_factor = {}
    direction_by_factor = {}
    for _start_myr, slab_factor, midline_factor in LOAD_RAMP:
        factor = (slab_factor, midline_factor)
        if factor in demand_by_factor:
            continue
        scaled_force = (
            source_base
            + (slab_factor - 1.0) * field.slab_pull_n
            + (midline_factor - 1.0) * midline_force
        )
        shell, demand_by_factor[factor], direction_by_factor[factor] = _opening_capacity_demand(
            plate, scaled_force
        )
        shell["target_demand_quantiles_n_per_m"] = np.quantile(
            demand_by_factor[factor], [0.1, 0.5, 0.9]
        ).tolist()
        shell_by_factor[f"slab{slab_factor:g}_midline{midline_factor:g}"] = shell
    thickness = np.maximum(
        plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m"), 1.0
    )
    areas = plate.node_areas_m2()
    domain = elevation_lines.effective_is_continental_from_codes(
        plate.collect("crust_type_code"), plate.crust_type == "continental"
    )
    if not domain.any():
        raise RuntimeError("selected bilateral slab-pull target has no continental capacity domain")
    capacity_fields, capacity_references = _capacity_fields(domain)
    whole_plate = np.ones_like(domain, dtype=bool)
    capacity_event_screens = {}
    for factor, demand in demand_by_factor.items():
        screen = {}
        for case_name, capacity in capacity_fields.items():
            whole_footprint = _connected_overcapacity(
                plate, demand, whole_plate, areas, capacity
            )
            continental_footprint = _connected_overcapacity(
                plate, demand, domain, areas, capacity
            )
            screen[case_name] = {
                "capacity_reference_n_per_m_by_crust_type": capacity_references[case_name],
                "whole_plate": {
                    **whole_footprint,
                    "qualifies_25pct_for_10myr_at_this_load_stage": bool(
                        whole_footprint["largest_connected_overcapacity_fraction"] >= 0.25
                    ),
                },
                "continental_domain": {
                    **continental_footprint,
                    "qualifies_25pct_for_10myr_at_this_load_stage": bool(
                        continental_footprint["largest_connected_overcapacity_fraction"] >= 0.25
                    ),
                },
            }
        capacity_event_screens[f"slab{factor[0]:g}_midline{factor[1]:g}"] = screen
    # The screen runs separate local damage histories against the same fixed source snapshot.
    outcomes = {}
    for case_name, capacity_n_per_m in capacity_fields.items():
        epsilon = np.zeros(plate.node_count(), dtype=float)
        failed_history = []
        connectivity_history = []
        split = None
        elapsed_myr = 0.0
        steps = int(round(MAX_HISTORY_YEARS / STEP_YEARS))
        for step in range(1, steps + 1):
            elapsed_myr = step * STEP_YEARS / 1.0e6
            slab_factor, midline_factor = LOAD_RAMP[0][1:]
            for start_myr, slab_scale, midline_scale in LOAD_RAMP:
                if elapsed_myr >= start_myr:
                    slab_factor, midline_factor = slab_scale, midline_scale
            stress_pa = demand_by_factor[(slab_factor, midline_factor)] / thickness
            extension_direction = direction_by_factor[(slab_factor, midline_factor)]
            yield_pa = capacity_n_per_m / thickness
            epsilon = _advance_damage(epsilon, stress_pa, yield_pa, STEP_YEARS * 365.25 * 24.0 * 3600.0)
            cracked = epsilon >= EPSILON_FAIL
            failed_history.append(int(cracked.sum()))
            if step == 1 or step % 5 == 0 or cracked.any():
                connectivity_history.append({
                    "elapsed_myr": elapsed_myr,
                    "effective_slab_sink_length_multiplier": slab_factor,
                    "midline_force_multiplier": midline_factor,
                    "cracked_cells": int(cracked.sum()),
                    **_connectivity_summary(
                        plate, cracked, extension_direction, stress_pa, yield_pa, areas
                    ),
                })
            if cracked.any():
                split = _edge_cut_candidate(
                    plate, cracked, extension_direction, stress_pa, yield_pa, areas
                )
            if split is not None:
                break
        outcomes[case_name] = {
            "capacity_reference_n_per_m_by_crust_type": capacity_references[case_name],
            "completed_split_candidate": split is not None,
            "split_time_myr": elapsed_myr if split is not None else None,
            "right_censored_myr": None if split is not None else MAX_HISTORY_YEARS / 1.0e6,
            "split": split,
            "cracked_cell_count_at_end": failed_history[-1],
            "peak_damage_strain": float(epsilon.max()),
            "area_fraction_at_failure": float(areas[epsilon >= EPSILON_FAIL].sum() / areas.sum()),
            "history_every_step_failed_cells": failed_history,
            "connectivity_history": connectivity_history,
            "final_connectivity": _connectivity_summary(
                plate, epsilon >= EPSILON_FAIL, extension_direction, stress_pa, yield_pa, areas
            ),
        }
        if split is not None:
            break
    report = {
        "title": "RFC 0010 Gate 3b coupled load-to-crack-to-split diagnostic",
        "status": "diagnostic-only; no production damage or split code is invoked",
        "seed": SEED,
        "node_density": DENSITY,
        "startup_myr": STARTUP_YEARS / 1.0e6,
        "step_myr": STEP_YEARS / 1.0e6,
        "max_history_myr": MAX_HISTORY_YEARS / 1.0e6,
        "selected_plate_id": int(plate.plate_id),
        "plate_area_km2": float(areas.sum() / 1.0e6),
        "continental_capacity_domain_area_km2": float(areas[domain].sum() / 1.0e6),
        "continental_capacity_domain_fraction": float(areas[domain].sum() / areas.sum()),
        "bilateral_precondition": candidate["side_groups"],
        "selection_candidates": candidates,
        "midline_source_p95_cm_per_yr_at_1x": native_p95,
        "midline_source_force_scale": float(source_scale),
        "load_ramp": [{
            "elapsed_myr": start,
            "effective_slab_sink_length_multiplier": slab_factor,
            "midline_force_multiplier": midline_factor,
        } for start, slab_factor, midline_factor in LOAD_RAMP],
        "slab_sink_length_km": float(torque.SUBDUCTION_LSINK_M / 1000.0),
        "active_slab_front_length_km": float(torque.boundary_node_length_m(
            inputs, shell_audit.line_spacing_rad(DENSITY)
        )[subducting].sum() / 1000.0),
        "shell_mesh_by_load_case": {
            factor: {key: shell[key] for key in (
                "target_representatives", "shell_node_count", "mesh_triangles",
                "production_plate_area_coverage", "accounting_pass", "force_transfer_residual_fraction",
                "torque_transfer_residual_fraction", "equilibrium_residual_fraction", "shell_energy_residual_fraction",
                "target_demand_quantiles_n_per_m",
            )}
            for factor, shell in shell_by_factor.items()
        },
        "capacity_event_screen_by_load_case": capacity_event_screens,
        "capacity_event_screen_rule": {
            "domains": ["whole selected plate", "eligible continental area on selected plate"],
            "largest_connected_overcapacity_fraction": 0.25,
            "persistence_myr": 10.0,
            "each_registered_load_stage_is_held_for_at_least_10_myr": True,
        },
        "capacity_assignment": {
            "continental_cells": "Gate 2 continental capacity anchors",
            "noncontinental_cells": "Gate 2 mature-oceanic capacity anchors",
            "oceanic_age_state": "not represented; mature-oceanic is the fixed scenario",
        },
        "damage_law": {
            "epsilon_onset": EPSILON_ONSET,
            "epsilon_fail": EPSILON_FAIL,
            "viscosity_pa_s": DAMAGE_VISCOSITY_PA_S,
            "residual_strength_fraction": RESIDUAL_STRENGTH_FRACTION,
            "provenance": "Gate 3a local-feedback diagnostic scenario values; not calibrated production parameters",
        },
        "split_rule": {
            "basis": "failed local cells must disconnect the exact production cell-adjacency graph",
            "minimum_daughter_area_fraction": MIN_DAUGHTER_AREA_FRACTION,
            "minimum_intact_area_fraction": MIN_INTACT_AREA_FRACTION,
        },
        "force_components_at_native_load_n": {
            "slab_pull": np.asarray(field.slab_pull_n).sum(axis=0).tolist(),
            "midline_current_basal_drag": np.asarray(midline_force).sum(axis=0).tolist(),
            "total_production_field": np.asarray(field.total_n).sum(axis=0).tolist(),
        },
        "capacity_cases": outcomes,
        "elapsed_wall_seconds": time.perf_counter() - started,
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.out)
    print(json.dumps({
        "capacity_cases": {
            key: {"completed_split_candidate": value["completed_split_candidate"],
                  "split_time_myr": value["split_time_myr"],
                  "cracked_cell_count_at_end": value["cracked_cell_count_at_end"]}
            for key, value in result["capacity_cases"].items()
        },
        "elapsed_wall_seconds": result["elapsed_wall_seconds"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
