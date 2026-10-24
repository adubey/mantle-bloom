#!/usr/bin/env python3
"""Replay one Gate 2 event and compare source-band loads under alternate rigid-mode projections.

The paired solve audits mutual contact-drag force/torque closure and compares applying that
pair interaction with a bonded/no-slip solve where the shared interface carries the reaction.
This is not a calibrated finite-slip law.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
      bin/debug/replay_gate2_coupled_event.py --seed 7 --years 40 \
      --target-myr 38 --plate 0 --kind split \
      --output analysis/issue487/gate2-coupled-event-seed7-38myr.json

The source and its nearest force-contributing neighbor are remeshed together on one clipped
curved mesh. The shared displacement nodes impose a bonded/no-slip interface. The force-included
and tied/no-pair-contact-drag cases are sensitivity diagnostics; neither defines a geological
finite-slip law or models faults.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import replay_gate2_capacity as replay
from app import lithosphere, thin_shell, torque
from app.elevation_lines import line_spacing_rad
import prototype_mincut_split as mincut
import validate_thin_shell_plate_meshes as mesh_validation


def _neighbor_for(source, world):
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    neighbors = source.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(source, neighbors, spacing, reach,
                                                getattr(world, "boundary_search_cache", None))
    valid = inputs.neighbor_plate_id[inputs.neighbor_plate_id >= 0]
    if not len(valid):
        raise ValueError("source plate has no force-contributing neighbor")
    ids, counts = np.unique(valid, return_counts=True)
    neighbor_id = int(ids[np.argmax(counts)])
    return next(plate for plate in neighbors if plate.plate_id == neighbor_id), int(counts.max())


def _pair_contact_drag(world, plate, neighbor):
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    neighbors = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(
        plate, neighbors, spacing, reach, getattr(world, "boundary_search_cache", None)
    )
    _convergent, _divergent, _transform, contested = torque.classify_boundary_nodes(
        plate, neighbors, inputs, reach, getattr(world, "boundary_search_cache", None)
    )
    collision = torque.continental_contact_mask(plate, inputs, contested)
    craton = torque.contact_craton_m(plate, neighbors, inputs)
    jam = torque.contact_jam_factor(plate, world, neighbors, inputs, reach)
    force = torque.contact_drag_force_vectors_n(plate, inputs, collision, spacing, craton, jam)
    return force * (inputs.neighbor_plate_id == neighbor.plate_id)[:, None]


def _project_rigid_rotation(points_xyz, areas_m2, force_n):
    """Mirror the solver's area-weighted rigid-rotation force projection."""
    points = np.asarray(points_xyz, dtype=float)
    areas = np.asarray(areas_m2, dtype=float)
    force = np.asarray(force_n, dtype=float)
    first, second = thin_shell._tangent_basis(points)
    tangent = force - np.sum(force * points, axis=1)[:, None] * points
    force_dofs = np.column_stack((np.sum(tangent * first, axis=1),
                                  np.sum(tangent * second, axis=1))).ravel()
    modes = np.empty((2 * len(points), 3))
    xyz_m = points * lithosphere.PLANET_RADIUS_M
    for axis in range(3):
        omega = np.eye(3)[axis]
        velocity = np.cross(np.broadcast_to(omega, xyz_m.shape), xyz_m)
        modes[:, axis] = np.column_stack((np.sum(velocity * first, axis=1),
                                          np.sum(velocity * second, axis=1))).ravel()
    weights = np.repeat(areas, 2)
    coeff, *_ = np.linalg.lstsq(modes * np.sqrt(weights)[:, None],
                                force_dofs / np.sqrt(weights), rcond=None)
    rigid_dofs = weights * (modes @ coeff)
    projected_dofs = force_dofs - rigid_dofs
    rigid = rigid_dofs.reshape(-1, 2)[:, :1] * first + rigid_dofs.reshape(-1, 2)[:, 1:] * second
    projected = projected_dofs.reshape(-1, 2)[:, :1] * first + projected_dofs.reshape(-1, 2)[:, 1:] * second
    return projected, rigid


def _solve_pair(world, source, cut_points, cut_normal):
    neighbor, nearest_node_count = _neighbor_for(source, world)
    source_field = replay._event_force_field(world, source)
    neighbor_field = replay._event_force_field(world, neighbor)
    source_points = np.asarray(source_field.points_xyz, dtype=float)
    neighbor_points = np.asarray(neighbor_field.points_xyz, dtype=float)
    points = np.vstack((source_points, neighbor_points))
    proxy = SimpleNamespace(
        plate_id=f"{source.plate_id}+{neighbor.plate_id}",
        contains_batch=lambda probes: source.contains_batch(probes) | neighbor.contains_batch(probes),
    )
    triangles, rejected, total = mesh_validation._candidate_triangles(proxy, points)
    triangles = replay._nondegenerate_triangles(points, triangles)
    if len(triangles) == 0:
        raise ValueError("joint clipped mesh has no nondegenerate triangles")
    used, triangles = replay._remap_triangle_nodes(triangles)
    component_sizes = thin_shell._mesh_component_sizes(len(used), triangles)
    if len(component_sizes) != 1:
        raise ValueError(f"joint mesh is not edge-connected: {component_sizes.tolist()}")
    force = np.vstack((source_field.total_n, neighbor_field.total_n))[used]
    area = np.r_[source_field.areas_m2, neighbor_field.areas_m2][used]
    hc = np.r_[source.collect("crustal_thickness_m"), neighbor.collect("crustal_thickness_m")][used]
    hm = np.r_[source.collect("mantle_lithosphere_thickness_m"), neighbor.collect("mantle_lithosphere_thickness_m")][used]
    capacity = np.full(len(used), 1e10)
    shell_inputs = dict(
        points_xyz=points[used],
        areas_m2=area,
        thickness_m=np.maximum(hc + hm, 1.0),
        young_modulus_pa=replay.EFFECTIVE_MODULUS_PA,
        extensional_capacity_resultant_n_per_m=capacity,
        mesh_triangles=triangles,
    )
    joined = thin_shell.solve_thin_shell(applied_force_n=force, **shell_inputs)

    # Audit whether treating the pair as one rigid body suppresses relative plate modes.
    source_local = used < len(source_points)
    source_projected, source_rigid = _project_rigid_rotation(
        points[used][source_local], area[source_local], force[source_local]
    )
    neighbor_projected, neighbor_rigid = _project_rigid_rotation(
        points[used][~source_local], area[~source_local], force[~source_local]
    )
    platewise_force = np.empty_like(force)
    platewise_force[source_local] = source_projected
    platewise_force[~source_local] = neighbor_projected
    platewise_rigid_force = np.zeros_like(force)
    platewise_rigid_force[source_local] = source_rigid
    platewise_rigid_force[~source_local] = neighbor_rigid
    platewise = thin_shell.solve_thin_shell(applied_force_n=platewise_force, **shell_inputs)

    # In a tied/no-slip interface model, contact exchange is carried by the interface
    # constraint reaction. Do not also apply the torque engine's contact drag as a
    # pair-external load in this variant.
    source_pair_contact = _pair_contact_drag(world, source, neighbor)
    neighbor_pair_contact = _pair_contact_drag(world, neighbor, source)
    pair_contact_force = np.vstack((source_pair_contact, neighbor_pair_contact))
    tied_interface = thin_shell.solve_thin_shell(
        applied_force_n=force - pair_contact_force[used], **shell_inputs
    )

    # Recreate the same isolated-source clipped mesh and load used by the matched replay.
    source_triangles, _rejected, _total = mesh_validation._candidate_triangles(source, source_points)
    source_triangles = replay._nondegenerate_triangles(source_points, source_triangles)
    source_triangles, _repair = replay._repair_connected_source_mesh(source, source_points, source_triangles)
    source_used, source_triangles = replay._remap_triangle_nodes(source_triangles)
    isolated = thin_shell.solve_thin_shell(
        source_points[source_used], source_field.areas_m2[source_used], source_field.total_n[source_used],
        thickness_m=np.maximum(source.collect("crustal_thickness_m")[source_used]
                                + source.collect("mantle_lithosphere_thickness_m")[source_used], 1.0),
        young_modulus_pa=replay.EFFECTIVE_MODULUS_PA,
        extensional_capacity_resultant_n_per_m=np.full(len(source_used), 1e10),
        mesh_triangles=source_triangles,
    )

    spacing = line_spacing_rad(world.node_density)
    band_chord = 2.0 * np.sin(0.5 * replay.merge_split.FAILED_RIFT_BAND_MULT * spacing)
    graph_i, graph_j, _ = mincut._graph(cut_points, spacing)
    flow_mask = cut_points @ cut_normal > 0.0
    band = mincut._boundary_band(cut_points, graph_i, graph_j, flow_mask, band_chord)
    coupled_lookup = {int(original): i for i, original in enumerate(used)}
    paired = [(int(node), coupled_lookup[int(node)]) for node in source_used if int(node) in coupled_lookup and band[int(node)]]
    if len(paired) < 0.95 * max(int(band.sum()), 1):
        raise ValueError("coupled mesh does not cover 95% of source cut-band nodes")
    isolated_lookup = {int(original): i for i, original in enumerate(source_used)}
    paired = [(isolated_lookup[node], coupled_lookup[node]) for node in band.nonzero()[0]
              if node in isolated_lookup and node in coupled_lookup]
    if not paired:
        raise ValueError("no matched source cut-band nodes in both meshes")
    solo_values = np.asarray([isolated.max_tensile_resultant_n_per_m[a] for a, _ in paired])
    coupled_values = np.asarray([joined.max_tensile_resultant_n_per_m[b] for _, b in paired])
    tied_values = np.asarray([tied_interface.max_tensile_resultant_n_per_m[b] for _, b in paired])
    platewise_values = np.asarray([platewise.max_tensile_resultant_n_per_m[b] for _, b in paired])
    scale_n_per_m = max(float(np.mean(solo_values)), 1.0)
    joined_difference = np.abs(coupled_values - solo_values)
    tied_difference = np.abs(tied_values - solo_values)
    zero_tension = solo_values <= 1.0
    pair_force = np.vstack((source_field.total_n, neighbor_field.total_n))
    pair_points = np.vstack((source_points, neighbor_points))
    raw_force = pair_force.sum(axis=0)
    raw_torque = np.cross(pair_points * lithosphere.PLANET_RADIUS_M, pair_force).sum(axis=0)
    contact_torque = np.cross(
        pair_points * lithosphere.PLANET_RADIUS_M, pair_contact_force
    ).sum(axis=0)
    contact_scale = max(float(np.linalg.norm(pair_contact_force, axis=1).sum()), 1.0)
    contact_torque_scale = max(
        float(np.linalg.norm(np.cross(pair_points * lithosphere.PLANET_RADIUS_M, pair_contact_force), axis=1).sum()),
        1.0,
    )
    force_components = {}
    for name in ("slab_pull_n", "slab_drag_n", "ridge_push_n", "basal_drag_n", "contact_drag_n"):
        component = np.vstack((getattr(source_field, name), getattr(neighbor_field, name)))
        force_components[name.removesuffix("_n")] = {
            "source_resultant_n": np.sum(getattr(source_field, name), axis=0).tolist(),
            "neighbor_resultant_n": np.sum(getattr(neighbor_field, name), axis=0).tolist(),
            "pair_resultant_n": component.sum(axis=0).tolist(),
        }
    return {
        "neighbor_plate_id": int(neighbor.plate_id),
        "neighbor_node_count": int(len(neighbor_points)),
        "source_nodes_nearest_to_neighbor": nearest_node_count,
        "joint_mesh_nodes": int(len(used)),
        "joint_mesh_triangles": int(len(triangles)),
        "joint_mesh_rejected_candidate_triangles": int(rejected),
        "joint_mesh_candidate_triangles": int(total),
        "source_cut_band_node_count": int(band.sum()),
        "matched_source_cut_band_node_count": len(paired),
        "matched_cut_band_coverage": len(paired) / max(int(band.sum()), 1),
        "raw_pair_force_n": raw_force.tolist(),
        "raw_pair_torque_nm": raw_torque.tolist(),
        "force_component_resultants": force_components,
        "mutual_pair_contact_force_l1_sum_n": float(np.linalg.norm(pair_contact_force, axis=1).sum()),
        "source_pair_contact_force_resultant_n": source_pair_contact.sum(axis=0).tolist(),
        "neighbor_pair_contact_force_resultant_n": neighbor_pair_contact.sum(axis=0).tolist(),
        "pair_contact_force_action_reaction_closure_ratio": float(np.linalg.norm(pair_contact_force.sum(axis=0)) / contact_scale),
        "pair_contact_torque_action_reaction_closure_ratio": float(np.linalg.norm(contact_torque) / contact_torque_scale),
        "joined_solver_stop_code": int(joined.solver_stop_code),
        "joined_equilibrium_residual_over_load": joined.equilibrium_residual_n / max(float(np.linalg.norm(joined.equilibrated_force_n)), 1.0),
        "joined_rigid_rotation_force_fraction": float(np.linalg.norm(joined.rigid_rotation_force_n) / max(np.linalg.norm(force), 1.0)),
        "joined_rigid_rotation_removed_torque_nm": np.cross(points[used] * lithosphere.PLANET_RADIUS_M, joined.rigid_rotation_force_n).sum(axis=0).tolist(),
        "sourcewise_rigid_rotation_force_fraction_of_raw_pair_load": float(np.linalg.norm(source_rigid) / max(float(np.linalg.norm(force)), 1.0)),
        "neighborwise_rigid_rotation_force_fraction_of_raw_pair_load": float(np.linalg.norm(neighbor_rigid) / max(float(np.linalg.norm(force)), 1.0)),
        "platewise_total_rigid_rotation_force_fraction_of_raw_pair_load": float(np.linalg.norm(platewise_rigid_force) / max(float(np.linalg.norm(force)), 1.0)),
        "pairwise_projection_removed_force_fraction_of_raw_pair_load": float(np.linalg.norm(joined.rigid_rotation_force_n) / max(float(np.linalg.norm(force)), 1.0)),
        "platewise_projection_pair_cleanup_fraction_of_raw_pair_load": float(np.linalg.norm(platewise.rigid_rotation_force_n) / max(float(np.linalg.norm(force)), 1.0)),
        "projection_policy_removed_force_difference_fraction_of_raw_pair_load": float(np.linalg.norm(platewise_rigid_force - joined.rigid_rotation_force_n) / max(float(np.linalg.norm(force)), 1.0)),
        "platewise_solver_stop_code": int(platewise.solver_stop_code),
        "platewise_equilibrium_residual_over_load": platewise.equilibrium_residual_n / max(float(np.linalg.norm(platewise.equilibrated_force_n)), 1.0),
        "isolated_solver_stop_code": int(isolated.solver_stop_code),
        "isolated_equilibrium_residual_over_load": isolated.equilibrium_residual_n / max(float(np.linalg.norm(isolated.equilibrated_force_n)), 1.0),
        "isolated_rigid_rotation_force_fraction": float(np.linalg.norm(isolated.rigid_rotation_force_n) / max(np.linalg.norm(source_field.total_n[source_used]), 1.0)),
        "isolated_source_cut_band_n_ext_mean_n_per_m": float(solo_values.mean()),
        "joined_source_cut_band_n_ext_mean_n_per_m": float(coupled_values.mean()),
        "platewise_projection_source_cut_band_n_ext_mean_n_per_m": float(platewise_values.mean()),
        "platewise_vs_pairwise_projection_mean_abs_n_ext_difference_fraction_of_pairwise_mean": float(
            np.mean(np.abs(platewise_values - coupled_values)) / max(float(np.mean(coupled_values)), 1.0)
        ),
        "platewise_vs_pairwise_projection_p95_abs_n_ext_difference_fraction_of_pairwise_mean": float(
            np.quantile(np.abs(platewise_values - coupled_values), 0.95) / max(float(np.mean(coupled_values)), 1.0)
        ),
        "tied_interface_without_applied_contact_drag_source_cut_band_n_ext_mean_n_per_m": float(tied_values.mean()),
        "comparison_scale_mean_isolated_n_ext_n_per_m": scale_n_per_m,
        "joined_mean_absolute_difference_fraction_of_isolated_band_mean": float(np.mean(joined_difference) / scale_n_per_m),
        "tied_mean_absolute_difference_fraction_of_isolated_band_mean": float(np.mean(tied_difference) / scale_n_per_m),
        "joined_vs_tied_mean_absolute_difference_fraction_of_joined_band_mean": float(
            np.mean(np.abs(tied_values - coupled_values)) / max(float(np.mean(coupled_values)), 1.0)
        ),
        "joined_p95_absolute_difference_fraction_of_isolated_band_mean": float(np.quantile(joined_difference, 0.95) / scale_n_per_m),
        "source_band_fraction_zero_to_positive_tension_joined": float(
            np.mean(zero_tension & (coupled_values > 1.0))
        ),
        "source_band_fraction_zero_to_positive_tension_tied": float(
            np.mean(zero_tension & (tied_values > 1.0))
        ),
        "interface_assumptions": {
            "joined_with_pair_contact_drag": "Shared displacement mesh plus each side's current contact drag for this same pair; equilibrium depends on those discretized tractions being action-reaction balanced.",
            "tied_without_pair_contact_drag": "Perfectly tied/no-slip displacement interface; mutual contact transfer is represented by the interface constraint reaction, not also by applied pair contact-drag forces. Contact loads from other neighbors remain applied.",
            "limitation": "The tied case is an idealized bonded endmember, not a calibrated tectonic contact law or an upper bound on extensional stress. Relative node positions and unequal front sampling may prevent exact force/torque cancellation.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, choices=(7, 11, 2024), required=True)
    parser.add_argument("--years", type=float, required=True)
    parser.add_argument("--target-myr", type=float, required=True)
    parser.add_argument("--plate", type=int, required=True)
    parser.add_argument("--kind", choices=("split", "failed_rift"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    original_score = replay._score_event
    captured = {}

    def score_event(world, plate, points, normal, kind, capacities, **kwargs):
        row = original_score(world, plate, points, normal, kind, capacities, **kwargs)
        if (row is not None and row.get("plate") == args.plate and row.get("kind") == args.kind
                and abs(row.get("myr", -1) - args.target_myr) < 1e-9):
            try:
                row["coupled_neighbor_diagnostic"] = _solve_pair(world, plate, points, normal)
            except Exception as error:
                row["coupled_neighbor_diagnostic"] = {"error": f"{type(error).__name__}: {error}"}
            captured.update(row["coupled_neighbor_diagnostic"])
        return row

    replay._score_event = score_event
    try:
        report = replay._run(args.seed, args.years * 1e6, args.output, capture_only=True)
    finally:
        replay._score_event = original_score
    report["coupled_event_target"] = {"myr": args.target_myr, "plate": args.plate, "kind": args.kind}
    report["coupled_event_found"] = bool(captured)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"seed": args.seed, "years": args.years, "target": report["coupled_event_target"], "found": report["coupled_event_found"], "diagnostic": captured}, indent=2))
    if not captured or "error" in captured:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
