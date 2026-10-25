#!/usr/bin/env python3
"""Finite-strain local-feedback prototype with nonlocal softening and rate loading.

The strip uses total-Lagrangian linear triangles, a compressible neo-Hookean
elastic response, multiplicative inelastic deformation, and a current-configuration
Gaussian nonlocal damage driver. This is a synthetic diagnostic, not production
planet physics.

Usage from the repository root:
    backend/.venv/bin/python bin/debug/prototype_gate3_finite_strain.py \
        --output analysis/issue487/gate3-finite-strain-prototype.json
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.sparse.linalg import spsolve

from prototype_gate3_local_feedback import (
    EPSILON_FAIL,
    EPSILON_ONSET,
    ETA_PA_S,
    LENGTH_M,
    POISSON_RATIO,
    RESIDUAL_FRACTION,
    SECONDS_PER_MYR,
    SEED_WIDTH_M,
    THICKNESS_M,
    WIDTH_M,
    YIELD_PA,
    build_mesh,
    damage_strength,
    principal_stress_and_direction,
)

YOUNG_MODULUS_PA = 70.0e9
LAME_MU_PA = YOUNG_MODULUS_PA / (2.0 * (1.0 + POISSON_RATIO))
LAME_LAMBDA_PA = (
    YOUNG_MODULUS_PA * POISSON_RATIO
    / ((1.0 + POISSON_RATIO) * (1.0 - 2.0 * POISSON_RATIO))
)
REFERENCE_RATE_PER_YEAR = 2.0e-8
REFERENCE_LENGTH_M = 25.0e3
MAX_LOG_EXTENSION = 0.30
HISTORY_INTERVAL_MYR = 0.1
NEWTON_TOLERANCE = 1.0e-8
NEWTON_MAX_ITERATIONS = 30


def reference_gradients(mesh):
    gradients = np.empty((len(mesh.triangles), 3, 2))
    areas = np.empty(len(mesh.triangles))
    for index, (b_matrix, triangle) in enumerate(zip(mesh.b_matrices, mesh.triangles)):
        gradients[index, :, 0] = b_matrix[0, (0, 2, 4)]
        gradients[index, :, 1] = b_matrix[1, (1, 3, 5)]
        points = mesh.coordinates[triangle]
        areas[index] = abs(
            (points[1, 0] - points[0, 0]) * (points[2, 1] - points[0, 1])
            - (points[2, 0] - points[0, 0]) * (points[1, 1] - points[0, 1])
        ) / 2.0
    return gradients, areas


def boundary_node_mask(mesh):
    edge_counts = Counter(
        tuple(sorted((int(first), int(second))))
        for triangle in mesh.triangles
        for first, second in ((triangle[0], triangle[1]), (triangle[1], triangle[2]), (triangle[2], triangle[0]))
    )
    mask = np.zeros(len(mesh.coordinates), dtype=bool)
    for (first, second), count in edge_counts.items():
        if count == 1:
            mask[first] = True
            mask[second] = True
    return mask


def nonlocal_operator(current_centers, current_areas, length_m):
    if length_m <= 0.0:
        raise ValueError("nonlocal length must be positive")
    tree = cKDTree(current_centers)
    pairs = tree.sparse_distance_matrix(
        tree, 3.0 * length_m, output_type="coo_matrix"
    )
    weights = np.exp(-0.5 * pairs.data**2 / length_m**2) * current_areas[pairs.col]
    row_totals = np.bincount(pairs.row, weights=weights, minlength=len(current_centers))
    operator = coo_matrix(
        (weights / row_totals[pairs.row], (pairs.row, pairs.col)),
        shape=(len(current_centers), len(current_centers)),
    ).tocsr()
    if np.max(np.abs(operator @ np.ones(len(current_centers)) - 1.0)) > 1.0e-12:
        raise ValueError("nonlocal operator failed its constant-field invariant")
    return operator


def element_response(u_nodes, gradients, fp):
    identity = np.eye(2)
    displacement_gradient = np.einsum("eai,eaj->eij", u_nodes, gradients)
    deformation = identity + displacement_gradient
    det_f = np.linalg.det(deformation)
    fp_inverse = np.linalg.inv(fp)
    elastic_f = deformation @ fp_inverse
    det_fe = np.linalg.det(elastic_f)
    if np.any(det_f <= 0.0) or np.any(det_fe <= 0.0):
        raise ValueError("deformation gradient lost positive orientation")

    inverse_elastic_transpose = np.linalg.inv(elastic_f).transpose(0, 2, 1)
    log_j = np.log(det_fe)
    elastic_pk = (
        LAME_MU_PA * (elastic_f - inverse_elastic_transpose)
        + LAME_LAMBDA_PA * log_j[:, None, None] * inverse_elastic_transpose
    )
    first_piola = elastic_pk @ fp_inverse.transpose(0, 2, 1)
    cauchy = first_piola @ deformation.transpose(0, 2, 1) / det_f[:, None, None]
    return deformation, elastic_f, fp_inverse, inverse_elastic_transpose, log_j, first_piola, cauchy


def assemble_internal_force(u, fp, mesh, gradients, areas, with_tangent):
    triangles = mesh.triangles
    u_nodes = u[triangles]
    (deformation, _elastic_f, fp_inverse, inv_fe_t, log_j,
     first_piola, cauchy) = element_response(u_nodes, gradients, fp)
    internal_element = THICKNESS_M * areas[:, None, None] * np.einsum(
        "eaj,eij->eai", gradients, first_piola
    )
    internal = np.zeros(2 * len(mesh.coordinates))
    dofs = np.stack((2 * triangles, 2 * triangles + 1), axis=2).reshape(len(triangles), 6)
    np.add.at(internal, dofs.ravel(), internal_element.reshape(-1))
    if not with_tangent:
        return internal, None, deformation, cauchy

    element_count = len(triangles)
    q = inv_fe_t
    inv_fe = q.transpose(0, 2, 1)
    a_tensor = np.empty((element_count, 2, 2, 2, 2))
    for displacement_axis in range(2):
        for reference_axis in range(2):
            d_fe = np.zeros_like(fp_inverse)
            d_fe[:, displacement_axis, :] = fp_inverse[:, reference_axis, :]
            d_log_j = np.trace(inv_fe @ d_fe, axis1=1, axis2=2)
            d_q = -q @ d_fe.transpose(0, 2, 1) @ q
            d_elastic_pk = (
                LAME_MU_PA * (d_fe - d_q)
                + LAME_LAMBDA_PA * (
                    d_log_j[:, None, None] * q + log_j[:, None, None] * d_q
                )
            )
            a_tensor[:, :, :, displacement_axis, reference_axis] = (
                d_elastic_pk @ fp_inverse.transpose(0, 2, 1)
            )
    element_tangent = THICKNESS_M * areas[:, None, None, None, None] * np.einsum(
        "eaI,eiIjJ,ebJ->eaibj", gradients, a_tensor, gradients, optimize=True
    )
    element_tangent = element_tangent.reshape(element_count, 6, 6)
    rows = np.repeat(dofs, 6, axis=1).ravel()
    columns = np.tile(dofs, (1, 6)).ravel()
    values = element_tangent.ravel()
    tangent = coo_matrix(
        (values, (rows, columns)), shape=(len(internal), len(internal))
    ).tocsr()
    return internal, tangent, deformation, cauchy


def solve_equilibrium(
    mesh, gradients, areas, fp, initial_u, log_extension, external_nodal_force=None,
    prescribed_boundary_displacement=None, additional_fixed_y_nodes=None,
    cohesive_stiffness=None,
):
    u = initial_u.copy()
    node_count = len(mesh.coordinates)
    if external_nodal_force is None:
        external_flat = np.zeros(2 * node_count)
    else:
        external_nodal_force = np.asarray(external_nodal_force, dtype=float)
        if external_nodal_force.shape != (node_count, 2):
            raise ValueError("external nodal force must have shape (node_count, 2)")
        external_flat = external_nodal_force.reshape(-1)
    constrained: dict[int, float] = {}
    if prescribed_boundary_displacement is None:
        for node in mesh.left_nodes:
            constrained[2 * node] = 0.0
        constrained[2 * mesh.bottom_left_node + 1] = 0.0
        right_displacement = LENGTH_M * math.expm1(log_extension)
        for node in mesh.right_nodes:
            constrained[2 * node] = right_displacement
        for node in additional_fixed_y_nodes or ():
            constrained[2 * node + 1] = 0.0
    else:
        prescribed_boundary_displacement = np.asarray(prescribed_boundary_displacement, dtype=float)
        if prescribed_boundary_displacement.shape != mesh.coordinates.shape:
            raise ValueError("prescribed boundary displacement must match mesh coordinate shape")
        boundary = boundary_node_mask(mesh)
        for node in np.flatnonzero(boundary):
            for axis in range(2):
                value = prescribed_boundary_displacement[node, axis]
                # NaN marks a free boundary component, allowing a cutout to leave
                # interior crack faces traction-free (#487).
                if np.isfinite(value):
                    constrained[2 * node + axis] = value
    constrained_dofs = np.asarray(sorted(constrained), dtype=int)
    constrained_values = np.asarray([constrained[dof] for dof in constrained_dofs])
    free = np.ones(2 * node_count, dtype=bool)
    free[constrained_dofs] = False
    u_flat = u.reshape(-1)
    u_flat[constrained_dofs] = constrained_values
    force_scale = max(YIELD_PA * THICKNESS_M * WIDTH_M, 1.0)
    cohesive_flat = None if cohesive_stiffness is None else cohesive_stiffness.tocsr()

    for iteration in range(1, NEWTON_MAX_ITERATIONS + 1):
        internal, tangent, deformation, cauchy = assemble_internal_force(
            u, fp, mesh, gradients, areas, with_tangent=True
        )
        if cohesive_flat is not None:
            internal = internal + cohesive_flat @ u.reshape(-1)
            tangent = tangent + cohesive_flat
        residual = (internal - external_flat)[free]
        residual_norm = float(np.linalg.norm(residual))
        relative_residual = residual_norm / force_scale
        if relative_residual <= NEWTON_TOLERANCE:
            return u, deformation, cauchy, iteration, relative_residual

        delta = np.zeros_like(u_flat)
        delta[free] = spsolve(tangent[free][:, free], -residual)
        step_length = 1.0
        accepted = False
        while step_length >= 1.0 / 128.0:
            trial_flat = u_flat + step_length * delta
            trial_flat[constrained_dofs] = constrained_values
            try:
                trial_internal, _, _, _ = assemble_internal_force(
                    trial_flat.reshape((-1, 2)), fp, mesh, gradients, areas,
                    with_tangent=False,
                )
                if cohesive_flat is not None:
                    trial_internal = trial_internal + cohesive_flat @ trial_flat
                trial_norm = float(
                    np.linalg.norm((trial_internal - external_flat)[free])
                )
            except (ValueError, np.linalg.LinAlgError):
                trial_norm = math.inf
            if trial_norm < residual_norm:
                u_flat = trial_flat
                accepted = True
                break
            step_length *= 0.5
        if not accepted:
            raise RuntimeError(
                f"Newton line search failed at iteration {iteration}; "
                f"relative residual {relative_residual:.3e}"
            )
        u = u_flat.reshape((-1, 2))
    raise RuntimeError(
        f"Newton solve exceeded {NEWTON_MAX_ITERATIONS} iterations; "
        f"relative residual {residual_norm / force_scale:.3e}"
    )


def update_plastic_deformation(fp, direction, increment):
    updated = fp.copy()
    active = increment > 0.0
    if not np.any(active):
        return updated
    normals = direction[active]
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    identity = np.eye(2)
    deviator = 2.0 * np.einsum("ei,ej->eij", normals, normals) - identity
    values = increment[active]
    # Isotropic in-plane plastic log extension/compression: det(Fp) stays one.
    exponential = (
        np.cosh(values)[:, None, None] * identity
        + np.sinh(values)[:, None, None] * deviator
    )
    updated[active] = exponential @ updated[active]
    return updated


def verify_finite_strain_mechanics():
    """Run constitutive, objectivity, small-strain, and tangent patch checks."""
    mesh = build_mesh(2, 1)
    gradients, areas = reference_gradients(mesh)
    triangles = mesh.triangles
    identity = np.eye(2)
    fp = np.repeat(identity[None, :, :], len(triangles), axis=0)

    def response_for_deformation(deformation_gradient, plastic=fp):
        displacement = mesh.coordinates @ (deformation_gradient - identity).T
        return element_response(displacement[triangles], gradients, plastic)[-1]

    zero_stress = response_for_deformation(identity)
    zero_stress_error_pa = float(np.max(np.abs(zero_stress)))

    stretch = np.array([[1.17, 0.0], [0.0, 0.93]])
    stretch_stress = response_for_deformation(stretch)
    jacobian = np.linalg.det(stretch)
    left_cauchy_green = stretch @ stretch.T
    expected_stress = (
        LAME_MU_PA * (left_cauchy_green - identity)
        + LAME_LAMBDA_PA * math.log(jacobian) * identity
    ) / jacobian
    affine_stretch_relative_error = float(
        np.linalg.norm(stretch_stress[0] - expected_stress)
        / max(np.linalg.norm(expected_stress), 1.0)
    )

    angle = 0.73
    rotation = np.array([
        [math.cos(angle), -math.sin(angle)],
        [math.sin(angle), math.cos(angle)],
    ])
    rotated_identity_stress = response_for_deformation(rotation)
    rotated_stretch_stress = response_for_deformation(rotation @ stretch)
    expected_rotated_stress = rotation @ stretch_stress[0] @ rotation.T
    elastic_rotation_objectivity_relative_error = float(
        np.linalg.norm(rotated_stretch_stress[0] - expected_rotated_stress)
        / max(np.linalg.norm(stretch_stress[0]), 1.0)
    )
    rigid_rotation_relative_stress = float(
        np.max(np.abs(rotated_identity_stress))
        / max(LAME_MU_PA, 1.0)
    )

    small_strain_errors = []
    strain_shape = np.array([[1.0, 0.4], [0.4, -0.3]])
    linear_stress_shape = (
        2.0 * LAME_MU_PA * strain_shape
        + LAME_LAMBDA_PA * np.trace(strain_shape) * identity
    )
    for magnitude in (1.0e-4, 1.0e-5):
        finite_stress = response_for_deformation(identity + magnitude * strain_shape)[0]
        linear_stress = magnitude * linear_stress_shape
        small_strain_errors.append(float(
            np.linalg.norm(finite_stress - linear_stress)
            / max(np.linalg.norm(linear_stress), 1.0)
        ))

    rng = np.random.default_rng(487)
    element_count = len(triangles)
    flow_directions = rng.normal(size=(element_count, 2))
    plastic = update_plastic_deformation(
        np.repeat(identity[None, :, :], element_count, axis=0),
        flow_directions,
        np.full(element_count, 0.03),
    )
    plastic_determinant_error = float(np.max(np.abs(np.linalg.det(plastic) - 1.0)))
    base_deformation = np.array([[1.015, 0.012], [-0.006, 0.992]])
    displacement = mesh.coordinates @ (base_deformation - identity).T
    displacement += rng.normal(scale=5.0, size=displacement.shape)
    _, tangent, _, _ = assemble_internal_force(
        displacement, plastic, mesh, gradients, areas, with_tangent=True
    )
    direction = rng.normal(size=displacement.shape)
    direction /= np.linalg.norm(direction)
    perturbation = 0.1
    force_plus = assemble_internal_force(
        displacement + perturbation * direction, plastic, mesh, gradients, areas,
        with_tangent=False,
    )[0]
    force_minus = assemble_internal_force(
        displacement - perturbation * direction, plastic, mesh, gradients, areas,
        with_tangent=False,
    )[0]
    finite_difference_tangent = (force_plus - force_minus) / (2.0 * perturbation)
    analytic_tangent = tangent @ direction.reshape(-1)
    tangent_relative_error = float(
        np.linalg.norm(finite_difference_tangent - analytic_tangent)
        / max(np.linalg.norm(finite_difference_tangent), 1.0)
    )

    damage_alpha = np.array([EPSILON_ONSET + 0.1])
    damage_stress = np.array([1.5 * YIELD_PA])
    damage_stress_tangent = np.array([2.0 * YIELD_PA])
    dt_over_eta = 0.1 / YIELD_PA
    damage_beta = (
        YIELD_PA * (1.0 - RESIDUAL_FRACTION)
        / (EPSILON_FAIL - EPSILON_ONSET)
    )
    initial_damage_strength = damage_strength(damage_alpha)
    expected_damage_increment = (
        dt_over_eta * (damage_stress - initial_damage_strength)
        / (1.0 + dt_over_eta * (damage_stress_tangent - damage_beta))
    )
    scalar_operator = coo_matrix(([1.0], ([0], [0])), shape=(1, 1)).tocsr()
    updated_alpha, _ = advance_damage_implicit(
        damage_alpha,
        np.zeros(1),
        damage_stress,
        damage_stress_tangent,
        dt_over_eta * ETA_PA_S,
        scalar_operator,
    )
    implicit_damage_update_relative_error = float(
        np.max(np.abs(
            updated_alpha - damage_alpha - expected_damage_increment
        ))
    )

    checks = {
        "undeformed_max_abs_cauchy_stress_pa": zero_stress_error_pa,
        "affine_stretch_relative_error": affine_stretch_relative_error,
        "rigid_rotation_relative_stress": rigid_rotation_relative_stress,
        "elastic_rotation_objectivity_relative_error": elastic_rotation_objectivity_relative_error,
        "small_strain_relative_errors_at_1e-4_and_1e-5": small_strain_errors,
        "plastic_determinant_max_abs_error": plastic_determinant_error,
        "consistent_tangent_directional_relative_error": tangent_relative_error,
        "implicit_damage_update_abs_error": implicit_damage_update_relative_error,
    }
    checks["passed"] = bool(
        zero_stress_error_pa < 1.0e-5
        and affine_stretch_relative_error < 1.0e-12
        and rigid_rotation_relative_stress < 1.0e-12
        and elastic_rotation_objectivity_relative_error < 1.0e-12
        and small_strain_errors[1] < small_strain_errors[0]
        and small_strain_errors[1] < 1.0e-4
        and plastic_determinant_error < 1.0e-12
        and tangent_relative_error < 1.0e-5
        and implicit_damage_update_relative_error < 1.0e-10
    )
    return checks


def projected_stress_tangent(u, fp, direction, stress, mesh, gradients):
    """Finite-difference tangent of tensile Cauchy stress along the flow direction."""
    perturbation = 1.0e-6
    fp_perturbed = update_plastic_deformation(
        fp, direction, np.full(len(fp), perturbation)
    )
    perturbed_stress = element_response(
        u[mesh.triangles], gradients, fp_perturbed
    )[-1]
    perturbed_maximum, _ = principal_stress_and_direction(perturbed_stress)
    tangent = (stress - perturbed_maximum) / perturbation
    return np.maximum(tangent, 1.0e6)


def advance_damage_implicit(alpha, post_failure, stress, tangent, duration_s, operator):
    """Implicit local Maxwell update with the nonlocal softening field iterated."""
    diagonal = operator.diagonal()
    beta = YIELD_PA * (1.0 - RESIDUAL_FRACTION) / (EPSILON_FAIL - EPSILON_ONSET)
    alpha_increment = np.zeros_like(alpha)
    dt_over_eta = duration_s / ETA_PA_S
    for _ in range(100):
        nonlocal_alpha = np.clip(operator @ (alpha + alpha_increment), 0.0, EPSILON_FAIL)
        next_increment = np.zeros_like(alpha_increment)
        active = alpha < EPSILON_FAIL
        indices = np.flatnonzero(active)
        without_self = nonlocal_alpha[indices] - diagonal[indices] * alpha_increment[indices]
        base_strength = damage_strength(without_self)
        softening_slope = np.where(
            (without_self > EPSILON_ONSET) & (without_self < EPSILON_FAIL), beta, 0.0
        )
        denominator = 1.0 + dt_over_eta * np.maximum(
            tangent[indices] - softening_slope * diagonal[indices], 1.0e6
        )
        candidate = dt_over_eta * np.maximum(stress[indices] - base_strength, 0.0) / denominator
        next_increment[indices] = np.minimum(
            EPSILON_FAIL - alpha[indices], candidate
        )
        if np.max(np.abs(next_increment - alpha_increment)) <= 1.0e-11:
            alpha_increment = next_increment
            break
        alpha_increment = next_increment
    else:
        raise RuntimeError("nonlocal implicit damage iteration did not converge")

    post_increment = np.zeros_like(post_failure)
    failed_before = alpha >= EPSILON_FAIL
    post_overstress = np.maximum(
        stress[failed_before] - RESIDUAL_FRACTION * YIELD_PA, 0.0
    )
    post_increment[failed_before] = (
        dt_over_eta * post_overstress
        / (1.0 + dt_over_eta * tangent[failed_before])
    )
    return alpha + alpha_increment, post_failure + post_increment


def failed_component_count(mesh, failed):
    edge_owner: dict[tuple[int, int], int] = {}
    rows: list[int] = []
    cols: list[int] = []
    for index, triangle in enumerate(mesh.triangles):
        if not failed[index]:
            continue
        for a, b in ((triangle[0], triangle[1]), (triangle[1], triangle[2]),
                     (triangle[2], triangle[0])):
            key = tuple(sorted((int(a), int(b))))
            if key in edge_owner:
                neighbor = edge_owner[key]
                rows.extend((index, neighbor))
                cols.extend((neighbor, index))
            else:
                edge_owner[key] = index
    graph = coo_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(failed), len(failed))
    ).tocsr()
    return connected_components(graph[failed][:, failed], directed=False)[0]


def initialize_case_fields(
    mesh, seed_level, seed_centers_x_m, seed_length_m,
    initial_plastic_deformation, initial_damage_field, diagnostic_regions,
):
    element_count = len(mesh.triangles)
    fp = np.repeat(np.eye(2)[None, :, :], element_count, axis=0)
    if initial_plastic_deformation is not None:
        initial_plastic_deformation = np.asarray(initial_plastic_deformation, dtype=float)
        if initial_plastic_deformation.shape != fp.shape:
            raise ValueError("initial plastic deformation must have shape (element_count, 2, 2)")
        fp = initial_plastic_deformation.copy()

    alpha = np.zeros(element_count)
    seed_centers = seed_centers_x_m or (LENGTH_M / 2.0,)
    seed_mask = np.zeros(element_count, dtype=bool)
    for center_x in seed_centers:
        seed_mask |= (
            (np.abs(mesh.element_centers[:, 1] - WIDTH_M / 2.0) <= SEED_WIDTH_M / 2.0)
            & (np.abs(mesh.element_centers[:, 0] - center_x) <= seed_length_m / 2.0)
        )
    if initial_damage_field is not None:
        initial_damage_field = np.asarray(initial_damage_field, dtype=float)
        if initial_damage_field.shape != (element_count,):
            raise ValueError("initial damage field must have shape (element_count,)")
        if np.any(initial_damage_field < 0.0) or np.any(initial_damage_field > EPSILON_FAIL):
            raise ValueError("initial damage field must be within [0, epsilon_fail]")
        alpha = initial_damage_field.copy()
        seed_mask = alpha > 0.0
    elif seed_level > 0.0:
        alpha[seed_mask] = seed_level

    regions = {}
    for name, region in (diagnostic_regions or {}).items():
        region = np.asarray(region, dtype=bool)
        if region.shape != (element_count,):
            raise ValueError(f"diagnostic region {name!r} must have shape (element_count,)")
        regions[name] = region
    return fp, alpha, seed_mask, seed_centers, regions


def support_constraint_mask(mesh, prescribed_boundary_displacement, additional_fixed_y_nodes):
    constrained = np.zeros(2 * len(mesh.coordinates), dtype=bool)
    if prescribed_boundary_displacement is not None:
        boundary = boundary_node_mask(mesh)
        prescribed = np.asarray(prescribed_boundary_displacement, dtype=float)
        if prescribed.shape != mesh.coordinates.shape:
            raise ValueError("prescribed boundary displacement must match mesh coordinate shape")
        constrained = (np.isfinite(prescribed) & boundary[:, None]).reshape(-1)
    else:
        constrained[2 * mesh.left_nodes] = True
        constrained[2 * mesh.bottom_left_node + 1] = True
        constrained[2 * mesh.right_nodes] = True
        if additional_fixed_y_nodes is not None:
            constrained[2 * np.asarray(additional_fixed_y_nodes, dtype=int) + 1] = True
    return constrained


def case_stop_reason(tensile_stress, process_zone_stress, deformation, fp,
                    regularize_tensile_driver):
    """Classify non-finite, over-limit or inverted states after equilibrium."""
    if not np.all(np.isfinite(tensile_stress)):
        return "non-finite Cauchy stress"
    if np.max(process_zone_stress) > 5.0 * YIELD_PA:
        return (
            "maximum process-zone averaged tensile stress exceeded 5x initial strength"
            if regularize_tensile_driver
            else "maximum tensile Cauchy stress exceeded 5x initial strength"
        )
    if np.any(np.linalg.det(deformation) <= 0.0) or np.any(np.linalg.det(fp) <= 0.0):
        return "deformation lost positive orientation"
    return None


def run_case(
    *, nx, ny, dt_myr, duration_myr, rate_per_year, length_m,
    seed_level=0.0, seed_centers_x_m=None, seed_length_m=50.0e3,
    external_nodal_force=None, external_load_metadata=None,
    prescribed_boundary_displacement=None, mesh=None,
    include_nodal_support_reaction=False,
    initial_plastic_deformation=None, initial_damage_field=None,
    diagnostic_regions=None, additional_fixed_y_nodes=None,
    regularize_tensile_driver=False,
    cohesive_stiffness=None,
):
    requested_log_extension = rate_per_year * duration_myr * 1.0e6
    if requested_log_extension > MAX_LOG_EXTENSION:
        raise ValueError(
            f"requested log extension {requested_log_extension:.4g} exceeds "
            f"prototype limit {MAX_LOG_EXTENSION}"
        )
    mesh = build_mesh(nx, ny) if mesh is None else mesh
    if external_nodal_force is None:
        external_nodal_force = np.zeros_like(mesh.coordinates)
    else:
        external_nodal_force = np.asarray(external_nodal_force, dtype=float)
        if external_nodal_force.shape != mesh.coordinates.shape:
            raise ValueError("external nodal force must match mesh coordinate shape")
    if prescribed_boundary_displacement is not None:
        prescribed_boundary_displacement = np.asarray(prescribed_boundary_displacement, dtype=float)
        if prescribed_boundary_displacement.shape != mesh.coordinates.shape:
            raise ValueError("prescribed boundary displacement must match mesh coordinate shape")
    gradients, reference_areas = reference_gradients(mesh)
    element_count = len(mesh.triangles)
    fp, alpha, seed_mask, seed_centers, diagnostic_regions = initialize_case_fields(
        mesh, seed_level, seed_centers_x_m, seed_length_m,
        initial_plastic_deformation, initial_damage_field, diagnostic_regions,
    )
    post_failure = np.zeros(element_count)

    u = np.zeros_like(mesh.coordinates)
    steps = math.ceil(duration_myr / dt_myr)
    actual_dt_myr = duration_myr / steps
    sample_every = max(1, round(HISTORY_INTERVAL_MYR / actual_dt_myr))
    history = []
    previous_stress = np.zeros(element_count)
    total_newton_iterations = 0
    stop_reason = "completed requested interval"
    completed_time = 0.0

    for step in range(steps):
        time_myr = (step + 1) * actual_dt_myr
        log_extension = rate_per_year * time_myr * 1.0e6
        boundary_displacement = None
        if prescribed_boundary_displacement is not None:
            boundary_displacement = prescribed_boundary_displacement.copy()
            boundary_displacement[:, 0] += (
                mesh.coordinates[:, 0] / LENGTH_M * LENGTH_M * math.expm1(log_extension)
            )
        alpha_start = alpha.copy()
        post_start = post_failure.copy()
        try:
            u, deformation, cauchy, iterations, residual = solve_equilibrium(
                mesh, gradients, reference_areas, fp, u, log_extension,
                external_nodal_force, boundary_displacement,
                additional_fixed_y_nodes, cohesive_stiffness,
            )
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            stop_reason = f"finite-strain equilibrium failed: {error}"
            break
        total_newton_iterations += iterations
        tensile_stress, direction = principal_stress_and_direction(cauchy)

        current_coordinates = mesh.coordinates + u
        current_centers = current_coordinates[mesh.triangles].mean(axis=1)
        det_f = np.linalg.det(deformation)
        current_areas = reference_areas * det_f
        operator = nonlocal_operator(current_centers, current_areas, length_m)
        tangent = projected_stress_tangent(
            u, fp, direction, tensile_stress, mesh, gradients
        )
        damage_stress = operator @ tensile_stress if regularize_tensile_driver else tensile_stress
        damage_tangent = operator @ tangent if regularize_tensile_driver else tangent
        alpha, post_failure = advance_damage_implicit(
            alpha_start, post_failure, damage_stress, damage_tangent,
            actual_dt_myr * SECONDS_PER_MYR, operator,
        )
        post_increment = post_failure - post_start
        plastic_increment = alpha - alpha_start + post_increment
        fp = update_plastic_deformation(fp, direction, plastic_increment)

        try:
            u, deformation, cauchy, iterations, residual = solve_equilibrium(
                mesh, gradients, reference_areas, fp, u, log_extension,
                external_nodal_force, boundary_displacement,
                additional_fixed_y_nodes, cohesive_stiffness,
            )
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            stop_reason = f"post-inelastic finite-strain equilibrium failed: {error}"
            break
        total_newton_iterations += iterations
        tensile_stress, _ = principal_stress_and_direction(cauchy)
        final_coordinates = mesh.coordinates + u
        final_centers = final_coordinates[mesh.triangles].mean(axis=1)
        final_areas = reference_areas * np.linalg.det(deformation)
        final_operator = nonlocal_operator(final_centers, final_areas, length_m)
        process_zone_stress = (
            final_operator @ tensile_stress if regularize_tensile_driver else tensile_stress
        )
        completed_time = time_myr
        newly_failed = (alpha_start < EPSILON_FAIL) & (alpha >= EPSILON_FAIL)
        step_stop_reason = case_stop_reason(
            tensile_stress, process_zone_stress, deformation, fp,
            regularize_tensile_driver,
        )
        is_stop = step_stop_reason is not None
        if is_stop:
            stop_reason = step_stop_reason

        if (step + 1) % sample_every == 0 or is_stop or np.any(newly_failed) or step + 1 == steps:
            final_coordinates = mesh.coordinates + u
            if prescribed_boundary_displacement is None:
                boundary_log_extension = math.log(
                    np.mean(final_coordinates[mesh.right_nodes, 0]) / LENGTH_M
                )
            else:
                # The conforming cutout has no structured right edge; measure the mean
                # boundary x-displacement gradient instead (#487).
                boundary = boundary_node_mask(mesh)
                if prescribed_boundary_displacement is not None:
                    prescribed = np.asarray(prescribed_boundary_displacement, dtype=float)
                    boundary &= np.any(np.isfinite(prescribed), axis=1)
                boundary_x = mesh.coordinates[boundary, 0]
                boundary_ux = u[boundary, 0]
                centered_x = boundary_x - np.mean(boundary_x)
                denominator = float(np.dot(centered_x, centered_x))
                extension_slope = (
                    float(np.dot(centered_x, boundary_ux - np.mean(boundary_ux)) / denominator)
                    if denominator > 0.0 else 0.0
                )
                boundary_log_extension = math.log1p(extension_slope) if extension_slope > -1.0 else float("nan")
            current_centers = final_coordinates[mesh.triangles].mean(axis=1)
            current_areas = reference_areas * np.linalg.det(deformation)
            operator_after = nonlocal_operator(current_centers, current_areas, length_m)
            history_entry = {
                "time_myr": time_myr,
                "prescribed_log_extension": log_extension,
                "measured_mean_log_extension": boundary_log_extension,
                "failed_elements": int(np.count_nonzero(alpha >= EPSILON_FAIL)),
                "newly_failed_elements": int(np.count_nonzero(newly_failed)),
                "max_damage_strain": float(np.max(alpha)),
                "max_nonlocal_damage_driver": float(np.max(np.clip(operator_after @ alpha, 0.0, EPSILON_FAIL))),
                "peak_max_principal_cauchy_stress_pa": float(np.max(tensile_stress)),
                "peak_process_zone_averaged_tensile_stress_pa": float(
                    np.max(process_zone_stress)
                ),
                "median_max_principal_cauchy_stress_pa": float(np.median(tensile_stress)),
                "peak_stress_change_pa": float(np.max(tensile_stress - previous_stress)),
                "minimum_det_deformation": float(np.min(np.linalg.det(deformation))),
                "minimum_det_plastic_deformation": float(np.min(np.linalg.det(fp))),
                "newton_iterations_this_step": iterations,
                "relative_equilibrium_residual": residual,
            }
            for region_name, region in diagnostic_regions.items():
                if np.any(region):
                    history_entry[f"{region_name}_max_damage_strain"] = float(np.max(alpha[region]))
                    history_entry[f"{region_name}_peak_tensile_stress_pa"] = float(
                        np.max(tensile_stress[region])
                    )
                    history_entry[
                        f"{region_name}_peak_process_zone_averaged_tensile_stress_pa"
                    ] = float(np.max(process_zone_stress[region]))
                    history_entry[f"{region_name}_mean_sigma_xx_pa"] = float(
                        np.mean(cauchy[region, 0, 0])
                    )
            history.append(history_entry)
        previous_stress = tensile_stress
        if is_stop:
            break

    failed = alpha >= EPSILON_FAIL
    failed_centers = mesh.element_centers[failed]
    seed_failed = int(np.count_nonzero(failed & seed_mask))
    outside_failed = int(np.count_nonzero(failed & ~seed_mask))
    external_flat = external_nodal_force.reshape(-1)
    final_internal, _, _, _ = assemble_internal_force(
        u, fp, mesh, gradients, reference_areas, with_tangent=False
    )
    if cohesive_stiffness is not None:
        final_internal = final_internal + cohesive_stiffness.tocsr() @ u.reshape(-1)
    support_reaction = final_internal - external_flat
    constrained = support_constraint_mask(
        mesh, prescribed_boundary_displacement, additional_fixed_y_nodes
    )
    reaction_flat = np.where(constrained, support_reaction, 0.0)
    final_coordinates = mesh.coordinates + u
    external_resultant = external_nodal_force.sum(axis=0)
    reaction_resultant = reaction_flat.reshape((-1, 2)).sum(axis=0)
    external_moment = float(np.sum(
        final_coordinates[:, 0] * external_nodal_force[:, 1]
        - final_coordinates[:, 1] * external_nodal_force[:, 0]
    ))
    reaction_nodes = reaction_flat.reshape((-1, 2))
    reaction_moment = float(np.sum(
        final_coordinates[:, 0] * reaction_nodes[:, 1]
        - final_coordinates[:, 1] * reaction_nodes[:, 0]
    ))
    force_scale = max(
        float(np.sum(np.linalg.norm(external_nodal_force, axis=1)))
        + float(np.sum(np.linalg.norm(reaction_nodes, axis=1))),
        1.0,
    )
    moment_scale = max(
        float(np.sum(np.abs(
            final_coordinates[:, 0] * external_nodal_force[:, 1]
            - final_coordinates[:, 1] * external_nodal_force[:, 0]
        )))
        + float(np.sum(np.abs(
            final_coordinates[:, 0] * reaction_nodes[:, 1]
            - final_coordinates[:, 1] * reaction_nodes[:, 0]
        ))),
        1.0,
    )
    result = {
        "nx": nx,
        "ny": ny,
        "element_spacing_x_km": float(np.ptp(mesh.coordinates[:, 0]) / nx / 1000.0),
        "element_spacing_y_km": float(np.ptp(mesh.coordinates[:, 1]) / ny / 1000.0),
        "nonlocal_length_m": length_m,
        "tensile_driver_regularized": regularize_tensile_driver,
        "dt_myr": actual_dt_myr,
        "requested_duration_myr": duration_myr,
        "simulated_duration_myr": completed_time,
        "extension_log_rate_per_year": rate_per_year,
        "prescribed_log_extension": rate_per_year * completed_time * 1.0e6,
        "initial_seed_level": seed_level,
        "initial_seed_elements": int(np.count_nonzero(seed_mask)) if (seed_level > 0.0 or initial_damage_field is not None) else 0,
        "initial_seed_area_km2": float(np.sum(reference_areas[seed_mask]) / 1.0e6) if (seed_level > 0.0 or initial_damage_field is not None) else 0.0,
        "initial_seed_centers_km": [center / 1000.0 for center in seed_centers] if seed_level > 0 else [],
        "max_damage_strain_seed": float(np.max(alpha[seed_mask])) if np.any(seed_mask) else 0.0,
        "max_damage_strain_outside_seed": float(np.max(alpha[~seed_mask])) if np.any(~seed_mask) else 0.0,
        "seed_failed_elements": seed_failed,
        "failed_outside_seed_elements": outside_failed,
        "failed_elements_total": int(np.count_nonzero(failed)),
        "failed_components": int(failed_component_count(mesh, failed)),
        "failed_x_span_km": [float(failed_centers[:, 0].min() / 1000.0), float(failed_centers[:, 0].max() / 1000.0)] if len(failed_centers) else None,
        "failed_y_span_km": [float(failed_centers[:, 1].min() / 1000.0), float(failed_centers[:, 1].max() / 1000.0)] if len(failed_centers) else None,
        "total_newton_iterations": total_newton_iterations,
        "stop_reason": stop_reason,
        "external_load_metadata": external_load_metadata,
        "boundary_loading_metadata": (
            "whole-plate displacement transfer plus prescribed extension-rate field"
            if prescribed_boundary_displacement is not None and rate_per_year != 0.0
            else "whole-plate displacement transfer"
            if prescribed_boundary_displacement is not None else None
        ),
        "external_resultant_force_n": external_resultant.tolist(),
        "support_reaction_resultant_force_n": reaction_resultant.tolist(),
        "force_balance_relative_residual": float(
            np.linalg.norm(external_resultant + reaction_resultant) / force_scale
        ),
        "external_moment_about_strip_origin_n_m": external_moment,
        "support_reaction_moment_about_strip_origin_n_m": reaction_moment,
        "moment_balance_relative_residual": abs(external_moment + reaction_moment) / moment_scale,
        "history": history,
    }
    if include_nodal_support_reaction:
        result["boundary_support_reaction_n"] = reaction_nodes.tolist()
    return result


def run_experiment(output_path=None):
    cases = []
    requested_duration_myr = 2.5

    def add_case(name, **parameters):
        result = run_case(**parameters)
        result["case"] = name
        cases.append(result)
        if output_path:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(json.dumps(
                experiment_metadata(cases), sort_keys=True, separators=(",", ":")
            ) + "\n", encoding="utf-8")
        print(f"completed {len(cases)} cases: {name}", flush=True)

    for rate in (1.0e-8, REFERENCE_RATE_PER_YEAR, 4.0e-8):
        duration = min(requested_duration_myr, MAX_LOG_EXTENSION / (rate * 1.0e6))
        add_case("intact_rate_control", nx=32, ny=10, dt_myr=0.02,
                 duration_myr=duration, rate_per_year=rate,
                 length_m=REFERENCE_LENGTH_M)
    add_case("seeded_extension_rate", nx=32, ny=10, dt_myr=0.02,
             duration_myr=requested_duration_myr,
             rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=REFERENCE_LENGTH_M,
             seed_level=1.15, seed_length_m=50.0e3)
    add_case("two_seed_linkage", nx=32, ny=10, dt_myr=0.02,
             duration_myr=requested_duration_myr,
             rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=REFERENCE_LENGTH_M,
             seed_level=1.15, seed_length_m=50.0e3,
             seed_centers_x_m=(140.0e3, 260.0e3))
    for dt in (0.04, 0.02, 0.01):
        add_case("seeded_timestep_refinement", nx=32, ny=10, dt_myr=dt,
                 duration_myr=requested_duration_myr,
                 rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=REFERENCE_LENGTH_M,
                 seed_level=1.15, seed_length_m=50.0e3)
    for nx, ny in ((16, 5), (24, 8)):
        add_case("seeded_mesh_refinement", nx=nx, ny=ny, dt_myr=0.02,
                 duration_myr=requested_duration_myr,
                 rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=REFERENCE_LENGTH_M,
                 seed_level=1.15, seed_length_m=50.0e3)
    for length in (12.5e3, 25.0e3, 50.0e3):
        add_case("nonlocal_length_sensitivity", nx=32, ny=10, dt_myr=0.02,
                 duration_myr=requested_duration_myr,
                 rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=length,
                 seed_level=1.15, seed_length_m=50.0e3)
    for seed_level in (1.15, 1.24):
        add_case("seed_proximity_sensitivity", nx=32, ny=10, dt_myr=0.02,
                 duration_myr=requested_duration_myr,
                 rate_per_year=REFERENCE_RATE_PER_YEAR, length_m=REFERENCE_LENGTH_M,
                 seed_level=seed_level, seed_length_m=50.0e3)
    return experiment_metadata(cases)


def experiment_metadata(cases):
    return {
        "title": "RFC 0010 Gate 3 finite-strain nonlocal softening and extension-rate prototype",
        "constitutive_basis": {
            "elastic_energy": "compressible neo-Hookean, plane strain, total-Lagrangian formulation",
            "deformation_split": "F = Fe Fp; elastic stress from Fe and inelastic flow updates Fp multiplicatively",
            "elastic_modulus_pa": YOUNG_MODULUS_PA,
            "poisson_ratio": POISSON_RATIO,
            "plastic_flow": "isochoric in-plane Hencky extension/compression aligned with maximum tensile Cauchy stress",
            "damage_law": {
                "epsilon_onset": EPSILON_ONSET,
                "epsilon_fail": EPSILON_FAIL,
                "eta_rift_pa_s": ETA_PA_S,
                "residual_strength_fraction": RESIDUAL_FRACTION,
                "initial_average_strength_pa": YIELD_PA,
            },
        },
        "regularization": {
            "type": "area-weighted Gaussian integral nonlocal average of accumulated damage strain",
            "reference_length_m": REFERENCE_LENGTH_M,
            "support_radius_in_lengths": 3.0,
            "configuration": "current deformed geometry with current element-area weights",
        },
        "boundary_condition": {
            "type": "prescribed constant logarithmic extension rate; left edge fixed and right edge ramps",
            "reference_rate_per_year": REFERENCE_RATE_PER_YEAR,
            "rate_sensitivity_per_year": [1.0e-8, 4.0e-8],
            "duration_cap": "prescribed log extension <= 0.30",
        },
        "prototype_assumptions": [
            "2-D plane-strain material with no thickness stretch or out-of-plane stress relaxation",
            "reference mesh is total-Lagrangian; finite deformation gradients and current coordinates are updated at each solve",
            "compressible neo-Hookean elasticity and multiplicative inelastic deformation are verification choices, not calibrated lithosphere rheology",
            "damage driving strain is nonlocal in the current configuration; the selected length is a design sensitivity, not a geological calibration",
            "maximum-principal Cauchy stress drives positive inelastic Hencky extension using the frozen Gate 3 damage law",
            "constitutive update is staggered: solve nonlinear equilibrium, use an implicit local Maxwell update with a finite-difference projected stress tangent, then solve equilibrium again at the same boundary displacement",
            "the uniform strip and initial weak patches are synthetic; no world state or production topology changes are made",
        ],
        "numerical_controls": {
            "newton_relative_residual_tolerance": NEWTON_TOLERANCE,
            "newton_max_iterations": NEWTON_MAX_ITERATIONS,
            "maximum_log_extension": MAX_LOG_EXTENSION,
            "maximum_tensile_stress_pa": 5.0 * YIELD_PA,
            "history_sampling_myr": HISTORY_INTERVAL_MYR,
        },
        "cases": cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output:
        run_experiment(args.output)
    else:
        print(json.dumps(run_experiment(), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
