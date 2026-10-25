#!/usr/bin/env python3
"""Prototype local stress redistribution from strain-softening inelastic extension.

This is a diagnostic-only, small-strain plane-stress strip. It couples the frozen
Gate 3 scalar damage law to an incremental elastic equilibrium solve through local
inelastic eigenstrain. It does not modify production rifting or plate topology.

Usage from the repository root:
    backend/.venv/bin/python bin/debug/prototype_gate3_local_feedback.py \
        --output analysis/issue487/gate3-local-feedback-prototype.json
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

SECONDS_PER_MYR = 1.0e6 * 365.25 * 24.0 * 3600.0
YOUNG_MODULUS_PA = 70.0e9
POISSON_RATIO = 0.25
THICKNESS_M = 100.0e3
WIDTH_M = 100.0e3
LENGTH_M = 400.0e3
YIELD_PA = 62.8e6
EPSILON_ONSET = 0.25
EPSILON_FAIL = 1.25
ETA_PA_S = 9.0e21
RESIDUAL_FRACTION = 1.0 / 3.0
SEED_WIDTH_M = 20.0e3


@dataclass(frozen=True)
class Mesh:
    coordinates: np.ndarray
    triangles: np.ndarray
    b_matrices: np.ndarray
    stiffness: object
    eigenstrain_load_matrices: np.ndarray
    left_nodes: np.ndarray
    right_nodes: np.ndarray
    bottom_left_node: int
    element_centers: np.ndarray
    nx: int
    ny: int


def plane_stress_matrix() -> np.ndarray:
    scale = YOUNG_MODULUS_PA / (1.0 - POISSON_RATIO**2)
    return scale * np.array([
        [1.0, POISSON_RATIO, 0.0],
        [POISSON_RATIO, 1.0, 0.0],
        [0.0, 0.0, (1.0 - POISSON_RATIO) / 2.0],
    ])


def triangle_b_matrix(coordinates: np.ndarray) -> tuple[np.ndarray, float]:
    jacobian = np.column_stack((coordinates[1] - coordinates[0], coordinates[2] - coordinates[0]))
    twice_area = float(np.linalg.det(jacobian))
    if twice_area <= 0.0:
        raise ValueError("mesh triangles must be counterclockwise and nondegenerate")
    inverse_jacobian = np.linalg.inv(jacobian)
    physical_gradients = inverse_jacobian.T
    gradients = np.column_stack((-physical_gradients.sum(axis=1), physical_gradients))
    b_matrix = np.zeros((3, 6))
    b_matrix[0, 0::2] = gradients[0]
    b_matrix[1, 1::2] = gradients[1]
    b_matrix[2, 0::2] = gradients[1]
    b_matrix[2, 1::2] = gradients[0]
    return b_matrix, twice_area / 2.0


def build_mesh(nx: int, ny: int) -> Mesh:
    x_values = np.linspace(0.0, LENGTH_M, nx + 1)
    y_values = np.linspace(0.0, WIDTH_M, ny + 1)
    coordinates = np.array([(x, y) for y in y_values for x in x_values], dtype=float)
    node_id = lambda i, j: j * (nx + 1) + i
    triangles = []
    for j in range(ny):
        for i in range(nx):
            lower_left = node_id(i, j)
            lower_right = node_id(i + 1, j)
            upper_left = node_id(i, j + 1)
            upper_right = node_id(i + 1, j + 1)
            triangles.extend(((lower_left, lower_right, upper_right),
                              (lower_left, upper_right, upper_left)))
    triangles_array = np.asarray(triangles, dtype=np.int64)
    constitutive = plane_stress_matrix()
    element_dofs = []
    b_matrices = []
    eigen_loads = []
    centers = []
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for triangle in triangles_array:
        b_matrix, area = triangle_b_matrix(coordinates[triangle])
        dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        stiffness = THICKNESS_M * area * b_matrix.T @ constitutive @ b_matrix
        load_matrix = THICKNESS_M * area * b_matrix.T @ constitutive
        rows.extend(np.repeat(dofs, 6).tolist())
        columns.extend(np.tile(dofs, 6).tolist())
        values.extend(stiffness.ravel().tolist())
        element_dofs.append(dofs)
        b_matrices.append(b_matrix)
        eigen_loads.append(load_matrix)
        centers.append(coordinates[triangle].mean(axis=0))
    stiffness_global = coo_matrix(
        (values, (rows, columns)), shape=(2 * len(coordinates), 2 * len(coordinates))
    ).tocsr()
    return Mesh(
        coordinates=coordinates,
        triangles=triangles_array,
        b_matrices=np.asarray(b_matrices),
        stiffness=stiffness_global,
        eigenstrain_load_matrices=np.asarray(eigen_loads),
        left_nodes=np.asarray([node_id(0, j) for j in range(ny + 1)], dtype=int),
        right_nodes=np.asarray([node_id(nx, j) for j in range(ny + 1)], dtype=int),
        bottom_left_node=node_id(0, 0),
        element_centers=np.asarray(centers),
        nx=nx,
        ny=ny,
    )


def damage_strength(epsilon: np.ndarray) -> np.ndarray:
    damage = np.clip(
        (epsilon - EPSILON_ONSET) / (EPSILON_FAIL - EPSILON_ONSET), 0.0, 1.0
    )
    return YIELD_PA * (1.0 - damage * (1.0 - RESIDUAL_FRACTION))


def advance_damage_scalar(
    epsilon: np.ndarray, stress_pa: np.ndarray, duration_s: float
) -> tuple[np.ndarray, np.ndarray]:
    """Exact scalar update at fixed maximum-principal stress.

    Returns the capped damage strain and additional post-failure plastic extension.
    Post-failure plastic extension is separate from the capped damage counter.
    """
    updated = epsilon.copy()
    post_failure = np.zeros_like(epsilon)
    for index, (old, stress) in enumerate(zip(epsilon, stress_pa)):
        if stress <= 0.0:
            continue
        remaining = duration_s
        current = float(old)
        if current < EPSILON_ONSET:
            overstress = stress - YIELD_PA
            if overstress <= 0.0:
                continue
            to_onset = (EPSILON_ONSET - current) * ETA_PA_S / overstress
            if to_onset >= remaining:
                updated[index] = current + overstress * remaining / ETA_PA_S
                continue
            current = EPSILON_ONSET
            remaining -= to_onset
        if current < EPSILON_FAIL and remaining > 0.0:
            span = EPSILON_FAIL - EPSILON_ONSET
            beta = YIELD_PA * (1.0 - RESIDUAL_FRACTION) / span
            x = max(0.0, current - EPSILON_ONSET)
            current_strength = YIELD_PA - beta * x
            overstress = stress - current_strength
            if overstress > 0.0:
                if beta == 0.0:
                    x_new = x + overstress * remaining / ETA_PA_S
                    time_to_failure = math.inf
                else:
                    offset = (stress - YIELD_PA) / beta
                    x_new = (x + offset) * math.exp(
                        min(beta * remaining / ETA_PA_S, 700.0)
                    ) - offset
                    if x_new >= span:
                        time_to_failure = ETA_PA_S / beta * math.log(
                            (span + offset) / (x + offset)
                        )
                    else:
                        time_to_failure = math.inf
                if x_new >= span:
                    remaining_after_failure = max(0.0, remaining - time_to_failure)
                    post_failure[index] = max(
                        0.0, stress - RESIDUAL_FRACTION * YIELD_PA
                    ) * remaining_after_failure / ETA_PA_S
                    current = EPSILON_FAIL
                else:
                    current = EPSILON_ONSET + max(0.0, x_new)
        elif current >= EPSILON_FAIL and remaining > 0.0:
            post_failure[index] = max(
                0.0, stress - RESIDUAL_FRACTION * YIELD_PA
            ) * remaining / ETA_PA_S
        updated[index] = min(EPSILON_FAIL, current)
    return updated, post_failure


def advance_displacement_hold(
    epsilon: np.ndarray,
    stress_pa: np.ndarray,
    directions: np.ndarray,
    duration_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate local viscoplastic relaxation at fixed total strain.

    The element tangent is projected onto the maximum-tension flow direction. This
    is an exact local relaxation update for that fixed strain, followed by a global
    equilibrium re-solve; spatially coupled accuracy is assessed by refinement.
    """
    updated = epsilon.copy()
    post_failure = np.zeros_like(epsilon)
    constitutive = plane_stress_matrix()
    beta = YIELD_PA * (1.0 - RESIDUAL_FRACTION) / (EPSILON_FAIL - EPSILON_ONSET)
    for index, (old, stress, direction) in enumerate(zip(epsilon, stress_pa, directions)):
        if stress <= 0.0:
            continue
        flow = np.array([direction[0] ** 2, direction[1] ** 2,
                         2.0 * direction[0] * direction[1]])
        tangent = float(flow @ constitutive @ flow)
        remaining = duration_s
        current = float(old)

        def relaxed_increment(overstress: float, modulus: float, time_s: float) -> float:
            if overstress <= 0.0 or time_s <= 0.0:
                return 0.0
            exponent = modulus * time_s / ETA_PA_S
            return overstress / modulus * (-math.expm1(-min(exponent, 700.0)))

        def time_for_increment(overstress: float, modulus: float, increment: float) -> float:
            fraction = modulus * increment / overstress
            if fraction >= 1.0:
                return math.inf
            return -ETA_PA_S / modulus * math.log1p(-fraction)

        if current < EPSILON_ONSET:
            overstress = stress - YIELD_PA
            to_onset_strain = EPSILON_ONSET - current
            attainable = overstress / tangent if overstress > 0.0 else 0.0
            if attainable <= to_onset_strain:
                updated[index] = current + relaxed_increment(overstress, tangent, remaining)
                continue
            time_to_onset = time_for_increment(overstress, tangent, to_onset_strain)
            if time_to_onset >= remaining:
                updated[index] = current + relaxed_increment(overstress, tangent, remaining)
                continue
            current = EPSILON_ONSET
            remaining -= time_to_onset
            stress -= tangent * to_onset_strain

        if current < EPSILON_FAIL:
            softening_modulus = tangent - beta
            if softening_modulus <= 0.0:
                raise ValueError("softening slope exceeds elastic tangent")
            strength = YIELD_PA - beta * (current - EPSILON_ONSET)
            overstress = stress - strength
            to_failure = EPSILON_FAIL - current
            attainable = overstress / softening_modulus if overstress > 0.0 else 0.0
            if attainable <= to_failure:
                current += relaxed_increment(overstress, softening_modulus, remaining)
                updated[index] = current
                continue
            time_to_failure = time_for_increment(
                overstress, softening_modulus, to_failure
            )
            if time_to_failure >= remaining:
                updated[index] = current + relaxed_increment(
                    overstress, softening_modulus, remaining
                )
                continue
            current = EPSILON_FAIL
            remaining -= time_to_failure
            stress -= tangent * to_failure

        residual_overstress = stress - RESIDUAL_FRACTION * YIELD_PA
        post_failure[index] = relaxed_increment(residual_overstress, tangent, remaining)
        updated[index] = min(EPSILON_FAIL, current)
    return updated, post_failure


def principal_stress_and_direction(stress: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh(stress)
    max_index = np.argmax(values, axis=1)
    maximum = values[np.arange(len(values)), max_index]
    direction = vectors[np.arange(len(values)), :, max_index]
    return maximum, direction


def solve_equilibrium(
    mesh: Mesh,
    plastic_strain: np.ndarray,
    load_ratio: float,
    boundary: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    node_count = len(mesh.coordinates)
    rhs = np.zeros(2 * node_count)
    for element_index, (matrix, triangle) in enumerate(
        zip(mesh.eigenstrain_load_matrices, mesh.triangles)
    ):
        dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        rhs[dofs] += matrix @ plastic_strain[element_index]

    constrained: dict[int, float] = {}
    for node in mesh.left_nodes:
        constrained[2 * node] = 0.0
    constrained[2 * mesh.bottom_left_node + 1] = 0.0
    target_stress = load_ratio * YIELD_PA
    if boundary == "displacement":
        imposed_strain = target_stress / YOUNG_MODULUS_PA
        displacement = imposed_strain * LENGTH_M
        for node in mesh.right_nodes:
            constrained[2 * node] = displacement
    elif boundary == "force":
        total_force = target_stress * THICKNESS_M * WIDTH_M
        edge_weights = np.ones(len(mesh.right_nodes), dtype=float)
        edge_weights[[0, -1]] = 0.5
        edge_weights /= mesh.ny
        rhs[2 * mesh.right_nodes] += total_force * edge_weights
    else:
        raise ValueError(f"unknown boundary control {boundary!r}")

    constrained_dofs = np.asarray(sorted(constrained), dtype=int)
    constrained_values = np.asarray([constrained[dof] for dof in constrained_dofs])
    free = np.ones(2 * node_count, dtype=bool)
    free[constrained_dofs] = False
    displacement = np.zeros(2 * node_count)
    displacement[constrained_dofs] = constrained_values
    free_rhs = rhs[free] - mesh.stiffness[free][:, constrained_dofs] @ constrained_values
    displacement[free] = spsolve(mesh.stiffness[free][:, free], free_rhs)

    constitutive = plane_stress_matrix()
    element_stress = np.empty((len(mesh.triangles), 3))
    for index, (b_matrix, triangle) in enumerate(zip(mesh.b_matrices, mesh.triangles)):
        dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        total_strain = b_matrix @ displacement[dofs]
        element_stress[index] = constitutive @ (total_strain - plastic_strain[index])
    stress_tensor = np.zeros((len(element_stress), 2, 2))
    stress_tensor[:, 0, 0] = element_stress[:, 0]
    stress_tensor[:, 1, 1] = element_stress[:, 1]
    stress_tensor[:, 0, 1] = element_stress[:, 2]
    stress_tensor[:, 1, 0] = element_stress[:, 2]
    maximum, direction = principal_stress_and_direction(stress_tensor)
    return displacement.reshape(node_count, 2), maximum, direction


def flow_strain_tensor(direction: np.ndarray, magnitude: np.ndarray) -> np.ndarray:
    nx, ny = direction[:, 0], direction[:, 1]
    return magnitude[:, None] * np.column_stack((nx * nx, ny * ny, 2.0 * nx * ny))


def run_case(
    *, nx: int, ny: int, dt_myr: float, duration_myr: float,
    load_ratio: float, boundary: str, seed_level: float = 0.0,
    seed_length_m: float = 0.0,
    seed_centers_x_m: tuple[float, ...] | None = None,
) -> dict[str, object]:
    mesh = build_mesh(nx, ny)
    eps = np.zeros(len(mesh.triangles))
    seed_centers = seed_centers_x_m or (LENGTH_M / 2.0,)
    if seed_level > 0.0:
        in_seed = np.zeros(len(mesh.triangles), dtype=bool)
        for center_x in seed_centers:
            in_seed |= (
                (np.abs(mesh.element_centers[:, 1] - WIDTH_M / 2.0) <= SEED_WIDTH_M / 2.0)
                & (np.abs(mesh.element_centers[:, 0] - center_x) <= seed_length_m / 2.0)
            )
        eps[in_seed] = seed_level
    plastic = np.zeros((len(mesh.triangles), 3))
    steps = math.ceil(duration_myr / dt_myr)
    actual_dt = duration_myr / steps
    sample_every = max(1, round(0.25 / actual_dt))
    snapshots = []
    previous_stress = np.zeros(len(mesh.triangles))
    stop_reason = "completed requested interval"
    for step in range(steps):
        dt_s = actual_dt * SECONDS_PER_MYR
        eps_start = eps.copy()
        plastic_start = plastic.copy()
        _, trial_stress, direction = solve_equilibrium(
            mesh, plastic_start, load_ratio, boundary
        )
        if boundary == "force":
            eps, post_failure = advance_damage_scalar(eps_start, trial_stress, dt_s)
        else:
            eps, post_failure = advance_displacement_hold(
                eps_start, trial_stress, direction, dt_s
            )
        plastic_increment = eps - eps_start + post_failure
        plastic = plastic_start + flow_strain_tensor(direction, plastic_increment)
        _, max_stress, _ = solve_equilibrium(mesh, plastic, load_ratio, boundary)
        newly_failed = (eps_start < EPSILON_FAIL) & (eps >= EPSILON_FAIL)
        snapshot = {
            "time_myr": (step + 1) * actual_dt,
            "failed_elements": int(np.count_nonzero(eps >= EPSILON_FAIL)),
            "newly_failed_elements": int(np.count_nonzero(newly_failed)),
            "max_damage_strain": float(np.max(eps)),
            "peak_principal_stress_pa": float(np.max(max_stress)),
            "median_principal_stress_pa": float(np.median(max_stress)),
            "peak_stress_change_pa": float(np.max(max_stress - previous_stress)),
            "max_abs_inelastic_eigenstrain": float(np.max(np.abs(plastic))),
        }
        previous_stress = max_stress
        if not np.all(np.isfinite(max_stress)):
            stop_reason = "non-finite stress"
            is_stop = True
        elif np.max(max_stress) > 5.0 * YIELD_PA:
            stop_reason = "stress exceeded 5x initial strength; local small-strain result is out of envelope"
            is_stop = True
        elif np.max(np.abs(plastic)) > 0.05:
            stop_reason = "inelastic eigenstrain exceeded 0.05; fixed-geometry small-strain result is out of envelope"
            is_stop = True
        else:
            is_stop = False
        if step % sample_every == 0 or is_stop or np.any(newly_failed):
            snapshots.append(snapshot)
        if is_stop:
            break

    failed = eps >= EPSILON_FAIL
    failed_centers = mesh.element_centers[failed]
    seed_mask = np.zeros(len(mesh.triangles), dtype=bool)
    for center_x in seed_centers:
        seed_mask |= (
            (np.abs(mesh.element_centers[:, 1] - WIDTH_M / 2.0) <= SEED_WIDTH_M / 2.0)
            & (np.abs(mesh.element_centers[:, 0] - center_x) <= seed_length_m / 2.0)
        )
    collar_failed = int(np.count_nonzero(failed & ~seed_mask))
    seed_failed = int(np.count_nonzero(failed & seed_mask))
    failed_x = failed_centers[:, 0] if len(failed_centers) else np.array([])
    return {
        "nx": nx,
        "ny": ny,
        "dt_myr": actual_dt,
        "duration_myr": duration_myr,
        "load_ratio_to_initial_strength": load_ratio,
        "boundary_control": "fixed_force" if boundary == "force" else "fixed_displacement_hold",
        "initial_seed_level": seed_level,
        "initial_seed_length_km": seed_length_m / 1000.0,
        "initial_seed_width_km": SEED_WIDTH_M / 1000.0 if seed_level > 0.0 else 0.0,
        "initial_seed_centers_km": [float(center / 1000.0) for center in seed_centers]
            if seed_level > 0.0 else [],
        "max_damage_strain_seed": float(np.max(eps[seed_mask])) if np.any(seed_mask) else 0.0,
        "max_damage_strain_outside_seed": float(np.max(eps[~seed_mask])) if np.any(~seed_mask) else 0.0,
        "failed_elements_total": int(np.count_nonzero(failed)),
        "failed_seed_elements": seed_failed,
        "failed_outside_seed_elements": collar_failed,
        "failed_x_span_km": (
            [float(failed_x.min() / 1000.0), float(failed_x.max() / 1000.0)]
            if len(failed_x) else None
        ),
        "peak_final_principal_stress_pa": snapshots[-1]["peak_principal_stress_pa"],
        "stop_reason": stop_reason,
        "simulated_duration_myr": snapshots[-1]["time_myr"],
        "integration": "exact local constant-stress update for force control; projected-tangent relaxation then equilibrium re-solve for fixed-displacement hold",
        "history": snapshots,
    }


def run_experiment() -> dict[str, object]:
    cases = []
    cases.append({
        "case": "subyield_tension",
        **run_case(nx=24, ny=6, dt_myr=0.05, duration_myr=8.0,
                   load_ratio=0.8, boundary="force"),
    })
    cases.append({
        "case": "compression",
        **run_case(nx=24, ny=6, dt_myr=0.05, duration_myr=8.0,
                   load_ratio=-0.8, boundary="force"),
    })
    # Uniform-load baseline: no seeded damage, 20% above initial strength.
    for boundary in ("force", "displacement"):
        cases.append({
            "case": "uniform_overstrength_baseline",
            **run_case(nx=24, ny=6, dt_myr=0.025, duration_myr=18.0,
                       load_ratio=1.2, boundary=boundary),
        })
    for dt in (0.0125, 0.00625):
        cases.append({
            "case": "uniform_overstrength_timestep_refinement",
            **run_case(nx=24, ny=6, dt_myr=dt, duration_myr=18.0,
                       load_ratio=1.2, boundary="force"),
        })
    # A finite, partially softened seed is loaded below intact strength. Any failure
    # outside it must arise from the equilibrium/inelastic feedback, not an imposed tip rule.
    for boundary in ("force", "displacement"):
        cases.append({
            "case": "seeded_subyield_tip_screen",
            **run_case(nx=32, ny=8, dt_myr=0.01, duration_myr=8.0,
                       load_ratio=0.8, boundary=boundary, seed_level=1.15,
                       seed_length_m=50.0e3),
        })
    cases.append({
        "case": "two_seed_linkage_screen",
        **run_case(nx=32, ny=8, dt_myr=0.01, duration_myr=8.0,
                   load_ratio=0.8, boundary="force", seed_level=1.15,
                   seed_length_m=50.0e3, seed_centers_x_m=(140.0e3, 260.0e3)),
    })
    # Timestep and mesh refinement at fixed force-controlled subyield loading.
    for nx, ny, dt in (
        (24, 6, 0.025),
        (24, 6, 0.0125),
        (24, 6, 0.00625),
        (48, 12, 0.0125),
    ):
        cases.append({
            "case": "seeded_refinement",
            **run_case(nx=nx, ny=ny, dt_myr=dt, duration_myr=8.0,
                       load_ratio=0.8, boundary="force", seed_level=1.15,
                       seed_length_m=50.0e3),
        })
    return {
        "title": "RFC 0010 Gate 3 local stress-redistribution prototype",
        "law": {
            "epsilon_onset": EPSILON_ONSET,
            "epsilon_fail": EPSILON_FAIL,
            "eta_rift_pa_s": ETA_PA_S,
            "residual_strength_fraction": RESIDUAL_FRACTION,
            "initial_average_strength_pa": YIELD_PA,
        },
        "prototype_assumptions": [
            "2-D small-strain plane-stress triangular finite elements with linear elasticity",
            "positive maximum-principal stress drives a viscoplastic extension eigenstrain",
            "damage strain follows the frozen two-threshold linear-softening law",
            "after full damage, additional plastic extension continues at residual strength while the damage counter stays capped",
            "yield strength does not scale Young's modulus",
            "force control holds remote traction; displacement control applies and holds the elastic displacement matching the initial traction",
            "the displacement-hold case is a one-time elastic loading followed by fixed end displacement, not a sustained extension-rate boundary condition",
            "seeded cells begin with epsilon_rift=1.15 but zero inelastic eigenstrain to isolate a pre-weakened local state",
            "time integration is staggered: solve equilibrium, update local strain using current stress, then re-solve equilibrium",
            "fixed geometry; no rift opening, fault friction, healing, or plate topology update",
        ],
        "diagnostic_envelope": {
            "max_principal_stress_pa": 5.0 * YIELD_PA,
            "max_abs_inelastic_eigenstrain": 0.05,
            "purpose": "stop outside the fixed-geometry small-strain prototype's declared range; these are numerical validity guards, not physical failure criteria",
            "history_sampling_myr": 0.25,
        },
        "material_scale": {
            "young_modulus_pa": YOUNG_MODULUS_PA,
            "poisson_ratio": POISSON_RATIO,
            "thickness_m": THICKNESS_M,
            "width_m": WIDTH_M,
            "length_m": LENGTH_M,
            "interpretation": "elastic regularization for a synthetic verification problem, not an inferred lithosphere property",
        },
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_experiment()
    encoded = json.dumps(result, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
