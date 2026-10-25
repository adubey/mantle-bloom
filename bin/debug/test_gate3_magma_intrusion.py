#!/usr/bin/env python3
"""Screen a localized pressurized dike against Gate 3 failure mechanics.

The diagnostic builds an explicit internal crack, calibrates its wall pressure to
an observed Afar aperture, checks elastic closure/work, then compares matched
magma-off/on damage runs. It does not modify production rifting or resolve magma
flow, solidification, or dike recurrence.

Usage from the repository root:
    backend/.venv/bin/python bin/debug/test_gate3_magma_intrusion.py \
        --output analysis/issue487/gate3-magma-intrusion-review.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import prototype_gate3_finite_strain as mechanics
import prototype_gate3_local_feedback as linear_mesh


MESH_NX = 32
MESH_NY = 20
DOMAIN_LENGTH_M = mechanics.LENGTH_M
DOMAIN_WIDTH_M = 130.0e3
DX_M = DOMAIN_LENGTH_M / MESH_NX
DY_M = DOMAIN_WIDTH_M / MESH_NY
DIKE_LENGTH_M = 65.0e3
DIKE_MEAN_APERTURE_M = 4.0
DIKE_MAX_APERTURE_M = 8.0
DIKE_HALF_LENGTH_M = DIKE_LENGTH_M / 2.0
DIKE_X_M = DOMAIN_LENGTH_M / 2.0
DIKE_Y_M = DOMAIN_WIDTH_M / 2.0
RATE_PER_YEAR = 4.0e-3 / DOMAIN_LENGTH_M
DURATION_MYR = 2.5
TIMESTEP_MYR = 0.02
SEED_LEVEL = 1.249
TIP_SEED_LENGTH_M = 12.5e3
TIP_PRESSURE_TAPER_M = 12.5e3


def build_cracked_mesh():
    """Split the mesh nodes on a finite dike, leaving shared nodes at its tips."""
    linear_mesh.WIDTH_M = DOMAIN_WIDTH_M
    mechanics.WIDTH_M = DOMAIN_WIDTH_M
    base = mechanics.build_mesh(MESH_NX, MESH_NY)
    x_index = MESH_NX // 2
    y_start = round((DIKE_Y_M - DIKE_HALF_LENGTH_M) / DY_M)
    y_end = round((DIKE_Y_M + DIKE_HALF_LENGTH_M) / DY_M)
    if not np.isclose(y_start * DY_M, DIKE_Y_M - DIKE_HALF_LENGTH_M) or not np.isclose(
        y_end * DY_M, DIKE_Y_M + DIKE_HALF_LENGTH_M
    ):
        raise ValueError("mesh does not place dike tips on mesh nodes")
    coordinates = base.coordinates.tolist()
    original_by_row = {j: j * (MESH_NX + 1) + x_index for j in range(y_start, y_end + 1)}
    right_face_by_row = dict(original_by_row)
    for row in range(y_start + 1, y_end):
        original = original_by_row[row]
        right_face_by_row[row] = len(coordinates)
        coordinates.append(base.coordinates[original].tolist())
    coordinates = np.asarray(coordinates, dtype=float)

    triangles = base.triangles.copy()
    for element_index, triangle in enumerate(triangles):
        center = base.element_centers[element_index]
        if center[0] <= DIKE_X_M or not (
            DIKE_Y_M - DIKE_HALF_LENGTH_M <= center[1] <= DIKE_Y_M + DIKE_HALF_LENGTH_M
        ):
            continue
        for local_node, node_id in enumerate(triangle):
            point = base.coordinates[node_id]
            if not np.isclose(point[0], DIKE_X_M):
                continue
            row = round(point[1] / DY_M)
            if row in right_face_by_row:
                triangles[element_index, local_node] = right_face_by_row[row]

    # Duplicating face-node IDs changes connectivity, not reference coordinates, so
    # the base element gradients remain valid for every triangle.
    b_matrices = base.b_matrices.copy()
    centers = coordinates[triangles].mean(axis=1)

    mesh = SimpleNamespace(
        coordinates=coordinates,
        triangles=triangles,
        b_matrices=b_matrices,
        left_nodes=base.left_nodes,
        right_nodes=base.right_nodes,
        bottom_left_node=base.bottom_left_node,
        bottom_right_node=MESH_NX,
        element_centers=centers,
        nx=MESH_NX,
        ny=MESH_NY,
    )
    return mesh, base, original_by_row, right_face_by_row, y_start, y_end


def pressure_force(mesh, pressure_pa, right_face_by_row, y_start, y_end):
    """Apply equal-and-opposite magma pressure to the two internal crack faces."""
    force = np.zeros_like(mesh.coordinates)
    def taper(y_m):
        tip_distance = min(
            y_m - (DIKE_Y_M - DIKE_HALF_LENGTH_M),
            DIKE_Y_M + DIKE_HALF_LENGTH_M - y_m,
        )
        return min(1.0, max(0.0, tip_distance / TIP_PRESSURE_TAPER_M))

    for row in range(y_start, y_end):
        left_nodes = (row * (MESH_NX + 1) + MESH_NX // 2,
                      (row + 1) * (MESH_NX + 1) + MESH_NX // 2)
        right_nodes = (right_face_by_row[row], right_face_by_row[row + 1])
        pressure_0 = pressure_pa * taper(row * DY_M)
        pressure_1 = pressure_pa * taper((row + 1) * DY_M)
        nodal_0 = mechanics.THICKNESS_M * DY_M * (2.0 * pressure_0 + pressure_1) / 6.0
        nodal_1 = mechanics.THICKNESS_M * DY_M * (pressure_0 + 2.0 * pressure_1) / 6.0
        for index, magnitude in enumerate((nodal_0, nodal_1)):
            force[left_nodes[index], 0] -= magnitude
            force[right_nodes[index], 0] += magnitude
    return force


def crack_opening(mesh, displacement, original_by_row, right_face_by_row, y_start, y_end):
    rows = np.arange(y_start, y_end + 1)
    aperture = np.array([
        displacement[right_face_by_row[int(row)], 0]
        - displacement[original_by_row[int(row)], 0]
        for row in rows
    ])
    y = rows * DY_M
    average = float(np.trapezoid(aperture, y) / DIKE_LENGTH_M)
    return aperture, average


def solve_pressure_case(mesh, pressure_pa, right_face_by_row, y_start, y_end):
    gradients, areas = mechanics.reference_gradients(mesh)
    identity_fp = np.repeat(np.eye(2)[None, :, :], len(mesh.triangles), axis=0)
    force = pressure_force(mesh, pressure_pa, right_face_by_row, y_start, y_end)
    displacement, deformation, cauchy, iterations, residual = mechanics.solve_equilibrium(
        mesh, gradients, areas, identity_fp, np.zeros_like(mesh.coordinates), 0.0, force,
        additional_fixed_y_nodes=(mesh.bottom_right_node,),
    )
    internal, _tangent, _deformation, _cauchy = mechanics.assemble_internal_force(
        displacement, identity_fp, mesh, gradients, areas, with_tangent=False
    )
    support = internal.reshape((-1, 2)) - force
    constrained = np.zeros(len(mesh.coordinates), dtype=bool)
    constrained[mesh.left_nodes] = True
    constrained[mesh.right_nodes] = True
    constrained[mesh.bottom_left_node] = True
    constrained[mesh.bottom_right_node] = True
    support[~constrained] = 0.0
    return {
        "gradients": gradients,
        "areas": areas,
        "force": force,
        "displacement": displacement,
        "deformation": deformation,
        "cauchy": cauchy,
        "iterations": iterations,
        "residual": residual,
        "support": support,
    }


def elastic_metrics(mesh, solution, original_by_row, right_face_by_row, y_start, y_end):
    aperture, mean_aperture = crack_opening(
        mesh, solution["displacement"], original_by_row, right_face_by_row, y_start, y_end
    )
    force = solution["force"]
    support = solution["support"]
    coordinates = mesh.coordinates + solution["displacement"]
    force_moment = float(np.sum(coordinates[:, 0] * force[:, 1] - coordinates[:, 1] * force[:, 0]))
    support_moment = float(np.sum(
        coordinates[:, 0] * support[:, 1] - coordinates[:, 1] * support[:, 0]
    ))
    force_norm = max(float(np.sum(np.linalg.norm(force, axis=1))), 1.0)
    moment_norm = max(float(np.sum(np.abs(
        coordinates[:, 0] * force[:, 1] - coordinates[:, 1] * force[:, 0]
    ))), 1.0)

    elastic_f = mechanics.element_response(
        solution["displacement"][mesh.triangles], solution["gradients"],
        np.repeat(np.eye(2)[None, :, :], len(mesh.triangles), axis=0),
    )[1]
    det_fe = np.linalg.det(elastic_f)
    log_j = np.log(det_fe)
    mu = mechanics.LAME_MU_PA
    lam = mechanics.LAME_LAMBDA_PA
    energy_density = (
        0.5 * mu * (np.einsum("eij,eij->e", elastic_f, elastic_f) - 2.0)
        - mu * log_j
        + 0.5 * lam * log_j**2
    )
    stored_energy = float(np.sum(
        mechanics.THICKNESS_M * solution["areas"] * energy_density
    ))
    pressure_work = 0.5 * float(np.sum(force * solution["displacement"]))

    principal, _direction = mechanics.principal_stress_and_direction(solution["cauchy"])
    centers = mesh.element_centers
    source_region = (
        (np.abs(centers[:, 0] - DIKE_X_M) <= 1.5 * DX_M)
        & (np.abs(centers[:, 1] - DIKE_Y_M) <= DIKE_HALF_LENGTH_M)
    )
    tip_region = (
        (np.abs(centers[:, 0] - DIKE_X_M) <= 2.0 * DX_M)
        & (np.minimum(
            np.abs(centers[:, 1] - (DIKE_Y_M - DIKE_HALF_LENGTH_M)),
            np.abs(centers[:, 1] - (DIKE_Y_M + DIKE_HALF_LENGTH_M)),
        ) <= 25.0e3)
    )

    def summary(region):
        return {
            "element_count": int(np.count_nonzero(region)),
            "mean_sigma_xx_pa": float(np.mean(solution["cauchy"][region, 0, 0])) if np.any(region) else None,
            "mean_max_principal_stress_pa": float(np.mean(principal[region])) if np.any(region) else None,
            "peak_max_principal_stress_pa": float(np.max(principal[region])) if np.any(region) else None,
        }

    source_resultant = np.sum(force, axis=0)
    support_resultant = np.sum(support, axis=0)
    return {
        "opening_m": {
            "mean": mean_aperture,
            "maximum": float(np.max(aperture)),
            "minimum": float(np.min(aperture)),
            "samples": aperture.tolist(),
        },
        "equilibrium_and_work": {
            "nonlinear_solver_relative_residual": solution["residual"],
            "newton_iterations": solution["iterations"],
            "source_force_resultant_n": source_resultant.tolist(),
            "support_reaction_resultant_n": support_resultant.tolist(),
            "relative_force_balance_residual": float(
                np.linalg.norm(source_resultant + support_resultant) / force_norm
            ),
            "source_moment_n_m": force_moment,
            "support_reaction_moment_n_m": support_moment,
            "relative_moment_balance_residual": float(
                abs(force_moment + support_moment) / moment_norm
            ),
            "elastic_pressure_work_j": pressure_work,
            "stored_elastic_energy_j": stored_energy,
            "relative_work_energy_difference": abs(pressure_work - stored_energy)
            / max(pressure_work, stored_energy, 1.0),
        },
        "stress_regions": {
            "host_rock_along_dike": summary(source_region),
            "dike_tips_and_collars": summary(tip_region),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    mesh, intact_mesh, original_by_row, right_face_by_row, y_start, y_end = build_cracked_mesh()
    calibration_pressure_pa = 1.0e6
    calibration = solve_pressure_case(
        mesh, calibration_pressure_pa, right_face_by_row, y_start, y_end
    )
    _opening, calibration_mean = crack_opening(
        mesh, calibration["displacement"], original_by_row, right_face_by_row, y_start, y_end
    )
    if calibration_mean <= 0.0:
        raise RuntimeError("positive dike-wall pressure did not open the crack")
    pressure_pa = calibration_pressure_pa * DIKE_MEAN_APERTURE_M / calibration_mean
    elastic = solve_pressure_case(mesh, pressure_pa, right_face_by_row, y_start, y_end)
    elastic_basis = elastic_metrics(
        mesh, elastic, original_by_row, right_face_by_row, y_start, y_end
    )

    magma_force = pressure_force(mesh, pressure_pa, right_face_by_row, y_start, y_end)

    cases = []
    for case_name, intrusion_on, seed_on in (
        ("tectonic_intact", False, False),
        ("tectonic_plus_dike_intact", True, False),
        ("tectonic_near_tip_seed", False, True),
        ("tectonic_plus_dike_near_tip_seed", True, True),
    ):
        case_mesh = mesh if intrusion_on else intact_mesh
        case_centers = case_mesh.element_centers
        case_source_region = (
            (np.abs(case_centers[:, 0] - DIKE_X_M) <= 1.5 * DX_M)
            & (np.abs(case_centers[:, 1] - DIKE_Y_M) <= DIKE_HALF_LENGTH_M)
        )
        case_tip_region = (
            (np.abs(case_centers[:, 0] - DIKE_X_M) <= 2.0 * DX_M)
            & (np.minimum(
                np.abs(case_centers[:, 1] - (DIKE_Y_M - DIKE_HALF_LENGTH_M)),
                np.abs(case_centers[:, 1] - (DIKE_Y_M + DIKE_HALF_LENGTH_M)),
            ) <= 25.0e3)
        )
        case_tip_seeds = (
            (np.abs(case_centers[:, 0] - DIKE_X_M) <= DX_M)
            & (
                (
                    (case_centers[:, 1] >= DIKE_Y_M + DIKE_HALF_LENGTH_M)
                    & (case_centers[:, 1] <= DIKE_Y_M + DIKE_HALF_LENGTH_M + TIP_SEED_LENGTH_M)
                )
                | (
                    (case_centers[:, 1] <= DIKE_Y_M - DIKE_HALF_LENGTH_M)
                    & (case_centers[:, 1] >= DIKE_Y_M - DIKE_HALF_LENGTH_M - TIP_SEED_LENGTH_M)
                )
            )
        )
        case_damage_seed = np.zeros(len(case_mesh.triangles))
        case_damage_seed[case_tip_seeds] = SEED_LEVEL
        case_regions = {
            "dike_source_flanks": case_source_region,
            "dike_tip_collar": case_tip_region,
            "far_field": ~(case_source_region | case_tip_region),
        }
        result = mechanics.run_case(
            nx=MESH_NX,
            ny=MESH_NY,
            dt_myr=TIMESTEP_MYR,
            duration_myr=DURATION_MYR,
            rate_per_year=RATE_PER_YEAR,
            length_m=mechanics.REFERENCE_LENGTH_M,
            seed_level=0.0,
            mesh=case_mesh,
            external_nodal_force=magma_force if intrusion_on else None,
            external_load_metadata=(
                {"type": "constant pressure on explicit dike faces", "pressure_pa": pressure_pa}
                if intrusion_on else None
            ),
            initial_damage_field=case_damage_seed if seed_on else None,
            diagnostic_regions=case_regions,
            # The same vertical pin is applied to both sides of the paired source
            # comparison so the dike load is the only changed case input.
            additional_fixed_y_nodes=(mesh.bottom_right_node,),
        )
        result.update({
            "case": case_name,
            "magma_intrusion_on": intrusion_on,
            "near_tip_seed_on": seed_on,
            "initial_seed_level": SEED_LEVEL if seed_on else 0.0,
            "initial_seed_elements": int(np.count_nonzero(case_tip_seeds)) if seed_on else 0,
            "initial_seed_geometry": "two 12.5 km collars immediately beyond dike tips",
            "dike_wall_pressure_pa": pressure_pa if intrusion_on else 0.0,
        })
        cases.append(result)
        print(f"completed {case_name}", flush=True)

    output = {
        "title": "RFC 0010 Gate 3 localized magma-intrusion source screen",
        "classification": "diagnostic-only finite pressurized crack; no production model change",
        "matched_boundary_conditions": (
            "all cases fix left/right horizontal displacement and bottom-left/bottom-right "
            "vertical displacement"
        ),
        "empirical_aperture_basis": {
            "reference": "Grandin et al. (2009), 2005 Manda Hararo-Dabbahu dike inversion",
            "dike_length_m": DIKE_LENGTH_M,
            "target_mean_aperture_m": DIKE_MEAN_APERTURE_M,
            "observed_maximum_aperture_bound_m": DIKE_MAX_APERTURE_M,
            "source_depth_range_m": [2.0e3, 8.0e3],
            "pressure_calibration": "wall pressure tapered to zero over the final 12.5 km at each crack tip, then scaled in the elastic mesh to match 4 m mean opening; peak opening is an independent check",
        },
        "mesh_and_material": {
            "mesh": f"{MESH_NX} x {MESH_NY} structured P1 triangles with duplicated internal dike-face nodes",
            "domain_length_m": DOMAIN_LENGTH_M,
            "domain_width_m": DOMAIN_WIDTH_M,
            "element_spacing_x_m": DX_M,
            "element_spacing_y_m": DY_M,
            "elastic_modulus_pa": mechanics.YOUNG_MODULUS_PA,
            "poisson_ratio": mechanics.POISSON_RATIO,
            "effective_thickness_m": mechanics.THICKNESS_M,
            "strength_pa": mechanics.YIELD_PA,
        },
        "calibrated_dike_wall_pressure_pa": pressure_pa,
        "elastic_dike_source_check": elastic_basis,
        "damage_test": {
            "tectonic_extension_rate_per_year": RATE_PER_YEAR,
            "equivalent_opening_rate_mm_per_year": RATE_PER_YEAR * DOMAIN_LENGTH_M * 1.0e3,
            "duration_myr": DURATION_MYR,
            "timestep_myr": TIMESTEP_MYR,
            "initial_near_tip_damage_strain": SEED_LEVEL,
            "matched_pairs": [
                ["tectonic_intact", "tectonic_plus_dike_intact"],
                ["tectonic_near_tip_seed", "tectonic_plus_dike_near_tip_seed"],
            ],
            "cases": cases,
        },
        "limits": [
            "The pressure is held constant during the diagnostic interval; the model does not resolve episode duration, magma transport, or post-intrusion solidification.",
            "The observed Afar dike constrains this diagnostic aperture scale only, not global dike pressure, event cadence, fracture toughness, or plume attribution.",
            "A positive result is a local-mechanics screen, not calibration or validation of geological rift propagation.",
        ],
    }
    serialized = json.dumps(output, sort_keys=True, separators=(",", ":")) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized)


if __name__ == "__main__":
    main()
