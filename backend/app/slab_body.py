"""Bounded, material-bearing 3-D slab ribbons for RFC 0010 Gate 3b (#699).

This is an experimental reduced-order body model. A ribbon carries consumed lithosphere
through a small set of down-dip stations, evolves its centerline under buoyancy, bending and
viscous mantle drag, and transmits the integrated buoyant tension at its plate attachment.
It does not solve mantle flow or represent a full 3-D continuum.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import geometry
from .lithosphere import PLANET_RADIUS_M

GRAVITY_M_S2 = 9.81
SECONDS_PER_YEAR = 365.25 * 24.0 * 60.0 * 60.0
POISSON_RATIO = 0.25


@dataclass
class SlabBody:
    """One attached slab ribbon in a subducting plate's co-rotating frame."""

    body_id: int
    plate_id: int
    trench_local_xyz: np.ndarray
    centerline_local_xyz_m: np.ndarray
    reference_centerline_local_xyz_m: np.ndarray
    station_volume_m3: np.ndarray
    density_anomaly_kg_m3: np.ndarray
    temperature_k: np.ndarray
    width_m: float
    thickness_m: float
    max_depth_m: float
    bending_rigidity_nm2: float
    source_cell_keys: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    attached: bool = True
    deep_reservoir_volume_m3: float = 0.0

    def __post_init__(self) -> None:
        self.trench_local_xyz = geometry.normalize(np.asarray(self.trench_local_xyz, dtype=float))
        self.centerline_local_xyz_m = np.asarray(self.centerline_local_xyz_m, dtype=float)
        self.reference_centerline_local_xyz_m = np.asarray(self.reference_centerline_local_xyz_m, dtype=float)
        self.station_volume_m3 = np.asarray(self.station_volume_m3, dtype=float)
        self.density_anomaly_kg_m3 = np.asarray(self.density_anomaly_kg_m3, dtype=float)
        self.temperature_k = np.asarray(self.temperature_k, dtype=float)
        self.source_cell_keys = np.asarray(self.source_cell_keys, dtype=np.int64)
        count = len(self.station_volume_m3)
        if count < 3 or self.centerline_local_xyz_m.shape != (count, 3):
            raise ValueError("slab body requires at least three 3-D stations and one volume per station")
        if self.reference_centerline_local_xyz_m.shape != (count, 3):
            raise ValueError("reference slab centerline must match station geometry")
        if self.density_anomaly_kg_m3.shape != (count,) or self.temperature_k.shape != (count,):
            raise ValueError("slab density and temperature arrays must match the station count")
        if not np.all(np.isfinite(self.centerline_local_xyz_m)) or np.any(self.station_volume_m3 <= 0.0):
            raise ValueError("slab geometry must be finite and station volumes positive")
        if np.any(self.density_anomaly_kg_m3 < 0.0) or np.any(self.temperature_k <= 0.0):
            raise ValueError("slab density anomaly must be nonnegative and temperatures positive")
        if min(self.width_m, self.thickness_m, self.max_depth_m, self.bending_rigidity_nm2) <= 0.0:
            raise ValueError("slab dimensions, depth limit, and bending rigidity must be positive")

    @property
    def consumed_volume_m3(self) -> float:
        """Material currently represented by the ribbon, including deep-reservoir stations."""
        return float(self.station_volume_m3.sum())

    @property
    def active_station_mask(self) -> np.ndarray:
        """Stations still inside the modeled mantle; deeper material stays accounted."""
        depth = PLANET_RADIUS_M - np.linalg.norm(self.centerline_local_xyz_m, axis=1)
        return depth <= self.max_depth_m + 1e-6

    def metrics(self) -> dict[str, float | int]:
        """Return physical geometry and material measures for gate reports."""
        points = self.centerline_local_xyz_m
        segments = np.diff(points, axis=0)
        lengths = np.linalg.norm(segments, axis=1)
        directions = geometry.normalize(segments)
        turns = np.arccos(np.clip(np.einsum("ij,ij->i", directions[:-1], directions[1:]), -1.0, 1.0))
        curvature = turns / np.maximum(0.5 * (lengths[:-1] + lengths[1:]), np.finfo(float).tiny)
        radial = self.trench_local_xyz
        first = directions[0]
        radial_component = float(np.dot(first, radial))
        dip = np.arctan2(-radial_component, np.linalg.norm(first - radial_component * radial))
        depth = PLANET_RADIUS_M - np.linalg.norm(points, axis=1)
        return {
            "station_count": int(len(points)),
            "length_m": float(lengths.sum()),
            "max_depth_m": float(depth.max(initial=0.0)),
            "dip_degrees": float(np.degrees(dip)),
            "mean_curvature_m_inv": float(curvature.mean()) if len(curvature) else 0.0,
            "max_curvature_m_inv": float(curvature.max(initial=0.0)),
            "width_m": float(self.width_m),
            "thickness_m": float(self.thickness_m),
            "consumed_volume_m3": self.consumed_volume_m3,
            "deep_reservoir_volume_m3": float(self.station_volume_m3[~self.active_station_mask].sum()),
        }


def create_slab_body(
    *,
    body_id: int,
    plate_id: int,
    plate_frame: np.ndarray,
    trench_world_xyz: np.ndarray,
    down_dip_world_xyz: np.ndarray,
    source_cell_keys: np.ndarray,
    width_m: float,
    thickness_m: float,
    max_depth_m: float,
    station_count: int,
    dip_degrees: float,
    density_anomaly_kg_m3: float | np.ndarray,
    surface_temperature_k: float,
    geothermal_gradient_k_per_m: float,
    bending_modulus_pa: float,
) -> SlabBody:
    """Create a bounded ribbon whose volume is supplied by consumed surface columns.

    The caller owns the transfer of `source_cell_keys` out of the surface mesh and must
    make their summed rock volume equal the returned ``consumed_volume_m3``. Dip is measured
    below horizontal; temperatures and density anomaly are explicit scenario inputs.
    """
    if station_count < 3 or not 0.0 < dip_degrees < 90.0:
        raise ValueError("slab ribbon needs at least three stations and a dip between 0 and 90 degrees")
    if not np.all(np.isfinite([width_m, thickness_m, max_depth_m, surface_temperature_k, geothermal_gradient_k_per_m, bending_modulus_pa])):
        raise ValueError("slab dimensions and material properties must be finite")
    frame = np.asarray(plate_frame, dtype=float)
    if frame.shape != (3, 3) or not np.all(np.isfinite(frame)):
        raise ValueError("plate_frame must be a finite 3x3 rotation matrix")
    trench_local = geometry.normalize(geometry.to_local(frame, np.asarray(trench_world_xyz, dtype=float)))
    inward_local = geometry.normalize(geometry.to_local(frame, np.asarray(down_dip_world_xyz, dtype=float)))
    inward_local -= np.dot(inward_local, trench_local) * trench_local
    inward_local = geometry.normalize(inward_local)
    depths = np.linspace(0.0, max_depth_m, station_count)
    horizontal = depths / np.tan(np.radians(dip_degrees))
    # Keep each station on its requested radial depth while moving the trenchward
    # surface projection inward by the dip-controlled distance. A straight Cartesian
    # offset makes the chord's radial norm shallower than its specified depth.
    radial_midpoints = PLANET_RADIUS_M - 0.5 * (depths[:-1] + depths[1:])
    angular_steps = np.diff(horizontal) / np.maximum(radial_midpoints, 1.0)
    angular_offset = np.concatenate(([0.0], np.cumsum(angular_steps)))
    radial_direction = (
        np.cos(angular_offset)[:, None] * trench_local
        + np.sin(angular_offset)[:, None] * inward_local
    )
    centerline = (PLANET_RADIUS_M - depths[:, None]) * radial_direction
    segment_lengths = np.linalg.norm(np.diff(centerline, axis=0), axis=1)
    station_lengths = np.empty(station_count)
    station_lengths[0] = 0.5 * segment_lengths[0]
    station_lengths[-1] = 0.5 * segment_lengths[-1]
    station_lengths[1:-1] = 0.5 * (segment_lengths[:-1] + segment_lengths[1:])
    volumes = station_lengths * width_m * thickness_m
    density = np.broadcast_to(np.asarray(density_anomaly_kg_m3, dtype=float), (station_count,)).copy()
    temperatures = surface_temperature_k + geothermal_gradient_k_per_m * depths
    rigidity = bending_modulus_pa * width_m * thickness_m**3 / 12.0
    return SlabBody(
        body_id=body_id,
        plate_id=plate_id,
        trench_local_xyz=trench_local,
        centerline_local_xyz_m=centerline,
        reference_centerline_local_xyz_m=centerline.copy(),
        station_volume_m3=volumes,
        density_anomaly_kg_m3=density,
        temperature_k=temperatures,
        width_m=width_m,
        thickness_m=thickness_m,
        max_depth_m=max_depth_m,
        bending_rigidity_nm2=rigidity,
        source_cell_keys=source_cell_keys,
    )


def attachment_force_local_n(body: SlabBody) -> tuple[np.ndarray, np.ndarray]:
    """Return equal-and-opposite plate and slab hinge forces in the plate-local frame."""
    if not body.attached:
        zero = np.zeros(3)
        return zero, zero.copy()
    active = body.active_station_mask
    tension_n = float(np.dot(body.station_volume_m3[active], body.density_anomaly_kg_m3[active])) * GRAVITY_M_S2
    first = body.centerline_local_xyz_m[1] - body.centerline_local_xyz_m[0]
    inward = first - np.dot(first, body.trench_local_xyz) * body.trench_local_xyz
    inward = geometry.normalize(inward)
    force_on_plate = tension_n * inward
    return force_on_plate, -force_on_plate


def apply_slab_load_to_plate_geometry(
    plate,
    bodies: list[SlabBody],
    years: float,
    *,
    effective_extensional_viscosity_pa_s: float,
    load_spread_degrees: float,
    extensional_capacity_resultant_n_per_m: np.ndarray | float,
    additional_force_local_n: np.ndarray | None = None,
) -> dict:
    """Transfer attached-ribbon forces through the shell and apply their strain to a plate.

    This controlled Gate 3b coupling uses a linear membrane response with an effective
    modulus ``viscosity / timestep``. It retains only the fitted homogeneous strain when
    updating sparse-quad geometry; the full local stress/failure field remains diagnostic
    output for the caller. Extensional capacity is supplied explicitly; the caller owns slab
    material transfer and local crack growth.
    """
    from . import thin_shell

    if years <= 0.0 or not np.isfinite(years):
        raise ValueError("years must be finite and positive")
    if effective_extensional_viscosity_pa_s <= 0.0 or load_spread_degrees <= 0.0:
        raise ValueError("extensional viscosity and load spread must be positive")
    if any(body.plate_id != plate.plate_id for body in bodies):
        raise ValueError("all slab bodies must attach to the supplied plate")
    if not bodies and additional_force_local_n is None:
        raise ValueError("an attached slab body or additional force field is required")

    nodes = plate.surface_nodes()
    points = np.asarray(nodes.local_xyz, dtype=float)
    areas = np.asarray(nodes.area_m2, dtype=float)
    forces = np.zeros_like(points)
    normal_reaction = np.zeros(3)
    angular_spread = np.radians(load_spread_degrees)
    total_force = np.zeros(3)
    source_torque = np.zeros(3)
    source_force_norm_sum = 0.0
    source_torque_norm_sum = 0.0
    for body in bodies:
        force_local, _reaction = attachment_force_local_n(body)
        total_force += force_local
        body_torque = np.cross(body.trench_local_xyz * PLANET_RADIUS_M, force_local)
        source_torque += body_torque
        source_force_norm_sum += float(np.linalg.norm(force_local))
        source_torque_norm_sum += float(np.linalg.norm(body_torque))
        cosine = np.clip(points @ body.trench_local_xyz, -1.0, 1.0)
        distance = np.arccos(cosine)
        weights = np.exp(-0.5 * (distance / angular_spread) ** 2) * areas
        weight_sum = float(weights.sum())
        if weight_sum <= 0.0:
            raise ValueError("slab attachment has no plate nodes in its load footprint")
        tangent_force = force_local[None, :] - (
            np.einsum("ij,j->i", points, force_local)[:, None] * points
        )
        weighted = (weights / weight_sum)[:, None]
        body_force = tangent_force * weighted
        body_normal_reaction = np.sum(
            (np.einsum("ij,j->i", points, force_local) * weighted[:, 0])[:, None] * points,
            axis=0,
        )
        # Spread the hinge load over the registered footprint while preserving the
        # attachment's exact resultant and moment. The unconstrained Gaussian transfers
        # force but shifts its line of action on a curved plate (#699).
        first, second = thin_shell._tangent_basis(points)
        basis = np.zeros((6, 2 * len(points)))
        basis[:3, 0::2] = first.T
        basis[:3, 1::2] = second.T
        basis[3:, 0::2] = np.cross(points, first).T
        basis[3:, 1::2] = np.cross(points, second).T
        desired_moment = np.cross(body.trench_local_xyz * PLANET_RADIUS_M, force_local)
        current_moment = np.cross(points * PLANET_RADIUS_M, body_force).sum(axis=0)
        target = np.concatenate((force_local - body_force.sum(axis=0) - body_normal_reaction,
                                 (desired_moment - current_moment) / PLANET_RADIUS_M))
        node_weight = np.maximum(weights / weight_sum, float(np.max(weights / weight_sum)) * 1e-4)
        dof_weight = np.repeat(node_weight, 2)
        normal_matrix = (basis * dof_weight[None, :]) @ basis.T
        lagrange = np.linalg.lstsq(normal_matrix, target, rcond=None)[0]
        correction = dof_weight * (basis.T @ lagrange)
        body_force += first * correction[0::2, None] + second * correction[1::2, None]
        forces += body_force
        normal_reaction += body_normal_reaction

    if additional_force_local_n is not None:
        extra = np.asarray(additional_force_local_n, dtype=float)
        if extra.shape != points.shape or not np.all(np.isfinite(extra)):
            raise ValueError("additional_force_local_n must be a finite (N, 3) force field")
        radial_force = np.einsum("ij,ij->i", points, extra)
        forces += extra - radial_force[:, None] * points
        total_force += extra.sum(axis=0)
        extra_torque = np.cross(points * PLANET_RADIUS_M, extra)
        source_torque += extra_torque.sum(axis=0)
        source_force_norm_sum += float(np.linalg.norm(extra, axis=1).sum())
        source_torque_norm_sum += float(np.linalg.norm(extra_torque, axis=1).sum())
        normal_reaction += np.sum(radial_force[:, None] * points, axis=0)

    thickness = np.maximum(
        plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m"),
        1.0,
    )
    timestep_s = years * SECONDS_PER_YEAR
    young_modulus_pa = effective_extensional_viscosity_pa_s / timestep_s
    capacity = np.broadcast_to(
        np.asarray(extensional_capacity_resultant_n_per_m, dtype=float), (len(points),)
    ).copy()
    shell = thin_shell.solve_thin_shell(
        points_xyz=points,
        areas_m2=areas,
        applied_force_n=forces,
        thickness_m=thickness,
        young_modulus_pa=young_modulus_pa,
        extensional_capacity_resultant_n_per_m=capacity,
    )
    strain = thin_shell.fit_homogeneous_strain_tensor(points, shell.displacement_m)
    max_log_strain = plate.apply_surface_strain(strain)
    load_work_j = float(np.sum(forces * shell.displacement_m))
    applied_resultant = forces.sum(axis=0)
    applied_torque = np.cross(points * PLANET_RADIUS_M, forces).sum(axis=0)

    # Independently integrate the linear membrane strain energy for the registered work
    # closure. The reduced solve reports the full nodal displacement and shell triangulation.
    first, second = thin_shell._tangent_basis(points)
    displacement_dofs = np.column_stack((
        np.sum(shell.displacement_m * first, axis=1),
        np.sum(shell.displacement_m * second, axis=1),
    )).ravel()
    modulus = np.broadcast_to(np.asarray(young_modulus_pa, dtype=float), (len(points),))
    stored_energy = 0.0
    for triangle in shell.triangles:
        _stiffness, strain_matrix, constitutive, element_area, _axis_x, _axis_y = thin_shell._triangle_element(
            points, triangle, first, second, thickness, modulus,
            POISSON_RATIO, PLANET_RADIUS_M,
        )
        element_dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        element_strain = strain_matrix @ displacement_dofs[element_dofs]
        stored_energy += 0.5 * float(element_strain @ constitutive @ element_strain) * element_area
    equilibrated_work = 0.5 * float(np.sum(shell.equilibrated_force_n * shell.displacement_m))
    work_residual = abs(equilibrated_work - stored_energy) / max(abs(equilibrated_work), abs(stored_energy), 1.0)

    stress = shell.stress_resultant_n_per_m
    stress_matrices = np.empty((len(stress), 2, 2), dtype=float)
    stress_matrices[:, 0, 0] = stress[:, 0]
    stress_matrices[:, 1, 1] = stress[:, 1]
    stress_matrices[:, 0, 1] = stress_matrices[:, 1, 0] = stress[:, 2]
    _principal_values, principal_vectors = np.linalg.eigh(stress_matrices)
    principal_extension = first * principal_vectors[:, 0, 1, None] + second * principal_vectors[:, 1, 1, None]
    return {
        "total_attachment_force_n": total_force,
        "source_force_norm_sum_n": float(source_force_norm_sum),
        "applied_force_n": forces,
        "applied_resultant_force_n": applied_resultant,
        "normal_reaction_force_n": normal_reaction,
        "force_transfer_error_n": applied_resultant + normal_reaction - total_force,
        "source_torque_nm": source_torque,
        "source_torque_norm_sum_nm": float(source_torque_norm_sum),
        "applied_torque_nm": applied_torque,
        "torque_transfer_error_nm": applied_torque - source_torque,
        "shell_displacement_m": shell.displacement_m,
        "stress_resultant_n_per_m": shell.stress_resultant_n_per_m,
        "failure_ratio": shell.failure_ratio,
        "strain_increment": strain,
        "max_log_strain": max_log_strain,
        "load_work_j": load_work_j,
        "stored_elastic_energy_j": float(stored_energy),
        "work_residual_fraction": float(work_residual),
        "max_tensile_resultant_n_per_m": shell.max_tensile_resultant_n_per_m,
        "principal_extension_xyz": principal_extension,
        "extensional_failure_ratio": shell.failure_ratio,
        "equilibrium_residual_n": shell.equilibrium_residual_n,
        "solver_stop_code": shell.solver_stop_code,
    }


def advance_slab_body(
    body: SlabBody,
    years: float,
    *,
    mantle_viscosity_pa_s: float,
    shear_zone_thickness_m: float,
    convergence_speed_m_per_s: float = 0.0,
) -> dict[str, float]:
    """Advance body bending under buoyancy with an implicit overdamped beam solve.

    The first station stays attached to the plate. Bending resists motion away from the
    initial centerline; mantle viscosity resists each station's velocity. The update is
    implicit in both resistances, so large geological timesteps do not destabilize the
    station chain.
    """
    if years <= 0.0 or not np.isfinite(years):
        raise ValueError("years must be finite and positive")
    if mantle_viscosity_pa_s <= 0.0 or shear_zone_thickness_m <= 0.0 or convergence_speed_m_per_s < 0.0:
        raise ValueError("mantle viscosity and shear-zone thickness must be positive")
    points = body.centerline_local_xyz_m
    count = len(points)
    dt_s = years * SECONDS_PER_YEAR
    spacing = float(np.mean(np.linalg.norm(np.diff(points, axis=0), axis=1)))
    second_difference = np.zeros((count - 2, count))
    rows = np.arange(count - 2)
    second_difference[rows, rows] = 1.0 / spacing**2
    second_difference[rows, rows + 1] = -2.0 / spacing**2
    second_difference[rows, rows + 2] = 1.0 / spacing**2
    stiffness = body.bending_rigidity_nm2 * spacing * second_difference.T @ second_difference
    drag = mantle_viscosity_pa_s * body.width_m * spacing / shear_zone_thickness_m
    drag_over_dt = drag / dt_s
    radius = np.linalg.norm(points, axis=1)
    radial = geometry.normalize(points)
    buoyancy = -(
        body.density_anomaly_kg_m3 * GRAVITY_M_S2 * body.station_volume_m3
    )[:, None] * radial
    active = body.active_station_mask
    buoyancy[~active] = 0.0
    free = slice(1, count)
    system = stiffness[free, free] + np.eye(count - 1) * drag_over_dt
    rhs = (
        drag_over_dt * points[free]
        + buoyancy[free]
        + (stiffness @ body.reference_centerline_local_xyz_m)[free]
        - stiffness[free, 0, None] * body.trench_local_xyz[None, :] * PLANET_RADIUS_M
    )
    updated = points.copy()
    updated[free] = np.linalg.solve(system, rhs)
    updated[0] = body.trench_local_xyz * PLANET_RADIUS_M
    beam_updated = updated.copy()

    # Subduction feeds the ribbon down-dip at the measured plate-convergence speed. The
    # attachment stays fixed while interior stations advect; material passing max_depth_m
    # remains in the body's deep-reservoir account instead of disappearing (#699).
    advected_distance = convergence_speed_m_per_s * dt_s
    if advected_distance > 0.0:
        radial = geometry.normalize(updated)
        first = updated[1] - updated[0]
        radial_component = float(np.dot(first, body.trench_local_xyz))
        horizontal = geometry.normalize(first - radial_component * body.trench_local_xyz)
        dip = float(np.arctan2(-radial_component, np.linalg.norm(first - radial_component * body.trench_local_xyz)))
        local_horizontal = horizontal[None, :] - (
            np.einsum("ij,j->i", radial, horizontal)[:, None] * radial
        )
        local_horizontal = geometry.normalize(local_horizontal)
        angular_step = np.cos(dip) * advected_distance / np.maximum(
            np.linalg.norm(updated, axis=1), 1.0
        )
        advected_radial = (
            np.cos(angular_step)[:, None] * radial
            + np.sin(angular_step)[:, None] * local_horizontal
        )
        advected_depth = PLANET_RADIUS_M - np.linalg.norm(updated, axis=1)
        advected_depth += np.sin(dip) * advected_distance
        updated = np.maximum(PLANET_RADIUS_M - advected_depth, 1.0)[:, None] * advected_radial
        updated[0] = body.trench_local_xyz * PLANET_RADIUS_M
    # The force balance below checks the implicit buoyancy/bending/drag relaxation. The
    # measured convergence advection is an imposed kinematic transport and is reported
    # separately, not folded into that relaxation residual.
    velocity = (beam_updated - points) / dt_s
    bending_force = -(stiffness @ (beam_updated - body.reference_centerline_local_xyz_m))
    drag_force = -drag * velocity
    external_work_j = float(np.sum(buoyancy * (beam_updated - points)))
    drag_dissipation_j = float(drag * np.sum(velocity**2) * dt_s)
    before_delta = points - body.reference_centerline_local_xyz_m
    after_delta = beam_updated - body.reference_centerline_local_xyz_m
    bend_energy_before_j = 0.5 * float(np.sum(before_delta * (stiffness @ before_delta)))
    bend_energy_after_j = 0.5 * float(np.sum(after_delta * (stiffness @ after_delta)))
    bending_work_j = float(np.sum((stiffness @ after_delta) * (beam_updated - points)))
    work_scale = max(abs(external_work_j), drag_dissipation_j, abs(bend_energy_after_j - bend_energy_before_j), 1.0)
    residual = abs(external_work_j - drag_dissipation_j - bending_work_j) / work_scale
    body.centerline_local_xyz_m = updated
    body.deep_reservoir_volume_m3 = float(body.station_volume_m3[~body.active_station_mask].sum())
    return {
        "external_buoyancy_work_j": external_work_j,
        "mantle_drag_dissipation_j": drag_dissipation_j,
        "bending_energy_before_j": bend_energy_before_j,
        "bending_energy_after_j": bend_energy_after_j,
        "bending_work_j": bending_work_j,
        "work_residual_fraction": residual,
        "max_station_speed_m_per_s": float(np.linalg.norm(velocity, axis=1).max(initial=0.0)),
        "max_station_displacement_m": float(np.linalg.norm(updated - points, axis=1).max(initial=0.0)),
        "advected_distance_m": float(advected_distance),
        "work_balance_excludes_kinematic_advection": True,
        "bending_force_resultant_n": float(np.linalg.norm(bending_force.sum(axis=0))),
        "mantle_drag_force_resultant_n": float(np.linalg.norm(drag_force.sum(axis=0))),
        "buoyancy_force_resultant_n": float(np.linalg.norm(buoyancy.sum(axis=0))),
        "station_radius_min_m": float(radius.min(initial=PLANET_RADIUS_M)),
    }
