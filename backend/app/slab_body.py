"""Bounded, material-bearing 3-D slab ribbons for RFC 0010 Gate 3b (#699).

This is an experimental reduced-order body model. A ribbon carries consumed lithosphere
through a small set of down-dip stations, evolves its centerline under buoyancy, bending and
viscous mantle drag, and transmits an axial tension balance at its plate attachment.
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
    last_attachment_force_local_n: np.ndarray = field(default_factory=lambda: np.zeros(3, dtype=float))
    source_crust_volume_m3: float = 0.0
    source_mantle_lithosphere_volume_m3: float = 0.0
    source_continental_material_m3: float = 0.0
    source_craton_crust_m3: float = 0.0
    # Slab material with no surface source in this screen: rock consumed before the save
    # (not retained by production) and length fed at the trench by convergence (#699).
    antecedent_volume_m3: float = 0.0
    kinematic_inflow_volume_m3: float = 0.0

    def __post_init__(self) -> None:
        self.trench_local_xyz = geometry.normalize(np.asarray(self.trench_local_xyz, dtype=float))
        self.centerline_local_xyz_m = np.asarray(self.centerline_local_xyz_m, dtype=float)
        self.reference_centerline_local_xyz_m = np.asarray(self.reference_centerline_local_xyz_m, dtype=float)
        self.station_volume_m3 = np.asarray(self.station_volume_m3, dtype=float)
        self.density_anomaly_kg_m3 = np.asarray(self.density_anomaly_kg_m3, dtype=float)
        self.temperature_k = np.asarray(self.temperature_k, dtype=float)
        self.source_cell_keys = np.asarray(self.source_cell_keys, dtype=np.int64)
        self.last_attachment_force_local_n = np.asarray(self.last_attachment_force_local_n, dtype=float)
        if self.last_attachment_force_local_n.shape != (3,) or not np.all(np.isfinite(self.last_attachment_force_local_n)):
            raise ValueError("last attachment force must be a finite 3-vector")
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
        if min(
            self.source_crust_volume_m3,
            self.source_mantle_lithosphere_volume_m3,
            self.source_continental_material_m3,
            self.source_craton_crust_m3,
            self.antecedent_volume_m3,
            self.kinematic_inflow_volume_m3,
        ) < 0.0:
            raise ValueError("slab provenance volumes must be nonnegative")

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
            "force_active_depth_m": float(self.max_depth_m),
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
    source_crust_volume_m3: float = 0.0,
    source_mantle_lithosphere_volume_m3: float = 0.0,
    source_continental_material_m3: float = 0.0,
    source_craton_crust_m3: float = 0.0,
) -> SlabBody:
    """Create a bounded ribbon whose volume is supplied by consumed surface columns.

    The caller owns the transfer of `source_cell_keys` out of the surface mesh and must
    make their summed rock volume plus ``antecedent_volume_m3`` equal the returned
    ``consumed_volume_m3``. Dip is measured
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
        source_crust_volume_m3=source_crust_volume_m3,
        source_mantle_lithosphere_volume_m3=source_mantle_lithosphere_volume_m3,
        source_continental_material_m3=source_continental_material_m3,
        source_craton_crust_m3=source_craton_crust_m3,
    )


def attachment_force_local_n(body: SlabBody) -> tuple[np.ndarray, np.ndarray]:
    """Return the last beam-equilibrium hinge force and its equal-and-opposite reaction.

    The force is set by :func:`advance_slab_body` from the ribbon's axial tension,
    projected into the plate tangent plane because the membrane cannot accept a radial
    load. Before the first advance, no solved hinge load is available.
    """
    if not body.attached:
        zero = np.zeros(3)
        return zero, zero.copy()
    force_on_plate = np.asarray(body.last_attachment_force_local_n, dtype=float).copy()
    return force_on_plate, -force_on_plate


def apply_slab_load_to_plate_geometry(
    plate,
    bodies: list[SlabBody],
    years: float,
    *,
    effective_extensional_viscosity_pa_s: float,
    load_spread_degrees: float,
    extensional_capacity_resultant_n_per_m: np.ndarray | float,
    shell_support_thickness_m: np.ndarray | None = None,
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
        # Truncate the Gaussian at three standard deviations. Keeping its infinite
        # tails with a global weight floor would apply hinge corrections to distant
        # plate nodes, turning a local attachment into a whole-plate load (#699).
        footprint = distance <= 3.0 * angular_spread
        weights = np.exp(-0.5 * (distance / angular_spread) ** 2) * areas
        weights[~footprint] = 0.0
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
        node_weight = weights / weight_sum
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

    if shell_support_thickness_m is None:
        thickness = plate.collect("crustal_thickness_m") + plate.collect("mantle_lithosphere_thickness_m")
    else:
        thickness = np.asarray(shell_support_thickness_m, dtype=float)
    if thickness.shape != (len(points),) or not np.all(np.isfinite(thickness)) or np.any(thickness <= 0.0):
        raise ValueError("shell support thickness must be finite, positive, and match plate nodes")
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
    shell_mesh_area_m2 = 0.0
    for triangle in shell.triangles:
        _stiffness, strain_matrix, constitutive, element_area, _axis_x, _axis_y = thin_shell._triangle_element(
            points, triangle, first, second, thickness, modulus,
            POISSON_RATIO, PLANET_RADIUS_M,
        )
        element_dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        element_strain = strain_matrix @ displacement_dofs[element_dofs]
        stored_energy += 0.5 * float(element_strain @ constitutive @ element_strain) * element_area
        shell_mesh_area_m2 += element_area
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
        "shell_mesh_area_m2": shell_mesh_area_m2,
        "shell_area_coverage_fraction": shell_mesh_area_m2 / max(float(areas.sum()), 1.0),
        "shell_minimum_thickness_m": float(thickness.min(initial=np.inf)),
        "work_residual_fraction": float(work_residual),
        "max_tensile_resultant_n_per_m": shell.max_tensile_resultant_n_per_m,
        "principal_extension_xyz": principal_extension,
        "extensional_failure_ratio": shell.failure_ratio,
        "equilibrium_residual_n": shell.equilibrium_residual_n,
        "solver_stop_code": shell.solver_stop_code,
    }


def _active_path(body: SlabBody) -> tuple[float, float]:
    """Return ribbon length and depth drop above the force-active depth.

    Segments crossing that depth are clipped by interpolation, so the axial budget does
    not jump when a station crosses it (#699).
    """
    points = body.centerline_local_xyz_m
    depth = PLANET_RADIUS_M - np.linalg.norm(points, axis=1)
    clipped = np.minimum(depth, body.max_depth_m)
    span = np.diff(depth)
    fraction = np.divide(np.diff(clipped), span, out=np.ones_like(span), where=np.abs(span) > 0.0)
    lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    return float(np.sum(lengths * fraction)), float(np.sum(np.diff(clipped)))


def _extended_path(path: np.ndarray, extension_m: float) -> np.ndarray:
    """Lengthen a polyline along its tip tangent and resample its stations uniformly."""
    tip = path[-1] + extension_m * geometry.normalize(path[-1] - path[-2])
    extended = np.vstack((path, tip))
    arc = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(extended, axis=0), axis=1))))
    stations = np.linspace(0.0, arc[-1], len(path))
    sampled = np.column_stack([np.interp(stations, arc, extended[:, axis]) for axis in range(3)])
    # Chords cut the polyline's corners, so restore the full length: a shrinking path would
    # lose ribbon the inflow ledger still counts (#699).
    directions = geometry.normalize(np.diff(sampled, axis=0))
    spacing = arc[-1] / (len(path) - 1)
    return path[0] + np.vstack((np.zeros(3), np.cumsum(directions * spacing, axis=0)))


def advance_slab_body(
    body: SlabBody,
    years: float,
    *,
    mantle_viscosity_pa_s: float,
    shear_zone_thickness_m: float,
    convergence_speed_m_per_s: float = 0.0,
    axial_drag_factor: float = 1.0,
) -> dict[str, float]:
    """Advance a ribbon and solve its hinge load from an axial force balance.

    The ribbon is inextensible down-dip. Its active weight along the path, less mantle
    drag on its sliding at the convergence speed, is the tension it can transmit. A
    negative balance is not a push on the plate: the slab cannot be sunk by its own weight
    at that speed, so the remainder is booked as driver force and the hinge carries none.
    Weight normal to the path and bending evolve the centerline through an implicit
    overdamped beam solve. The earlier chain of free nodes had no axial member: its weight
    went to each node's own drag, and its hinge load was stretch against the reference
    line, a factor of 1e-7 of the weight (#699).

    Convergence lengthens the tip rather than stretching the first segment. The inflow is
    booked in `kinematic_inflow_volume_m3` because this screen does not move the plate into
    the trench. `axial_drag_factor` scales only the axial drag, to bracket gross slab pull.
    """
    if years <= 0.0 or not np.isfinite(years):
        raise ValueError("years must be finite and positive")
    if mantle_viscosity_pa_s <= 0.0 or shear_zone_thickness_m <= 0.0 or convergence_speed_m_per_s < 0.0:
        raise ValueError("mantle viscosity and shear-zone thickness must be positive")
    if axial_drag_factor < 0.0:
        raise ValueError("axial drag factor must be nonnegative")
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
    segment_directions = geometry.normalize(np.diff(points, axis=0))
    tangents = np.empty_like(points)
    tangents[0] = segment_directions[0]
    tangents[-1] = segment_directions[-1]
    tangents[1:-1] = geometry.normalize(segment_directions[:-1] + segment_directions[1:])
    normal_buoyancy = buoyancy - np.einsum("ij,ij->i", buoyancy, tangents)[:, None] * tangents
    system = stiffness[free, free] + np.eye(count - 1) * drag_over_dt
    rhs = (
        drag_over_dt * points[free]
        + normal_buoyancy[free]
        + (stiffness @ body.reference_centerline_local_xyz_m)[free]
        - stiffness[free, 0, None] * body.trench_local_xyz[None, :] * PLANET_RADIUS_M
    )
    updated = points.copy()
    updated[free] = np.linalg.solve(system, rhs)
    updated[0] = body.trench_local_xyz * PLANET_RADIUS_M
    velocity = (updated - points) / dt_s
    displacement = updated - points
    bending_force = -(stiffness @ (updated - body.reference_centerline_local_xyz_m))
    drag_force = -drag * velocity
    beam_work_j = float(np.sum(normal_buoyancy * displacement))
    drag_dissipation_j = float(drag * np.sum(velocity**2) * dt_s)
    before_delta = points - body.reference_centerline_local_xyz_m
    after_delta = updated - body.reference_centerline_local_xyz_m
    bend_energy_before_j = 0.5 * float(np.sum(before_delta * (stiffness @ before_delta)))
    bend_energy_after_j = 0.5 * float(np.sum(after_delta * (stiffness @ after_delta)))
    bending_work_j = float(np.sum((stiffness @ after_delta) * displacement))
    beam_scale = max(abs(beam_work_j), drag_dissipation_j, abs(bending_work_j), 1.0)
    beam_residual = abs(beam_work_j - drag_dissipation_j - bending_work_j) / beam_scale
    beam_force_residual = (normal_buoyancy + bending_force + drag_force)[1:]

    path_length = float(np.sum(np.linalg.norm(np.diff(points, axis=0), axis=1)))
    line_weight_n_per_m = float(
        np.dot(body.density_anomaly_kg_m3, body.station_volume_m3) * GRAVITY_M_S2 / path_length
    )
    active_length_m, active_depth_m = _active_path(body)
    axial_buoyancy_n = line_weight_n_per_m * active_depth_m
    axial_drag_n = (
        axial_drag_factor * mantle_viscosity_pa_s * body.width_m * active_length_m
        * convergence_speed_m_per_s / shear_zone_thickness_m
    )
    axial_balance_n = axial_buoyancy_n - axial_drag_n
    tension_n = max(axial_balance_n, 0.0) if body.attached else 0.0
    driver_force_n = max(-axial_balance_n, 0.0) if body.attached else 0.0
    distance_m = convergence_speed_m_per_s * dt_s
    axial_buoyancy_work_j = axial_buoyancy_n * distance_m
    axial_drag_work_j = axial_drag_n * distance_m
    hinge_work_j = tension_n * distance_m
    driver_work_j = driver_force_n * distance_m
    axial_scale = max(axial_buoyancy_work_j, axial_drag_work_j, 1.0)
    if body.attached:
        axial_residual = abs(
            axial_buoyancy_work_j + driver_work_j - axial_drag_work_j - hinge_work_j
        ) / axial_scale
    else:
        axial_residual = 0.0

    # The membrane takes only the tension's tangent part. The root shear goes into plate
    # flexure and the mantle reaction at the bend: adding it clamps the hinge, and a clamp
    # holds a hanging slab vertically, cancelling the in-plane pull (#699).
    radial_hinge = body.trench_local_xyz
    hinge_load = tension_n * tangents[0]
    hinge_radial_n = float(np.dot(hinge_load, radial_hinge))
    tangent_force = hinge_load - hinge_radial_n * radial_hinge
    body.last_attachment_force_local_n = tangent_force if body.attached else np.zeros(3)
    # The beam solve is linear, so a sinking path would lengthen without limit. Restore
    # each segment's length from the hinge outward to keep the ribbon inextensible (#699).
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    solved_directions = geometry.normalize(np.diff(updated, axis=0))
    updated = updated[0] + np.vstack((
        np.zeros(3), np.cumsum(solved_directions * segment_lengths[:, None], axis=0)
    ))
    if distance_m > 0.0:
        volume_before_m3 = body.consumed_volume_m3
        inflow_m3 = distance_m * volume_before_m3 / path_length
        extended = _extended_path(updated, distance_m)
        body.reference_centerline_local_xyz_m = _extended_path(
            body.reference_centerline_local_xyz_m, distance_m
        )
        share = np.linalg.norm(np.diff(extended, axis=0), axis=1)
        station_length = np.empty(count)
        station_length[0] = 0.5 * share[0]
        station_length[-1] = 0.5 * share[-1]
        station_length[1:-1] = 0.5 * (share[:-1] + share[1:])
        # Normalize to the booked total: beam sag also lengthens the path, and that is
        # not new material (#699).
        body.station_volume_m3 = station_length * (volume_before_m3 + inflow_m3) / station_length.sum()
        body.kinematic_inflow_volume_m3 += inflow_m3
        updated = extended
    body.centerline_local_xyz_m = updated
    body.deep_reservoir_volume_m3 = float(body.station_volume_m3[~body.active_station_mask].sum())
    return {
        "beam_normal_buoyancy_work_j": beam_work_j,
        "mantle_drag_dissipation_j": drag_dissipation_j,
        "bending_energy_before_j": bend_energy_before_j,
        "bending_energy_after_j": bend_energy_after_j,
        "bending_work_j": bending_work_j,
        "beam_work_residual_fraction": beam_residual,
        "axial_buoyancy_force_n": axial_buoyancy_n,
        "axial_mantle_drag_force_n": axial_drag_n,
        "hinge_tension_n": tension_n,
        "plate_driver_force_n": driver_force_n,
        "axial_buoyancy_work_j": axial_buoyancy_work_j,
        "axial_mantle_drag_work_j": axial_drag_work_j,
        "hinge_work_j": hinge_work_j,
        "plate_driver_work_j": driver_work_j,
        "axial_work_residual_fraction": axial_residual,
        "work_residual_fraction": max(beam_residual, axial_residual),
        "max_station_speed_m_per_s": float(np.linalg.norm(velocity, axis=1).max(initial=0.0)),
        "max_station_displacement_m": float(np.linalg.norm(displacement, axis=1).max(initial=0.0)),
        "advected_distance_m": float(distance_m),
        "active_length_m": active_length_m,
        "active_depth_m": active_depth_m,
        "attachment_force_local_n": body.last_attachment_force_local_n.copy(),
        "attachment_reaction_local_n": -body.last_attachment_force_local_n,
        "hinge_radial_support_n": hinge_radial_n,
        "hinge_bending_shear_n": float(np.linalg.norm(bending_force[0])),
        "beam_force_residual_n": float(np.linalg.norm(beam_force_residual)),
        "station_radius_min_m": float(radius.min(initial=PLANET_RADIUS_M)),
    }
