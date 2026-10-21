"""Reduced thin-shell membrane diagnostic for per-node tectonic loads.

This is an offline diagnostic, not a production plate-motion or deformation solver. It
projects applied loads onto the shell tangent plane, removes the part that can accelerate a
rigid plate rotation, solves linear membrane equilibrium on a triangular mesh, and reports
incremental depth-integrated stress resultants. A caller supplies the column's additional
extensional capacity in N/m; this module does not invent a geotherm or lithosphere strength law.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsmr
from scipy.spatial import ConvexHull

from . import geometry, lithosphere


@dataclass(frozen=True)
class ThinShellResult:
    """Nodal incremental extension and capacity fields from membrane equilibrium."""

    displacement_m: np.ndarray  # (N, 3), tangent displacement
    stress_resultant_n_per_m: np.ndarray  # (N, 3): xx, yy, xy in each node's tangent basis
    max_tensile_resultant_n_per_m: np.ndarray  # (N,), compatibility alias for max extension increment
    tensile_strength_resultant_n_per_m: np.ndarray  # (N,), compatibility alias for extensional capacity
    failure_ratio: np.ndarray  # (N,), compatibility alias for extensional failure ratio
    equilibrated_force_n: np.ndarray  # (N, 3), applied tangent force after rigid-mode removal
    rigid_rotation_force_n: np.ndarray  # (N, 3), removed rigid-rotation load component
    fixed_node_mask: np.ndarray  # (N,), nodes constrained to zero tangent displacement
    triangles: np.ndarray  # (T, 3), shell elements
    solver_stop_code: int
    equilibrium_residual_n: float

    @property
    def max_incremental_extensional_resultant_n_per_m(self) -> np.ndarray:
        """Largest positive principal stress increment, measured from zero load."""
        return self.max_tensile_resultant_n_per_m

    @property
    def extensional_capacity_resultant_n_per_m(self) -> np.ndarray:
        """Integrated additional strength to normal-fault yield for the chosen K0."""
        return self.tensile_strength_resultant_n_per_m

    @property
    def extensional_failure_ratio(self) -> np.ndarray:
        """Incremental extensional resultant divided by additional column strength."""
        return self.failure_ratio


def _tangent_basis(points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return stable orthonormal tangent axes at each unit-sphere point."""
    reference = np.zeros_like(points)
    use_z = np.abs(points[:, 2]) < 0.9
    reference[use_z, 2] = 1.0
    reference[~use_z, 0] = 1.0
    first = geometry.normalize(np.cross(points, reference))
    second = np.cross(points, first)
    return first, second


def _spherical_surface_triangles(points: np.ndarray) -> np.ndarray:
    """Keep outward-facing spherical facets and discard the convex-hull closure below a cap."""
    hull = ConvexHull(points)
    centroids = geometry.normalize(points[hull.simplices].mean(axis=1))
    radial_alignment = np.sum(hull.equations[:, :3] * centroids, axis=1)
    triangles = hull.simplices[radial_alignment > 0.5]
    if len(triangles) == 0:
        raise ValueError("node set has no outward-facing spherical hull facets")
    return np.asarray(triangles, dtype=np.int64)


def _mesh_component_labels(triangles: np.ndarray) -> np.ndarray:
    """Label triangle components joined through complete shared edges."""
    edges = np.sort(
        np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]])),
        axis=1,
    )
    edge_order = np.lexsort((edges[:, 1], edges[:, 0]))
    sorted_edges = edges[edge_order]
    edge_triangles = np.tile(np.arange(len(triangles)), 3)[edge_order]
    same_edge = np.all(sorted_edges[1:] == sorted_edges[:-1], axis=1)
    edge_groups = np.split(edge_triangles, np.flatnonzero(~same_edge) + 1)
    triangle_neighbors = [np.unique(group) for group in edge_groups if len(np.unique(group)) > 1]
    if triangle_neighbors:
        pairs = np.column_stack((
            np.concatenate([group[:1].repeat(len(group) - 1) for group in triangle_neighbors]),
            np.concatenate([group[1:] for group in triangle_neighbors]),
        ))
        rows = np.r_[pairs[:, 0], pairs[:, 1]]
        cols = np.r_[pairs[:, 1], pairs[:, 0]]
    else:
        rows = cols = np.empty(0, dtype=np.int64)
    graph = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(triangles), len(triangles))).tocsr()
    _, labels = connected_components(graph, directed=False)
    return labels


def _mesh_component_sizes(node_count: int, triangles: np.ndarray) -> np.ndarray:
    """Return triangle counts for edge-connected components, including unused nodes."""
    triangle_component_sizes = np.bincount(_mesh_component_labels(triangles))
    used_node_count = len(np.unique(triangles))
    unused_node_count = node_count - used_node_count
    if unused_node_count:
        return np.r_[np.sort(triangle_component_sizes)[::-1], np.zeros(unused_node_count, dtype=np.int64)]
    return np.sort(triangle_component_sizes)[::-1]


def _require_connected_mesh(node_count: int, triangles: np.ndarray) -> None:
    """Reject meshes not joined into one component by shared triangle edges."""
    component_sizes = _mesh_component_sizes(node_count, triangles)
    if len(component_sizes) > 1:
        raise ValueError(
            "thin-shell mesh has disconnected element components "
            f"with triangle counts {component_sizes.tolist()}"
        )


def _triangle_element(points_xyz: np.ndarray, triangle: np.ndarray, first: np.ndarray, second: np.ndarray,
                      thickness_m: np.ndarray, young_modulus_pa: np.ndarray,
                      poisson_ratio: float, radius_m: float) -> tuple[
                          np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, np.ndarray
                      ]:
    """Return global-node element stiffness, strain matrix, constitutive matrix and area."""
    ids = triangle
    radial = geometry.normalize(points_xyz[ids].sum(axis=0))
    edge = points_xyz[ids[1]] - points_xyz[ids[0]]
    axis_x = geometry.normalize(edge - np.dot(edge, radial) * radial)
    axis_y = np.cross(radial, axis_x)
    xy = np.column_stack((points_xyz[ids] @ axis_x, points_xyz[ids] @ axis_y)) * radius_m
    x, y = xy[:, 0], xy[:, 1]
    twice_area = (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0])
    area = abs(twice_area) * 0.5
    if area <= np.finfo(float).eps * radius_m**2:
        raise ValueError("thin-shell mesh contains a degenerate triangle")
    b = np.array([
        [y[1] - y[2], 0.0, y[2] - y[0], 0.0, y[0] - y[1], 0.0],
        [0.0, x[2] - x[1], 0.0, x[0] - x[2], 0.0, x[1] - x[0]],
        [x[2] - x[1], y[1] - y[2], x[0] - x[2], y[2] - y[0], x[1] - x[0], y[0] - y[1]],
    ]) / twice_area
    e_h = float(np.mean(thickness_m[ids]))
    e_pa = float(np.mean(young_modulus_pa[ids]))
    scale = e_pa * e_h / (1.0 - poisson_ratio**2)
    constitutive = scale * np.array([
        [1.0, poisson_ratio, 0.0],
        [poisson_ratio, 1.0, 0.0],
        [0.0, 0.0, (1.0 - poisson_ratio) * 0.5],
    ])
    transform = np.zeros((6, 6))
    for local, node in enumerate(ids):
        transform[2 * local, 2 * local:2 * local + 2] = (np.dot(axis_x, first[node]), np.dot(axis_x, second[node]))
        transform[2 * local + 1, 2 * local:2 * local + 2] = (np.dot(axis_y, first[node]), np.dot(axis_y, second[node]))
    stiffness = area * transform.T @ b.T @ constitutive @ b @ transform
    return stiffness, b @ transform, constitutive, area, axis_x, axis_y


def solve_thin_shell(
    points_xyz: np.ndarray,
    areas_m2: np.ndarray,
    applied_force_n: np.ndarray,
    thickness_m: np.ndarray | float,
    young_modulus_pa: np.ndarray | float,
    tensile_strength_resultant_n_per_m: np.ndarray | None = None,
    *,
    extensional_capacity_resultant_n_per_m: np.ndarray | None = None,
    poisson_ratio: float = 0.25,
    radius_m: float = lithosphere.PLANET_RADIUS_M,
    tolerance: float = 1e-9,
    max_iterations: int | None = None,
    fixed_node_mask: np.ndarray | None = None,
    mesh_triangles: np.ndarray | None = None,
) -> ThinShellResult:
    """Solve small-strain membrane equilibrium on a spherical node set.

    `applied_force_n` is an external force vector per node. Its radial component is supported
    by shell-normal pressure and is excluded from in-plane equilibrium. Rigid-rotation load
    is removed because plate-scale torque is handled by the existing motion solve. The
    remaining tangential load is transferred through linear triangular membrane elements.
    By default all edges are free. `fixed_node_mask` optionally clamps selected nodes to zero
    tangent displacement in the plate-motion frame; reactions at those nodes are supplied by
    the constraint and are excluded from the free-node residual.
    The membrane tensor is the elastic stress increment from the zero-load reference, with
    tension positive. Its largest positive principal value is the incremental extensional
    resultant to compare with the column's integrated *additional* strength to normal-fault
    yield. That capacity already incorporates the column's assumed background horizontal stress
    (`K0`). This comparison assumes vertical compression remains the major principal stress and
    does not represent a total prestressed lithosphere state. Supply the capacity as
    `extensional_capacity_resultant_n_per_m`; the old `tensile_strength_resultant_n_per_m`
    argument remains a compatibility alias. Each element uses plane-stress elasticity with
    local effective thickness and modulus.

    By default, outward-facing convex-hull facets define a first-order mesh. Pass explicit
    `mesh_triangles` to preserve a caller's plate outline, including concave boundaries. This is
    a screening diagnostic: a supplied mesh must be element-connected, and disconnected meshes
    are rejected because each component adds rigid modes not handled by the plate-wide projection.
    The solve does not model faults, bending, plasticity, changing topology, or a thermal strength
    envelope.
    """
    points = geometry.normalize(np.asarray(points_xyz, dtype=float))
    areas = np.asarray(areas_m2, dtype=float)
    force = np.asarray(applied_force_n, dtype=float)
    n = len(points)
    thickness = np.broadcast_to(np.asarray(thickness_m, dtype=float), (n,)).copy()
    modulus = np.broadcast_to(np.asarray(young_modulus_pa, dtype=float), (n,)).copy()
    if (tensile_strength_resultant_n_per_m is None) == (extensional_capacity_resultant_n_per_m is None):
        raise ValueError("supply exactly one extensional capacity resultant")
    capacity_input = (
        extensional_capacity_resultant_n_per_m
        if extensional_capacity_resultant_n_per_m is not None
        else tensile_strength_resultant_n_per_m
    )
    strength = np.asarray(capacity_input, dtype=float)
    fixed_nodes = np.zeros(n, dtype=bool) if fixed_node_mask is None else np.asarray(fixed_node_mask, dtype=bool)
    if points.shape != (n, 3) or force.shape != (n, 3) or areas.shape != (n,) or strength.shape != (n,):
        raise ValueError("points, forces, areas, and strength must have matching node counts")
    if fixed_nodes.shape != (n,):
        raise ValueError("fixed_node_mask must have one value per node")
    if np.all(fixed_nodes):
        raise ValueError("thin-shell solve requires at least one unconstrained node")
    if n < 4:
        raise ValueError("thin-shell solve requires at least four non-coplanar nodes")
    if (np.any(~np.isfinite(points)) or np.any(~np.isfinite(force)) or np.any(~np.isfinite(areas))
            or np.any(~np.isfinite(thickness)) or np.any(~np.isfinite(modulus))
            or np.any(~np.isfinite(strength))):
        raise ValueError("thin-shell inputs must be finite")
    if np.any(areas <= 0.0) or np.any(thickness <= 0.0) or np.any(modulus <= 0.0) or np.any(strength <= 0.0):
        raise ValueError("areas, thickness, elastic modulus, and tensile strength must be positive")
    if not 0.0 <= poisson_ratio < 0.5 or radius_m <= 0.0 or tolerance <= 0.0:
        raise ValueError("invalid Poisson ratio, radius, or solver tolerance")

    if mesh_triangles is None:
        triangles = _spherical_surface_triangles(points)
    else:
        raw_triangles = np.asarray(mesh_triangles)
        if raw_triangles.ndim != 2 or raw_triangles.shape[1] != 3 or len(raw_triangles) == 0:
            raise ValueError("mesh_triangles must have shape (T, 3) with at least one element")
        if not np.issubdtype(raw_triangles.dtype, np.integer):
            raise ValueError("mesh_triangles must contain integer node indices")
        triangles = raw_triangles.astype(np.int64, copy=False)
        if np.any(triangles < 0) or np.any(triangles >= n):
            raise ValueError("mesh_triangles contains an out-of-range node index")
        if np.any(np.diff(np.sort(triangles, axis=1), axis=1) == 0):
            raise ValueError("mesh_triangles cannot repeat a node within an element")
    # A plate may have holes, but its membrane elements must form one connected shell.
    # Clipping candidate triangles can sever narrow plate sections, leaving extra rigid
    # modes that the plate-wide torque projection does not remove. Fail closed on that mesh.
    _require_connected_mesh(n, triangles)
    first, second = _tangent_basis(points)
    tangent_force = force - np.sum(force * points, axis=1)[:, None] * points
    force_dofs = np.column_stack((np.sum(tangent_force * first, axis=1), np.sum(tangent_force * second, axis=1))).ravel()

    # Rigid rotations span the membrane stiffness nullspace. Remove their load component
    # using node area as the force-distribution weight, so the remaining torque is zero.
    rotation_modes = np.empty((2 * n, 3))
    xyz_m = points * radius_m
    for axis in range(3):
        omega = np.eye(3)[axis]
        velocity = np.cross(np.broadcast_to(omega, xyz_m.shape), xyz_m)
        rotation_modes[:, axis] = np.column_stack((np.sum(velocity * first, axis=1), np.sum(velocity * second, axis=1))).ravel()
    rotation_weights = np.repeat(areas, 2)
    weighted_modes = rotation_modes * np.sqrt(rotation_weights)[:, None]
    weighted_force = force_dofs / np.sqrt(rotation_weights)
    rotation_coefficients, *_ = np.linalg.lstsq(weighted_modes, weighted_force, rcond=None)
    rigid_dofs = rotation_weights * (rotation_modes @ rotation_coefficients)
    equilibrated_dofs = force_dofs - rigid_dofs

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    element_data = []
    for triangle in triangles:
        stiffness, strain_matrix, constitutive, element_area, axis_x, axis_y = _triangle_element(
            points, triangle, first, second, thickness, modulus, poisson_ratio, radius_m
        )
        dofs = np.array([dof for node in triangle for dof in (2 * node, 2 * node + 1)])
        rows.extend(np.repeat(dofs, 6).tolist())
        cols.extend(np.tile(dofs, 6).tolist())
        values.extend(stiffness.ravel().tolist())
        element_data.append((triangle, dofs, strain_matrix, constitutive, element_area, axis_x, axis_y))
    stiffness = coo_matrix((values, (rows, cols)), shape=(2 * n, 2 * n)).tocsr()
    limit = max_iterations if max_iterations is not None else max(1000, 20 * n)
    free_dofs = ~np.repeat(fixed_nodes, 2)
    reduced_stiffness = stiffness[free_dofs][:, free_dofs]
    diagonal = reduced_stiffness.diagonal()
    if np.any(diagonal <= 0.0):
        raise ValueError("mesh contains unconstrained degrees of freedom with no element stiffness")
    scaling = 1.0 / np.sqrt(diagonal)
    scaled_stiffness = reduced_stiffness.multiply(scaling[:, None]).multiply(scaling[None, :]).tocsr()
    scaled_force = equilibrated_dofs[free_dofs] * scaling
    solve = lsmr(scaled_stiffness, scaled_force, atol=tolerance, btol=tolerance, maxiter=limit)
    displacement_dofs = np.zeros(2 * n)
    displacement_dofs[free_dofs] = scaling * solve[0]
    stop_code = int(solve[1])
    if stop_code not in (0, 1, 2, 4, 5):
        raise RuntimeError(f"thin-shell equilibrium solver did not converge (stop code {stop_code})")

    nodal_stress = np.zeros((n, 3))
    nodal_weight = np.zeros(n)
    for triangle, dofs, strain_matrix, constitutive, element_area, axis_x, axis_y in element_data:
        resultant = constitutive @ strain_matrix @ displacement_dofs[dofs]
        stress_element = np.array([[resultant[0], resultant[2]], [resultant[2], resultant[1]]])
        for node in triangle:
            rotation = np.array([
                [np.dot(first[node], axis_x), np.dot(first[node], axis_y)],
                [np.dot(second[node], axis_x), np.dot(second[node], axis_y)],
            ])
            stress_node = rotation @ stress_element @ rotation.T
            nodal_stress[node] += np.array([stress_node[0, 0], stress_node[1, 1], stress_node[0, 1]]) * element_area
            nodal_weight[node] += element_area
    nodal_stress /= np.maximum(nodal_weight[:, None], np.finfo(float).tiny)
    mean = nodal_stress[:, :2].sum(axis=1) * 0.5
    radius = np.sqrt(((nodal_stress[:, 0] - nodal_stress[:, 1]) * 0.5) ** 2 + nodal_stress[:, 2] ** 2)
    max_tensile = np.maximum(mean + radius, 0.0)
    failure_ratio = max_tensile / strength
    displacement_xyz = displacement_dofs.reshape(n, 2, 1) * np.stack((first, second), axis=1)
    displacement_xyz = displacement_xyz.sum(axis=1)
    rigid_force = rigid_dofs.reshape(n, 2)
    rigid_force_xyz = rigid_force[:, :1] * first + rigid_force[:, 1:] * second
    equilibrated_force_xyz = equilibrated_dofs.reshape(n, 2)
    equilibrated_force_xyz = equilibrated_force_xyz[:, :1] * first + equilibrated_force_xyz[:, 1:] * second
    return ThinShellResult(
        displacement_m=displacement_xyz,
        stress_resultant_n_per_m=nodal_stress,
        max_tensile_resultant_n_per_m=max_tensile,
        tensile_strength_resultant_n_per_m=strength,
        failure_ratio=failure_ratio,
        equilibrated_force_n=equilibrated_force_xyz,
        rigid_rotation_force_n=rigid_force_xyz,
        fixed_node_mask=fixed_nodes,
        triangles=triangles,
        solver_stop_code=stop_code,
        equilibrium_residual_n=float(np.linalg.norm((stiffness @ displacement_dofs - equilibrated_dofs)[free_dofs])),
    )


def solve_force_field(
    force_field,
    thickness_m: np.ndarray | float,
    young_modulus_pa: np.ndarray | float,
    tensile_strength_resultant_n_per_m: np.ndarray | None = None,
    *,
    extensional_capacity_resultant_n_per_m: np.ndarray | None = None,
    **kwargs,
) -> ThinShellResult:
    """Transfer a `torque.PerNodeForceField` through the membrane equilibrium diagnostic."""
    return solve_thin_shell(
        points_xyz=force_field.points_xyz,
        areas_m2=force_field.areas_m2,
        applied_force_n=force_field.total_n,
        thickness_m=thickness_m,
        young_modulus_pa=young_modulus_pa,
        tensile_strength_resultant_n_per_m=tensile_strength_resultant_n_per_m,
        extensional_capacity_resultant_n_per_m=extensional_capacity_resultant_n_per_m,
        **kwargs,
    )
