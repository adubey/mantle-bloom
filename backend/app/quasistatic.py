"""Reduced simultaneous quasi-static plate solve with bounded contact forces.

This module contains the small algebraic pilot for RFC 0012 G2. It accepts assembled
whole-plate drag/torque terms and contact response maps; geometry and constitutive laws
remain owned by their process-specific callers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GRAVITY_M_PER_S2 = 9.81


def column_gpe_anomaly_n_per_m(
    layer_thickness_m: np.ndarray,
    layer_density_kg_per_m3: np.ndarray,
    compensation_depth_m: float,
    *,
    asthenosphere_density_kg_per_m3: float,
    gravity_m_per_s2: float = GRAVITY_M_PER_S2,
) -> float:
    """Column gravitational-potential-energy anomaly relative to asthenosphere (N/m).

    Layers are ordered from the surface downward. Potential is referenced to one shared
    compensation plane at ``compensation_depth_m``; the reference column contains
    asthenosphere at every depth. This is an integrated column contrast, before projecting
    it onto a boundary normal or integrating it along a ridge/collision segment.
    """
    thickness = np.asarray(layer_thickness_m, dtype=float)
    density = np.asarray(layer_density_kg_per_m3, dtype=float)
    if thickness.ndim != 1 or density.shape != thickness.shape:
        raise ValueError("layer thickness and density arrays must be matching vectors")
    if np.any(thickness < 0.0) or not np.all(np.isfinite(thickness)) or not np.all(np.isfinite(density)):
        raise ValueError("column layers must contain finite values and nonnegative thicknesses")
    if not np.isfinite(compensation_depth_m) or compensation_depth_m < float(thickness.sum()):
        raise ValueError("compensation depth must be at or below the full column")
    top = np.concatenate(([0.0], np.cumsum(thickness)[:-1]))
    bottom = top + thickness
    layer_integral = compensation_depth_m * thickness - 0.5 * (bottom**2 - top**2)
    return float(
        gravity_m_per_s2
        * np.sum((density - asthenosphere_density_kg_per_m3) * layer_integral)
    )


@dataclass(frozen=True)
class ContactSolveResult:
    omega: np.ndarray
    contact_force_n_per_m: np.ndarray
    closing_rate: np.ndarray
    residual_torque_nm: np.ndarray
    converged: bool
    sweeps: int


def contact_force_at_capacity(force: np.ndarray, capacity: np.ndarray) -> np.ndarray:
    """Return the solver's shared tolerance-based capacity classification."""
    force = np.asarray(force, dtype=float)
    capacity = np.asarray(capacity, dtype=float)
    if force.shape != capacity.shape:
        raise ValueError("force and capacity must have matching shapes")
    tolerance = max(1e-12, float(capacity.max(initial=0.0)) * 1e-12)
    return force >= capacity - tolerance


def solve_quasistatic_contacts(
    drag_matrix: np.ndarray,
    driving_torque_nm: np.ndarray,
    contact_torque_per_n_per_m: np.ndarray | None = None,
    closing_rate_per_omega: np.ndarray | None = None,
    contact_capacity_n_per_m: np.ndarray | None = None,
    *,
    closing_rate_offset: np.ndarray | None = None,
    rate_tolerance: float = 1e-12,
    max_sweeps: int = 1000,
) -> ContactSolveResult:
    """Solve ``K Ω = τ + B f`` with bounded non-tensile contact forces.

    ``B`` maps a contact line load (N/m) to plate torques (N·m); ``A`` maps plate angular
    rates (rad/yr) to contact-normal closing rates in the caller's chosen rate units. For
    every contact, ``0 <= f <= capacity`` and a sub-capacity force holds the contact at
    non-closing rate. A contact at capacity may continue closing and is marked yielded.

    The Moore-Penrose solve uses the minimum-norm rate when the physical resistance leaves
    a free mode. In particular, zero load and zero resistance returns zero motion; an
    unbalanced torque remains visible in ``residual_torque_nm`` instead of being converted
    into an invented rate.
    """
    matrix = np.asarray(drag_matrix, dtype=float)
    torque = np.asarray(driving_torque_nm, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("drag_matrix must be square")
    if torque.shape != (matrix.shape[0],):
        raise ValueError("driving_torque_nm must match drag_matrix")
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(torque)):
        raise ValueError("motion system must contain finite values")

    if contact_torque_per_n_per_m is None:
        if closing_rate_per_omega is not None or contact_capacity_n_per_m is not None:
            raise ValueError("contact maps and capacities must be supplied together")
        contact_map = np.zeros((len(torque), 0), dtype=float)
        closure_map = np.zeros((0, len(torque)), dtype=float)
        capacity = np.zeros(0, dtype=float)
    else:
        contact_map = np.asarray(contact_torque_per_n_per_m, dtype=float)
        closure_map = np.asarray(closing_rate_per_omega, dtype=float)
        capacity = np.asarray(contact_capacity_n_per_m, dtype=float)
        if contact_map.ndim != 2 or contact_map.shape[0] != len(torque):
            raise ValueError("contact torque map must have one row per motion degree of freedom")
        count = contact_map.shape[1]
        if closure_map.shape != (count, len(torque)) or capacity.shape != (count,):
            raise ValueError("contact maps and capacities have incompatible shapes")
        if np.any(capacity < 0) or not np.all(np.isfinite(contact_map)):
            raise ValueError("contact maps must be finite and capacities nonnegative")

    if not len(capacity):
        omega = np.linalg.lstsq(matrix, torque, rcond=None)[0]
        residual = matrix @ omega - torque
        return ContactSolveResult(omega, capacity, np.zeros(0), residual, True, 0)

    offset = (
        np.zeros(len(capacity), dtype=float)
        if closing_rate_offset is None
        else np.asarray(closing_rate_offset, dtype=float)
    )
    if offset.shape != (len(capacity),) or not np.all(np.isfinite(offset)):
        raise ValueError("closing_rate_offset must contain one finite value per contact")

    # A single response factorization keeps every projected contact update on the same
    # pre-step drag system. Least squares also gives the prescribed zero-rate free mode.
    response = np.linalg.lstsq(matrix, contact_map, rcond=None)[0]
    free_omega = np.linalg.lstsq(matrix, torque, rcond=None)[0]
    compliance = closure_map @ response
    omega = free_omega.copy()
    force = np.zeros(len(capacity), dtype=float)
    converged = False

    for sweep in range(1, max_sweeps + 1):
        for index in range(len(force)):
            rate = float(closure_map[index] @ omega + offset[index])
            slope = float(compliance[index, index])
            # Compliance is a rate per line force and is naturally far below machine
            # epsilon in SI units (m/s per N/m). An absolute epsilon here silently disables
            # every physical contact update even when its sign is well resolved.
            if slope < 0.0:
                old_force = force[index]
                new_force = np.clip(old_force + rate / -slope, 0.0, capacity[index])
                force[index] = new_force
                omega += response[:, index] * (new_force - old_force)
        closing = closure_map @ omega + offset
        force_tolerance = max(1e-12, np.max(capacity) * 1e-12)
        at_zero = force <= force_tolerance
        at_cap = contact_force_at_capacity(force, capacity)
        violation = np.where(at_zero, np.maximum(closing, 0.0), 0.0)
        violation = np.where(at_cap, np.maximum(-closing, 0.0), violation)
        interior = ~(at_zero | at_cap)
        violation = np.where(interior, np.abs(closing), violation)
        if not len(violation) or float(np.max(violation)) <= rate_tolerance:
            converged = True
            break

    closing = closure_map @ omega + offset
    rhs = torque + contact_map @ force
    residual = matrix @ omega - rhs
    return ContactSolveResult(omega, force, closing, residual, converged, sweep)
