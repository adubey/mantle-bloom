"""Offline prototype of RFC 0018's damage-band rifting (2026-11-14 revision), for its diagnosis gates.

Not production code: D0 and D1 use it to calibrate and tune the method before #846 implements it.

This module holds the parts D0 needs:

* ``LoadParts`` / ``load_parts``: the cheap directional membrane-load tensor, built from the per-node force
  field the torque engine already computes (``torque.gather_per_node_force_field``). Basal drag enters as a
  signed traction difference along each adjacency, scaled by ``L_drag``. Slab forces (pull plus slab drag)
  enter as tension along their force direction, and ridge push and contact drag as compression, each spread
  inward by an exponential kernel of length ``L_pull``. Tensors are 2x2 symmetric, stored as (xx, yy, xy)
  in each node's ``thin_shell._tangent_basis``, the basis the thin-shell diagnostic reports in.
* ``principal``: largest principal value and its direction.
* ``ligament_splits``: the integrated path-normal breakthrough rule.

Usage: imported by ``calibrate_rift_load_proxy.py`` (D0); see that script.
"""

from dataclasses import dataclass

import numpy as np

from app import lithosphere, thin_shell

RADIUS_M = lithosphere.PLANET_RADIUS_M


@dataclass
class LoadParts:
    """Separable pieces of the proxy: ``drag_unit`` scales linearly with L_drag; the kernel parts depend on L_pull."""

    drag_unit: np.ndarray  # (N, 3) tensor for L_drag = 1 m
    sources: dict  # name -> (positions (S, 3), unit directions (S, 3), forces N (S,), sign +1 tension / -1 compression)


def tangent_dyad(direction: np.ndarray, first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """(xx, yy, xy) of the unit dyad along ``direction`` (projected into each node's tangent plane)."""
    a, b = np.sum(direction * first, axis=1), np.sum(direction * second, axis=1)
    norm = np.maximum(np.hypot(a, b), 1e-15)
    a, b = a / norm, b / norm
    return np.column_stack([a * a, b * b, a * b])


def drag_unit_tensor(points, areas, i, j, basal_drag_n, first, second) -> np.ndarray:
    """Mean over each node's adjacencies of Δτ (t̂⊗t̂), with t̂ the i→j centre line; per metre of L_drag."""
    tau = basal_drag_n / np.maximum(areas, 1.0)[:, None]  # Pa
    out = np.zeros((len(points), 3))
    degree = np.zeros(len(points))
    for a, b in ((i, j), (j, i)):  # each node sees the pair oriented away from itself; Δτ·t̂ is symmetric
        chord = points[b] - points[a]
        dtau = np.sum((tau[b] - tau[a]) * chord, axis=1) / np.maximum(np.linalg.norm(chord, axis=1), 1e-15)
        dyad = tangent_dyad(chord, first[a], second[a])
        np.add.at(out, a, dtau[:, None] * dyad)
        np.add.at(degree, a, 1.0)
    return out / np.maximum(degree, 1.0)[:, None]


def _force_sources(points, force_n, sign):
    magnitude = np.linalg.norm(force_n, axis=1)
    keep = magnitude > 0.0
    return points[keep], force_n[keep] / magnitude[keep, None], magnitude[keep], sign


def load_parts(points, areas, i, j, field) -> LoadParts:
    first, second = thin_shell._tangent_basis(points)
    sources = {
        "slab": _force_sources(points, field.slab_pull_n + field.slab_drag_n, 1.0),
        "ridge": _force_sources(points, field.ridge_push_n, -1.0),
        "contact": _force_sources(points, field.contact_drag_n, -1.0),
    }
    return LoadParts(drag_unit_tensor(points, areas, i, j, field.basal_drag_n, first, second), sources)


def kernel_tensor(points, sources, l_pull_m: float) -> np.ndarray:
    """Σ_s sign F_s exp(-d/L)/(2L) (ê_s⊗ê_s): a line load per unit length summed over its sources, so a long straight
    front of line load f gives about f at the front, decaying inward over L."""
    first, second = thin_shell._tangent_basis(points)
    positions, directions, forces, sign = sources
    out = np.zeros((len(points), 3))
    for start in range(0, len(positions), 256):  # chunk sources to bound memory
        p, e, f = positions[start:start + 256], directions[start:start + 256], forces[start:start + 256]
        d = RADIUS_M * np.arccos(np.clip(points @ p.T, -1.0, 1.0))  # (N, S)
        w = sign * f[None, :] * np.exp(-d / l_pull_m) / (2.0 * l_pull_m)
        for k in range(len(p)):
            out += w[:, k, None] * tangent_dyad(np.broadcast_to(e[k], points.shape), first, second)
    return out


def assemble(parts: LoadParts, kernels: dict, l_drag_m: float) -> np.ndarray:
    total = parts.drag_unit * l_drag_m
    for value in kernels.values():
        total = total + value
    return total


def balanced_tangent_force(points, areas, force_n) -> np.ndarray:
    """Tangent force with its area-weighted rigid-rotation part removed, as `thin_shell.solve_thin_shell` does,
    so the remainder has zero net torque and a line cut through the plate sees only the part that stresses it."""
    tangent = force_n - np.sum(force_n * points, axis=1)[:, None] * points
    xyz = points * RADIUS_M
    modes = np.stack([np.cross(axis, xyz) for axis in np.eye(3)], axis=2)  # (N, 3 xyz, 3 rotation axes)
    root = np.sqrt(np.maximum(areas, 1.0))
    design = (modes * root[:, None, None]).reshape(-1, 3)
    target = (tangent / root[:, None]).ravel()
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    return tangent - areas[:, None] * (modes @ coefficients)


CUT_ANGLES = np.linspace(0.0, np.pi, 6, endpoint=False)
_CUT_DESIGN = np.column_stack([np.cos(CUT_ANGLES) ** 2, np.sin(CUT_ANGLES) ** 2, 2 * np.cos(CUT_ANGLES) * np.sin(CUT_ANGLES)])
_CUT_INVERSE = np.linalg.pinv(_CUT_DESIGN)


def _line_resultants(coords, force_2d, areas, spacing_m, angle) -> np.ndarray:
    """Pull across the straight line through each cell normal to `angle`: the balanced force on its far side along
    the normal, divided by the plate's width there (cell area within half a spacing of the line, over the spacing)."""
    normal = np.array([np.cos(angle), np.sin(angle)])
    position, push = coords @ normal, force_2d @ normal
    order = np.argsort(position, kind="stable")
    ranked = position[order]
    beyond = np.empty(len(position))
    beyond[order] = np.cumsum(push[order][::-1])[::-1] - 0.5 * push[order]  # strictly beyond, plus half its own cell
    cumulative = np.concatenate([[0.0], np.cumsum(areas[order])])
    upper = np.searchsorted(ranked, position + 0.5 * spacing_m, side="right")
    lower = np.searchsorted(ranked, position - 0.5 * spacing_m, side="left")
    width = np.maximum((cumulative[upper] - cumulative[lower]) / spacing_m, spacing_m)
    return beyond / width


def cut_tensor(points, areas, force_n, spacing_m) -> np.ndarray:
    """Proxy C: the membrane tensor fitted to line-cut resultants in six directions, in each node's tangent basis."""
    balanced = balanced_tangent_force(points, areas, force_n)
    centre = points.mean(axis=0)
    centre /= np.linalg.norm(centre)
    e_u, e_v = (axis[0] for axis in thin_shell._tangent_basis(centre[None, :]))
    coords = RADIUS_M * np.column_stack([points @ e_u, points @ e_v])
    force_2d = np.column_stack([balanced @ e_u, balanced @ e_v])
    resultants = np.column_stack([_line_resultants(coords, force_2d, areas, spacing_m, a) for a in CUT_ANGLES])
    frame = resultants @ _CUT_INVERSE.T  # (N, 3): xx, yy, xy in the (e_u, e_v) frame
    first, second = thin_shell._tangent_basis(points)
    r11, r12, r21, r22 = first @ e_u, first @ e_v, second @ e_u, second @ e_v  # rotation from frame to node basis
    xx, yy, xy = frame[:, 0], frame[:, 1], frame[:, 2]
    return np.column_stack([
        r11 * r11 * xx + 2 * r11 * r12 * xy + r12 * r12 * yy,
        r21 * r21 * xx + 2 * r21 * r22 * xy + r22 * r22 * yy,
        r11 * r21 * xx + (r11 * r22 + r12 * r21) * xy + r12 * r22 * yy,
    ])


def principal(tensor: np.ndarray):
    """Largest principal value, smallest principal value and the angle (rad, in the tangent basis) of the largest."""
    xx, yy, xy = tensor[:, 0], tensor[:, 1], tensor[:, 2]
    mean, radius = 0.5 * (xx + yy), np.hypot(0.5 * (xx - yy), xy)
    return mean + radius, mean - radius, 0.5 * np.arctan2(2.0 * xy, xx - yy)


def normal_component(tensor: np.ndarray, normal_ab: np.ndarray) -> np.ndarray:
    """n̂ᵀ𝐍n̂ for in-plane unit normals given as (a, b) in the tangent basis."""
    a, b = normal_ab[:, 0], normal_ab[:, 1]
    return tensor[:, 0] * a * a + tensor[:, 1] * b * b + 2.0 * tensor[:, 2] * a * b


def ligament_splits(tensor, normal_ab, ds, intact, n_eff) -> bool:
    """Σ_path max(0, n̂ᵀ𝐍n̂) ds ≥ Σ_intact N_eff ds; a fully broken path splits while it carries path-normal load."""
    pulling = np.maximum(0.0, normal_component(tensor, normal_ab))
    load = float(np.sum(pulling * ds))
    strength = float(np.sum((n_eff * ds)[intact]))
    return load > 0.0 and load >= strength
