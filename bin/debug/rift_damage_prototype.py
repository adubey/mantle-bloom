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
