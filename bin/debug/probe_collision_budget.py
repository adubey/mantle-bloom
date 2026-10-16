#!/usr/bin/env python3
"""Collision-resolution gate (issue #452): where overlap builds at continental contacts, and
what force would stop it. Part of the exploration #449; no model change.

Replays one seed from generation (`world.generate_world(seed)`, defaults: the #403 controls)
and taps the step in this process only; nothing is fed back. Every step, for every contact
pair in scope (continental on both sides, or craton within `--reach-km` of the contact, as
`probe_craton_contacts.Probe` reads them):

- **G1, overlap.** Area of each plate's nodes inside the other (the `contested` containment
  test, which `plates.compute_node_overlap` shares), penetration depth (each overlapping
  node's distance to its own plate's nearest non-overlapping node, in km and cells), and the
  double-counted volume D: both columns' Hc over the overlap ground (half of each side's sum
  of own plus facing column, so the ground is counted once). Each overlapping node is
  **craton**-blocked (a retreat candidate holding craton, or craton on both sides: RFC
  0008's jam test), **cap**-blocked (its column or the facing one within CAP_TOL_M of
  `SUTURE_ACCRETION_MAX_HC_M`), or **deformable**. D's change from one prepass to the next
  is booked node by node, by (plate, cell key), to the node's class.
- **G2/G3, forces.** At the shift each plate solves `(I + gK) w = I w_old + g (tau + b)`
  (`torque.integrate_omega`). The probe keeps each plate's system and its contact drag per
  neighbour, and solves the pair jointly (6x6), with the drag's coupling to the other
  plate's w moved to the left-hand side. A force f per metre of trench, normal to the
  contact and pushing the plates apart, enters as a torque on each side's edge nodes
  (`CONTACT_EDGE_SPACINGS`). The **stall force** is the f that holds the edge's mean normal
  closing at zero. Basal drag, slab pull and drag, ridge push and the plates' other contacts
  stay in each system as the model has them, so the plate's drive needs no apportioning
  between its boundaries by hand. The slab-pull share is the drop in stall force with both
  plates' slab-pull torques removed (#315). `mantle.clamp_rate` is left out.
  The **push** is the gravitational potential energy contrast between the contact band's
  continental columns (both sides' edge nodes, `isostatic_elevation`'s own Airy columns with
  Hm, unclipped) and the reference continental column, integrated to a common compensation
  depth. `push_airy` is #449's crust-only Airy figure, for comparison. H* is the uniform band
  Hc, at the band's mean Hm, whose push equals the stall force. The steady closing with the
  push added keeps only the normal component of the contact drag (#449 review, point 4).
- **G4, budget.** Per pair-step: crust driven in (edge nodes' closing x dt x length share x
  Hc on both sides, in D's units), D's change over the shift, D removed with cells the
  deform pass consumed (retreat), band thickening over deform, and the deform pass's sinks
  (`no_outlet_subducted_m3`, `ceiling_overflow_no_outlet_m3`, `collision_subducted_m3`) and
  `suture_donated_m3`, each plate's split over its contacts by contact nodes. Craton sinks
  near contacts come from `probe_craton_contacts`' rows, written alongside.
- **G5 and checks.** At each merge: the absorbed plate's area and Hc over third plates, on
  their continental and on their oceanic crust, with the overlap's age
  (`overlap_onset_years`); Hc and craton the merge lost; |w_a - w_b|. Splits of plates in a
  jammed pair. Pair symmetry: the two sides' model contact-drag normal forces.

`--capture-dir` saves each merge with RFC 0001's signature (`probe_craton_clips.tap_capture`)
for `measure_merge_placement.py`. `--untapped` runs the same loop with no taps and writes only
the end craton volume, for the check that the probe doesn't change the run.

Run seeds one at a time (each holds a full world, ~1.5 GB):

    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed 326494039 \\
        --out /tmp/i452/seed326494039.jsonl --capture-dir ~/mantle-bloom-captures/issue452
    backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed 326494039 --untapped \\
        --out /tmp/i452/untapped-326494039.jsonl
    backend/.venv/bin/python bin/debug/probe_collision_budget.py \\
        --summarize /tmp/i452/seed*.jsonl --json analysis/issue452/gate.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app import collision_polarity, cratons, geometry, lithosphere, merge_split, quad_tectonics, torque  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import MAX_ELEVATION_M, effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M  # noqa: E402

import probe_craton_clips  # noqa: E402
import probe_craton_contacts  # noqa: E402

R = lithosphere.PLANET_RADIUS_M
SPY = torque.SECONDS_PER_YEAR
KM3 = 1e9
CAP_M = SUTURE_ACCRETION_MAX_HC_M
CAP_TOL_M = 50.0
COMPENSATION_DEPTH_M = 400_000.0  # below sea level, under every column's lithosphere
CLASSES = ("craton", "cap", "deform")
CRATON, CAP, DEFORM = 0, 1, 2
OROGENY_ACCOUNTS = ("no_outlet_subducted_m3", "ceiling_overflow_no_outlet_m3", "suture_donated_m3")
BUDGET_ACCOUNTS = OROGENY_ACCOUNTS + ("collision_subducted_m3",)
PRE_MERGE_WINDOW_MYR = 30.0
STRENGTH_N_PER_M = (1e12, 1e13)  # published integrated lithospheric strength (#452 G3)


def _pair(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _key(pair: tuple[int, int]) -> str:
    return f"{pair[0]}-{pair[1]}"


def _r(x, digits: int = 4):
    return None if x is None else float(f"{x:.{digits}g}")


# --- columns and push ----------------------------------------------------------------------


def unclipped_elevation(hc: np.ndarray, hm: np.ndarray, rho_c) -> np.ndarray:
    """`lithosphere.isostatic_elevation` without its final clip to the elevation bounds."""
    rho_a = lithosphere.RHO_ASTHENOSPHERE
    bracket = hc * (1.0 - rho_c / rho_a) + hm * (1.0 - lithosphere.RHO_LITHOSPHERE_MANTLE / rho_a) + lithosphere.ISOSTATIC_REFERENCE_OFFSET_M
    return np.where(bracket <= 0.0, bracket * rho_a / (rho_a - lithosphere.RHO_WATER), bracket)


def column_gpe(hc, hm, rho_c) -> np.ndarray:
    """Vertically integrated lithostatic pressure (N/m) from each column's top down to
    COMPENSATION_DEPTH_M below sea level: water if submerged, crust, mantle lithosphere, then
    asthenosphere. The difference between two columns is the horizontal force per unit length
    their contrast exerts (the columns balance at depth, so the isostatic offset cancels)."""
    hc, hm = np.asarray(hc, dtype=float), np.asarray(hm, dtype=float)
    e = unclipped_elevation(hc, hm, rho_c)
    water = np.clip(-e, 0.0, None)
    asthenosphere = COMPENSATION_DEPTH_M + e - hc - hm
    g = lithosphere.GRAVITY_M_S2
    layers = ((water, lithosphere.RHO_WATER), (hc, rho_c), (hm, lithosphere.RHO_LITHOSPHERE_MANTLE), (asthenosphere, lithosphere.RHO_ASTHENOSPHERE))
    p_top, total = np.zeros_like(hc), np.zeros_like(hc)
    for h, rho in layers:
        total += p_top * h + 0.5 * rho * g * h**2
        p_top = p_top + rho * g * h
    return total


REFERENCE_GPE = float(
    column_gpe([lithosphere.REFERENCE_HC_CONTINENTAL_M], [lithosphere.REFERENCE_HM_CONTINENTAL_M], lithosphere.RHO_CONTINENTAL_CRUST)[0]
)


def push_n_per_m(hc, hm) -> float:
    """Mean GPE contrast of continental columns against the reference continental column."""
    if len(hc) == 0:
        return 0.0
    return float(np.mean(column_gpe(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST)) - REFERENCE_GPE)


def airy_push_n_per_m(hc_m: float) -> float:
    """#449's crust-only Airy figure, 1/2 rho_c g (1 - rho_c/rho_a) (H^2 - H0^2)."""
    rho_c, rho_a = lithosphere.RHO_CONTINENTAL_CRUST, lithosphere.RHO_ASTHENOSPHERE
    h0 = lithosphere.REFERENCE_HC_CONTINENTAL_M
    return 0.5 * rho_c * lithosphere.GRAVITY_M_S2 * (1.0 - rho_c / rho_a) * (hc_m**2 - h0**2)


def h_star_km(stall: float, hm_m: float) -> float | None:
    """Uniform band Hc (km), at Hm `hm_m`, whose push equals `stall`. None when even the cap's
    push falls short (H* above the cap); the reference Hc when a reference band's would do."""
    h0 = lithosphere.REFERENCE_HC_CONTINENTAL_M
    if stall <= push_n_per_m([h0], [hm_m]):
        return h0 / 1e3
    if push_n_per_m([CAP_M], [hm_m]) < stall:
        return None
    lo, hi = h0, CAP_M
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if push_n_per_m([mid], [hm_m]) < stall else (lo, mid)
    return 0.5 * (lo + hi) / 1e3


# --- shift-time systems --------------------------------------------------------------------


class ShiftSystem:
    """One plate's shift solve this step, with its contact drag per neighbour."""

    def __init__(self, plate, inputs, coefficient: np.ndarray, spacing: float) -> None:
        self.plate_id = plate.plate_id
        self.inputs = inputs
        self.coefficient = coefficient
        self.areas = plate.accounting_areas_m2(spacing)
        self.c = coefficient * self.areas * R**2 / SPY
        self.slab = np.zeros(3)
        self.lhs = self.rhs = self.omega = None

    def facing(self, neighbour_id: int) -> np.ndarray:
        return self.inputs.neighbor_plate_id == neighbour_id

    def drag_k(self, neighbour_id: int, normal_only: bool = False) -> np.ndarray:
        """K of the contact drag from the nodes facing `neighbour_id`; its b is K @ w_neighbour
        (`torque.collision_drag_coefficients`). With `normal_only`, the drag of the normal
        closing alone: sum c m m^T, m = p x n."""
        sel = self.facing(neighbour_id) & (self.c > 0.0)
        if not np.any(sel):
            return np.zeros((3, 3))
        p, c = self.inputs.own_points[sel], self.c[sel]
        if normal_only:
            m = np.cross(p, self.inputs.direction_to_neighbor[sel])
            return np.einsum("n,ni,nj->ij", c, m, m)
        return c.sum() * np.eye(3) - np.einsum("n,ni,nj->ij", c, p, p)

    def neighbour_omega(self, neighbour_id: int) -> np.ndarray:
        sel = self.facing(neighbour_id)
        return self.inputs.neighbor_omega[sel][0] if np.any(sel) else np.zeros(3)

    def edge(self, neighbour_id: int, spacing: float) -> np.ndarray:
        return self.facing(neighbour_id) & (self.inputs.dist_to_neighbor <= collision_polarity.CONTACT_EDGE_SPACINGS * spacing)

    def drag_normal_force_n(self, neighbour_id: int, omega_other: np.ndarray) -> float:
        """The model's contact drag on this plate from `neighbour_id`, normal component (N)."""
        sel = self.facing(neighbour_id) & (self.coefficient > 0.0)
        if not np.any(sel) or self.omega is None:
            return 0.0
        p, n = self.inputs.own_points[sel], self.inputs.direction_to_neighbor[sel]
        closing = np.einsum("ni,ni->n", np.cross(self.omega - omega_other, p), n) * R / SPY
        return float(np.sum(self.coefficient[sel] * self.areas[sel] * closing))


def _edge_lever(system: ShiftSystem, sel: np.ndarray, length_m: float) -> tuple[np.ndarray, np.ndarray]:
    """(torque per unit f in N/m spread over the side's edge nodes, sum of p x n)."""
    p, n = system.inputs.own_points[sel], system.inputs.direction_to_neighbor[sel]
    return np.cross(p * R, -n).sum(axis=0) * (length_m / len(p)), np.cross(p, n).sum(axis=0)


def _joint_solve(a: ShiftSystem, b: ShiftSystem, g: float, rhs_a: np.ndarray, rhs_b: np.ndarray, normal_only: bool) -> np.ndarray:
    """Both plates' w from their shift systems, with the pair's contact drag coupling them
    (kept whole, or only its normal part)."""
    k_ab, k_ba = a.drag_k(b.plate_id), b.drag_k(a.plate_id)
    lhs_a, lhs_b, couple_ab, couple_ba = a.lhs, b.lhs, k_ab, k_ba
    if normal_only:
        couple_ab, couple_ba = a.drag_k(b.plate_id, True), b.drag_k(a.plate_id, True)
        lhs_a, lhs_b = lhs_a + g * (couple_ab - k_ab), lhs_b + g * (couple_ba - k_ba)
    big = np.block([[lhs_a, -g * couple_ab], [-g * couple_ba, lhs_b]])
    return np.linalg.solve(big, np.concatenate([rhs_a, rhs_b]))


def pair_forces(a: ShiftSystem, b: ShiftSystem, years: float, spacing: float) -> dict | None:
    """Stall force, its slab-pull share, closing response and the model drag on each side."""
    edge_a, edge_b = a.edge(b.plate_id, spacing), b.edge(a.plate_id, spacing)
    if not np.any(edge_a) or not np.any(edge_b) or a.lhs is None or b.lhs is None:
        return None
    na, nb = int(edge_a.sum()), int(edge_b.sum())
    length = 0.5 * (na + nb) * spacing * R
    g = years * SPY**2
    t_a, m_a = _edge_lever(a, edge_a, length)
    t_b, m_b = _edge_lever(b, edge_b, length)
    m = (m_a - m_b) * R * 100.0 / (na + nb)  # mean normal closing in cm/yr = m . (w_a - w_b)
    # Each plate's rhs less the pair's drag b, which the coupling term carries instead.
    base_a = a.rhs - g * a.drag_k(b.plate_id) @ a.neighbour_omega(b.plate_id)
    base_b = b.rhs - g * b.drag_k(a.plate_id) @ b.neighbour_omega(a.plate_id)

    def closing(x: np.ndarray) -> float:
        return float(m @ (x[:3] - x[3:]))

    c0 = closing(_joint_solve(a, b, g, base_a, base_b, False))
    s = closing(_joint_solve(a, b, g, g * t_a, g * t_b, False))
    c0_noslab = closing(_joint_solve(a, b, g, base_a - g * a.slab, base_b - g * b.slab, False))
    c0_normal = closing(_joint_solve(a, b, g, base_a, base_b, True))
    s_normal = closing(_joint_solve(a, b, g, g * t_a, g * t_b, True))
    fa, fb = a.drag_normal_force_n(b.plate_id, b.omega), b.drag_normal_force_n(a.plate_id, a.omega)
    return {
        "len_km": length / 1e3,
        "close_now_cm": float(m @ (a.omega - b.omega)),
        "close_joint_cm": c0,
        "stall": -c0 / s if s < 0.0 else None,
        "stall_noslab": -c0_noslab / s if s < 0.0 else None,
        "normal": (c0_normal, s_normal),
        "drag_a": fa / length,
        "drag_b": fb / length,
    }


# --- per-plate node reads ------------------------------------------------------------------


class PlateRead:
    """One plate's contact geometry at a moment in the step."""

    def __init__(self, world, plate, plates: list, spacing: float, reach: float, with_classes: bool) -> None:
        cache = world.boundary_search_cache
        neighbours = cache.neighbours(plate, plates, reach) if cache is not None else plate.get_neighbours(plates, reach)
        self.plate = plate
        self.inputs = inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
        _, _, _, self.contested = torque.classify_boundary_nodes(plate, neighbours, inputs, reach, cache)
        self.keys = np.asarray(plate.cell_keys).copy()
        self.areas = plate.accounting_areas_m2(spacing)
        self.hc, self.hm, self.craton = inputs.own_hc, inputs.own_hm, inputs.own_craton_m
        self.owner = inputs.neighbor_plate_id
        self.band = inputs.dist_to_neighbor <= reach
        self.edge = inputs.dist_to_neighbor <= collision_polarity.CONTACT_EDGE_SPACINGS * spacing
        self.continental = effective_is_continental_from_codes(inputs.own_crust_type_codes, plate.crust_type == "continental")
        self.other_hc, self.other_craton = np.zeros(len(self.hc)), np.zeros(len(self.hc))
        for neighbour in neighbours:
            mine = self.owner == neighbour.plate_id
            if np.any(mine):
                idx = inputs.neighbor_node_index[mine]
                self.other_hc[mine] = neighbour.collect("crustal_thickness_m")[idx]
                self.other_craton[mine] = neighbour.collect("craton_crust_m")[idx]
        self.d = 0.5 * self.areas * (self.hc + self.other_hc)
        if with_classes:
            self._classify(world, spacing)

    def _classify(self, world, spacing: float) -> None:
        candidate = collision_polarity.retreat_candidates(world.collision_polarity_frame, self.plate, self.contested)
        craton = (self.craton > 0.0) & (candidate | (self.other_craton > 0.0))
        cap = (self.hc >= CAP_M - CAP_TOL_M) | (self.other_hc >= CAP_M - CAP_TOL_M)
        self.cls = np.where(craton, CRATON, np.where(cap, CAP, DEFORM))
        points = self.inputs.own_points
        free = ~self.contested
        self.depth = np.zeros(len(points))
        if np.any(self.contested) and np.any(free):
            dist, _ = cKDTree(points[free]).query(points[self.contested])
            self.depth[self.contested] = 2.0 * np.arcsin(np.clip(dist / 2.0, 0.0, 1.0)) * R
        self.closing_m_per_s = torque.boundary_closing_rate_m_per_s(self.plate, self.inputs)

    def contact_nodes(self) -> dict[int, int]:
        """Contact nodes (band or contested) per neighbour, the weights a plate's deform-pass
        sinks are split by."""
        contact = (self.band | self.contested) & (self.owner >= 0)
        ids, counts = np.unique(self.owner[contact], return_counts=True)
        return {int(i): int(c) for i, c in zip(ids, counts)}


def _side_overlap(side: PlateRead, other_id: int, years: float, length_m: float) -> dict:
    """One side's overlap, D by class, penetration and crust driven in, facing `other_id`."""
    ov = side.contested & (side.owner == other_id)
    edge = side.edge & (side.owner == other_id)
    n_edge = int(edge.sum())
    driven = np.clip(side.closing_m_per_s[edge], 0.0, None) * SPY * years * (length_m / max(n_edge, 1)) * side.hc[edge]
    return {
        "ov": ov,
        "area": float(side.areas[ov].sum()),
        "d_cls": [float(side.d[ov & (side.cls == k)].sum()) for k in range(3)],
        "depth": side.depth[ov],
        "edge_cls": [int((edge & (side.cls == k)).sum()) for k in range(3)],
        "in": float(driven.sum()),
        "nodes": {(side.plate.plate_id, int(k)): (float(d), int(c)) for k, d, c in zip(side.keys[ov], side.d[ov], side.cls[ov])},
        "band_keys": side.keys[side.band & (side.owner == other_id)],
        "band_hc": side.hc[side.band & (side.owner == other_id)],
        "edge_cont": edge & side.continental,
    }


def _growth_by_class(now: dict, before: dict) -> list[float]:
    """D's change since the last prepass, node by node, to each node's class (a node gone
    since is booked to the class it had)."""
    out = [0.0, 0.0, 0.0]
    for key, (d, cls) in now.items():
        out[cls] += d - before.get(key, (0.0, cls))[0]
    for key, (d, cls) in before.items():
        if key not in now:
            out[cls] -= d
    return out


# --- the probe -----------------------------------------------------------------------------


class Probe:
    def __init__(self, world, contacts: probe_craton_contacts.Probe) -> None:
        self.world = world
        self.contacts = contacts
        self.spacing = line_spacing_rad(world.node_density)
        self.reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * self.spacing
        self.systems: dict[int, ShiftSystem] = {}
        self.slab: dict[int, np.ndarray] = {}
        self.prev_nodes: dict[str, dict] = {}
        self.pre: dict[str, dict] = {}
        self.weights: dict[int, dict[int, int]] = {}
        self.deform_budget: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
        self.record: dict | None = None
        self.rows: list[dict] = []
        self.ledgers = self._ledgers()

    # --- shift taps

    def on_slab_pull(self, plate, value: np.ndarray) -> None:
        self.slab[plate.plate_id] = np.asarray(value, dtype=float)

    def on_collision_drag(self, plate, inputs, mask, craton_m, jam) -> None:
        coefficient = torque.contact_drag_pa_s_per_m(plate, inputs, mask, self.spacing, craton_m, jam)
        system = ShiftSystem(plate, inputs, coefficient, self.spacing)
        system.slab = self.slab.pop(plate.plate_id, np.zeros(3))
        self.systems[plate.plate_id] = system

    def on_integrate(self, plate, explicit, drag_b, drag_k, inertia, years, result) -> None:
        system = self.systems.get(plate.plate_id)
        if system is None:
            return
        g = years * SPY**2
        system.lhs = inertia + g * drag_k
        system.rhs = inertia @ plate.omega + g * (explicit + drag_b)
        system.omega = np.asarray(result, dtype=float)

    # --- prepass

    def on_prepass(self, years: float) -> None:
        world = self.world
        self.years = years
        plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
        scope = {k: row for k, row in self.contacts.record["pairs"].items() if row["cc"] or sum(row["reach_km2"]) > 0.0}
        ids = {int(pid) for key in scope for pid in key.split("-")}
        reads = {p.plate_id: PlateRead(world, p, plates, self.spacing, self.reach, True) for p in plates if p.plate_id in ids}
        self.weights = {pid: read.contact_nodes() for pid, read in reads.items()}
        self.record = {"myr": world.elapsed_years / 1e6, "pairs": {}, "merges": [], "splits": []}
        self.pre = {}
        nodes_now = {}
        for key, row in scope.items():
            a, b = (int(pid) for pid in key.split("-"))
            out, nodes = self._pair_row(reads[a], reads[b], row)
            out["dD_cls"] = [round(x / KM3, 2) for x in _growth_by_class(nodes, self.prev_nodes.get(key, {}))]
            nodes_now[key] = nodes
            self.record["pairs"][key] = out
        self.prev_nodes = nodes_now
        self.systems = {}
        self.deform_budget.clear()

    def _pair_row(self, ra: PlateRead, rb: PlateRead, contact_row: dict) -> tuple[dict, dict]:
        a, b = ra.plate.plate_id, rb.plate.plate_id
        length = contact_row["len_km"] * 1e3
        sa, sb = _side_overlap(ra, b, self.years, length), _side_overlap(rb, a, self.years, length)
        self.pre[_key((a, b))] = {"a": sa, "b": sb}
        d_cls = [x + y for x, y in zip(sa["d_cls"], sb["d_cls"])]
        depth = np.concatenate([sa["depth"], sb["depth"]])
        edge_cls = [x + y for x, y in zip(sa["edge_cls"], sb["edge_cls"])]
        weights = d_cls if sum(d_cls) > 0.0 else edge_cls
        row = {
            "jammed": contact_row["jammed"],
            "cc": contact_row["cc"],
            "len_km": contact_row["len_km"],
            "close_cm": contact_row["close_cm"],
            "ov_km2": round(0.5 * (sa["area"] + sb["area"]) / 1e6, 1),
            "D_km3": round(sum(d_cls) / KM3, 2),
            "D_cls_km3": [round(x / KM3, 2) for x in d_cls],
            "class": CLASSES[int(np.argmax(weights))] if sum(weights) > 0 else "none",
            "depth_km": [round(float(depth.mean()) / 1e3, 1), round(float(depth.max()) / 1e3, 1)] if len(depth) else [0.0, 0.0],
            "in_km3": round((sa["in"] + sb["in"]) / KM3, 2),
        }
        if row["ov_km2"] > 0.0:
            row["craton_ov_km3"] = round(
                (probe_craton_clips.craton_overlap_m3(self.world, ra.plate, rb.plate) + probe_craton_clips.craton_overlap_m3(self.world, rb.plate, ra.plate)) / KM3, 1
            )
        row.update(self._force_row(ra, rb, sa, sb))
        return row, {**sa["nodes"], **sb["nodes"]}

    def _force_row(self, ra: PlateRead, rb: PlateRead, sa: dict, sb: dict) -> dict:
        hc = np.concatenate([ra.hc[sa["edge_cont"]], rb.hc[sb["edge_cont"]]])
        hm = np.concatenate([ra.hm[sa["edge_cont"]], rb.hm[sb["edge_cont"]]])
        out = {}
        if len(hc):
            push = push_n_per_m(hc, hm)
            out.update(
                band_hc_km=round(float(hc.mean()) / 1e3, 1), band_hm_km=round(float(hm.mean()) / 1e3, 1), push=_r(push),
                push_airy=_r(float(np.mean(airy_push_n_per_m(hc)))),
                clip_share=round(float(np.mean(unclipped_elevation(hc, hm, lithosphere.RHO_CONTINENTAL_CRUST) > MAX_ELEVATION_M)), 3),
            )
        a, b = self.systems.get(ra.plate.plate_id), self.systems.get(rb.plate.plate_id)
        forces = pair_forces(a, b, self.years, self.spacing) if a is not None and b is not None else None
        if forces is None:
            return out
        stall = forces["stall"]
        c_n, s_n = forces["normal"]
        out.update(
            close_now_cm=round(forces["close_now_cm"], 3), close_joint_cm=round(forces["close_joint_cm"], 3),
            stall=_r(stall), stall_noslab=_r(forces["stall_noslab"]),
            drag=[_r(forces["drag_a"]), _r(forces["drag_b"])],
        )
        if stall is not None and stall > 0.0 and len(hc):
            out["h_star_km"] = _r(h_star_km(stall, float(hm.mean())))
            out["close_push_cm"] = round(max(c_n + out["push"] * s_n, 0.0), 3)
        return out

    # --- deform taps

    def on_deform(self, plate, run) -> None:
        before = self._ledgers()
        run()
        after = self._ledgers()
        for account in BUDGET_ACCOUNTS:
            self.deform_budget[plate.plate_id][account] += after[account] - before[account]

    def on_finish_deform(self) -> None:
        if self.record is None:
            return
        world = self.world
        plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
        ids = {int(pid) for key in self.record["pairs"] for pid in key.split("-")}
        by_id = {p.plate_id: p for p in plates}
        reads = {pid: PlateRead(world, by_id[pid], plates, self.spacing, self.reach, False) for pid in ids if pid in by_id}
        for key, row in self.record["pairs"].items():
            a, b = (int(pid) for pid in key.split("-"))
            pre = self.pre[key]
            row.update(self._deform_row(reads.get(a), reads.get(b), pre["a"], pre["b"]))
            for pid, other in ((a, b), (b, a)):
                weights = self.weights.get(pid, {})
                share = weights.get(other, 0) / max(sum(weights.values()), 1)
                for account in BUDGET_ACCOUNTS:
                    row.setdefault("sinks_km3", {}).setdefault(account, 0.0)
                    row["sinks_km3"][account] = round(row["sinks_km3"][account] + share * self.deform_budget[pid][account] / KM3, 2)

    def _deform_row(self, ra: PlateRead | None, rb: PlateRead | None, sa: dict, sb: dict) -> dict:
        d_after, consumed, thick = 0.0, 0.0, 0.0
        for read, side, other in ((ra, sa, rb), (rb, sb, ra)):
            if read is None:
                consumed += sum(d for d, _ in side["nodes"].values())
                continue
            other_id = other.plate.plate_id if other is not None else -2
            ov = read.contested & (read.owner == other_id)
            d_after += float(read.d[ov].sum())
            keys = set(read.keys.tolist())
            consumed += sum(d for (_, k), (d, _) in side["nodes"].items() if k not in keys)
            order = np.argsort(read.keys)
            pos = np.clip(np.searchsorted(read.keys[order], side["band_keys"]), 0, len(order) - 1)
            found = read.keys[order[pos]] == side["band_keys"]
            idx = order[pos[found]]
            thick += float(((read.hc[idx] - side["band_hc"][found]) * read.areas[idx]).sum())
        return {"D_after_km3": round(d_after / KM3, 2), "consumed_D_km3": round(consumed / KM3, 2), "thicken_km3": round(thick / KM3, 2)}

    # --- topology taps

    def on_merge(self, keep, absorb, run) -> None:
        world = self.world
        rel = np.linalg.norm(np.asarray(keep.omega) - np.asarray(absorb.omega)) * 180.0 / np.pi * 1e6  # deg/Myr
        third = self._third_plate_overflow(keep, absorb)
        hc_before = _hc_volume(keep) + _hc_volume(absorb)
        before = self._ledgers()
        run()
        after = self._ledgers()
        survivor = next(p for p in world.plates if p.plate_id == keep.plate_id)
        self.record["merges"].append({
            "pair": _key(_pair(keep.plate_id, absorb.plate_id)), "keep": keep.plate_id, "absorb": absorb.plate_id,
            "rel_deg_per_myr": round(rel, 4), "hc_lost_km3": round((hc_before - _hc_volume(survivor)) / KM3, 1),
            "merge_capped_km3": round((after["merge_capped_m3"] - before["merge_capped_m3"]) / KM3, 2),
            "material_clipped_km3": round((after["collision_subducted_m3"] - before["collision_subducted_m3"]) / KM3, 2),
            **third,
        })

    def _third_plate_overflow(self, keep, absorb) -> dict:
        """The absorbed plate's area and Hc over the keeper and over third plates (their
        continental and oceanic crust), and how long it had sat there."""
        points, _ = absorb.all_points_and_elevation()
        areas, hc = absorb.accounting_areas_m2(self.spacing), absorb.collect("crustal_thickness_m")
        onset = np.asarray(absorb.collect("overlap_onset_years"), dtype=float)
        age = np.where(onset > 0.0, self.world.elapsed_years - onset, 0.0) / 1e6
        out = {}
        covered = np.zeros(len(points), dtype=bool)
        for other in self.world.plates:
            if other.plate_id == absorb.plate_id or other.node_count() == 0:
                continue
            centre, radius = geometry.bounding_sphere(other.all_points_and_elevation()[0])
            near = np.flatnonzero(~covered & (geometry.angular_distance(points, centre) <= radius + 2.0 * self.spacing))
            if not len(near):
                continue
            inside = near[other.contains_batch(points[near])]
            if not len(inside):
                continue
            covered[inside] = True
            if other.plate_id == keep.plate_id:
                _add_cover(out, "on_keep", inside, areas, hc, age)
                continue
            _, idx = other.get_node_kdtree().query(points[inside])
            codes = other.collect("crust_type_code")[idx]
            continental = effective_is_continental_from_codes(codes, other.crust_type == "continental")
            _add_cover(out, "third_cont", inside[continental], areas, hc, age)
            _add_cover(out, "third_ocean", inside[~continental], areas, hc, age)
        return out

    def on_split(self, plate, result) -> None:
        if result is None or self.record is None:
            return
        pid = plate.plate_id
        pairs = [row for key, row in self.record["pairs"].items() if pid in (int(x) for x in key.split("-"))]
        self.record["splits"].append({
            "plate": pid, "jammed": any(r["jammed"] for r in pairs), "overlapping": any(r["ov_km2"] > 0.0 for r in pairs),
        })

    def on_end_topology(self) -> None:
        if self.record is None:
            return
        now = self._ledgers()
        self.record["step_km3"] = {k: round((now[k] - self.ledgers[k]) / KM3, 3) for k in now}
        self.ledgers = now
        self.rows.append(self.record)
        self.record = None

    def _ledgers(self) -> dict[str, float]:
        world = self.world
        budget = getattr(world, "orogenic_relief_budget", None) or {}
        out = {k: float(budget.get(k, 0.0)) for k in OROGENY_ACCOUNTS}
        material = getattr(world, "continental_material_ledger", None) or {}
        out["collision_subducted_m3"] = float(material.get("collision_subducted_m3", 0.0))
        ledger = getattr(world, "craton_ledger", None) or {}
        out["merge_capped_m3"] = float(ledger.get("merge_capped_m3", 0.0))
        out["craton_sinks_m3"] = cratons.booked_sinks_m3(world)
        return out


def _hc_volume(plate) -> float:
    return float(plate.collect("crustal_thickness_m") @ plate.node_areas_m2())


def _add_cover(out: dict, name: str, idx: np.ndarray, areas: np.ndarray, hc: np.ndarray, age: np.ndarray) -> None:
    if not len(idx):
        return
    entry = out.setdefault(name, {"km2": 0.0, "hc_km3": 0.0, "age_myr_max": 0.0, "_age_area": 0.0})
    entry["km2"] = round(entry["km2"] + float(areas[idx].sum()) / 1e6, 1)
    entry["hc_km3"] = round(entry["hc_km3"] + float(hc[idx] @ areas[idx]) / KM3, 1)
    entry["age_myr_max"] = round(max(entry["age_myr_max"], float(age[idx].max())), 1)
    entry["_age_area"] += float(age[idx] @ areas[idx]) / 1e6
    entry["age_myr_mean"] = round(entry["_age_area"] / entry["km2"], 1) if entry["km2"] else 0.0


# --- installing the taps -------------------------------------------------------------------


def install(probe: Probe) -> None:
    """Wrap the step's functions (after `probe_craton_contacts.install` and the capture tap,
    so each wrapper here sees the call first)."""
    slab_pull_torque = torque.slab_pull_torque
    collision_drag_coefficients = torque.collision_drag_coefficients
    integrate_omega = torque.integrate_omega
    observe_contacts = collision_polarity.observe_contacts
    finish_deform_pass = collision_polarity.finish_deform_pass
    end_topology = collision_polarity.end_topology
    deform = quad_tectonics.deform
    merge_with_audit = merge_split._merge_plates_with_craton_audit
    maybe_split_plate = merge_split.maybe_split_plate

    def tapped_slab_pull_torque(plate, inputs, mask, spacing_rad):
        value = slab_pull_torque(plate, inputs, mask, spacing_rad)
        probe.on_slab_pull(plate, value)
        return value

    def tapped_collision_drag_coefficients(plate, inputs, mask, spacing_rad, craton_m=None, jam=None):
        probe.on_collision_drag(plate, inputs, mask, craton_m, jam)
        return collision_drag_coefficients(plate, inputs, mask, spacing_rad, craton_m, jam)

    def tapped_integrate_omega(plate, explicit_torque, drag_b, drag_k, inertia_tensor, years):
        result = integrate_omega(plate, explicit_torque, drag_b, drag_k, inertia_tensor, years)
        probe.on_integrate(plate, explicit_torque, drag_b, drag_k, inertia_tensor, years, result)
        return result

    def tapped_observe_contacts(world, years):
        frame = observe_contacts(world, years)
        probe.on_prepass(years)
        return frame

    def tapped_deform(plate, world, other_plates, years, max_distance):
        probe.on_deform(plate, lambda: deform(plate, world, other_plates, years, max_distance))

    def tapped_finish_deform_pass(world):
        probe.on_finish_deform()
        return finish_deform_pass(world)

    def tapped_end_topology(world, frames_before):
        end_topology(world, frames_before)
        probe.on_end_topology()

    def tapped_merge_with_audit(world, id_keep, id_absorb):
        keep = next(p for p in world.plates if p.plate_id == id_keep)
        absorb = next(p for p in world.plates if p.plate_id == id_absorb)
        probe.on_merge(keep, absorb, lambda: merge_with_audit(world, id_keep, id_absorb))

    def tapped_maybe_split_plate(world, plate):
        result = maybe_split_plate(world, plate)
        probe.on_split(plate, result)
        return result

    torque.slab_pull_torque = tapped_slab_pull_torque
    torque.collision_drag_coefficients = tapped_collision_drag_coefficients
    torque.integrate_omega = tapped_integrate_omega
    collision_polarity.observe_contacts = tapped_observe_contacts
    collision_polarity.finish_deform_pass = tapped_finish_deform_pass
    collision_polarity.end_topology = tapped_end_topology
    quad_tectonics.deform = tapped_deform
    merge_split._merge_plates_with_craton_audit = tapped_merge_with_audit
    merge_split.maybe_split_plate = tapped_maybe_split_plate


def run(args) -> None:
    world = world_mod.generate_world(args.seed)
    start_live = cratons.live_volume_m3(world)
    header = {"seed": args.seed, "step_myr": args.step / 1e6, "reach_km": args.reach_km, "start_craton_m3": start_live,
              "commit": probe_craton_clips._commit(), "untapped": args.untapped,
              "cell_km": line_spacing_rad(world.node_density) * R / 1e3}
    probe = contacts = None
    if not args.untapped:
        contacts = probe_craton_contacts.Probe(world, args.reach_km)
        probe_craton_contacts.install(contacts)
        if args.capture_dir is not None:
            probe_craton_clips.tap_capture(args.capture_dir, 0.01, args.max_captures)
        probe = Probe(world, contacts)
        install(probe)
    started = time.time()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as handle:
        handle.write(json.dumps(header) + "\n")
        while world.elapsed_years < args.years - 1.0:
            world_mod.step_world(world, min(args.step, args.years - world.elapsed_years))
            if probe is not None:
                _write_rows(handle, probe, contacts)
                print(f"seed {args.seed} {world.elapsed_years / 1e6:.0f} Myr, {time.time() - started:.0f} s", flush=True)
        if contacts is not None:
            contacts.finish()
            _write_rows(handle, probe, contacts)
        handle.write(json.dumps({"end_craton_m3": cratons.live_volume_m3(world), "ledger": dict(world.craton_ledger), "wall_s": time.time() - started}) + "\n")


def _write_rows(handle, probe: Probe, contacts) -> None:
    for row in probe.rows:
        handle.write(json.dumps({"kind": "pairs", **row}) + "\n")
    for row in contacts.rows:
        handle.write(json.dumps({"kind": "contacts", **row}) + "\n")
    handle.flush()
    probe.rows, contacts.rows = [], []


# --- summary -------------------------------------------------------------------------------


def _pct(values: list[float], qs=(10, 50, 90)) -> dict | None:
    return {f"p{q}": _r(float(np.percentile(values, q)), 3) for q in qs} if values else None


def _share(part: float, whole: float) -> float | None:
    return round(part / whole, 4) if whole else None


def _load(path: Path) -> tuple[dict, list[dict], list[dict], dict]:
    lines = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    header = lines[0]
    tail = lines[-1] if "end_craton_m3" in lines[-1] else {}
    pairs = [r for r in lines[1:] if r.get("kind") == "pairs"]
    contacts = [r for r in lines[1:] if r.get("kind") == "contacts"]
    return header, pairs, contacts, tail


def _merge_windows(rows: list[dict]) -> dict[str, float]:
    """Pair -> myr of its merge (the last, if a pair id pair merged more than once)."""
    out = {}
    for row in rows:
        for merge in row["merges"]:
            out[merge["pair"]] = row["myr"]
    return out


def _g1(rows: list[dict], cell_km: float) -> dict:
    """G1 over every pair in scope, and over continental pairs alone."""
    out = {"all": _g1_shares(rows, False), "cc": _g1_shares(rows, True)}
    depth = [p["depth_km"] for row in rows for p in row["pairs"].values() if p["ov_km2"] > 0.0]
    out["depth_mean_cells"] = _pct([d[0] / cell_km for d in depth])
    out["depth_max_cells"] = _pct([d[1] / cell_km for d in depth])
    out["depth_mean_km"] = _pct([d[0] for d in depth])
    out["cell_km"] = round(cell_km, 1)
    return out


def _g1_shares(rows: list[dict], cc_only: bool) -> dict:
    """Shares of D's growth by node class: all pair-steps, and the 30 Myr before each merge."""
    merged = _merge_windows(rows)
    net, gross, window = np.zeros(3), np.zeros(3), np.zeros(3)
    for row in rows:
        for key, pair in row["pairs"].items():
            if cc_only and not pair["cc"]:
                continue
            growth = np.array(pair["dD_cls"])
            net += growth
            gross += np.clip(growth, 0.0, None)
            t = merged.get(key)
            if t is not None and t - PRE_MERGE_WINDOW_MYR <= row["myr"] <= t:
                window += growth
    out = {}
    for name, v in (("net", net), ("gross", gross), ("pre_merge_30myr", window)):
        total = float(v.sum())
        out[name] = {"total_km3": round(total, 1), **{c: _share(float(v[i]), total) for i, c in enumerate(CLASSES)}}
        out[name]["blocked"] = _share(float(v[CRATON] + v[CAP]), total)
    return out


def _pair_steps(rows: list[dict], cls: tuple[str, ...]) -> list[dict]:
    return [p for row in rows for p in row["pairs"].values() if p["class"] in cls and p.get("stall") is not None and p["stall"] > 0.0]


def _g2(rows: list[dict]) -> dict:
    """At deformable contacts: H* against the cap, push against stall, closing with the push."""
    steps = _pair_steps(rows, ("deform",))
    with_push = [p for p in steps if "h_star_km" in p]
    h_star = [p["h_star_km"] for p in with_push if p["h_star_km"] is not None]
    tibet = [p for p in steps if "band_hc_km" in p and 65.0 <= p["band_hc_km"] <= 75.0]
    slab = [1.0 - p["stall_noslab"] / p["stall"] for p in steps if p.get("stall_noslab") is not None]
    return {
        "pair_steps_closing": len(steps),
        "jammed_pair_steps": sum(p["jammed"] for p in steps),
        "stall_n_per_m": _pct([p["stall"] for p in steps]),
        "h_star_km": _pct(h_star),
        "share_h_star_below_cap": _share(len(h_star), len(with_push)),
        "share_push_exceeds_stall": _share(sum(p["push"] >= p["stall"] for p in with_push), len(with_push)),
        "close_now_cm": _pct([p["close_now_cm"] for p in steps]),
        "close_with_push_cm": _pct([p["close_push_cm"] for p in with_push]),
        "share_stopped_by_push": _share(sum(p["close_push_cm"] <= 0.0 for p in with_push), len(with_push)),
        "band_hc_km": _pct([p["band_hc_km"] for p in with_push]),
        "band_hm_km": _pct([p["band_hm_km"] for p in with_push]),
        "share_band_clipped": _share(sum(p["clip_share"] > 0.0 for p in with_push), len(with_push)),
        "slab_pull_share_of_stall": _pct(slab),
        "tibet_like": {"n": len(tibet), "push": _pct([p["push"] for p in tibet if "push" in p]), "stall": _pct([p["stall"] for p in tibet]),
                       "push_airy": _pct([p["push_airy"] for p in tibet if "push_airy" in p])},
    }


def _g3(rows: list[dict], step_myr: float) -> dict:
    """At blocked contacts: the stall force against published lithospheric strength, and how
    fast their overlap grows."""
    out = {}
    lo, hi = STRENGTH_N_PER_M
    for cls in ("craton", "cap"):
        steps = _pair_steps(rows, (cls,))
        stall = [p["stall"] for p in steps]
        out[cls] = {
            "pair_steps_closing": len(steps),
            "stall_n_per_m": _pct(stall),
            "share_below_1e13": _share(sum(s <= hi for s in stall), len(stall)),
            "share_in_1e12_1e13": _share(sum(lo <= s <= hi for s in stall), len(stall)),
            "share_below_1e12": _share(sum(s < lo for s in stall), len(stall)),
            "D_growth_km3_per_myr": _pct([sum(p["dD_cls"]) / step_myr for p in steps]),
            "overlap_km2": _pct([p["ov_km2"] for p in steps]),
        }
    return out


def _g4(rows: list[dict], contacts: list[dict]) -> dict:
    """Contacts' volume budget and the cap's costs before merges against merges' losses."""
    sums = defaultdict(float)
    for row in rows:
        for p in row["pairs"].values():
            sums["in"] += p["in_km3"]
            sums["D_growth"] += sum(p["dD_cls"])
            sums["consumed_D"] += p.get("consumed_D_km3", 0.0)
            sums["thicken"] += p.get("thicken_km3", 0.0)
            sums["D_change_over_deform"] += p.get("D_after_km3", p["D_km3"]) - p["D_km3"]
            for account, value in p.get("sinks_km3", {}).items():
                sums[f"contact_{account}"] += value
        for account, value in row.get("step_km3", {}).items():
            sums[f"world_{account}"] += value
        for merge in row["merges"]:
            sums["merge_hc_lost"] += merge["hc_lost_km3"]
            sums["merge_capped"] += merge["merge_capped_km3"]
            sums["merge_material_clipped"] += merge["material_clipped_km3"]
    # D's growth over the shift: this prepass against the last deform pass's end.
    for before, row in zip(rows, rows[1:]):
        for key, p in row["pairs"].items():
            if key in before["pairs"] and "D_after_km3" in before["pairs"][key]:
                sums["D_growth_over_shift"] += p["D_km3"] - before["pairs"][key]["D_after_km3"]
    fused = sum(f["lost_m3"] for c in contacts for f in c.get("fusions", [])) / KM3
    near = sum(sum(c.get("loss", {}).values()) for c in contacts) / KM3
    sums["craton_loss_near_contacts_excl_fusion"] = near - fused
    world_cap = sums["world_no_outlet_subducted_m3"] + sums["world_ceiling_overflow_no_outlet_m3"]
    return {
        **{k: round(v, 1) for k, v in sorted(sums.items())},
        "premerge_cap_hc_over_merge_hc": _share(world_cap, sums["merge_hc_lost"]),
        "premerge_craton_near_contacts_over_merge_capped": _share(sums["craton_loss_near_contacts_excl_fusion"], sums["merge_capped"]),
        "share_of_in_left_as_overlap": _share(sums["D_growth"], sums["in"]),
        "share_of_shift_growth_left_as_overlap": _share(sums["D_growth"], sums["D_growth_over_shift"]),
        "share_of_shift_growth_consumed": _share(sums["consumed_D"], sums["D_growth_over_shift"]),
    }


def _jam_episode(rows: list[dict], key: str, t_merge: float) -> dict | None:
    """The pair's values at the first step of the jammed run ending at its merge."""
    series = [(row["myr"], row["pairs"][key]) for row in rows if key in row["pairs"] and row["myr"] <= t_merge]
    if not series or not series[-1][1]["jammed"]:
        return None
    i = len(series) - 1
    while i > 0 and series[i - 1][1]["jammed"]:
        i -= 1
    start, end = series[i], series[-1]
    return {
        "onset_myr": start[0], "jammed_myr": round(t_merge - start[0], 1),
        "ov_km2": [start[1]["ov_km2"], end[1]["ov_km2"]], "D_km3": [start[1]["D_km3"], end[1]["D_km3"]],
        "craton_ov_km3": [start[1].get("craton_ov_km3", 0.0), end[1].get("craton_ov_km3", 0.0)],
    }


def _merges(rows: list[dict]) -> list[dict]:
    out = []
    for row in rows:
        for merge in row["merges"]:
            pair = row["pairs"].get(merge["pair"], {})
            out.append({
                "myr": row["myr"], **{k: v for k, v in merge.items() if k != "pair"}, "pair": merge["pair"],
                "close_cm": pair.get("close_now_cm", pair.get("close_cm")), "jammed": pair.get("jammed"),
                "jam_episode": _jam_episode(rows, merge["pair"], row["myr"]),
            })
    for merge in out:
        for name in ("on_keep", "third_cont", "third_ocean"):
            if name in merge:
                merge[name].pop("_age_area", None)
    return out


def _checks(rows: list[dict], merges: list[dict]) -> dict:
    mismatch = []
    for row in rows:
        for p in row["pairs"].values():
            if p["jammed"] and p.get("drag") and None not in p["drag"]:
                fa, fb = (abs(x) for x in p["drag"])
                if max(fa, fb) > 0.0:
                    mismatch.append(abs(fa - fb) / max(fa, fb))
    splits = [s for row in rows for s in row["splits"]]
    third = {name: sum(m.get(name, {}).get("hc_km3", 0.0) for m in merges) for name in ("third_cont", "third_ocean")}
    return {
        "pair_drag_mismatch_at_jams": _pct(mismatch),
        "splits": len(splits),
        "splits_of_jammed_plates": sum(s["jammed"] for s in splits),
        "splits_of_overlapping_plates": sum(s["overlapping"] for s in splits),
        "merge_rel_rotation_deg_per_myr": _pct([m["rel_deg_per_myr"] for m in merges]),
        "third_plate_hc_km3": third,
        "third_plate_ocean_share": _share(third["third_ocean"], third["third_cont"] + third["third_ocean"]),
        "close_joint_vs_now_cm_abs_err": _pct([abs(p["close_joint_cm"] - p["close_now_cm"]) for row in rows for p in row["pairs"].values() if "close_joint_cm" in p]),
    }


def summarize_seed(path: Path) -> dict:
    header, rows, contacts, tail = _load(path)
    merges = _merges(rows)
    return {
        "seed": header["seed"], "step_myr": header["step_myr"], "commit": header.get("commit"),
        "myr": rows[-1]["myr"] if rows else 0.0,
        "start_craton_km3": round(header["start_craton_m3"] / KM3, 1),
        "end_craton_km3": round(tail.get("end_craton_m3", float("nan")) / KM3, 1),
        "merges": len(merges), "merges_booking_merge_capped": sum(m["merge_capped_km3"] > 0.0 for m in merges),
        "G1": _g1(rows, header.get("cell_km", 100.0)), "G2": _g2(rows), "G3": _g3(rows, header["step_myr"]), "G4": _g4(rows, contacts),
        "checks": _checks(rows, merges),
        "craton_ledger_km3": {k: round(v / KM3, 1) for k, v in tail.get("ledger", {}).items() if v},
        # Per-merge rows only for merges that lost craton to the cap, to keep the JSON compact.
        "merge_events": [m for m in merges if m["merge_capped_km3"] > 0.0],
    }


def summarize_pooled(paths: list[Path]) -> dict:
    """G1-G3, the checks and the merges over several seeds' pair-steps together, each pair
    key prefixed by its seed so merge windows don't collide."""
    rows, merges, seeds, step = [], [], [], None
    for path in paths:
        header, seed_rows, _, _ = _load(path)
        seeds.append(header["seed"])
        step = header["step_myr"]
        tagged = []
        for row in seed_rows:
            tag = f"{header['seed']}:"
            tagged.append({**row, "pairs": {tag + k: v for k, v in row["pairs"].items()},
                           "merges": [{**m, "pair": tag + m["pair"]} for m in row["merges"]]})
        rows.extend(tagged)
        merges.extend(_merges(tagged))
    return {"seeds": seeds, "step_myr": step, "G1": _g1(rows, _load(paths[0])[0].get("cell_km", 100.0)), "G2": _g2(rows),
            "G3": _g3(rows, step), "checks": _checks(rows, merges), "merges": _merge_digest(merges)}


def _merge_digest(merges: list[dict]) -> dict:
    """Merges that lose craton to the cap: jammed at the merge, overlap at jam onset against
    at the merge (#402's ceiling), and the third-plate overflow's age."""
    capped = [m for m in merges if m["merge_capped_km3"] > 0.0]
    episodes = [m["jam_episode"] for m in capped if m["jam_episode"]]
    onset = sum(e["D_km3"][0] for e in episodes)
    at_merge = sum(e["D_km3"][1] for e in episodes)
    return {
        "n": len(merges), "capped": len(capped), "capped_jammed_at_merge": len(episodes),
        "capped_jammed_myr": _pct([e["jammed_myr"] for e in episodes]),
        "capped_D_at_jam_onset_over_at_merge": _share(onset, at_merge),
        "capped_craton_ov_at_jam_onset_over_at_merge": _share(sum(e["craton_ov_km3"][0] for e in episodes), sum(e["craton_ov_km3"][1] for e in episodes)),
        "capped_rel_rotation_deg_per_myr": _pct([m["rel_deg_per_myr"] for m in capped]),
        "capped_close_cm": _pct([m["close_cm"] for m in capped if m["close_cm"] is not None]),
        "capped_third_ocean_age_myr_mean": _pct([m["third_ocean"]["age_myr_mean"] for m in capped if "third_ocean" in m]),
        "capped_third_cont_age_myr_mean": _pct([m["third_cont"]["age_myr_mean"] for m in capped if "third_cont" in m]),
        "capped_on_keep_age_myr_mean": _pct([m["on_keep"]["age_myr_mean"] for m in capped if "on_keep" in m]),
    }


def summarize(paths: list[Path], json_out: Path | None, pools: list[list[Path]] | None = None) -> None:
    seeds = [summarize_seed(path) for path in paths]
    for s in seeds:
        print(json.dumps({k: v for k, v in s.items() if k != "merge_events"}, indent=1))
    out = {"seeds": seeds}
    if pools:
        out["pooled"] = [summarize_pooled(pool) for pool in pools]
        print(json.dumps(out["pooled"], indent=1))
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps(out, indent=1) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--reach-km", type=float, default=500.0)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--capture-dir", type=Path, help="save each merge with RFC 0001's signature here")
    parser.add_argument("--max-captures", type=int, default=6)
    parser.add_argument("--untapped", action="store_true", help="no taps: only the end craton volume")
    parser.add_argument("--summarize", type=Path, nargs="+")
    parser.add_argument("--json", type=Path, help="with --summarize: write the per-seed summary here")
    parser.add_argument("--pool", type=Path, nargs="+", action="append", help="with --summarize: also pool these runs' pair-steps (repeatable)")
    args = parser.parse_args()
    if args.summarize:
        summarize(args.summarize, args.json, args.pool)
        return
    if args.out is None or args.seed is None:
        parser.error("--seed and --out are required unless --summarize is given")
    run(args)


if __name__ == "__main__":
    main()
