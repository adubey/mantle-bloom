#!/usr/bin/env python3
"""Craton jams at plate contacts: the measurement gate of RFC 0008 (issue #395).

Replays one seed from generation at fixed steps and taps the step in this process only, so
no saved state changes. Every step, for every pair of plates in contact (any crust type):

- **contact**: length (the edge row, nodes within 1.5 spacings of the other plate, averaged
  over both sides), length-weighted normal closing rate, and craton area within
  `--reach-km` (RFC 0008's `CRATON_BRAKE_REACH_KM`, 500 km) of the contact on each side;
- **jam**: a contact node is jammed when it is a retreat candidate holding craton, so
  `cratons.retreat_allowed` blocks it, or when it and the neighbour's nearest node are both
  cratonic. The candidates are read before the no-craton gate: on an active polarized front
  the frozen `run_eligible` lower-plate cells, elsewhere `boundary_context`'s `shrinkable`.
  A pair is jammed when at least `MERGE_MIN_CONTACT_NODES` of its nodes are and it isn't
  diverging (closing >= -1 cm/yr); f is its jammed share of contact nodes;
- **merge state**: `collision_progress` before and after `update_collision_progress`, the
  pairs that rolled and passed, drops of a pair still in contact and jammed, each fusion's
  kind (roll or forced), craton lost in it, and the craton-on-craton overlap volume at fusion
  (the absorbed plate's cratonic cells inside the keeper, at the lesser of the two cratons);
- **craton loss near contacts**: the per-cell drop in `craton_crust_m` from the prepass to
  the topology pass (deform, faults) and from the end of topology to the next prepass
  (erosion, the end-of-step audit), cells matched by (plate id, cell key), each booked to the
  nearest contact within reach. Fusion losses are booked to the fused pair; other topology
  losses (splits, failed rifts, removals) to no contact. The step's ledger sinks ride along
  as a check on the per-cell sums;
- **buoyant pull**: the share of slab-pull torque magnitude (`torque.slab_pull_force_n`, as
  applied, so after solution 1's buoyancy, #400) from continental-coded or cratonic
  subducting nodes;
- **plates**: plate count, area-weighted mean surface speed, and the largest plate's share
  of the surface (RFC 0008's bounds on slowing);
- **contact force**: per unit trench length, the collision brake as applied
  (`torque.contact_drag_pa_s_per_m` wherever it is non-zero: contested nodes, #399, jammed
  contacts, #400, and contacts with craton in brake reach, #423), and the gate's estimate of a jam drag at the probe's craton-jammed
  contacts, `CONTACT_RESISTANCE_N_PER_M x j x v / v_ref`, with j those nodes' mean jam factor.

The summary reports, per seed: J (craton destruction within reach of a jammed contact whose
pair didn't fuse within 50 Myr, as a share of all destruction), dropped collision time on
jammed pairs, the buoyant-pull share, closing rates at converging contacts with and without
craton in reach, the jammed share of collision contact length, and the contact forces.

Run seeds one at a time (each replay holds a full world):

    backend/.venv/bin/python bin/debug/probe_craton_contacts.py --seed 7 --years 1e9 \\
        --out /tmp/rfc0008-seed7.jsonl
    backend/.venv/bin/python bin/debug/probe_craton_contacts.py \\
        --summarize /tmp/rfc0008-seed*.jsonl --json analysis/rfc0008/gate.json
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

from app import collision_polarity, cratons, lithosphere, merge_split, quad_tectonics, torque  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402
from app.lithosphere_plate import CONTINENTAL_CONTESTED_RETREAT_MIN_RUN  # noqa: E402

R = lithosphere.PLANET_RADIUS_M
CM_PER_YR = 0.01 / torque.SECONDS_PER_YEAR  # m/s
TRANSFORM_CM = 1.0
EDGE_SPACINGS = 1.5
FUSION_WINDOW_MYR = 50.0


def _pair(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def _key(pair: tuple[int, int]) -> str:
    return f"{pair[0]}-{pair[1]}"


def retreat_candidates(world, plate, others: list, contested: np.ndarray) -> np.ndarray:
    """`boundary_context`'s `shrinkable` before `_retreat`'s no-craton filter, with the
    craton-free `retreatable` of an active polarized front replaced by its `run_eligible`."""
    if plate.crust_type != "continental":
        shrinkable = contested.copy()
    else:
        shrinkable = quad_tectonics.components_of_at_least(plate, contested, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN)
    frame = world.collision_polarity_frame
    masks = frame.masks.get(plate.plate_id) if frame is not None else None
    if masks is None or len(masks.front_id) != len(contested) or not np.array_equal(masks.node_keys, plate.cell_keys):
        return shrinkable
    live = {p.plate_id for p in others}
    active = []
    for front_id, (lower_id, upper_id) in frame.polarity.items():
        if plate.plate_id not in (lower_id, upper_id):
            continue
        other_id = upper_id if lower_id == plate.plate_id else lower_id
        other = frame.masks.get(other_id)
        if other_id not in live or other is None:
            continue
        own_front, other_front = masks.front_id == front_id, other.front_id == front_id
        if not np.any(own_front & masks.continental) or not np.any(other_front & other.continental):
            continue
        lower = masks if lower_id == plate.plate_id else other
        if np.any(lower.retreatable & (lower.front_id == front_id) & lower.continental):
            active.append(front_id)
    if not active:
        return shrinkable
    nodes = np.isin(masks.front_id, active)
    return (shrinkable & ~nodes) | (masks.run_eligible & masks.lower & nodes & masks.continental)


class Probe:
    def __init__(self, world, reach_km: float) -> None:
        self.world = world
        self.spacing = line_spacing_rad(world.node_density)
        self.reach_rad = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * self.spacing
        self.craton_reach_rad = reach_km * 1e3 / R
        self.ds_m = self.spacing * R
        self.record: dict | None = None
        self.pending: dict | None = None
        self.snapshot: dict | None = None
        self.contact_tree = None
        self.contact_owner = np.zeros(0, dtype=np.int64)
        self.contact_pairs: list[tuple[int, int]] = []
        self.rows: list[dict] = []
        self._new_step_scratch()

    def _new_step_scratch(self) -> None:
        self.pull = {"total": 0.0, "buoyant": 0.0, "cratonic": 0.0, "plates": []}
        self.brake: dict[tuple[int, int], list[float]] = defaultdict(list)
        self.forced: tuple[int, int] | None = None
        self.fusions: list[dict] = []
        self.progress: dict = {}
        self.sinks_mark = cratons.booked_sinks_m3(self.world)

    # --- shift-time taps ----------------------------------------------------------------

    def on_slab_pull(self, plate, inputs, mask) -> None:
        if not np.any(mask):
            return
        pts = inputs.own_points[mask]
        force = torque.slab_pull_force_n(plate, inputs, mask, self.spacing)
        torque_mag = np.linalg.norm(np.cross(pts * R, force[:, None] * inputs.direction_to_neighbor[mask]), axis=1)
        continental = effective_is_continental_from_codes(inputs.own_crust_type_codes[mask], plate.crust_type == "continental")
        cratonic = plate.collect("craton_crust_m")[mask] > 0.0
        total = float(torque_mag.sum())
        buoyant = float(torque_mag[continental | cratonic].sum())
        self.pull["total"] += total
        self.pull["buoyant"] += buoyant
        self.pull["cratonic"] += float(torque_mag[cratonic].sum())
        self.pull["plates"].append(round(buoyant / total, 4) if total > 0.0 else 0.0)

    def on_collision_drag(self, plate, inputs, collision_mask, craton, jam) -> None:
        """The brake as applied per node, `torque.contact_drag_pa_s_per_m` (with the jam, #400)
        times the node's area and real relative speed, summed per neighbour."""
        coefficient = torque.contact_drag_pa_s_per_m(plate, inputs, collision_mask, self.spacing, craton, jam)
        mask = coefficient > 0.0
        if not np.any(mask):
            return
        coefficient = coefficient[mask]
        pts = inputs.own_points[mask]
        relative = np.cross(plate.omega - inputs.neighbor_omega[mask], pts) * R / torque.SECONDS_PER_YEAR
        force = coefficient * plate.accounting_areas_m2(self.spacing)[mask] * np.linalg.norm(relative, axis=1)
        owners = inputs.neighbor_plate_id[mask]
        for other in np.unique(owners):
            self.brake[_pair(plate.plate_id, int(other))].append(float(force[owners == other].sum()))

    # --- prepass: contacts and jams ------------------------------------------------------

    def on_prepass(self) -> None:
        world = self.world
        if self.snapshot is not None:
            # Erosion and the end-of-step audit since the last topology pass belong to the
            # previous step's contacts.
            self._book_diff(self.pending, "after_topology")
            self._flush()
        self.sinks_mark = cratons.booked_sinks_m3(world)
        self.record = {"myr": world.elapsed_years / 1e6, "pairs": {}, "loss": defaultdict(float), "loss_far_m3": 0.0}
        cache = world.boundary_search_cache
        plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
        sides: dict[tuple[int, int], dict[int, dict]] = defaultdict(dict)
        cratons_by_plate = {}
        for plate in plates:
            neighbours = cache.neighbours(plate, plates, self.reach_rad)
            inputs = torque.gather_boundary_force_inputs(plate, neighbours, self.spacing, self.reach_rad, cache)
            _, _, _, contested = torque.classify_boundary_nodes(plate, neighbours, inputs, self.reach_rad, cache)
            craton = plate.collect("craton_crust_m")
            cratons_by_plate[plate.plate_id] = (inputs.own_points[craton > 0.0], plate.accounting_areas_m2(self.spacing)[craton > 0.0])
            if not neighbours:
                continue
            band = inputs.dist_to_neighbor <= self.reach_rad
            contact = band | contested
            closing = torque.boundary_closing_rate_m_per_s(plate, inputs) / CM_PER_YR
            continental = effective_is_continental_from_codes(inputs.own_crust_type_codes, plate.crust_type == "continental")
            other_craton = np.zeros(len(craton))
            for neighbour in neighbours:
                mine = inputs.neighbor_plate_id == neighbour.plate_id
                if np.any(mine):
                    other_craton[mine] = neighbour.collect("craton_crust_m")[inputs.neighbor_node_index[mine]]
            candidate = retreat_candidates(world, plate, neighbours, contested)
            jammed = contact & (((candidate & (craton > 0.0))) | ((craton > 0.0) & (other_craton > 0.0)))
            edge = inputs.dist_to_neighbor <= EDGE_SPACINGS * self.spacing
            hc_cont = np.where(continental, inputs.own_hc, 0.0)
            for neighbour in neighbours:
                mine = contact & (inputs.neighbor_plate_id == neighbour.plate_id)
                if not np.any(mine):
                    continue
                edge_mine = mine & edge
                sides[_pair(plate.plate_id, neighbour.plate_id)][plate.plate_id] = {
                    "points": inputs.own_points[mine],
                    "edge": int(edge_mine.sum()),
                    "close_sum": float(closing[edge_mine].sum()),
                    "contact": int(mine.sum()),
                    "contested": int((mine & contested).sum()),
                    "continental": float(continental[mine].mean()),
                    "jam": int((mine & jammed).sum()),
                    "jam_j": float(torque.jam_factor(hc_cont[mine & jammed]).sum()),
                    "nodes": plate.node_count(),
                }
        self._build_pairs(sides, cratons_by_plate)
        self.record["plates"] = self._plate_stats(plates)
        self.snapshot = self._snapshot()

    def _plate_stats(self, plates: list) -> dict:
        """Plate count, area-weighted mean surface speed and the largest plate's share of the
        surface, at the omegas the plates moved with last step: RFC 0008's bounds on slowing."""
        areas = [float(plate.accounting_areas_m2(self.spacing).sum()) for plate in plates]
        moving = 0.0
        for plate, area in zip(plates, areas):
            points, _ = plate.all_points_and_elevation()
            speed = np.linalg.norm(np.cross(plate.omega, points), axis=1) * R * 100.0  # cm/yr
            moving += float(np.dot(speed, plate.accounting_areas_m2(self.spacing)))
        total = sum(areas)
        return {
            "n": len(plates),
            "speed_cm": round(moving / total, 4) if total else 0.0,
            "largest": round(max(areas) / total, 4) if total else 0.0,
        }

    def _build_pairs(self, sides: dict, cratons_by_plate: dict) -> None:
        points, owners = [], []
        self.contact_pairs = []
        for pair, by_side in sorted(sides.items()):
            if len(by_side) < 2:
                continue  # one-sided reach (a sliver); not a contact both plates feel
            a, b = (by_side[pid] for pid in pair)
            edge = a["edge"] + b["edge"]
            length_m = 0.5 * edge * self.ds_m
            close = (a["close_sum"] + b["close_sum"]) / edge if edge else 0.0
            contact_pts = np.concatenate([a["points"], b["points"]])
            tree = cKDTree(contact_pts)
            reach = []
            for pid in pair:
                pts, areas = cratons_by_plate[pid]
                if len(pts) == 0:
                    reach.append(0.0)
                    continue
                dist, _ = tree.query(pts, distance_upper_bound=self.craton_reach_rad)
                reach.append(float(areas[np.isfinite(dist)].sum()) / 1e6)
            jam = a["jam"] + b["jam"]
            contact = a["contact"] + b["contact"]
            jammed = jam >= merge_split.MERGE_MIN_CONTACT_NODES and close >= -TRANSFORM_CM
            j_mean = (a["jam_j"] + b["jam_j"]) / jam if jam else 0.0
            brake = self.brake.get(pair, [])
            row = {
                "len_km": round(length_m / 1e3, 1),
                "close_cm": round(close, 3),
                "cc": bool(a["continental"] >= 0.5 and b["continental"] >= 0.5),
                "reach_km2": [round(x) for x in reach],
                "contact": contact,
                "jam": jam,
                "jammed": bool(jammed),
                "f": round(jam / contact, 4) if contact else 0.0,
                "overlap": round(max(a["contested"], b["contested"]) / max(1, min(a["nodes"], b["nodes"])), 4),
                "j": round(j_mean, 3),
            }
            if length_m > 0.0 and brake:
                row["brake_n_per_m"] = float(np.mean(brake)) / length_m
            if jammed:
                row["proposed_n_per_m"] = torque.CONTACT_RESISTANCE_N_PER_M * j_mean * max(close, 0.0) / torque.CONTACT_RESISTANCE_REF_CM_PER_YR
            self.record["pairs"][_key(pair)] = row
            self.contact_pairs.append(pair)
            points.append(contact_pts)
            owners.append(np.full(len(contact_pts), len(self.contact_pairs) - 1))
        self.contact_tree = cKDTree(np.concatenate(points)) if points else None
        self.contact_owner = np.concatenate(owners) if owners else np.zeros(0, dtype=np.int64)

    # --- craton loss by cell -------------------------------------------------------------

    def _snapshot(self) -> dict:
        out = {}
        for plate in self.world.plates:
            if plate.node_count() == 0:
                continue
            craton = plate.collect("craton_crust_m")
            cratonic = craton > 0.0
            if not np.any(cratonic):
                continue
            points, _ = plate.all_points_and_elevation()
            out[plate.plate_id] = (plate.cell_keys[cratonic].copy(), craton[cratonic], plate.accounting_areas_m2(self.spacing)[cratonic], points[cratonic])
        return out

    def _book_diff(self, record: dict | None, phase: str) -> None:
        if record is None or self.snapshot is None:
            return
        live = {p.plate_id: p for p in self.world.plates}
        drops, where = [], []
        for pid, (keys, craton, areas, points) in self.snapshot.items():
            plate = live.get(pid)
            after = np.zeros(len(keys))
            if plate is not None and plate.node_count() > 0:
                order = np.argsort(plate.cell_keys)
                sorted_keys = plate.cell_keys[order]
                pos = np.clip(np.searchsorted(sorted_keys, keys), 0, len(sorted_keys) - 1)
                found = sorted_keys[pos] == keys
                after[found] = plate.collect("craton_crust_m")[order[pos[found]]]
            drop = np.maximum(craton - after, 0.0) * areas
            lost = drop > 0.0
            drops.append(drop[lost])
            where.append(points[lost])
        drop = np.concatenate(drops) if drops else np.zeros(0)
        record[f"cell_loss_{phase}_m3"] = float(drop.sum())
        if len(drop) == 0:
            return
        owner = np.full(len(drop), -1)
        if self.contact_tree is not None:
            dist, idx = self.contact_tree.query(np.concatenate(where), distance_upper_bound=self.craton_reach_rad)
            near = np.isfinite(dist)
            owner[near] = self.contact_owner[idx[near]]
        record["loss_far_m3"] += float(drop[owner < 0].sum())
        for i in np.unique(owner[owner >= 0]):
            record["loss"][_key(self.contact_pairs[i])] += float(drop[owner == i].sum())

    # --- topology taps -------------------------------------------------------------------

    def on_begin_topology(self) -> None:
        self._book_diff(self.record, "deform")
        self.record["sinks_deform_m3"] = cratons.booked_sinks_m3(self.world) - self.sinks_mark
        self.topology_mark = cratons.booked_sinks_m3(self.world)

    def on_progress(self, before: dict, after: dict, ready: list, rolled: list) -> None:
        self.progress = {"before": before, "after": after, "ready": ready, "rolled": rolled}

    def on_merge(self, keep, absorb, run) -> None:
        world = self.world
        overlap = 0.0
        a_craton = absorb.collect("craton_crust_m")
        if np.any(a_craton > 0.0) and np.any(keep.collect("craton_crust_m") > 0.0):
            pts, _ = absorb.all_points_and_elevation()
            cratonic = np.flatnonzero(a_craton > 0.0)
            inside = keep.contains_batch(pts[cratonic])
            if np.any(inside):
                _, idx = keep.get_node_kdtree().query(pts[cratonic[inside]])
                stacked = np.minimum(a_craton[cratonic[inside]], keep.collect("craton_crust_m")[idx])
                overlap = float(stacked @ absorb.accounting_areas_m2(self.spacing)[cratonic[inside]])
        before = cratons.plate_volume_m3(world, keep) + cratons.plate_volume_m3(world, absorb)
        pair = _pair(keep.plate_id, absorb.plate_id)
        run()
        survivor = next(p for p in world.plates if p.plate_id == keep.plate_id)
        lost = max(before - cratons.plate_volume_m3(world, survivor), 0.0)
        self.fusions.append({"pair": _key(pair), "kind": "forced" if self.forced == pair else "roll", "overlap_m3": overlap, "lost_m3": lost})
        self.record["loss"][_key(pair)] += lost

    def on_end_topology(self) -> None:
        sinks = cratons.booked_sinks_m3(self.world) - self.topology_mark
        fused = sum(f["lost_m3"] for f in self.fusions)
        self.record["sinks_topology_m3"] = sinks
        self.record["loss_far_m3"] += max(sinks - fused, 0.0)
        self.record["fusions"] = self.fusions
        progress = self.progress
        merged = {f["pair"] for f in self.fusions}
        dropped = []
        for pair, years in progress.get("before", {}).items():
            key = _key(pair)
            if pair in progress["after"] or key in merged:
                continue
            row = self.record["pairs"].get(key)
            dropped.append({"pair": key, "myr": years / 1e6, "jammed": bool(row and row["jammed"])})
        self.record["progress"] = {
            "tracked": len(progress.get("after", {})),
            "rolled": [_key(p) for p in progress.get("rolled", [])],
            "passed": [_key(p) for p in progress.get("ready", [])],
            "dropped": dropped,
        }
        self.record["pull"] = {k: (round(v, 6) if isinstance(v, float) else v) for k, v in self.pull.items()}
        self.record["sinks_after_mark"] = cratons.booked_sinks_m3(self.world)
        self.pending, self.record = self.record, None
        self.snapshot = self._snapshot()
        self._new_step_scratch()

    def _flush(self) -> None:
        record = self.pending
        if record is None:
            return
        record["sinks_after_topology_m3"] = cratons.booked_sinks_m3(self.world) - record.pop("sinks_after_mark")
        record["loss"] = dict(record["loss"])
        self.rows.append(record)
        self.pending = None

    def finish(self) -> None:
        self._book_diff(self.pending, "after_topology")
        self._flush()


def install(probe: Probe) -> None:
    observe_contacts = collision_polarity.observe_contacts
    begin_topology = collision_polarity.begin_topology
    end_topology = collision_polarity.end_topology
    slab_pull_torque = torque.slab_pull_torque
    collision_drag_coefficients = torque.collision_drag_coefficients
    update_collision_progress = merge_split.update_collision_progress
    pop_ready_forced_merge = merge_split.pop_ready_forced_merge
    merge_with_audit = merge_split._merge_plates_with_craton_audit

    def tapped_observe_contacts(world, years):
        frame = observe_contacts(world, years)
        probe.on_prepass()
        return frame

    def tapped_begin_topology(world):
        probe.on_begin_topology()
        return begin_topology(world)

    def tapped_end_topology(world, frames_before):
        end_topology(world, frames_before)
        probe.on_end_topology()

    def tapped_slab_pull_torque(plate, inputs, mask, spacing_rad):
        probe.on_slab_pull(plate, inputs, mask)
        return slab_pull_torque(plate, inputs, mask, spacing_rad)

    def tapped_collision_drag_coefficients(plate, inputs, mask, spacing_rad, craton=None, jam=None):
        probe.on_collision_drag(plate, inputs, mask, craton, jam)
        return collision_drag_coefficients(plate, inputs, mask, spacing_rad, craton, jam)

    def tapped_update_collision_progress(world, years):
        before = dict(world.collision_progress)
        ready = update_collision_progress(world, years)
        after = dict(world.collision_progress)
        rolled = [p for p, acc in after.items() if acc >= merge_split._collision_threshold_years(world.seed, p)]
        probe.on_progress(before, after, list(ready), rolled)
        return ready

    def tapped_pop_ready_forced_merge(world):
        forced = pop_ready_forced_merge(world)
        probe.forced = _pair(*forced) if forced is not None else None
        return forced

    def tapped_merge_with_audit(world, id_keep, id_absorb):
        keep = next(p for p in world.plates if p.plate_id == id_keep)
        absorb = next(p for p in world.plates if p.plate_id == id_absorb)
        probe.on_merge(keep, absorb, lambda: merge_with_audit(world, id_keep, id_absorb))

    collision_polarity.observe_contacts = tapped_observe_contacts
    collision_polarity.begin_topology = tapped_begin_topology
    collision_polarity.end_topology = tapped_end_topology
    torque.slab_pull_torque = tapped_slab_pull_torque
    torque.collision_drag_coefficients = tapped_collision_drag_coefficients
    merge_split.update_collision_progress = tapped_update_collision_progress
    merge_split.pop_ready_forced_merge = tapped_pop_ready_forced_merge
    merge_split._merge_plates_with_craton_audit = tapped_merge_with_audit


def run(seed: int, years: float, step: float, reach_km: float, out: Path) -> None:
    world = world_mod.generate_world(seed)
    probe = Probe(world, reach_km)
    install(probe)
    start_live = cratons.live_volume_m3(world)
    started = time.time()
    with out.open("w") as handle:
        handle.write(json.dumps({"seed": seed, "step_myr": step / 1e6, "reach_km": reach_km, "start_craton_m3": start_live}) + "\n")
        while world.elapsed_years < years - 1.0:
            world_mod.step_world(world, min(step, years - world.elapsed_years))
            for row in probe.rows:
                handle.write(json.dumps(row) + "\n")
            handle.flush()
            if probe.rows:
                print(f"seed {seed} {probe.rows[-1]['myr']:.0f} Myr, {len(probe.rows[-1]['pairs'])} contacts, {time.time() - started:.0f} s", flush=True)
            probe.rows = []
        probe.finish()
        for row in probe.rows:
            handle.write(json.dumps(row) + "\n")
        handle.write(json.dumps({"end_craton_m3": cratons.live_volume_m3(world), "ledger": dict(world.craton_ledger), "wall_s": time.time() - started}) + "\n")


# --- summary -------------------------------------------------------------------------------


def _weighted_mean(values: list[float], weights: list[float]) -> float | None:
    total = float(np.sum(weights)) if weights else 0.0
    return float(np.dot(values, weights) / total) if total > 0.0 else None


def _loss_tally(rows: list[dict]) -> dict:
    """Craton destruction: in all, near any contact, near a jammed contact, and near a jammed
    contact whose pair didn't fuse within FUSION_WINDOW_MYR (J's numerator)."""
    fused_at: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        for fusion in row.get("fusions", []):
            fused_at[fusion["pair"]].append(row["myr"])
    out = {"total": 0.0, "near_contact": 0.0, "near_jammed": 0.0, "jammed_unfused": 0.0, "sinks": 0.0}
    for row in rows:
        loss = row.get("loss", {})
        out["total"] += row.get("loss_far_m3", 0.0) + sum(loss.values())
        out["sinks"] += sum(row.get(f"sinks_{phase}_m3", 0.0) for phase in ("deform", "topology", "after_topology"))
        for key, value in loss.items():
            out["near_contact"] += value
            pair = row["pairs"].get(key)
            if not (pair and pair["jammed"]):
                continue
            out["near_jammed"] += value
            if not any(row["myr"] <= t <= row["myr"] + FUSION_WINDOW_MYR for t in fused_at.get(key, [])):
                out["jammed_unfused"] += value
    return out


def _contact_tally(rows: list[dict]) -> dict:
    """Closing rates at converging contacts by crust and craton in reach, the jammed share of
    continental collision length, and contact forces at jammed pairs."""
    close = {k: ([], []) for k in ("cc_craton", "cc_clear", "all_craton", "all_clear")}
    out = {"close": close, "collision_len": 0.0, "jammed_len": 0.0, "brake": [], "proposed": [], "jammed_pair_steps": 0}
    for row in rows:
        for pair in row["pairs"].values():
            out["jammed_pair_steps"] += pair["jammed"]
            if pair["jammed"]:
                if "brake_n_per_m" in pair:
                    out["brake"].append(pair["brake_n_per_m"])
                if pair.get("proposed_n_per_m", 0.0) > 0.0:
                    out["proposed"].append(pair["proposed_n_per_m"])
            if pair["close_cm"] <= TRANSFORM_CM:
                continue
            reach = "craton" if sum(pair["reach_km2"]) > 0.0 else "clear"
            for group in ("cc", "all") if pair["cc"] else ("all",):
                close[f"{group}_{reach}"][0].append(pair["close_cm"])
                close[f"{group}_{reach}"][1].append(pair["len_km"])
            if pair["cc"]:
                out["collision_len"] += pair["len_km"]
                out["jammed_len"] += pair["len_km"] * pair["f"] if pair["jammed"] else 0.0
    return out


def _progress_and_pull(rows: list[dict]) -> dict:
    out = {"dropped": 0, "dropped_jammed": 0, "dropped_jammed_myr": 0.0, "pull": {"total": 0.0, "buoyant": 0.0, "cratonic": 0.0}, "plate_shares": []}
    for row in rows:
        for drop in row.get("progress", {}).get("dropped", []):
            out["dropped"] += 1
            if drop["jammed"]:
                out["dropped_jammed"] += 1
                out["dropped_jammed_myr"] += drop["myr"]
        pull = row.get("pull", {})
        for key in out["pull"]:
            out["pull"][key] += pull.get(key, 0.0)
        out["plate_shares"].extend(pull.get("plates", []))
    return out


def _plate_means(rows: list[dict]) -> dict:
    stats = [row["plates"] for row in rows if "plates" in row]
    if not stats:
        return {}
    return {
        "count_mean": round(float(np.mean([s["n"] for s in stats])), 2),
        "speed_cm_per_yr_mean": round(float(np.mean([s["speed_cm"] for s in stats])), 3),
        "largest_share_mean": round(float(np.mean([s["largest"] for s in stats])), 4),
        "largest_share_max": round(float(np.max([s["largest"] for s in stats])), 4),
    }


def summarize_seed(path: Path) -> dict:
    lines = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    header, rows = lines[0], [r for r in lines[1:] if "myr" in r]
    tail = lines[-1] if "end_craton_m3" in lines[-1] else {}
    fusions = [fusion for row in rows for fusion in row.get("fusions", [])]
    losses, contacts, other = _loss_tally(rows), _contact_tally(rows), _progress_and_pull(rows)
    total = losses["total"]
    pull, plate_shares = other["pull"], other["plate_shares"]
    brake, proposed = contacts["brake"], contacts["proposed"]

    def share(x: float) -> float | None:
        return round(x / total, 4) if total > 0.0 else None

    def pct(values: list[float], q: float) -> float | None:
        return float(np.percentile(values, q)) if values else None

    rates = {k: _weighted_mean(*v) for k, v in contacts["close"].items()}
    return {
        "seed": header["seed"],
        "myr": rows[-1]["myr"] + header["step_myr"] if rows else 0.0,
        "start_craton_km3": header["start_craton_m3"] / 1e9,
        "end_craton_km3": tail.get("end_craton_m3", float("nan")) / 1e9,
        "destroyed_km3": total / 1e9,
        "ledger_sinks_km3": losses["sinks"] / 1e9,
        "share_near_contact": share(losses["near_contact"]),
        "share_near_jammed_contact": share(losses["near_jammed"]),
        "J": share(losses["jammed_unfused"]),
        "jammed_pair_steps": contacts["jammed_pair_steps"],
        "collision_time_drops": other["dropped"],
        "collision_time_drops_jammed": other["dropped_jammed"],
        "collision_time_dropped_jammed_myr": round(other["dropped_jammed_myr"], 1),
        "buoyant_pull_share": round(pull["buoyant"] / pull["total"], 4) if pull["total"] else None,
        "cratonic_pull_share": round(pull["cratonic"] / pull["total"], 4) if pull["total"] else None,
        "plate_steps_buoyant_pull_over_10pct": round(float(np.mean(np.array(plate_shares) > 0.1)), 4) if plate_shares else None,
        "closing_cm_per_yr": {k: (round(v, 3) if v is not None else None) for k, v in rates.items()},
        "jammed_share_of_collision_length": round(contacts["jammed_len"] / contacts["collision_len"], 4) if contacts["collision_len"] else None,
        "brake_n_per_m_at_jams": {"p10": pct(brake, 10), "p50": pct(brake, 50), "p90": pct(brake, 90)},
        "proposed_n_per_m_at_jams": {"p10": pct(proposed, 10), "p50": pct(proposed, 50), "p90": pct(proposed, 90)},
        "plates": _plate_means(rows),
        "fusions": {
            "roll": sum(f["kind"] == "roll" for f in fusions),
            "forced": sum(f["kind"] == "forced" for f in fusions),
            "craton_overlap_km3": round(sum(f["overlap_m3"] for f in fusions) / 1e9, 1),
            "craton_lost_km3": round(sum(f["lost_m3"] for f in fusions) / 1e9, 1),
        },
    }


def summarize(paths: list[Path], json_out: Path | None) -> None:
    seeds = [summarize_seed(path) for path in paths]
    for s in seeds:
        rates = s["closing_cm_per_yr"]
        print(
            f"seed {s['seed']} to {s['myr']:.0f} Myr: destroyed {s['destroyed_km3']:.3g} km^3 (ledger sinks {s['ledger_sinks_km3']:.3g}); "
            f"near contact {s['share_near_contact']}, near jammed {s['share_near_jammed_contact']}, J {s['J']}\n"
            f"  collision time drops {s['collision_time_drops']} ({s['collision_time_drops_jammed']} jammed, "
            f"{s['collision_time_dropped_jammed_myr']} Myr lost); jammed pair-steps {s['jammed_pair_steps']}\n"
            f"  buoyant pull share {s['buoyant_pull_share']} (cratonic {s['cratonic_pull_share']}); "
            f"jammed share of collision length {s['jammed_share_of_collision_length']}\n"
            f"  closing cm/yr, continental pairs: craton {rates['cc_craton']} / clear {rates['cc_clear']}; "
            f"all: craton {rates['all_craton']} / clear {rates['all_clear']}\n"
            f"  brake at jams {s['brake_n_per_m_at_jams']}; proposed {s['proposed_n_per_m_at_jams']}; fusions {s['fusions']}\n"
            f"  plates {s['plates']}"
        )
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(json.dumps({"seeds": seeds}, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--reach-km", type=float, default=500.0)
    parser.add_argument(
        "--brake-reach-km", type=float, help="override torque.CRATON_BRAKE_REACH_KM for this replay (RFC 0008 solution 2's sweep, #423)"
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--summarize", type=Path, nargs="+")
    parser.add_argument("--json", type=Path, help="with --summarize: write the per-seed summary here")
    args = parser.parse_args()
    if args.summarize:
        summarize(args.summarize, args.json)
        return
    if args.out is None:
        parser.error("--out is required unless --summarize is given")
    if args.brake_reach_km is not None:
        torque.CRATON_BRAKE_REACH_KM = args.brake_reach_km
    run(args.seed, args.years, args.step, args.reach_km, args.out)


if __name__ == "__main__":
    main()
