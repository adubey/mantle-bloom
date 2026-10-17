#!/usr/bin/env python3
"""Which force-proxy rules bind (RFC 0012 gate G1, issue #465). No model change.

Replays one seed from generation (`world.generate_world(seed)`, default controls) and taps the
rules RFC 0012's Problem table lists, in this process only; nothing is fed back:

- **Rate clamp.** Every `mantle.clamp_rate` call, by caller (`torque.integrate_omega`,
  `merge_split.maybe_split_plate`, `quad_merge.merge`): the share clipped at the 0.5 cm/yr floor and at the
  15 cm/yr cap.
- **Motion solve.** Per `integrate_omega` call: the slowest relaxation time of
  `(I + g K) w = I w_old + g (tau + b)`, the largest eigenvalue of `K^-1 I` (in seconds, to
  set beside the RFC's illustrative 7e-7 s), the step over that time, and how far the step's
  omega is from the quasi-static `K^-1 (tau + b)` before the clamp.
- **Merges.** Each fusion by path: the contact timer and probability roll
  (`update_collision_progress`) or the deep-overlap backstop (`pop_ready_forced_merge`), with
  the years accumulated against the pair's timer, and the contact's normal closing rate and
  |w_a - w_b| R at the moment of merging. Steps on which a pair past its timer failed the roll.
- **Splits.** Each `maybe_split_plate` call by the gate that decided it (size, age, flow-fit
  residual, later geometric gates, the 0.55 roll, split), and for those past the residual
  gate, whether plate size or the accumulated stress proxy set the relaxation.
- **Convergent yield.** Every `rheology.convergent_strain` call: convergent nodes, the share
  past yield, and the share within 5% of the rate law's saturation; `shortening.ring_cascade`
  demand absorbed, lost to decay, and unrouted. Each 10 Myr: continental nodes at 99% of the
  Hc or Hm cap.
- **Order dependence** (`--order-myr` checkpoints): on copies of the world, each plate's
  omega after shifting every plate in `--order-shuffles` shuffled orders, against the model's
  own order.
- **Band rows** at the same checkpoints: subducting and divergent band nodes per distinct
  nearest neighbour node (one per spacing of boundary), the factor
  `probe_boundary_load_normalization.py` finds slab pull and ridge push scaled by.

`--profile` instead times one density's step phases (`--profile-steps` steps after one warm-up),
for gate G1's cost row.

Usage:
    backend/.venv/bin/python bin/debug/audit_force_rules.py --seed 326494039 --myr 500 \\
        --out analysis/rfc0012/rules_326494039.json
    backend/.venv/bin/python bin/debug/audit_force_rules.py --seed 326494039 --profile \\
        --density 4 --out analysis/rfc0012/profile_d4.json
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import resource
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

SOURCE_ROOT = Path(os.environ.get("MB_SOURCE_ROOT", Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(SOURCE_ROOT / "backend"))

from app import lithosphere, mantle, merge_split, plates, rheology, shortening, torque  # noqa: E402
from app import stats  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402

R = lithosphere.PLANET_RADIUS_M
CM_PER_RAD_YR = R * 100.0  # rad/yr -> cm/yr at the surface
SATURATION = (1.0 - np.tan(rheology.INTERNAL_FRICTION_ANGLE_RAD)) / np.tan(rheology.INTERNAL_FRICTION_ANGLE_RAD)


class Tap:
    def __init__(self, world):
        self.world = world
        self.clamp = defaultdict(Counter)
        self.solve = []  # per integrate_omega call
        self.merges = []
        self.failed_rolls = 0
        self.split_gates = Counter()
        self.split_relax = Counter()
        self.yield_counts = Counter()
        self.cascade = Counter()
        self.caps = []
        self.order = []
        self.band_rows = []
        self._merge_path = None

    # -- rate clamp -----------------------------------------------------------------------
    def install(self) -> None:
        tap = self
        orig_clamp = mantle.clamp_rate

        def clamp_rate(omega, min_rate=mantle.MIN_PLATE_RATE, max_rate=mantle.MAX_PLATE_RATE):
            frame = sys._getframe(1)
            caller = f"{Path(frame.f_code.co_filename).stem}.{frame.f_code.co_name}"
            rate = float(np.linalg.norm(omega))
            c = tap.clamp[caller]
            c["calls"] += 1
            c["floor"] += rate < min_rate
            c["cap"] += rate > max_rate
            return orig_clamp(omega, min_rate, max_rate)

        mantle.clamp_rate = clamp_rate

        orig_integrate = torque.integrate_omega

        def integrate_omega(plate, explicit_torque, drag_b, drag_k, inertia_tensor, years):
            g = years * torque.SECONDS_PER_YEAR**2
            lhs = inertia_tensor + g * drag_k
            raw = np.linalg.lstsq(lhs, inertia_tensor @ plate.omega + g * (explicit_torque + drag_b), rcond=None)[0]
            row = {"rate_cm_yr": float(np.linalg.norm(raw)) * CM_PER_RAD_YR}
            try:
                relax = np.linalg.eigvals(np.linalg.solve(drag_k, inertia_tensor)).real
                qs = np.linalg.solve(drag_k, explicit_torque + drag_b)
                # K^-1 I is in yr^2 * (s/yr)^-2 = years / SECONDS_PER_YEAR^2 -- see lhs above.
                t_yr = float(relax.max()) / torque.SECONDS_PER_YEAR**2
                row["relax_s"] = t_yr * torque.SECONDS_PER_YEAR
                row["step_over_relax"] = years / t_yr if t_yr > 0 else None
                row["qs_rel_diff"] = float(np.linalg.norm(raw - qs) / max(np.linalg.norm(qs), 1e-30))
            except np.linalg.LinAlgError:
                row["singular_k"] = True
            tap.solve.append(row)
            return orig_integrate(plate, explicit_torque, drag_b, drag_k, inertia_tensor, years)

        torque.integrate_omega = integrate_omega

        # -- merges ---------------------------------------------------------------------
        orig_progress = merge_split.update_collision_progress

        def update_collision_progress(world, years):
            ready = orig_progress(world, years)
            past = [p for p, acc in world.collision_progress.items() if acc >= merge_split._collision_threshold_years(world.seed, p)]
            tap.failed_rolls += len([p for p in past if p not in ready])
            tap._merge_path = ("timer_roll", ready[0]) if ready else None
            return ready

        merge_split.update_collision_progress = update_collision_progress

        orig_pop = merge_split.pop_ready_forced_merge

        def pop_ready_forced_merge(world):
            acc = dict(world.overlap_progress)
            pair = orig_pop(world)
            if pair is not None:
                key = tuple(sorted(pair))
                tap._merge_path = ("forced_overlap", pair, acc.get(key, 0.0))
            return pair

        merge_split.pop_ready_forced_merge = pop_ready_forced_merge

        orig_merge = merge_split._merge_plates_with_craton_audit

        def merge(world, id_keep, id_absorb):
            path = tap._merge_path[0] if tap._merge_path else "other"
            row = {"myr": world.elapsed_years / 1e6, "path": path, **contact_rate(world, id_keep, id_absorb)}
            pair = (min(id_keep, id_absorb), max(id_keep, id_absorb))
            if path == "timer_roll":
                row["threshold_myr"] = merge_split._collision_threshold_years(world.seed, pair) / 1e6
            elif path == "forced_overlap":
                row["overlap_myr"] = tap._merge_path[2] / 1e6
            tap.merges.append(row)
            tap._merge_path = None
            return orig_merge(world, id_keep, id_absorb)

        merge_split._merge_plates_with_craton_audit = merge

        # -- splits ---------------------------------------------------------------------
        orig_split = merge_split.maybe_split_plate

        def maybe_split_plate(world, plate):
            gate, relax_by = split_gate(world, plate)
            # The log is trimmed to MAX_EVENT_LOG_LENGTH, so compare its last entry, not its length.
            last_before = world.events[-1] if world.events else None
            result = orig_split(world, plate)
            if gate == "past_residual":
                tap.split_relax[relax_by] += 1
                if result is not None:
                    gate = "split"
                elif world.events and world.events[-1] is not last_before and "rift failed" in world.events[-1][1]:
                    gate = "failed_roll"
                else:
                    gate = "geometric_gate"
            tap.split_gates[gate] += 1
            return result

        merge_split.maybe_split_plate = maybe_split_plate

        # -- convergent yield -------------------------------------------------------------
        orig_strain = rheology.convergent_strain

        def convergent_strain(closing_rate_m_per_s, years_myr, fault_factor, strength=1.0):
            closing = np.asarray(closing_rate_m_per_s)
            convergent = closing > 0.0
            rate = np.clip(rheology.plastic_strain_rate_per_myr(closing), 0.0, None)
            norm = rate / rheology.PLASTIC_THICKENING_RATE_PER_MYR_PER_YIELD_EXCESS
            tap.yield_counts["convergent"] += int(np.count_nonzero(convergent))
            tap.yield_counts["past_yield"] += int(np.count_nonzero(rate > 0.0))
            tap.yield_counts["near_saturation"] += int(np.count_nonzero(norm >= 0.95 * SATURATION))
            return orig_strain(closing_rate_m_per_s, years_myr, fault_factor, strength)

        rheology.convergent_strain = convergent_strain

        orig_cascade = shortening.ring_cascade

        def ring_cascade(problem):
            result = orig_cascade(problem)
            tap.cascade["demand_m2"] += float(np.sum(problem.demand_m2))
            tap.cascade["absorbed_m2"] += float(np.sum(result.absorbed_m2))
            tap.cascade["decay_m2"] += float(np.sum(result.returned_decay_m2))
            tap.cascade["unrouted_m2"] += float(np.sum(result.returned_unrouted_m2))
            return result

        shortening.ring_cascade = ring_cascade

    def record_caps(self) -> None:
        hc_cap, hm_cap = lithosphere.MAX_CRUSTAL_THICKNESS_M, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
        n = hc_n = hm_n = 0
        for plate in self.world.plates:
            if plate.node_count() == 0:
                continue
            cont = plates.node_is_continental(plate)
            n += int(np.count_nonzero(cont))
            hc_n += int(np.count_nonzero(cont & (plate.collect("crustal_thickness_m") >= 0.99 * hc_cap)))
            hm_n += int(np.count_nonzero(cont & (plate.collect("mantle_lithosphere_thickness_m") >= 0.99 * hm_cap)))
        self.caps.append({"myr": self.world.elapsed_years / 1e6, "continental_nodes": n, "hc_at_cap": hc_n, "hm_at_cap": hm_n})


def contact_rate(world, id_a: int, id_b: int) -> dict:
    """The a-b contact's mean normal closing rate over a's band nodes facing b, and |w_a - w_b| R."""
    by_id = {p.plate_id: p for p in world.plates}
    a, b = by_id.get(id_a), by_id.get(id_b)
    if a is None or b is None or a.node_count() == 0 or b.node_count() == 0:
        return {}
    spacing = line_spacing_rad(world.node_density)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    inputs = torque.gather_boundary_force_inputs(a, [b], spacing, reach)
    band = inputs.dist_to_neighbor <= reach
    row = {"relative_omega_cm_yr": float(np.linalg.norm(a.omega - b.omega)) * CM_PER_RAD_YR, "band_nodes": int(band.sum())}
    if np.any(band):
        closing = torque.closing_rad_per_yr(a, inputs)[band] * CM_PER_RAD_YR
        row["closing_mean_cm_yr"] = float(np.mean(closing))
        row["closing_abs_median_cm_yr"] = float(np.median(np.abs(closing)))
    return row


def split_gate(world, plate) -> tuple[str, str | None]:
    """`maybe_split_plate`'s first gates, read without changing anything: which one rejects
    this plate, or `past_residual`, with which relaxation signal was larger."""
    if plate.node_count() < max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density)):
        return "size", None
    if plate.age_years < merge_split.SPLIT_MIN_AGE_YEARS:
        return "age", None
    points, _ = plate.all_points_and_elevation()
    velocities = mantle.flow_at(points, world.mantle_centers)
    _, radius = merge_split.geometry.bounding_sphere(points)
    size_frac = min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD)
    stress_frac = min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT)
    relax = max(size_frac, stress_frac)
    if merge_split._fit_residual_rms(points, velocities, plate.omega) < merge_split.SPLIT_RMS_RESIDUAL_THRESHOLD * (1.0 - relax):
        return "residual", None
    by = "certain" if relax >= 1.0 else ("size" if size_frac >= stress_frac else "stress")
    return "past_residual", by


def order_dependence(world, years: float, shuffles: int, rng) -> dict:
    """Every plate's omega after one shift pass in the model's order and in shuffled orders,
    each from its own copy of the world. Returns the median and max over plates of the
    largest relative omega change any shuffle caused."""

    def shifted(order_ids):
        w = copy.deepcopy(world)
        by_id = {p.plate_id: p for p in w.plates}
        w.plates = [by_id[i] for i in order_ids]
        for plate in w.plates:
            plate.shift(w, years)
        return {p.plate_id: np.array(p.omega) for p in w.plates}

    ids = [p.plate_id for p in world.plates]
    base = shifted(ids)
    worst = dict.fromkeys(ids, 0.0)
    for _ in range(shuffles):
        other = shifted(list(rng.permutation(ids)))
        for i in ids:
            worst[i] = max(worst[i], float(np.linalg.norm(other[i] - base[i]) / max(np.linalg.norm(base[i]), 1e-30)))
    values = np.array(list(worst.values()))
    return {"myr": world.elapsed_years / 1e6, "plates": len(ids), "median_rel": float(np.median(values)), "max_rel": float(values.max())}


def band_rows(world) -> dict:
    """Subducting and divergent band nodes per distinct nearest neighbour node, over all plates."""
    spacing = line_spacing_rad(world.node_density)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    counts = Counter()
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        neighbours = plate.get_neighbours([p for p in world.plates if p.plate_id != plate.plate_id])
        if not neighbours:
            continue
        inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach)
        _c, divergent, _t, _x = torque.classify_boundary_nodes(plate, neighbours, inputs, reach)
        for name, mask in (("subducting", torque.subducting_boundary_mask(plate, inputs, reach)), ("divergent", divergent)):
            counts[f"{name}_nodes"] += int(mask.sum())
            fronts = set(zip(inputs.neighbor_plate_id[mask].tolist(), inputs.neighbor_node_index[mask].tolist()))
            counts[f"{name}_front_nodes"] += len(fronts)
    return {
        "myr": world.elapsed_years / 1e6,
        **counts,
        **{f"{k}_rows": counts[f"{k}_nodes"] / max(1, counts[f"{k}_front_nodes"]) for k in ("subducting", "divergent")},
    }


def summarize(tap: Tap) -> dict:
    solve = tap.solve
    relax = np.array([r["relax_s"] for r in solve if "relax_s" in r])
    ratio = np.array([r["step_over_relax"] for r in solve if r.get("step_over_relax")])
    qs = np.array([r["qs_rel_diff"] for r in solve if "qs_rel_diff" in r])
    rates = np.array([r["rate_cm_yr"] for r in solve])
    merges_by = Counter(m["path"] for m in tap.merges)
    yc = tap.yield_counts
    last_cap = tap.caps[-1] if tap.caps else {}
    continental_hc, continental_hm = [], []
    for plate in tap.world.plates:
        continental = plates.node_is_continental(plate)
        continental_hc.extend(plate.collect("crustal_thickness_m")[continental])
        continental_hm.extend(plate.collect("mantle_lithosphere_thickness_m")[continental])

    def percentiles(values):
        if not values:
            return None
        return {str(p): float(np.percentile(values, p)) for p in (10, 50, 90, 99)}

    return {
        "clamp": {
            caller: {**c, "floor_share": c["floor"] / c["calls"], "cap_share": c["cap"] / c["calls"]} for caller, c in tap.clamp.items()
        },
        "solve": {
            "calls": len(solve),
            "singular_k": sum(1 for r in solve if r.get("singular_k")),
            "relax_s_median": float(np.median(relax)) if len(relax) else None,
            "relax_s_max": float(relax.max()) if len(relax) else None,
            "step_over_relax_min": float(ratio.min()) if len(ratio) else None,
            "qs_rel_diff_max": float(qs.max()) if len(qs) else None,
            "rate_cm_yr_pctl_5_50_95_99": np.percentile(rates, [5, 50, 95, 99]).tolist() if len(rates) else None,
        },
        "merges": {"by_path": dict(merges_by), "failed_roll_pair_steps": tap.failed_rolls, "events": tap.merges},
        "splits": {"by_gate": dict(tap.split_gates), "relaxed_by": dict(tap.split_relax)},
        "yield": {
            **yc,
            "past_yield_share": yc["past_yield"] / max(1, yc["convergent"]),
            "near_saturation_share": yc["near_saturation"] / max(1, yc["convergent"]),
            "cascade": dict(tap.cascade),
        },
        "caps": {
            "series": tap.caps,
            "hc_at_cap_share_mean": float(np.mean([c["hc_at_cap"] / max(1, c["continental_nodes"]) for c in tap.caps])) if tap.caps else None,
            "hm_at_cap_share_mean": float(np.mean([c["hm_at_cap"] / max(1, c["continental_nodes"]) for c in tap.caps])) if tap.caps else None,
            "last": last_cap,
        },
        "order": tap.order,
        "band_rows": tap.band_rows,
        "final_state": {
            "land_fraction": float(stats.compute_stats(tap.world)["land_fraction"]),
            "continental_hc_m_p10_p50_p90_p99": percentiles(continental_hc),
            "continental_hm_m_p10_p50_p90_p99": percentiles(continental_hm),
            "continental_nodes": len(continental_hc),
        },
    }


def run_audit(args) -> dict:
    world = world_mod.generate_world(args.seed, node_density=args.density)
    tap = Tap(world)
    tap.install()
    step = args.step_myr * 1e6
    order_at = {round(m * 1e6) for m in args.order_myr}
    rng = np.random.default_rng(args.seed)
    t0 = time.perf_counter()
    tap.record_caps()
    while world.elapsed_years < args.myr * 1e6 - 1.0:
        world_mod.step_world(world, step)
        now = round(world.elapsed_years)
        if now % round(10e6) == 0:
            tap.record_caps()
            print(f"  {now / 1e6:.0f} Myr  plates={len(world.plates)}  {time.perf_counter() - t0:.0f} s", flush=True)
        if now in order_at:
            tap.band_rows.append(band_rows(world))
            tap.order.append(order_dependence(world, step, args.order_shuffles, rng))
    return {"seconds": time.perf_counter() - t0, **summarize(tap)}


def run_profile(args) -> dict:
    """Wall time of each step phase over `--profile-steps` steps, after one warm-up step."""
    world = world_mod.generate_world(args.seed, node_density=args.density)
    totals = Counter()
    targets = [
        (torque, "gather_boundary_force_inputs", "gather"),
        (torque, "classify_boundary_nodes", "classify"),
        (torque, "integrate_omega", "solve"),
        (torque, "apply_omega_and_rotate", "rotate"),
        (torque, "slab_pull_torque", "forces"),
        (torque, "ridge_push_torque", "forces"),
        (torque, "basal_drag_coefficients", "forces"),
        (torque, "collision_drag_coefficients", "forces"),
        (torque, "contact_jam_factor", "forces"),
        (lithosphere, "moment_of_inertia_tensor", "forces"),
        (merge_split, "apply_topology_changes", "topology"),
    ]
    active = {"on": False}

    def timed(fn, label):
        def wrapper(*a, **kw):
            if not active["on"]:
                return fn(*a, **kw)
            t = time.perf_counter()
            try:
                return fn(*a, **kw)
            finally:
                totals[label] += time.perf_counter() - t

        return wrapper

    for module, name, label in targets:
        if hasattr(module, name):
            setattr(module, name, timed(getattr(module, name), label))
    plate_cls = type(world.plates[0])
    plate_cls.shift = timed(plate_cls.shift, "shift_total")
    plate_cls.deform = timed(plate_cls.deform, "deform_total")

    step = args.step_myr * 1e6
    world_mod.step_world(world, step)  # warm-up: numba compiles, caches fill
    steps = []
    active["on"] = True
    for _ in range(args.profile_steps):
        t = time.perf_counter()
        world_mod.step_world(world, step)
        steps.append(time.perf_counter() - t)
    nodes = sum(p.node_count() for p in world.plates)
    return {
        "density": args.density,
        "nodes": nodes,
        "plates": len(world.plates),
        "step_s": steps,
        "median_step_s": float(np.median(steps)),
        "phase_s_per_step": {k: v / args.profile_steps for k, v in totals.items()},
        "python": sys.version.split()[0],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--myr", type=float, default=500.0)
    ap.add_argument("--step-myr", type=float, default=2.0)
    ap.add_argument("--density", type=float, default=1.0)
    ap.add_argument("--order-myr", type=float, nargs="*", default=[100.0, 300.0, 500.0])
    ap.add_argument("--order-shuffles", type=int, default=3)
    ap.add_argument("--profile", action="store_true")
    ap.add_argument("--profile-steps", type=int, default=5)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=SOURCE_ROOT, capture_output=True, text=True).stdout.strip()
    body = run_profile(args) if args.profile else run_audit(args)
    result = {"issue": 468, "source_issue": 465, "commit": commit, "args": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}, **body}
    result["max_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
