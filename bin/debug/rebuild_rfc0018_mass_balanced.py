"""Rebuild RFC 0018's mass-balanced candidate generator and breakup prior from recorded evidence (#848).

The code that G1/G2 froze is lost. This rebuilds it from its written spec and checks every
candidate against the recorded probe output on the same 20 G0 snapshots.

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/rebuild_rfc0018_mass_balanced.py \
        --snapshots ~/mantle-bloom-captures/rfc0018-g0-snapshots \
        --oracle ~/mantle-bloom-captures/rfc0018-issue803/mass-balanced-probe.json \
        --g0-script ~/mantle-bloom-captures/issue848/measure_rfc0018_g0.py \
        --out analysis/issue848/mass-balanced-rebuild.json [--seeds N N ...]

``--g0-script`` is ``bin/debug/measure_rfc0018_g0.py`` from origin/codex/issue-797-evidence-rfc0018-g0,
loaded by path; its resistance, trace, component and opening helpers are reused unchanged.

Rebuilt rule (each point is checked against the oracle by this script):
  * Fiedler vector of the weighted Laplacian with conductance 1/edge_cost (the G0/G2 mapping).
  * ``count = min(12, edges crossing the two flow clusters with a nonzero Fiedler gap)`` targets,
    ``linspace(0.15, 0.50, count)``, each a share of plate mass for the low-Fiedler side.
  * The low side is the nodes up to the FIRST one where cumulative mass along ascending Fiedler value
    reaches the target; the high side is the other daughter.
  * Every refinement iteration (up to 8) refits the resistance to the current side, re-solves the
    Fiedler vector, and re-applies the same mass-quantile cut. Identical partitions are dropped.
Prior (analysis/issue797/README.md): ``p_i = O_eff / (O_eff + O_ref * R_i / R_ref)`` with
``O_eff = max(O, 0) * S / S_ref``; ``q`` proportional to the nonnegative momentum score at ``mu = 1``
(uniform if all are zero); ``P_break = sum q_i p_i``; structurally invalid candidates have ``p_i = 0``.
"""

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from scipy.cluster.vq import kmeans2  # noqa: E402

from app import mantle, merge_split, persistence  # noqa: E402

TARGET_LO, TARGET_HI, MAX_CANDIDATES, MAX_FIT_ITERATIONS = 0.15, 0.50, 12, 8
O_REF_CM_YR, R_REF, S_REF_RAD_YR = 0.8734616232700081, 61.057432008023106, 6.465517316573429e-09


def load_g0(path: Path):
    spec = importlib.util.spec_from_file_location("g0_candidate_search", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def low_side(fiedler: np.ndarray, mass: np.ndarray, target: float) -> np.ndarray:
    """Nodes up to the first whose cumulative mass along ascending Fiedler value reaches target * total."""
    order = np.argsort(fiedler, kind="stable")
    cumulative = np.cumsum(mass[order])
    k = min(int(np.searchsorted(cumulative, target * cumulative[-1], side="left")), len(cumulative) - 2)
    low = np.zeros(len(fiedler), bool)
    low[order[: k + 1]] = True
    return low


def refine(solve, side, mass, target, cap, shortcut=False):
    """The frozen refinement loop: refit, re-solve and re-apply the mass quantile until unchanged or ``cap`` passes.

    Each pass is a fixed map from one partition to the next (the solver starts from a fixed ``v0``), so once a
    partition recurs the rest of the sequence repeats with the same period. With ``shortcut`` the loop stops at
    the first recurrence and returns the partition the capped loop would have ended on, from the period alone.
    The loop also records the cycle it found and how often the Fiedler vector's arbitrary sign flipped between
    passes, which moves the low-mass side to the other end of the plate."""
    states, first_seen = [side], {side.tobytes(): 0}
    info = {"fit_passes": 0, "converged": False, "cycle_start": None, "cycle_period": None,
            "cycle_detected_pass": None, "cycle_disagreement": None, "fiedler_sign_flips": 0}
    previous = None
    for t in range(cap):
        fiedler = solve(side)
        info["fit_passes"] = t + 1
        if previous is not None and float(fiedler @ previous) < 0:
            info["fiedler_sign_flips"] += 1
        previous = fiedler
        following = ~low_side(fiedler, mass, target)
        start = first_seen.get(following.tobytes())
        if start is not None and info["cycle_period"] is None:
            period = t + 1 - start
            info.update(cycle_start=start, cycle_period=period, cycle_detected_pass=t + 1,
                        cycle_disagreement=float(np.mean(states[start] != states[start + 1])) if period > 1 else 0.0)
        if np.array_equal(following, side):  # the frozen loop's own stop: the partition did not change
            info["converged"] = True
            return side, info
        if start is not None and shortcut:  # periodic from ``start``: the capped loop ends on this state
            return states[start + (cap - start) % (t + 1 - start)], info
        first_seen.setdefault(following.tobytes(), t + 1)
        states.append(following)
        side = following
    return side, info


def propose(g0, world, plate, max_fit_iterations: int = MAX_FIT_ITERATIONS, cycle_shortcut: bool = False) -> list[dict]:
    """The mass-balanced candidates for one plate, each with the same fields as ``g0.choose_cut`` variants.

    ``max_fit_iterations`` caps the refinement passes (the frozen method uses 8; 0 keeps the initial cut).
    ``cycle_shortcut`` stops a candidate's loop at its first repeated partition (see ``refine``). Each candidate
    also reports ``fit_passes`` (solves in the loop), ``converged`` (the loop stopped early) and its cycle."""
    points, _ = plate.all_points_and_elevation()
    flow = mantle.flow_at(points, world.mantle_centers)
    _, labels = kmeans2(flow, k=2, minit="++", seed=plate.plate_id)
    if len(np.unique(labels)) < 2:
        raise ValueError("mantle-flow clustering produced only one regime")
    cluster_omega = np.asarray([mantle.fit_euler_pole(points[labels == k], flow[labels == k]) for k in (0, 1)])
    i, j, lengths = g0._shared_edges(plate)
    node_resistance, areal_mass, _ = g0._node_resistance(world, plate, points, flow, list(cluster_omega))
    median_edge_length = max(float(np.median(lengths)), 1.0)

    def cost_for(side):
        resistance, _, _ = g0._node_resistance(world, plate, points, flow, cluster_omega, side)
        return lengths / median_edge_length * 0.5 * (resistance[i] + resistance[j])

    initial_cost = lengths / median_edge_length * 0.5 * (node_resistance[i] + node_resistance[j])
    mass = areal_mass * plate.node_areas_m2()
    fiedler = g0._fiedler_vector(len(points), i, j, initial_cost)
    eligible = np.flatnonzero((labels[i] != labels[j]) & (np.abs(fiedler[i] - fiedler[j]) > 1e-12))
    count = min(MAX_CANDIDATES, len(eligible))
    if not count:
        raise ValueError("no adjacency edge crosses both the flow clusters and the spectral partition")
    targets = np.linspace(TARGET_LO, TARGET_HI, count) if count > 1 else np.array([TARGET_LO])
    spacing_rad = g0.line_spacing_rad(world.node_density)
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    _, radius = merge_split.geometry.bounding_sphere(points)
    relax = max(min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    pole_threshold = merge_split.SPLIT_MIN_POLE_SEPARATION * (1.0 - relax)
    boundary_nodes = plate.exposed_sides().any(axis=1)
    candidates, seen = [], set()
    for target in targets:
        side, loop = refine(lambda s: g0._fiedler_vector(len(points), i, j, cost_for(s)),
                            ~low_side(fiedler, mass, target), mass, target, max_fit_iterations, cycle_shortcut)
        passes, converged = loop["fit_passes"], loop["converged"]
        if not side.any() or side.all():
            continue
        key = min(side.tobytes(), (~side).tobytes())
        if key in seen:
            continue
        seen.add(key)
        cost = cost_for(side)
        cut = side[i] != side[j]
        opening = g0._opening(points, flow, i, j, lengths, side, mass, plate.omega)
        if opening is None:
            continue
        components = g0._largest_component_share(len(points), i, j, side)
        trace_count, trace_share, main_touches, _ = g0._trace_stats(points, i, j, lengths, cut, boundary_nodes, spacing_rad)
        omega_a = mantle.fit_euler_pole(points[side], flow[side])
        omega_b = mantle.fit_euler_pole(points[~side], flow[~side])
        pole_separation = float(np.linalg.norm(omega_a - omega_b))
        connected = min(components) >= g0.MIN_COMPONENT_SHARE
        trace_ok = trace_share >= 0.90 and main_touches >= 2
        size_ok = int(side.sum()) >= min_nodes and int((~side).sum()) >= min_nodes
        cut_nodes = np.unique(np.concatenate([i[cut], j[cut]]))
        candidates.append({
            **loop, "target": float(target), "fit_passes": passes, "converged": converged, "daughter_node_counts": [int(side.sum()), int((~side).sum())],
            "cut_cost": float(cost[cut].sum()), "admissible_connected": bool(connected),
            "single_edge_to_edge_trace": bool(trace_ok), "daughter_size_ok": bool(size_ok),
            "structurally_valid": bool(connected and trace_ok and size_ok),
            "pole_separation_rad_per_yr": pole_separation, "pole_separation_ok": pole_separation >= pole_threshold,
            "boundary_touch_nodes": int(boundary_nodes[cut_nodes].sum()), "cut_trace_components": int(trace_count),
            "opening": opening, "side_mask": side,
        })
    return candidates


def breakup_prior(candidates: list[dict], mu_key: str = "1.0") -> dict:
    """Frozen O/R/S prior and proposal weights (analysis/issue797/README.md); fills q, p per candidate."""
    scores = np.array([max(float(c["opening"]["momentum_scores_by_mu"][mu_key]), 0.0) for c in candidates])
    q = scores / scores.sum() if scores.sum() > 0 else np.full(len(candidates), 1.0 / max(len(candidates), 1))
    for c, weight in zip(candidates, q):
        o_eff = max(float(c["opening"]["net_opening_cm_per_yr"]), 0.0) * c["pole_separation_rad_per_yr"] / S_REF_RAD_YR
        denom = o_eff + O_REF_CM_YR * c["cut_cost"] / R_REF
        c["q"] = float(weight)
        c["p"] = float(o_eff / denom) if c["structurally_valid"] and denom > 0 else 0.0
    return {"P_break": float(sum(c["q"] * c["p"] for c in candidates)), "candidates": len(candidates)}


def compare(candidates: list[dict], recorded: list[dict]) -> dict:
    """Per-candidate agreement with the recorded probe, keyed by position (the order is the target order)."""
    fields = {"partition": 0, "structural": 0, "momentum": 0, "net_opening": 0}
    for mine, ref in zip(candidates, recorded):
        fields["partition"] += int(mine["daughter_node_counts"] == ref["daughter_node_counts"]
                                   and abs(mine["cut_cost"] - ref["cut_cost"]) < 1e-6 * max(1.0, abs(ref["cut_cost"])))
        ref_valid = bool(ref["admissible_connected"] and ref["single_edge_to_edge_trace"] and ref["daughter_size_ok"])
        fields["structural"] += int(mine["structurally_valid"] == ref_valid)
        a, b = mine["opening"], ref["opening"]
        fields["momentum"] += int(np.isclose(a["momentum_scores_by_mu"]["1.0"], b["momentum_scores_by_mu"]["1.0"], rtol=1e-6))
        fields["net_opening"] += int(np.isclose(a["net_opening_cm_per_yr"], b["net_opening_cm_per_yr"], rtol=1e-6, atol=1e-12))
    return fields


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--snapshots", required=True, type=Path)
    parser.add_argument("--oracle", required=True, type=Path)
    parser.add_argument("--g0-script", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--seeds", type=int, nargs="*")
    args = parser.parse_args()
    g0 = load_g0(args.g0_script.expanduser())
    oracle = json.loads(args.oracle.expanduser().read_text())
    started = time.perf_counter()
    totals = {"opportunities": 0, "candidates": 0, "recorded_candidates": 0, "valid": 0, "recorded_valid": 0,
              "opps_with_valid": 0, "opps_positive_p_break": 0, "agreement": {"partition": 0, "structural": 0, "momentum": 0, "net_opening": 0}}
    p_breaks, rows = [], []
    for world_row in oracle["worlds"]:
        if args.seeds and int(world_row["seed"]) not in args.seeds:
            continue
        with open(args.snapshots.expanduser() / f"g0-snapshot-seed-{world_row['seed']}.bloom", "rb") as handle:
            world = persistence.load_world_file(handle)
        by_id = {int(p.plate_id): p for p in world.plates}
        for opp in world_row["opportunities"]:
            candidates = propose(g0, world, by_id[int(opp["plate_id"])])
            prior = breakup_prior(candidates)
            agreement = compare(candidates, opp["candidates"])
            valid = sum(c["structurally_valid"] for c in candidates)
            recorded_valid = sum(bool(c["admissible_connected"] and c["single_edge_to_edge_trace"] and c["daughter_size_ok"])
                                 for c in opp["candidates"])
            totals["opportunities"] += 1
            totals["candidates"] += len(candidates)
            totals["recorded_candidates"] += len(opp["candidates"])
            totals["valid"] += valid
            totals["recorded_valid"] += recorded_valid
            totals["opps_with_valid"] += int(valid > 0)
            totals["opps_positive_p_break"] += int(prior["P_break"] > 0)
            for key, value in agreement.items():
                totals["agreement"][key] += value
            p_breaks.append(prior["P_break"])
            rows.append({"seed": int(world_row["seed"]), "plate_id": int(opp["plate_id"]), "candidates": len(candidates),
                         "recorded": len(opp["candidates"]), "valid": valid, "P_break": prior["P_break"], **agreement})
            print(json.dumps(rows[-1]), flush=True)
    totals["mean_P_break"] = float(np.mean(p_breaks)) if p_breaks else None
    totals["runtime_s"] = round(time.perf_counter() - started, 1)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({"totals": totals, "opportunities": rows}, indent=1) + "\n")
    print(json.dumps(totals, indent=1))


if __name__ == "__main__":
    main()
