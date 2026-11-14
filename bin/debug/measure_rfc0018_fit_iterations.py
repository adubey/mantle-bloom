"""How much do refinement passes 2 to 8 change RFC 0018's candidates? (#848)

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/measure_rfc0018_fit_iterations.py \
        --g0-script <measure_rfc0018_g0.py> --out analysis/issue848/fit-iterations.json \
        [--snapshots <rfc0018-g0-snapshots> --seeds N ...] [--saves <a.mbworld> <b.mbworld> ...] \
        [--caps 0 1 2 3 4 cycle]

The frozen candidate method refits resistance and re-solves the Fiedler vector up to 8 times per candidate
(one sparse solve per pass), so a decision costs up to 1 + 12 x 8 = 97 solves. This reruns the rebuilt
generator (``rebuild_rfc0018_mass_balanced.propose``) on every split opportunity with the pass cap lowered
and compares each cap against the 8-pass reference: share of identical partitions, node disagreement, the
structural-validity flag, cut cost, and each opportunity's ``P_break``. Candidates are aligned by their
mass target because a lower cap can change which candidates survive the duplicate drop. The ``cycle`` arm
keeps the cap of 8 but stops each loop at its first repeated partition (``rebuild_rfc0018_mass_balanced.refine``);
if the solve is deterministic it must reproduce the reference exactly, so any disagreement there is a finding.

``--snapshots`` takes the 20 G0 snapshots (opportunities as at generation, 2 Myr); ``--saves`` takes
checkpoints of mature worlds, measured as they stand (no stepping).
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "backend")]

import rebuild_rfc0018_mass_balanced as rebuild  # noqa: E402

persistence = rebuild.persistence

REFERENCE = rebuild.MAX_FIT_ITERATIONS


CYCLE_ARM = "cycle"  # the full cap, but each loop stops at its first repeated partition


def run_cap(g0, world, plate, cap, shortcut=False):
    started = time.perf_counter()
    candidates = rebuild.propose(g0, world, plate, max_fit_iterations=cap, cycle_shortcut=shortcut)
    prior = rebuild.breakup_prior(candidates)
    return candidates, prior["P_break"], time.perf_counter() - started


def compare(reference, other):
    """Per-candidate agreement between a capped run and the reference, aligned by mass target."""
    by_target = {round(c["target"], 9): c for c in other}
    rows = []
    for ref in reference:
        match = by_target.get(round(ref["target"], 9))
        if match is None:
            rows.append({"present": False})
            continue
        a, b = ref["side_mask"], match["side_mask"]
        disagree = float(min(np.mean(a != b), np.mean(a == b)))
        rows.append({"present": True, "identical": bool(disagree == 0.0), "node_disagreement": disagree,
                     "valid_equal": ref["structurally_valid"] == match["structurally_valid"],
                     "cost_rel_diff": abs(ref["cut_cost"] - match["cut_cost"]) / max(abs(ref["cut_cost"]), 1e-12)})
    return rows


def measure_world(g0, world, caps, label, rows):
    for plate in world.plates:
        if not g0._is_opportunity(world, plate):
            continue
        try:
            reference, p_ref, t_ref = run_cap(g0, world, plate, REFERENCE)
        except ValueError:
            continue
        row = {"world": label, "plate_id": int(plate.plate_id), "nodes": int(plate.node_count()),
               "reference": {"p_break": p_ref, "seconds": t_ref, "valid": sum(c["structurally_valid"] for c in reference),
                             "candidates": len(reference), "converged": sum(c["converged"] for c in reference),
                             "passes": [c["fit_passes"] for c in reference],
                             "cycles": [{k: c[k] for k in CYCLE_FIELDS} for c in reference]},
               "caps": {}}
        for arm in caps:
            cap, shortcut = (REFERENCE, True) if arm == CYCLE_ARM else (int(arm), False)
            capped, p_cap, t_cap = run_cap(g0, world, plate, cap, shortcut)
            row["caps"][arm] = {"p_break": p_cap, "seconds": t_cap, "valid": sum(c["structurally_valid"] for c in capped),
                                     "candidates": len(capped), "compare": compare(reference, capped),
                                     "passes": sum(c["fit_passes"] for c in capped)}
        rows.append(row)
        print(json.dumps({"world": label, "plate": row["plate_id"], "nodes": row["nodes"]}), flush=True)


CYCLE_FIELDS = ["cycle_start", "cycle_period", "cycle_detected_pass", "cycle_disagreement", "fiedler_sign_flips"]


def histogram(values) -> dict:
    keys, counts = np.unique([str(v) for v in values], return_counts=True)
    return {k: int(n) for k, n in zip(keys, counts)}


def cycle_summary(cycles) -> dict:
    """How the reference loops end: the period and start of the first repeat, and Fiedler sign flips."""
    looping = [c for c in cycles if c["cycle_period"] and c["cycle_period"] > 1]
    return {"candidates": len(cycles), "period_histogram": histogram([c["cycle_period"] for c in cycles]),
            "cycle_start_histogram": histogram([c["cycle_start"] for c in cycles]),
            "detected_pass_histogram": histogram([c["cycle_detected_pass"] for c in cycles]),
            "sign_flips_histogram": histogram([c["fiedler_sign_flips"] for c in cycles]),
            "cycle_disagreement_mean": float(np.mean([c["cycle_disagreement"] for c in looping])) if looping else None,
            "cycle_disagreement_median": float(np.median([c["cycle_disagreement"] for c in looping])) if looping else None}


def summarize(rows, caps):
    out = {"opportunities": len(rows), "reference_cap": REFERENCE}
    cands = [c for r in rows for c in range(r["reference"]["candidates"])]
    out["reference"] = {"candidates": len(cands),
                        "converged_share": float(sum(r["reference"]["converged"] for r in rows) / len(cands)),
                        "passes_histogram": {str(k): int(v) for k, v in zip(*np.unique(
                            [p for r in rows for p in r["reference"]["passes"]], return_counts=True))},
                        "seconds": float(sum(r["reference"]["seconds"] for r in rows)),
                        "p_break_total": float(sum(r["reference"]["p_break"] for r in rows)),
                        "valid_total": int(sum(r["reference"]["valid"] for r in rows))}
    out["caps"] = {}
    out["reference"]["cycles"] = cycle_summary([c for r in rows for c in r["reference"]["cycles"]])
    for arm in caps:
        comp = [c for r in rows for c in r["caps"][arm]["compare"]]
        present = [c for c in comp if c["present"]]
        p_ref = np.array([r["reference"]["p_break"] for r in rows])
        p_cap = np.array([r["caps"][arm]["p_break"] for r in rows])
        v_ref = np.array([r["reference"]["valid"] for r in rows])
        v_cap = np.array([r["caps"][arm]["valid"] for r in rows])
        disagreement = np.array([c["node_disagreement"] for c in present])
        out["caps"][arm] = {
            "candidates_missing": len(comp) - len(present),
            "identical_partition_share": float(np.mean([c["identical"] for c in present])),
            "node_disagreement_mean": float(disagreement.mean()), "node_disagreement_p95": float(np.percentile(disagreement, 95)),
            "node_disagreement_max": float(disagreement.max()),
            "valid_flag_agreement": float(np.mean([c["valid_equal"] for c in present])),
            "valid_count_equal_opportunities": float(np.mean(v_ref == v_cap)),
            "p_break_total": float(p_cap.sum()), "p_break_total_rel_error": float(abs(p_cap.sum() - p_ref.sum()) / max(p_ref.sum(), 1e-12)),
            "p_break_max_abs_diff": float(np.abs(p_cap - p_ref).max()),
            "positive_p_break_agreement": float(np.mean((p_ref > 0) == (p_cap > 0))),
            "seconds": float(sum(r["caps"][arm]["seconds"] for r in rows)),
            "solves_per_opportunity": float(np.mean([r["caps"][arm]["passes"] + 1 for r in rows])),
        }
    out["reference"]["solves_per_opportunity"] = float(np.mean([sum(r["reference"]["passes"]) + 1 for r in rows]))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--g0-script", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--snapshots", type=Path)
    parser.add_argument("--seeds", type=int, nargs="*")
    parser.add_argument("--saves", type=Path, nargs="*", default=[])
    parser.add_argument("--caps", nargs="*", default=["0", "1", "2", "3", "4"],
                        help=f"pass caps to compare with the reference, or '{CYCLE_ARM}' for the repeat shortcut")
    args = parser.parse_args()
    g0 = rebuild.load_g0(args.g0_script.expanduser())
    rows: list[dict] = []
    paths = [(f"snapshot-{s}", args.snapshots.expanduser() / f"g0-snapshot-seed-{s}.bloom") for s in (args.seeds or [])]
    paths += [(f"{p.parent.name}/{p.name}", p.expanduser()) for p in args.saves]
    for label, path in paths:
        with open(path, "rb") as handle:
            world = persistence.load_world_file(handle)
        measure_world(g0, world, args.caps, label, rows)
    summary = summarize(rows, args.caps)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
