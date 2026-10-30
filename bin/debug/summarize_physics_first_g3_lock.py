"""Summarize consecutive merge-contact lock histories from a G3 trace JSONL."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


EPSILONS_CM_PER_YR = (0.25, 0.5, 1.0, 2.0, 5.0, 10.0)
TAUS_MYR = (1.0, 10.0, 25.0, 50.0, 75.0, 100.0)


def summarize(path: Path) -> dict:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    steps = [row for row in rows if row.get("kind") == "step"]
    summary = next((row for row in rows if row.get("kind") == "summary"), {})
    if not steps:
        raise ValueError("trace contains no step rows")
    step_myr = float(summary.get("step_myr") or (steps[1]["elapsed_myr"] - steps[0]["elapsed_myr"]))
    active = {eps: {} for eps in EPSILONS_CM_PER_YR}
    episodes = {eps: [] for eps in EPSILONS_CM_PER_YR}
    merge_qualifications = {(eps, tau): 0 for eps in EPSILONS_CM_PER_YR for tau in TAUS_MYR}
    merge_total = 0

    for row in steps:
        elapsed = float(row["elapsed_myr"])
        present = {}
        for detail in row.get("collision_search", {}).get("pair_details", []):
            if int(detail.get("close_contact_points", 0)) < 4:
                continue
            rates = detail.get("closing_abs_cm_per_yr")
            if rates is not None:
                present[tuple(sorted(detail["pair"]))] = float(rates["max"])

        for eps, tracks in active.items():
            for pair, (started, _last_seen) in list(tracks.items()):
                if pair not in present or present[pair] >= eps:
                    episodes[eps].append(elapsed - started)
                    del tracks[pair]
            for pair, rate in present.items():
                if rate < eps:
                    started, _last_seen = tracks.get(pair, (elapsed, elapsed))
                    tracks[pair] = (started, elapsed)

        if row.get("actual_merge") is not None:
            merge_total += 1
            pair = tuple(sorted((row["actual_merge"]["keep_plate_id"], row["actual_merge"]["absorbed_plate_id"])))
            for eps, tracks in active.items():
                track = tracks.get(pair)
                streak = elapsed - track[0] + step_myr if track else 0.0
                for tau in TAUS_MYR:
                    merge_qualifications[(eps, tau)] += int(streak >= tau)

    last_elapsed = float(steps[-1]["elapsed_myr"])
    for eps, tracks in active.items():
        episodes[eps].extend(last_elapsed - started + step_myr for started, _last in tracks.values())

    return {
        "source_trace": path.name,
        "raw_trace_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "model_commit": summary.get("commit") or steps[0].get("commit"),
        "seed": summary.get("seed") or steps[0].get("seed"),
        "steps": len(steps),
        "step_myr": step_myr,
        "trace_complete": bool(summary),
        "elapsed_myr": last_elapsed,
        "actual_merges": merge_total,
        "runtime_seconds": summary.get("runtime_seconds"),
        "warning_note": "The replay emitted glacier depth-cap warnings; see the analysis report.",
        "continuous_lock_episode_sensitivity": [
            {
                "epsilon_cm_per_yr": eps,
                "longest_episode_myr": max(episodes[eps], default=0.0),
                "episodes_meeting_tau": {
                    str(tau): sum(duration >= tau for duration in episodes[eps])
                    for tau in TAUS_MYR
                },
            }
            for eps in EPSILONS_CM_PER_YR
        ],
        "merge_events_meeting_preceding_continuous_lock_by_epsilon": {
            str(eps): {
                "qualifying_tau_myr": [
                    tau for tau in TAUS_MYR if merge_qualifications[(eps, tau)] > 0
                ],
                "event_count": max(merge_qualifications[(eps, tau)] for tau in TAUS_MYR),
            }
            for eps in EPSILONS_CM_PER_YR
        },
        "method_note": (
            "rates are maximum absolute nodal closing rates among at least four close continental "
            "contacts; absence of qualifying close nodes breaks the lock interval. Complete "
            "traces include all configured steps; interrupted traces provide lower bounds."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trace", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--glacier-cap-warnings", type=int, default=0)
    args = parser.parse_args()
    result = summarize(args.trace)
    result["glacier_max_depth_warnings"] = args.glacier_cap_warnings
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.out:
        args.out.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
