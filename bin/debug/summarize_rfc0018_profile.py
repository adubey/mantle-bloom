"""Summarize profile_rfc0018_candidate_search.py output into compact JSON (#848).

Usage (repo root):
    backend/.venv/bin/python bin/debug/summarize_rfc0018_profile.py \
        ~/mantle-bloom-captures/issue848/profile.jsonl analysis/issue848/profile-summary.json

Per mapping it reports total and median wall time, sparse-solve counts, each helper's
inclusive share of time, how many variants used all 8 fit iterations, and, for windows of 1, 3
and 5 steps, the share of opportunities and of search time on plates that ran
``quad_tectonics._open_rift`` within that window (the rift-candidate marker proxy).
"""

import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

WINDOWS = (1, 3, 5)


def summarize(rows: list[dict]) -> dict:
    out = {"rows": len(rows), "measured_steps": max(r["step"] for r in rows), "mappings": {}}
    for mapping in sorted({r["mapping"] for r in rows}):
        group = [r for r in rows if r["mapping"] == mapping]
        ok = [r for r in group if r["status"] == "ok"]
        total = sum(r["total_s"] for r in group)
        seconds = defaultdict(float)
        for r in group:
            for name, value in r["seconds"].items():
                seconds[name] += value
        variants = sum(r["variants"] for r in ok)
        capped = sum(r["variants_at_iteration_cap"] for r in ok)
        windows = {}
        for k in WINDOWS:
            near = [r for r in group if r["steps_since_open_rift"] is not None and r["steps_since_open_rift"] <= k]
            windows[str(k)] = {
                "opportunities": len(near),
                "opportunity_share": round(len(near) / len(group), 4),
                "time_share": round(sum(r["total_s"] for r in near) / total, 4),
            }
        out["mappings"][mapping] = {
            "opportunities": len(group),
            "not_ok": len(group) - len(ok),
            "total_s": round(total, 1),
            "median_s": round(statistics.median(r["total_s"] for r in group), 3),
            "sparse_solves": sum(r["sparse_solves"] for r in group),
            "median_solves": statistics.median(r["sparse_solves"] for r in group),
            "median_nodes": statistics.median(r["nodes"] for r in group),
            "max_nodes": max(r["nodes"] for r in group),
            "inclusive_time_share": {k: round(v / total, 4) for k, v in sorted(seconds.items(), key=lambda x: -x[1])},
            "variants": variants,
            "variants_at_iteration_cap": capped,
            "iteration_cap_share": round(capped / variants, 4),
            "usable_variants": sum(r["usable_variants"] for r in ok),
            "selected_accepted": sum(1 for r in ok if r["selected_accepted"]),
            "extension_window_steps": windows,
        }
    return out


def main() -> None:
    source, target = (Path(a).expanduser() for a in sys.argv[1:3])
    rows = [json.loads(line) for line in source.read_text().splitlines() if line.strip()]
    summary = summarize(rows)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
