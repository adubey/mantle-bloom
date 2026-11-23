"""Summarise `big_plate_rift_probe.py` outputs: which rift gate holds big plates back, and how craton share relates (#995).

    python3 -I bin/debug/rift_state_analysis.py --dir analysis/issue995 --out analysis/issue995/rift-state-analysis.json
"""

import argparse
import collections
import glob
import json
from pathlib import Path

import numpy as np
from scipy import stats

BIG = 0.10  # share of the planet's area


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    out = {}
    for setting in ("b123-brake-off", "brake-o0.25-c1"):
        files = sorted(glob.glob(str(args.dir / f"rift-state-{setting}-*.json")))
        gates = {"big": collections.Counter(), "small": collections.Counter()}
        rows, rejected, tests, splits = [], collections.Counter(), 0, 0
        for f in files:
            d = json.load(open(f))
            tests += d["counters"]["break_tests"]
            splits += len(d["splits"])
            rejected.update(d["counters"]["rejected_because"])
            for s in d["samples"]:
                for p in s["plates"]:
                    gates["big" if p["area_share"] >= BIG else "small"][p["gate"]] += 1
                    rows.append(p)
        has_line = [r for r in rows if "broken_share" in r]
        a = np.array([r["area_share"] for r in has_line]); b = np.array([r["broken_share"] for r in has_line])
        c = np.array([r["craton_share"] for r in has_line]); age = np.array([r["age_myr"] for r in has_line])
        allc = np.array([r["craton_share"] for r in rows]); alla = np.array([r["area_share"] for r in rows])
        big = a >= BIG
        out[setting] = {
            "runs": len(files), "splits": splits, "break_tests": tests, "rejected_because": dict(rejected),
            "gates_big": dict(gates["big"]), "gates_small": dict(gates["small"]),
            "broken_share_mean_big": round(float(b[big].mean()), 3), "broken_share_mean_small": round(float(b[~big].mean()), 3),
            "craton_share_mean_big": round(float(allc[alla >= BIG].mean()), 3), "craton_share_mean_small": round(float(allc[alla < BIG].mean()), 3),
            "spearman_broken_vs_area": round(float(stats.spearmanr(a, b)[0]), 2),
            "spearman_broken_vs_craton": round(float(stats.spearmanr(c, b)[0]), 2),
            "spearman_broken_vs_age": round(float(stats.spearmanr(age, b)[0]), 2),
            "spearman_area_vs_craton": round(float(stats.spearmanr(a, c)[0]), 2),
            "passes_share_by_craton_tercile": [
                round(float(np.mean([r["gate"] == "passes" for r in has_line if lo <= r["craton_share"] <= hi])), 3)
                for lo, hi in ((0, 0.2), (0.2, 0.5), (0.5, 1.0))],
        }
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
