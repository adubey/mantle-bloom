"""Does plate count or the largest plate's share explain the fall in mean plate speed over time? (#995)

Reads the ``slowdown-<setting>-<seed>.json`` files that ``plate_contact_graph.py record`` writes and prints, per
setting, within-run (each run's own mean removed) correlations of the area-weighted NNR mean speed with time, plate
count and the largest plate's area share, plus per-seed speed trends and a pooled per-plate speed-by-area table.

    python3 -I bin/debug/slowdown_analysis.py --dir analysis/issue995 --out analysis/issue995/slowdown-analysis.json
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy import stats


def corr(a, b) -> float:
    return round(float(np.corrcoef(a, b)[0, 1]), 2)


def residual(y, x):
    return y - np.polyval(np.polyfit(x, y, 1), x)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = {}
    for setting in ("b123-brake-off", "brake-o0.25-c1"):
        T, V, B, N, per_plate, trends = [], [], [], [], [], []
        files = sorted(glob.glob(str(args.dir / f"slowdown-{setting}-*.json")))
        for f in files:
            q = []
            for s in json.load(open(f))["samples"]:
                a = np.array([v["area_km2"] for v in s["plates"].values()])
                sp = np.array([v["speed_nnr_cm_per_yr"] for v in s["plates"].values()])
                q.append((s["myr"], s["mean_speed_nnr_cm_per_yr"], a.max() / a.sum(), len(a)))
                per_plate += list(zip(a, sp))
            q = np.array(q)
            trends.append(float(np.mean(q[-5:, 1]) - np.mean(q[:5, 1])))  # last 50 Myr minus first 50 Myr
            q = q - q.mean(0)
            T += list(q[:, 0]); V += list(q[:, 1]); B += list(q[:, 2]); N += list(q[:, 3])
        T, V, B, N = map(np.array, (T, V, B, N))
        area, spd = (np.array(x) for x in zip(*per_plate))
        quart = np.quantile(area, [0, 0.25, 0.5, 0.75, 1])
        trends = np.array(trends)
        result[setting] = {
            "runs": len(files),
            "speed_vs_time": corr(V, T), "largest_share_vs_time": corr(B, T), "speed_vs_largest_share": corr(V, B),
            "speed_vs_time_given_largest_share": corr(residual(V, B), residual(T, B)),
            "speed_vs_count_given_time": corr(residual(V, T), residual(N, T)),
            "count_vs_time": corr(N, T),
            "per_plate_speed_vs_area_spearman": round(float(stats.spearmanr(area, spd)[0]), 2),
            "speed_by_area_quartile": [round(float(spd[(area >= lo) & (area <= hi)].mean()), 1) for lo, hi in zip(quart[:-1], quart[1:])],
            "seeds_slower_last_vs_first_50myr": int((trends < 0).sum()),
            "median_change_cm_per_yr": round(float(np.median(trends)), 2),
            "mean_change_cm_per_yr": round(float(trends.mean()), 2),
        }
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
