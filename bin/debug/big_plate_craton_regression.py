"""Separate craton share from plate size in the big-plate break-up statistic (#995).

Logistic regression of "splits within the window" on area share and craton share (standardised), over big plate-samples
pooled across seeds of the given arms, with a seed-level bootstrap for the coefficients. Craton share is correlated with
size, so the pooled craton effect in `breakup-stat.json` cannot be read as causal on its own.

    python3 -I bin/debug/big_plate_craton_regression.py --dir analysis/issue995 --out analysis/issue995/craton-regression.json \\
        brake-o0.25-c1 brake-o0.25-c1-heal250 brake-o0.25-c1-kcap07
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from big_plate_breakup_stat import arm_runs


def records(directory: Path, setting: str, big: float, window: float) -> list[tuple]:
    out = []
    for run in arm_runs(directory, setting):
        splits: dict[int, list[float]] = {}
        for s in run["splits"]:
            splits.setdefault(s["plate"], []).append(s["myr"])
        end = run["samples"][-1]["myr"]
        for sample in run["samples"]:
            t = sample["myr"]
            if t + window > end:
                continue
            for p in sample["plates"]:
                if p["area_share"] >= big:
                    hit = any(t < x <= t + window for x in splits.get(p["plate"], []))
                    out.append((run["seed"], p["area_share"], p["craton_share"], float(hit)))
    return out


def fit(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    X = np.column_stack([np.ones(len(x)), x])

    def nll(b: np.ndarray) -> float:
        z = X @ b
        return float(np.sum(np.logaddexp(0, z) - y * z) + 1e-3 * np.sum(b[1:] ** 2))

    return minimize(nll, np.zeros(X.shape[1]), method="BFGS").x


def analyse(rows: list[tuple], rng: np.random.Generator) -> dict:
    arr = np.array([(r[1], r[2], r[3]) for r in rows])
    seeds = np.array([r[0] for r in rows])
    mu, sd = arr[:, :2].mean(0), arr[:, :2].std(0)
    x, y = (arr[:, :2] - mu) / sd, arr[:, 2]
    full = fit(x, y)
    alone = fit(x[:, 1:], y)
    uniq = np.unique(seeds)
    boots = []
    for _ in range(300):
        pick = rng.choice(uniq, len(uniq))
        idx = np.concatenate([np.where(seeds == s)[0] for s in pick])
        boots.append(fit(x[idx], y[idx]))
    lo, hi = np.quantile(boots, [0.025, 0.975], axis=0)
    return {"samples": len(y), "corr_area_craton": round(float(np.corrcoef(x.T)[0, 1]), 2),
            "craton_alone_coef": round(float(alone[1]), 2),
            "coef_per_sd": {"area": round(float(full[1]), 2), "craton": round(float(full[2]), 2)},
            "bootstrap_95": {"area": [round(float(lo[1]), 2), round(float(hi[1]), 2)],
                             "craton": [round(float(lo[2]), 2), round(float(hi[2]), 2)]}}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--big", type=float, default=0.10)
    ap.add_argument("--window", type=float, default=20.0)
    ap.add_argument("settings", nargs="+")
    args = ap.parse_args()
    rng = np.random.default_rng(995)
    out = {s: analyse(records(args.dir, s, args.big, args.window), rng) for s in args.settings}
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
