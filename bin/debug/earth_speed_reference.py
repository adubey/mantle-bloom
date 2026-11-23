"""Earth reference for plate relative speed per boundary class and per plate-area pair (#998, step 1).

Reads the PB2002 steps file (Bird 2003, one row per boundary step, with the relative speed, the divergence rate
and the strike-slip rate in mm/yr) and the plate area shares in `analysis/issue958/earth-plate-speed.json`.

Classes (fixed before any model sample was scored, so the model is not tuned to them): by obliquity of the
relative motion to the boundary normal, *sliding* if |normal rate| < 0.5 x speed (motion within 30 degrees of
the boundary), otherwise *converging* (closing) or *diverging* (opening) by the sign of the normal rate. Bird's
own codes (SUB, OSR, ...) are reported as a cross-check.

Two reference distributions, both length-weighted (a step counts by its length in km):
  1. per class: the weighted quantiles of step relative speed, and a lognormal fit (weighted mean and sd of
     ln speed). Uncertainty is a bootstrap over *plate pairs* (steps along one contact are not independent).
  2. per contact (plate pair with total length >= MIN_CONTACT_KM): ln(contact speed) regressed on
     ln(area of the larger plate) and ln(area of the smaller plate), with a bootstrap over contacts.

Usage (python3 -I is fine; numpy only):
    python3 bin/debug/earth_speed_reference.py --steps <dir>/PB2002_steps.dat.txt \\
        --plates analysis/issue958/earth-plate-speed.json --out analysis/issue998/earth-reference.json
"""

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

MIN_CONTACT_KM = 500.0  # as bin/debug/plate_contact_graph.py (#995)
SLIDING_OBLIQUITY = 0.5  # |normal| / speed below this is sliding
BOOTSTRAP = 2000
QUANTILES = (0.1, 0.25, 0.5, 0.75, 0.9)
CODE = re.compile(r":?([A-Z]{3})\*?\s*$")


def classify(speed: float, closing: float) -> str:
    if abs(closing) < SLIDING_OBLIQUITY * speed:
        return "sliding"
    return "converging" if closing > 0 else "diverging"


def read_steps(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        t = line.split()
        if len(t) < 14:
            continue
        pair = t[1].lstrip(":")
        for sep in ("\\", "/", "-"):
            if sep in pair:
                a, b = pair.split(sep)
                break
        else:
            continue
        speed = float(t[8]) / 10.0
        closing = -float(t[10]) / 10.0  # file's divergence is positive for opening
        m = CODE.search(line)
        rows.append({"pair": tuple(sorted((a, b))), "km": float(t[6]), "speed": speed, "closing": closing,
                     "cls": classify(speed, closing), "bird": m.group(1) if m else "?"})
    return rows


def wquantile(x: np.ndarray, w: np.ndarray, qs) -> list[float]:
    order = np.argsort(x)
    x, w = x[order], w[order]
    cum = (np.cumsum(w) - 0.5 * w) / w.sum()
    return [float(np.interp(q, cum, x)) for q in qs]


def lognormal(x: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    ln = np.log(np.maximum(x, 0.05))  # PB2002 speeds are quoted to 0.1 cm/yr; floor the zeros
    mu = np.average(ln, weights=w)
    return float(mu), float(np.sqrt(np.average((ln - mu) ** 2, weights=w)))


def class_block(rows: list[dict], rng: np.random.Generator) -> dict:
    out = {}
    for cls in ("converging", "diverging", "sliding", "all"):
        sel = [r for r in rows if cls == "all" or r["cls"] == cls]
        x = np.array([r["speed"] for r in sel])
        w = np.array([r["km"] for r in sel])
        pairs = sorted({r["pair"] for r in sel})
        by_pair = {p: [r for r in sel if r["pair"] == p] for p in pairs}
        boot = []
        for _ in range(BOOTSTRAP):
            pick = [by_pair[pairs[i]] for i in rng.integers(0, len(pairs), len(pairs))]
            flat = [r for grp in pick for r in grp]
            xs = np.array([r["speed"] for r in flat])
            ws = np.array([r["km"] for r in flat])
            boot.append([*lognormal(xs, ws), *wquantile(xs, ws, (0.5,)), float((ws[xs < 2].sum() / ws.sum())), float((ws[xs > 8].sum() / ws.sum()))])
        boot = np.array(boot)
        mu, sigma = lognormal(x, w)
        ci = lambda col: [float(v) for v in np.percentile(boot[:, col], [2.5, 97.5])]  # noqa: E731
        out[cls] = {
            "n_steps": len(sel), "n_pairs": len(pairs), "length_km": float(w.sum()),
            "length_share": float(w.sum() / sum(r["km"] for r in rows)),
            "quantiles_cm_per_yr": dict(zip(map(str, QUANTILES), wquantile(x, w, QUANTILES))),
            "share_under_2": float(w[x < 2].sum() / w.sum()), "share_under_4": float(w[x < 4].sum() / w.sum()),
            "share_over_8": float(w[x > 8].sum() / w.sum()),
            "lognormal": {"mu": mu, "sigma": sigma, "mu_ci95": ci(0), "sigma_ci95": ci(1)},
            "median_ci95": ci(2), "share_under_2_ci95": ci(3), "share_over_8_ci95": ci(4),
        }
    return out


def contact_block(rows: list[dict], areas: dict, rng: np.random.Generator) -> dict:
    agg = defaultdict(lambda: [0.0, 0.0])
    for r in rows:
        agg[r["pair"]][0] += r["km"]
        agg[r["pair"]][1] += r["km"] * r["speed"]
    contacts = [(p, km, sp / km) for p, (km, sp) in agg.items() if km >= MIN_CONTACT_KM and p[0] in areas and p[1] in areas]
    big = np.array([max(areas[p[0]]["area_share"], areas[p[1]]["area_share"]) for p, _, _ in contacts])
    small = np.array([min(areas[p[0]]["area_share"], areas[p[1]]["area_share"]) for p, _, _ in contacts])
    y = np.log(np.maximum([s for _, _, s in contacts], 0.05))
    X = np.column_stack([np.ones(len(y)), np.log(big), np.log(small)])

    def fit(idx):
        beta, *_ = np.linalg.lstsq(X[idx], y[idx], rcond=None)
        resid = y[idx] - X[idx] @ beta
        return beta, float(resid.std(ddof=3))

    beta, sd = fit(np.arange(len(y)))
    boots = np.array([fit(rng.integers(0, len(y), len(y)))[0] for _ in range(BOOTSTRAP)])
    r2 = 1 - ((y - X @ beta) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    total = sum(km for _, km, _ in contacts)
    speeds = np.array([s for _, _, s in contacts])
    lengths = np.array([km for _, km, _ in contacts])
    return {
        "n_contacts": len(contacts), "contact_length_km": total,
        "model": "ln(speed cm/yr) = b0 + b1 ln(larger area share) + b2 ln(smaller area share) + N(0, sd)",
        "coefficients": {"b0": float(beta[0]), "b1_larger": float(beta[1]), "b2_smaller": float(beta[2])},
        "coefficients_ci95": {k: [float(v) for v in np.percentile(boots[:, i], [2.5, 97.5])] for i, k in enumerate(("b0", "b1_larger", "b2_smaller"))},
        "residual_sd_ln": sd, "r2": float(r2),
        "length_weighted_contact_speed_quantiles": dict(zip(map(str, QUANTILES), wquantile(speeds, lengths, QUANTILES))),
        "contacts": [{"a": p[0], "b": p[1], "km": km, "speed": s, "area_larger": float(bg), "area_smaller": float(sm)}
                     for (p, km, s), bg, sm in zip(contacts, big, small)],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--steps", type=Path, required=True)
    ap.add_argument("--plates", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=998)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    rows = read_steps(args.steps)
    areas = json.loads(args.plates.read_text())["plates"]
    bird = defaultdict(float)
    for r in rows:
        bird[r["bird"]] += r["km"]
    result = {
        "source": "PB2002 steps (Bird 2003), plate areas from analysis/issue958/earth-plate-speed.json",
        "sliding_obliquity": SLIDING_OBLIQUITY, "min_contact_km": MIN_CONTACT_KM, "bootstrap": BOOTSTRAP,
        "by_class": class_block(rows, rng),
        "bird_code_length_km": dict(bird),
        "by_contact_area": contact_block(rows, areas, rng),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for cls, b in result["by_class"].items():
        print(f"{cls:11s} steps {b['n_steps']:5d} pairs {b['n_pairs']:3d} share {b['length_share']:.2f} median {b['quantiles_cm_per_yr']['0.5']:.2f} "
              f"<2 {b['share_under_2']:.2f} >8 {b['share_over_8']:.2f} lognormal mu {b['lognormal']['mu']:.2f} sigma {b['lognormal']['sigma']:.2f}")
    c = result["by_contact_area"]
    print(c["n_contacts"], "contacts", c["coefficients"], "R2", round(c["r2"], 2), "sd", round(c["residual_sd_ln"], 2))


if __name__ == "__main__":
    main()
