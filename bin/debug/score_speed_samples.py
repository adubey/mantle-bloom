"""Score recorded model plate-speed samples against the PB2002 reference (#998, step 3).

Implements the criterion fixed in `analysis/issue998/README.md` section 2: the length-weighted mean joint
log-likelihood of boundary class and relative speed under Earth's per-class lognormals, rejected below the 5th
percentile of the same score over random 50-pair subsets of Earth's plate pairs.

    python3 bin/debug/score_speed_samples.py --reference analysis/issue998/earth-reference.json \\
        --steps <captures>/pb2002/PB2002_steps.dat.txt --runs <captures>/b123-brake-off \\
        --out analysis/issue998/scores-b123-brake-off.json

`--runs` is a directory of `<seed>.json` files from `plate_speed_samples.py`. numpy only.
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import earth_speed_reference as earth  # noqa: E402

CLASSES = ("converging", "diverging", "sliding")
SUBSET_PAIRS = 50
SUBSETS = 2000
REJECT_QUANTILE = 0.05
FLOOR = 0.05


def make_scorer(reference: dict):
    fits = {c: reference["by_class"][c]["lognormal"] for c in CLASSES}
    log_pi = {c: np.log(reference["by_class"][c]["length_share"]) for c in CLASSES}

    def node_loglik(speed: np.ndarray, cls: np.ndarray) -> np.ndarray:
        ln = np.log(np.maximum(speed, FLOOR))
        out = np.empty(len(ln))
        for c in CLASSES:
            m = cls == c
            mu, sigma = fits[c]["mu"], fits[c]["sigma"]
            out[m] = log_pi[c] - 0.5 * np.log(2 * np.pi * sigma**2) - (ln[m] - mu) ** 2 / (2 * sigma**2) - ln[m]
        return out

    return node_loglik


def classes(speed: np.ndarray, closing: np.ndarray) -> np.ndarray:
    return np.where(np.abs(closing) < earth.SLIDING_OBLIQUITY * speed, "sliding", np.where(closing > 0, "converging", "diverging"))


def score_nodes(node_loglik, km: np.ndarray, speed: np.ndarray, closing: np.ndarray) -> dict:
    cls = classes(speed, closing)
    ll = node_loglik(speed, cls)
    total = km.sum()
    return {"L": float((km * ll).sum() / total), "ll_node": ll, "cls": cls}


def earth_threshold(rows: list[dict], node_loglik, rng: np.random.Generator) -> dict:
    by_pair = defaultdict(list)
    for r in rows:
        by_pair[r["pair"]].append(r)
    pairs = list(by_pair)
    scores = []
    for _ in range(SUBSETS):
        pick = [pairs[i] for i in rng.choice(len(pairs), SUBSET_PAIRS, replace=False)]
        flat = [r for p in pick for r in by_pair[p]]
        km = np.array([r["km"] for r in flat])
        speed = np.array([r["speed"] for r in flat])
        closing = np.array([r["closing"] for r in flat])
        scores.append(score_nodes(node_loglik, km, speed, closing)["L"])
    scores = np.array(scores)
    return {"L_star": float(np.quantile(scores, REJECT_QUANTILE)), "earth_L_median": float(np.median(scores)),
            "earth_L_quantiles": {str(q): float(np.quantile(scores, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}}


def class_stats(km, speed, cls) -> dict:
    total = km.sum()
    out = {}
    for c in CLASSES:
        m = cls == c
        w = km[m]
        out[c] = {"share": float(w.sum() / total),
                  "under_2": float(w[speed[m] < 2].sum() / w.sum()) if w.sum() else None,
                  "over_8": float(w[speed[m] > 8].sum() / w.sum()) if w.sum() else None,
                  "median": float(np.interp(0.5, (np.cumsum(w[np.argsort(speed[m])]) - 0.5 * w[np.argsort(speed[m])]) / w.sum(), np.sort(speed[m]))) if w.sum() else None}
    return out


def score_sample(sample: dict, node_loglik, seed: int, threshold: float) -> tuple[dict, list[dict]]:
    km, speed, closing, owner = [], [], [], []
    for ci, contact in enumerate(sample["contacts"]):
        for w, v, c in contact["nodes"]:
            km.append(w), speed.append(v), closing.append(c), owner.append(ci)
    km, speed, closing, owner = map(np.array, (km, speed, closing, owner))
    scored = score_nodes(node_loglik, km, speed, closing)
    areas = sample["areas_km2"]
    total_area = sum(areas.values())
    contacts = []
    for ci, contact in enumerate(sample["contacts"]):
        m = owner == ci
        if not m.any() or km[m].sum() < earth.MIN_CONTACT_KM:
            continue
        a, b = areas[str(contact["a"])] / total_area, areas[str(contact["b"])] / total_area
        contacts.append({"seed": seed, "myr": sample["myr"], "km": float(km[m].sum()), "speed": float((km[m] * speed[m]).sum() / km[m].sum()),
                         "L": float((km[m] * scored["ll_node"][m]).sum() / km[m].sum()), "larger": max(a, b), "smaller": min(a, b)})
    row = {"seed": seed, "myr": sample["myr"], "L": scored["L"], "rejected": scored["L"] < threshold, "plates": len(areas), "contacts": len(contacts),
           "boundary_km": float(km.sum()), "classes": class_stats(km, speed, scored["cls"]),
           "median_contact_speed": float(np.median([c["speed"] for c in contacts])) if contacts else None}
    return row, contacts


def terciles(contacts: list[dict]) -> list[dict]:
    """Mean per-contact L by (larger, smaller) plate-area tercile, terciles of the pooled model contacts."""
    big = np.array([c["larger"] for c in contacts])
    small = np.array([c["smaller"] for c in contacts])
    eb, es = np.quantile(big, [1 / 3, 2 / 3]), np.quantile(small, [1 / 3, 2 / 3])
    rows = []
    for i, nb in enumerate(("small", "mid", "large")):
        for j, ns in enumerate(("small", "mid", "large")):
            sel = [c for c, b, s in zip(contacts, big, small) if np.searchsorted(eb, b) == i and np.searchsorted(es, s) == j]
            if sel:
                km = np.array([c["km"] for c in sel])
                rows.append({"larger_plate": nb, "smaller_plate": ns, "contacts": len(sel), "length_share": float(km.sum() / sum(c["km"] for c in contacts)),
                             "mean_L": float(np.average([c["L"] for c in sel], weights=km)),
                             "median_speed": float(np.median([c["speed"] for c in sel]))})
    return [{"larger_edges": [float(x) for x in eb], "smaller_edges": [float(x) for x in es]}, *rows]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--steps", type=Path, required=True)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--seed", type=int, default=998)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    node_loglik = make_scorer(json.loads(args.reference.read_text()))
    threshold = earth_threshold(earth.read_steps(args.steps), node_loglik, rng)
    rows, contacts = [], []
    for path in sorted(args.runs.glob("*.json")):
        run = json.loads(path.read_text())
        for sample in run["samples"]:
            row, cs = score_sample(sample, node_loglik, run["seed"], threshold["L_star"])
            rows.append(row)
            contacts += cs
    times = sorted({r["myr"] for r in rows})
    by_time = [{"myr": t, "samples": sum(r["myr"] == t for r in rows), "accepted": sum(r["myr"] == t and not r["rejected"] for r in rows),
                "acceptance_rate": float(np.mean([not r["rejected"] for r in rows if r["myr"] == t])),
                "L_median": float(np.median([r["L"] for r in rows if r["myr"] == t])),
                "median_contact_speed": float(np.median([r["median_contact_speed"] for r in rows if r["myr"] == t]))} for t in times]
    result = {"runs": str(args.runs.name), "threshold": threshold, "subset_pairs": SUBSET_PAIRS, "subsets": SUBSETS,
              "n_seeds": len({r["seed"] for r in rows}), "by_time": by_time, "samples": rows,
              "area_terciles": {str(t): terciles([c for c in contacts if c["myr"] == t]) for t in times if any(c["myr"] == t for c in contacts)}}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    print("L* =", round(threshold["L_star"], 3), "Earth median", round(threshold["earth_L_median"], 3))
    for b in by_time:
        print(f"{b['myr']:5.0f} Myr  accepted {b['accepted']:2d}/{b['samples']}  L median {b['L_median']:.2f}  median contact speed {b['median_contact_speed']:.2f}")


if __name__ == "__main__":
    main()
