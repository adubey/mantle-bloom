"""Attribute the #998 speed rejections to driving-force and layout features (step 4).

Reads the per-run samples (`plate_speed_samples.py`, with `plate_torque_breakdown.sample` rows) and the scores
(`score_speed_samples.py`), and reports, per sample time:
  * pooled features, one value per world: driving-torque shares (slab pull, ridge push, mantle flow; norms, summed
    over plates), the area share of plates with a subducting node, the pulled share of boundary length, the share
    of boundary length opening past the gate, plate count, mean and largest plate area share, area-weighted speed,
    and the speed with mantle flow alone (`counterfactual_speed_cm_per_yr.mantle_flow_only`);
  * the paired change of each feature from t=0 (accepted) to each later time, with a bootstrap interval over seeds;
  * across seeds, the rank correlation of each feature with the score L, and for a time with both outcomes the
    standardized difference (Cohen's d) between accepted and rejected;
  * per-plate, at each time: rank correlation of plate speed with each torque term per unit area.

    python3 bin/debug/attribute_speed_rejection.py --runs <captures>/b123-brake-off \\
        --scores analysis/issue998/scores-b123-brake-off.json --out analysis/issue998/attribution-b123-brake-off.json
"""

import argparse
import json
from pathlib import Path

import numpy as np

TERMS = ("slab_pull", "ridge_push", "mantle_flow", "basal_drag", "slab_drag", "contact_brake")
DRIVING = ("slab_pull", "ridge_push", "mantle_flow")
RESISTING = ("basal_drag", "slab_drag", "contact_brake")
BOOT = 2000


def ranks(x: np.ndarray) -> np.ndarray:
    order = np.argsort(x, kind="stable")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    for v in np.unique(x):  # average ties
        m = x == v
        r[m] = r[m].mean()
    return r


def spearman(x, y) -> float | None:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4 or x.std() == 0 or y.std() == 0:
        return None
    return float(np.corrcoef(ranks(x), ranks(y))[0, 1])


def world_features(sample: dict) -> dict:
    plates = sample["breakdown"]["plates"]
    area = np.array([p["area_km2"] for p in plates])
    total = area.sum()
    norm = {t: sum(p["torque_norm_nm"].get(t, 0.0) for p in plates) for t in TERMS}
    drive = sum(norm[t] for t in DRIVING) or 1.0
    band = sum(p["boundary_km"]["band_total"] for p in plates) or 1.0
    speed = np.array([p["speed_cm_per_yr"]["all"] for p in plates])
    flow = np.array([p["counterfactual_speed_cm_per_yr"]["mantle_flow_only"]["all"] for p in plates])
    out = {f"driving_share_{t}": norm[t] / drive for t in DRIVING}
    # Resisting terms are signed along each plate's motion axis (`torque_along_axis_nm`), not in `torque_norm_nm`.
    along = {t: sum(p["torque_along_axis_nm"].get(t, 0.0) for p in plates) for t in TERMS}
    gross = sum(along[t] for t in DRIVING) or 1.0
    for t in ("basal_drag", "slab_drag", "contact_brake"):
        out[f"resisting_over_driving_{t}"] = -along[t] / gross
    out["slab_pull_cancellation_kept"] = along["slab_pull"] / norm["slab_pull"] if norm["slab_pull"] else 0.0
    out["slab_bearing_area_share"] = float(sum(a for a, p in zip(area, plates) if p["g0"]["has_subducting_node"]) / total)
    out["slab_attached_area_share_10pct"] = float(sum(a for a, p in zip(area, plates) if p["g0"]["consumed_10myr_km2"] >= 0.1 * a) / total)
    out["pulled_boundary_share"] = sum(p["boundary_km"]["down_going_closing_pulled"] for p in plates) / band
    out["opening_past_gate_share"] = sum(p["boundary_km"]["opening_past_gate"] for p in plates) / band
    out["plates"] = len(plates)
    out["mean_plate_area_share"] = float(1.0 / len(plates))
    out["largest_plate_area_share"] = float(area.max() / total)
    out["speed_area_weighted_cm_per_yr"] = float((speed * area).sum() / total)
    out["flow_only_speed_cm_per_yr"] = float((flow * area).sum() / total)
    out["speed_over_flow_only"] = out["speed_area_weighted_cm_per_yr"] / out["flow_only_speed_cm_per_yr"]
    return out


def boot_ci(values: np.ndarray, rng: np.random.Generator) -> list[float]:
    means = [values[rng.integers(0, len(values), len(values))].mean() for _ in range(BOOT)]
    return [float(v) for v in np.percentile(means, [2.5, 97.5])]


def plate_term(p: dict, t: str) -> float:
    return p["torque_along_axis_nm"].get(t, 0.0) if t in RESISTING else p["torque_norm_nm"].get(t, 0.0)


def per_plate_correlations(samples_by_seed: dict, myr: float) -> dict:
    rows = []
    for sample in samples_by_seed.values():
        if sample["myr"] == myr:
            rows += sample["breakdown"]["plates"]
    speed = [p["speed_cm_per_yr"]["all"] for p in rows]
    area = [p["area_km2"] for p in rows]
    out = {"plates": len(rows)}
    for t in TERMS:
        out[f"rho_speed_vs_{t}_per_area"] = spearman(speed, [plate_term(p, t) / a for p, a in zip(rows, area)])
    out["rho_speed_vs_area"] = spearman(speed, area)
    out["rho_speed_vs_slab_bearing"] = spearman(speed, [float(p["g0"]["has_subducting_node"]) for p in rows])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", type=Path, required=True)
    ap.add_argument("--scores", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    rng = np.random.default_rng(998)
    scores = {(s["seed"], s["myr"]): s for s in json.loads(args.scores.read_text())["samples"]}
    feats: dict[tuple, dict] = {}
    by_time_samples: dict[float, dict] = {}
    for path in sorted(args.runs.glob("*.json")):
        run = json.loads(path.read_text())
        for sample in run["samples"]:
            feats[(run["seed"], sample["myr"])] = world_features(sample)
            by_time_samples.setdefault(sample["myr"], {})[run["seed"]] = sample
    times = sorted({t for _, t in feats})
    seeds = sorted({s for s, _ in feats})
    names = list(next(iter(feats.values())))
    result = {"seeds": len(seeds), "by_time": {}, "paired_change_from_t0": {}, "across_seeds": {}, "per_plate": {}}
    for t in times:
        result["by_time"][str(t)] = {n: {"mean": float(np.mean([feats[(s, t)][n] for s in seeds])), "min": float(np.min([feats[(s, t)][n] for s in seeds])),
                                         "max": float(np.max([feats[(s, t)][n] for s in seeds]))} for n in names}
        result["per_plate"][str(t)] = per_plate_correlations(by_time_samples[t], t)
    for t in times[1:]:
        result["paired_change_from_t0"][str(t)] = {}
        for n in names:
            d = np.array([feats[(s, t)][n] - feats[(s, 0.0)][n] for s in seeds])
            result["paired_change_from_t0"][str(t)][n] = {"mean": float(d.mean()), "ci95": boot_ci(d, rng), "seeds_increasing": int((d > 0).sum())}
    for t in times:
        L = [scores[(s, t)]["L"] for s in seeds]
        rej = np.array([scores[(s, t)]["rejected"] for s in seeds])
        block = {}
        for n in names:
            x = np.array([feats[(s, t)][n] for s in seeds])
            entry = {"rho_with_L": spearman(x, L)}
            if 0 < rej.sum() < len(rej) and min(rej.sum(), (~rej).sum()) >= 2:
                a, r = x[~rej], x[rej]
                pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(r) - 1) * r.var(ddof=1)) / (len(a) + len(r) - 2))
                entry["cohens_d_rejected_minus_accepted"] = float((r.mean() - a.mean()) / pooled) if pooled else None
            block[n] = entry
        result["across_seeds"][str(t)] = {"accepted": int((~rej).sum()), "rejected": int(rej.sum()), "features": block}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    key = ["driving_share_slab_pull", "driving_share_ridge_push", "driving_share_mantle_flow", "slab_bearing_area_share", "pulled_boundary_share",
           "plates", "speed_area_weighted_cm_per_yr", "flow_only_speed_cm_per_yr"]
    print("mean over seeds by time:")
    print("Myr  " + "  ".join(f"{k[:18]:>18s}" for k in key))
    for t in times:
        print(f"{t:4.0f} " + "  ".join(f"{result['by_time'][str(t)][k]['mean']:18.3f}" for k in key))


if __name__ == "__main__":
    main()
