#!/usr/bin/env python3
"""RFC 0015 gate G1 on the march (#766): apply the preregistered rules to the probe output.

Reads `<dir>/<save>/g1.jsonl`, written by the gate branch's `replay_g0.py --g1-myr` (its
`bin/debug/g1_march.py`), and `<dir>/<save>.json` for the replay's checkpoints. Per save it
reports the volume-weighted τ distribution and G1's pass rule, the substep limit from τ/4, the
convergence check by cap, turnover, frictional heat beside basal melt, and lake water under ice.
See `analysis/issue766/README.md` for the rules.

    backend/.venv/bin/python bin/debug/summarize_g1_march.py <dir> --out analysis/issue766/g1.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

STEP_LIMIT_FRACTION = 0.1
FINE_WINDOW_YR = 10_000.0
LIMIT_ROUNDING_YR = (1000.0, 500.0, 250.0, 125.0)
CONVERGENCE_LIMIT = 0.01
VANISHED_FRACTION = 0.01


def body_tau(body: dict) -> float:
    """The median of the four perturbations' e-folding times; a perturbation that never
    reaches 1/e counts as beyond the fine window."""
    taus = [c["tau_efold_yr"] if c["tau_efold_yr"] is not None else np.inf for c in body["fine"].values()]
    return float(np.median(taus))


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order]) / weights.sum()
    return float(values[order][np.searchsorted(cumulative, q)])


def summarise_save(records: list[dict], checkpoints: list[dict]) -> dict:
    bodies = [(r, b) for r in records for b in r["bodies"] if b["body"] != "small"]
    tau = np.array([body_tau(b) for _, b in bodies])
    volume = np.array([b["opening_km3"] for _, b in bodies])
    step_years = records[0]["step_years"]
    threshold = STEP_LIMIT_FRACTION * step_years
    over = float(volume[tau > threshold].sum() / volume.sum())
    p5 = weighted_quantile(tau, volume, 0.05)
    limit = next((r for r in LIMIT_ROUNDING_YR if r <= p5 / 4.0), LIMIT_ROUNDING_YR[-1])
    caps = list(bodies[0][1]["convergence_gap"])
    convergence = {}
    for cap in caps:
        larger = cap.split("->")[0]
        # A body the larger cap melts away has no meaningful ratio: count it as failing, apart.
        vanished = sum(1 for _, b in bodies if cap in b["convergence_gap"]
                       and b["closing_by_cap_km3"][larger] < VANISHED_FRACTION * max(b["opening_km3"], 1e-9))
        gaps = np.array([abs(b["convergence_gap"][cap]) for _, b in bodies if cap in b["convergence_gap"]
                         and b["closing_by_cap_km3"][larger] >= VANISHED_FRACTION * max(b["opening_km3"], 1e-9)])
        convergence[cap] = {"passes": bool((gaps <= CONVERGENCE_LIMIT).all()) and vanished == 0,
                            "max_gap_pct": round(100 * float(gaps.max()), 2) if len(gaps) else None,
                            "bodies_over_1pct": int((gaps > CONVERGENCE_LIMIT).sum()), "bodies_vanished": vanished,
                            "bodies": int(len(gaps)) + vanished,
                            "volume_share_over_1pct": round(sum(b["opening_km3"] for _, b in bodies if cap in b["convergence_gap"]
                                                                and abs(b["convergence_gap"][cap]) > CONVERGENCE_LIMIT) / volume.sum(), 4)}
    largest_passing = next((int(c.split("->")[0]) for c in caps if convergence[c]["passes"]), None)
    turnover = np.array([b["turnover_yr"] or np.nan for _, b in bodies], dtype=float)
    finite = np.isfinite(tau)
    probed = {r["target_myr"] for r in records}
    unprobed = [c["step_seconds_mean"] for c in checkpoints[1:] if round(c["elapsed_myr"], 2) not in probed]
    return {
        "probed_steps_myr": sorted(probed),
        "bodies_measured": len(bodies),
        "tau_yr_weighted": {f"p{int(q * 100)}": weighted_quantile(tau, volume, q) for q in (0.05, 0.5, 0.95)},
        "tau_yr_range": [float(tau[finite].min()), float(tau[finite].max())] if finite.any() else None,
        "bodies_tau_beyond_window": int((~finite).sum()),
        "volume_share_tau_over_threshold": round(over, 4),
        "pass_rule": "pass" if over <= 0.05 else ("falsified" if over > 0.20 else "fail"),
        "substep_limit_yr": limit,
        "convergence": convergence,
        "largest_cap_converging_yr": largest_passing,
        "turnover_yr_weighted_median": weighted_quantile(np.nan_to_num(turnover, nan=np.inf), volume, 0.5),
        "frictional_over_basal": [round(r["frictional_melt_km3"] / r["basal_melt_km3"], 3) for r in records if r["frictional_melt_km3"] and r["basal_melt_km3"]],
        "lake_under_ice_share_max": max((c.get("lake_under_ice_share") or 0.0) for c in checkpoints),
        "ice_volume_mkm3": [round(checkpoints[0]["ice_volume_km3"] / 1e6, 2), round(max(c["ice_volume_km3"] for c in checkpoints) / 1e6, 2), round(checkpoints[-1]["ice_volume_km3"] / 1e6, 2)],
        "unprobed_step_seconds_median": round(float(np.median(unprobed)), 1),
        "probe_seconds": [r["probe_seconds"] for r in records],
        "per_body": [{"step_myr": r["target_myr"], "body": b["body"], "lat": b["lat"], "opening_km3": b["opening_km3"], "tau_yr": body_tau(b),
                      "perturbation_tau_yr": {k: c["tau_efold_yr"] for k, c in b["fine"].items()},
                      "perturbation_closing_ratio": {k: c["closing_ratio"] for k, c in b["fine"].items()},
                      "turnover_yr": b["turnover_yr"], "convergence_gap": b["convergence_gap"]} for r, b in bodies],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", type=Path)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    report = {}
    for probes in sorted(args.runs.glob("*/g1.jsonl")):
        save = probes.parent.name
        records = [json.loads(line) for line in probes.read_text().splitlines() if line.strip()]
        checkpoints = json.loads((args.runs / f"{save}.json").read_text())["checkpoints"]
        report[save] = summarise_save(records, checkpoints)
    args.out.write_text(json.dumps(report, indent=1) + "\n")
    for save, row in report.items():
        print(save, json.dumps({k: v for k, v in row.items() if k != "per_body"}))


if __name__ == "__main__":
    main()
