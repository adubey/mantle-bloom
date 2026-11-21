"""Round 2 analysis for #1001: score a craton_budget_trajectory.py replay against the pre-registration.

Usage: analyse_round2.py RUNS_ROOT TAG SUBDIR
  RUNS_ROOT/SUBDIR/seed<N>-<TAG>.jsonl   round 2 replay (both seeds)
  RUNS_ROOT/runs/seed<N>-{on,off}.jsonl  round 1 arms (off arm resumes the on arm's 190 Myr save)
Prints JSON. Speed is the 210-500 Myr mean of the all-node area-weighted speed; SE is within-run only.
"""

import json, math, sys
from pathlib import Path

S = Path(sys.argv[1])
tag = sys.argv[2]
rdir = S / sys.argv[3]


def load(p):
    return [json.loads(l) for l in open(p) if l.strip()]


def by(rows):
    return {round(r["elapsed_myr"]): r for r in rows}


def speed(b, key="all"):
    xs = [b[t][f"speed_{key}_cm_per_yr"] for t in range(210, 501, 10)]
    m = sum(xs) / len(xs)
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return m, sd / math.sqrt(len(xs))


out = {}
pooled_dec = pooled_a0 = 0.0
for seed in (326494039, 111065103):
    r2 = by(load(rdir / f"seed{seed}-{tag}.jsonl"))
    on = by(load(S / "runs" / f"seed{seed}-on.jsonl"))
    offr = by(load(S / "runs" / f"seed{seed}-off.jsonl"))
    off = {**{k: v for k, v in on.items() if k <= 190}, **offr}
    a0 = r2[0]["craton_area_mkm2"]
    v0 = r2[0]["craton_volume_mkm3"]
    e = r2[500]
    ac = e.get("area_counters_mkm2", {})
    formed = ac.get("formed_area_m2")
    dec = ac.get("decratonized_area_m2")
    res = {
        "area0_mkm2": a0,
        "area500_mkm2": e["craton_area_mkm2"],
        "formed_area_pct_per_gyr": formed / a0 * 100 / 0.5
        if formed is not None
        else None,
        "decratonized_area_pct_per_gyr": dec / a0 * 100 / 0.5
        if dec is not None
        else None,
        "net_area_pct_per_gyr": (e["craton_area_mkm2"] / a0 - 1) * 100 / 0.5,
        "net_volume_pct_per_gyr": (e["craton_volume_mkm3"] / v0 - 1) * 100 / 0.5,
        "share": {
            t: round(r2[t]["craton_share_of_genuine_continental"], 4)
            for t in (0, 200, 300, 400, 500)
        },
        "share_on": {
            t: round(on[t]["craton_share_of_genuine_continental"], 4)
            for t in (200, 500)
        },
        "share_off": {
            t: round(off[t]["craton_share_of_genuine_continental"], 4)
            for t in (200, 500)
        },
        "area500_vs_on_pct": (e["craton_area_mkm2"] / on[500]["craton_area_mkm2"] - 1)
        * 100,
        "area500_on": on[500]["craton_area_mkm2"],
        "max_abs_balance_error_mkm3": max(
            abs(r["balance_error_mkm3"]) for r in r2.values()
        ),
        "peak_rss_gb": max(r["peak_rss_gb"] for r in r2.values()),
        "cores": {t: r2[t]["craton_core_count"] for t in (0, 200, 500)},
    }
    if dec is not None:
        pooled_dec += dec
        pooled_a0 += a0
    if "formation_credit_earned_m2" in ac:
        resid = [
            abs(
                r["area_counters_mkm2"]["formation_credit_earned_m2"]
                - r["area_counters_mkm2"]["formed_area_m2"]
                - r["formation_credit_mkm2"]
                - r["area_counters_mkm2"]["formation_credit_discarded_m2"]
            )
            for r in r2.values()
        ]
        res["credit"] = {
            "earned": ac["formation_credit_earned_m2"],
            "formed": formed,
            "banked": e["formation_credit_mkm2"],
            "discarded": ac["formation_credit_discarded_m2"],
            "max_abs_residual_mkm2": max(resid),
            "formed_le_earned_all": all(
                r["area_counters_mkm2"]["formed_area_m2"]
                <= r["area_counters_mkm2"]["formation_credit_earned_m2"] + 1e-9
                for r in r2.values()
            ),
        }
    first_formed = next(
        (
            t
            for t in sorted(r2)
            if r2[t].get("area_counters_mkm2", {}).get("formed_area_m2", 0) > 0
        ),
        None,
    )
    res["first_formation_myr_checkpoint"] = first_formed
    for arm, b in (("on", on), ("off", off), ("r2", r2)):
        m, se = speed(b)
        res[f"speed_{arm}"] = {"mean": m, "se": se}
    for ref in ("on", "off"):
        d = res["speed_r2"]["mean"] - res[f"speed_{ref}"]["mean"]
        se = math.hypot(res["speed_r2"]["se"], res[f"speed_{ref}"]["se"])
        res[f"speed_r2_minus_{ref}"] = {
            "diff": d,
            "z": d / se,
            "pct": d / res[f"speed_{ref}"]["mean"] * 100,
        }
    res["speed_traj"] = {
        t: (
            round(on[t]["speed_all_cm_per_yr"], 3),
            round(r2[t]["speed_all_cm_per_yr"], 3),
        )
        for t in range(0, 501, 50)
    }
    out[str(seed)] = res
out["pooled_decratonized_pct_per_gyr"] = (
    pooled_dec / pooled_a0 * 100 / 0.5 if pooled_a0 else None
)
print(json.dumps(out, indent=1))
