#!/usr/bin/env python3
"""#842 round 3: score the probe runs against the preregistered decision rule.

Reads every `probe_822_march_cost.py` JSON in a directory and applies
`analysis/issue842/preregistration.md`: each criterion is reported against its old and its
new (gating) threshold, and the verdict is computed under both, so the final evaluation can
say whether keeping the old thresholds would change it.

File names carry the role, since the probe JSON doesn't:

- `<seed>-<arm>-<myr>.json`: an E1/E2 (diagnosis) or E5 (held-out) replay, 40 x 0.25 Myr
  or 10 x 1 Myr. `<arm>` is `d2` or `trend-1k`; `<myr>` is `0.25` or `1`.
- `e3-<step>-<arm>-<scale>.json`: one E3 member; `<arm>` is `trend-1k`, `d2` or `d2-10k`
  (the positive control at 10 kyr substeps), `<scale>` the opening-ice factor.
- `e4-<seed>-sub<years>-<myr>.json`: one E4 `d2` convergence run.
- `e6-<name>-<arm>.json`: a report-only timing run (#822's monster and CI-gate worlds).

    backend/.venv/bin/python bin/debug/score_842_round3.py <capture-dir> \\
        --json analysis/issue842/round3.json
"""

from __future__ import annotations

import argparse
import json
import statistics
from collections import defaultdict
from pathlib import Path

DIAGNOSIS = ("997271774", "649658380", "374464662")
HELD_OUT = ("447276134", "628977350", "295569319")
ICE_AGE_SAVES = ("374464662", "447276134")
ARMS = ("d2", "trend-1k")

LEDGER_TOL = 1.0e-3  # 0.1% of the step's accumulation
SWING_LIMIT = 0.25
GAP_OLD_HELD_OUT, GAP_OLD_DIAGNOSIS, GAP_NEW = 0.10, 0.20, 0.20
RUNTIME_LIMIT = 1.15
DROP_FRACTION = 0.88  # a member drops when it closes below 88% of the unperturbed d2 control
SUBSTEP_MOVE_LIMIT = 0.20
REPLAY_STEPS = 10
LAND_TOL = 1.0e-3  # #789's land-ledger criterion; RFC 0015 D8's meltwater ledgers aren't on main yet
# The old runtime threshold was +15% against the pre-march solve, which main no longer has.
# #789 measured it on the diagnosis saves at 0.25 Myr steps (median s per step), so the old
# check is estimated against those figures: another day and machine load, not like-for-like.
PRE_MARCH_MEDIAN_S = {"997271774": 7.0, "649658380": 6.0, "374464662": 8.0}


def _load(directory: Path) -> dict[str, dict]:
    return {p.stem: json.loads(p.read_text()) for p in sorted(directory.glob("*.json"))}


def _ledger(runs: dict, seed: str, myr: str) -> dict:
    """Gate 1 on the first REPLAY_STEPS steps: old (absolute) and new (or no worse than d2)."""
    d2 = runs[f"{seed}-d2-{myr}"]["per_step"][:REPLAY_STEPS]
    trend = runs[f"{seed}-trend-1k-{myr}"]["per_step"][:REPLAY_STEPS]
    over = {arm: sum(s["ice_closure"] > LEDGER_TOL for s in steps) for arm, steps in (("d2", d2), ("trend-1k", trend))}
    new_fail = [
        t["step"] for t, d in zip(trend, d2)
        if t["ice_closure"] > LEDGER_TOL and t["ice_closure"] > d["ice_closure"]
    ]
    melt_ok = all(s["land_closure"] <= LAND_TOL for s in trend)
    return {
        "steps": len(trend),
        "steps_over_0.1pct": over,
        "worst": {arm: max(s["ice_closure"] for s in steps) for arm, steps in (("d2", d2), ("trend-1k", trend))},
        "new_failing_steps": new_fail,
        "land_ledger_ok": melt_ok,
        "old_pass": over["trend-1k"] == 0 and melt_ok,
        "new_pass": not new_fail and melt_ok,
    }


def _swing(run: dict) -> float:
    """Largest step-to-step change in ice volume after the first 2 Myr."""
    steps, myr = run["per_step"], run["myr_per_step"]
    vols = [s["ice_closing_km3"] for s in steps]
    after = [i for i in range(1, len(vols)) if i * myr >= 2.0]
    return max((abs(vols[i] - vols[i - 1]) / max(vols[i - 1], 1e-9) for i in after), default=0.0)


def _gap(runs: dict, seed: str, arm: str, field: str) -> float | None:
    fine, coarse = runs.get(f"{seed}-{arm}-0.25"), runs.get(f"{seed}-{arm}-1")
    if not fine or not coarse:
        return None
    a, b = fine["per_step"][-1][field], coarse["per_step"][-1][field]
    if a is None or b is None:
        return None
    return abs(a - b) / max(a, 1e-9)  # over the 0.25 Myr run's value, as #765 reported it


def _replays(runs: dict, seeds: tuple[str, ...], held_out: bool) -> dict:
    out = {}
    for seed in seeds:
        rec: dict = {"ledger": {}, "swing": {}, "gap": {}, "agreement_vs_d2": {}}
        for myr in ("0.25", "1"):
            if f"{seed}-trend-1k-{myr}" in runs and f"{seed}-d2-{myr}" in runs:
                rec["ledger"][myr] = _ledger(runs, seed, myr)
                d2_end = runs[f"{seed}-d2-{myr}"]["per_step"][-1]["ice_closing_km3"]
                tr_end = runs[f"{seed}-trend-1k-{myr}"]["per_step"][-1]["ice_closing_km3"]
                rec["agreement_vs_d2"][myr] = (tr_end - d2_end) / max(d2_end, 1e-9)
            for arm in ARMS:
                if f"{seed}-{arm}-{myr}" in runs and seed not in ICE_AGE_SAVES:
                    rec["swing"][f"{arm}-{myr}"] = _swing(runs[f"{seed}-{arm}-{myr}"])
        for arm in ARMS:
            rec["gap"][arm] = {f: _gap(runs, seed, arm, f) for f in ("ice_closing_km3", "ice_area_km2")}
        gap = max((g for g in rec["gap"]["trend-1k"].values() if g is not None), default=None)
        old_gap = GAP_OLD_HELD_OUT if held_out else GAP_OLD_DIAGNOSIS
        swing_ok = all(v <= SWING_LIMIT for k, v in rec["swing"].items() if k.startswith("trend-1k"))
        rec["gates"] = {
            "ledger": {g: all(v[g] for v in rec["ledger"].values()) for g in ("old_pass", "new_pass")},
            "swing": swing_ok,
            "gap": {"value": gap, "old_limit": old_gap, "old_pass": gap is not None and gap < old_gap,
                    "new_pass": gap is not None and gap < GAP_NEW},
        }
        out[seed] = rec
    return out


def _e3(runs: dict) -> dict:
    members: dict = defaultdict(lambda: defaultdict(dict))
    for name, run in runs.items():
        if name.startswith("e3-") and not name.startswith("e3-rerun-"):
            _, step, rest = name.split("-", 2)
            arm, scale = rest.rsplit("-", 1)
            s = run["per_step"][0]
            members[step][arm][float(scale)] = (s["ice_closing_km3"], (s["march"] or {}).get("worst_substep_move"))
    out = {}
    for step, arms in members.items():
        control = arms.get("d2", {}).get(1.0, (None,))[0]
        rec = {"control_km3": control}
        for arm, by_scale in arms.items():
            perturbed = {k: v for k, v in by_scale.items() if k != 1.0}
            drops = None if control is None else sum(v[0] < DROP_FRACTION * control for v in perturbed.values())
            rec[arm] = {"members": len(perturbed), "drops": drops,
                        "closing_km3": [round(v[0], 3) for _, v in sorted(by_scale.items())],
                        "worst_substep_move": max((v[1] or 0.0) for v in by_scale.values())}
        positive = rec.get("d2-10k", {}).get("drops")
        rec["drop_step_confirmed"] = positive is not None and positive >= 8
        trend, d2 = rec.get("trend-1k", {}).get("drops"), rec.get("d2", {}).get("drops")
        rec["pass"] = trend is not None and d2 is not None and trend <= d2 and (d2 > 0 or trend == 0)
        out[step] = rec
    return out


def _worst_substep(runs: dict) -> dict:
    worst = {}
    for name, run in runs.items():
        if "-trend-1k-" in name and not name.startswith(("e3-", "e4-", "e6-")):
            worst[name] = max(((s["march"] or {}).get("worst_substep_move") or 0.0) for s in run["per_step"])
    return {"by_run": worst, "pass": all(v <= SUBSTEP_MOVE_LIMIT for v in worst.values())}


def _e4(runs: dict) -> dict:
    out: dict = defaultdict(dict)
    for name, run in runs.items():
        if name.startswith("e4-"):
            _, seed, sub, myr = name.split("-")
            out[f"{seed}-{myr}"][int(sub.removeprefix("sub"))] = run["per_step"][0]["ice_closing_km3"]
    report = {}
    for key, by_sub in out.items():
        subs = sorted(by_sub, reverse=True)
        report[key] = {
            "closing_km3": {s: round(by_sub[s], 3) for s in subs},
            "successive_change": {f"{a}->{b}": abs(by_sub[b] - by_sub[a]) / max(by_sub[a], 1e-9) for a, b in zip(subs, subs[1:])},
        }
    return report


def _e6(runs: dict) -> dict:
    out = {}
    for seed in DIAGNOSIS:
        for myr in ("0.25", "1"):
            arms = {arm: runs.get(f"{seed}-{arm}-{myr}") for arm in ARMS}
            if all(arms.values()):
                med = {arm: statistics.median(s["wall_s"] for s in r["per_step"]) for arm, r in arms.items()}
                out[f"{seed}-{myr}"] = {"median_wall_s": med, "ratio": med["trend-1k"] / med["d2"]}
    report_only = {name: [s["wall_s"] for s in run["per_step"]] for name, run in runs.items() if name.startswith("e6-")}
    old = {seed: out[f"{seed}-0.25"]["median_wall_s"]["trend-1k"] / base
           for seed, base in PRE_MARCH_MEDIAN_S.items() if f"{seed}-0.25" in out}
    return {"by_save": out, "report_only": report_only,
            "old_estimate_vs_pre_march": old,
            "old_pass": bool(old) and all(v <= RUNTIME_LIMIT for v in old.values()),
            "pass": bool(out) and all(v["ratio"] <= RUNTIME_LIMIT for v in out.values())}


def _verdict(e1, e5, e3, substep, e6, which: str) -> dict:
    replays = list(e1.values()) + list(e5.values())
    conditions = {
        "1_ledger_and_swings": all(r["gates"]["ledger"][f"{which}_pass"] and r["gates"]["swing"] for r in replays),
        "2_step_gap": all(r["gates"]["gap"][f"{which}_pass"] for r in replays),
        # An ensemble gates only on a step its 10 kyr positive control shows is a drop step.
        "3_drops": substep["pass"] and bool(e3) and all(r["pass"] for r in e3.values() if r["drop_step_confirmed"]),
        "4_runtime": e6["pass"] if which == "new" else e6["old_pass"],
    }
    return {"conditions": conditions, "pass": all(conditions.values())}


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    args.add_argument("directory", type=Path)
    args.add_argument("--json", type=Path)
    ns = args.parse_args()
    runs = _load(ns.directory)
    e1, e5 = _replays(runs, DIAGNOSIS, False), _replays(runs, HELD_OUT, True)
    e3, substep, e4, e6 = _e3(runs), _worst_substep(runs), _e4(runs), _e6(runs)
    result = {
        "E1_E2_diagnosis": e1, "E5_held_out": e5, "E3_ensemble": e3, "E3_worst_substep": substep,
        "E4_convergence": e4, "E6_runtime": e6,
        "verdict_new_thresholds": _verdict(e1, e5, e3, substep, e6, "new"),
        "verdict_old_thresholds": _verdict(e1, e5, e3, substep, e6, "old"),
    }
    text = json.dumps(result, indent=1, default=str) + "\n"
    if ns.json:
        ns.json.write_text(text)
    for key in ("verdict_new_thresholds", "verdict_old_thresholds"):
        print(key, json.dumps(result[key]))


if __name__ == "__main__":
    main()
