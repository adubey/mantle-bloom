#!/usr/bin/env python3
"""#842: score the earlier round 3 (2026-11-08/09, commit a81c4ccc) under its own rule.

That run was preregistered in `analysis/issue842/README.md` on branch
`claude/issue-842-adaptive-march-default` and never scored (see this branch's
`analysis/issue842/preregistration.md`, amendment 5). This script applies that rule as
written. Where the text leaves a choice open, it was settled here before any of the
run's numbers were read:

- Rule 2 asks for median *march* time, but `replay_789_march.py`'s E1/E5 records hold only
  whole-step time (`wall_s`). Whole-step time stands in. That is conservative for the flip,
  since non-march work dilutes the ratio toward 1.
- Rule 2 pairs steps by index: a step counts when `d2` marched more than 30 substeps on it.
- Rule 3's E4 exception applies to seed 997271774 only, the one E1 save with an E4 run:
  `d2`'s |250 − 125 yr| closing-ice difference on E4's step 3, over the 250 yr value.
- Rule 4's gap is over the 0.25 Myr run's value (as #765 reported), for volume and area,
  and the larger gates. The gap over the 1 Myr value is reported beside it.
- Rule 5 counts flips over all three ensemble worlds together. Per-world counts are
  reported too. The no-move-over-20% part applies to `trend-adapt-1k` members (the rule is
  "all of these must hold for `trend-adapt-1k`"). The other arms' moves are reported.

    backend/.venv/bin/python bin/debug/score_842_round3_nov8.py <issue-842 capture dir> \\
        --json analysis/issue842/round3_nov8.json
"""

from __future__ import annotations

import argparse
import json
import statistics
from pathlib import Path

E1 = {"s997": "997271774", "s649": "649658380", "s374": "374464662"}
E5 = {"h447": "447276134", "h628": "628977350", "h295": "295569319"}
TREND = "trend-adapt-1k"


def _steps(directory: Path, name: str) -> list[dict]:
    return json.loads((directory / f"{name}.json").read_text())["per_step"]


def _rule1(d2: list[dict], tr: list[dict]) -> dict:
    fails = [t["step"] for t, d in zip(tr, d2) if t["ice_closure"] > max(1e-3, d["ice_closure"] + 1e-4)]
    land = [t["step"] for t in tr if t["land_closure"] > 1e-3]
    return {"glacier_failing_steps": fails, "land_failing_steps": land,
            "worst_glacier": {"d2": max(s["ice_closure"] for s in d2), TREND: max(s["ice_closure"] for s in tr)},
            "pass": not fails and not land}


def _rule2(d2: list[dict], tr: list[dict]) -> dict:
    slow = [i for i, d in enumerate(d2) if ((d["march"] or {}).get("substeps") or 0) > 30]
    if slow:
        ratio = statistics.median(tr[i]["wall_s"] for i in slow) / statistics.median(d2[i]["wall_s"] for i in slow)
        limit = 0.6
    else:
        ratio = statistics.median(s["wall_s"] for s in tr) / statistics.median(s["wall_s"] for s in d2)
        limit = 1.1
    return {"slow_steps": slow, "whole_step_time_ratio": ratio, "limit": limit, "pass": ratio <= limit}


def _e4_self_difference(directory: Path) -> float | None:
    try:
        a = _steps(directory, "e4-s997-d2-250")[3]["ice_closing_km3"]
        b = _steps(directory, "e4-s997-d2-125")[3]["ice_closing_km3"]
    except (FileNotFoundError, IndexError):
        return None
    return abs(a - b) / max(a, 1e-9)


def _rule3(d2: list[dict], tr: list[dict], seed: str, self_diff: float | None) -> dict:
    diff = abs(tr[-1]["ice_closing_km3"] - d2[-1]["ice_closing_km3"]) / max(d2[-1]["ice_closing_km3"], 1e-9)
    limit = 0.02
    if seed == "997271774" and self_diff is not None and self_diff > limit:
        limit = self_diff
    return {"ice_vs_d2": diff, "limit": limit, "pass": diff <= limit}


def _gap(fine: list[dict], coarse: list[dict], field: str) -> tuple[float, float]:
    a, b = fine[-1][field], coarse[-1][field]
    return abs(a - b) / max(a, 1e-9), abs(a - b) / max(b, 1e-9)


def _rule4(directory: Path) -> dict:
    out = {}
    for arm in ("d2", TREND):
        fine = _steps(directory, f"e2-s997-{arm}-025")
        coarse = _steps(directory, f"e2-s997-{arm}-1")
        out[arm] = {f: _gap(fine, coarse, f) for f in ("ice_closing_km3", "ice_area_km2")}
    gate = {arm: max(v[0] for v in out[arm].values()) for arm in out}
    return {"gaps_over_fine_and_coarse": out, "gap": gate,
            "pass": gate[TREND] <= 0.10 or gate[TREND] <= gate["d2"]}


def _rule5(directory: Path) -> dict:
    worlds, flips, moves = {}, {"d2": 0, TREND: 0, "trend-adapt": 0}, {}
    for world in ("s997", "s649", "gen30"):
        worlds[world] = {}
        for arm in flips:
            path = directory / f"e3-{world}-{arm}.json"
            if not path.exists():
                continue
            members = json.loads(path.read_text())["ensemble"]
            base = next(m for m in members if m["member"] == 0)["ice_closing_km3"]
            flipped = [m["member"] for m in members
                       if abs(m["ice_closing_km3"] - base) / max(base, 1e-9) > 0.02]
            worst = max(((m["march"] or {}).get("max_substep_change") or 0.0) for m in members)
            worlds[world][arm] = {"members": len(members), "flipped": flipped, "worst_substep_move": worst}
            flips[arm] += len(flipped)
            moves[arm] = max(moves.get(arm, 0.0), worst)
    return {"by_world": worlds, "flips": flips, "worst_substep_move": moves,
            "pass": moves.get(TREND, 1.0) <= 0.20 and flips[TREND] <= flips["d2"] + 2}


def main() -> None:
    args = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    args.add_argument("directory", type=Path)
    args.add_argument("--json", type=Path)
    ns = args.parse_args()
    self_diff = _e4_self_difference(ns.directory)
    saves = {}
    for prefix, labels in (("e1", E1), ("e5", E5)):
        for label, seed in labels.items():
            d2 = _steps(ns.directory, f"{prefix}-{label}-d2")
            tr = _steps(ns.directory, f"{prefix}-{label}-{TREND}")
            saves[seed] = {"set": prefix, "rule1": _rule1(d2, tr), "rule2": _rule2(d2, tr),
                           "rule3": _rule3(d2, tr, seed, self_diff)}
    rule4, rule5 = _rule4(ns.directory), _rule5(ns.directory)
    conditions = {f"rule{k}": all(s[f"rule{k}"]["pass"] for s in saves.values()) for k in (1, 2, 3)}
    conditions["rule4"], conditions["rule5"] = rule4["pass"], rule5["pass"]
    averaged = _amendment6(ns.directory, saves, rule4, rule5)
    result = {"tree": "a81c4ccc", "e4_d2_self_difference_997_step3": self_diff, "saves": saves,
              "rule4_step_gap_997": rule4, "rule5_ensembles": rule5,
              "verdict_as_written": {"conditions": conditions, "pass": all(conditions.values())},
              "amendment6_averaged": averaged}
    text = json.dumps(result, indent=1) + "\n"
    if ns.json:
        ns.json.write_text(text)
    print("as written", json.dumps(result["verdict_as_written"]))
    print("amendment 6", json.dumps(averaged["verdict"]))


def _amendment6(directory: Path, saves: dict, rule4: dict, rule5: dict) -> dict:
    """The owner's 2026-11-13 averaged rules, set after this run's results were seen.

    Rule 1: mean glacier closure over all 60 checked steps <= max(1e-3, d2's mean + 1e-4); the
    land ledger stays per step. Rule 2: the mean of the six per-save ratios <= 0.6. Rule 5
    split into 5a (no substep move > 20%) and 5b (flips <= d2's + 2). Rules 3 and 4 unchanged.
    """
    closures: dict = {"d2": [], TREND: []}
    per_save = {}
    for prefix, labels in (("e1", E1), ("e5", E5)):
        for label, seed in labels.items():
            for arm in closures:
                closures[arm] += [s["ice_closure"] for s in _steps(directory, f"{prefix}-{label}-{arm}")]
            per_save[seed] = {arm: statistics.mean(closures[arm][-10:]) for arm in closures}
    mean = {arm: statistics.mean(v) for arm, v in closures.items()}
    limit = max(1e-3, mean["d2"] + 1e-4)
    land_ok = all(not s["rule1"]["land_failing_steps"] for s in saves.values())
    ratio = statistics.mean(s["rule2"]["whole_step_time_ratio"] for s in saves.values())
    conditions = {
        "rule1_mean": mean[TREND] <= limit and land_ok,
        "rule2_mean": ratio <= 0.6,
        "rule3": all(s["rule3"]["pass"] for s in saves.values()),
        "rule4": rule4["pass"],
        "rule5a_substep_move": rule5["worst_substep_move"].get(TREND, 1.0) <= 0.20,
        "rule5b_flips": rule5["flips"][TREND] <= rule5["flips"]["d2"] + 2,
    }
    return {"rule1": {"mean_closure": mean, "limit": limit, "per_save_mean_report_only": per_save},
            "rule2": {"mean_ratio": ratio, "limit": 0.6},
            "verdict": {"conditions": conditions, "pass": all(conditions.values())}}


if __name__ == "__main__":
    main()
