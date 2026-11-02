#!/usr/bin/env python3
"""RFC 0016 G1b: score paired replay arms against #680's owner-selected materiality rule.

Reads `replay_rfc0016_river_water.py` outputs named `<window>-<arm>.json`, one directory on the
100 kyr production schedule and one with the registered 10 kyr connectivity trace. A window is
material for an arm when its registered event's `|ΔV|` changes by at least 10% from the control,
or when the fixed anchor's wet/dry class differs from the control at two consecutive 10 kyr
outputs. An arm is material when that holds in at least two of the three windows. The event `ΔV`
is the target body's volume against its matched volume at the previous step (#588's overlap rule),
because A's anchor is dry the step before its event although its lake is present (#751).

Example:

    python3 bin/debug/score_rfc0016_g1b.py --production <captures>/step100kyr-matched \\
        --trace <captures>/step10kyr --out analysis/rfc0016/g1b-summary.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

WINDOWS = ("a", "d", "e")
# The output whose change from the previous one is the registered event: the reversal two steps
# after each checkpoint, as #680 placed it.
EVENT_STEP = {"a": 2, "d": 2, "e": 2}
COUNTED = (
    "central", "capacity-150", "capacity-1000", "depletion-0.3", "depletion-0.7",
    "demand-0.85", "demand-1.15", "split-low", "split-high",
)
REPORTED = ("partition-match", "coupling-off", "sham-up", "sham-down")
EVENT_THRESHOLD = 0.10
CONSECUTIVE_OUTPUTS = 2


def event_delta_km3(steps: list[dict], window: str) -> float | None:
    """The registered event's `ΔV`, or None when the anchor's body has no match at the step before."""
    target = steps[EVENT_STEP[window] - 1]["target"]
    previous = target.get("matched_previous_km3")
    if not target["anchor_wet"] or previous is None:
        return None
    return target["volume_km3"] - previous


def longest_mismatch(control: list[dict], arm: list[dict]) -> tuple[int, list[int]]:
    """The longest run of consecutive outputs whose anchor class differs, and where they differ."""
    mismatched = [
        index + 1 for index, (c, a) in enumerate(zip(control, arm, strict=True))
        if c["target"]["anchor_wet"] != a["target"]["anchor_wet"]
    ]
    longest = run = 0
    for index, output in enumerate(mismatched):
        run = run + 1 if index and output == mismatched[index - 1] + 1 else 1
        longest = max(longest, run)
    return longest, mismatched


def ledger(steps: list[dict]) -> dict:
    stores = [s["land_store"] for s in steps if "max_identity_mm" in s.get("land_store", {})]
    return {
        "max_normalized_land_residual": max(s["normalized_land_residual"] for s in steps),
        "max_normalized_lake_residual": max(s["normalized_lake_residual"] for s in steps),
        "stranded_overflow_m3": sum(s["land_rain_only"]["stranded_overflow_m3"] for s in steps),
        "max_store_identity_mm": max((s["max_identity_mm"] for s in stores), default=None),
    }


def describe(change: float | None) -> str:
    return "lost" if change is None else f"{100.0 * change:+.2f}%"


def score_window(window: str, arm: str, production: Path, trace: Path) -> dict:
    load = lambda folder, name: json.loads((folder / f"{window}-{name}.json").read_text())["steps"]  # noqa: E731
    control, treated = load(production, "control"), load(production, arm)
    control_dv, treated_dv = event_delta_km3(control, window), event_delta_km3(treated, window)
    change = None if treated_dv is None else abs(treated_dv) / abs(control_dv) - 1.0
    longest, mismatched = longest_mismatch(load(trace, "control"), load(trace, arm))
    stores = [s["land_store"] for s in treated if "excess_ratio" in s.get("land_store", {})]
    return {
        "control_event_dv_km3": control_dv,
        "event_dv_km3": treated_dv,
        "event_change": change,
        # The rule scores class changes on the 10 kyr trace, so a production event the treatment
        # loses is flagged here, not scored.
        "event_lost": change is None,
        "event_material": change is not None and abs(change) >= EVENT_THRESHOLD,
        "target_volume_km3": [s["target"].get("volume_km3") for s in treated],
        "trace_mismatched_outputs": mismatched,
        "trace_longest_mismatch": longest,
        "trace_material": longest >= CONSECUTIVE_OUTPUTS,
        "excess_ratio": [s["excess_ratio"] for s in stores] or None,
        "ledger": ledger(treated),
        "trace_ledger": ledger(load(trace, arm)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--production", type=Path, required=True, help="100 kyr arm outputs")
    ap.add_argument("--trace", type=Path, required=True, help="10 kyr arm outputs")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    arms = {}
    for arm in COUNTED + REPORTED:
        windows = {w: score_window(w, arm, args.production, args.trace) for w in WINDOWS}
        material = [w for w, s in windows.items() if s["event_material"] or s["trace_material"]]
        arms[arm] = {"counted": arm in COUNTED, "material_windows": material, "material": len(material) >= 2,
                     "windows": windows}
        print(arm, "material" if len(material) >= 2 else "not material", material, " ".join(
            f"{w}:{describe(s['event_change'])}/run{s['trace_longest_mismatch']}" for w, s in windows.items()
        ))
    controls = {w: ledger(json.loads((args.production / f"{w}-control.json").read_text())["steps"]) for w in WINDOWS}
    summary = {
        "rule": {"event_threshold": EVENT_THRESHOLD, "consecutive_outputs": CONSECUTIVE_OUTPUTS,
                 "windows_needed": 2, "event_step": EVENT_STEP},
        "control_ledger": controls,
        "package_material": any(arms[a]["material"] for a in COUNTED),
        "arms": arms,
    }
    args.out.write_text(json.dumps(summary, indent=1) + "\n")


if __name__ == "__main__":
    main()
