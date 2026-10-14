#!/usr/bin/env python3
"""Compact re-baseline numbers for RFCs 0001-0007 (issue #403).

Combines, per seed, `compare_cratons.py`'s checkpoint rows (topology events and plate
counts) with `ledger_history.summarize` on the run's final save (`compare_cratons.py
--save-dir`), and adds what the RFCs quote that neither prints:

- craton sinks per Gyr as a share of the volume at the window start (ledger_history's table);
- merge-sized `delaminated_m3` steps: how many samples (one per step) booked more than 0.1%,
  0.5% and 1% of the live craton, and the largest (RFC 0001; RFC 0008's gate found fusions
  carry this account, so a sample's jump is all the merges in that step);
- `formed_m3` per 500 Myr window, and the craton share of genuine continental area at each
  500 Myr checkpoint (RFC 0004);
- splits, merges and failed rifts per Gyr, and plate count at each checkpoint;
- continental-material fluxes per window (RFCs 0005-0007), with each sink account's rate and
  the live inventory at each window's ends.

Usage:
    backend/.venv/bin/python bin/debug/rebaseline_summary.py --rows /tmp/runs/*.jsonl \
        --saves /tmp/runs/*.mbworld --baseline "origin/main 2562f63" --out analysis/issue403/baseline.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from ledger_history import CONTINENTAL_SINKS, CRATON_SINKS, FLUX_WINDOWS_MYR, summarize  # noqa: E402

from app import persistence  # noqa: E402

JUMP_SHARES = (0.001, 0.005, 0.01)
FORMATION_WINDOW_MYR = 500.0


def _round(value, digits: int = 4):
    if isinstance(value, float):
        return round(value, digits)
    if isinstance(value, dict):
        return {key: _round(item, digits) for key, item in value.items()}
    if isinstance(value, list):
        return [_round(item, digits) for item in value]
    return value


def history_figures(history: list[dict]) -> dict:
    t = np.array([row["elapsed_years"] for row in history]) / 1e6
    volume = np.array([row["craton_volume_km3"] for row in history])
    ledger = {key: np.array([row["craton_ledger_km3"].get(key, 0.0) for row in history]) for key in ("delaminated_m3", "formed_m3")}
    jump = np.diff(ledger["delaminated_m3"])
    share = np.divide(jump, volume[:-1], out=np.zeros_like(jump), where=volume[:-1] > 0.0)
    largest = int(np.argmax(share)) if len(share) else 0
    windows = []
    for lo in np.arange(0.0, t[-1], FORMATION_WINDOW_MYR):
        i = int(np.searchsorted(t, lo))
        j = min(int(np.searchsorted(t, lo + FORMATION_WINDOW_MYR)), len(t) - 1)
        windows.append({"window_myr": [float(t[i]), float(t[j])], "formed_mkm3": (ledger["formed_m3"][j] - ledger["formed_m3"][i]) / 1e6})
    shares = []
    for at in np.arange(0.0, t[-1] + 1.0, FORMATION_WINDOW_MYR):
        k = min(int(np.searchsorted(t, at)), len(t) - 1)
        shares.append({"myr": float(t[k]), "craton_share_of_genuine_continental": history[k]["craton_continental_fraction"]})
    material = [row["continental_material_ledger"] for row in history]
    sinks = []
    for lo, hi in FLUX_WINDOWS_MYR:
        i, j = min(int(np.searchsorted(t, lo)), len(t) - 1), min(int(np.searchsorted(t, hi)), len(t) - 1)
        years = (t[j] - t[i]) * 1e6
        if years <= 0.0:
            continue
        sinks.append({
            "window_myr": [float(t[i]), float(t[j])],
            "inventory_mkm3": [material[i]["surface_continental_derived_m3"] / 1e15, material[j]["surface_continental_derived_m3"] / 1e15],
            **{f"{key.removesuffix('_m3')}_km3_per_yr": (material[j][key] - material[i][key]) / 1e9 / years for key in CONTINENTAL_SINKS},
        })
    return {
        "samples": len(history),
        "initial_continental_mkm3": material[0]["initial_continental_m3"] / 1e15,
        "continental_sinks_per_window": sinks,
        "sample_spacing_myr": float(np.median(np.diff(t))) if len(t) > 1 else 0.0,
        "delaminated_samples_over_share_of_live": {f"{s:.1%}": int(np.count_nonzero(share > s)) for s in JUMP_SHARES},
        "delaminated_largest_sample": {"myr": float(t[largest + 1]), "mkm3": float(jump[largest]) / 1e6, "share_of_live": float(share[largest])},
        "formed_per_window": windows,
        "craton_share_at_checkpoints": shares,
    }


def row_figures(rows: list[dict]) -> dict:
    rows = sorted(rows, key=lambda row: row["myr"])
    last = rows[-1]
    gyr = last["myr"] / 1000.0
    by_gyr = defaultdict(dict)
    for row in rows:
        if round(row["myr"]) % 1000 == 0:
            by_gyr[round(row["myr"]) // 1000] = row["events"]
    events_per_gyr = []
    for g in sorted(by_gyr)[1:]:
        before, after = by_gyr[g - 1], by_gyr[g]
        events_per_gyr.append({"gyr": g - 1, **{key: after.get(key, 0) - before.get(key, 0) for key in ("split", "merged", "failed_rift")}})
    plates = [row["plate_count"] for row in rows]
    return {
        "myr": last["myr"],
        "events_total": last["events"],
        "events_per_gyr": events_per_gyr or [{"gyr": 0, **{key: last["events"].get(key, 0) / gyr for key in ("split", "merged", "failed_rift")}}],
        "plate_count": {"start": plates[0], "end": plates[-1], "min": min(plates), "max": max(plates), "mean": float(np.mean(plates))},
        "plate_count_by_checkpoint": {str(round(row["myr"])): row["plate_count"] for row in rows},
        "craton_balance_error_km3": last["craton_balance_error_km3"],
        "genuine_continental_area_mkm2": {str(round(row["myr"])): row["genuine_continental_area_km2"] / 1e6 for row in rows},
        "land_fraction_node": {str(round(row["myr"])): row["land_fraction_node"] for row in rows},
        "rift_cuts": last.get("rift_cuts", {}),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rows", type=Path, nargs="+", required=True, help="compare_cratons.py --out files")
    parser.add_argument("--saves", type=Path, nargs="+", required=True, help="compare_cratons.py --save-dir worlds, named seed<S>-...")
    parser.add_argument("--baseline", required=True, help='the commit measured, e.g. "origin/main 2562f63"')
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = defaultdict(list)
    for path in args.rows:
        for line in path.read_text().splitlines():
            if line.strip():
                row = json.loads(line)
                rows[row["seed"]].append(row)
    seeds = {}
    for save in args.saves:
        seed = int(save.name.removeprefix("seed").split("-")[0])
        history = persistence.load_world_bytes(save.read_bytes()).stats_history
        ledger = summarize(history)
        seeds[str(seed)] = {
            **row_figures(rows[seed]),
            **history_figures(history),
            "craton_loss_per_gyr": ledger["craton_loss_per_gyr"],
            "craton_checkpoints": ledger["checkpoints"],
            "craton_largest_jumps": {key: ledger["craton_largest_jumps"][key] for key in (*CRATON_SINKS, "formed_m3")},
            "continental_fluxes": ledger["continental_fluxes"],
        }
        print(f"seed {seed}: {len(history)} samples", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(_round({"baseline": args.baseline, "seeds": seeds}), indent=1) + "\n")


if __name__ == "__main__":
    main()
