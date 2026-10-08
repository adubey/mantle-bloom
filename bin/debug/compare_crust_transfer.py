#!/usr/bin/env python3
"""Issue #320: compare `attribute_hm_ratchet.py` replays with and without crust transfer.

Each replay JSON is one run of the same save and span. Reports, per run, at the final
checkpoint: land, continental and emergent continental area, land hypsometry, continental-
material inventory, collision loss of continental material (total and per Myr, by account),
the continental Hm-cap fraction, runtime, and -- where the engine has it -- crust-transfer
fronts, volumes and per-front/per-step cost.

    python3 bin/debug/compare_crust_transfer.py \\
        main=/path/to/base.json transfer=/path/to/transfer.json \\
        --baseline314 analysis/issue314/cascade-seed997271774.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

LOSS_ACCOUNTS = (
    "delaminated_lower_crust_m3",
    "collision_subducted_m3",
    "collision_lower_crust_subducted_m3",
    "deeply_subducted_m3",
)


def _delta(first: dict, last: dict, key: str) -> float:
    return float(last.get(key, 0.0)) - float(first.get(key, 0.0))


def summarize(path: Path) -> dict:
    run = json.loads(path.read_text())
    first, last = run["checkpoints"][0], run["checkpoints"][-1]
    m0, m1 = first["metrics"], last["metrics"]
    span_myr = m1["elapsed_myr"] - m0["elapsed_myr"]
    out = {
        "commit": run.get("commit", "?")[:10],
        "span_myr": span_myr,
        "land_area_Mkm2": m1["land_area_km2"] / 1e6,
        "continental_area_Mkm2": m1["continental_area_km2"] / 1e6,
        "emergent_continental_Mkm2": m1["emergent_continental_area_km2"] / 1e6,
        "land_elevation_p10_p50_p90_p99_m": [round(m1["land_elevation_m"][f"p{q}"]) for q in (10, 50, 90, 99)],
        "continental_hc_p50_p90_p99_km": [round(m1["continental_hc_km"][f"p{q}"], 2) for q in (50, 90, 99)],
        "continental_hm_cap_pct": 100.0 * m1["continental"]["hm_cap_fraction"],
        "continental_hm_cap_change_pp": 100.0 * (m1["continental"]["hm_cap_fraction"] - m0["continental"]["hm_cap_fraction"]),
        "runtime_s": run.get("runtime_seconds"),
        "peak_memory_mb": run.get("peak_memory_mb"),
    }
    if "crust_accounts" in first and "crust_accounts" in last:
        a0, a1 = first["crust_accounts"], last["crust_accounts"]
        live0 = a0["continental_material_live_km3"]
        out["continental_material_Mkm3"] = a1["continental_material_live_km3"] / 1e6
        out["continental_material_change_pct"] = 100.0 * (a1["continental_material_live_km3"] - live0) / live0
        losses = {k: _delta(a0["continental_ledger_km3"], a1["continental_ledger_km3"], k) for k in LOSS_ACCOUNTS}
        collision = sum(v for k, v in losses.items() if k != "deeply_subducted_m3")
        out["collision_loss_km3_by_account"] = {k: round(v) for k, v in losses.items()}
        out["collision_loss_pct_inventory_per_myr"] = 100.0 * collision / live0 / span_myr
        budget0, budget1 = a0["orogeny_budget_km3"], a1["orogeny_budget_km3"]
        donated = _delta(budget0, budget1, "suture_donated_m3")
        out["suture_donated_pct_inventory_per_myr"] = 100.0 * donated / live0 / span_myr
        out["suture_budget_km3"] = {
            k: round(_delta(budget0, budget1, k))
            for k in (
                "suture_donated_m3", "suture_scraped_m3", "suture_underthrust_m3", "suture_underthrust_placed_m3",
                "suture_lower_crust_subducted_m3", "upper_plate_placed_m3", "overrider_placed_m3",
                "no_outlet_subducted_m3", "delamination_completed_m3",
            )
            if k in budget1
        }
        out["ledger_errors_km3"] = {
            k: a1.get(k) for k in ("continental_balance_error_km3", "craton_balance_error_km3", "mobile_cover_balance_error_km3")
        }
        transfer = a1.get("suture_transfer")
        if transfer:
            by_step = transfer.get("by_step", {})
            per_step = sorted(row["seconds"] for row in by_step.values())
            out["transfer"] = {
                "fronts": transfer.get("fronts"),
                "seconds": transfer.get("seconds"),
                "mean_front_ms": 1e3 * transfer["seconds"] / max(transfer["fronts"], 1),
                "max_front_ms": 1e3 * transfer.get("max_front_seconds", 0.0),
                "steps_with_fronts": len(by_step),
                "median_step_ms": 1e3 * per_step[len(per_step) // 2] if per_step else 0.0,
                "max_step_ms": 1e3 * per_step[-1] if per_step else 0.0,
            }
    steps = run.get("steps", [])
    if steps:
        newly = 0.0
        for step in steps:
            row = step.get("newly_capped_by_phase_km2", {}).get("boundary_retreat", {}).get("continental_node", {})
            newly += row.get("newly_capped_area_km2", 0.0)
        out["boundary_retreat_newly_capped_continental_km2"] = newly
        times = sorted(step["runtime_seconds"] for step in steps)
        out["median_step_s"] = times[len(times) // 2]
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="label=path to an attribute_hm_ratchet.py JSON")
    ap.add_argument("--baseline314", type=Path, help="issue #314's checked-in replay (metrics only)")
    args = ap.parse_args()
    rows = {}
    if args.baseline314:
        rows["#314"] = summarize(args.baseline314)
    for spec in args.runs:
        label, _, path = spec.partition("=")
        rows[label] = summarize(Path(path))
    print(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
