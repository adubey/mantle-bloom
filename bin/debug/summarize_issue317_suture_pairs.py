#!/usr/bin/env python3
"""Render the recorded issue #317 suture Hm totals by donor/neighbor pair."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "analysis/issue317/cascade-after-316-seed997271774.json"
DEFAULT_OUTPUT = ROOT / "analysis/issue317/suture-pair-breakdown-after-316.md"


def _millions(value: float) -> str:
    return f"{value / 1e6:,.3f}"


def render(data: dict) -> str:
    checkpoint = data["checkpoints"][-1]
    budget = checkpoint["hm_suture_budget_km3"]
    by_pair = budget["by_pair"]
    bilateral: dict[tuple[int, int], dict[str, dict]] = {}
    multi_candidate: list[tuple[str, dict]] = []

    for key, row in by_pair.items():
        donor_text, receivers_text = key.split("->", 1)
        receivers = [int(value) for value in receivers_text.split("|")]
        donor = int(donor_text)
        if len(receivers) != 1:
            multi_candidate.append((key, row))
            continue
        receiver = receivers[0]
        pair = tuple(sorted((donor, receiver)))
        direction = f"{donor}->{receiver}"
        bilateral.setdefault(pair, {})[direction] = row

    recorded_fronts = sum(row["fronts"] for row in by_pair.values())
    if recorded_fronts != budget["fronts"]:
        raise ValueError(f"pair groups count {recorded_fronts} fronts; total is {budget['fronts']}")
    for field in ("donor_hm_km3", "placed_hm_km3", "unplaced_hm_km3"):
        grouped = sum(row[field] for row in by_pair.values())
        if not math.isclose(grouped, budget[field], rel_tol=1e-12, abs_tol=1e-3):
            raise ValueError(f"pair groups sum to {grouped} {field}; total is {budget[field]}")

    lines = [
        "# Issue #317 suture Hm by plate pair",
        "",
        f"Replay commit `{data['commit']}`, final checkpoint at {checkpoint['metrics']['elapsed_myr']:.1f} Myr.",
        "All volumes are exact-area Hm volumes, in million km³ (M km³).",
        "",
        f"The final ledger contains {budget['fronts']:,} connected donor fronts across "
        f"{len(bilateral)} single-neighbor bilateral pairs and {len(multi_candidate)} "
        "multi-neighbor candidate groups. Each direction is donor plate → neighbor plate.",
        "",
        "For a single-neighbor pair, both donor directions are shown separately. A front "
        "whose adjacent cells identify multiple candidate plates remains grouped under its "
        "full candidate set; its volume is not divided among candidates. The saved ledger "
        "aggregates fronts by this key, so it reports front counts and grouped volumes rather "
        "than individual spatial front IDs.",
        "",
        "## Single-neighbor bilateral pairs",
        "",
        "| Donor → neighbor | Fronts | Donor Hm | Placed on survivors | Unplaced / delaminated |",
        "|---|---:|---:|---:|---:|",
    ]
    directed_rows = []
    for (plate_a, plate_b), directions in bilateral.items():
        for donor, receiver in ((plate_a, plate_b), (plate_b, plate_a)):
            key = f"{donor}->{receiver}"
            row = directions.get(key)
            if row is not None:
                directed_rows.append((row["donor_hm_km3"], key, row))
    for _, key, row in sorted(directed_rows, reverse=True):
        lines.append(
            f"| {key} | {row['fronts']:,} | {_millions(row['donor_hm_km3'])} | "
            f"{_millions(row['placed_hm_km3'])} | {_millions(row['unplaced_hm_km3'])} |"
        )

    lines.extend(
        [
            "",
            "## Multi-neighbor candidate groups",
            "",
            "| Donor → candidate neighbors | Fronts | Donor Hm | Placed on survivors | Unplaced / delaminated |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for key, row in sorted(multi_candidate, key=lambda item: item[1]["donor_hm_km3"], reverse=True):
        lines.append(
            f"| {key} | {row['fronts']:,} | {_millions(row['donor_hm_km3'])} | "
            f"{_millions(row['placed_hm_km3'])} | {_millions(row['unplaced_hm_km3'])} |"
        )

    lines.extend(
        [
            "",
            "## Reconciliation",
            "",
            "| Measure | By pair/candidate groups | Recorded total |",
            "|---|---:|---:|",
            f"| Fronts | {recorded_fronts:,} | {budget['fronts']:,} |",
            f"| Donor Hm (M km³) | {_millions(sum(r['donor_hm_km3'] for r in by_pair.values()))} | {_millions(budget['donor_hm_km3'])} |",
            f"| Placed Hm (M km³) | {_millions(sum(r['placed_hm_km3'] for r in by_pair.values()))} | {_millions(budget['placed_hm_km3'])} |",
            f"| Unplaced / delaminated Hm (M km³) | {_millions(sum(r['unplaced_hm_km3'] for r in by_pair.values()))} | {_millions(budget['unplaced_hm_km3'])} |",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    report = render(data)
    args.output.write_text(report)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
