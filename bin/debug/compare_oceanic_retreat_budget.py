#!/usr/bin/env python3
"""Before/after volume comparison for GitHub issue #216's oceanic self-plate retreat budget fix
(see OCEANIC_SELF_RETREAT_BUDGET_MULTIPLIER in lithosphere_plate.py). That fix caps an oceanic
self-plate's own ordinary subduction retreat -- both end retreat and the interior-subduction
carve -- by however much this same step's decompression-melting creation is adding across the
whole plate, the same shape of cap issue #177 already applied to a continental plate's
oceanic-override retreat against arc-magmatic creation. Before the fix, oceanic self-plate
retreat was bounded only by geometric/count caps (n_distance_cap/max_extend_nodes), with
nothing tying it to any compensating creation -- the issue's own phase-budget measurements
found this uncompensated loss was the dominant net Hc/Hm driver, once the self-canceling
stretch/regularize bookkeeping pair was excluded.

Run this script twice against two different git states of lithosphere_plate.py (before the fix
and after), each with a distinct --label, then pass both --out files to --compare to print the
side-by-side table. Reuses SWEEP_SEEDS (the same 10 seeds issue #171's own sweep harness uses,
drawn from real reported-problem saves) and STEP_YEARS/NODE_DENSITY conventions from
sweep_lib.py so this fits the project's existing sweep infrastructure rather than inventing a
new one, even though it isn't a parameter-multiplier sweep (there's no tunable knob for this --
it's a structural code change) and so doesn't go through run_sweep.py itself.

Usage (run from the repo root; each invocation reads whatever git state HEAD currently is):
    backend/.venv/bin/python bin/debug/compare_oceanic_retreat_budget.py --label before \
        --out bin/debug/results/oceanic_retreat_budget_before.jsonl
    # ... git checkout the fix, or apply it ...
    backend/.venv/bin/python bin/debug/compare_oceanic_retreat_budget.py --label after \
        --out bin/debug/results/oceanic_retreat_budget_after.jsonl
    backend/.venv/bin/python bin/debug/compare_oceanic_retreat_budget.py --compare \
        bin/debug/results/oceanic_retreat_budget_before.jsonl \
        bin/debug/results/oceanic_retreat_budget_after.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import SWEEP_SEEDS  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import stats  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

NODE_DENSITY = 1.0
STEP_YEARS = 10_000_000
CHECKPOINTS_MY = (10, 30, 60)


def run_one(seed: int, label: str) -> list[dict]:
    world = generate_world(seed=seed, node_density=NODE_DENSITY)

    rows = []
    years_done = 0.0
    for checkpoint_my in CHECKPOINTS_MY:
        target = checkpoint_my * 1_000_000
        while years_done < target:
            step = min(STEP_YEARS, target - years_done)
            step_world(world, years=step)
            years_done += step

        total_hc = 0.0
        total_hm = 0.0
        total_nodes = 0
        max_plate_nodes = 0
        oceanic_nodes = 0
        oceanic_hc = 0.0
        oceanic_hm = 0.0
        for plate in world.plates:
            if not plate.lines:
                continue
            hc = plate.collect("crustal_thickness_m")
            hm = plate.collect("mantle_lithosphere_thickness_m")
            n = len(hc)
            total_hc += float(hc.sum())
            total_hm += float(hm.sum())
            total_nodes += n
            max_plate_nodes = max(max_plate_nodes, n)
            if plate.crust_type == "oceanic":
                oceanic_nodes += n
                oceanic_hc += float(hc.sum())
                oceanic_hm += float(hm.sum())

        land_fraction = stats.compute_stats(world)["land_fraction"]

        rows.append({
            "label": label,
            "seed": seed,
            "years_my": years_done / 1_000_000,
            "total_nodes": total_nodes,
            "max_plate_nodes": max_plate_nodes,
            "num_plates": len(world.plates),
            "mean_hc_m": total_hc / total_nodes if total_nodes else None,
            "mean_hm_m": total_hm / total_nodes if total_nodes else None,
            "sum_hc_m": total_hc,
            "sum_hm_m": total_hm,
            "oceanic_mean_hc_m": oceanic_hc / oceanic_nodes if oceanic_nodes else None,
            "oceanic_mean_hm_m": oceanic_hm / oceanic_nodes if oceanic_nodes else None,
            "land_fraction": land_fraction,
        })

    return rows


def run_sweep(label: str, out_path: Path) -> None:
    print(f"{len(SWEEP_SEEDS)} seed(s) to run, label={label!r}, checkpoints at {CHECKPOINTS_MY} My.", file=sys.stderr)
    all_rows: list[dict] = []
    t0 = time.perf_counter()
    for n, seed in enumerate(SWEEP_SEEDS, start=1):
        try:
            rows = run_one(seed, label)
        except Exception as exc:  # noqa: BLE001 -- keep the sweep going past one bad seed
            print(f"FAILED seed={seed}: {exc!r}", file=sys.stderr)
            continue
        all_rows.extend(rows)
        elapsed = time.perf_counter() - t0
        print(f"[{n}/{len(SWEEP_SEEDS)}] seed={seed} done ({elapsed:.0f}s elapsed)", file=sys.stderr)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        for row in all_rows:
            f.write(json.dumps(row) + "\n")
    print(f"\nWrote {len(all_rows)} row(s) to {out_path}", file=sys.stderr)


def _load(path: Path) -> list[dict]:
    with path.open() as f:
        return [json.loads(line) for line in f if line.strip()]


def compare(before_path: Path, after_path: Path) -> None:
    before = {(r["seed"], r["years_my"]): r for r in _load(before_path)}
    after = {(r["seed"], r["years_my"]): r for r in _load(after_path)}
    keys = sorted(set(before) & set(after))
    if not keys:
        print("No matching (seed, years_my) rows between the two files.")
        return

    print(f"{'seed':>12s} {'yrs My':>7s} {'sum Hc before':>15s} {'sum Hc after':>15s} {'d%':>7s}  "
          f"{'sum Hm before':>15s} {'sum Hm after':>15s} {'d%':>7s}  "
          f"{'land% before':>12s} {'land% after':>11s}  {'nodes before':>12s} {'nodes after':>11s}")
    for seed, years_my in keys:
        b, a = before[(seed, years_my)], after[(seed, years_my)]
        dhc = 100.0 * (a["sum_hc_m"] - b["sum_hc_m"]) / b["sum_hc_m"] if b["sum_hc_m"] else float("nan")
        dhm = 100.0 * (a["sum_hm_m"] - b["sum_hm_m"]) / b["sum_hm_m"] if b["sum_hm_m"] else float("nan")
        print(f"{seed:>12d} {years_my:>7.0f} {b['sum_hc_m']:>15,.0f} {a['sum_hc_m']:>15,.0f} {dhc:>6.1f}%  "
              f"{b['sum_hm_m']:>15,.0f} {a['sum_hm_m']:>15,.0f} {dhm:>6.1f}%  "
              f"{100*b['land_fraction']:>11.2f}% {100*a['land_fraction']:>10.2f}%  "
              f"{b['total_nodes']:>12d} {a['total_nodes']:>11d}")

    # Summary at the final checkpoint only -- the headline "did this help" numbers.
    final_years = max(years_my for _, years_my in keys)
    final_keys = [k for k in keys if k[1] == final_years]
    b_hc = np.array([before[k]["sum_hc_m"] for k in final_keys])
    a_hc = np.array([after[k]["sum_hc_m"] for k in final_keys])
    b_hm = np.array([before[k]["sum_hm_m"] for k in final_keys])
    a_hm = np.array([after[k]["sum_hm_m"] for k in final_keys])
    b_land = np.array([before[k]["land_fraction"] for k in final_keys])
    a_land = np.array([after[k]["land_fraction"] for k in final_keys])
    b_nodes = np.array([before[k]["total_nodes"] for k in final_keys])
    a_nodes = np.array([after[k]["total_nodes"] for k in final_keys])
    print(f"\nAt {final_years:.0f} My, across {len(final_keys)} seeds:")
    print(f"  mean sum Hc:  before {b_hc.mean():,.0f}  after {a_hc.mean():,.0f}  "
          f"({100*(a_hc.mean()-b_hc.mean())/b_hc.mean():+.2f}%)")
    print(f"  mean sum Hm:  before {b_hm.mean():,.0f}  after {a_hm.mean():,.0f}  "
          f"({100*(a_hm.mean()-b_hm.mean())/b_hm.mean():+.2f}%)")
    print(f"  mean land_fraction:  before {b_land.mean():.4f}  after {a_land.mean():.4f}  "
          f"(delta {a_land.mean()-b_land.mean():+.4f})")
    print(f"  mean total_nodes:  before {b_nodes.mean():.0f}  after {a_nodes.mean():.0f}  "
          f"({100*(a_nodes.mean()-b_nodes.mean())/b_nodes.mean():+.2f}%)")
    print(f"  max total_nodes (any seed):  before {b_nodes.max():.0f}  after {a_nodes.max():.0f}"
          "  -- watch for runaway growth in the 'after' column specifically")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--label", help="Tag for this run's rows (e.g. 'before' or 'after').")
    parser.add_argument("--out", type=Path, help="Where to write this run's JSONL rows.")
    parser.add_argument("--compare", nargs=2, metavar=("BEFORE_JSONL", "AFTER_JSONL"),
                         help="Skip running; just print the before/after comparison table.")
    args = parser.parse_args()

    if args.compare:
        compare(Path(args.compare[0]), Path(args.compare[1]))
        return

    if not args.label or not args.out:
        parser.error("--label and --out are required unless --compare is given")
    run_sweep(args.label, args.out)


if __name__ == "__main__":
    main()
