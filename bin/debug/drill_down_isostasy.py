#!/usr/bin/env python3
"""Drill-down for issue #189 ("isostasy formula reads out MAX_ELEVATION_M for land whose Hc is
well under its own cap"). #189 pooled #180's rerun (10 sweep_lib.SWEEP_SEEDS seeds, 120 My,
collision_uplift_multiplier in {1.0, 0.25}) and found that of the land nodes pinned at
elevation_lines.MAX_ELEVATION_M (9000m), ~67-69% have crustal_thickness_m under half of its own
cap (lithosphere.MAX_CRUSTAL_THICKNESS_M) -- too far from that cap for the pin to be explained by
`apply_convergent_deformation`'s own Hc/Hm overflow path.

This reruns the same (seed, collision_uplift_multiplier) grid and, for every pinned land node,
independently recomputes `lithosphere.isostatic_elevation(hc, hm, rho_c)` from that node's own
recorded Hc/Hm/rho_c -- i.e. what the node's *own* Airy-isostasy column actually supports, with
the same clip applied so it's an apples-to-apples comparison against the stored `elevation`. If
that recomputed value agrees with the stored 9000m, the formula itself (given the Hc/Hm on
record) really does read out the cap -- pointing at the formula/calibration, or at anomalous
Hm/rho_c. If it *disagrees* (comes out well under 9000), the 9000m stamped on the line isn't
something this node's own crustal column backs at all -- i.e. "isostatic debt" (the same concept
bin/debug/drill_down_volcanism.py's `_total_isostatic_debt_km3` already measures in aggregate for
issue #173), meaning some other code path wrote `elevation` directly without a matching Hc/Hm
change, pointing at a bypass bug elsewhere rather than at isostatic_elevation itself.

Usage: backend/.venv/bin/python bin/debug/drill_down_isostasy.py
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import SWEEP_SEEDS  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import lithosphere  # noqa: E402
from app.elevation_lines import ELEV_CHANGE_LABELS, MAX_ELEVATION_M  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

NODE_DENSITY = 1.0
STEP_YEARS = 10_000_000
TOTAL_YEARS = 120_000_000
COLLISION_UPLIFT_MULTIPLIERS = (1.0, 0.25)

HC_CAP = lithosphere.MAX_CRUSTAL_THICKNESS_M
PIN_EPS_M = 1e-6  # np.clip's own output is exactly MAX_ELEVATION_M, no float slop expected


def run_one(seed: int, collision_uplift_multiplier: float) -> list[dict]:
    world = generate_world(seed=seed, node_density=NODE_DENSITY)
    world.collision_uplift_multiplier = collision_uplift_multiplier

    years_done = 0.0
    while years_done < TOTAL_YEARS:
        step = min(STEP_YEARS, TOTAL_YEARS - years_done)
        step_world(world, years=step)
        years_done += step

    rows = []
    for plate in world.plates:
        if not plate.lines:
            continue
        elevation = plate.collect("elevation")
        land = elevation > world.sea_level_m
        pinned = land & (elevation >= MAX_ELEVATION_M - PIN_EPS_M)
        if not np.any(pinned):
            continue

        hc = plate.collect("crustal_thickness_m")[pinned]
        hm = plate.collect("mantle_lithosphere_thickness_m")[pinned]
        crust_type_code = plate.collect("crust_type_code")[pinned]
        reason = plate.collect("elev_change_reason")[pinned]
        rho_c = lithosphere.node_crust_density(crust_type_code, plate.crust_type)
        recomputed = lithosphere.isostatic_elevation(hc, hm, rho_c)

        for i in range(hc.size):
            rows.append({
                "seed": seed,
                "collision_uplift_multiplier": collision_uplift_multiplier,
                "plate_crust_type": plate.crust_type,
                "hc_m": float(hc[i]),
                "hm_m": float(hm[i]),
                "hc_cap_ratio": float(hc[i] / HC_CAP),
                "rho_c": float(rho_c[i]) if np.ndim(rho_c) else float(rho_c),
                "elev_change_reason": int(reason[i]),
                "recomputed_elevation_m": float(recomputed[i]),
                "matches_stored": bool(recomputed[i] >= MAX_ELEVATION_M - PIN_EPS_M),
            })
    return rows


def main() -> None:
    jobs = [(seed, m) for seed in SWEEP_SEEDS for m in COLLISION_UPLIFT_MULTIPLIERS]
    print(f"{len(jobs)} job(s) to run.", file=sys.stderr)

    all_rows: list[dict] = []
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(run_one, seed, m): (seed, m) for seed, m in jobs}
        for n, future in enumerate(as_completed(futures), start=1):
            seed, m = futures[future]
            try:
                rows = future.result()
            except Exception as exc:  # noqa: BLE001 -- keep the sweep going past one bad job
                print(f"FAILED seed={seed} x{m}: {exc!r}", file=sys.stderr)
                continue
            all_rows.extend(rows)
            elapsed = time.perf_counter() - t0
            print(f"[{n}/{len(jobs)}] seed={seed} x{m} done, {len(rows)} pinned land node(s) "
                  f"({elapsed:.0f}s elapsed)", file=sys.stderr)

    out_dir = Path(__file__).resolve().parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "drill_down_isostasy.jsonl"
    with out_path.open("w") as f:
        for row in all_rows:
            f.write(json.dumps(row) + "\n")
    print(f"\nWrote {len(all_rows)} row(s) to {out_path}", file=sys.stderr)

    if not all_rows:
        print("No pinned land nodes found at all -- nothing to analyze.")
        return

    low_hc = [r for r in all_rows if r["hc_cap_ratio"] < 0.50]
    overflow = [r for r in all_rows if r["hc_cap_ratio"] >= 0.99]
    mid = [r for r in all_rows if 0.50 <= r["hc_cap_ratio"] < 0.99]

    print(f"\nTotal pinned land nodes: {len(all_rows)}")
    print(f"  Hc >= 99% of cap (overflow):        {len(overflow)} ({100 * len(overflow) / len(all_rows):.1f}%)")
    print(f"  Hc in [50%, 99%) of cap:             {len(mid)} ({100 * len(mid) / len(all_rows):.1f}%)")
    print(f"  Hc < 50% of cap:                     {len(low_hc)} ({100 * len(low_hc) / len(all_rows):.1f}%)")

    if low_hc:
        matches = [r for r in low_hc if r["matches_stored"]]
        debt = [r for r in low_hc if not r["matches_stored"]]
        print(f"\nOf the Hc<50%-of-cap pinned group ({len(low_hc)} nodes):")
        print(f"  Own Hc/Hm/rho_c formula recomputes to >=9000m (formula-driven): "
              f"{len(matches)} ({100 * len(matches) / len(low_hc):.1f}%)")
        print(f"  Own Hc/Hm/rho_c formula recomputes BELOW 9000m (isostatic debt -- some other "
              f"path wrote elevation without backing Hc/Hm): {len(debt)} ({100 * len(debt) / len(low_hc):.1f}%)")

        if debt:
            gaps = np.array([MAX_ELEVATION_M - r["recomputed_elevation_m"] for r in debt])
            print(f"\n  Debt group: recomputed elevation vs stored 9000m -- "
                  f"mean gap {gaps.mean():.0f}m, median {np.median(gaps):.0f}m, max {gaps.max():.0f}m")
            reasons = Counter(r["elev_change_reason"] for r in debt)
            print("  Debt group elev_change_reason breakdown:")
            for code, count in reasons.most_common():
                label = ELEV_CHANGE_LABELS[code] if 0 <= code < len(ELEV_CHANGE_LABELS) else f"unknown({code})"
                print(f"    {label:35s} {count:6d} ({100 * count / len(debt):.1f}%)")

        if matches:
            hm_vals = np.array([r["hm_m"] for r in matches])
            hc_vals = np.array([r["hc_m"] for r in matches])
            rho_vals = np.array([r["rho_c"] for r in matches])
            print(f"\n  Formula-driven group ({len(matches)} nodes): "
                  f"Hc mean/min/max = {hc_vals.mean():.0f}/{hc_vals.min():.0f}/{hc_vals.max():.0f}m, "
                  f"Hm mean/min/max = {hm_vals.mean():.0f}/{hm_vals.min():.0f}/{hm_vals.max():.0f}m, "
                  f"rho_c values present: {sorted(set(rho_vals.tolist()))}")
            reasons = Counter(r["elev_change_reason"] for r in matches)
            print("  Formula-driven group elev_change_reason breakdown:")
            for code, count in reasons.most_common():
                label = ELEV_CHANGE_LABELS[code] if 0 <= code < len(ELEV_CHANGE_LABELS) else f"unknown({code})"
                print(f"    {label:35s} {count:6d} ({100 * count / len(matches):.1f}%)")


if __name__ == "__main__":
    main()
