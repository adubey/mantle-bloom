#!/usr/bin/env python3
"""Measure the cratons-first layout's craton, belt and continent-size figures that the
constants in `continent_layout.py` cite in their comments (#470, #793).

For each `CONTINENT_SIZE_CONCENTRATION` value, builds the layout for seeds 1..N and the
three calibration seeds (7, 11, 2024) at land fraction 0.29, and reports over the N seeds:
craton count; share of cratons within `CRATON_EDGE_NEAR_KM` of the continental edge; belt
and core share of continental area; share of nuclei that joined a craton; craton areas
(5th percentile, share under 1 M km^2, largest, cratons of 10 M km^2 or more per world);
and each continent's share of continental area (smallest, largest). A layout takes about
1 s at node density 1.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/continent_layout_measures.py
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/continent_layout_measures.py \
        --seeds 40 --concentration 2 32 --out analysis/issue793/layout-measures.json
"""

import argparse
import json

import numpy as np

from app import continent_layout as cl
from app.lithosphere_plate import cratons_first_layout

CALIBRATION_SEEDS = (7, 11, 2024)
M_KM2 = 1e6


def _measure(seed: int) -> dict:
    layout = cratons_first_layout(seed, land_fraction=0.29).layout
    diag = layout.diagnostics()
    area = layout.lattice.area_km2
    sizes = np.array([area[layout.core == c].sum() for c in np.unique(layout.core[layout.is_core])]) / M_KM2
    continents = np.array([c["area_km2"] for c in diag["continents"]])
    shares = continents[continents > 0] / continents.sum()
    return {
        "seed": seed,
        "cratons": diag["cratons"],
        "near_edge": diag["cratons_near_edge_share"],
        "belt": diag["belt_share_of_continent"],
        "core": diag["core_share_of_continent"],
        "joined": diag["joined_nucleus_share"],
        "sizes": sizes,
        "smallest_continent": float(shares.min()),
        "largest_continent": float(shares.max()),
    }


def _summary(rows: list[dict]) -> dict:
    sizes = np.concatenate([r["sizes"] for r in rows])
    column = lambda key: np.array([r[key] for r in rows], dtype=float)  # noqa: E731
    return {
        "seeds": len(rows),
        "cratons_min_max": [int(column("cratons").min()), int(column("cratons").max())],
        "near_edge_median_min_max": [round(float(f(column("near_edge"))), 3) for f in (np.median, np.min, np.max)],
        "belt_share_min_max": [round(float(column("belt").min()), 3), round(float(column("belt").max()), 3)],
        "core_share_mean_min_max": [round(float(f(column("core"))), 3) for f in (np.mean, np.min, np.max)],
        "joined_nucleus_share_mean": round(float(column("joined").mean()), 3),
        "cratons_of_10m_km2_median_per_world": float(np.median([int((r["sizes"] >= 10).sum()) for r in rows])),
        "craton_area_p5_m_km2": round(float(np.percentile(sizes, 5)), 2),
        "craton_share_under_1m_km2": round(float((sizes < 1).mean()), 3),
        "craton_area_max_m_km2": round(float(sizes.max()), 1),
        "smallest_continent_share_min": round(float(column("smallest_continent").min()), 3),
        "largest_continent_share_max": round(float(column("largest_continent").max()), 3),
        "seeds_outside_5pct_40pct": int(((column("smallest_continent") < 0.05) | (column("largest_continent") > 0.40)).sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=int, default=40, help="seeds 1..N (default 40)")
    parser.add_argument("--concentration", type=float, nargs="+", default=[cl.CONTINENT_SIZE_CONCENTRATION])
    parser.add_argument("--out", help="write the summaries as JSON")
    args = parser.parse_args()
    result = {}
    for concentration in args.concentration:
        cl.CONTINENT_SIZE_CONCENTRATION = concentration
        rows = [_measure(seed) for seed in range(1, args.seeds + 1)]
        result[f"concentration_{concentration:g}"] = {
            **_summary(rows),
            "core_share_on_calibration_seeds": {str(s): round(_measure(s)["core"], 3) for s in CALIBRATION_SEEDS},
        }
        print(f"concentration {concentration:g}: {json.dumps(result[f'concentration_{concentration:g}'])}", flush=True)
    if args.out:
        with open(args.out, "w") as handle:
            json.dump(result, handle, indent=1)
            handle.write("\n")


if __name__ == "__main__":
    main()
