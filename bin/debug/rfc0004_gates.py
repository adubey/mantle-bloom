#!/usr/bin/env python3
"""RFC 0004 formation and craton-share gates, from #403's checked-in replay JSON (issue #409).

Needs no replay or save: reads `analysis/issue403/tuning.json` and `long.json` (seeds 7, 11,
2024 to 1 Gyr and seed 456862899 to 3 Gyr at `2562f63`), plus the cratonic area inside the
genuine set at each of their checkpoints, from `analysis/issue409/craton_genuine_area.jsonl`
(a replay of the same runs by `craton_genuine_area.py`, checked here against #403's genuine
and total craton areas). It reports, per run:

- formation timing: the single-sample pulse at CRATON_FORMATION_MYR, when quiet interior
  that generation left unseeded (between CRATON_FORMATION_MARGIN_KM and
  CRATON_SEED_MARGIN_KM in) matures together, and formation after it;
- absolute inventories at each 100 Myr checkpoint: genuine continental area, craton area
  (all of it, and the part inside the genuine set), non-cratonic genuine area (the host),
  craton volume and cumulative `formed_m3`;
- an upper bound on the area that formation converted (formed volume over the 28 km host
  floor) against the net loss of host. Host gained in the meantime only adds to the loss,
  so the share of the net loss that the bound can't explain is a floor on non-formation
  loss;
- for each 500 Myr window, the original qualifier (genuine area >= 25% of initial, tested
  at the window start and throughout) and the revised one (non-cratonic genuine area at
  the window start >= 10% of its initial value), plus post-pulse formation in the window.

With --generation it also generates each seed's world and measures the quiet host that
seeding leaves uncratonised with its clock at 0: the band the 200 Myr pulse comes from.

Usage:
    python3 bin/debug/rfc0004_gates.py [--json analysis/issue409/gates.json]
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/rfc0004_gates.py --generation \
        --json analysis/issue409/gates.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCES = (ROOT / "analysis/issue403/tuning.json", ROOT / "analysis/issue403/long.json")
GENUINE_CRATON = ROOT / "analysis/issue409/craton_genuine_area.jsonl"
# The replay must reproduce #403's trajectory: genuine and total craton areas agree to this.
MATCH_RTOL = 1e-6
# Mirrors cratons.CRATON_FORMATION_MYR and CRATON_HOST_MIN_HC_M at 2562f63; kept literal so
# the script reads the published JSON without importing a model that may have moved since.
PULSE_MYR = 200.0
HOST_MIN_HC_KM = 28.0
WINDOW_MYR = 500.0
GENUINE_FLOOR = 0.25
NON_CRATONIC_FLOOR = 0.10


def _checkpoint(row: dict, replay: dict) -> dict:
    genuine = row["genuine_continental_area_km2"] / 1e6
    craton = row["craton_area_km2"] / 1e6
    for key in ("genuine_continental_area_km2", "craton_area_km2"):
        if abs(replay[key] - row[key]) > MATCH_RTOL * max(row[key], 1.0):
            raise SystemExit(f"replay diverges from #403 at {row['myr']} Myr: {key} {replay[key]} != {row[key]}")
    craton_genuine = replay["craton_genuine_area_km2"] / 1e6
    return {
        "myr": row["myr"],
        "genuine_mkm2": genuine,
        "craton_mkm2": craton,
        "craton_genuine_mkm2": craton_genuine,
        "non_cratonic_mkm2": genuine - craton_genuine,
        "craton_share": row["craton_continental_fraction"],
        "craton_mkm3": row["craton_volume_km3"] / 1e6,
        "formed_mkm3": row["craton_ledger_km3"]["formed_m3"] / 1e6,
        "continental_inventory_mkm3": row["continental_material_ledger"]["surface_continental_derived_m3"] / 1e15,
    }


def _pulse_mkm3(run: dict) -> float:
    jumps = run["ledger_history"]["craton_largest_jumps"].get("formed_m3", [])
    return sum(j["mkm3"] for j in jumps if j["myr"] == PULSE_MYR)


def _at(points: list[dict], myr: float) -> dict | None:
    return next((p for p in points if p["myr"] == myr), None)


def _windows(points: list[dict], pulse: float) -> list[dict]:
    first = points[0]
    end_myr = points[-1]["myr"]
    rows = []
    start = 0.0
    while start < end_myr:
        stop = start + WINDOW_MYR
        lo, hi = _at(points, start), _at(points, stop)
        inside = [p for p in points if start <= p["myr"] <= stop]
        formed = hi["formed_mkm3"] - lo["formed_mkm3"]
        if start <= PULSE_MYR < stop:
            formed -= pulse
        rows.append({
            "window_myr": [start, stop],
            "genuine_start_frac": lo["genuine_mkm2"] / first["genuine_mkm2"],
            "genuine_min_frac": min(p["genuine_mkm2"] for p in inside) / first["genuine_mkm2"],
            "non_cratonic_start_mkm2": lo["non_cratonic_mkm2"],
            "non_cratonic_start_frac": lo["non_cratonic_mkm2"] / first["non_cratonic_mkm2"],
            "qualifies_genuine_at_start": lo["genuine_mkm2"] >= GENUINE_FLOOR * first["genuine_mkm2"],
            "qualifies_genuine_throughout": all(p["genuine_mkm2"] >= GENUINE_FLOOR * first["genuine_mkm2"] for p in inside),
            "qualifies_non_cratonic_at_start": lo["non_cratonic_mkm2"] >= NON_CRATONIC_FLOOR * first["non_cratonic_mkm2"],
            "post_pulse_formed_mkm3": formed,
        })
        start = stop
    return rows


def _conversion(points: list[dict], myr: float) -> dict:
    """Formation can have turned at most formed / 28 km of host into craton by `myr`."""
    first, last = points[0], _at(points, myr)
    return {
        "myr": myr,
        "non_cratonic_lost_mkm2": first["non_cratonic_mkm2"] - last["non_cratonic_mkm2"],
        "converted_max_mkm2": last["formed_mkm3"] / HOST_MIN_HC_KM,
        "craton_area_retained": last["craton_mkm2"] / first["craton_mkm2"],
        "craton_genuine_area_retained": last["craton_genuine_mkm2"] / first["craton_genuine_mkm2"],
        "non_cratonic_retained": last["non_cratonic_mkm2"] / first["non_cratonic_mkm2"],
    }


def summarise(run: dict, replay: list[dict]) -> dict:
    points = [_checkpoint(row, rep) for row, rep in zip(run["checkpoints"], replay, strict=True)]
    pulse = _pulse_mkm3(run)
    total = points[-1]["formed_mkm3"]
    return {
        "seed": run["seed"],
        "formation": {
            "pulse_at_200_myr_mkm3": pulse,
            "post_pulse_mkm3": total - pulse,
            "total_mkm3": total,
            "last_positive_sample_myr": run["formation"]["last_positive_sample_myr"],
            "positive_samples": run["formation"]["positive_samples"],
        },
        "conversion_by_500_myr": _conversion(points, 500.0),
        "windows": _windows(points, pulse),
        "checkpoints": points,
    }


def generation_band(seed: int) -> dict:
    """Quiet, unseeded host in a freshly generated world: it all starts its clock at 0."""
    import numpy as np
    from app import cratons, world
    from app.elevation_lines import line_spacing_rad

    generated = world.generate_world(seed)
    hops = cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, line_spacing_rad(generated.node_density))
    volume = area = clock = 0.0
    for plate in generated.plates:
        host = cratons._genuine_continental(plate)
        interior = cratons._interior_hops(plate, host, hops)
        if interior is None:
            continue
        hc = plate.collect("crustal_thickness_m")
        band = host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & (interior >= hops) & (plate.collect("craton_crust_m") <= 0.0)
        areas = cratons._areas(generated, plate)
        volume += float(np.dot(hc[band], areas[band]))
        area += float(areas[band].sum())
        clock = max(clock, float(plate.collect("stable_continental_myr")[band].max(initial=0.0)))
    return {"quiet_unseeded_mkm3": volume / 1e15, "quiet_unseeded_mkm2": area / 1e12, "max_clock_myr": clock}


def _print(result: dict) -> None:
    for run in result["runs"]:
        f, c = run["formation"], run["conversion_by_500_myr"]
        print(f"seed {run['seed']}: formed {f['total_mkm3']:.3f} M km^3 = pulse {f['pulse_at_200_myr_mkm3']:.3f}"
              f" + after {f['post_pulse_mkm3']:.3f}; last positive {f['last_positive_sample_myr']:.0f} Myr")
        if "generation_band" in run:
            g = run["generation_band"]
            print(f"  generation: quiet unseeded host {g['quiet_unseeded_mkm2']:.2f} M km^2, {g['quiet_unseeded_mkm3']:.1f}"
                  f" M km^3, clocks <= {g['max_clock_myr']:.0f} Myr; pulse is"
                  f" {f['pulse_at_200_myr_mkm3'] / g['quiet_unseeded_mkm3']:.1%} of it")
        print(f"  by 500 Myr: non-cratonic lost {c['non_cratonic_lost_mkm2']:.2f} M km^2, formation converted"
              f" <= {c['converted_max_mkm2']:.2f} ({c['converted_max_mkm2'] / c['non_cratonic_lost_mkm2']:.1%});"
              f" retained craton {c['craton_area_retained']:.1%} (in genuine {c['craton_genuine_area_retained']:.1%}),"
              f" non-cratonic {c['non_cratonic_retained']:.1%}")
        for w in run["windows"]:
            print(f"  {w['window_myr'][0]:5.0f}-{w['window_myr'][1]:5.0f}: genuine start {w['genuine_start_frac']:6.1%}"
                  f" min {w['genuine_min_frac']:6.1%}  non-cratonic start {w['non_cratonic_start_mkm2']:6.2f}"
                  f" ({w['non_cratonic_start_frac']:6.1%})  qualifies start/throughout/revised"
                  f" {w['qualifies_genuine_at_start']:d}/{w['qualifies_genuine_throughout']:d}"
                  f"/{w['qualifies_non_cratonic_at_start']:d}  post-pulse formed {w['post_pulse_formed_mkm3']:.3f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=Path, help="write the summary here")
    parser.add_argument("--generation", action="store_true", help="also measure each seed's unseeded quiet band (needs the backend)")
    args = parser.parse_args()
    runs = []
    baseline = None
    lines = GENUINE_CRATON.read_text().splitlines()
    replays = {run["seed"]: run["checkpoints"] for run in map(json.loads, filter(None, lines))}
    for path in SOURCES:
        data = json.loads(path.read_text())
        baseline = data["baseline_commit"]
        runs.extend(summarise(run, replays[run["seed"]]) for run in data["runs"])
    if args.generation:
        for run in runs:
            run["generation_band"] = generation_band(run["seed"])
    result = {
        "issue": 409,
        "source": [str(p.relative_to(ROOT)) for p in (*SOURCES, GENUINE_CRATON)],
        "baseline_commit": baseline,
        "window_myr": WINDOW_MYR,
        "genuine_floor": GENUINE_FLOOR,
        "non_cratonic_floor": NON_CRATONIC_FLOOR,
        "runs": runs,
    }
    _print(result)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
