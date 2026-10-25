#!/usr/bin/env python3
"""Plume-breakup comparator (issue #582).

Measures, per seed, the fraction of continental plume-head arrivals followed by
a split of the host plate within 20 Myr, beside two plume-off references:

- matched: fraction of the same scheduled head times whose following 20 Myr
  window contains any split in a plume-off twin run (multiplier 0, identical
  schedule thinned away, same seed/controls/steps);
- background: splits per 500 Myr in the twin run (expect ~26/28/36 on seeds
  7/11/2024 from #407; deviation flags baseline drift, not plume physics).

Method follows bin/debug/attribute_reset_sources.py (#407): default
`world.generate_world(seed)`, cratons on, 5 Myr steps to 500 Myr. Every tap
wraps a function, records, delegates, and is restored afterwards; taps never
mutate world state. Split counting matches #407 exactly: `maybe_split_plate`
successes only (defragmentation "fragmented" landmasses excluded on both sides).

A head whose host plate is absorbed by a merge or consumed by defragmentation
inside its 20 Myr window is censored (reported, excluded from the denominator);
a host that splits keeps its id on daughter A, so hits are never censored.

Usage:
    backend/.venv/bin/python bin/debug/compare_plume_breakup.py \
        --out analysis/issue582/breakup.json
    backend/.venv/bin/python bin/debug/compare_plume_breakup.py --probe-only \
        --out /tmp/breakup-probe.json

One world runs at a time (replays are RAM-heavy). No saves are written; only the
compact JSON and analysis/issue582/README.md are reviewed.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
from app import merge_split, plumes  # noqa: E402
from app import plates as plates_mod  # noqa: E402
from app import world as world_mod  # noqa: E402

SCRIPT_NAME = Path(__file__).name
WINDOW_YEARS = 20e6
BACKGROUND_SPLITS = {7: 26, 11: 28, 2024: 36}  # per 500 Myr, analysis/issue407/README.md


class BreakupTaps:
    """Record head arrivals, splits, and host-id disappearances per run."""

    def __init__(self) -> None:
        self.heads: list[dict] = []  # {born, host, continental}
        self.splits: list[dict] = []  # {time, parent}
        self.vanished: list[dict] = []  # {time, plate} merged/defragged away
        self._saved: list[tuple[object, str, object]] = []

    def _wrap(self, obj, name, make_wrapper) -> None:
        original = getattr(obj, name)
        setattr(obj, name, make_wrapper(original))
        self._saved.append((obj, name, original))

    def install(self) -> None:
        taps = self

        def apply_wrapper(original):
            def apply(world, years, *args, **kwargs):
                t1 = world.elapsed_years
                out = original(world, years, *args, **kwargs)
                for plume in world.mantle_plumes:
                    if plume.is_tail_only or not (t1 - years < plume.born_years <= t1):
                        continue
                    plates = [p for p in world.plates if p.node_count() > 0]
                    continental, host = plumes._center_host(plates, plume.position)
                    taps.heads.append(
                        {"born": float(plume.born_years), "host": host, "continental": bool(continental)}
                    )
                return out

            return apply

        def split_wrapper(original):
            def maybe_split_plate(world, plate, *args, **kwargs):
                out = original(world, plate, *args, **kwargs)
                if out is not None:
                    taps.splits.append({"time": float(world.elapsed_years), "parent": plate.plate_id})
                return out

            return maybe_split_plate

        def merge_wrapper(original):
            def _merge(world, keep, absorb, *args, **kwargs):
                out = original(world, keep, absorb, *args, **kwargs)
                taps.vanished.append({"time": float(world.elapsed_years), "plate": absorb})
                return out

            return _merge

        def defrag_wrapper(original):
            def defragment(self, next_id, connect_radius_rad, min_fragment_nodes, world):
                old = self.plate_id
                out = original(self, next_id, connect_radius_rad, min_fragment_nodes, world)
                if out is not None:
                    replacements, _, _ = out
                    if all(r.plate_id != old for r in replacements):
                        taps.vanished.append({"time": float(world.elapsed_years), "plate": old})
                return out

            return defragment

        self._wrap(plumes, "apply", apply_wrapper)
        self._wrap(merge_split, "maybe_split_plate", split_wrapper)
        self._wrap(merge_split, "_merge_plates_with_craton_audit", merge_wrapper)
        self._wrap(plates_mod.Plate, "defragment", defrag_wrapper)

    def uninstall(self) -> None:
        for obj, name, original in reversed(self._saved):
            setattr(obj, name, original)
        self._saved.clear()


def _run(seed: int, years: float, step: float, plume_on: bool) -> dict:
    """One 500 Myr world; returns tap records plus wall time."""
    taps = BreakupTaps()
    taps.install()
    try:
        world = world_mod.generate_world(seed)
        world.plume_multiplier = 1.0 if plume_on else 0.0
        t0 = time.time()
        while world.elapsed_years < years - 1.0:
            world_mod.step_world(world, min(step, years - world.elapsed_years))
        wall = time.time() - t0
        return {
            "seed": seed,
            "plume_on": plume_on,
            "wall_s": round(wall, 1),
            "heads": taps.heads,
            "splits": taps.splits,
            "vanished": taps.vanished,
        }
    finally:
        taps.uninstall()


def _analyze(on: dict, off: dict) -> dict:
    """Head-following split fraction vs matched plume-off windows."""
    heads = [h for h in on["heads"] if h["continental"]]
    vanish_at = {}
    for v in on["vanished"]:
        vanish_at.setdefault(v["plate"], v["time"])
    hits = censored = 0
    for h in heads:
        lo, hi = h["born"], h["born"] + WINDOW_YEARS
        if any(s["parent"] == h["host"] and lo < s["time"] <= hi for s in on["splits"]):
            hits += 1
        elif h["host"] in vanish_at and lo < vanish_at[h["host"]] <= hi:
            censored += 1
    matched = 0
    for h in heads:
        lo, hi = h["born"], h["born"] + WINDOW_YEARS
        if any(lo < s["time"] <= hi for s in off["splits"]):
            matched += 1
    return {
        "seed": on["seed"],
        "continental_heads": len(heads),
        "oceanic_heads": sum(1 for h in on["heads"] if not h["continental"]),
        "head_followed_by_host_split": hits,
        "censored_windows": censored,
        "head_split_fraction": round(hits / max(1, len(heads) - censored), 4),
        "matched_off_windows_with_any_split": matched,
        "matched_off_fraction": round(matched / max(1, len(heads)), 4),
        "on_splits_per_500myr": len(on["splits"]),
        "off_splits_per_500myr": len(off["splits"]),
        "off_background_expected": BACKGROUND_SPLITS.get(on["seed"]),
        "on_wall_s": on["wall_s"],
        "off_wall_s": off["wall_s"],
    }


def _provenance() -> dict:
    root = Path(__file__).resolve().parents[2]
    def git(*args):
        try:
            return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, timeout=30).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return "unavailable"

    return {
        "script": SCRIPT_NAME,
        "git_head": git("rev-parse", "HEAD"),
        "git_status_short": git("status", "--short"),
        "python": platform.python_version(),
        "controls": {"steps": "default generate_world, cratons on", "plume_multiplier_on": 1.0, "plume_multiplier_off": 0.0},
        "knobs": {
            "PLUME_CRATON_WEAKENING": plumes.PLUME_CRATON_WEAKENING,
            "PLUME_HEAD_STRESS": plumes.PLUME_HEAD_STRESS,
            "PLUME_HEAT_MYR": plumes.PLUME_HEAT_MYR,
            "STRESS_CERTAIN_RIFT": merge_split.STRESS_CERTAIN_RIFT,
            "RIFT_SUCCESS_PROBABILITY": merge_split.RIFT_SUCCESS_PROBABILITY,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True, help="JSON output path")
    parser.add_argument("--seeds", default="7,11,2024", help="comma-separated seeds")
    parser.add_argument("--myr", type=float, default=500.0, help="horizon per run in Myr")
    parser.add_argument("--step-myr", type=float, default=5.0, help="step length in Myr")
    parser.add_argument("--probe-only", action="store_true", help="one short seed-7 run to verify taps and JSON shape")
    args = parser.parse_args()

    seeds = [7] if args.probe_only else [int(s) for s in args.seeds.split(",")]
    horizon = 20e6 if args.probe_only else args.myr * 1e6
    step = args.step_myr * 1e6
    results = []
    for seed in seeds:
        on = _run(seed, horizon, step, True)
        print(json.dumps({"seed": seed, "plume_on": True, "heads": len(on["heads"]), "splits": len(on["splits"]), "wall_s": on["wall_s"]}), flush=True)
        if args.probe_only:
            results.append({"seed": seed, "probe_on": on})
            continue
        off = _run(seed, horizon, step, False)
        print(json.dumps({"seed": seed, "plume_on": False, "heads": len(off["heads"]), "splits": len(off["splits"]), "wall_s": off["wall_s"]}), flush=True)
        results.append(_analyze(on, off))
    out = {"provenance": _provenance(), "window_years": WINDOW_YEARS, "results": results}
    Path(args.out).write_text(json.dumps(out, indent=2))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
