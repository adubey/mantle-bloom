#!/usr/bin/env python3
"""Compare RFC 0014 lake storage at 10 kyr and 100 kyr under one frozen forcing.

Loads a saved world, captures the production lake forcing during one 10 kyr step, then
solves the same terrain, prior liquid/ice state, climate rates, and temperature inputs
at both step lengths. This isolates the lake timestep response from tectonic and climate
feedback. RFC 0014 currently keeps the endpoint optimization disabled, so every measured
body is classified as prognostic-C fallback.

Run one save at a time, never alongside another replay, probe, or test:

    backend/.venv/bin/python bin/debug/measure_frozen_lake_step_size.py \
        <held-out-save.mbworld> --out <scratch>/step-size.json

The compact result records the save checksum, commit, aggregate lake volume and matched
body-level liquid-plus-ice differences. Saves and raw replay outputs stay out of the repo.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))

from app import persistence, world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402


def _args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def _components(depth, ice, area, neighbors):
    wet = depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    groups = hydrology.lake_components(wet, neighbors)
    return [
        {
            "nodes": nodes,
            "area_km2": float(area[nodes].sum() / 1e6),
            "volume_m3": float(np.dot(depth[nodes] + ice[nodes], area[nodes])),
        }
        for nodes in groups
    ]


def _matched_differences(left, right):
    candidates = []
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            overlap = np.intersect1d(a["nodes"], b["nodes"], assume_unique=True).size
            if overlap:
                candidates.append((overlap, i, j))
    candidates.sort(reverse=True)
    used_left, used_right, matches = set(), set(), []
    for overlap, i, j in candidates:
        if i in used_left or j in used_right:
            continue
        smaller = min(len(left[i]["nodes"]), len(right[j]["nodes"]))
        if overlap < max(1, smaller // 2):
            continue
        used_left.add(i)
        used_right.add(j)
        a, b = left[i], right[j]
        if min(a["area_km2"], b["area_km2"]) < 10_000.0:
            continue
        relative = abs(a["volume_m3"] - b["volume_m3"]) / max(a["volume_m3"], b["volume_m3"], 1.0)
        matches.append({
            "area_10k_km2": a["area_km2"],
            "area_100k_km2": b["area_km2"],
            "volume_10k_m3": a["volume_m3"],
            "volume_100k_m3": b["volume_m3"],
            "relative_difference": relative,
            "pass_5pct_fallback_bound": bool(relative <= 0.05),
        })
    return matches


def main() -> None:
    args = _args()
    raw = args.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    world.ice_age_period_years = 0.0
    captured: dict = {}
    original_solve = lakes.solve_lakes_prognostic
    original_ice = hydrology.advance_lake_ice

    def capture_solve(forest, elevation, previous_depth, years, climate, hypsometry, **kwargs):
        captured["solve"] = copy.deepcopy(
            (forest, elevation, previous_depth, years, climate, hypsometry, kwargs.get("is_ocean"))
        )
        return original_solve(forest, elevation, previous_depth, years, climate, hypsometry, **kwargs)

    def capture_ice(prior_ice, liquid_depth, area, mean_c, amplitude_c, ice_cover_fraction, years):
        captured["phase"] = copy.deepcopy(
            (prior_ice, area, mean_c, amplitude_c, ice_cover_fraction)
        )
        return original_ice(prior_ice, liquid_depth, area, mean_c, amplitude_c, ice_cover_fraction, years)

    lakes.solve_lakes_prognostic = capture_solve
    hydrology.advance_lake_ice = capture_ice
    started = time.perf_counter()
    world_mod.step_world(world, 10_000.0)
    capture_runtime_s = time.perf_counter() - started
    lakes.solve_lakes_prognostic = original_solve
    hydrology.advance_lake_ice = original_ice
    if "solve" not in captured or "phase" not in captured:
        raise RuntimeError("failed to capture production lake and ice forcing")

    fields = world.hydrology_cache
    neighbors = np.asarray(fields.neighbor_idx).copy()
    source_forest, elevation, previous_depth, years, climate, hypsometry, is_ocean = captured["solve"]
    prior_ice, area, mean_c, amplitude_c, ice_cover = captured["phase"]
    outputs = {}
    components = {}
    for label, years in (("10k", 10_000.0), ("100k", 100_000.0)):
        result = original_solve(
            copy.deepcopy(source_forest), elevation, previous_depth, years, climate,
            hypsometry, is_ocean=is_ocean,
        )
        ice, liquid, _ice_budget = original_ice(
            prior_ice, result.lake_depth.copy(), area, mean_c, amplitude_c, ice_cover, years,
        )
        outputs[label] = {
            "lake_volume_m3": float(np.dot(liquid, hypsometry.area_m2)),
            "ice_volume_m3": float(np.dot(ice, hypsometry.area_m2)),
            "combined_volume_m3": float(np.dot(liquid + ice, hypsometry.area_m2)),
            "body_count": int(len(hydrology.lake_components(liquid > hydrology.LAKE_MIN_VISIBLE_DEPTH_M, neighbors))),
        }
        components[label] = _components(liquid, ice, hypsometry.area_m2, neighbors)

    matched = _matched_differences(components["10k"], components["100k"])
    result = {
        "issue": 673,
        "rfc": "0014",
        "backend_commit": __import__("subprocess").run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip(),
        "save_sha256": hashlib.sha256(raw).hexdigest(),
        "ice_age_control": "never",
        "forcing_source": "captured once from a 10 kyr production step, then frozen across both solves",
        "capture_step_years": 10_000.0,
        "capture_step_runtime_s": capture_runtime_s,
        "solver": "prognostic C for all bodies; endpoint optimization disabled",
        "results": outputs,
        "matched_bodies_over_10000_km2": matched,
        "matched_large_body_count": len(matched),
        "fallback_5pct_pass_count": sum(item["pass_5pct_fallback_bound"] for item in matched),
        "fallback_5pct_fail_count": sum(not item["pass_5pct_fallback_bound"] for item in matched),
        "largest_relative_difference": max((item["relative_difference"] for item in matched), default=None),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
