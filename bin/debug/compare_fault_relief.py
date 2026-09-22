#!/usr/bin/env python3
"""Paired long-run Hc/Hm sweep for issue #213.

Uses sweep_lib's fresh-world generation, 10 My steps, and first ten reported-problem
seeds, plus twenty deterministic seeds. Each seed runs the pre-#213 additive relief
and the current conserved relief from the same initial state. Output is resumable
JSONL, with a compact JSON summary and dependency-free SVG time-series plot.

Usage:
    backend/.venv/bin/python bin/debug/compare_fault_relief.py --workers 4
    backend/.venv/bin/python bin/debug/compare_fault_relief.py --seeds 829071382 --checkpoints 0,10,20
"""

from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import (  # noqa: E402
    BASELINE_AVG_RATE_CM_YR, BASELINE_MAX_RATE_CM_YR, NODE_DENSITY,
    STEP_YEARS, SWEEP_SEEDS, _apply_rotation_rate_overrides,
)
from legacy_fault_relief import apply_plate_fault_relief as legacy_relief  # noqa: E402
from plot_fault_relief import render as render_png  # noqa: E402
from app import faults, lithosphere  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402
import numpy as np  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = THIS_DIR / "results" / "issue213_hc_hm.jsonl"


def default_seeds() -> tuple[int, ...]:
    rng = random.Random(213)
    extra = []
    while len(extra) < 20:
        seed = rng.randrange(1, 1_000_000_000)
        if seed not in SWEEP_SEEDS and seed not in extra:
            extra.append(seed)
    return (*SWEEP_SEEDS, *extra)


def measure(world) -> dict:
    hc = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates])
    hm = np.concatenate([p.collect("mantle_lithosphere_thickness_m") for p in world.plates])
    elevation = np.concatenate([p.collect("elevation") for p in world.plates])
    land = elevation > world.sea_level_m
    result = {"node_count": int(hc.size), "land_node_count": int(land.sum())}
    for name, values, cap in (
        ("hc", hc, lithosphere.MAX_CRUSTAL_THICKNESS_M),
        ("hm", hm, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M),
    ):
        result.update({
            f"{name}_mean_m": float(values.mean()),
            f"{name}_std_m": float(values.std()),
            f"{name}_p95_m": float(np.percentile(values, 95)),
            f"{name}_max_m": float(values.max()),
            f"{name}_at_max_fraction": float(np.mean(values >= cap - 1e-6)),
            f"land_{name}_mean_m": float(values[land].mean()) if np.any(land) else None,
        })
    result["land_above_3km_fraction"] = float(np.mean(elevation[land] > 3000)) if np.any(land) else None
    result["land_hc_70km_fraction"] = float(np.mean(hc[land] >= 70000)) if np.any(land) else None
    return result


def run_seed(seed: int, node_density: float, checkpoints: tuple[int, ...]) -> list[dict]:
    rows = []
    original = faults._apply_plate_fault_relief
    try:
        for approach, relief in (("old_additive", legacy_relief), ("conserved", original)):
            faults._apply_plate_fault_relief = relief
            _apply_rotation_rate_overrides(BASELINE_AVG_RATE_CM_YR, BASELINE_MAX_RATE_CM_YR)
            world = generate_world(seed=seed, node_density=node_density)
            years_done = 0
            for checkpoint in checkpoints:
                while years_done < checkpoint:
                    step = min(STEP_YEARS, checkpoint - years_done)
                    step_world(world, years=step)
                    years_done += step
                rows.append({
                    "seed": seed, "approach": approach, "checkpoint_years": checkpoint,
                    "node_density": node_density, "step_years": STEP_YEARS,
                    **measure(world),
                })
    finally:
        faults._apply_plate_fault_relief = original
    return rows


METRICS = ("hc_mean_m", "hc_p95_m", "hc_at_max_fraction", "hm_mean_m", "hm_p95_m", "hm_at_max_fraction")


def summarize(rows: list[dict], seeds: list[int], checkpoints: tuple[int, ...]) -> dict:
    by_key = {(r["seed"], r["approach"], r["checkpoint_years"]): r for r in rows}
    points = []
    for checkpoint in checkpoints:
        paired = [seed for seed in seeds if all((seed, approach, checkpoint) in by_key for approach in ("old_additive", "conserved"))]
        point = {"checkpoint_years": checkpoint, "paired_seeds": len(paired), "metrics": {}}
        for metric in METRICS:
            old = [by_key[(s, "old_additive", checkpoint)][metric] for s in paired]
            new = [by_key[(s, "conserved", checkpoint)][metric] for s in paired]
            delta = [b - a for a, b in zip(old, new)]
            point["metrics"][metric] = {
                "old_mean": statistics.fmean(old) if old else None,
                "new_mean": statistics.fmean(new) if new else None,
                "paired_delta_mean": statistics.fmean(delta) if delta else None,
                "paired_delta_se": statistics.stdev(delta) / math.sqrt(len(delta)) if len(delta) > 1 else None,
            }
        points.append(point)
    return {"seed_count_requested": len(seeds), "seeds": seeds, "approaches": ["old_additive", "conserved"],
            "step_years": STEP_YEARS, "checkpoints": points}


def write_svg(summary: dict, path: Path) -> None:
    panels = (("hc_mean_m", "Mean Hc (km)", 0.001), ("hc_p95_m", "95th percentile Hc (km)", 0.001),
              ("hm_mean_m", "Mean Hm (km)", 0.001), ("hm_p95_m", "95th percentile Hm (km)", 0.001),
              ("hc_at_max_fraction", "Hc at cap (%)", 100), ("hm_at_max_fraction", "Hm at cap (%)", 100))
    points = [p for p in summary["checkpoints"] if p["paired_seeds"]]
    if not points:
        return
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="960" viewBox="0 0 1050 960">',
             '<rect width="1050" height="960" fill="#101827"/>',
             '<text x="28" y="33" fill="#e8edf5" font-family="sans-serif" font-size="20">Issue #213: paired fault relief sweep</text>',
             '<text x="28" y="55" fill="#aebbd0" font-family="sans-serif" font-size="12">Blue: old additive   Orange: conserved   Each point is the paired-seed mean</text>']
    ages = [p["checkpoint_years"] / 1e6 for p in points]
    age_max = max(ages) or 1
    for i, (metric, label, factor) in enumerate(panels):
        col, row = i % 2, i // 2
        x0, y0 = 58 + col * 520, 100 + row * 280
        w, h = 430, 195
        values = [p["metrics"][metric][key] * factor for p in points for key in ("old_mean", "new_mean")]
        lo, hi = min(values), max(values)
        pad = max((hi - lo) * 0.12, 0.01)
        lo -= pad; hi += pad
        parts.append(f'<text x="{x0}" y="{y0-18}" fill="#e8edf5" font-family="sans-serif" font-size="15">{label}</text>')
        parts.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#172338" stroke="#52627a"/>')
        parts.append(f'<text x="{x0}" y="{y0+h+17}" fill="#aebbd0" font-family="sans-serif" font-size="11">0</text>')
        parts.append(f'<text x="{x0+w-34}" y="{y0+h+17}" fill="#aebbd0" font-family="sans-serif" font-size="11">{age_max:g} My</text>')
        parts.append(f'<text x="{x0+3}" y="{y0+13}" fill="#aebbd0" font-family="sans-serif" font-size="11">{hi:.1f}</text>')
        parts.append(f'<text x="{x0+3}" y="{y0+h-5}" fill="#aebbd0" font-family="sans-serif" font-size="11">{lo:.1f}</text>')
        for key, color in (("old_mean", "#66b5ff"), ("new_mean", "#ffae68")):
            coords = []
            for age, point in zip(ages, points):
                value = point["metrics"][metric][key] * factor
                coords.append(f'{x0+w*age/age_max:.1f},{y0+h*(hi-value)/(hi-lo):.1f}')
            parts.append(f'<polyline points="{" ".join(coords)}" fill="none" stroke="{color}" stroke-width="2.5"/>')
    parts.append('</svg>')
    path.write_text(''.join(parts))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", help="comma-separated seed IDs (default: 30 fixed seeds)")
    parser.add_argument("--checkpoints", default="0,30,60,90,120,150,180,210,240,270,300,320", help="Myr, comma-separated")
    parser.add_argument("--node-density", type=float, default=NODE_DENSITY)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    seeds = list(dict.fromkeys(int(s) for s in args.seeds.split(","))) if args.seeds else list(default_seeds())
    checkpoints = tuple(sorted({int(float(m) * 1e6) for m in args.checkpoints.split(",")}))
    if not seeds or not checkpoints or checkpoints[0] < 0 or args.workers < 1:
        parser.error("provide seeds, nonnegative checkpoints, and at least one worker")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    rows = [json.loads(line) for line in args.out.read_text().splitlines() if line.strip()] if args.out.exists() else []
    for row in rows:
        if row["node_density"] != args.node_density or row["step_years"] != STEP_YEARS:
            parser.error("output already contains a different node density or step size; choose another --out")
    complete = {(r["seed"], r["approach"], r["checkpoint_years"]) for r in rows}
    needed = {(s, approach, cp) for s in seeds for approach in ("old_additive", "conserved") for cp in checkpoints}
    pending = [s for s in seeds if any((s, approach, cp) not in complete for approach in ("old_additive", "conserved") for cp in checkpoints)]
    print(f"{len(pending)} seeds pending of {len(seeds)}; {len(checkpoints)} checkpoints; output {args.out}", flush=True)
    start = time.perf_counter()
    with args.out.open("a") as out, ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_seed, s, args.node_density, checkpoints): s for s in pending}
        for i, future in enumerate(as_completed(futures), 1):
            seed = futures[future]
            try:
                records = future.result()
            except Exception as exc:
                print(f"FAILED seed {seed}: {exc!r}", file=sys.stderr, flush=True)
                continue
            for row in records:
                key = (row["seed"], row["approach"], row["checkpoint_years"])
                if key not in complete:
                    out.write(json.dumps(row) + "\n")
                    rows.append(row)
                    complete.add(key)
            out.flush()
            print(f"[{i}/{len(pending)}] seed {seed} ({time.perf_counter()-start:.0f}s)", flush=True)
    summary = summarize(rows, seeds, checkpoints)
    summary["node_density"] = args.node_density
    summary_path = args.out.with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2))
    svg_path = args.out.with_suffix(".svg")
    write_svg(summary, svg_path)
    png_path = args.out.with_suffix(".png")
    render_png(summary, png_path)
    print(f"Wrote {summary_path}, {svg_path}, and {png_path}; {len(needed - complete)} checkpoint rows still missing.")


if __name__ == "__main__":
    main()
