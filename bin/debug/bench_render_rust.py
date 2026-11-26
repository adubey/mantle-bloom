"""Benchmark sequential Python, sequential Rust, and pipelined parallel Rust rendering (#1039, #1035).

Runs each benchmark arm in a fresh subprocess so `ru_maxrss` peak RSS is uncontaminated by
prior runs, and measures:
  1. Single-frame render wall-clock time (cold and warm) and peak RSS for `combined`, `biome`,
     `elevation`, and `globe` (`combined`) at 2200x1222 on Density 1 (`1/1/1`) and Density 4 (`4/4/2`).
  2. Animation wall-clock time, per-frame time, and peak RSS across `steps_per_frame` (1, 10)
     comparing sequential Python, sequential Rust, pipelined parallel Rust (2 render threads),
     and constrained parallel Rust (1 render thread).

Usage (from the repository or worktree root, after checking `compute_budget.py`):

    backend/.venv/bin/python bin/debug/bench_render_rust.py --densities 1,4
    backend/.venv/bin/python bin/debug/bench_render_rust.py --densities 1 --json /tmp/bench.json
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "backend"))


def _peak_rss_mb() -> float:
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    if sys.platform == "darwin":
        return raw / (1024.0 * 1024.0)
    return raw / 1024.0


def _build_world(seed: int, density: int, warmup_steps: int = 2):
    from app.world import generate_world, step_world

    fluid = 1 if density == 1 else 2
    world = generate_world(
        seed=seed,
        node_density=float(density),
        climate_density=float(density),
        fluid_density=float(fluid),
    )
    for _ in range(warmup_steps):
        step_world(world, 100_000.0)
    return world


def _run_single_worker(seed: int, density: int, width: int, height: int) -> dict:
    from app import _render_rust, render_image

    world = _build_world(seed, density)
    cases = [
        ("behrmann", "combined"),
        ("behrmann", "biome"),
        ("behrmann", "elevation"),
        ("globe", "combined"),
    ]
    results = {}
    for proj, view in cases:
        t0 = time.perf_counter()
        png0 = render_image.render_png(world, proj, view, width, height)
        cold_s = time.perf_counter() - t0
        warm_times = []
        for _ in range(3):
            t1 = time.perf_counter()
            png1 = render_image.render_png(world, proj, view, width, height)
            warm_times.append(time.perf_counter() - t1)
            assert len(png1) == len(png0) or len(png1) > 0
        results[f"{proj}/{view}"] = {
            "cold_s": round(cold_s, 4),
            "warm_mean_s": round(sum(warm_times) / len(warm_times), 4),
            "warm_min_s": round(min(warm_times), 4),
        }
    return {
        "rust_available": _render_rust.RUST_AVAILABLE,
        "render_threads": _render_rust.render_pool_threads(),
        "density": density,
        "peak_rss_mb": round(_peak_rss_mb(), 1),
        "views": results,
    }


def _run_anim_worker(
    seed: int,
    density: int,
    width: int,
    height: int,
    steps_per_frame: int,
    num_frames: int,
    parallel: bool,
) -> dict:
    from app import _render_rust, render_image

    world = _build_world(seed, density)
    t0 = time.perf_counter()
    msgs = list(
        render_image.stream_animation_mp4(
            world,
            "behrmann",
            "combined",
            width,
            height,
            None,
            100_000.0,
            steps_per_frame,
            num_frames,
            parallel_render=parallel,
        )
    )
    total_s = time.perf_counter() - t0
    prog = [m for m in msgs if m[0] == "progress"]
    done = msgs[-1]
    assert done[0] == "done" and done[3] is None
    return {
        "rust_available": _render_rust.RUST_AVAILABLE,
        "render_threads": _render_rust.render_pool_threads(),
        "density": density,
        "steps_per_frame": steps_per_frame,
        "num_frames": len(prog),
        "parallel": parallel,
        "total_s": round(total_s, 4),
        "per_frame_s": round(total_s / max(len(prog), 1), 4),
        "peak_rss_mb": round(_peak_rss_mb(), 1),
    }


def _spawn_arm(args_list: list[str], extra_env: dict[str, str]) -> dict:
    env = os.environ.copy()
    env.update(extra_env)
    cmd = [sys.executable, str(Path(__file__).resolve()), *args_list]
    proc = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True, check=True)
    return json.loads(proc.stdout.strip().splitlines()[-1])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--densities", default="1,4")
    ap.add_argument("--width", type=int, default=2200)
    ap.add_argument("--height", type=int, default=1222)
    ap.add_argument("--steps-per-frame", default="1,10")
    ap.add_argument("--num-frames", type=int, default=4)
    ap.add_argument("--json", type=Path, default=None)
    # Internal worker flags for isolated RSS measurement.
    ap.add_argument("--worker", choices=["single", "anim"], default=None)
    ap.add_argument("--worker-density", type=int, default=1)
    ap.add_argument("--worker-spf", type=int, default=1)
    ap.add_argument("--worker-parallel", action="store_true")
    args = ap.parse_args()

    if args.worker == "single":
        print(json.dumps(_run_single_worker(args.seed, args.worker_density, args.width, args.height)))
        return 0
    if args.worker == "anim":
        print(
            json.dumps(
                _run_anim_worker(
                    args.seed,
                    args.worker_density,
                    args.width,
                    args.height,
                    args.worker_spf,
                    args.num_frames,
                    args.worker_parallel,
                )
            )
        )
        return 0

    densities = [int(d.strip()) for d in args.densities.split(",") if d.strip()]
    spfs = [int(s.strip()) for s in args.steps_per_frame.split(",") if s.strip()]

    report: dict = {"width": args.width, "height": args.height, "seed": args.seed, "single": {}, "anim": {}}
    for d in densities:
        py_res = _spawn_arm(
            ["--worker", "single", "--seed", str(args.seed), "--worker-density", str(d), "--width", str(args.width), "--height", str(args.height)],
            {"MANTLE_BLOOM_DISABLE_RUST": "1"},
        )
        rs_res = _spawn_arm(
            ["--worker", "single", "--seed", str(args.seed), "--worker-density", str(d), "--width", str(args.width), "--height", str(args.height)],
            {"MANTLE_BLOOM_DISABLE_RUST": "0", "MANTLE_BLOOM_RENDER_THREADS": "2"},
        )
        report["single"][f"d{d}"] = {"python": py_res, "rust_2t": rs_res}

        anim_rows = []
        for spf in spfs:
            arms = [
                ("py_seq", False, {"MANTLE_BLOOM_DISABLE_RUST": "1"}),
                ("rust_seq_2t", False, {"MANTLE_BLOOM_DISABLE_RUST": "0", "MANTLE_BLOOM_RENDER_THREADS": "2"}),
                ("rust_par_2t", True, {"MANTLE_BLOOM_DISABLE_RUST": "0", "MANTLE_BLOOM_RENDER_THREADS": "2"}),
                ("rust_par_1t", True, {"MANTLE_BLOOM_DISABLE_RUST": "0", "MANTLE_BLOOM_RENDER_THREADS": "1"}),
            ]
            for label, par, env in arms:
                cmd_args = [
                    "--worker", "anim",
                    "--seed", str(args.seed),
                    "--worker-density", str(d),
                    "--width", str(args.width),
                    "--height", str(args.height),
                    "--worker-spf", str(spf),
                    "--num-frames", str(args.num_frames),
                ]
                if par:
                    cmd_args.append("--worker-parallel")
                res = _spawn_arm(cmd_args, env)
                res["label"] = label
                anim_rows.append(res)
        report["anim"][f"d{d}"] = anim_rows

    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
