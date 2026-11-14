"""Profile RFC 0018's candidate search on a mature world, and how often it runs where it can't pay off (#848).

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/profile_rfc0018_candidate_search.py \
        --save ~/mantle-bloom-captures/issue848/world_100myr.mbworld \
        --g0-script /path/to/measure_rfc0018_g0.py \
        --out ~/mantle-bloom-captures/issue848/profile.jsonl \
        [--warmup-steps 5] [--steps 15] [--mappings inverse,direct]

``--g0-script`` is ``bin/debug/measure_rfc0018_g0.py`` from origin/codex/issue-797-evidence-rfc0018-g0
(the frozen G0/G2 candidate method, not on main); it is loaded by path so it isn't copied here.
Its ``choose_cut`` is called unchanged; only timers and counters are wrapped around its helpers.
``inverse`` is that script's ``conductance = 1 / edge_cost``; ``direct`` swaps in
``conductance = edge_cost`` (the RFC 0018 revision, #850), changing nothing else.

The world steps with the production great-circle splitter (a baseline-arm history) and
``simulate_climate_biomes`` off. After ``--warmup-steps`` unmeasured steps, every split
opportunity (the g0 script's ``_is_opportunity``) is measured from the same pre-split state,
as the g0 runner does. Each row records wall time per helper, sparse-solve counts, how many
variants used all 8 fit iterations (the cap, or a convergence on the eighth), and when the
plate last ran ``quad_tectonics._open_rift``. That last field is the rift-candidate marker proxy:
``steps_since_open_rift`` is None if the plate hasn't been stretched since warm-up began.

Timers are inclusive, so helper times overlap (``_fiedler_vector`` contains ``eigsh``).
Profiling overhead stays out of any acceptance timing (RFC 0018 G3): use these numbers to
divide cost, not to state a runtime ratio.
"""

import argparse
import importlib.util
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import merge_split, persistence, quad_tectonics  # noqa: E402
from app import world as world_mod  # noqa: E402

TIMED = ["_fiedler_vector", "_node_resistance", "_opening", "_trace_stats", "_largest_component_share",
         "_shared_edges", "eigsh", "kmeans2", "_fit_euler_pole"]


def load_g0(path: Path):
    spec = importlib.util.spec_from_file_location("g0_candidate_search", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Meter:
    """Inclusive wall time and call count per wrapped function, reset per measurement."""

    def __init__(self):
        self.seconds = defaultdict(float)
        self.calls = defaultdict(int)
        self.solve_sizes = []

    def wrap(self, name, fn):
        def timed(*args, **kwargs):
            started = time.perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                self.seconds[name] += time.perf_counter() - started
                self.calls[name] += 1
                if name == "_fiedler_vector":
                    self.solve_sizes.append(int(args[0]))

        return timed

    def reset(self):
        self.seconds.clear()
        self.calls.clear()
        self.solve_sizes.clear()


def install(g0, meter, mapping):
    """Wrap g0's helpers in place; returns a restore function."""
    saved = {}
    original_fiedler = g0._fiedler_vector

    def direct_fiedler(n, i, j, edge_cost):
        # The inverse code builds conductance = 1 / edge_cost; passing 1 / edge_cost makes it edge_cost.
        return original_fiedler(n, i, j, 1.0 / np.maximum(edge_cost, 1e-10))

    pairs = {
        "_fiedler_vector": direct_fiedler if mapping == "direct" else original_fiedler,
        "_node_resistance": g0._node_resistance, "_opening": g0._opening, "_trace_stats": g0._trace_stats,
        "_largest_component_share": g0._largest_component_share, "_shared_edges": g0._shared_edges,
        "eigsh": g0.eigsh, "kmeans2": g0.kmeans2, "_fit_euler_pole": g0.mantle.fit_euler_pole,
    }
    for name, fn in pairs.items():
        target = g0.mantle if name == "_fit_euler_pole" else g0
        attr = "fit_euler_pole" if name == "_fit_euler_pole" else name
        saved[(id(target), attr)] = (target, attr, getattr(target, attr))
        setattr(target, attr, meter.wrap(name, fn))

    def restore():
        for target, attr, fn in saved.values():
            setattr(target, attr, fn)

    return restore


def measure(g0, meter, world, plate, mapping):
    meter.reset()
    restore = install(g0, meter, mapping)
    started = time.perf_counter()
    row = {"mapping": mapping, "status": "ok"}
    try:
        _, chosen, variants, _ = g0.choose_cut(world, plate)
        row["variants"] = len(variants)
        row["fit_iterations"] = [int(v["candidate_fit_iterations"]) for v in variants]
        row["variants_at_iteration_cap"] = sum(1 for v in variants if v["candidate_fit_iterations"] >= 8)
        row["usable_variants"] = sum(1 for v in variants if v["usable_candidate"])
        row["selected_accepted"] = bool(chosen.get("accepted"))
    except Exception as exc:  # the g0 method raises ValueError when no candidate exists
        row["status"] = f"{type(exc).__name__}: {exc}"
    finally:
        restore()
    row["total_s"] = time.perf_counter() - started
    row["seconds"] = {k: round(v, 4) for k, v in meter.seconds.items()}
    row["calls"] = dict(meter.calls)
    row["sparse_solves"] = meter.calls["_fiedler_vector"]
    row["solve_nodes"] = meter.solve_sizes[0] if meter.solve_sizes else None
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--save", required=True, type=Path)
    parser.add_argument("--g0-script", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--warmup-steps", type=int, default=5)
    parser.add_argument("--steps", type=int, default=15)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--mappings", default="inverse,direct")
    args = parser.parse_args()
    mappings = [m for m in args.mappings.split(",") if m]
    g0 = load_g0(args.g0_script.expanduser())
    with open(args.save.expanduser(), "rb") as handle:
        world = persistence.load_world_file(handle)
    world.simulate_climate_biomes = False
    meter = Meter()
    step_index = 0
    last_open_rift = {}  # plate_id -> step index of its latest _open_rift call
    original_open_rift = quad_tectonics._open_rift
    original_split = merge_split.maybe_split_plate
    out = args.out.expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    handle = open(out, "w")
    diagnosed = {"step": None}

    def open_rift_spy(plate, *a, **kw):
        last_open_rift[plate.plate_id] = step_index
        return original_open_rift(plate, *a, **kw)

    def split_spy(current_world, plate, step_years=0.0):
        if step_index > args.warmup_steps and diagnosed["step"] != step_index:
            diagnosed["step"] = step_index
            for probe in list(current_world.plates):
                if not g0._is_opportunity(current_world, probe):
                    continue
                since = last_open_rift.get(probe.plate_id)
                base = {
                    "step": step_index - args.warmup_steps, "elapsed_myr": current_world.elapsed_years / 1e6,
                    "plate_id": int(probe.plate_id), "nodes": int(probe.node_count()),
                    "steps_since_open_rift": None if since is None else step_index - since,
                }
                for mapping in mappings:
                    handle.write(json.dumps({**base, **measure(g0, meter, current_world, probe, mapping)}) + "\n")
                    handle.flush()
        return original_split(current_world, plate, step_years)

    quad_tectonics._open_rift = open_rift_spy
    merge_split.maybe_split_plate = split_spy
    try:
        for _ in range(args.warmup_steps + args.steps):
            step_index += 1
            world_mod.step_world(world, args.step_myr * 1e6)
            print(json.dumps({"step": step_index, "plates": len(world.plates)}), flush=True)
    finally:
        quad_tectonics._open_rift = original_open_rift
        merge_split.maybe_split_plate = original_split
        handle.close()


if __name__ == "__main__":
    main()
