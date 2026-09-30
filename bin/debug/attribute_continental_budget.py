#!/usr/bin/env python3
"""Issue #249: attribute the continental-crust budget to step phases on line vs quad worlds.

The parity campaign found quad worlds losing continental crust (area-weighted continental Hc
volume) much faster than line worlds from the same seed. This steps one world per surface
and, around every top-level `step_world` phase (the same buckets as `surface_parity`'s
timers, each plate's `deform` call included), measures the whole world's

- continental Hc volume, sum(Hc * area) over nodes whose effective crust type is continental;
- continental area, the same sum without Hc -- its change separates material that changed
  type or was created/destroyed from thickness change;
- total Hc volume.

Areas are each node's actual area (`SurfaceNodes.area_m2`), as in the parity harness. Finer
buckets inside `deform` and `topology` (`SUB_PHASES`) are reported alongside, nested.

    cd backend
    .venv/bin/python ../bin/debug/attribute_continental_budget.py --seed 1 --steps 120 --out ../analysis/issue249-campaign/continental-budget-seed1.json
"""

from __future__ import annotations

import argparse
import functools
import inspect
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import elevation_lines, lithosphere_plate, merge_split, quad_tectonics, surface_parity as sp, world as world_mod  # noqa: E402

# Finer buckets inside `deform` and `topology`, measured independently of the top-level ones
# (so they nest inside them). Column physics is shared by both engines; boundary moves are not.
SUB_PHASES = (
    (lithosphere_plate, "deform_columns", "columns (lines)"),
    (quad_tectonics, "deform_columns", "columns (quad)"),
    (quad_tectonics, "_retreat", "quad retreat"),
    (quad_tectonics, "_advance", "quad advance"),
    (merge_split, "remove_defunct_plates", "remove defunct plates"),
    (merge_split, "defragment_plates", "defragment"),
    (merge_split, "relattice_continental_plates", "relattice"),
    (merge_split, "merge_plates", "merge"),
    (merge_split, "maybe_split_plate", "split"),
)


# Node areas depend only on a plate's local node set, so they are cached per (plate object,
# topology revision); line plates estimate them with a KD-tree per call otherwise.
_AREAS: dict[int, tuple[object, int, np.ndarray]] = {}


def _areas(plate) -> np.ndarray:
    cached = _AREAS.get(id(plate))
    if cached is not None and cached[0] is plate and cached[1] == plate.topology_revision:
        return cached[2]
    areas = plate.surface_nodes().area_m2
    _AREAS[id(plate)] = (plate, plate.topology_revision, areas)
    return areas


def budget(world) -> np.ndarray:
    """(continental Hc volume km^3, continental area km^2, total Hc volume km^3)."""
    out = np.zeros(3)
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        hc = plate.collect("crustal_thickness_m")
        area = _areas(plate)
        continental = elevation_lines.effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        volume = hc * area
        out += [volume[continental].sum() / 1e9, area[continental].sum() / 1e6, volume.sum() / 1e9]
    return out


class Attribution:
    def __init__(self, world_ref: list) -> None:
        self.world_ref = world_ref
        self.delta: dict[str, np.ndarray] = defaultdict(lambda: np.zeros(3))
        self.depth = 0

    def wrap(self, original, bucket: str):
        attribution = self

        def measured(call):
            if attribution.depth or attribution.world_ref[0] is None:
                return call()
            attribution.depth += 1
            before = budget(attribution.world_ref[0])
            try:
                return call()
            finally:
                attribution.delta[bucket] += budget(attribution.world_ref[0]) - before
                attribution.depth -= 1

        if inspect.isgeneratorfunction(original):

            @functools.wraps(original)
            def generator(*args, **kwargs):
                return (yield from measured(lambda: list(original(*args, **kwargs))))

            return generator

        @functools.wraps(original)
        def plain(*args, **kwargs):
            return measured(lambda: original(*args, **kwargs))

        return plain


def run(seed: int, surface: str, steps: int, density: float, every: int) -> dict:
    world_ref = [None]
    attribution, sub = Attribution(world_ref), Attribution(world_ref)
    patches = []
    for recorder, phases in ((attribution, [p for p in sp._PHASES if p[2] != "record_stats"]), (sub, SUB_PHASES)):
        for owner, attribute, bucket in phases:
            original = owner.__dict__[attribute] if isinstance(owner, type) else getattr(owner, attribute)
            patches.append((owner, attribute, original))
            setattr(owner, attribute, recorder.wrap(getattr(owner, attribute), bucket))
    try:
        world = world_mod.generate_world(seed=seed, node_density=density, surface=surface)
        world_ref[0] = world
        start = budget(world)
        rows = []
        for step in range(1, steps + 1):
            before_step = budget(world)
            world_mod.step_world(world, 1e6)
            attribution.delta["step_total"] += budget(world) - before_step
            if step % every == 0 or step == steps:
                rows.append(
                    {
                        "step": step,
                        "budget": budget(world).tolist(),
                        "cumulative": {k: v.tolist() for k, v in sorted(attribution.delta.items())},
                        "cumulative_sub": {k: v.tolist() for k, v in sorted(sub.delta.items())},
                    }
                )
                shares = {k: round(float(v[0] / start[0]), 4) for k, v in sorted({**attribution.delta, **sub.delta}.items()) if abs(v[0]) > 1e-4 * start[0]}
                print(surface, step, shares, flush=True)
    finally:
        for owner, attribute, original in reversed(patches):
            setattr(owner, attribute, original)
    return {"seed": seed, "surface": surface, "density": density, "start": start.tolist(), "rows": rows}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--steps", type=int, default=120)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--every", type=int, default=10)
    parser.add_argument("--surfaces", default="lines,quad")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    report = {surface: run(args.seed, surface, args.steps, args.density, args.every) for surface in args.surfaces.split(",")}
    if args.out:
        args.out.write_text(json.dumps(sp.sig(report), indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
