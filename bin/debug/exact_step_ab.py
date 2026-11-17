"""Do two versions of the code step a world to exactly the same state? (#147)

For changes that must not move results (a faster kernel, a cache). Run `step` once per
version, one after the other, then `compare` the two exact saves:

    git archive origin/main backend/app | tar -x -C <dir-outside-repo>/base
    backend/.venv/bin/python bin/debug/exact_step_ab.py step --backend <dir>/base/backend \
        --load ~/Downloads/<save>.mbworld --steps 3 --out <dir>/old.mbworld
    backend/.venv/bin/python bin/debug/exact_step_ab.py step \
        --load ~/Downloads/<save>.mbworld --steps 3 --out <dir>/new.mbworld
    backend/.venv/bin/python bin/debug/exact_step_ab.py compare <dir>/old.mbworld <dir>/new.mbworld

`step` loads `--load` (or generates `--seed` at `--density` when omitted), steps it
`--steps` times by `--step-years`, prints each step's wall time and sea level, and writes
an exact save to `--out`. `--backend` picks the `backend/` directory whose `app` package
runs (default: this checkout's), so the old arm needs no second virtualenv.

`compare` walks the two manifests together and compares every array and scalar. Fields
are matched by name: the numbering of a save's `.npy` members follows dict order, which
varies from one process to the next, so comparing members by file name reports hundreds
of false differences. Wall-clock timers stored in the world (`*_seconds`) are skipped.
Exits 1 when anything else differs.

Follow the compute-budget skill first, and run the arms back to back on an otherwise idle
machine if the step times matter.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
import zipfile
from pathlib import Path

import numpy as np


def _step(args: argparse.Namespace) -> int:
    backend = args.backend or Path(__file__).resolve().parents[2] / "backend"
    sys.path.insert(0, str(backend))
    import app
    from app import persistence
    from app.world import generate_world, step_world

    if args.load:
        with args.load.expanduser().open("rb") as handle:
            world = persistence.load_world_file(handle)
    else:
        world = generate_world(
            seed=args.seed, node_density=args.density, climate_density=args.density, fluid_density=args.density
        )
    print(f"app: {Path(app.__file__).parent}")
    for index in range(args.steps):
        start = time.perf_counter()
        step_world(world, args.step_years)
        print(f"step {index}: {time.perf_counter() - start:.1f} s, sea level {float(world.sea_level_m)!r}", flush=True)
    args.out.write_bytes(persistence.save_exact_bytes(world))
    return 0


def _by_name(node):
    """A manifest dict stored as `[key, value]` pairs, as a real dict, so the walk pairs
    entries by key rather than by position."""
    if isinstance(node, dict) and set(node) == {"__dict__"}:
        pairs = node["__dict__"]
        if all(isinstance(pair, list) and len(pair) == 2 for pair in pairs):
            return {json.dumps(key, sort_keys=True): value for key, value in pairs}
    return node


def _array(archive: zipfile.ZipFile, ref: str) -> np.ndarray:
    name = next(n for n in (f"a/{ref}.npy", f"{ref}.npy", ref) if n in archive.NameToInfo)
    return np.load(io.BytesIO(archive.read(name)), allow_pickle=True)


class _Comparison:
    def __init__(self, old: zipfile.ZipFile, new: zipfile.ZipFile, show: int):
        self.old, self.new, self.show = old, new, show
        self.arrays = self.scalars = self.skipped = self.differing = 0

    def _differs(self, path: str, why: str) -> None:
        self.differing += 1
        if self.differing <= self.show:
            print(f"DIFF {path}: {why}")

    def _arrays(self, a: dict, b: dict, path: str) -> None:
        x, y = _array(self.old, a["__npy__"]), _array(self.new, b["__npy__"])
        self.arrays += 1
        if x.shape != y.shape or x.dtype != y.dtype:
            self._differs(path, f"{x.dtype}{x.shape} vs {y.dtype}{y.shape}")
        elif not np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"):
            self._differs(path, f"{int(np.count_nonzero(x != y))} of {x.size} elements")

    def walk(self, a, b, path: str = "") -> None:
        a, b = _by_name(a), _by_name(b)
        if isinstance(a, dict) and isinstance(b, dict):
            if set(a) == {"__npy__"} and set(b) == {"__npy__"}:
                self._arrays(a, b, path)
            elif set(a) != set(b):
                self._differs(path, f"keys differ: {sorted(set(a) ^ set(b))}")
            else:
                for key in a:
                    self.walk(a[key], b[key], f"{path}/{key}")
        elif isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                self._differs(path, f"length {len(a)} vs {len(b)}")
            else:
                for index, (x, y) in enumerate(zip(a, b)):
                    self.walk(x, y, f"{path}[{index}]")
        elif path.rstrip('"').endswith("_seconds"):
            self.skipped += 1
        else:
            self.scalars += 1
            both_nan = isinstance(a, float) and isinstance(b, float) and a != a and b != b
            if a != b and not both_nan:
                self._differs(path, f"{a!r} vs {b!r}")


def _compare(args: argparse.Namespace) -> int:
    with zipfile.ZipFile(args.old) as old, zipfile.ZipFile(args.new) as new:
        comparison = _Comparison(old, new, args.show)
        comparison.walk(json.loads(old.read("manifest.json")), json.loads(new.read("manifest.json")))
        members = (len(old.namelist()), len(new.namelist()))
    print(
        f"{comparison.arrays} arrays and {comparison.scalars} scalars compared, {comparison.differing} differ "
        f"({comparison.skipped} wall-clock timers skipped; members {members[0]} vs {members[1]})"
    )
    return 1 if comparison.differing or members[0] != members[1] else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="command", required=True)
    step = sub.add_parser("step", help="step a world and write an exact save")
    step.add_argument("--backend", type=Path, help="backend/ directory whose app package to run")
    step.add_argument("--load", type=Path, help="saved world to start from")
    step.add_argument("--seed", type=int, default=0)
    step.add_argument("--density", type=float, default=1.0, help="node, climate and fluid density of a generated world")
    step.add_argument("--steps", type=int, default=3)
    step.add_argument("--step-years", type=float, default=100_000)
    step.add_argument("--out", type=Path, required=True)
    step.set_defaults(run=_step)
    compare = sub.add_parser("compare", help="compare two exact saves field by field")
    compare.add_argument("old", type=Path)
    compare.add_argument("new", type=Path)
    compare.add_argument("--show", type=int, default=25, help="differences to print")
    compare.set_defaults(run=_compare)
    args = ap.parse_args()
    return args.run(args)


if __name__ == "__main__":
    raise SystemExit(main())
