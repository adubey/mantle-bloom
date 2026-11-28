"""Where does `cKDTree` time go, per call site, during stepping and rendering? (#1042)

Swaps an instrumented `cKDTree` subclass into every imported `app.*` module (each binds its own
`cKDTree` name at import time) and into `scipy.spatial`, and times the shared
`app.spatial_index.SpatialIndex` (rows prefixed `si.`) where it exists, then runs `--steps` `step_world` calls
after `--warmup` unmeasured ones, followed by `--renders` `combined` map renders after one more
step. Every tree construction is charged to the call site that built it, and every query to the
call site that issued it, so a cached tree queried from five places shows up as five query rows.

Reported per call site (`file:line:function`), split into the step and render scenarios:

  * builds: count, total ms, mean points per tree;
  * `query` (k-NN), `query_ball_point`, `query_pairs` and others: count, total ms, query points,
    and for `query` the largest `k`, whether a finite `distance_upper_bound` was given, and
    whether the result indices were used (always assumed: we cannot see that from here).

This is the measurement that picks which call sites the shared Rust spatial index (#1042, step 5
of #1046) should take over first. It records timings only; it does not change any result.

Usage (from the repository root, with this worktree's own venv; follow the compute-budget skill):

    backend/.venv/bin/python bin/debug/profile_kdtree_sites.py --seed 0 --steps 3 --renders 2 \
        --json /path/outside/repo/kdtree-young.json
    backend/.venv/bin/python bin/debug/profile_kdtree_sites.py --load <save>.bloom --steps 3

Options: --seed, --node-density/--climate-density/--fluid-density (UI defaults 4/4/2),
--step-years (default 100000), --warmup (default 1), --steps (default 3), --renders (default 2),
--load (saved world instead of a generated one), --wind-model (default diagnostic), --top (rows
printed per scenario, default 25), --json (write every row).
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import scipy.spatial  # noqa: E402
from scipy.spatial import cKDTree as _ScipyKDTree  # noqa: E402

_BACKEND = str(Path(__file__).resolve().parents[2] / "backend")
_STATS: dict[tuple[str, str, str], dict] = defaultdict(
    lambda: {"count": 0, "ms": 0.0, "points": 0, "max_k": 0, "bounded": 0, "ball_return_length": 0}
)
_SCENARIO = ["setup"]


def _site() -> str:
    """The innermost `app/` frame outside this module, as `path:line:function`."""
    frame = sys._getframe(2)
    while frame is not None:
        name = frame.f_code.co_filename
        if name.startswith(_BACKEND) and "/app/" in name:
            rel = name.split("/backend/", 1)[1]
            return f"{rel}:{frame.f_lineno}:{frame.f_code.co_name}"
        frame = frame.f_back
    return "<outside app>"


def _record(kind: str, site: str, ms: float, points: int, **extra) -> None:
    row = _STATS[(_SCENARIO[0], site, kind)]
    row["count"] += 1
    row["ms"] += ms
    row["points"] += points
    row["max_k"] = max(row["max_k"], int(extra.get("k", 0)))
    row["bounded"] += int(bool(extra.get("bounded", False)))
    row["ball_return_length"] += int(bool(extra.get("return_length", False)))


def _n_queries(x) -> int:
    arr = np.asarray(x)
    return 1 if arr.ndim == 1 else int(np.prod(arr.shape[:-1]))


class InstrumentedKDTree(_ScipyKDTree):
    def __init__(self, data, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        super().__init__(data, *args, **kwargs)
        _record("build", site, 1e3 * (time.perf_counter() - start), int(self.n))

    def query(self, x, k=1, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().query(x, k, *args, **kwargs)
        kmax = int(np.max(k)) if np.ndim(k) else int(k)
        bounded = np.isfinite(kwargs.get("distance_upper_bound", args[2] if len(args) > 2 else np.inf))
        _record("query", site, 1e3 * (time.perf_counter() - start), _n_queries(x), k=kmax, bounded=bounded)
        return out

    def query_ball_point(self, x, r, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().query_ball_point(x, r, *args, **kwargs)
        _record("query_ball_point", site, 1e3 * (time.perf_counter() - start), _n_queries(x),
                return_length=kwargs.get("return_length", False))
        return out

    def query_pairs(self, r, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().query_pairs(r, *args, **kwargs)
        _record("query_pairs", site, 1e3 * (time.perf_counter() - start), int(self.n))
        return out

    def query_ball_tree(self, other, r, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().query_ball_tree(other, r, *args, **kwargs)
        _record("query_ball_tree", site, 1e3 * (time.perf_counter() - start), int(other.n))
        return out

    def count_neighbors(self, other, r, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().count_neighbors(other, r, *args, **kwargs)
        _record("count_neighbors", site, 1e3 * (time.perf_counter() - start), int(other.n))
        return out

    def sparse_distance_matrix(self, other, max_distance, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = super().sparse_distance_matrix(other, max_distance, *args, **kwargs)
        _record("sparse_distance_matrix", site, 1e3 * (time.perf_counter() - start), int(other.n))
        return out


def _instrument_spatial_index() -> bool:
    """Time `app.spatial_index.SpatialIndex` (#1042) the same way, as `si.*` rows. Patching the
    class's own methods reaches every module that imported the name."""
    try:
        from app.spatial_index import SpatialIndex
    except ImportError:
        return False
    original = {name: getattr(SpatialIndex, name) for name in ("__init__", "query", "query_ball_point", "query_ball_point_csr", "query_pairs")}

    def init(self, points, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        original["__init__"](self, points, *args, **kwargs)
        _record("si.build", site, 1e3 * (time.perf_counter() - start), int(self.n))

    def query(self, x, k=1, *args, **kwargs):
        site = _site()
        start = time.perf_counter()
        out = original["query"](self, x, k, *args, **kwargs)
        bounded = np.isfinite(kwargs.get("distance_upper_bound", args[2] if len(args) > 2 else np.inf))
        _record("si.query", site, 1e3 * (time.perf_counter() - start), _n_queries(x), k=int(k), bounded=bounded)
        return out

    def timed(name):
        def method(self, x, *args, **kwargs):
            site = _site()
            start = time.perf_counter()
            out = original[name](self, x, *args, **kwargs)
            points = int(self.n) if name == "query_pairs" else _n_queries(x)
            _record(f"si.{name}", site, 1e3 * (time.perf_counter() - start), points,
                    return_length=kwargs.get("return_length", False))
            return out
        return method

    SpatialIndex.__init__ = init
    SpatialIndex.query = query
    for name in ("query_ball_point", "query_ball_point_csr", "query_pairs"):
        setattr(SpatialIndex, name, timed(name))
    return True


def install() -> int:
    """Rebind every imported `app.*` module's `cKDTree` name (except `app.spatial_index`'s own,
    which `SpatialIndex` rows already cover) and time `SpatialIndex`; return how many modules
    were patched."""
    import importlib
    import pkgutil

    import app

    for info in pkgutil.walk_packages(app.__path__, "app."):
        try:
            importlib.import_module(info.name)
        except Exception:  # noqa: BLE001 -- optional modules (tests-only deps) may not import
            continue
    scipy.spatial.cKDTree = InstrumentedKDTree
    patched = 0
    for name, module in list(sys.modules.items()):
        if name == "app.spatial_index":
            continue
        if (name == "app" or name.startswith("app.")) and getattr(module, "cKDTree", None) is _ScipyKDTree:
            module.cKDTree = InstrumentedKDTree
            patched += 1
    _instrument_spatial_index()
    return patched


def _rows(scenario: str) -> list[dict]:
    rows = []
    for (scen, site, kind), row in _STATS.items():
        if scen == scenario:
            rows.append({"site": site, "kind": kind, **row})
    rows.sort(key=lambda r: -r["ms"])
    return rows


def _print(scenario: str, rows: list[dict], divisor: float, top: int) -> None:
    total = sum(r["ms"] for r in rows) / divisor
    build = sum(r["ms"] for r in rows if r["kind"].endswith("build")) / divisor
    native = sum(r["ms"] for r in rows if r["kind"].startswith("si.")) / divisor
    print(f"\n== {scenario}: {total:.1f} ms in spatial queries per call ({build:.1f} build, {total - build:.1f} query; "
          f"{native:.1f} in SpatialIndex, {total - native:.1f} in bare cKDTree)")
    print(f"{'ms/call':>9} {'count':>7} {'points':>11} {'k':>3} {'bnd':>4}  kind               site")
    for r in rows[:top]:
        print(
            f"{r['ms'] / divisor:9.1f} {r['count'] / divisor:7.1f} {r['points'] / max(r['count'], 1):11.0f} "
            f"{r['max_k']:3d} {r['bounded']:4d}  {r['kind']:<18} {r['site']}"
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for flag, kind, default in (
        ("--seed", int, 0), ("--warmup", int, 1), ("--steps", int, 3), ("--renders", int, 2),
        ("--node-density", float, 4.0), ("--climate-density", float, 4.0), ("--fluid-density", float, 2.0),
        ("--step-years", float, 100_000.0),
    ):
        ap.add_argument(flag, type=kind, default=default)
    ap.add_argument("--load", type=Path)
    ap.add_argument("--wind-model", choices=("diagnostic", "cfd"), default="diagnostic")
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    patched = install()
    from app import persistence, render_image, world as world_mod

    if args.load:
        with args.load.open("rb") as handle:
            world = persistence.load_world_file(handle)
    else:
        world = world_mod.generate_world(
            args.seed, node_density=args.node_density, climate_density=args.climate_density,
            fluid_density=args.fluid_density,
        )
    world.wind_model = args.wind_model
    print(f"patched {patched} modules; {len(world.plates)} plates, {world.elapsed_years / 1e6:.1f} Myr", file=sys.stderr)

    for _ in range(args.warmup):
        world_mod.step_world(world, args.step_years)
    _SCENARIO[0] = "step"
    step_s = []
    for _ in range(args.steps):
        start = time.perf_counter()
        world_mod.step_world(world, args.step_years)
        step_s.append(time.perf_counter() - start)
    _SCENARIO[0] = "setup"
    world_mod.step_world(world, args.step_years)  # a step drops the render caches, as in the UI
    _SCENARIO[0] = "render"
    render_s = []
    for _ in range(args.renders):
        start = time.perf_counter()
        render_image.render_png(world, "eckert4", "combined", 1100, 611)
        render_s.append(time.perf_counter() - start)

    out = {"args": {k: str(v) for k, v in vars(args).items()}, "step_seconds": step_s, "render_seconds": render_s}
    for scenario, divisor in (("step", max(args.steps, 1)), ("render", max(args.renders, 1))):
        rows = _rows(scenario)
        _print(scenario, rows, float(divisor), args.top)
        out[scenario] = {"divisor": divisor, "rows": rows}
    print("\nstep s: " + ", ".join(f"{s:.2f}" for s in step_s) + "; render s: " + ", ".join(f"{s:.2f}" for s in render_s))
    if args.json:
        args.json.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
