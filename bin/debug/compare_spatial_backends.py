"""Does `SpatialIndex` answer every query a real step makes exactly as `cKDTree` would? (#1042)

Wraps `app.spatial_index.SpatialIndex`'s query methods so each call is also answered by a plain
`scipy.spatial.cKDTree` over the same points, then steps a world `--steps` times and reports, per
call site (`file:line:function`), how many calls and query rows disagreed and why:

  * `query`: rows whose indices differ, split into rows where every differing slot is an exact
    distance tie (same distance, different point: the documented tie-breaking difference) and
    rows where the distances themselves differ (a real discrepancy);
  * `query_ball_point` / `query_ball_point_csr`: rows whose sorted index sets differ;
  * `query_pairs`: calls whose pair sets differ.

The step itself runs on `SpatialIndex`'s answers, so the world evolves as on the branch.

Usage (from the repository root, with this worktree's own venv; follow the compute-budget skill):

    backend/.venv/bin/python bin/debug/compare_spatial_backends.py --load <save>.bloom --steps 1
    backend/.venv/bin/python bin/debug/compare_spatial_backends.py --seed 0 --steps 1 --json out.json

Options: --seed, --load, --steps (default 1), --step-years (default 100000), --wind-model
(default diagnostic), --reference (`ckdtree`, default: a plain `cKDTree`, so documented tie
differences show up as `ties`; `fallback`: `SpatialIndex`'s own scipy fallback, which must
agree exactly, ties and ordering included), --json (write every row).
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import spatial_index  # noqa: E402
from app.spatial_index import SpatialIndex  # noqa: E402

_BACKEND = str(Path(__file__).resolve().parents[2] / "backend")
STATS: dict[tuple[str, str], dict] = defaultdict(lambda: {"calls": 0, "rows": 0, "bad_calls": 0, "tie_rows": 0, "distance_rows": 0, "set_rows": 0})


def _site() -> str:
    frame = sys._getframe(2)
    while frame is not None:
        name = frame.f_code.co_filename
        if name.startswith(_BACKEND) and "/app/" in name and not name.endswith("spatial_index.py"):
            return f"{name.split('/backend/', 1)[1]}:{frame.f_lineno}:{frame.f_code.co_name}"
        frame = frame.f_back
    return "<outside app>"


def _compare_query(row, got, want):
    (gd, gi), (wd, wi) = got, want
    gd, gi, wd, wi = (np.atleast_1d(np.asarray(a)) for a in (gd, gi, wd, wi))
    gd, gi, wd, wi = (a.reshape(len(np.atleast_2d(a)) if a.ndim > 1 else a.size, -1) for a in (gd, gi, wd, wi))
    differ = np.any(gi != wi, axis=1)
    row["rows"] += len(gi)
    if not np.any(differ):
        return
    row["bad_calls"] += 1
    distance = np.any(gd != wd, axis=1)
    row["distance_rows"] += int(np.count_nonzero(differ & distance))
    row["tie_rows"] += int(np.count_nonzero(differ & ~distance))


def _fallback_index(points):
    """A `SpatialIndex` on its scipy fallback, over the same points."""
    saved = spatial_index._KdIndex
    spatial_index._KdIndex = None
    try:
        return SpatialIndex(points)
    finally:
        spatial_index._KdIndex = saved


def install(reference: str = "ckdtree") -> None:
    original = {name: getattr(SpatialIndex, name) for name in ("query", "query_ball_point", "query_ball_point_csr", "query_pairs")}
    if reference == "fallback":
        _install_fallback_reference(original)
        return

    def query(self, x, k=1, *args, **kwargs):
        got = original["query"](self, x, k, *args, **kwargs)
        kwargs.pop("workers", None)
        want = cKDTree(self.data).query(x, k, *args, **kwargs)
        row = STATS[(_site(), "query")]
        row["calls"] += 1
        _compare_query(row, got, want)
        return got

    def ball(self, x, r, *args, **kwargs):
        got = original["query_ball_point"](self, x, r, *args, **kwargs)
        kwargs.pop("workers", None)
        want = cKDTree(self.data).query_ball_point(x, r, *args, **kwargs)
        row = STATS[(_site(), "query_ball_point")]
        row["calls"] += 1
        if kwargs.get("return_length"):
            g, w = np.atleast_1d(got), np.atleast_1d(want)
            bad = int(np.count_nonzero(g != w))
        else:
            g = [got] if isinstance(got, list) else list(got)
            w = [want] if isinstance(want, list) else list(want)
            bad = sum(sorted(a) != sorted(b) for a, b in zip(g, w))
        row["rows"] += len(g)
        row["set_rows"] += bad
        row["bad_calls"] += int(bad > 0)
        return got

    def csr(self, x, r):
        indptr, indices = original["query_ball_point_csr"](self, x, r)
        rows = np.asarray(x, dtype=float).reshape(-1, 3)
        radii = np.broadcast_to(np.asarray(r, dtype=float), (len(rows),))
        tree = cKDTree(self.data)
        row = STATS[(_site(), "query_ball_point_csr")]
        row["calls"] += 1
        row["rows"] += len(rows)
        bad = sum(sorted(tree.query_ball_point(rows[i], radii[i])) != indices[indptr[i]:indptr[i + 1]].tolist() for i in range(len(rows)))
        row["set_rows"] += bad
        row["bad_calls"] += int(bad > 0)
        return indptr, indices

    def pairs(self, r, *args, **kwargs):
        got = original["query_pairs"](self, r, *args, **kwargs)
        row = STATS[(_site(), "query_pairs")]
        row["calls"] += 1
        want = cKDTree(self.data).query_pairs(r)
        g = got if isinstance(got, set) else set(map(tuple, np.asarray(got).tolist()))
        row["bad_calls"] += int(g != want)
        return got

    SpatialIndex.query = query
    SpatialIndex.query_ball_point = ball
    SpatialIndex.query_ball_point_csr = csr
    SpatialIndex.query_pairs = pairs


def _same(a, b) -> bool:
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, set) or isinstance(b, set):
        return a == b
    if isinstance(a, np.ndarray) and a.dtype == object:
        return isinstance(b, np.ndarray) and a.shape == b.shape and all(list(x) == list(y) for x, y in zip(a.ravel(), b.ravel()))
    return np.array_equal(np.asarray(a), np.asarray(b)) and np.asarray(a).dtype == np.asarray(b).dtype


def _install_fallback_reference(original) -> None:
    """Answer every call on the native index and on a fallback index too; count calls whose
    results differ in any way (values, order, dtype or shape)."""
    for name, fn in original.items():
        def method(self, *args, _name=name, _fn=fn, **kwargs):
            got = _fn(self, *args, **kwargs)
            want = _fn(_fallback_index(self.data), *args, **kwargs)
            row = STATS[(_site(), _name)]
            row["calls"] += 1
            if not _same(got, want):
                row["bad_calls"] += 1
                if row["bad_calls"] <= 3:
                    print(f"MISMATCH {_name} at {_site()}: args={[np.shape(a) for a in args]} kwargs={kwargs}", file=sys.stderr)
            return got
        setattr(SpatialIndex, name, method)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--load", type=Path)
    ap.add_argument("--steps", type=int, default=1)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--wind-model", choices=("diagnostic", "cfd"), default="diagnostic")
    ap.add_argument("--reference", choices=("ckdtree", "fallback"), default="ckdtree")
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()
    if not spatial_index.RUST_AVAILABLE:
        raise SystemExit("mantle_bloom_rust.KdIndex unavailable: nothing to compare")
    from app import persistence, world as world_mod

    if args.load:
        with args.load.open("rb") as handle:
            world = persistence.load_world_file(handle)
    else:
        world = world_mod.generate_world(args.seed)
    world.wind_model = args.wind_model
    install(args.reference)
    for _ in range(args.steps):
        world_mod.step_world(world, args.step_years)

    rows = [{"site": site, "kind": kind, **row} for (site, kind), row in STATS.items()]
    rows.sort(key=lambda r: (-(r["tie_rows"] + r["distance_rows"] + r["set_rows"]), r["site"]))
    print(f"{'calls':>7} {'rows':>9} {'badcalls':>8} {'ties':>6} {'dist':>6} {'sets':>6}  kind                 site")
    for r in rows:
        print(f"{r['calls']:7d} {r['rows']:9d} {r['bad_calls']:8d} {r['tie_rows']:6d} {r['distance_rows']:6d} {r['set_rows']:6d}  {r['kind']:<20} {r['site']}")
    if args.json:
        args.json.write_text(json.dumps(rows, indent=1))


if __name__ == "__main__":
    main()
