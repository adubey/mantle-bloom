#!/usr/bin/env python3
"""Interactive explorer for the cratons-first layout (RFC 0009, #685): change the layout knobs
and see batches of worlds side by side.

Serves one page (layout_explorer.html, beside this script) on localhost. Each batch is a list
of jobs (a seed plus constant overrides): N seeds at fixed settings, one knob swept across
columns (one column of the same N seeds per value), or one knob swept down a single column at
a fixed seed. Every job runs `lithosphere_plate.cratons_first_layout` (continents, cratons,
belts, land mask and plate partition; no plates or relief) in a process pool, with the
requested module constants patched for that job only. Each map is the layout lattice
resampled onto an equirectangular grid (nearest lattice point per pixel), streamed back to
the page as it finishes, with the layout's own diagnostics.

A constant changes the rng draws after it, so the same seed under different settings is a
different world from the first knob that moves onward, not the same world perturbed.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/layout_explorer.py --port 8765
"""

from __future__ import annotations

import argparse
import base64
import json
import threading
import time
import traceback
from concurrent.futures import CancelledError, ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import import_module
from pathlib import Path

import numpy as np

PAGE = Path(__file__).with_name("layout_explorer.html")
MAP_WIDTH, MAP_HEIGHT = 480, 240
MAX_JOBS = 40

# (module, constant, kind, group, label). Kinds: "float", "int", "frange" and "irange" (a
# (low, high) pair). Every constant is read from its module at call time, so patching the
# module attribute is enough.
KNOBS = [
    ("continent_layout", "CONTINENT_COUNT_RANGE", "irange", "Continents", "Continent count (min, max)"),
    ("continent_layout", "CONTINENT_SPACING_KM", "float", "Continents", "Poisson-disc centre spacing, km"),
    ("continent_layout", "CONTINENT_SIZE_CONCENTRATION", "float", "Continents", "Dirichlet size concentration"),
    ("continent_layout", "CONTINENT_OUTLINE_NOISE_KM", "float", "Continents", "Outline noise, km"),
    ("continent_layout", "CONTINENT_GROWTH_SPEED_SIGMA", "float", "Continents", "Growth-speed noise sigma"),
    ("continent_layout", "CONTINENT_CENTRE_PULL", "float", "Continents", "Pull toward continent centre"),
    ("continent_layout", "CRATON_SHARE_OF_CONTINENT", "float", "Cratons", "Craton share of continent"),
    ("continent_layout", "CRATON_SHARE_MAX", "float", "Cratons", "Craton share cap"),
    ("continent_layout", "CRATON_CLUSTER_CONCENTRATION", "float", "Cratons", "CRP concentration (alpha)"),
    ("continent_layout", "NUCLEUS_JOIN_STEP_RANGE", "frange", "Cratons", "Join step, x craton reach"),
    ("continent_layout", "NUCLEUS_STEP_RANGE", "frange", "Cratons", "New-craton step, x core spacing"),
    ("continent_layout", "CRATON_SIZE_SCALE_RANGE", "frange", "Cratons", "Craton size scale"),
    ("continent_layout", "CRATON_CORE_RADIUS_KM", "float", "Cratons", "Core radius, km"),
    ("continent_layout", "CRATON_CORE_SPACING_KM", "float", "Cratons", "Core spacing, km"),
    ("continent_layout", "CRATON_MAX_AREA_KM2", "float", "Cratons", "Max composite area, km²"),
    ("continent_layout", "CRATON_CORE_OUTLINE_NOISE_KM", "float", "Cratons", "Core outline noise, km"),
    ("continent_layout", "CRATON_BELT_WIDTH_KM", "float", "Belts", "Belt width, km"),
    ("continent_layout", "BELT_CLOSURE_KM", "float", "Belts", "Belt closure radius, km"),
    ("continent_layout", "BELT_LENS_BELT_WIDTHS", "float", "Belts", "Belt lens, x belt width"),
    ("plate_partition", "POINT_SEED_COUNT", "int", "Plates", "Point sites"),
    ("plate_partition", "PARTITION_CANDIDATES", "int", "Plates", "Candidate partitions"),
    ("plate_partition", "MAX_PARTITION_ATTEMPTS", "int", "Plates", "Max partition draws"),
    ("plate_partition", "REFINE_MAX_SWEEPS", "int", "Plates", "Refinement sweeps"),
    ("plate_partition", "PLATE_AREA_BAND", "frange", "Plates", "Refine area band"),
    ("plate_partition", "PLATE_COST_BELT", "float", "Plates", "Boundary cost: belt"),
    ("plate_partition", "PLATE_COST_MARGIN", "float", "Plates", "Boundary cost: margin"),
    ("plate_partition", "PLATE_COST_CORE", "float", "Plates", "Boundary cost: core (follows margin)"),
    ("plate_partition", "PLATE_COST_LAND_LINK", "float", "Plates", "Boundary cost: land link"),
    ("lithosphere_plate", "MIN_AUTO_PLATES", "int", "Plates", "Auto plate count min"),
    ("lithosphere_plate", "MAX_AUTO_PLATES", "int", "Plates", "Auto plate count max"),
]
KNOB_BY_NAME = {k[1]: k for k in KNOBS}
# Knobs the backend sets from another knob at import (`PLATE_COST_CORE = PLATE_COST_MARGIN`).
# Overriding the source alone moves the derived knob too, as editing the source line would.
DERIVED_KNOBS = {"PLATE_COST_CORE": "PLATE_COST_MARGIN"}
# (low, high) knob pairs whose effective values must stay ordered.
ORDERED_KNOBS = [("MIN_AUTO_PLATES", "MAX_AUTO_PLATES")]

# Each worker's import-time knob values, captured before its first override.
_worker_defaults: dict[str, object] = {}


def _module(name: str):
    return import_module(f"app.{name}")


def _jsonable(value):
    return list(value) if isinstance(value, tuple) else value


def knob_defaults() -> list[dict]:
    """Each knob with its value on this checkout, for the page's form."""
    return [
        {"name": name, "module": module, "kind": kind, "group": group, "label": label, "default": _jsonable(getattr(_module(module), name))}
        for module, name, kind, group, label in KNOBS
    ]


def _coerce(name: str, value):
    kind = KNOB_BY_NAME[name][2]
    if kind == "float":
        return float(value)
    if kind == "int":
        return int(value)
    low, high = (float(v) for v in value) if kind == "frange" else (int(v) for v in value)
    if low > high:
        raise ValueError(f"{name}: min {low} is above max {high}")
    return (low, high)


def _patch(overrides: dict) -> None:
    """Apply `overrides` to the worker's modules, first restoring every knob to its import
    value: workers are reused across jobs."""
    for module, name, *_ in KNOBS:
        mod = _module(module)
        _worker_defaults.setdefault(name, getattr(mod, name))
        setattr(mod, name, _worker_defaults[name])
    for name, value in overrides.items():
        setattr(_module(KNOB_BY_NAME[name][0]), name, _coerce(name, value))
    for derived, source in DERIVED_KNOBS.items():
        if source in overrides and derived not in overrides:
            setattr(_module(KNOB_BY_NAME[derived][0]), derived, _coerce(source, overrides[source]))


def _check_order(overrides: dict, defaults: dict) -> None:
    for low, high in ORDERED_KNOBS:
        lo, hi = overrides.get(low, defaults[low]), overrides.get(high, defaults[high])
        if lo > hi:
            raise ValueError(f"{low} = {lo} is above {high} = {hi}")


def _pixel_xyz() -> np.ndarray:
    lon = (np.arange(MAP_WIDTH) + 0.5) / MAP_WIDTH * 2 * np.pi - np.pi
    lat = np.pi / 2 - (np.arange(MAP_HEIGHT) + 0.5) / MAP_HEIGHT * np.pi
    lon, lat = np.meshgrid(lon, lat)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1).reshape(-1, 3)


def _lonlat(xyz: np.ndarray) -> list[list[float]]:
    xyz = np.asarray(xyz).reshape(-1, 3)
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    lat = np.degrees(np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0)))
    return np.round(np.stack([lon, lat], axis=1), 2).tolist()


def _b64(array: np.ndarray) -> str:
    return base64.b64encode(np.ascontiguousarray(array).tobytes()).decode("ascii")


def _maps(layout) -> dict:
    """The layout's per-pixel maps, base64 for the page. `kind` is 0 ocean, 1 margin, 2 belt,
    3 core, with bit 2 set on land; craton and continent ids are shifted so 0 is none."""
    nearest = layout.nearest(_pixel_xyz())
    kind = np.zeros(layout.lattice.size, dtype=np.uint8)
    kind[layout.continent >= 0] = 1
    kind[layout.belt] = 2
    kind[layout.core >= 0] = 3
    kind |= layout.land.astype(np.uint8) << 2
    return {
        "width": MAP_WIDTH,
        "height": MAP_HEIGHT,
        "kind": _b64(kind[nearest]),
        "plate": _b64(layout.plate[nearest].astype(np.uint16)),
        "craton": _b64((layout.core[nearest] + 1).astype(np.uint16)),
        "continent": _b64((layout.continent[nearest] + 1).astype(np.uint16)),
        "nuclei": _lonlat(layout.nuclei),
        "nucleus_craton": layout.nucleus_craton.tolist(),
        "sites": _lonlat(layout.sites) if len(layout.sites) else [],
    }


def _stats(layout) -> dict:
    area = layout.lattice.area_km2
    diagnostics = layout.diagnostics()
    diagnostics.pop("continents", None)
    return {
        **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in diagnostics.items()},
        "continents": int(len(layout.continent_area_km2)),
        "plates": int(layout.plate.max()) + 1,
        "land_fraction": round(float(area[layout.land].sum() / area.sum()), 4),
        "J_start": round(layout.refine_costs[0]),
        "J_end": round(layout.refine_costs[-1]),
        "lattice_points": int(layout.lattice.size),
    }


def run_job(job: dict) -> dict:
    """One layout: `job` has `index`, `seed`, `world` (num_plates, continental_fraction,
    land_fraction, node_density) and `overrides` (knob name -> value)."""
    from app import lithosphere_plate

    started = time.time()
    world = job["world"]
    result = {"index": job["index"], "seed": job["seed"], "overrides": job["overrides"]}
    try:
        _patch(job["overrides"])
        layout = lithosphere_plate.cratons_first_layout(
            job["seed"],
            num_plates=world.get("num_plates"),
            continental_fraction=world.get("continental_fraction"),
            land_fraction=world.get("land_fraction"),
            node_density=float(world.get("node_density") or 1.0),
        ).layout
    except Exception as exc:  # any knob can break the layout; the page shows it on the tile
        result.update(ok=False, error=f"{type(exc).__name__}: {exc}", trace=traceback.format_exc(limit=4), seconds=time.time() - started)
        return result
    result.update(ok=True, seconds=time.time() - started, stats=_stats(layout), **_maps(layout))
    return result


class Handler(BaseHTTPRequestHandler):
    pool: ProcessPoolExecutor
    workers: int
    defaults: list[dict]
    # Guards `pool`, `active` and `stopped` across request threads.
    lock = threading.RLock()
    # Batch token -> its futures, so /api/stop can cancel a batch's queued jobs.
    active: dict[str, list] = {}
    # Tokens stopped before their batch arrived; the batch is refused when it does.
    stopped: set[str] = set()

    def log_message(self, format, *args):  # noqa: A002 -- BaseHTTPRequestHandler's signature
        pass

    def do_GET(self):
        if self.path not in ("/", "/index.html"):
            return self._send(404, "text/plain", b"not found")
        body = PAGE.read_text().replace("/*__KNOBS__*/[]", json.dumps(self.defaults)).encode()
        return self._send(200, "text/html; charset=utf-8", body)

    def do_POST(self):
        try:
            request = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        except (TypeError, ValueError) as exc:  # no Content-Length, or malformed JSON
            return self._send(400, "text/plain", f"bad request: {exc}".encode())
        if not isinstance(request, dict):
            return self._send(400, "text/plain", b"bad request: expected a JSON object")
        token = str(request.get("batch"))
        if self.path == "/api/stop":
            return self._stop(token)
        if self.path != "/api/generate":
            return self._send(404, "text/plain", b"not found")
        defaults = {d["name"]: d["default"] for d in self.defaults}
        try:
            jobs = []
            for i, job in enumerate(request["jobs"][:MAX_JOBS]):
                overrides = {k: _jsonable(_coerce(k, v)) for k, v in job.get("overrides", {}).items()}
                _check_order(overrides, defaults)
                jobs.append({"index": i, "seed": int(job["seed"]), "world": request.get("world", {}), "overrides": overrides})
        except (KeyError, ValueError, TypeError) as exc:
            return self._send(400, "text/plain", str(exc).encode())
        return self._stream(token, jobs)

    def _stop(self, token: str) -> None:
        """Cancel a batch's queued jobs; jobs a worker has already taken still finish. A stop
        that beats its batch here is remembered, and the batch is refused on arrival."""
        with self.lock:
            futures = self.active.pop(token, None)
            if futures is None:
                self.stopped.add(token)
        for future in futures or []:
            future.cancel()
        return self._send(200, "text/plain", b"stopped")

    def _submit(self, jobs: list[dict]) -> tuple[ProcessPoolExecutor, list]:
        with self.lock:
            pool = Handler.pool
            try:
                return pool, [pool.submit(run_job, job) for job in jobs]
            except BrokenProcessPool:
                self._replace_pool(pool)
                return Handler.pool, [Handler.pool.submit(run_job, job) for job in jobs]

    def _replace_pool(self, broken: ProcessPoolExecutor) -> None:
        """A worker died (killed, out of memory): every later submit to that pool fails, so
        start a new one, once, however many request threads notice."""
        with self.lock:
            if Handler.pool is broken:
                Handler.pool = ProcessPoolExecutor(max_workers=self.workers)
        broken.shutdown(wait=False, cancel_futures=True)

    def _stream(self, token: str, jobs: list[dict]) -> None:
        """Run `jobs` and write each result as one NDJSON line as it finishes. A failure the
        page can't infer from a closed stream is sent as a final `{"fatal": ...}` line."""
        self.send_response(200)
        self.send_header("Content-Type", "application/x-ndjson")
        self.end_headers()
        with self.lock:
            if token in self.stopped:
                self.stopped.discard(token)
                return
            pool, futures = self._submit(jobs)
            self.active[token] = futures
        try:
            for future in as_completed(futures):
                if not future.cancelled():
                    self._write(future.result())
        except (BrokenPipeError, ConnectionResetError, CancelledError):
            pass  # the page closed the stream or stopped the batch
        except Exception as exc:  # a dead worker, or a result that won't serialise
            if isinstance(exc, BrokenProcessPool):
                self._replace_pool(pool)
            try:
                self._write({"fatal": f"{type(exc).__name__}: {exc}"})
            except OSError:
                pass
        finally:
            for future in futures:
                future.cancel()
            with self.lock:
                self.active.pop(token, None)

    def _write(self, record: dict) -> None:
        self.wfile.write((json.dumps(record) + "\n").encode())
        self.wfile.flush()

    def _send(self, code: int, content_type: str, body: bytes) -> None:
        self.send_response(code)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--port", type=int, default=8765, help="never 8000 or 5173, the user's own instance")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    Handler.defaults = knob_defaults()
    Handler.workers = args.workers
    Handler.pool = ProcessPoolExecutor(max_workers=args.workers)
    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"Layout explorer on http://127.0.0.1:{args.port}/", flush=True)
    try:
        server.serve_forever()
    finally:
        Handler.pool.shutdown(cancel_futures=True)


if __name__ == "__main__":
    main()
