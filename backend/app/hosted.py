"""Lightweight Railway entry point. Only its child process imports the simulator."""

from __future__ import annotations

import asyncio
import base64
import binascii
import logging
import os
import re
import tempfile
from contextlib import asynccontextmanager
from pathlib import Path

import httpx
from fastapi import FastAPI, HTTPException, Request, Response
from fastapi.openapi.docs import get_redoc_html, get_swagger_ui_html, get_swagger_ui_oauth2_redirect_html
from fastapi.openapi.utils import get_openapi
from fastapi.responses import JSONResponse, RedirectResponse
from starlette.routing import Route

from . import auth, checkpoint, public_view
from .frontend_files import mount_frontend
from .hosted_worker import HostedWorker, WorkerUnavailable

_logger = logging.getLogger(__name__)
_WORLD_PATH = "/public/worlds/{world_id}"
_WORLD_ID = re.compile(r"[0-9a-f]{32}")
_IDLE_SECONDS = 90.0
_RECOVERY_RETRY_SECONDS = 5.0
_STATE_POLL_INTERVAL = 0.5
_NO_STORE = {"Cache-Control": "no-store"}


def _root() -> checkpoint.CheckpointStore:
    store = checkpoint.CheckpointStore.from_env()
    if store is None:
        raise RuntimeError("hosted mode requires a checkpoint volume")
    return store


def _records(root: checkpoint.CheckpointStore) -> list[dict]:
    try:
        records = root.load_index() or []
    except (OSError, checkpoint.RegistryIndexError) as exc:
        raise HTTPException(503, "the world index is unavailable") from exc
    records.sort(key=lambda row: row["updated_at"], reverse=True)
    return records


def _idle_seconds() -> float:
    message = "MANTLE_BLOOM_IDLE_SECONDS must be positive and finite"
    try:
        seconds = float(os.environ.get("MANTLE_BLOOM_IDLE_SECONDS", _IDLE_SECONDS))
    except ValueError as exc:
        raise ValueError(message) from exc
    if not 0 < seconds < float("inf"):
        raise ValueError(message)
    return seconds


@asynccontextmanager
async def _lifespan(app: FastAPI):
    root = _root()
    root.check_writable()
    # Shared with every child in this container when OAuth signing is unconfigured.
    public_view.share_fallback_secret()
    idle_seconds = _idle_seconds()
    records = _records(root)
    stores = [root.for_world(row["world_id"]) for row in records[:checkpoint.max_worlds()]]
    for store in stores:
        store.check_writable()
    with tempfile.TemporaryDirectory(prefix="mb-worker-") as directory:
        worker = HostedWorker(Path(directory) / "app.sock", idle_seconds)
        app.state.worker = worker
        app.state.warm_task = None
        app.state.public_warm_after = 0.0
        app.state.openapi_lock = asyncio.Lock()
        app.state.openapi = None
        supervisor = asyncio.create_task(_supervise_worker(root, worker))
        try:
            yield
        finally:
            tasks = [supervisor]
            if app.state.warm_task is not None:
                tasks.append(app.state.warm_task)
            for task in tasks:
                task.cancel()
            await asyncio.gather(*tasks, return_exceptions=True)
            await worker.close()


def _saved_run_to_resume(root: checkpoint.CheckpointStore) -> bool:
    records = _records(root)
    stores = [root.for_world(row["world_id"]) for row in records[:checkpoint.max_worlds()]]
    if not records:
        stores = [root]
    return any(store.world_path.is_file() and store.load_status().running for store in stores)


async def _supervise_worker(root: checkpoint.CheckpointStore, worker: HostedWorker) -> None:
    retry_at = 0.0
    clock = asyncio.get_running_loop().time
    while True:
        if worker.ready:
            await worker.retire_if_idle()
        elif clock() >= retry_at:
            try:
                if _saved_run_to_resume(root):
                    await worker.prepare()
            except (OSError, WorkerUnavailable, httpx.RequestError, HTTPException):
                _logger.exception("could not resume saved Play sessions")
            retry_at = clock() + _RECOVERY_RETRY_SECONDS
        await asyncio.sleep(min(worker.idle_seconds, _RECOVERY_RETRY_SECONDS))


async def _warm(app: FastAPI, world_id: str) -> None:
    worker = app.state.worker
    try:
        await worker.warm(world_id)
    except (OSError, WorkerUnavailable, httpx.RequestError):
        _logger.exception("could not warm the simulation")
    finally:
        app.state.public_warm_after = asyncio.get_running_loop().time() + 2 * worker.idle_seconds
        app.state.warm_task = None


def _public_record(world_id: str, *, include_running: bool = False) -> tuple[dict, checkpoint.CheckpointStore]:
    root = _root()
    for record in _records(root)[:checkpoint.max_worlds()]:
        if record["world_id"] == world_id:
            store = root.for_world(world_id)
            if record["is_public"] or (include_running and store.load_status().running):
                return record, store
            break
    raise HTTPException(404, "no such public world")


def _summary(store: checkpoint.CheckpointStore) -> dict | None:
    saved = store.load_hosted_summary()
    return None if saved is None else _viewer_summary(saved)


def _viewer_summary(saved: dict) -> dict:
    return {key: saved[key] for key in ("revision", "seed", "elapsed_years")}


def _presentation(record: dict) -> tuple[str, list[float]]:
    view = record.get("viewer_view")
    return (
        view if view in public_view.GUEST_VIEWS else public_view.GUEST_DEFAULT_VIEW,
        record.get("viewer_rotation") or list(public_view.GUEST_IDENTITY_ROTATION),
    )


def _state(record: dict, world: dict, request: Request, guest_demo: bool = False, *, running: bool = False) -> dict:
    view, rotation = _presentation(record)
    return {
        "name": record["name"], "world": world, "running": running,
        "is_default": record.get("is_default", False), "view": view, "rotation": rotation,
        "steps_remaining": public_view.steps_remaining(request, auth.load_config(request), guest_demo),
        "step_limit": public_view.GUEST_STEP_LIMIT, "step_years": public_view.GUEST_STEP_YEARS,
        "step_progress": None, "snapshot_first": True,
    }


def _snapshot(
    store: checkpoint.CheckpointStore, record: dict, *, revision: str | None = None,
) -> tuple[dict, str] | None:
    generation = store.world_generation()
    world = _summary(store)
    if generation is None or world is None or (revision is not None and world["revision"] != revision):
        return None
    entry = store.load_viewer_image()
    if entry is None or entry.get("generation") != generation:
        return None
    view, rotation = _presentation(record)
    expected = {
        "seed": world["seed"], "elapsed_years": world["elapsed_years"], "view": view, "rotation": rotation,
        "width": public_view.GUEST_DEFAULT_WIDTH, "height": public_view.GUEST_DEFAULT_HEIGHT,
    }
    if entry["key"] != expected or store.world_generation() != generation:
        return None
    return world, entry["image_base64"]


async def _active(request: Request) -> Response | None:
    return await request.app.state.worker.forward(request, wake=False, activity=False)


def _passive(request: Request) -> bool:
    path = auth._route_path(request.scope)[1].rstrip("/")
    return request.method == "GET" and (path == "/worlds" or path.endswith("/run"))


async def _proxy(request: Request) -> Response:
    worker = request.app.state.worker
    return await worker.forward(request, activity=not _passive(request))


def _admin_records(root: checkpoint.CheckpointStore) -> dict:
    worlds = []
    cap = checkpoint.max_worlds()
    for index, row in enumerate(_records(root)):
        store = root.for_world(row["world_id"])
        load_error = None
        summary = None
        status = checkpoint.RunStatus()
        if index < cap:
            status = store.load_status()
            summary = store.load_hosted_summary()
            if summary is not None:
                summary["replay_logs"] = []
        else:
            load_error = checkpoint.world_limit_error(cap)
        worlds.append({
            **row, "load_error": load_error, "summary": summary,
            "run": {"running": status.running, "step_years": status.step_years, "stepping": False, "error": None},
        })
    return {"worlds": worlds}


async def _public_records(request: Request) -> Response:
    active = await _active(request)
    if active is not None:
        return active
    root = _root()
    records = _records(root)
    worlds = []
    for row in records[:checkpoint.max_worlds()]:
        store = root.for_world(row["world_id"])
        running = store.load_status().running
        if not (row["is_public"] or running):
            continue
        summary = _summary(store)
        if summary is None and store.world_path.is_file():
            # A checkpoint from an older version needs a one-time metadata upgrade.
            return await _proxy(request)
        if summary is not None:
            worlds.append({
                **summary, "world_id": row["world_id"], "name": row["name"], "running": running,
                "is_default": row.get("is_default", False),
            })
    if not records and root.world_path.is_file():
        return await _proxy(request)
    return JSONResponse({"worlds": worlds}, headers=_NO_STORE)


async def _world_state(
    world_id: str, request: Request, after: str | None = None, guest_demo: bool = False,
) -> Response:
    # Hold the long poll here, with only short worker reads. A read-only viewer
    # can stay open past the idle timeout without pinning the simulation.
    clock = asyncio.get_running_loop().time
    started = clock()
    deadline = started + public_view.GUEST_POLL_WAIT_S
    first = None
    previous = None
    state = None
    while True:
        worker = request.app.state.worker
        source = (worker.generation, worker.ready)
        # Only the worker writes these files. Reuse a sleeping response until
        # a worker starts or exits, even if it fails before becoming ready.
        if state is None or source[1] or source != previous:
            state = await _current_state(world_id, request, guest_demo)
        previous = source
        if state is None:
            return await _proxy(request)
        if first is None:
            first = state
        # A repeated poll with the same revision and progress should still
        # yield time to the step, even if the client immediately polls again.
        now = clock()
        if (state["world"]["revision"] != after or state != first or now >= deadline
                or (state["step_progress"] is not None and now - started >= _STATE_POLL_INTERVAL)):
            return JSONResponse(state, headers=_NO_STORE)
        if await request.is_disconnected():
            return Response(status_code=204)
        await asyncio.sleep(_STATE_POLL_INTERVAL)


async def _live_world_state(world_id: str, request: Request) -> httpx.Response | None:
    if _WORLD_ID.fullmatch(world_id) is None:
        raise HTTPException(404, "no such public world")
    return await request.app.state.worker.read(
        f"/public/worlds/{world_id}/world", request,
        {
            "guest_demo": request.query_params.get("guest_demo", "false"),
            "following": "after" in request.query_params,
        },
    )


async def _current_state(world_id: str, request: Request, guest_demo: bool) -> dict | None:
    active = await _live_world_state(world_id, request)
    if active is not None:
        if active.status_code == 404:
            raise HTTPException(404, "no such public world")
        if not active.is_success:
            raise HTTPException(503, "public world state is unavailable")
        return {**active.json(), "snapshot_first": True}
    record, store = _public_record(world_id, include_running=True)
    world = _summary(store)
    if world is None:
        if store.world_path.is_file():
            return None
        raise HTTPException(404, "no such public world")
    return _state(record, world, request, guest_demo, running=store.load_status().running)


async def _preview(world_id: str, request: Request, revision: str, view: str, rotation: str) -> Response:
    record, store = _public_record(world_id)
    saved_view, saved_rotation = _presentation(record)
    try:
        matching_rotation = [float(value) for value in rotation.split(",")] == saved_rotation
    except ValueError:
        matching_rotation = False
    if view != saved_view or not matching_rotation:
        raise HTTPException(404, "no such public preview")
    snapshot = _snapshot(store, record, revision=revision)
    if snapshot is None:
        raise HTTPException(404, "no such public preview")
    try:
        png = base64.b64decode(snapshot[1], validate=True)
    except binascii.Error as exc:
        raise HTTPException(404, "no such public preview") from exc
    if not png.startswith(b"\x89PNG\r\n\x1a\n"):
        raise HTTPException(404, "no such public preview")
    # Recheck visibility after reading the file. The same request must not
    # authorize a preview against an older index after an unpublish.
    active = await _live_world_state(world_id, request)
    if active is not None and not _matches_active_preview(active, revision, saved_view, saved_rotation):
        raise HTTPException(404, "no such public preview")
    current, current_store = _public_record(world_id)
    if _presentation(current) != (saved_view, saved_rotation) or _summary(current_store) != snapshot[0]:
        raise HTTPException(404, "no such public preview")
    return Response(png, media_type="image/png", headers=_NO_STORE)


def _matches_active_preview(response: httpx.Response, revision: str, view: str, rotation: list[float]) -> bool:
    if not response.is_success:
        return False
    state = response.json()
    return ((state.get("world") or {}).get("revision"), state["view"], state["rotation"]) == (revision, view, rotation)


async def _warm_world(world_id: str, request: Request) -> Response:
    _, store = _public_record(world_id)
    if not store.world_path.is_file():
        raise HTTPException(404, "no such public world")
    now = asyncio.get_running_loop().time()
    # Apply the cooldown across worlds so rotating IDs cannot keep the worker
    # resident. Admin and step requests do not go through this hint.
    if request.app.state.warm_task is None and now >= request.app.state.public_warm_after:
        request.app.state.warm_task = asyncio.create_task(_warm(request.app, world_id))
    return Response(status_code=202)


async def _step(world_id: str, request: Request, guest_demo: bool = False) -> Response:
    worker = request.app.state.worker
    if not worker.running:
        _, store = _public_record(world_id, include_running=True)
        if not store.world_path.is_file():
            raise HTTPException(404, "no such public world")
        # Exhausted and stale requests can be answered from committed metadata.
        # Bound the body this process inspects. Request.stream() can forward
        # the cached body to the worker if an actual step is needed.
        world = _summary(store)
        size = request.headers.get("content-length", "")
        if world is not None and size.isdigit() and len(size) <= 4 and int(size) <= 1024:
            body = await request.body()
            try:
                payload = public_view.PublicStepRequest.model_validate_json(body)
            except ValueError:
                return await _proxy(request)
            if worker.running:
                return await _proxy(request)
            remaining = public_view.steps_remaining(request, auth.load_config(request), guest_demo)
            if payload.revision != world["revision"]:
                return JSONResponse({"world": world, "computed": False, "steps_remaining": remaining})
            if remaining == 0:
                raise HTTPException(429, "step allowance used up")
    return await _proxy(request)


async def _render(
    world_id: str, request: Request, view: str = "elevation", width: int = 1100, height: int = 611,
) -> Response:
    if view not in public_view.GUEST_VIEWS:
        raise HTTPException(400, "unknown public view")
    limit = public_view.GUEST_MAX_RENDER_PX
    if not (1 <= width <= limit and 1 <= height <= limit):
        raise HTTPException(400, f"width/height must be in [1, {limit}]")
    if not request.app.state.worker.running:
        record, store = _public_record(world_id, include_running=True)
        if not store.world_path.is_file():
            raise HTTPException(404, "no such public world")
        if (width, height) == (public_view.GUEST_DEFAULT_WIDTH, public_view.GUEST_DEFAULT_HEIGHT):
            snapshot = _snapshot(store, record)
            if snapshot is not None:
                world, image = snapshot
                return JSONResponse({**world, "image_base64": image}, headers=_NO_STORE)
    return await _proxy(request)


async def _admin_list(request: Request) -> Response:
    active = await _active(request)
    if active is not None:
        return active
    root = _root()
    payload = _admin_records(root)
    if any(
        row["summary"] is None and row["load_error"] is None and root.for_world(row["world_id"]).world_path.is_file()
        for row in payload["worlds"]
    ) or (not payload["worlds"] and root.world_path.is_file()):
        return await _proxy(request)
    return JSONResponse(payload, headers=_NO_STORE)


async def _admin_run(request: Request) -> Response:
    world_id = request.path_params["world_id"]
    if _WORLD_ID.fullmatch(world_id) is None:
        raise HTTPException(404, "unknown world")
    # An active run already vetoes sleep. Preserve its existing event-driven
    # poll and use a disk response once the worker is gone.
    active = await _active(request)
    if active is not None:
        return active
    world = next((row for row in _admin_records(_root())["worlds"] if row["world_id"] == world_id), None)
    if world is None:
        raise HTTPException(404, "unknown world")
    status = world["run"]
    saved = world["summary"]
    if saved is None and world["load_error"] is None and _root().for_world(world_id).world_path.is_file():
        return await _proxy(request)
    if status["running"]:
        return await _proxy(request)
    return JSONResponse({
        **status, "run_version": 0, "phase": None,
        "world": None if saved is None else _viewer_summary(saved),
    }, headers=_NO_STORE)


def _health() -> dict:
    root = _root()
    try:
        root.check_writable()
    except OSError as exc:
        raise HTTPException(503, "the checkpoint volume is not writable") from exc
    return {"status": "ok"}


def create_app(dist: str | None = None) -> FastAPI:
    app = FastAPI(title="Mantle-Bloom", lifespan=_lifespan, openapi_url=None)
    app.add_middleware(auth.AuthMiddleware, allow_local=False)
    app.include_router(auth.router)
    app.add_api_route("/healthz", _health)
    app.add_api_route("/openapi.json", _openapi, include_in_schema=False)
    app.add_api_route("/docs", _swagger, include_in_schema=False)
    app.add_api_route("/docs/oauth2-redirect", get_swagger_ui_oauth2_redirect_html, include_in_schema=False)
    app.add_api_route("/redoc", _redoc, include_in_schema=False)
    app.add_api_route("/demo", _demo)
    app.add_api_route("/public/worlds", _public_records)
    app.add_api_route(_WORLD_PATH + "/world", _world_state)
    app.add_api_route(_WORLD_PATH + "/render", _render)
    app.add_api_route(_WORLD_PATH + "/preview.png", _preview)
    app.add_api_route(_WORLD_PATH + "/warm", _warm_world, methods=["POST"])
    app.add_api_route(_WORLD_PATH + "/step", _step, methods=["POST"])
    app.add_api_route("/worlds", _admin_list)
    app.add_api_route("/_worker/{path:path}", _private_route, include_in_schema=False)
    app.routes.append(Route("/worlds/{world_id}/run", _admin_run, methods=["GET"]))
    for path in ("/worlds", "/worlds/{path:path}", "/debug/{path:path}"):
        app.routes.append(Route(path, _proxy, methods=["GET", "POST", "PATCH", "DELETE", "PUT", "OPTIONS"]))
    if dist:
        mount_frontend(app, dist)
    return app


async def _openapi(request: Request) -> dict:
    app = request.app
    async with app.state.openapi_lock:
        if app.state.openapi is None:
            try:
                worker_schema = await app.state.worker.openapi()
            except (OSError, WorkerUnavailable, httpx.HTTPError) as exc:
                raise HTTPException(503, "the API schema is temporarily unavailable") from exc
            # The simulator owns admin operations and their request/response
            # models. Add only public operations implemented by this process.
            public_schema = get_openapi(title=app.title, version=app.version, routes=app.routes)
            for path, operations in public_schema["paths"].items():
                worker_schema["paths"].setdefault(path, operations)
            schemas = worker_schema.setdefault("components", {}).setdefault("schemas", {})
            for name, schema in public_schema.get("components", {}).get("schemas", {}).items():
                schemas.setdefault(name, schema)
            app.state.openapi = worker_schema
    return app.state.openapi


def _swagger(request: Request) -> Response:
    root = request.scope.get("root_path", "")
    return get_swagger_ui_html(
        openapi_url=f"{root}/openapi.json", title="Mantle-Bloom API",
        oauth2_redirect_url=f"{root}/docs/oauth2-redirect",
    )


def _redoc(request: Request) -> Response:
    return get_redoc_html(openapi_url=f'{request.scope.get("root_path", "")}/openapi.json', title="Mantle-Bloom API")


def _demo() -> RedirectResponse:
    return RedirectResponse("/login?demo=true", status_code=302)


def _private_route(path: str) -> None:
    raise HTTPException(404, "not found")


app = create_app(os.environ.get("MANTLE_BLOOM_FRONTEND_DIST"))
