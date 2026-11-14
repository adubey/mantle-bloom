"""Own the hosted simulation process and stream requests over a private Unix socket."""

from __future__ import annotations

import asyncio
import logging
import os
import sys
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from http.cookiejar import CookieJar, DefaultCookiePolicy
from pathlib import Path
from typing import AsyncContextManager

import httpx
from anyio import CancelScope
from fastapi import Request, Response
from fastapi.responses import JSONResponse, StreamingResponse
from starlette.background import BackgroundTask

_logger = logging.getLogger(__name__)
_HOP_HEADERS = {
    b"connection", b"keep-alive", b"proxy-authenticate", b"proxy-authorization",
    b"te", b"trailer", b"transfer-encoding", b"upgrade",
}


class WorkerUnavailable(RuntimeError):
    pass


class HostedWorker:
    def __init__(self, socket: Path, idle_seconds: float) -> None:
        self.socket = socket
        self.idle_seconds = idle_seconds
        self.process: asyncio.subprocess.Process | None = None
        self.generation = 0
        self._ready = False
        self._startup = asyncio.Lock()
        self._requests = 0
        self._last_activity = 0.0
        self._client = httpx.AsyncClient(
            base_url="http://simulation",
            transport=httpx.AsyncHTTPTransport(uds=str(socket)),
            timeout=httpx.Timeout(None, connect=2),
            # Public responses set host-only cookies. Forward them to the
            # caller without putting them in the shared client cookie jar.
            cookies=CookieJar(policy=DefaultCookiePolicy(allowed_domains=())),
        )

    @property
    def running(self) -> bool:
        return self.process is not None and self.process.returncode is None

    @property
    def ready(self) -> bool:
        return self.running and self._ready

    def touch(self) -> None:
        self._last_activity = asyncio.get_running_loop().time()

    @asynccontextmanager
    async def lease(self, *, wake: bool, activity: bool) -> AsyncIterator[bool]:
        # Requests can arrive while an idle check awaits the worker. Count them
        # before any await so the check cannot retire their worker underneath them.
        self._requests += 1
        try:
            if wake:
                await self._ensure_ready()
            yield self.ready
        finally:
            if activity:
                self.touch()
            self._requests -= 1

    async def _ensure_ready(self) -> None:
        async with self._startup:
            if self.running and self._ready:
                return
            if not self.running:
                self._ready = False
                self.socket.unlink(missing_ok=True)
                self.process = await asyncio.create_subprocess_exec(
                    sys.executable, "-m", "uvicorn", "app.main:app",
                    "--uds", str(self.socket), "--no-proxy-headers", "--no-access-log",
                    "--timeout-graceful-shutdown", "5",
                    env={**os.environ, "MANTLE_BLOOM_HOSTED_WORKER": "1", "MANTLE_BLOOM_LOCAL_NO_LOGIN": "0"},
                )
                self.generation += 1
            try:
                async with asyncio.timeout(120):
                    while self.running:
                        try:
                            response = await self._client.get("/healthz", timeout=2)
                        except httpx.RequestError:
                            await asyncio.sleep(0.05)
                            continue
                        if response.is_success:
                            self._ready = True
                            return
                        break
            except TimeoutError:
                _logger.error("simulation worker startup timed out")
            await self._stop()
            raise WorkerUnavailable("the simulation could not start")

    async def _stop(self) -> None:
        self._ready = False
        if self.running:
            process = self.process
            try:
                process.terminate()
            except ProcessLookupError:
                pass
            try:
                await asyncio.wait_for(process.wait(), 7)
            except TimeoutError:
                try:
                    process.kill()
                except ProcessLookupError:
                    pass
                await process.wait()
        self.process = None

    async def close(self) -> None:
        async with self._startup:
            await self._stop()
        await self._client.aclose()

    async def retire_if_idle(self) -> None:
        if not self.running or not self._ready or self._requests or not self._idle_timeout_elapsed():
            return
        async with self._startup:
            try:
                response = await self._client.get("/_worker/idle", timeout=2)
                idle = response.is_success and response.json().get("idle")
            except (httpx.RequestError, ValueError):
                return
            if idle and not self._requests and self._idle_timeout_elapsed():
                await self._stop()

    def _idle_timeout_elapsed(self) -> bool:
        return asyncio.get_running_loop().time() - self._last_activity >= self.idle_seconds

    async def warm(self, world_id: str) -> None:
        was_ready = self.ready
        async with self.lease(wake=True, activity=False):
            response = await self._client.post(f"/_worker/warm/{world_id}", timeout=120)
            if response.status_code >= 500:
                raise WorkerUnavailable("the checkpoint could not be opened")
            if not was_ready or (response.is_success and response.json().get("warmed")):
                self.touch()

    async def prepare(self) -> None:
        async with self.lease(wake=True, activity=True):
            pass

    async def openapi(self) -> dict:
        async with self.lease(wake=True, activity=True):
            response = await self._client.get("/openapi.json", timeout=10)
            response.raise_for_status()
            return response.json()

    async def read(self, path: str, request: Request, params: dict | None = None) -> httpx.Response | None:
        async with self.lease(wake=False, activity=False) as ready:
            if not ready:
                return None
            try:
                return await self._client.get(
                    path,
                    headers=[(k, v) for k, v in request.headers.raw if k.lower() in (b"cookie", b"host")],
                    params=params,
                    timeout=2,
                )
            except httpx.RequestError:
                if not self.ready:
                    return None
                # While a live worker owns visibility changes, a timeout
                # cannot authorize serving older metadata from disk.
                return httpx.Response(503, json={"detail": "The simulation is unavailable. Try again shortly."})

    async def forward(self, request: Request, *, wake: bool = True, activity: bool = True) -> Response | None:
        lease = self.lease(wake=wake, activity=activity)
        entered = False
        try:
            ready = await lease.__aenter__()
            entered = True
            if not ready:
                await lease.__aexit__(None, None, None)
                return None
            upstream = await self._client.send(_upstream_request(request), stream=True)
        except (WorkerUnavailable, httpx.RequestError, OSError):
            if entered:
                await lease.__aexit__(None, None, None)
            if not wake and not self.ready:
                return None
            return JSONResponse({"detail": "The simulation is waking up. Try again shortly."}, 503)
        except BaseException:
            if entered:
                await lease.__aexit__(None, None, None)
            raise

        return self._stream_response(upstream, lease)

    @staticmethod
    def _stream_response(upstream: httpx.Response, lease: AsyncContextManager[bool]) -> Response:
        released = False

        async def release():
            nonlocal released
            if released:
                return
            released = True
            with CancelScope(shield=True):
                try:
                    await upstream.aclose()
                finally:
                    await lease.__aexit__(None, None, None)

        async def body():
            try:
                async for chunk in upstream.aiter_raw():
                    yield chunk
            finally:
                await release()

        response = StreamingResponse(
            body(), status_code=upstream.status_code, background=BackgroundTask(release),
        )
        for key, value in upstream.headers.raw:
            if key.lower() not in _HOP_HEADERS:
                response.headers.append(key.decode("latin-1"), value.decode("latin-1"))
        return response


def _upstream_request(request: Request) -> httpx.Request:
    # Preserve the original Host, cookies and Origin. The worker trusts neither a
    # loopback address nor a forwarding header as authentication.
    excluded = _HOP_HEADERS | {
        token.strip().lower().encode("latin-1") for token in request.headers.get("connection", "").split(",")
    }
    headers = [(k, v) for k, v in request.headers.raw if k.lower() not in excluded]
    url = httpx.URL("http://simulation").copy_with(
        raw_path=request.scope["raw_path"] + (
            b"?" + request.scope["query_string"] if request.scope["query_string"] else b""
        ),
    )
    # Uvicorn exposes the socket worker at the root, like the public deployment.
    content = request.stream() if request.method not in ("GET", "HEAD", "OPTIONS") else None
    return httpx.Request(request.method, url, headers=headers, content=content)
