"""A ceiling on request body size for every route (#920).

FastAPI reads a JSON body into memory before it validates it, and nothing upstream
bounded how much it would read: a logged-out visitor could post an arbitrarily large body
to `/public/worlds/{world_id}/step` and the worker buffered all of it. The cap is applied
here, ahead of routing, so it also covers a route added later.

- `/public/*` takes `MAX_PUBLIC_BODY_BYTES`: its only body is a step's `{"revision": ...}`.
- The two save-upload routes are exempt. They stream to disk under their own, much larger
  `MAX_LOAD_WORLD_BYTES` cap (`main._spool_upload_to_tempfile`).
- Everything else takes `MAX_BODY_BYTES`, sized for the largest legitimate JSON body: a
  generate request carrying a base64 sketch PNG.

A declared `Content-Length` over the cap is refused with `413` before any of the body is
read. A body with no declared length (chunked) is counted as it streams: the route receives
up to the cap, then the read fails and the answer is `413`.
"""

from __future__ import annotations

import re

from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse

from .auth import _route_path

MAX_PUBLIC_BODY_BYTES = 4 * 1024
MAX_BODY_BYTES = 32 * 1024 * 1024

_PUBLIC_PREFIX = "/public/"
# POST /worlds/import and POST /worlds/{world_id}/load.
_UPLOAD_PATH = re.compile(r"/worlds/(import|[^/]+/load)")


class _BodyTooLarge(HTTPException):
    """Raised from `receive` when a streamed body passes its cap. An `HTTPException` so
    FastAPI's body reader, which turns any other error into a `400`, lets it through as
    the `413` it is."""

    def __init__(self, limit: int) -> None:
        super().__init__(status_code=413, detail=_detail(limit))


def _detail(limit: int) -> str:
    return f"request body exceeds the {limit}-byte limit"


def limit_for(method: str, path: str) -> int | None:
    """The body cap for a request, or None for a route that enforces its own."""
    if method == "POST" and _UPLOAD_PATH.fullmatch(path):
        return None
    return MAX_PUBLIC_BODY_BYTES if path.startswith(_PUBLIC_PREFIX) else MAX_BODY_BYTES


def _declared_length(scope) -> int | None:
    for name, value in scope.get("headers") or ():
        if name == b"content-length":
            try:
                return int(value)
            except ValueError:
                return None
    return None


class BodyLimitMiddleware:
    """Pure ASGI, like `auth.AuthMiddleware`, so streamed responses pass through untouched."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        limit = limit_for(scope["method"], _route_path(scope)[1])
        if limit is None:
            await self.app(scope, receive, send)
            return
        refusal = JSONResponse({"detail": _detail(limit)}, status_code=413)
        declared = _declared_length(scope)
        if declared is not None and declared > limit:
            await refusal(scope, receive, send)
            return

        received = 0
        started = False

        async def counted_receive():
            nonlocal received
            message = await receive()
            if message["type"] == "http.request":
                received += len(message.get("body", b""))
                if received > limit:
                    raise _BodyTooLarge(limit)
            return message

        async def tracked_send(message):
            nonlocal started
            if message["type"] == "http.response.start":
                started = True
            await send(message)

        try:
            await self.app(scope, counted_receive, tracked_send)
        except Exception as exc:
            # Usually the app's own exception handling has already answered 413 and nothing
            # reaches here. A route that had begun answering can't be given a new status.
            if started or not _is_overrun(exc):
                raise
            await refusal(scope, receive, send)


def _is_overrun(exc: BaseException) -> bool:
    """Whether `exc` is the cap's own exception, bare or grouped: a streaming response
    reads `receive` inside a task group, which wraps what it raises in an ExceptionGroup."""
    if isinstance(exc, _BodyTooLarge):
        return True
    grouped = getattr(exc, "exceptions", None)
    return bool(grouped) and all(_is_overrun(inner) for inner in grouped)


def install(app: FastAPI) -> None:
    app.add_middleware(BodyLimitMiddleware)
