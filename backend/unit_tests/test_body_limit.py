"""Request body cap (#920): refused before a route buffers it, declared length or not."""

import asyncio

import pytest
from fastapi import FastAPI, Request
from fastapi.testclient import TestClient
from pydantic import BaseModel

from app import body_limit, hosted, main, security_headers
from app.body_limit import MAX_BODY_BYTES, MAX_PUBLIC_BODY_BYTES


class _Step(BaseModel):
    revision: str


def _client(root_path: str = "") -> tuple[TestClient, list[int]]:
    app = FastAPI(root_path=root_path)
    body_limit.install(app)
    security_headers.install(app)
    seen: list[int] = []

    @app.post("/public/worlds/{world_id}/step")
    def _step(world_id: str, req: _Step) -> dict:
        return {"revision": req.revision}

    @app.post("/worlds/import")
    @app.post("/worlds/{world_id}/load")
    @app.post("/worlds/{world_id}/generate")
    async def _count(request: Request) -> dict:
        size = 0
        async for chunk in request.stream():
            size += len(chunk)
            seen.append(size)
        return {"size": size}

    return TestClient(app), seen


def test_small_public_body_passes():
    client, _ = _client()
    response = client.post("/public/worlds/abc/step", json={"revision": "r1"})
    assert response.status_code == 200 and response.json() == {"revision": "r1"}


def test_declared_oversize_public_body_is_refused_with_headers():
    client, _ = _client()
    response = client.post("/public/worlds/abc/step", content=b" " * (MAX_PUBLIC_BODY_BYTES + 1))
    assert response.status_code == 413
    assert response.headers["x-content-type-options"] == "nosniff"


async def _raw_post(app, path: str, messages: list[bytes], headers: list[tuple[bytes, bytes]] = ()) -> tuple[int, dict, int]:
    """POST `messages` as separate ASGI body messages with no Content-Length, the way a
    chunked upload arrives (the test client hands a body over as one message). Returns the
    status, the response headers and how many body messages the app asked for."""
    pending = list(messages)
    asked = 0
    sent: list[dict] = []

    async def receive():
        nonlocal asked
        if not pending:
            return {"type": "http.disconnect"}
        asked += 1
        body = pending.pop(0)
        return {"type": "http.request", "body": body, "more_body": bool(pending)}

    async def send(message):
        sent.append(message)

    scope = {
        "type": "http", "asgi": {"version": "3.0"}, "http_version": "1.1", "method": "POST", "path": path,
        "raw_path": path.encode(), "root_path": "", "query_string": b"", "scheme": "http",
        "headers": [(b"host", b"testserver"), (b"content-type", b"application/json"), *headers],
        "client": ("203.0.113.5", 1234), "server": ("testserver", 80),
    }
    await app(scope, receive, send)
    start = next(m for m in sent if m["type"] == "http.response.start")
    return start["status"], {k.decode(): v.decode() for k, v in start["headers"]}, asked


def test_chunked_oversize_body_stops_at_the_cap_with_413():
    client, seen = _client()
    chunk = b" " * (1024 * 1024)
    total = MAX_BODY_BYTES // len(chunk) + 8
    status, _, asked = asyncio.run(_raw_post(client.app, "/worlds/abc/generate", [chunk] * total))
    assert status == 413
    assert max(seen) <= MAX_BODY_BYTES
    assert asked < total


def test_chunked_oversize_json_body_is_413_not_a_parse_error():
    client, _ = _client()
    chunk = b" " * 1024
    status, headers, asked = asyncio.run(_raw_post(client.app, "/public/worlds/abc/step", [chunk] * 10))
    assert status == 413
    assert headers["x-content-type-options"] == "nosniff"
    assert asked <= MAX_PUBLIC_BODY_BYTES // len(chunk) + 1


def test_malformed_content_length_falls_back_to_counting():
    client, _ = _client()
    chunk = b" " * 1024
    for declared in (b"-5", b"abc"):
        status, _, _ = asyncio.run(_raw_post(
            client.app, "/public/worlds/abc/step", [chunk] * 10, [(b"content-length", declared)],
        ))
        assert status == 413


def test_simulator_413_carries_cors_headers_for_the_dev_frontend():
    origin = next(iter(main._frontend_origins))
    response = TestClient(main.app).post(
        "/public/worlds/abc/step", content=b" " * (MAX_PUBLIC_BODY_BYTES + 1), headers={"Origin": origin},
    )
    assert response.status_code == 413
    assert response.headers["access-control-allow-origin"] == origin
    assert response.headers["x-content-type-options"] == "nosniff"


def test_hosted_refuses_a_declared_oversize_body_before_the_auth_guard():
    with pytest.MonkeyPatch.context() as patch:
        patch.setenv("MANTLE_BLOOM_LOCAL_NO_LOGIN", "0")
        app = hosted.create_app()
        response = TestClient(app).post("/worlds/abc/generate", content=b" " * (MAX_BODY_BYTES + 1))
    assert response.status_code == 413
    assert response.headers["x-frame-options"] == "DENY"


@pytest.mark.parametrize("path", ["/worlds/import", "/worlds/abc/load"])
def test_upload_routes_keep_their_own_larger_cap(path):
    client, _ = _client()
    response = client.post(path, content=b" " * (MAX_BODY_BYTES + 1))
    assert response.status_code == 200 and response.json() == {"size": MAX_BODY_BYTES + 1}


def test_root_path_is_stripped_before_matching():
    client, _ = _client(root_path="/mb")
    response = client.post("/mb/public/worlds/abc/step", content=b" " * (MAX_PUBLIC_BODY_BYTES + 1))
    assert response.status_code == 413


def test_limit_for():
    assert body_limit.limit_for("POST", "/worlds/import") is None
    assert body_limit.limit_for("POST", "/worlds/abc/load") is None
    assert body_limit.limit_for("GET", "/worlds/abc/load") == MAX_BODY_BYTES
    assert body_limit.limit_for("POST", "/worlds/abc/load/x") == MAX_BODY_BYTES
    assert body_limit.limit_for("POST", "/public/worlds/abc/step") == MAX_PUBLIC_BODY_BYTES
    assert body_limit.limit_for("POST", "/worlds/abc/generate") == MAX_BODY_BYTES


def test_both_entry_points_install_the_cap():
    for app in (main.app, hosted.create_app()):
        assert any(m.cls is body_limit.BodyLimitMiddleware for m in app.user_middleware)
