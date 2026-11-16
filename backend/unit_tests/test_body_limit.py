"""Request body cap (#920): refused before a route buffers it, declared length or not."""

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


def _chunks(total: int, size: int = 1024):
    for _ in range(total // size):
        yield b" " * size


def test_small_public_body_passes():
    client, _ = _client()
    response = client.post("/public/worlds/abc/step", json={"revision": "r1"})
    assert response.status_code == 200 and response.json() == {"revision": "r1"}


def test_declared_oversize_public_body_is_refused_with_headers():
    client, _ = _client()
    response = client.post("/public/worlds/abc/step", content=b" " * (MAX_PUBLIC_BODY_BYTES + 1))
    assert response.status_code == 413
    assert response.headers["x-content-type-options"] == "nosniff"


def test_chunked_oversize_body_is_cut_off_at_the_cap():
    client, seen = _client()
    response = client.post("/worlds/abc/generate", content=_chunks(MAX_BODY_BYTES + 64 * 1024, 64 * 1024))
    assert response.status_code == 413
    # The route never received past the cap (the test client hands the body over in one
    # message, so here it received nothing at all).
    assert max(seen, default=0) <= MAX_BODY_BYTES


def test_chunked_oversize_public_body_never_reaches_validation():
    client, _ = _client()
    response = client.post("/public/worlds/abc/step", content=_chunks(MAX_PUBLIC_BODY_BYTES * 4))
    # FastAPI reports a failed body read as 400; either way the body was not buffered.
    assert response.status_code in (400, 413)


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
