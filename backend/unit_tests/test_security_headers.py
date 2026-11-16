"""Browser security headers (#555 F5): present on every response, refusals included."""

from fastapi import FastAPI
from fastapi.responses import PlainTextResponse
from fastapi.testclient import TestClient

from app import hosted, security_headers
from app.security_headers import CONTENT_SECURITY_POLICY


def _client(root_path: str = "") -> TestClient:
    app = FastAPI(root_path=root_path)
    security_headers.install(app)

    @app.get("/boom")
    def _boom() -> PlainTextResponse:
        raise RuntimeError("unhandled")

    @app.get("/plain")
    def _plain() -> PlainTextResponse:
        return PlainTextResponse("ok")

    @app.get("/own-policy")
    def _own_policy() -> PlainTextResponse:
        return PlainTextResponse("ok", headers={"Content-Security-Policy": "default-src 'none'"})

    return TestClient(app)


def test_every_response_carries_the_headers():
    for path in ("/plain", "/missing"):
        headers = _client().get(path).headers
        assert headers["x-content-type-options"] == "nosniff"
        assert headers["x-frame-options"] == "DENY"
        assert headers["referrer-policy"] == "same-origin"
        assert headers["content-security-policy"] == CONTENT_SECURITY_POLICY


def test_policy_blocks_inline_and_foreign_scripts_and_framing():
    directives = dict(part.split(" ", 1) for part in CONTENT_SECURITY_POLICY.split("; "))
    assert directives["script-src"] == "'self'"
    assert directives["default-src"] == "'self'"
    assert directives["frame-ancestors"] == "'none'"
    assert directives["object-src"] == "'none'"


def test_hsts_only_for_https_requests():
    client = _client()
    assert "strict-transport-security" not in client.get("/plain").headers
    proxied = client.get("/plain", headers={"X-Forwarded-Proto": "https"}).headers
    assert proxied["strict-transport-security"] == security_headers.STRICT_TRANSPORT_SECURITY
    direct = TestClient(client.app, base_url="https://testserver").get("/plain").headers
    assert direct["strict-transport-security"] == security_headers.STRICT_TRANSPORT_SECURITY


def test_a_route_set_header_is_kept_not_duplicated():
    response = _client().get("/own-policy")
    assert response.headers.get_list("content-security-policy") == ["default-src 'none'"]
    assert response.headers["x-content-type-options"] == "nosniff"


def test_api_docs_keep_only_the_framing_rule():
    for client, prefix in ((_client(), ""), (_client("/prefix"), "/prefix")):
        for path in ("/docs", "/redoc"):
            response = client.get(prefix + path)
            assert response.status_code == 200
            assert response.headers["content-security-policy"] == "frame-ancestors 'none'"
            assert response.headers["x-frame-options"] == "DENY"


def test_hosted_auth_refusals_carry_the_headers(monkeypatch):
    for name in ("CLIENT_ID", "CLIENT_SECRET", "OWNER_ID", "SESSION_SECRET", "PUBLIC_URL"):
        monkeypatch.delenv(f"MANTLE_BLOOM_GITHUB_{name}", raising=False)
        monkeypatch.delenv(f"MANTLE_BLOOM_{name}", raising=False)
    # No lifespan: the guard refuses before any route needs the checkpoint volume.
    client = TestClient(hosted.create_app(), follow_redirects=False)
    refused = client.get("/worlds")
    assert refused.status_code == 401
    assert refused.headers["content-security-policy"] == CONTENT_SECURITY_POLICY
    assert refused.headers["x-frame-options"] == "DENY"
    login = client.get("/login")
    assert login.status_code == 200
    assert login.headers["content-security-policy"] == CONTENT_SECURITY_POLICY


def test_unhandled_exception_500_carries_the_headers():
    # Starlette builds this response outside every user middleware (PR #919 review).
    response = TestClient(_client().app, raise_server_exceptions=False).get("/boom")
    assert response.status_code == 500
    assert response.headers["content-security-policy"] == CONTENT_SECURITY_POLICY
    assert response.headers["x-frame-options"] == "DENY"
    assert response.headers["x-content-type-options"] == "nosniff"
