import re
import time
import urllib.parse
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from app import auth, main
from app.main import app

OWNER_ID = 4242
SECRET = "s" * 40


@pytest.fixture(autouse=True)
def configured(monkeypatch):
    """A fully configured GitHub login; tests that need it missing unset it themselves."""
    monkeypatch.setenv("MANTLE_BLOOM_GITHUB_CLIENT_ID", "client-id")
    monkeypatch.setenv("MANTLE_BLOOM_GITHUB_CLIENT_SECRET", "client-secret")
    monkeypatch.setenv("MANTLE_BLOOM_GITHUB_OWNER_ID", str(OWNER_ID))
    monkeypatch.setenv("MANTLE_BLOOM_SESSION_SECRET", SECRET)
    monkeypatch.delenv("MANTLE_BLOOM_PUBLIC_URL", raising=False)
    monkeypatch.delenv("MANTLE_BLOOM_LOCAL_NO_LOGIN", raising=False)


@pytest.fixture(autouse=True)
def _fresh_app_state():
    """One empty world record (ID `test-world`) and no pending animation stop,
    before and after: the CSRF tests post its animate/stop, which would
    otherwise leak a stop signal into later tests."""
    from app.world_registry import WorldSession

    for session in main._registry.list():
        if session.runner is not None:
            session.runner.close()
        main._registry.remove(session.world_id)
    session = main._wire_session(main._registry.adopt(WorldSession(world_id="test-world")))
    session.animation_stop_event.clear()
    yield
    session.world = None
    session.animation_stop_event.clear()


@pytest.fixture
def remote():
    """A remote browser: not loopback, and https so the Secure cookies round-trip."""
    return TestClient(app, base_url="https://testserver", follow_redirects=False)


def _local(**kwargs) -> TestClient:
    return TestClient(app, base_url="http://localhost", client=("127.0.0.1", 50000), follow_redirects=False, **kwargs)


def _fake_github(monkeypatch, user: dict, token: dict | None = None) -> list:
    calls = []

    def fake(req):
        calls.append(req)
        if req.full_url == auth.GITHUB_TOKEN_URL:
            return {"access_token": "gho_x"} if token is None else token
        assert req.full_url == auth.GITHUB_USER_URL
        assert req.get_header("Authorization") == "Bearer gho_x"
        return user

    monkeypatch.setattr(auth, "_http_json", fake)
    return calls


def _start(client: TestClient) -> str:
    resp = client.get("/auth/github")
    assert resp.status_code == 303
    return urllib.parse.parse_qs(urllib.parse.urlsplit(resp.headers["location"]).query)["state"][0]


def _sign_in(client: TestClient, monkeypatch, user_id: int = OWNER_ID):
    _fake_github(monkeypatch, {"id": user_id, "login": "owner"})
    state = _start(client)
    return client.get("/auth/github/callback", params={"code": "abc", "state": state})


def test_token_round_trips_only_with_same_secret_and_purpose():
    token = auth.sign_token({"uid": 1}, b"k" * 32, "session", 60)
    assert auth.read_token(token, b"k" * 32, "session")["uid"] == 1
    assert auth.read_token(token, b"j" * 32, "session") is None
    assert auth.read_token(token, b"k" * 32, "oauth-state") is None
    body, mac = token.split(".")
    assert auth.read_token(f"{body}x.{mac}", b"k" * 32, "session") is None
    assert auth.read_token("garbage", b"k" * 32, "session") is None
    assert auth.read_token(f"{body}.\u00e9", b"k" * 32, "session") is None


def test_expired_token_is_rejected(monkeypatch):
    token = auth.sign_token({"uid": 1}, b"k" * 32, "session", 60)
    monkeypatch.setattr(time, "time", lambda: 10**12)
    assert auth.read_token(token, b"k" * 32, "session") is None


def test_logged_out_root_redirects_to_view(remote):
    resp = remote.get("/")
    assert resp.status_code == 302
    assert resp.headers["location"] == "/view"


@pytest.mark.parametrize("method", ["HEAD", "PUT", "DELETE", "PATCH"])
def test_world_routes_reject_every_method(remote, method):
    assert remote.request(method, "/worlds/test-world/summary").status_code == 401


def test_prefixed_app_still_guards_world_routes():
    """Under a root_path the router strips the prefix, so the guard must too (PR #377 review)."""
    client = TestClient(app, base_url="https://testserver", root_path="/mb", follow_redirects=False)
    assert client.get("/mb/worlds/test-world/summary").status_code == 401
    assert client.post("/mb/worlds/test-world/animate/stop", headers={"Origin": "https://testserver"}).status_code == 401
    assert client.get("/mb/").headers["location"] == "/mb/view"


def test_world_routes_reject_unauthenticated_requests(remote):
    assert remote.get("/worlds/test-world/summary").status_code == 401
    assert remote.post("/worlds/test-world/step", json={"years": 1}).status_code == 401
    assert remote.post("/worlds/test-world/generate", json={}).status_code == 401
    assert remote.post("/worlds/test-world/load", content=b"x").status_code == 401
    # /debug/scenarios moved out from under /worlds/ but stays owner-only (#563 review).
    assert remote.get("/debug/scenarios").status_code == 401


def test_view_and_login_stay_public(remote):
    assert remote.get("/login").status_code == 200
    assert remote.get("/auth/session").json() == {
        "admin": False, "local": False, "login": None, "login_available": True,
    }


def test_railway_health_check_is_public(remote):
    # Railway's deploy health check (.railway/railway.ts) carries no session (#376).
    config = (Path(__file__).resolve().parents[2] / ".railway" / "railway.ts").read_text()
    match = re.search(r'healthcheckPath:\s*"([^"]+)"', config)
    assert match, "no healthcheckPath in .railway/railway.ts"
    assert remote.get(match.group(1)).json() == {"status": "ok"}


def test_unconfigured_server_fails_closed(monkeypatch, remote):
    for name in ("CLIENT_ID", "CLIENT_SECRET", "OWNER_ID"):
        monkeypatch.delenv(f"MANTLE_BLOOM_GITHUB_{name}", raising=False)
    assert remote.get("/worlds/test-world/summary").status_code == 401
    page = remote.get("/login")
    assert "isn&#x27;t configured" in page.text and 'href="/auth/github"' not in page.text
    assert remote.get("/auth/github").headers["location"] == "/login?error=unconfigured"


def test_unconfigured_login_page_tells_the_operator_what_to_set(monkeypatch, remote):
    # The admin's only symptom on a fresh Railway deploy is a /view whose Sign in
    # entry lands here, so the page names the setup instead of dead-ending (#626).
    for name in ("CLIENT_ID", "CLIENT_SECRET", "OWNER_ID"):
        monkeypatch.delenv(f"MANTLE_BLOOM_GITHUB_{name}", raising=False)
    page = remote.get("/login")
    assert "isn&#x27;t configured" in page.text
    for name in (
        "MANTLE_BLOOM_GITHUB_CLIENT_ID",
        "MANTLE_BLOOM_GITHUB_CLIENT_SECRET",
        "MANTLE_BLOOM_GITHUB_OWNER_ID",
        "MANTLE_BLOOM_SESSION_SECRET",
        "MANTLE_BLOOM_PUBLIC_URL",
    ):
        assert name in page.text
    assert "/auth/github/callback" in page.text
    assert 'href="/auth/github"' not in page.text


def test_short_session_secret_leaves_login_unconfigured(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_SESSION_SECRET", "short")
    assert not auth.load_config().login_configured


def test_local_requests_skip_login():
    with _local() as client:
        assert client.get("/worlds/test-world/summary").status_code == 404  # no world yet, but let through
        assert client.get("/auth/session").json()["admin"] is True


@pytest.mark.parametrize("headers", [
    {"X-Forwarded-For": "203.0.113.9"}, {"Forwarded": "for=203.0.113.9"}, {"Host": "mantle.example.com"},
    # A proxy that rewrites Host to localhost may forward only these (PR #377 review).
    {"X-Forwarded-Proto": "https"}, {"X-Forwarded-Port": "443"}, {"Via": "1.1 caddy"},
    {"CF-Connecting-IP": "203.0.113.9"},
])
def test_proxied_loopback_requests_need_login(headers):
    with _local() as client:
        assert client.get("/worlds/test-world/summary", headers=headers).status_code == 401


def test_local_bypass_can_be_disabled(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_LOCAL_NO_LOGIN", "0")
    with _local() as client:
        assert client.get("/worlds/test-world/summary").status_code == 401


def test_login_start_redirects_to_github_with_state_and_pkce(remote):
    resp = remote.get("/auth/github")
    url = urllib.parse.urlsplit(resp.headers["location"])
    query = urllib.parse.parse_qs(url.query)
    assert f"{url.scheme}://{url.netloc}{url.path}" == auth.GITHUB_AUTHORIZE_URL
    assert query["client_id"] == ["client-id"]
    assert query["redirect_uri"] == ["https://testserver/auth/github/callback"]
    assert query["code_challenge_method"] == ["S256"]
    cookie = resp.cookies[auth.STATE_COOKIE]
    assert auth.read_token(cookie, SECRET.encode(), "oauth-state")["state"] == query["state"][0]
    set_cookie = resp.headers["set-cookie"].lower()
    assert "httponly" in set_cookie and "secure" in set_cookie and "samesite=lax" in set_cookie


def test_public_url_sets_redirect_uri(monkeypatch, remote):
    monkeypatch.setenv("MANTLE_BLOOM_PUBLIC_URL", "https://mantle.example.com/")
    query = urllib.parse.parse_qs(urllib.parse.urlsplit(remote.get("/auth/github").headers["location"]).query)
    assert query["redirect_uri"] == ["https://mantle.example.com/auth/github/callback"]


def test_owner_signs_in_and_reaches_the_full_app(monkeypatch, remote):
    calls = _fake_github(monkeypatch, {"id": OWNER_ID, "login": "owner"})
    state = _start(remote)
    resp = remote.get("/auth/github/callback", params={"code": "abc", "state": state})
    assert resp.status_code == 303 and resp.headers["location"] == "/"
    token_form = urllib.parse.parse_qs(calls[0].data.decode())
    assert token_form["code"] == ["abc"] and token_form["code_verifier"][0]
    assert remote.get("/worlds/test-world/summary").status_code == 404  # authorized; just no world yet
    assert remote.get("/").status_code != 302
    assert remote.get("/auth/session").json() == {
        "admin": True, "local": False, "login": "owner", "login_available": True,
    }
    assert remote.get("/login").headers["location"] == "/"


def test_other_github_account_gets_no_session(monkeypatch, remote):
    resp = _sign_in(remote, monkeypatch, user_id=OWNER_ID + 1)
    assert resp.headers["location"] == "/login?error=forbidden"
    assert auth.SESSION_COOKIE not in remote.cookies
    assert remote.get("/worlds/test-world/summary").status_code == 401


@pytest.mark.parametrize("params, error", [
    ({"error": "access_denied", "state": "x"}, "denied"),
    ({"code": "abc", "state": "wrong"}, "state"),
    ({"code": "abc"}, "state"),
    ({"code": "abc", "state": "\u00e9"}, "state"),
])
def test_callback_failures_redirect_to_login(monkeypatch, remote, params, error):
    _fake_github(monkeypatch, {"id": OWNER_ID})
    _start(remote)
    resp = remote.get("/auth/github/callback", params=params)
    assert resp.headers["location"] == f"/login?error={error}"
    assert auth.SESSION_COOKIE not in remote.cookies


def test_callback_without_state_cookie_fails(monkeypatch, remote):
    _fake_github(monkeypatch, {"id": OWNER_ID})
    resp = remote.get("/auth/github/callback", params={"code": "abc", "state": "anything"})
    assert resp.headers["location"] == "/login?error=state"


def test_callback_reports_github_failures(monkeypatch, remote):
    _fake_github(monkeypatch, {"id": OWNER_ID}, token={"error": "bad_verification_code"})
    state = _start(remote)
    resp = remote.get("/auth/github/callback", params={"code": "abc", "state": state})
    assert resp.headers["location"] == "/login?error=github"


def test_login_page_shows_only_known_error_messages(remote):
    assert auth.LOGIN_ERRORS["forbidden"].replace("'", "&#x27;") in remote.get("/login?error=forbidden").text
    assert "<script>" not in remote.get("/login", params={"error": "<script>"}).text


def test_forged_or_stale_sessions_are_rejected(monkeypatch, remote):
    remote.cookies.set(auth.SESSION_COOKIE, auth.sign_token({"uid": OWNER_ID}, b"x" * 40, "session", 60))
    assert remote.get("/worlds/test-world/summary").status_code == 401
    # A session minted for a previous owner ID stops working when the owner changes.
    remote.cookies.set(auth.SESSION_COOKIE, auth.sign_token({"uid": OWNER_ID}, SECRET.encode(), "session", 60))
    monkeypatch.setenv("MANTLE_BLOOM_GITHUB_OWNER_ID", str(OWNER_ID + 1))
    assert remote.get("/worlds/test-world/summary").status_code == 401


def test_logout_clears_the_session(monkeypatch, remote):
    _sign_in(remote, monkeypatch)
    resp = remote.post("/logout", headers={"Origin": "https://testserver"})
    assert resp.status_code == 303 and resp.headers["location"] == "/view"
    assert auth.SESSION_COOKIE not in remote.cookies
    assert remote.get("/worlds/test-world/summary").status_code == 401


def test_cross_site_writes_are_refused(monkeypatch, remote):
    _sign_in(remote, monkeypatch)
    evil = {"Origin": "https://evil.example"}
    assert remote.post("/worlds/test-world/animate/stop", headers=evil).status_code == 403
    assert remote.post("/logout", headers=evil).status_code == 403
    assert remote.post("/worlds/test-world/animate/stop", headers={"Origin": "null"}).status_code == 403
    assert remote.post("/worlds/test-world/animate/stop", headers={"Sec-Fetch-Site": "cross-site"}).status_code == 403
    assert remote.post("/worlds/test-world/animate/stop", headers={"Origin": "https://testserver"}).status_code == 200


def test_dev_frontend_origin_may_write_locally():
    with _local() as client:
        resp = client.post("/worlds/test-world/animate/stop", headers={"Origin": f"http://localhost:{main._frontend_port}"})
        assert resp.status_code == 200
        assert client.post("/worlds/test-world/animate/stop", headers={"Origin": "https://evil.example"}).status_code == 403


def test_public_url_origin_may_write_behind_a_host_rewriting_proxy(monkeypatch, remote):
    """A proxy that sends Host: 127.0.0.1:<port> still lets the owner's own pages write."""
    _sign_in(remote, monkeypatch)
    origin = {"Origin": "https://mantle.example.com"}
    assert remote.post("/worlds/test-world/animate/stop", headers=origin).status_code == 403
    monkeypatch.setenv("MANTLE_BLOOM_PUBLIC_URL", "https://mantle.example.com/")
    assert remote.post("/worlds/test-world/animate/stop", headers=origin).status_code == 200


def test_dev_frontend_origin_may_not_write_remotely(monkeypatch, remote):
    _sign_in(remote, monkeypatch)
    dev = {"Origin": f"http://localhost:{main._frontend_port}"}
    assert remote.post("/worlds/test-world/animate/stop", headers=dev).status_code == 403


def test_tokens_minted_for_another_purpose_are_not_sessions(remote):
    """The guest allowance and OAuth state cookies share the session's secret; neither
    may pass as a session even when it names the owner (#920)."""
    for purpose in ("guest-steps", "oauth-state"):
        token = auth.sign_token({"uid": OWNER_ID, "used": 0}, SECRET.encode(), purpose, 60)
        remote.cookies.set(auth.SESSION_COOKIE, token)
        assert remote.get("/worlds/test-world/summary").status_code == 401, purpose
        assert remote.get("/auth/session").json()["admin"] is False


def _signed(body: bytes) -> str:
    encoded = auth._b64(body)
    return f"{encoded}.{auth._mac(SECRET.encode(), 'session', encoded)}"


@pytest.mark.parametrize("cookie", [
    "", ".", "a.b.c", "not-a-token", "%00.%00", "ZZZ.ZZZ",
    _signed(b"not json"), _signed(b"[]"), _signed(b'"owner"'),
    _signed(b'{"uid":4242}'), _signed(b'{"uid":4242,"exp":"never"}'), _signed(b'{"uid":4242,"exp":1e99}'),
    _signed(b'{"uid":"4242","exp":99999999999}'), _signed(b'{"uid":true,"exp":99999999999}'),
])
def test_malformed_session_cookies_are_refused_not_errors(cookie):
    """A correctly signed body that isn't a well-formed session is still no session, and
    nothing a client puts in the cookie can cause a 500 (#920)."""
    client = TestClient(app, base_url="https://testserver", follow_redirects=False, raise_server_exceptions=True)
    headers = {"Cookie": f"{auth.SESSION_COOKIE}={cookie}"}
    assert client.get("/worlds/test-world/summary", headers=headers).status_code == 401
    assert client.get("/auth/session", headers=headers).json()["admin"] is False


def test_copied_session_outlives_logout_until_it_expires(monkeypatch, remote):
    """Accepted for a single-owner service (#920 S10): sessions are stateless, so logout
    clears the browser's cookie but can't revoke a copy. Rotating the secret does."""
    _sign_in(remote, monkeypatch)
    copied = remote.cookies.get(auth.SESSION_COOKIE)
    remote.post("/logout", headers={"Origin": "https://testserver"})
    other = TestClient(app, base_url="https://testserver", follow_redirects=False)
    other.cookies.set(auth.SESSION_COOKIE, copied)
    assert other.get("/worlds/test-world/summary").status_code != 401
    monkeypatch.setenv("MANTLE_BLOOM_SESSION_SECRET", "r" * 40)
    assert other.get("/worlds/test-world/summary").status_code == 401


def _all_routes(routes) -> list:
    """Every route with methods, including those behind `include_router`."""
    found = []
    for route in routes:
        inner = getattr(route, "original_router", None)
        if inner is not None:
            found.extend(_all_routes(inner.routes))
        elif hasattr(route, "methods"):
            found.append(route)
    return found


# Every route the simulator app answers without a session. A new route must be added here
# on purpose, or the test below fails: nothing becomes public by accident (#920).
_PUBLIC_ROUTES = {
    "/login", "/auth/github", "/auth/github/callback", "/logout", "/auth/session", "/healthz",
    "/demo", "/demo/", "/openapi.json", "/docs", "/docs/oauth2-redirect", "/redoc",
    "/public/worlds", "/public/worlds/{world_id}/world", "/public/worlds/{world_id}/render",
    "/public/worlds/{world_id}/step",
}


def test_every_route_is_owner_only_or_listed_as_public(remote):
    routes = _all_routes(app.routes)
    guarded = [route for route in routes if auth._is_world_route(route.path)]
    assert {route.path for route in routes} - {route.path for route in guarded} == _PUBLIC_ROUTES
    assert len(guarded) > 40
    for route in guarded:
        path = route.path.replace("{world_id}", "test-world")
        for method in route.methods - {"HEAD", "OPTIONS"}:
            assert remote.request(method, path).status_code == 401, (method, path)


def test_a_finished_sign_in_cannot_be_replayed_in_the_same_browser(monkeypatch, remote):
    """The callback clears the state cookie, so a second callback with the same `state`
    and `code` starts no new session (#920)."""
    _fake_github(monkeypatch, {"id": OWNER_ID, "login": "owner"})
    state = _start(remote)
    params = {"code": "abc", "state": state}
    assert remote.get("/auth/github/callback", params=params).headers["location"] == "/"
    assert auth.STATE_COOKIE not in remote.cookies
    remote.cookies.delete(auth.SESSION_COOKIE)
    again = remote.get("/auth/github/callback", params=params)
    assert again.headers["location"] == "/login?error=state"
    assert auth.SESSION_COOKIE not in remote.cookies
