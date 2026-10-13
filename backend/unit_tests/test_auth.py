import time
import urllib.parse

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
    """No world and no pending animation stop, before and after: the CSRF tests post
    /world/animate/stop, which would otherwise leak a stop signal into later tests."""
    main._state["world"] = None
    main._animation_stop_event.clear()
    yield
    main._state["world"] = None
    main._animation_stop_event.clear()


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


def test_expired_token_is_rejected(monkeypatch):
    token = auth.sign_token({"uid": 1}, b"k" * 32, "session", 60)
    monkeypatch.setattr(time, "time", lambda: 10**12)
    assert auth.read_token(token, b"k" * 32, "session") is None


def test_logged_out_root_redirects_to_view(remote):
    resp = remote.get("/")
    assert resp.status_code == 302
    assert resp.headers["location"] == "/view"


def test_world_routes_reject_unauthenticated_requests(remote):
    assert remote.get("/world/summary").status_code == 401
    assert remote.post("/world/step", json={"years": 1}).status_code == 401
    assert remote.post("/world/generate", json={}).status_code == 401
    assert remote.post("/world/load", content=b"x").status_code == 401


def test_view_and_login_stay_public(remote):
    assert remote.get("/login").status_code == 200
    assert remote.get("/auth/session").json() == {
        "admin": False, "local": False, "login": None, "login_available": True,
    }


def test_unconfigured_server_fails_closed(monkeypatch, remote):
    for name in ("CLIENT_ID", "CLIENT_SECRET", "OWNER_ID"):
        monkeypatch.delenv(f"MANTLE_BLOOM_GITHUB_{name}", raising=False)
    assert remote.get("/world/summary").status_code == 401
    page = remote.get("/login")
    assert "isn&#x27;t configured" in page.text and "/auth/github" not in page.text
    assert remote.get("/auth/github").headers["location"] == "/login?error=unconfigured"


def test_short_session_secret_leaves_login_unconfigured(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_SESSION_SECRET", "short")
    assert not auth.load_config().login_configured


def test_local_requests_skip_login():
    with _local() as client:
        assert client.get("/world/summary").status_code == 404  # no world yet, but let through
        assert client.get("/auth/session").json()["admin"] is True


@pytest.mark.parametrize("headers", [{"X-Forwarded-For": "203.0.113.9"}, {"Forwarded": "for=203.0.113.9"},
                                     {"Host": "mantle.example.com"}])
def test_proxied_loopback_requests_need_login(headers):
    with _local() as client:
        assert client.get("/world/summary", headers=headers).status_code == 401


def test_local_bypass_can_be_disabled(monkeypatch):
    monkeypatch.setenv("MANTLE_BLOOM_LOCAL_NO_LOGIN", "0")
    with _local() as client:
        assert client.get("/world/summary").status_code == 401


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
    assert remote.get("/world/summary").status_code == 404  # authorized; just no world yet
    assert remote.get("/").status_code != 302
    assert remote.get("/auth/session").json() == {
        "admin": True, "local": False, "login": "owner", "login_available": True,
    }
    assert remote.get("/login").headers["location"] == "/"


def test_other_github_account_gets_no_session(monkeypatch, remote):
    resp = _sign_in(remote, monkeypatch, user_id=OWNER_ID + 1)
    assert resp.headers["location"] == "/login?error=forbidden"
    assert auth.SESSION_COOKIE not in remote.cookies
    assert remote.get("/world/summary").status_code == 401


@pytest.mark.parametrize("params, error", [
    ({"error": "access_denied", "state": "x"}, "denied"),
    ({"code": "abc", "state": "wrong"}, "state"),
    ({"code": "abc"}, "state"),
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
    assert remote.get("/world/summary").status_code == 401
    # A session minted for a previous owner ID stops working when the owner changes.
    remote.cookies.set(auth.SESSION_COOKIE, auth.sign_token({"uid": OWNER_ID}, SECRET.encode(), "session", 60))
    monkeypatch.setenv("MANTLE_BLOOM_GITHUB_OWNER_ID", str(OWNER_ID + 1))
    assert remote.get("/world/summary").status_code == 401


def test_logout_clears_the_session(monkeypatch, remote):
    _sign_in(remote, monkeypatch)
    resp = remote.post("/logout", headers={"Origin": "https://testserver"})
    assert resp.status_code == 303 and resp.headers["location"] == "/view"
    assert auth.SESSION_COOKIE not in remote.cookies
    assert remote.get("/world/summary").status_code == 401


def test_cross_site_writes_are_refused(monkeypatch, remote):
    _sign_in(remote, monkeypatch)
    evil = {"Origin": "https://evil.example"}
    assert remote.post("/world/animate/stop", headers=evil).status_code == 403
    assert remote.post("/logout", headers=evil).status_code == 403
    assert remote.post("/world/animate/stop", headers={"Origin": "null"}).status_code == 403
    assert remote.post("/world/animate/stop", headers={"Sec-Fetch-Site": "cross-site"}).status_code == 403
    assert remote.post("/world/animate/stop", headers={"Origin": "https://testserver"}).status_code == 200


def test_dev_frontend_origin_may_write_locally():
    with _local() as client:
        resp = client.post("/world/animate/stop", headers={"Origin": f"http://localhost:{main._frontend_port}"})
        assert resp.status_code == 200
        assert client.post("/world/animate/stop", headers={"Origin": "https://evil.example"}).status_code == 403
