"""Single-owner GitHub login guarding the full app (#372).

The full app at `/`, every `/worlds/{world_id}/*` route, and `/debug/scenarios` are for
one person: the GitHub account whose numeric user ID is in
`MANTLE_BLOOM_GITHUB_OWNER_ID`. The ID is stable where the login name
isn't -- a renamed or recycled login can't inherit admin. Everyone else is sent to `/view`.

Sessions are stateless: an HMAC-signed cookie carrying the owner's ID and an expiry, so there's
no session store to add. Rotating `MANTLE_BLOOM_SESSION_SECRET` signs everyone out. CSRF is
refused by origin: a state-changing request whose `Origin` isn't this server's own (or the
dev frontend's) is rejected before it reaches a route. The session cookie is also
`SameSite=Lax`, so a cross-site form post never carries it in the first place.

Requests made on this machine -- a loopback client, a loopback `Host`, and no proxy forwarding
headers -- skip login (the #372 thread), so `bin/restart.sh` and the desktop build keep working
with no OAuth app configured. A reverse proxy on the same host connects from loopback too, so
the `Host` and forwarding-header checks are what keep its traffic out of the bypass. Set
`MANTLE_BLOOM_LOCAL_NO_LOGIN=0` behind a proxy that rewrites `Host` and strips those headers.

Without a complete OAuth configuration, the guard fails closed: remote requests get no admin
access, and `/login` says login isn't configured.
"""

from __future__ import annotations

import base64
import hashlib
import hmac
import html
import ipaddress
import json
import os
import secrets
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass

from fastapi import APIRouter, Request
from fastapi.responses import HTMLResponse, JSONResponse, RedirectResponse, Response

SESSION_COOKIE = "mb_session"
STATE_COOKIE = "mb_oauth_state"
SESSION_MAX_AGE_S = 14 * 24 * 3600
# The authorize-and-return round trip happens in seconds; ten minutes is generous.
STATE_MAX_AGE_S = 10 * 60
# Shorter secrets make the HMAC guessable offline from a single cookie.
MIN_SESSION_SECRET_BYTES = 32

GITHUB_AUTHORIZE_URL = "https://github.com/login/oauth/authorize"
GITHUB_TOKEN_URL = "https://github.com/login/oauth/access_token"
GITHUB_USER_URL = "https://api.github.com/user"
GITHUB_TIMEOUT_S = 10

_SAFE_METHODS = frozenset({"GET", "HEAD", "OPTIONS"})
_LOOPBACK_HOSTNAMES = frozenset({"localhost", "127.0.0.1", "::1"})
# Any of these on a loopback request means a proxy relayed it from elsewhere. Every
# `X-Forwarded-*` header counts, not just -For/-Host: a proxy that rewrites Host to localhost
# may forward only -Proto or -Port (PR #377 review).
_PROXY_HEADERS = frozenset({
    "forwarded", "via", "x-real-ip", "x-client-ip", "true-client-ip", "cf-connecting-ip",
    "fastly-client-ip", "x-original-forwarded-for",
})
_PROXY_HEADER_PREFIX = "x-forwarded"

# /login?error=<code> only ever shows one of these, so the page never reflects input.
LOGIN_ERRORS = {
    "denied": "GitHub sign-in was cancelled.",
    "state": "That sign-in attempt expired or was already used. Please try again.",
    "github": "Couldn't reach GitHub to finish signing in. Please try again.",
    "forbidden": "That GitHub account isn't allowed to administer this server.",
    "unconfigured": "GitHub login isn't configured on this server.",
}


@dataclass(frozen=True)
class AuthConfig:
    client_id: str
    client_secret: str
    owner_id: int | None
    session_secret: bytes
    public_url: str
    local_no_login: bool

    @property
    def login_configured(self) -> bool:
        return bool(
            self.client_id
            and self.client_secret
            and self.owner_id is not None
            and len(self.session_secret) >= MIN_SESSION_SECRET_BYTES
        )


def load_config() -> AuthConfig:
    """Read the auth settings from the environment on every call, so tests and a restarted
    process see changes without a reload."""
    owner = os.environ.get("MANTLE_BLOOM_GITHUB_OWNER_ID", "").strip()
    return AuthConfig(
        client_id=os.environ.get("MANTLE_BLOOM_GITHUB_CLIENT_ID", "").strip(),
        client_secret=os.environ.get("MANTLE_BLOOM_GITHUB_CLIENT_SECRET", "").strip(),
        owner_id=int(owner) if owner.isdigit() else None,
        session_secret=os.environ.get("MANTLE_BLOOM_SESSION_SECRET", "").encode(),
        public_url=os.environ.get("MANTLE_BLOOM_PUBLIC_URL", "").strip().rstrip("/"),
        local_no_login=os.environ.get("MANTLE_BLOOM_LOCAL_NO_LOGIN", "1").strip() != "0",
    )


def _b64(raw: bytes) -> str:
    return base64.urlsafe_b64encode(raw).rstrip(b"=").decode()


def _unb64(text: str) -> bytes:
    return base64.urlsafe_b64decode(text + "=" * (-len(text) % 4))


def _same(a: str, b: str) -> bool:
    # compare_digest raises on non-ASCII str; a crafted cookie or query value must just fail.
    return hmac.compare_digest(a.encode(), b.encode())


def _mac(secret: bytes, purpose: str, body: str) -> str:
    # `purpose` keeps a token minted for one cookie from validating as the other.
    return _b64(hmac.new(secret, f"{purpose}.{body}".encode(), hashlib.sha256).digest())


def sign_token(payload: dict, secret: bytes, purpose: str, max_age_s: int) -> str:
    body = _b64(json.dumps({**payload, "exp": int(time.time()) + max_age_s}, separators=(",", ":")).encode())
    return f"{body}.{_mac(secret, purpose, body)}"


def read_token(token: str | None, secret: bytes, purpose: str) -> dict | None:
    """The payload of a token `sign_token` made with the same secret and purpose, or None if
    it's missing, forged, malformed or expired."""
    if not token or not secret or token.count(".") != 1:
        return None
    body, mac = token.split(".")
    if not _same(mac, _mac(secret, purpose, body)):
        return None
    try:
        payload = json.loads(_unb64(body))
    except ValueError:
        return None
    if not isinstance(payload, dict) or not isinstance(payload.get("exp"), int) or payload["exp"] < time.time():
        return None
    return payload


def is_local_request(request: Request, config: AuthConfig) -> bool:
    if not config.local_no_login or request.client is None:
        return False
    try:
        if not ipaddress.ip_address(request.client.host).is_loopback:
            return False
    except ValueError:
        return False
    if request.url.hostname not in _LOOPBACK_HOSTNAMES:
        return False
    return not any(_is_proxy_header(name) for name in request.headers.keys())


def _is_proxy_header(name: str) -> bool:
    # Starlette lowercases header names.
    return name in _PROXY_HEADERS or name.startswith(_PROXY_HEADER_PREFIX)


def session_user(request: Request, config: AuthConfig) -> dict | None:
    """The signed-in owner's session payload, or None. The ID is checked against the current
    owner, not just the signature, so changing the owner ID revokes old sessions."""
    if not config.login_configured:
        return None
    payload = read_token(request.cookies.get(SESSION_COOKIE), config.session_secret, "session")
    if payload is None or payload.get("uid") != config.owner_id:
        return None
    return payload


def is_admin(request: Request, config: AuthConfig) -> bool:
    return is_local_request(request, config) or session_user(request, config) is not None


def _origin_of(url: str) -> str:
    parts = urllib.parse.urlsplit(url)
    return f"{parts.scheme}://{parts.netloc}"


def is_cross_site_write(request: Request, config: AuthConfig, extra_origins: frozenset[str]) -> bool:
    """Whether a state-changing request came from another site's page. Browsers attach
    `Origin` to every POST, so a write without one is a non-browser client (curl, a script),
    which can't ride a victim's cookie anyway. `Sec-Fetch-Site` covers older browsers that
    omit `Origin` on same-origin posts.

    Same-site means the `Host` the request arrived with, or `MANTLE_BLOOM_PUBLIC_URL` for a
    proxy that rewrites `Host`. The dev frontend's origins count only for local requests, so
    a page on the owner's own localhost can't post to a remote deployment (PR #377 review)."""
    if request.method in _SAFE_METHODS:
        return False
    origin = request.headers.get("origin")
    if origin is None:
        return request.headers.get("sec-fetch-site") == "cross-site"
    if config.public_url and origin == _origin_of(config.public_url):
        return False
    if origin in extra_origins and is_local_request(request, config):
        return False
    return urllib.parse.urlsplit(origin).netloc != request.headers.get("host")


def _route_path(scope) -> tuple[str, str]:
    """(root_path, the path the router matches). Starlette keeps root_path inside
    scope["path"] (uvicorn --root-path, a mount) and strips it when routing, so the guard must
    strip it too or /<prefix>/worlds/{world_id}/... slips past (PR #377 review)."""
    path, root = scope["path"], scope.get("root_path", "")
    if root and (path == root or path.startswith(root + "/")):
        return root, path[len(root):] or "/"
    return root, path


def _is_world_route(path: str) -> bool:
    # /debug/scenarios is guarded too: it moved out from under /worlds/ but is
    # still an owner-route, and the guard below also carries the cross-site
    # write check for it.
    return path == "/worlds" or path.startswith("/worlds/") or path.startswith("/debug/")


def gate(request: Request, config: AuthConfig, extra_origins: frozenset[str]) -> Response | None:
    """The response that refuses `request`, or None to let it through."""
    if is_cross_site_write(request, config, extra_origins):
        return JSONResponse({"detail": "Cross-site request refused"}, status_code=403)
    root, path = _route_path(request.scope)
    if path == "/" and request.method in ("GET", "HEAD") and not is_admin(request, config):
        return RedirectResponse(f"{root}/view", status_code=302)
    if _is_world_route(path) and not is_admin(request, config):
        return JSONResponse({"detail": "Login required"}, status_code=401)
    return None


class AuthMiddleware:
    """Pure ASGI rather than BaseHTTPMiddleware, so the streamed /worlds/{world_id}/generate,
    /worlds/{world_id}/step and /worlds/{world_id}/animate responses pass through untouched."""

    def __init__(self, app, extra_origins: frozenset[str] = frozenset()):
        self.app = app
        self.extra_origins = extra_origins

    async def __call__(self, scope, receive, send):
        if scope["type"] == "http":
            refusal = gate(Request(scope), load_config(), self.extra_origins)
            if refusal is not None:
                await refusal(scope, receive, send)
                return
        await self.app(scope, receive, send)


class GitHubError(Exception):
    pass


def _http_json(req: urllib.request.Request) -> dict:
    try:
        with urllib.request.urlopen(req, timeout=GITHUB_TIMEOUT_S) as resp:
            payload = json.load(resp)
    except (urllib.error.URLError, TimeoutError, ValueError) as exc:
        raise GitHubError(str(exc)) from exc
    if not isinstance(payload, dict):
        raise GitHubError("unexpected GitHub response")
    return payload


def fetch_github_user(config: AuthConfig, code: str, verifier: str, redirect_uri: str) -> dict:
    """Trade an authorization code for the GitHub user it was issued to."""
    form = urllib.parse.urlencode({
        "client_id": config.client_id,
        "client_secret": config.client_secret,
        "code": code,
        "redirect_uri": redirect_uri,
        "code_verifier": verifier,
    }).encode()
    token = _http_json(urllib.request.Request(
        GITHUB_TOKEN_URL, data=form, headers={"Accept": "application/json"}, method="POST",
    ))
    access_token = token.get("access_token")
    if not isinstance(access_token, str) or not access_token:
        raise GitHubError(f"token exchange failed: {token.get('error', 'no access_token')}")
    return _http_json(urllib.request.Request(GITHUB_USER_URL, headers={
        "Accept": "application/vnd.github+json",
        "Authorization": f"Bearer {access_token}",
        "User-Agent": "mantle-bloom",
    }))


def _set_cookie(response: Response, name: str, value: str, max_age: int, path: str = "/") -> None:
    response.set_cookie(name, value, max_age=max_age, path=path, secure=True, httponly=True, samesite="lax")


def _delete_cookie(response: Response, name: str, path: str = "/") -> None:
    response.delete_cookie(name, path=path, secure=True, httponly=True, samesite="lax")


def _callback_url(request: Request, config: AuthConfig) -> str:
    # Behind a TLS-terminating proxy, request.base_url says http://; MANTLE_BLOOM_PUBLIC_URL
    # pins the URL that's registered with the GitHub OAuth app.
    base = config.public_url or str(request.base_url).rstrip("/")
    return f"{base}/auth/github/callback"


def _login_failure(code: str) -> RedirectResponse:
    response = RedirectResponse(f"/login?error={code}", status_code=303)
    _delete_cookie(response, STATE_COOKIE, path="/auth/github")
    return response


_LOGIN_PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sign in - Mantle-Bloom</title>
<style>
body {{ margin: 0; min-height: 100vh; display: grid; place-items: center; background: #0f1117;
  color: #e6e8ef; font-family: system-ui, sans-serif; }}
main {{ width: min(360px, calc(100vw - 32px)); text-align: center; }}
h1 {{ font-size: 22px; margin: 0 0 4px; }}
p {{ opacity: 0.75; }}
.error {{ opacity: 1; color: #ff9b9b; }}
a.button {{ display: inline-block; margin: 16px 0; padding: 10px 18px; border-radius: 6px;
  background: #e6e8ef; color: #0f1117; text-decoration: none; font-weight: 600; }}
a {{ color: #9db8ff; }}
</style></head>
<body><main>
<h1>Mantle-Bloom</h1>
<p>Sign in to run the full simulator.</p>
{error}{action}{help}
<p><a href="/view">Back to the viewer</a></p>
</main></body></html>
"""


router = APIRouter()


# Operator help for the unconfigured login page: variable names only, no values --
# they're already public in docs/hosting.md. Without this, the admin's only symptom
# on a fresh Railway deploy is a /view with no Sign in entry (#626).
_UNCONFIGURED_HELP_HTML = """<p>If you run this server: set MANTLE_BLOOM_GITHUB_CLIENT_ID,
MANTLE_BLOOM_GITHUB_CLIENT_SECRET, MANTLE_BLOOM_GITHUB_OWNER_ID, MANTLE_BLOOM_SESSION_SECRET
and MANTLE_BLOOM_PUBLIC_URL, and register &lt;public-url&gt;/auth/github/callback with your
GitHub OAuth app (see docs/hosting.md). The /view header's Sign in link lands here until
login is configured.</p>"""


@router.get("/login", include_in_schema=False)
def login_page(request: Request, error: str | None = None, demo: bool = False) -> Response:
    config = load_config()
    if is_admin(request, config) and not demo:
        return RedirectResponse("/", status_code=303)
    if not config.login_configured:
        error = "unconfigured"
    message = LOGIN_ERRORS.get(error or "")
    error_html = f'<p class="error" role="alert">{html.escape(message)}</p>' if message else ""
    action = '<a class="button" href="/auth/github">Sign in with GitHub</a>' if config.login_configured else ""
    help_html = _UNCONFIGURED_HELP_HTML if not config.login_configured else ""
    return HTMLResponse(_LOGIN_PAGE.format(error=error_html, action=action, help=help_html))


@router.get("/auth/github", include_in_schema=False)
def github_start(request: Request) -> Response:
    config = load_config()
    if not config.login_configured:
        return RedirectResponse("/login?error=unconfigured", status_code=303)
    state = secrets.token_urlsafe(24)
    # PKCE, so a leaked authorization code is useless without this cookie's verifier.
    verifier = secrets.token_urlsafe(48)
    challenge = _b64(hashlib.sha256(verifier.encode()).digest())
    query = urllib.parse.urlencode({
        "client_id": config.client_id,
        "redirect_uri": _callback_url(request, config),
        "state": state,
        "scope": "",
        "allow_signup": "false",
        "code_challenge": challenge,
        "code_challenge_method": "S256",
    })
    response = RedirectResponse(f"{GITHUB_AUTHORIZE_URL}?{query}", status_code=303)
    token = sign_token({"state": state, "verifier": verifier}, config.session_secret, "oauth-state", STATE_MAX_AGE_S)
    _set_cookie(response, STATE_COOKIE, token, STATE_MAX_AGE_S, path="/auth/github")
    return response


@router.get("/auth/github/callback", include_in_schema=False)
def github_callback(
    request: Request, code: str | None = None, state: str | None = None, error: str | None = None
) -> Response:
    config = load_config()
    if not config.login_configured:
        return _login_failure("unconfigured")
    if error:
        return _login_failure("denied")
    expected = read_token(request.cookies.get(STATE_COOKIE), config.session_secret, "oauth-state")
    if not code or not state or expected is None or not _same(state, str(expected.get("state"))):
        return _login_failure("state")
    try:
        user = fetch_github_user(config, code, str(expected.get("verifier")), _callback_url(request, config))
    except GitHubError:
        return _login_failure("github")
    if user.get("id") != config.owner_id:
        return _login_failure("forbidden")
    response = RedirectResponse("/", status_code=303)
    _delete_cookie(response, STATE_COOKIE, path="/auth/github")
    session = sign_token({"uid": config.owner_id, "login": str(user.get("login", ""))},
                         config.session_secret, "session", SESSION_MAX_AGE_S)
    _set_cookie(response, SESSION_COOKIE, session, SESSION_MAX_AGE_S)
    return response


@router.post("/logout", include_in_schema=False)
def logout() -> Response:
    # AuthMiddleware has already refused a cross-site post here.
    response = RedirectResponse("/view", status_code=303)
    _delete_cookie(response, SESSION_COOKIE)
    return response


@router.get("/auth/session")
def auth_session(request: Request) -> dict:
    config = load_config()
    session = session_user(request, config)
    local = is_local_request(request, config)
    return {
        "admin": local or session is not None,
        "local": local,
        "login": session.get("login") if session else None,
        "login_available": config.login_configured,
    }
