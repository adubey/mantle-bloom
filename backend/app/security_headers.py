"""Browser security headers on every HTTP response (#555 F5, tracked in #556).

Defense in depth for a same-origin app: the page, its bundle and the API share one origin
and load nothing from anywhere else, so the policy can be strict.

- `X-Content-Type-Options: nosniff` -- a response is only ever what its content type says.
- `X-Frame-Options: DENY` and CSP `frame-ancestors 'none'` -- no other site can frame the
  app and trick a signed-in owner into clicking through it.
- `Referrer-Policy: same-origin` -- outbound links carry no world IDs or paths.
- `Content-Security-Policy` -- scripts only from this origin, never inline. Styles allow
  inline because React writes `style=` attributes and the login page carries a `<style>`
  block. Images and media allow `data:` and `blob:`, which is how rendered maps and saved
  animations reach the page.
- `Strict-Transport-Security` -- only on a response to an HTTPS request, so plain-HTTP
  localhost is never pinned to HTTPS.

The interactive API docs (`/docs`, `/redoc`) load their bundle from a CDN and start it from
an inline script, so they get every header but the content policy's resource rules.

A header the route (or, in the hosted app, the proxied worker's copy of this middleware)
already set is left alone.
"""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse

from .auth import _route_path

CONTENT_SECURITY_POLICY = "; ".join((
    "default-src 'self'",
    "script-src 'self'",
    "style-src 'self' 'unsafe-inline'",
    "img-src 'self' data: blob:",
    "media-src 'self' blob:",
    "font-src 'self' data:",
    "connect-src 'self'",
    "object-src 'none'",
    "base-uri 'self'",
    "form-action 'self'",
    "frame-ancestors 'none'",
))
# Swagger UI and ReDoc: third-party bundle and inline bootstrap, so only the framing rule.
DOCS_CONTENT_SECURITY_POLICY = "frame-ancestors 'none'"
_DOCS_PATHS = frozenset({"/docs", "/docs/oauth2-redirect", "/redoc"})

# One year, this host only: nothing here knows whether subdomains serve HTTPS.
STRICT_TRANSPORT_SECURITY = "max-age=31536000"

_ALWAYS = (
    (b"x-content-type-options", b"nosniff"),
    (b"x-frame-options", b"DENY"),
    (b"referrer-policy", b"same-origin"),
)


def _is_https(scope) -> bool:
    """Whether the browser's request was HTTPS. Behind a TLS-terminating proxy the scope
    says `http`, so the proxy's `X-Forwarded-Proto` counts too. Nothing is trusted on its
    word here: a forged header only adds HSTS to a plain-HTTP response, which browsers
    ignore."""
    if scope.get("scheme") == "https":
        return True
    forwarded = dict(scope.get("headers") or ()).get(b"x-forwarded-proto", b"")
    return forwarded.split(b",")[0].strip().lower() == b"https"


def _headers_for(scope) -> list[tuple[bytes, bytes]]:
    _, path = _route_path(scope)
    policy = DOCS_CONTENT_SECURITY_POLICY if path in _DOCS_PATHS else CONTENT_SECURITY_POLICY
    headers = [*_ALWAYS, (b"content-security-policy", policy.encode())]
    if _is_https(scope):
        headers.append((b"strict-transport-security", STRICT_TRANSPORT_SECURITY.encode()))
    return headers


class SecurityHeadersMiddleware:
    """Pure ASGI, like `auth.AuthMiddleware`, so streamed responses pass through untouched.
    Add it last (outermost) so refusals from the auth guard and CORS carry the headers too."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        extra = _headers_for(scope)

        async def send_with_headers(message):
            if message["type"] == "http.response.start":
                headers = list(message.get("headers") or ())
                present = {name.lower() for name, _ in headers}
                headers.extend(pair for pair in extra if pair[0] not in present)
                message = {**message, "headers": headers}
            await send(message)

        await self.app(scope, receive, send_with_headers)


async def _server_error(request: Request, _exc: Exception) -> JSONResponse:
    """The 500 for an unhandled exception. Starlette builds it in `ServerErrorMiddleware`,
    which wraps every user middleware, so `SecurityHeadersMiddleware` never sees it; the
    handler adds the headers itself. The exception is still re-raised after this response,
    so the server logs the traceback as before."""
    return JSONResponse(
        {"detail": "Internal Server Error"}, status_code=500,
        headers={name.decode(): value.decode() for name, value in _headers_for(request.scope)},
    )


def install(app: FastAPI) -> None:
    """Add the headers to every response of `app`, unhandled-exception 500s included."""
    app.add_middleware(SecurityHeadersMiddleware)
    app.add_exception_handler(Exception, _server_error)
