"""The logged-out viewer's world access and its per-cookie step allowance (#373).

Anyone can open `/view` and look at the shared world, but only through the narrow `/public/*`
routes in main.py: a biome or elevation render in the Eckert IV projection, and a fixed-size
step. Everything that creates, saves, loads or reconfigures a world stays behind the admin
guard (auth.py).

A step is the one expensive thing a guest can ask for, so each guest cookie may trigger at most
`GUEST_STEP_LIMIT` step computations. The count lives in the cookie itself, HMAC-signed so it
can't be edited down; clearing the cookie resets it, which the issue accepts -- this is a usage
limit, not an anti-abuse control.

Only a request that actually runs the simulation is charged. Each step request names the world
revision its viewer is looking at. If the world has already moved past it (another viewer, or
the owner, stepped since), or another guest's step is already running, the request is answered
with the resulting world without computing anything, and costs nothing.
"""

from __future__ import annotations

import secrets
import threading
from dataclasses import dataclass
from typing import Callable

from fastapi import Request, Response

from . import auth

GUEST_STEP_LIMIT = 30
# Fixed server-side, so a guest can't buy a much longer (costlier) step with one allowance unit.
# The full app's default "Years per step".
GUEST_STEP_YEARS = 100_000.0
GUEST_PROJECTION = "eckert4"
GUEST_VIEWS = ("elevation", "biome")
# The viewer renders at 2x its 1100x611 display size; a little headroom past that, but well
# under the admin cap, since a guest render holds the world lock.
GUEST_MAX_RENDER_PX = 2400
# A guest waits this long for the owner's step or animation to release the world, then gets a
# 503. Blocking without a limit would let a crowd of viewers tie up every worker thread for the
# length of a long animation.
GUEST_LOCK_WAIT_S = 10.0
# How long a guest's step waits on another guest's step that's already running.
GUEST_INFLIGHT_WAIT_S = 300.0

GUEST_COOKIE = "mb_guest"
GUEST_COOKIE_MAX_AGE_S = 365 * 24 * 3600

# Used when no session secret is configured (e.g. a local run). A restart then resets every
# guest's allowance, which is no worse than a guest clearing their cookie.
_FALLBACK_SECRET = secrets.token_bytes(32)


def world_revision(world) -> str:
    """Identifies the world state a viewer is looking at. Seed plus elapsed time changes on
    every step and on almost every generate or load; a rare collision only means a guest is
    charged for a step they really did ask for, never that two computations are merged."""
    return f"{world.seed}:{float(world.elapsed_years)!r}"


def public_summary(world) -> dict:
    """The little a viewer needs: the world summary is mostly inspector data the viewer
    doesn't show."""
    return {"revision": world_revision(world), "seed": world.seed, "elapsed_years": world.elapsed_years}


def _secret(config: auth.AuthConfig) -> bytes:
    if len(config.session_secret) >= auth.MIN_SESSION_SECRET_BYTES:
        return config.session_secret
    return _FALLBACK_SECRET


def steps_used(request: Request, config: auth.AuthConfig) -> int:
    payload = auth.read_token(request.cookies.get(GUEST_COOKIE), _secret(config), "guest-steps")
    used = payload.get("used") if payload else None
    return used if isinstance(used, int) and used >= 0 else 0


def set_steps_used(response: Response, used: int, config: auth.AuthConfig) -> None:
    token = auth.sign_token({"used": used}, _secret(config), "guest-steps", GUEST_COOKIE_MAX_AGE_S)
    # Not `Secure`: this cookie grants nothing, and a viewer reached over plain HTTP should
    # still keep their count rather than silently getting a fresh allowance on every step.
    response.set_cookie(
        GUEST_COOKIE, token, max_age=GUEST_COOKIE_MAX_AGE_S, path="/", httponly=True, samesite="lax"
    )


@dataclass(frozen=True)
class StepOutcome:
    # "computed": this request ran the step. "shared": the world had already moved past the
    # viewer's revision, so it was served without computing. "exhausted": no allowance left
    # for a step that would need computing. "busy": the owner holds the world, or the step
    # we waited on didn't finish in time.
    status: str


class GuestStepCoordinator:
    """Runs at most one guest step at a time and lets concurrent requests share it."""

    def __init__(self) -> None:
        self._mutex = threading.Lock()
        self._inflight: threading.Event | None = None

    def request_step(
        self,
        seen_revision: str,
        may_compute: bool,
        current_world: Callable[[], object],
        world_lock: threading.Lock,
        compute: Callable[[object], None],
        inflight_wait_s: float = GUEST_INFLIGHT_WAIT_S,
    ) -> StepOutcome:
        while True:
            with self._mutex:
                inflight = self._inflight
                if inflight is None:
                    world = current_world()
                    if world_revision(world) != seen_revision:
                        return StepOutcome("shared")
                    if not may_compute:
                        return StepOutcome("exhausted")
                    # Never block here: the owner's step or animation can hold the lock for
                    # minutes, and we're inside `_mutex`.
                    if not world_lock.acquire(blocking=False):
                        return StepOutcome("busy")
                    done = self._inflight = threading.Event()
                    break
            # Another guest's step is running. Wait for it, then decide again: normally the
            # world has moved past `seen_revision` and this request is served for free.
            if not inflight.wait(inflight_wait_s):
                return StepOutcome("busy")
        try:
            compute(world)
        finally:
            world_lock.release()
            with self._mutex:
                self._inflight = None
            done.set()
        return StepOutcome("computed")
