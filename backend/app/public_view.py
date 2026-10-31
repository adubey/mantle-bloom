"""The logged-out viewer's world access and its per-cookie step allowance (#373).

Anyone can open `/view` and look at the shared world, but only through the narrow `/public/*`
routes in main.py: a map render in the Eckert IV projection, and a fixed-size
step. Everything that creates, saves, loads or reconfigures a world stays behind the admin
guard (auth.py).

A step is the one expensive thing a guest can ask for, so each guest cookie may trigger at most
`GUEST_STEP_LIMIT` step computations. The count lives in the cookie itself, HMAC-signed so it
can't be edited down; clearing the cookie resets it, which the issue accepts -- this is a usage
limit, not an anti-abuse control.

Only a request that actually runs the simulation is charged. Each step request names the world
revision its viewer is looking at. If the world has already moved past it (another viewer, or
the owner, stepped since), or another guest's step is already running, the request is answered
with the resulting world without computing anything, and costs nothing. A request refused
because the owner holds the world waits up to `GUEST_OWNER_STEP_WAIT_S` for their step to
land and shares that too: viewers asking for the same advance get one result however it was
computed (#653).

Every viewer shares the one world, so the work of following it is shared too (#374): each
change to the world publishes a new revision through `WorldFeed`, which viewers long-poll
instead of each asking for steps or renders of their own, and `RenderCache` renders each view
once per revision however many viewers ask for it. A computing guest step also publishes its
live (fraction, phase) to the session's progress slot, waking the long poll without a new
revision, which is what drives /view's wait progress bar mid-step; each animation frame
publishes its revision the same way, so viewers watch the owner's recording step by step
instead of nothing until it ends (#653).
"""

from __future__ import annotations

import asyncio
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
GUEST_VIEWS = ("elevation", "biome", "combined")
# The image /view itself asks for (Viewer.tsx's fixed combined view at
# MAX_DISPLAY_WIDTH * RENDER_SCALE): the one PNG the server pre-renders into the
# RenderCache on every published change, so the first viewer after a step hits
# cache instantly (#658). Kept as server constants -- not negotiated per request --
# so there is exactly one eager image per revision.
GUEST_DEFAULT_VIEW = "combined"
GUEST_DEFAULT_WIDTH = 2200
GUEST_DEFAULT_HEIGHT = 1222
# The viewer renders at 2x its 1100x611 display size; a little headroom past that, but well
# under the admin cap, since a guest render holds the world lock.
GUEST_MAX_RENDER_PX = 2400
# A guest waits this long for the owner's step or animation to release the world, then gets a
# 503. Blocking without a limit would let a crowd of viewers tie up every worker thread for the
# length of a long animation.
GUEST_LOCK_WAIT_S = 10.0
# How long a guest's step waits on another guest's step that's already running.
GUEST_INFLIGHT_WAIT_S = 300.0
# How long a guest refused as busy waits for the owner's step or animation frame to
# land and shares it, before failing with 503 (#653). Owner steps are the same
# advance the viewer asked for, so sharing beats refusing; past this window the
# request fails rather than pin a worker thread for the whole owner run.
GUEST_OWNER_STEP_WAIT_S = 60.0
# Poll granularity of that wait: the feed revision is a cheap mutex read.
GUEST_OWNER_STEP_POLL_S = 0.5
# How many guest requests may queue behind an already-running guest step. Each waiter
# pins a worker thread for up to GUEST_INFLIGHT_WAIT_S, and this route is public: without
# a cap, one burst during a slow step can occupy the whole AnyIO worker pool and starve
# the owner's own synchronous routes. Excess waiters get a 503 at once (the viewer's
# poll backs off and retries). Sized to leave headroom in the default 40-thread pool.
GUEST_MAX_WAITERS = 8
# How long GET /public/world?after=... holds a viewer's request open waiting for a new
# revision. Kept under the idle timeout of common proxies, so the request ends on our terms.
GUEST_POLL_WAIT_S = 25.0
# Renders kept for the current revision. The viewer asks for one combined render size; the cap
# stops a client cycling through sizes from holding many large PNGs.
GUEST_RENDER_CACHE_ENTRIES = 6

GUEST_COOKIE = "mb_guest"
GUEST_COOKIE_MAX_AGE_S = 365 * 24 * 3600

# Used when no session secret is configured (e.g. a local run). A restart then resets every
# guest's allowance, which is no worse than a guest clearing their cookie.
_FALLBACK_SECRET = secrets.token_bytes(32)


class WorldFeed:
    """The shared world's revision, and the viewer-facing summary of it, published after every
    change to the world (#374).

    The revision is a counter, bumped on every change, behind a per-process epoch so a restart
    can't hand out a revision a viewer already saw for a different world. The summary is taken
    at publish time, while the caller still holds the world lock, so readers get a consistent
    revision/elapsed-time pair without taking the lock, which a step may hold for minutes.
    """

    def __init__(self) -> None:
        self._mutex = threading.Lock()
        self._epoch = secrets.token_hex(4)
        self._count = 0
        self._summary: dict | None = None
        self._name: str | None = None
        self._waiters: list[tuple[asyncio.AbstractEventLoop, asyncio.Future]] = []

    def revision(self) -> str:
        with self._mutex:
            return self._revision()

    def _revision(self) -> str:
        return f"{self._epoch}-{self._count}"

    def summary(self) -> dict | None:
        """{"revision", "seed", "elapsed_years"} of the current world, or None before one exists."""
        with self._mutex:
            return self._summary

    def name(self) -> str | None:
        """The name viewers see. Published here rather than read off the session so it moves
        with the summary: a regenerate renames the session before its new world exists, and
        the new name must not appear beside the old world's age (#711)."""
        with self._mutex:
            return self._name

    def set_name(self, name: str) -> None:
        """Rename the world for viewers, waking them: a rename has no new revision (#711)."""
        with self._mutex:
            self._name = name
        self.wake()

    def publish(self, world, name: str | None = None) -> None:
        """Record a change to `world` (None: there's no world), and to its `name` if given,
        and wake every waiting viewer."""
        with self._mutex:
            self._count += 1
            if name is not None:
                self._name = name
            revision = self._revision()
            self._summary = None if world is None else {
                "revision": revision, "seed": world.seed, "elapsed_years": float(world.elapsed_years)
            }
            waiters, self._waiters = self._waiters, []
        _wake_all(waiters)

    def wake(self) -> None:
        """Wake every waiting request without changing the revision, for a change to something
        a `wait_past` caller's `changed` checks (the run status, #375)."""
        with self._mutex:
            waiters, self._waiters = self._waiters, []
        _wake_all(waiters)

    def has_waiters(self) -> bool:
        """A viewer is following this world right now. Gates the viewer preview
        pre-render (#658): with an audience the render pays for itself on the next
        poll; without one it is pure cost on every change. Best-effort by nature --
        a waiter arriving mid-render just renders on demand as before."""
        with self._mutex:
            return bool(self._waiters)

    async def wait_past(self, seen: str, timeout_s: float, changed: Callable[[], bool] | None = None) -> None:
        """Return once the revision differs from `seen` or `changed()` is true, or after
        `timeout_s`. Whatever `changed` reads must be updated before the `wake` that announces
        it: both run under the same mutex, so no update falls between the check and the wait."""
        loop = asyncio.get_running_loop()
        future = loop.create_future()
        with self._mutex:
            if self._revision() != seen or (changed is not None and changed()):
                return
            self._waiters.append((loop, future))
        try:
            await asyncio.wait_for(future, timeout_s)
        except asyncio.TimeoutError:
            pass
        finally:
            with self._mutex:
                if (loop, future) in self._waiters:
                    self._waiters.remove((loop, future))


def _wake_all(waiters: list[tuple[asyncio.AbstractEventLoop, asyncio.Future]]) -> None:
    for loop, future in waiters:
        try:
            loop.call_soon_threadsafe(_resolve, future)
        except RuntimeError:
            pass  # That request's event loop has closed; there's no one left to wake.


def _resolve(future: asyncio.Future) -> None:
    if not future.done():
        future.set_result(None)


class RenderCache:
    """Renders of the current revision, keyed by the viewer's (view, width, height); the
    projection and orientation are fixed for viewers. An entry for an older revision is never
    served: storing one for a newer revision drops them all, and `clear` frees them as soon as
    the world changes."""

    def __init__(self, max_entries: int = GUEST_RENDER_CACHE_ENTRIES) -> None:
        self._mutex = threading.Lock()
        self._max_entries = max_entries
        self._revision: str | None = None
        self._entries: dict[tuple, dict] = {}

    def get(self, revision: str, key: tuple) -> dict | None:
        with self._mutex:
            return self._entries.get(key) if revision == self._revision else None

    def put(self, revision: str, key: tuple, entry: dict) -> None:
        with self._mutex:
            if revision != self._revision:
                self._revision, self._entries = revision, {}
            self._entries[key] = entry
            while len(self._entries) > self._max_entries:
                del self._entries[next(iter(self._entries))]

    def clear(self) -> None:
        with self._mutex:
            self._revision, self._entries = None, {}


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
    """Runs at most one guest step at a time and lets concurrent requests share it. `compute`
    must publish the new revision before it returns, while it still holds the world lock."""

    def __init__(self) -> None:
        self._mutex = threading.Lock()
        self._inflight: threading.Event | None = None
        self._waiter_slots = threading.BoundedSemaphore(GUEST_MAX_WAITERS)

    def busy(self) -> bool:
        """A guest step is computing right now. Lets delete refuse while one is
        in flight instead of detaching the session underneath it."""
        with self._mutex:
            return self._inflight is not None

    def request_step(
        self,
        seen_revision: str,
        may_compute: bool,
        current_revision: Callable[[], str],
        world_lock: threading.Lock,
        compute: Callable[[], None],
        inflight_wait_s: float = GUEST_INFLIGHT_WAIT_S,
    ) -> StepOutcome:
        while True:
            with self._mutex:
                inflight = self._inflight
                if inflight is None:
                    if current_revision() != seen_revision:
                        return StepOutcome("shared")
                    if not may_compute:
                        return StepOutcome("exhausted")
                    # Never block here: the owner's step or animation can hold the lock for
                    # minutes, and we're inside `_mutex`.
                    if not world_lock.acquire(blocking=False):
                        return StepOutcome("busy")
                    # The owner may have changed the world and let go of it since the check
                    # above; don't charge this guest for stepping a world they haven't seen.
                    if current_revision() != seen_revision:
                        world_lock.release()
                        return StepOutcome("shared")
                    done = self._inflight = threading.Event()
                    break
            # Another guest's step is running. Wait for it, then decide again: normally the
            # world has moved past `seen_revision` and this request is served for free.
            # The wait takes a bounded slot first: an unbounded crowd of waiters would
            # pin a worker thread each for up to inflight_wait_s and starve the owner's
            # routes of threads.
            if not self._waiter_slots.acquire(blocking=False):
                return StepOutcome("busy")
            try:
                if not inflight.wait(inflight_wait_s):
                    return StepOutcome("busy")
            finally:
                self._waiter_slots.release()
        try:
            compute()
        finally:
            world_lock.release()
            with self._mutex:
                self._inflight = None
            done.set()
        return StepOutcome("computed")
