"""Server-side background animation state for one world (#1033).

Decouples recording execution from the initiating POST /worlds/{world_id}/animate
HTTP stream so refreshing the browser does not stop the animation or lose the
resulting MP4. Both the initiating POST stream and reconnected GET
/worlds/{world_id}/animate/stream clients subscribe to the same NDJSON stream.
"""

from __future__ import annotations

import asyncio
from collections import deque
import threading
import uuid
from typing import AsyncIterator, Iterator

# Cap unread preview/heartbeat lines buffered per slow subscriber so a stalled
# reader cannot accumulate hundreds of base64 PNGs in memory while the terminal
# done/error line is always delivered reliably (#1033).
MAX_SUBSCRIBER_BACKLOG = 8


class _Subscriber:
    """Bounded-backlog stream channel for one HTTP subscriber.

    Supports both synchronous reads (`next_line`) and cancellation-aware async
    reads (`next_line_async`) so disconnecting while the worker is silent
    immediately releases the ASGI task without holding a threadpool slot.
    """

    def __init__(self, maxlen: int = MAX_SUBSCRIBER_BACKLOG) -> None:
        self._cond = threading.Condition()
        self._previews: deque[str] = deque(maxlen=maxlen)
        self._terminal: str | None = None
        self._closed: bool = False
        self._async_waiter: tuple[asyncio.AbstractEventLoop, asyncio.Event] | None = None

    def _wake_waiters_locked(self) -> None:
        self._cond.notify()
        waiter = self._async_waiter
        if waiter is not None:
            loop, event = waiter
            if not loop.is_closed():
                loop.call_soon_threadsafe(event.set)

    def push_preview(self, line: str) -> None:
        with self._cond:
            if self._closed:
                return
            self._previews.append(line)
            self._wake_waiters_locked()

    def push_terminal(self, line: str | None) -> None:
        with self._cond:
            if self._closed:
                return
            if line is not None:
                self._terminal = line
            self._closed = True
            self._wake_waiters_locked()

    def _pop_ready_locked(self) -> tuple[bool, str | None]:
        if self._previews:
            return True, self._previews.popleft()
        if self._terminal is not None:
            line = self._terminal
            self._terminal = None
            return True, line
        if self._closed:
            return True, None
        return False, None

    def next_line(self) -> str | None:
        with self._cond:
            while True:
                ready, line = self._pop_ready_locked()
                if ready:
                    return line
                self._cond.wait()

    async def next_line_async(self) -> str | None:
        loop = asyncio.get_running_loop()
        while True:
            with self._cond:
                ready, line = self._pop_ready_locked()
                if ready:
                    return line
                event = asyncio.Event()
                self._async_waiter = (loop, event)
            try:
                await event.wait()
            finally:
                with self._cond:
                    if self._async_waiter is not None and self._async_waiter[1] is event:
                        self._async_waiter = None


class AnimationJob:
    """Thread-safe pub/sub state for a single world's background animation."""

    def __init__(
        self,
        *,
        total: int,
        unbounded: bool,
        projection: str,
        view: str,
        elapsed_years: float = 0.0,
    ) -> None:
        self.job_id: str = uuid.uuid4().hex
        self.total: int = total
        self.unbounded: bool = unbounded
        self.projection: str = projection
        self.view: str = view
        # Per-job stop event initialized before publication so an early Stop
        # can never be erased by startup clearing a shared session event (#1033).
        self.stop_event: threading.Event = threading.Event()
        self._lock = threading.Lock()
        self._active: bool = True
        self._stopping: bool = False
        self._frame: int = 0
        self._elapsed_years: float = elapsed_years
        self._seconds_per_step: float | None = None
        self._stopped_early: bool = False
        self._error: str | None = None
        self._has_result: bool = False
        # Latest progress NDJSON line so a reconnecting tab immediately gets
        # the current frame PNG, stats, frame counter, and elapsed years.
        self._latest_progress_line: str | None = None
        # Terminal "done" or "error" NDJSON line, kept until dismissed so a tab
        # refreshing after completion still surfaces SaveAnimationModal.
        self._terminal_line: str | None = None
        self._subscribers: list[_Subscriber] = []

    def is_active(self) -> bool:
        with self._lock:
            return self._active

    def has_result(self) -> bool:
        with self._lock:
            return self._has_result

    def subscriber_count(self) -> int:
        with self._lock:
            return len(self._subscribers)

    def frame(self) -> int:
        with self._lock:
            return self._frame

    def seconds_per_step(self) -> float | None:
        with self._lock:
            return self._seconds_per_step

    def set_seconds_per_step(self, rate: float) -> None:
        with self._lock:
            self._seconds_per_step = rate

    def summary(self) -> dict:
        """Compact, lock-free-from-world status payload for GET /worlds/{world_id}/run."""
        with self._lock:
            return {
                "job_id": self.job_id,
                "active": self._active,
                "stopping": self._stopping,
                "frame": self._frame,
                "total": self.total,
                "unbounded": self.unbounded,
                "elapsed_years": self._elapsed_years,
                "seconds_per_step": self._seconds_per_step,
                "has_result": self._has_result,
                "stopped_early": self._stopped_early,
                "error": self._error,
                "projection": self.projection,
                "view": self.view,
            }

    def request_stop(self) -> None:
        with self._lock:
            if self._active:
                self._stopping = True
                self.stop_event.set()

    def subscribe(self, *, replay_latest: bool = False) -> _Subscriber:
        """Register a stream consumer.

        `replay_latest=False` is for the initiating POST /animate stream (registered
        before the worker starts so it receives every frame in order).
        `replay_latest=True` is for a reconnecting GET /animate/stream client: if the
        job is still running it immediately receives the latest progress snapshot,
        and if the job already finished it immediately receives the terminal line.
        """
        sub = _Subscriber()
        with self._lock:
            if self._terminal_line is not None:
                sub.push_terminal(self._terminal_line)
                return sub
            if not self._active:
                sub.push_terminal(None)
                return sub
            if replay_latest and self._latest_progress_line is not None:
                sub.push_preview(self._latest_progress_line)
            self._subscribers.append(sub)
        return sub

    def unsubscribe(self, sub: _Subscriber) -> None:
        with self._lock:
            if sub in self._subscribers:
                self._subscribers.remove(sub)

    def iter_stream(self, sub: _Subscriber) -> Iterator[str]:
        """Yield NDJSON lines synchronously for unit tests."""
        try:
            while True:
                line = sub.next_line()
                if line is None:
                    break
                yield line
        finally:
            self.unsubscribe(sub)

    async def iter_stream_async(self, sub: _Subscriber) -> AsyncIterator[str]:
        """Yield NDJSON lines asynchronously so client disconnects release immediately."""
        try:
            while True:
                line = await sub.next_line_async()
                if line is None:
                    break
                yield line
        finally:
            self.unsubscribe(sub)

    def publish_progress(
        self,
        line: str,
        *,
        frame: int,
        elapsed_years: float,
        seconds_per_step: float | None = None,
    ) -> None:
        with self._lock:
            self._frame = frame
            self._elapsed_years = elapsed_years
            if seconds_per_step is not None:
                self._seconds_per_step = seconds_per_step
            self._latest_progress_line = line
            subscribers = list(self._subscribers)
        for sub in subscribers:
            sub.push_preview(line)

    def publish_preview(
        self,
        line: str,
        *,
        elapsed_years: float,
        seconds_per_step: float | None = None,
    ) -> None:
        """Publish an intermediate-step preview PNG without advancing the encoded frame count.

        Drops any stale preview older than the latest published elapsed_years (#1033).
        """
        with self._lock:
            if not self._active or elapsed_years < self._elapsed_years:
                return
            self._elapsed_years = elapsed_years
            if seconds_per_step is not None:
                self._seconds_per_step = seconds_per_step
            self._latest_progress_line = line
            subscribers = list(self._subscribers)
        for sub in subscribers:
            sub.push_preview(line)

    def publish_heartbeat(
        self,
        line: str,
        *,
        elapsed_years: float | None = None,
        seconds_per_step: float | None = None,
    ) -> None:
        with self._lock:
            if elapsed_years is not None:
                self._elapsed_years = elapsed_years
            if seconds_per_step is not None:
                self._seconds_per_step = seconds_per_step
            subscribers = list(self._subscribers)
        for sub in subscribers:
            sub.push_preview(line)

    def finish_done(
        self,
        line: str,
        *,
        stopped_early: bool,
        error: str | None,
        elapsed_years: float,
    ) -> None:
        with self._lock:
            self._active = False
            self._stopping = False
            self._stopped_early = stopped_early
            self._error = error
            self._elapsed_years = elapsed_years
            self._has_result = True
            self._latest_progress_line = None
            self._terminal_line = line
            subscribers = list(self._subscribers)
            self._subscribers.clear()
        for sub in subscribers:
            sub.push_terminal(line)

    def finish_error(self, line: str, *, detail: str) -> None:
        with self._lock:
            self._active = False
            self._stopping = False
            self._error = detail
            self._has_result = False
            self._latest_progress_line = None
            self._terminal_line = line
            subscribers = list(self._subscribers)
            self._subscribers.clear()
        for sub in subscribers:
            sub.push_terminal(line)

    def finish_empty(self) -> None:
        """Close out a job that bailed before producing a result (e.g. deleted world)."""
        with self._lock:
            self._active = False
            self._stopping = False
            self._latest_progress_line = None
            subscribers = list(self._subscribers)
            self._subscribers.clear()
        for sub in subscribers:
            sub.push_terminal(None)
