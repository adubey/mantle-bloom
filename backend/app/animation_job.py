"""Server-side background animation state for one world (#1033).

Decouples recording execution from the initiating POST /worlds/{world_id}/animate
HTTP stream so refreshing the browser does not stop the animation or lose the
resulting MP4. Both the initiating POST stream and reconnected GET
/worlds/{world_id}/animate/stream clients subscribe to the same NDJSON stream.
"""

from __future__ import annotations

import queue
import threading
import uuid
from typing import Iterator


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
        self._lock = threading.Lock()
        self._active: bool = True
        self._stopping: bool = False
        self._frame: int = 0
        self._elapsed_years: float = elapsed_years
        self._stopped_early: bool = False
        self._error: str | None = None
        self._has_result: bool = False
        # Latest progress NDJSON line so a reconnecting tab immediately gets
        # the current frame PNG, stats, frame counter, and elapsed years.
        self._latest_progress_line: str | None = None
        # Terminal "done" or "error" NDJSON line, kept until dismissed so a tab
        # refreshing after completion still surfaces SaveAnimationModal.
        self._terminal_line: str | None = None
        self._subscribers: list[queue.Queue[str | None]] = []

    def is_active(self) -> bool:
        with self._lock:
            return self._active

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
                "has_result": self._has_result,
                "stopped_early": self._stopped_early,
                "error": self._error,
                "projection": self.projection,
                "view": self.view,
            }

    def mark_stopping(self) -> None:
        with self._lock:
            if self._active:
                self._stopping = True

    def subscribe(self, *, replay_latest: bool = False) -> queue.Queue[str | None]:
        """Register a stream consumer.

        `replay_latest=False` is for the initiating POST /animate stream (registered
        before the worker starts so it receives every frame in order).
        `replay_latest=True` is for a reconnecting GET /animate/stream client: if the
        job is still running it immediately receives the latest progress snapshot,
        and if the job already finished it immediately receives the terminal line.
        """
        sub: queue.Queue[str | None] = queue.Queue()
        with self._lock:
            if self._terminal_line is not None:
                sub.put(self._terminal_line)
                sub.put(None)
                return sub
            if not self._active:
                sub.put(None)
                return sub
            if replay_latest and self._latest_progress_line is not None:
                sub.put(self._latest_progress_line)
            self._subscribers.append(sub)
        return sub

    def unsubscribe(self, sub: queue.Queue[str | None]) -> None:
        with self._lock:
            if sub in self._subscribers:
                self._subscribers.remove(sub)

    def iter_stream(self, sub: queue.Queue[str | None]) -> Iterator[str]:
        """Yield NDJSON lines for one HTTP stream subscriber without affecting the job."""
        try:
            while True:
                line = sub.get()
                if line is None:
                    break
                yield line
        finally:
            self.unsubscribe(sub)

    def publish_progress(self, line: str, *, frame: int, elapsed_years: float) -> None:
        with self._lock:
            self._frame = frame
            self._elapsed_years = elapsed_years
            self._latest_progress_line = line
            for sub in self._subscribers:
                sub.put(line)

    def publish_heartbeat(self, line: str, *, elapsed_years: float | None = None) -> None:
        with self._lock:
            if elapsed_years is not None:
                self._elapsed_years = elapsed_years
            for sub in self._subscribers:
                sub.put(line)

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
            sub.put(line)
            sub.put(None)

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
            sub.put(line)
            sub.put(None)

    def finish_empty(self) -> None:
        """Close out a job that bailed before producing a result (e.g. deleted world)."""
        with self._lock:
            self._active = False
            self._stopping = False
            self._latest_progress_line = None
            subscribers = list(self._subscribers)
            self._subscribers.clear()
        for sub in subscribers:
            sub.put(None)
