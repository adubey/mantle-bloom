"""The shared world's run/stopped status, and the server-side loop that advances it (#375).

Play in the full app used to be a client loop of `POST /world/step` calls, so the world only
advanced while the owner's tab stayed open, and a restarted server had no idea it had been
running. Now Play and Stop set a status on the server (`POST /world/run`, admin only), and a
background thread steps the world for as long as it's running. The status is checkpointed with
the world (checkpoint.py), so a hosted instance that sleeps or restarts picks up where it left
off: still running, or still stopped.
"""

from __future__ import annotations

import logging
import threading
from typing import Callable

from .checkpoint import RunStatus

# Pause between steps. The step holds the world lock throughout, so this is the window in
# which renders, stats and viewers' requests get the world (Python's lock isn't fair: without a
# pause the loop could take it straight back). Matches the full app's old client-side cadence.
RUN_INTERVAL_S = 0.4

_logger = logging.getLogger(__name__)

# step_once(step_years, should_step) runs one step of the shared world and returns False if
# there's nothing to run. It calls should_step() once it holds the world, and skips the step
# if that's False: the owner may have pressed Stop while the loop was waiting for the world.
StepOnce = Callable[[float, Callable[[], bool]], bool]


class WorldRunner:
    def __init__(
        self,
        step_once: StepOnce,
        on_status: Callable[[RunStatus], None] | None = None,
        interval_s: float = RUN_INTERVAL_S,
    ) -> None:
        self._step_once = step_once
        self._on_status = on_status
        self._interval_s = interval_s
        self._cond = threading.Condition()
        # Serializes status changes with their persistence, so the saved status is the last one set.
        self._status_mutex = threading.Lock()
        self._status = RunStatus()
        self._error: str | None = None
        self._closed = False
        self._thread: threading.Thread | None = None

    def status(self) -> RunStatus:
        with self._cond:
            return self._status

    def error(self) -> str | None:
        """Why the loop last stopped on its own, cleared by the next status change."""
        with self._cond:
            return self._error

    def set_status(self, status: RunStatus, persist: bool = True) -> None:
        """Start, stop or resize the run. `persist=False` is for restoring a status that was
        just read from the checkpoint."""
        with self._status_mutex:
            with self._cond:
                self._status = status
                self._error = None
                if status.running and not self._closed:
                    self._ensure_thread()
                self._cond.notify_all()
            if persist:
                self._persist(status)

    def close(self, timeout_s: float = 1.0) -> None:
        """Stop the loop at shutdown without changing the (saved) status, so a restart resumes
        it. Doesn't wait out a step in progress beyond `timeout_s`: the checkpoint written
        before that step is intact either way."""
        with self._cond:
            self._closed = True
            self._cond.notify_all()
            thread = self._thread
        if thread is not None:
            thread.join(timeout_s)

    def _wants_step(self) -> bool:
        with self._cond:
            return self._status.running and not self._closed

    def _ensure_thread(self) -> None:
        if self._thread is None or not self._thread.is_alive():
            self._thread = threading.Thread(target=self._loop, name="world-runner", daemon=True)
            self._thread.start()

    def _loop(self) -> None:
        while True:
            with self._cond:
                self._cond.wait_for(lambda: self._closed or self._status.running)
                if self._closed:
                    return
                step_years = self._status.step_years
            try:
                stepped = self._step_once(step_years, self._wants_step)
            except Exception as exc:  # noqa: BLE001 - a failing step stops the run rather than the thread
                _logger.exception("run: the step failed; stopping")
                self._stop_on_its_own(f"step failed: {type(exc).__name__}: {exc}")
                continue
            if not stepped:
                _logger.warning("run: there's no world to run; stopping")
                self._stop_on_its_own("there's no world to run")
                continue
            with self._cond:
                self._cond.wait_for(lambda: self._closed or not self._status.running, self._interval_s)

    def _stop_on_its_own(self, reason: str) -> None:
        with self._status_mutex:
            with self._cond:
                if not self._status.running:
                    return  # Stopped meanwhile; that status is already saved.
                self._status = RunStatus(running=False, step_years=self._status.step_years)
                self._error = reason
                status = self._status
            # Saved as stopped: a step that fails here would most likely fail again after a
            # restart, and a restart loop of failures is worse than waiting for the owner.
            self._persist(status)

    def _persist(self, status: RunStatus) -> None:
        if self._on_status is None:
            return
        try:
            self._on_status(status)
        except Exception:  # noqa: BLE001 - the in-memory status still holds; a write can be retried by the next change
            _logger.exception("run: couldn't save the run status")
