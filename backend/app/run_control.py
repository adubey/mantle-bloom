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


class StatusNotSaved(RuntimeError):
    """The status changed in memory but couldn't be written to the checkpoint. It's retried
    with the next world checkpoint (`retry_unsaved_status`)."""


class WorldRunner:
    def __init__(
        self,
        step_once: StepOnce,
        on_status: Callable[[RunStatus], None] | None = None,
        on_change: Callable[[], None] | None = None,
        interval_s: float = RUN_INTERVAL_S,
    ) -> None:
        self._step_once = step_once
        # Saves a status to the checkpoint.
        self._on_status = on_status
        # Called after `version` moves, to wake anyone long-polling the run state.
        self._on_change = on_change
        self._interval_s = interval_s
        self._cond = threading.Condition()
        # Serializes status changes with their persistence, so the saved status is the last one set.
        self._status_mutex = threading.Lock()
        self._status = RunStatus()
        self._error: str | None = None
        self._stepping = False
        self._unsaved = False
        # Bumped on every change to the status or error, and when a step finishes after Stop,
        # so a long poll can tell that something other than the world changed (Stop pressed in
        # another tab, or the last step landing). Not on every step of a run: those publish a
        # new world revision anyway, and a wake per step would also wake every viewer's poll.
        self._version = 0
        # Bumped by close(). A loop exits once its generation is out of date, so a later
        # start (the app's lifespan can run more than once in a process, as tests do) gets a
        # fresh loop rather than a runner that's closed for good.
        self._generation = 0
        self._thread: threading.Thread | None = None

    def status(self) -> RunStatus:
        with self._cond:
            return self._status

    def error(self) -> str | None:
        """Why the loop last stopped on its own, cleared by the next status change."""
        with self._cond:
            return self._error

    def stepping(self) -> bool:
        """Whether the loop is taking a step. It can still be after Stop: a step under way
        finishes, and the world moves once more."""
        with self._cond:
            return self._stepping

    def version(self) -> int:
        with self._cond:
            return self._version

    def set_status(self, status: RunStatus, persist: bool = True) -> None:
        """Start, stop or resize the run. `persist=False` is for restoring a status that was
        just read from the checkpoint. Raises StatusNotSaved if the status took effect but
        couldn't be saved, so the owner hears that a restart may not honor it."""
        with self._status_mutex:
            with self._cond:
                self._status = status
                self._error = None
                if status.running:
                    self._ensure_thread()
                self._bump()
            self._changed()
            if persist and not self._persist(status):
                raise StatusNotSaved("the run status couldn't be saved; a restart may not keep it")

    def retry_unsaved_status(self) -> None:
        """Save the status again if its last save failed. Called with each world checkpoint,
        so a volume that comes back catches the status up too."""
        with self._status_mutex:
            with self._cond:
                if not self._unsaved:
                    return
                status = self._status
            self._persist(status)

    def close(self, timeout_s: float = 1.0) -> None:
        """Stop the loop at shutdown without changing the (saved) status, so a restart resumes
        it. Doesn't wait out a step in progress beyond `timeout_s`: the checkpoint written
        before that step is intact either way."""
        with self._cond:
            self._generation += 1
            self._cond.notify_all()
            thread, self._thread = self._thread, None
        if thread is not None:
            thread.join(timeout_s)

    def _bump(self) -> None:
        # The caller holds `_cond`.
        self._version += 1
        self._cond.notify_all()

    def _changed(self) -> None:
        if self._on_change is not None:
            self._on_change()

    def _ensure_thread(self) -> None:
        if self._thread is None or not self._thread.is_alive():
            self._thread = threading.Thread(target=self._loop, args=(self._generation,), name="world-runner", daemon=True)
            self._thread.start()

    def _loop(self, generation: int) -> None:
        def closed() -> bool:
            return self._generation != generation

        def wants_step() -> bool:
            with self._cond:
                return self._status.running and not closed()

        while True:
            with self._cond:
                self._cond.wait_for(lambda: closed() or self._status.running)
                if closed():
                    return
                status = self._status
                self._stepping = True
            try:
                stepped = self._step_once(status.step_years, wants_step)
            except Exception as exc:  # noqa: BLE001 - a failing step stops the run rather than the thread
                _logger.exception("run: the step failed; stopping")
                self._stop_on_its_own(status, f"step failed: {type(exc).__name__}: {exc}")
                continue
            finally:
                with self._cond:
                    self._stepping = False
                    stopped = not self._status.running
                    if stopped:
                        self._bump()
                if stopped:
                    self._changed()
            if not stepped:
                _logger.warning("run: there's no world to run; stopping")
                self._stop_on_its_own(status, "there's no world to run")
                continue
            with self._cond:
                self._cond.wait_for(lambda: closed() or not self._status.running, self._interval_s)

    def _stop_on_its_own(self, stepped: RunStatus, reason: str) -> None:
        with self._status_mutex:
            with self._cond:
                # Only the status that step ran under: if the owner pressed Stop, or Play
                # again with a new step size, meanwhile, theirs is newer and stands.
                if self._status is not stepped:
                    return
                self._status = RunStatus(running=False, step_years=stepped.step_years)
                self._error = reason
                status = self._status
                self._bump()
            self._changed()
            # Saved as stopped: a step that fails here would most likely fail again after a
            # restart, and a restart loop of failures is worse than waiting for the owner.
            self._persist(status)

    def _persist(self, status: RunStatus) -> bool:
        """Save `status`; the caller holds `_status_mutex`. A failure is logged and marked for
        `retry_unsaved_status`."""
        if self._on_status is None:
            return True
        try:
            self._on_status(status)
        except Exception:  # noqa: BLE001 - the in-memory status holds; retried with the next world checkpoint
            _logger.exception("run: couldn't save the run status")
            saved = False
        else:
            saved = True
        with self._cond:
            self._unsaved = not saved
        return saved
