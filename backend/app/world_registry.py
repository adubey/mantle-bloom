"""Several independently running worlds in one process (#563).

Until now the server held exactly one world in `main._state`, addressed by
unqualified `/world/*` routes. This module is the registry that replaces that
singleton: each world gets a `WorldSession` carrying everything that used to be
a process-wide global for the one world -- its `World`, state lock, revision
feed, public render cache, background runner, animation stop event, guest-step
coordinator, replay log and checkpoint store.

World identity and management metadata (name, visibility, timestamps) live
here, outside the pickled `World`, so saves do not carry them and the save
format does not change.

Heavy simulation work still cannot run concurrently in one process: Numba's
`parallel=True` kernels are not safe to enter from two threads at once (see
main.py's lock comment). `FairComputeGate` therefore serializes compute across
worlds, first-come-first-served, so independently running worlds progress over
time without starving each other.
"""

from __future__ import annotations

import threading
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Callable, Iterator

from . import public_view, world_names

if TYPE_CHECKING:  # pragma: no cover - typing only, avoids an import cycle
    from .animation_job import AnimationJob
    from .checkpoint import CheckpointStore
    from .replay_log import ReplayLog
    from .run_control import WorldRunner
    from .world import World


class UnknownWorldError(KeyError):
    """No record with that world ID."""


class WorldLimitReached(RuntimeError):
    """Creating another record would exceed the configured world cap."""


class FairComputeGate:
    """FIFO serializer for simulation compute across all worlds.

    Tickets are taken in arrival order and served in that order, so a world
    whose runner asks for steps continuously cannot starve another world's
    runner or a manual step. Deliberately not reentrant: a holder must finish
    before taking another ticket.
    """

    def __init__(self) -> None:
        self._cond = threading.Condition()
        self._next_ticket = 0
        self._serving = 0
        self._cancelled: set[int] = set()

    def acquire(self, timeout: float | None = None) -> bool:
        """Take the next FIFO ticket. With `timeout`, give up (and cancel the
        ticket) after that many seconds; `timeout=0` is a non-blocking try."""
        deadline = None if timeout is None else time.monotonic() + timeout
        with self._cond:
            ticket = self._next_ticket
            self._next_ticket += 1
            while ticket != self._serving:
                if deadline is not None:
                    remaining = deadline - time.monotonic()
                    if remaining <= 0:
                        self._cancelled.add(ticket)
                        return False
                    self._cond.wait(remaining)
                else:
                    self._cond.wait()
            return True

    def release(self) -> None:
        with self._cond:
            self._serving += 1
            while self._serving in self._cancelled:
                self._cancelled.discard(self._serving)
                self._serving += 1
            self._cond.notify_all()

    @contextmanager
    def hold(self) -> Iterator[None]:
        self.acquire()
        try:
            yield
        finally:
            self.release()

@dataclass
class WorldSession:
    """One world's server-side state. `record` fields are the registry's
    metadata; everything else is the per-world replacement for a former
    process-wide singleton in main.py.

    `deleted` marks a session removed from the registry: a runner thread
    parked in the compute gate at delete time stays parked there, and must
    find nothing to do (and nowhere to checkpoint) when it wakes.

    `lock` is a plain Lock, never reentrant: see its comment below."""

    world_id: str
    name: str = "Untitled world"
    is_public: bool = False
    # The one world /view opens on (#755): at most one session carries this at
    # a time (WorldRegistry.set_default enforces it under the registry lock).
    is_default: bool = False
    # The shared /view presentation the owner decides for every viewer of this
    # world (#755): the map view, and the orientation as 9 row-major floats
    # (None is the identity). Guests cannot override either per request.
    viewer_view: str = public_view.GUEST_DEFAULT_VIEW
    viewer_rotation: list[float] | None = None
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)
    world: World | None = None
    deferred_summary: dict | None = None
    checkpoint_dirty: bool = False
    # Why the checkpointed world could not be restored at startup, if it could
    # not. Shown on the record rather than silently starting empty.
    load_error: str | None = None
    # A plain Lock, like the old process-wide `_world_lock`: a client
    # disconnect can finalize a streaming generator on a different thread
    # than the one that acquired it, and only Lock tolerates that release.
    # Never acquire it twice on one thread.
    lock: threading.Lock = field(default_factory=threading.Lock)
    feed: public_view.WorldFeed = field(default_factory=public_view.WorldFeed)
    public_renders: public_view.RenderCache = field(default_factory=public_view.RenderCache)
    guest_steps: public_view.GuestStepCoordinator = field(default_factory=public_view.GuestStepCoordinator)
    animation_stop_event: threading.Event = field(default_factory=threading.Event)
    # Active or finished-uncollected background animation job (#1033).
    animation_job: AnimationJob | None = None
    # Control edits submitted while an animation is running, applied at the next
    # frame boundary once the current frame finishes rendering (#1033).
    _accepting_pending_controls: bool = field(default=False, repr=False)
    _pending_controls: dict[str, object] = field(default_factory=dict, repr=False)
    _controls_snapshot: dict | None = field(default=None, repr=False)
    # Set by main.py once the session exists: the runner's callbacks need
    # main's step/checkpoint helpers, which cannot be imported here.
    runner: WorldRunner | None = None
    replay_log: ReplayLog | None = None
    checkpoint_store: CheckpointStore | None = None
    deleted: bool = field(default=False, repr=False)
    _operation_busy: bool = field(default=False, repr=False)
    _busy_lock: threading.Lock = field(default_factory=threading.Lock, repr=False)
    # A guest step's live progress for /view watchers: (fraction 0-1, phase label),
    # or None when no guest step is computing. Written by the computing thread while
    # it holds the session lock, so it gets its own mutex -- readers (the public
    # world poll) must never park on the session lock mid-step (#653).
    _guest_step_progress: tuple[float, str] | None = field(default=None, repr=False)
    _guest_step_progress_lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    # A counter bumped on every change to what /view shows while a step computes (a span's
    # label, a guest step's line, either clearing). Viewers send back the one they last saw, so a
    # change that lands between two long polls is still answered instead of becoming the new
    # baseline (#801).
    _progress_seq: int = field(default=0, repr=False)

    def bump_progress(self) -> int:
        with self._guest_step_progress_lock:
            self._progress_seq += 1
            return self._progress_seq

    def progress_seq(self) -> int:
        with self._guest_step_progress_lock:
            return self._progress_seq

    def set_guest_step_progress(self, progress: tuple[float, str] | None) -> None:
        """Publish a guest step's current position (or clear it when done). A new
        position wakes the world's viewers so a long poll answers for progress, not
        just for a new revision (see main.py's public_world). Clearing does not wake:
        the step's publish lands microseconds later and wakes for the new revision --
        and waking here would pop the waiter list just before `_world_changed`
        captures its pre-render audience (#658)."""
        with self._guest_step_progress_lock:
            self._guest_step_progress = progress
            self._progress_seq += 1
        if progress is not None:
            self.feed.wake()

    def guest_step_progress(self) -> tuple[float, str] | None:
        """The latest published guest-step position, for the public world poll."""
        with self._guest_step_progress_lock:
            return self._guest_step_progress

    def touch(self) -> None:
        self.updated_at = time.time()

    def begin_operation(self) -> bool:
        """Mark a long manual operation (generate/step/animate/load) as running
        on this world. False when one is already running -- the caller answers
        503, preserving the old same-world conflict behavior."""
        with self._busy_lock:
            if self._operation_busy:
                return False
            self._operation_busy = True
            return True

    def end_operation(self) -> None:
        with self._busy_lock:
            self._operation_busy = False

    def operation_busy(self) -> bool:
        with self._busy_lock:
            return self._operation_busy

    def get_animation_job(self) -> AnimationJob | None:
        with self._busy_lock:
            return self.animation_job

    def has_pending_animation_result(self) -> bool:
        """True while a completed animation's MP4 is still waiting to be saved or discarded."""
        with self._busy_lock:
            job = self.animation_job
            return job is not None and not job.is_active() and job.has_result()

    def begin_animation_operation(self, job: AnimationJob) -> str | None:
        """Atomically reserve the world and publish `job` under `_busy_lock` (#1033).

        Returns None on success, `"busy"` if another operation is already running,
        `"pending_result"` if a finished MP4 is still awaiting Save or Discard, or
        `"no_world"` if no world has been generated yet.
        """
        with self._busy_lock:
            if self._operation_busy:
                return "busy"
            existing = self.animation_job
            if existing is not None and not existing.is_active() and existing.has_result():
                return "pending_result"
            if self.world is None and self.deferred_summary is None:
                return "no_world"
            self.animation_stop_event.clear()
            self._operation_busy = True
            self._accepting_pending_controls = True
            self.animation_job = job
            return None

    def finish_animation_operation(self, complete_job: Callable[[], None]) -> None:
        """Publish the job's terminal state and release `_operation_busy` atomically (#1033).

        Prevents a gap where `_operation_busy` is False before `job.has_result()` is True,
        which could let a second POST /animate clobber the finished MP4 or let
        `_session_can_sleep` retire the worker.
        """
        with self._busy_lock:
            self._accepting_pending_controls = False
            complete_job()
            self._operation_busy = False

    def clear_finished_animation_job(self, job_id: str | None = None) -> None:
        """Atomically drop a completed animation job (and its buffered MP4) once saved,
        discarded, or replaced by a fresh world. Never clears an active job or races
        with a newly published one (#1033)."""
        with self._busy_lock:
            job = self.animation_job
            if job is None or job.is_active():
                return
            if job_id is not None and job.job_id != job_id:
                return
            self.animation_job = None

    def try_queue_pending_controls(self, updates: dict[str, object]) -> tuple[bool, dict[str, object]]:
        """Atomically enqueue control changes only while an animation is accepting them (#1033).

        Returns `(True, merged_pending)` if queued, or `(False, {})` once the worker has
        closed the queue for finalization so the caller applies them under the normal lock.
        """
        with self._busy_lock:
            if not self._accepting_pending_controls:
                return False, {}
            self._pending_controls.update(updates)
            return True, dict(self._pending_controls)

    def get_pending_controls(self) -> dict[str, object]:
        with self._busy_lock:
            return dict(self._pending_controls)

    def drain_pending_controls(self, *, close_queue: bool = False) -> dict[str, object]:
        with self._busy_lock:
            if close_queue:
                self._accepting_pending_controls = False
            pending = dict(self._pending_controls)
            self._pending_controls.clear()
            return pending

    def clear_pending_controls(self) -> None:
        with self._busy_lock:
            self._pending_controls.clear()

    def set_controls_snapshot(self, state: dict) -> None:
        with self._busy_lock:
            self._controls_snapshot = dict(state)

    def get_controls_snapshot(self) -> dict | None:
        with self._busy_lock:
            return None if self._controls_snapshot is None else dict(self._controls_snapshot)

    def record(self) -> dict:
        return {
            "world_id": self.world_id,
            "name": self.name,
            "is_public": self.is_public,
            "is_default": self.is_default,
            "viewer_view": self.viewer_view,
            "viewer_rotation": self.viewer_rotation,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
        }


class WorldRegistry:
    """All world sessions in this process, keyed by world ID."""

    def __init__(self, max_worlds: int) -> None:
        if max_worlds < 1:
            raise ValueError("max_worlds must be at least 1")
        self.max_worlds = max_worlds
        self._lock = threading.Lock()
        self._sessions: dict[str, WorldSession] = {}

    def __len__(self) -> int:
        with self._lock:
            return len(self._sessions)

    def create(self, name: str | None = None, is_public: bool = False, *, enforce_limit: bool = True) -> WorldSession:
        with self._lock:
            if enforce_limit and len(self._sessions) >= self.max_worlds:
                raise WorldLimitReached(
                    f"world limit reached ({self.max_worlds}); delete a world first"
                )
            session = WorldSession(
                world_id=uuid.uuid4().hex,
                name=name or world_names.default_world_name({s.name for s in self._sessions.values()}),
                is_public=is_public,
            )
            self._sessions[session.world_id] = session
            return session

    def adopt(self, session: WorldSession) -> WorldSession:
        """Register a session built elsewhere (checkpoint restore). Existing
        checkpointed worlds are always restored, even past the cap: the cap
        blocks creating new records, never recovering saved ones."""
        with self._lock:
            if session.world_id in self._sessions:
                raise ValueError(f"duplicate world_id {session.world_id!r}")
            self._sessions[session.world_id] = session
            return session

    def set_default(self, world_id: str | None) -> list[WorldSession]:
        """Make `world_id` the one default world (None: no default, #755).

        The check and the clearing of every other session happen under the
        registry lock, so two owner tabs setting different worlds at once can
        never leave two defaults. Returns the sessions whose flag changed."""
        with self._lock:
            sessions = list(self._sessions.values())
            if world_id is not None and world_id not in self._sessions:
                raise UnknownWorldError(world_id)
            changed = []
            for session in sessions:
                wanted = session.world_id == world_id
                if session.is_default != wanted:
                    session.is_default = wanted
                    changed.append(session)
            return changed

    def default(self) -> WorldSession | None:
        """The default world, if one is set (#755)."""
        with self._lock:
            for session in self._sessions.values():
                if session.is_default:
                    return session
            return None

    def get(self, world_id: str) -> WorldSession:
        with self._lock:
            session = self._sessions.get(world_id)
        if session is None:
            raise UnknownWorldError(world_id)
        return session

    def find(self, world_id: str) -> WorldSession | None:
        with self._lock:
            return self._sessions.get(world_id)

    def list(self) -> list[WorldSession]:
        with self._lock:
            return sorted(self._sessions.values(), key=lambda s: s.created_at)

    def remove(self, world_id: str) -> WorldSession:
        with self._lock:
            session = self._sessions.pop(world_id, None)
        if session is None:
            raise UnknownWorldError(world_id)
        return session
