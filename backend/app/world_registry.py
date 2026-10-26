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
from typing import TYPE_CHECKING, Iterator

from . import public_view, world_names

if TYPE_CHECKING:  # pragma: no cover - typing only, avoids an import cycle
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
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)
    world: World | None = None
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
    # Set by main.py once the session exists: the runner's callbacks need
    # main's step/checkpoint helpers, which cannot be imported here.
    runner: WorldRunner | None = None
    replay_log: ReplayLog | None = None
    checkpoint_store: CheckpointStore | None = None
    deleted: bool = field(default=False, repr=False)
    _operation_busy: bool = field(default=False, repr=False)
    _busy_lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

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

    def record(self) -> dict:
        return {
            "world_id": self.world_id,
            "name": self.name,
            "is_public": self.is_public,
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

    def create(self, name: str | None = None, is_public: bool = False) -> WorldSession:
        with self._lock:
            if len(self._sessions) >= self.max_worlds:
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
