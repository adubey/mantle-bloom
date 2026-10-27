"""The hosted world's checkpoint on persistent storage (#375).

A hosted instance (Railway) sleeps when idle and restarts on every deploy, and both throw away
the one in-memory world. So after every completed change to the world, main.py writes it here,
along with the owner's run/stopped status (run_control.py), and loads both back at startup:
the world resumes from exactly the saved point, without simulating the time it spent asleep,
and keeps running only if it was running before.

Nothing depends on a hook before sleep or shutdown, which may never run. Each file is instead
replaced atomically (written to a temporary file in the same directory, flushed to disk, then
renamed over the old one), so a process killed partway through a write leaves the previous
checkpoint intact. The world and the status are separate files because the status changes on
its own (Play/Stop), and rewriting a 50+ MB world to flip one flag would be wasteful; each file
is self-consistent, and any pairing of the two is a valid state to resume.

The world file is an ordinary `.mbworld` save (persistence.py), so it follows
docs/save-compatibility.md and can be downloaded and opened in the full app or bin/debug.

Checkpointing is off unless `MANTLE_BLOOM_CHECKPOINT_DIR` (or, failing that, Railway's
`RAILWAY_VOLUME_MOUNT_PATH`) names a directory. A local run keeps no checkpoint.
"""

from __future__ import annotations

import json
import logging
import os
import re
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path

from . import persistence
from .world import World

CHECKPOINT_DIR_ENV = "MANTLE_BLOOM_CHECKPOINT_DIR"
# Railway sets this on a service with a volume attached, so the volume works with no other
# configuration.
RAILWAY_VOLUME_ENV = "RAILWAY_VOLUME_MOUNT_PATH"

WORLD_FILE = "world.mbworld"
STATUS_FILE = "run_status.json"
STATUS_VERSION = 1
# The multi-world registry index (#563): metadata for every world whose
# checkpoint lives under `WORLDS_DIR/<world_id>/` in the same root.
WORLDS_INDEX_FILE = "worlds.json"
WORLDS_DIR = "worlds"
INDEX_VERSION = 1


class RegistryIndexError(ValueError):
    """The registry index exists but cannot be trusted. Startup must fail its
    health check rather than silently start with an empty registry."""
# The full app's default "Years per step" (App.tsx), for a status written before any Play.
DEFAULT_STEP_YEARS = 100_000.0
# Prefix of in-progress writes. Any left behind were cut off by a killed process.
_TEMP_PREFIX = ".checkpoint-"

_logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RunStatus:
    running: bool = False
    step_years: float = DEFAULT_STEP_YEARS


def valid_step_years(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and 0 < value < float("inf")


class CheckpointStore:
    def __init__(self, directory: str | os.PathLike[str]) -> None:
        self.directory = Path(directory)
        # Writers each hold the world lock already; this keeps a status write (which doesn't)
        # from interleaving with a world write in the same directory listing.
        self._mutex = threading.Lock()

    @classmethod
    def from_env(cls) -> CheckpointStore | None:
        directory = os.environ.get(CHECKPOINT_DIR_ENV, "").strip() or os.environ.get(RAILWAY_VOLUME_ENV, "").strip()
        return cls(directory) if directory else None

    # Registry IDs are uuid4 hex; anything else reaching a per-world path would
    # be a tampered index turning rmtree/replay-log setup into a path traversal.
    _WORLD_ID_PATTERN = re.compile(r"[0-9a-f]{32}")

    @staticmethod
    def _is_timestamp(value: object) -> bool:
        return isinstance(value, (int, float)) and not isinstance(value, bool)

    def for_world(self, world_id: str) -> CheckpointStore:
        """The store for one world's checkpoint subdirectory (#563)."""
        if not self._WORLD_ID_PATTERN.fullmatch(world_id):
            raise RegistryIndexError(f"refusing suspicious world_id {world_id!r}")
        return CheckpointStore(self.directory / WORLDS_DIR / world_id)

    @property
    def index_path(self) -> Path:
        return self.directory / WORLDS_INDEX_FILE

    def save_index(self, records: list[dict]) -> None:
        """Atomically replace the registry index. Called on record create,
        rename, publish/unpublish, import and delete -- not on every step."""
        data = json.dumps({"version": INDEX_VERSION, "worlds": records})
        self._replace(self.index_path, data.encode())

    def load_index(self) -> list[dict] | None:
        """The registry index records, or None when no index exists yet (a
        pre-multi-world checkpoint directory). Raises RegistryIndexError when
        an index exists but is unreadable or malformed."""
        path = self.index_path
        if not path.exists():
            return None
        try:
            payload = json.loads(path.read_text())
            if payload.get("version") != INDEX_VERSION or not isinstance(payload.get("worlds"), list):
                raise ValueError(f"unexpected contents {payload!r}")
            for record in payload["worlds"]:
                if (
                    not isinstance(record, dict)
                    or not isinstance(record.get("world_id"), str)
                    or not self._WORLD_ID_PATTERN.fullmatch(record["world_id"])
                    # Timestamps are read with float() at restore: a non-numeric
                    # one would raise outside any try and fail startup for every
                    # world, so reject it here where the index fails closed.
                    # (bool is an int subclass; a True timestamp is nonsense too.)
                    or not isinstance(record.get("name"), str)
                    or not isinstance(record.get("is_public"), bool)
                    or not self._is_timestamp(record.get("created_at"))
                    or not self._is_timestamp(record.get("updated_at"))
                ):
                    raise ValueError(f"unexpected record {record!r}")
        except Exception as exc:  # noqa: BLE001 - any bad index must fail closed
            raise RegistryIndexError(f"can't read registry index {path}: {exc}") from exc
        return payload["worlds"]

    @property
    def world_path(self) -> Path:
        return self.directory / WORLD_FILE

    @property
    def status_path(self) -> Path:
        return self.directory / STATUS_FILE

    def save_world(self, world: World) -> None:
        """Replace the world checkpoint. The caller holds the world lock, so the world can't
        change while it's serialized."""
        started = time.monotonic()
        data = persistence.save_world_bytes(world)
        self._replace(self.world_path, data)
        _logger.debug(
            "checkpointed world at %.0f years (%.1f MB, %.2f s)",
            world.elapsed_years, len(data) / 1e6, time.monotonic() - started,
        )

    def save_status(self, status: RunStatus) -> None:
        data = json.dumps({"version": STATUS_VERSION, "running": status.running, "step_years": status.step_years})
        self._replace(self.status_path, data.encode())

    def load(self) -> tuple[World | None, RunStatus]:
        """The checkpointed world (None if there isn't a usable one) and run status. Never
        raises for a missing or unreadable checkpoint: the server must still start, so it logs
        why and starts without that part. An unreadable world file is renamed aside rather
        than left to be overwritten by the next checkpoint, so it can still be examined.

        Raises OSError if the directory can't be created or written: a volume that can't hold
        a checkpoint must fail the deploy's health check (#376), not the first step's write."""
        self.directory.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryFile(prefix=_TEMP_PREFIX, dir=self.directory):
            pass
        self._remove_partial_writes()
        world = self._load_world()
        status = self._load_status()
        if status.running and world is None:
            _logger.warning("checkpoint: the saved status is running, but there's no world to run; staying stopped")
            status = RunStatus(running=False, step_years=status.step_years)
            # Saved too, or the next world the owner makes would start running at the next restart.
            self.save_status(status)
        return world, status

    def _load_world(self) -> World | None:
        path = self.world_path
        if not path.exists():
            # A warning, so a volume mounted at the wrong path shows in the deploy log.
            _logger.warning("checkpoint: no saved world in %s; starting without one", self.directory)
            return None
        try:
            # Stream from disk: a checkpoint is a whole world, and holding its bytes
            # beside the world it builds doubles restore-time RAM (#648).
            with open(path, "rb") as checkpoint_file:
                world = persistence.load_world_file(checkpoint_file)
        except Exception as exc:  # noqa: BLE001 - any unreadable file must leave the server startable
            aside = path.with_name(f"{path.name}.invalid-{int(time.time())}")
            _logger.error(
                "checkpoint: can't load %s (%s: %s); starting without a world and moving it to %s",
                path, type(exc).__name__, exc, aside.name,
            )
            try:
                path.replace(aside)
            except OSError:
                _logger.exception("checkpoint: couldn't move the unreadable world aside")
            return None
        _logger.info("checkpoint: restored world seed %s at %.0f years from %s", world.seed, world.elapsed_years, path)
        return world

    def _load_status(self) -> RunStatus:
        path = self.status_path
        if not path.exists():
            return RunStatus()
        try:
            payload = json.loads(path.read_text())
            running, step_years = payload["running"], payload["step_years"]
            if payload.get("version") != STATUS_VERSION or not isinstance(running, bool) or not valid_step_years(step_years):
                raise ValueError(f"unexpected contents {payload!r}")
        except Exception as exc:  # noqa: BLE001 - same as _load_world
            _logger.error("checkpoint: can't read %s (%s: %s); staying stopped", path, type(exc).__name__, exc)
            return RunStatus()
        _logger.info("checkpoint: restored run status %s (%.0f years per step)", "running" if running else "stopped", step_years)
        return RunStatus(running=running, step_years=float(step_years))

    def _replace(self, path: Path, data: bytes) -> None:
        with self._mutex:
            self.directory.mkdir(parents=True, exist_ok=True)
            fd, temp_name = tempfile.mkstemp(prefix=_TEMP_PREFIX, dir=self.directory)
            try:
                with os.fdopen(fd, "wb") as handle:
                    handle.write(data)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temp_name, path)
            except BaseException:
                Path(temp_name).unlink(missing_ok=True)
                raise
            _fsync_directory(self.directory)

    def _remove_partial_writes(self) -> None:
        for leftover in self.directory.glob(f"{_TEMP_PREFIX}*"):
            _logger.warning("checkpoint: removing %s, a write cut off by a stopped process", leftover.name)
            leftover.unlink(missing_ok=True)


def _fsync_directory(directory: Path) -> None:
    """Make the rename itself durable: until the directory entry reaches the disk, a power
    loss can still bring back the old file."""
    try:
        fd = os.open(directory, os.O_RDONLY)
    except OSError:
        return  # Not supported on this platform/filesystem; the rename is still atomic.
    try:
        os.fsync(fd)
    except OSError:
        pass
    finally:
        os.close(fd)
