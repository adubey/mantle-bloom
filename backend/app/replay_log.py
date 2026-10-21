"""Replay logs: save the world before a step whose health goes askew, for replay later (#442).

The RFC diagnosis gates hunt for rare events (a merge that clips craton, a ledger that stops
closing) by replaying a seed for hundreds of Myr, which takes hours and may not reproduce the
event under a newer build. The same events turn up in ordinary runs. While the server advances
the shared world (main.py's /world/step, /world/animate, the run loop and guest steps), this
module measures a few health metrics around each step and, when a trigger fires, writes a
record that a debug script or an RFC analysis can replay:

- `before.mbworld`: the world just before the step, an ordinary save (persistence.py) that
  loads in the app and in bin/debug.
- `record.json`: which triggers fired and their values, the step length, the world's scalar
  controls, the git commit, the health metrics before and after, the step's wall time, and a
  fingerprint of the world after the step (`fingerprint`).

Every random stream a step draws is seeded from `(world.seed, elapsed_years, ...)`, so the
pre-step world and the step length reproduce the step exactly (`replay` checks this against
the fingerprint, so a record says when it doesn't). The fingerprint stands in for a second
60 MB pickle of the post-step world.

Keeping the pre-step world means pickling it before every step, then dropping the bytes unless
a trigger fires: about 45 ms against a 1.8 s step at node density 1 (2.6%, #442's measurement).

Triggers are named functions of a `StepOutcome` in `TRIGGERS`. An RFC adds its own signature
with `TRIGGERS["name"] = fn`, without touching the run loop. Code inside a step can attach inputs
that only exist there with `annotate(key, value)`; they ride along in the record, though the
pre-step world must still be enough to regenerate them.

Logging is off unless `MANTLE_BLOOM_REPLAY_LOG_DIR` names a directory (`bin/restart.sh
--replay-log`). Records are rate-limited per trigger, capped by total size (oldest removed
first), and written by a background thread so viewers don't wait on a 60 MB write.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import queue
import shutil
import subprocess
import tempfile
import threading
import time
import traceback
from collections import deque
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterator

import numpy as np

from . import continental_ledger, cratons, persistence
from .elevation_lines import line_spacing_rad
from .world import World, step_world_progress

REPLAY_LOG_DIR_ENV = "MANTLE_BLOOM_REPLAY_LOG_DIR"
RECORD_VERSION = 1
BEFORE_FILE = "before.mbworld"
RECORD_FILE = "record.json"
# Total size of all records before the oldest go. One record is about one save (50-70 MB at
# node density 1), so this keeps a few dozen.
DEFAULT_MAX_BYTES = 2_000_000_000
# A trigger that has fired stays quiet for this many logged steps: a world that stays askew
# would otherwise write a record every step.
DEFAULT_MIN_STEPS_BETWEEN = 20
# Messages kept for the ?deb Event Console (main.py's summary).
MAX_RECENT_EVENTS = 50
_TEMP_PREFIX = ".replay-log-"

# Starting thresholds. Loose on purpose (#442): tighten them from what fires in practice.
# A single sink account destroying more than this share of the step's starting live craton.
# RFC 0001's merge-clip signature is 1%.
CRATON_STEP_LOSS_FRACTION = 0.01
# Change in cratons.balance_error_m3 over one step, relative to live craton volume. The change,
# not the level: an old save whose ledger never closed would otherwise fire on every step.
CRATON_BALANCE_TOLERANCE = 1e-9
# Relative change of genuine continental area in one step.
CONTINENTAL_AREA_JUMP_FRACTION = 0.05
# Absolute change of the land fraction in one step.
LAND_FRACTION_JUMP = 0.05
# Relative change of the plate count in one step.
PLATE_COUNT_JUMP_FRACTION = 0.25

_logger = logging.getLogger(__name__)


# --- Health metrics ----------------------------------------------------------------------------


@dataclass(frozen=True)
class Health:
    plate_count: int
    craton_live_m3: float
    craton_ledger: dict[str, float]
    craton_balance_error_m3: float
    # Why continental_ledger.assert_closed failed, or None when the ledger closes.
    continental_ledger_error: str | None
    genuine_continental_area_m2: float
    # Area-weighted share of nodes above sea level. Not stats.py's lake-aware land fraction,
    # which needs the climate grid; this one costs one pass over the plates.
    land_fraction: float


def measure(world: World) -> Health:
    cratons.ensure_initialized(world)
    spacing_rad = line_spacing_rad(world.node_density)
    land_m2 = total_m2 = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        _, elevation = plate.all_points_and_elevation()
        area_m2 = plate.accounting_areas_m2(spacing_rad)
        land_m2 += float(area_m2[elevation > world.sea_level_m].sum())
        total_m2 += float(area_m2.sum())
    try:
        continental_ledger.assert_closed(world)
        ledger_error = None
    except AssertionError as exc:
        ledger_error = str(exc)
    return Health(
        plate_count=len(world.plates),
        craton_live_m3=cratons.live_volume_m3(world),
        craton_ledger={key: float(value) for key, value in world.craton_ledger.items()},
        craton_balance_error_m3=cratons.balance_error_m3(world),
        continental_ledger_error=ledger_error,
        genuine_continental_area_m2=cratons.genuine_continental_area_m2(world),
        land_fraction=land_m2 / total_m2 if total_m2 > 0.0 else 0.0,
    )


# --- Triggers ----------------------------------------------------------------------------------


@dataclass
class StepOutcome:
    years: float
    before: Health
    # None when the step raised.
    after: Health | None
    error: BaseException | None = None
    # Inputs attached from inside the step with `annotate`.
    annotations: dict = field(default_factory=dict)


# A trigger returns the values that make it fire (JSON-serializable, for the record), or None.
Trigger = Callable[[StepOutcome], "dict | None"]


def _relative_change(before: float, after: float) -> float:
    if before == 0.0:
        return 0.0 if after == 0.0 else float("inf")
    return abs(after - before) / abs(before)


def _craton_step_loss(outcome: StepOutcome) -> dict | None:
    if outcome.after is None or outcome.before.craton_live_m3 <= 0.0:
        return None
    losses = {
        key: (outcome.after.craton_ledger.get(key, 0.0) - outcome.before.craton_ledger.get(key, 0.0))
        / outcome.before.craton_live_m3
        for key in cratons.CRATON_SINK_ACCOUNTS
    }
    over = {key: share for key, share in losses.items() if share > CRATON_STEP_LOSS_FRACTION}
    return {"share_of_live_craton": over, "threshold": CRATON_STEP_LOSS_FRACTION} if over else None


def _craton_ledger_balance(outcome: StepOutcome) -> dict | None:
    if outcome.after is None:
        return None
    change = outcome.after.craton_balance_error_m3 - outcome.before.craton_balance_error_m3
    live = max(outcome.after.craton_live_m3, outcome.before.craton_live_m3)
    if abs(change) <= CRATON_BALANCE_TOLERANCE * live:
        return None
    return {"balance_error_change_m3": change, "live_m3": live, "threshold": CRATON_BALANCE_TOLERANCE}


def _continental_ledger_open(outcome: StepOutcome) -> dict | None:
    # Fires on the step that opens the ledger, not on every step after it.
    if outcome.after is None or outcome.before.continental_ledger_error is not None:
        return None
    if outcome.after.continental_ledger_error is None:
        return None
    return {"error": outcome.after.continental_ledger_error}


def _jump(name: str, before: float, after: float, change: float, threshold: float) -> dict | None:
    if change <= threshold:
        return None
    return {f"{name}_before": before, f"{name}_after": after, "change": change, "threshold": threshold}


def _continental_area_jump(outcome: StepOutcome) -> dict | None:
    if outcome.after is None:
        return None
    before, after = outcome.before.genuine_continental_area_m2, outcome.after.genuine_continental_area_m2
    return _jump("area_m2", before, after, _relative_change(before, after), CONTINENTAL_AREA_JUMP_FRACTION)


def _land_fraction_jump(outcome: StepOutcome) -> dict | None:
    if outcome.after is None:
        return None
    before, after = outcome.before.land_fraction, outcome.after.land_fraction
    return _jump("land_fraction", before, after, abs(after - before), LAND_FRACTION_JUMP)


def _plate_count_jump(outcome: StepOutcome) -> dict | None:
    if outcome.after is None:
        return None
    before, after = outcome.before.plate_count, outcome.after.plate_count
    return _jump("plates", before, after, _relative_change(before, after), PLATE_COUNT_JUMP_FRACTION)


def _step_exception(outcome: StepOutcome) -> dict | None:
    if outcome.error is None:
        return None
    return {
        "type": type(outcome.error).__name__,
        "message": str(outcome.error),
        "traceback": "".join(traceback.format_exception(outcome.error)),
    }


# Every trigger a logged step is checked against, by name. An RFC adds its own signature with
# `TRIGGERS["name"] = fn`.
TRIGGERS: dict[str, Trigger] = {
    "craton_step_loss": _craton_step_loss,
    "craton_ledger_balance": _craton_ledger_balance,
    "continental_ledger_open": _continental_ledger_open,
    "continental_area_jump": _continental_area_jump,
    "land_fraction_jump": _land_fraction_jump,
    "plate_count_jump": _plate_count_jump,
    "step_exception": _step_exception,
}


# --- Annotations from inside a step ------------------------------------------------------------

_active = threading.local()


def annotate(key: str, value: object) -> None:
    """Attach `value` (JSON-serializable) to the record of the step under way, if it's being
    logged; a no-op otherwise, so simulation code can call it unconditionally."""
    annotations = getattr(_active, "annotations", None)
    if annotations is not None:
        annotations[key] = value


# --- Fingerprint of a world --------------------------------------------------------------------

# Dict keys and attributes that time the step itself; a replay never reproduces them.
_WALL_CLOCK_SUFFIX = "_seconds"


def _feed(h, value: object, seen: set[int]) -> None:
    if value is None or isinstance(value, (bool, int, float, complex, str, bytes, np.generic)):
        h.update(f"{type(value).__name__}:{value!r};".encode())
        return
    if isinstance(value, np.ndarray):
        h.update(f"nd:{value.dtype}:{value.shape};".encode())
        if value.dtype.hasobject:
            for item in value.ravel():
                _feed(h, item, seen)
        else:
            h.update(np.ascontiguousarray(value).tobytes())
        return
    # Only the objects on the current path, to stop at a cycle: which other objects a world
    # shares (an interned empty tuple, say) doesn't survive a save and load.
    if id(value) in seen:
        h.update(b"cycle;")
        return
    seen.add(id(value))
    try:
        _feed_container(h, value, seen)
    finally:
        seen.discard(id(value))


def _feed_container(h, value: object, seen: set[int]) -> None:
    if isinstance(value, dict):
        items = [(k, v) for k, v in value.items() if not (isinstance(k, str) and k.endswith(_WALL_CLOCK_SUFFIX))]
        h.update(f"dict:{len(items)};".encode())
        for key, item in items:
            _feed(h, key, seen)
            _feed(h, item, seen)
    elif isinstance(value, (list, tuple)):
        h.update(f"{type(value).__name__}:{len(value)};".encode())
        for item in value:
            _feed(h, item, seen)
    elif isinstance(value, (set, frozenset)):
        h.update(f"set:{sorted(repr(item) for item in value)!r};".encode())
    else:
        _feed_object(h, value, seen)


def _feed_object(h, value: object, seen: set[int]) -> None:
    h.update(f"obj:{type(value).__module__}.{type(value).__qualname__};".encode())
    if is_dataclass(value):
        names = [f.name for f in fields(value)]
    elif hasattr(value, "__dict__"):
        # Sorted: a loaded plate's __setstate__ sets its attributes in another order.
        names = sorted(vars(value))
    else:
        # A k-d tree or other extension object: a cache rebuilt from the state around it.
        return
    for name in names:
        if name.endswith(_WALL_CLOCK_SUFFIX):
            continue
        h.update(f".{name}=".encode())
        _feed(h, getattr(value, name, None), seen)


def fingerprint(world: World) -> dict[str, str]:
    """A digest per `World` attribute, so a replay that doesn't match says which parts differ.
    Covers everything a save carries except wall-clock timings."""
    result = {}
    for name, value in vars(world).items():
        h = hashlib.sha256()
        _feed(h, value, set())
        result[name] = h.hexdigest()[:16]
    return result


# --- Records -----------------------------------------------------------------------------------


@dataclass(frozen=True)
class _PendingRecord:
    name: str
    before_bytes: bytes
    sidecar: dict


def _controls(world: World) -> dict:
    """The world's scalar fields: the /world/controls settings, flags and densities. The
    pre-step save carries them too; these are for reading the record without loading it."""
    return {
        name: value for name, value in vars(world).items()
        if isinstance(value, (bool, int, float, str)) and not isinstance(value, np.generic)
    }


def _git_state() -> dict | None:
    root = Path(__file__).resolve().parents[2]
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True, timeout=5,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=no"], cwd=root, capture_output=True, text=True,
            check=True, timeout=5,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None  # A packaged build (bin/package.sh) has no checkout.
    return {"commit": commit, "dirty": bool(status.strip())}


def _record_name(seed: int, elapsed_years: float, fired: dict, now: float) -> str:
    stamp = datetime.fromtimestamp(now, timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    return f"{stamp}-seed{seed}-{elapsed_years / 1e6:.2f}myr-{'+'.join(sorted(fired))}"


class ReplayLog:
    def __init__(
        self,
        directory: str | os.PathLike[str],
        *,
        max_bytes: int = DEFAULT_MAX_BYTES,
        min_steps_between: int = DEFAULT_MIN_STEPS_BETWEEN,
        triggers: dict[str, Trigger] | None = None,
        step_progress: Callable[[World, float], Iterator[float]] = step_world_progress,
    ) -> None:
        self.directory = Path(directory)
        self.max_bytes = max_bytes
        self.min_steps_between = min_steps_between
        # TRIGGERS itself by default, so one registered after construction still counts.
        self.triggers = TRIGGERS if triggers is None else triggers
        self._step_progress = step_progress
        self._git = _git_state()
        self._steps = 0
        self._last_fired: dict[str, int] = {}
        self._recent: deque[dict] = deque(maxlen=MAX_RECENT_EVENTS)
        self._mutex = threading.Lock()
        self._queue: queue.Queue[_PendingRecord] = queue.Queue()
        self._writer: threading.Thread | None = None
        self.directory.mkdir(parents=True, exist_ok=True)
        for leftover in self.directory.glob(f"{_TEMP_PREFIX}*"):
            shutil.rmtree(leftover, ignore_errors=True)

    @classmethod
    def from_env(cls) -> ReplayLog | None:
        directory = os.environ.get(REPLAY_LOG_DIR_ENV, "").strip()
        return cls(Path(directory).expanduser()) if directory else None

    def step(self, world: World, years: float) -> Iterator[float]:
        """`step_world_progress(world, years)`, logging a record if a trigger fires. Never
        lets a failure of its own disturb the step."""
        try:
            before_bytes = persistence.save_world_bytes(world)
            before = measure(world)
        except Exception:  # noqa: BLE001 - the step goes ahead unlogged
            _logger.exception("replay log: couldn't snapshot the world; stepping without a log")
            yield from self._step_progress(world, years)
            return
        snapshot = (world.seed, world.elapsed_years, world.steps_taken)
        _active.annotations = annotations = {}
        started = time.monotonic()
        try:
            yield from self._step_progress(world, years)
        except Exception as exc:
            _active.annotations = None
            outcome = StepOutcome(years, before, None, exc, annotations)
            self._after_step(world, snapshot, before_bytes, outcome, time.monotonic() - started)
            raise
        _active.annotations = None
        seconds = time.monotonic() - started
        try:
            outcome = StepOutcome(years, before, measure(world), None, annotations)
        except Exception:  # noqa: BLE001
            _logger.exception("replay log: couldn't measure the stepped world")
            return
        self._after_step(world, snapshot, before_bytes, outcome, seconds)

    def _after_step(self, world: World, snapshot: tuple, before_bytes: bytes, outcome: StepOutcome, seconds: float) -> None:
        try:
            self._steps += 1
            fired = self._fired(outcome)
            if fired:
                self._enqueue(world, snapshot, before_bytes, outcome, fired, seconds)
        except Exception:  # noqa: BLE001
            _logger.exception("replay log: couldn't check or record the step")

    def _fired(self, outcome: StepOutcome) -> dict:
        fired = {}
        for name, check in self.triggers.items():
            try:
                values = check(outcome)
            except Exception as exc:  # noqa: BLE001 - one broken trigger mustn't hide the others
                _logger.exception("replay log: trigger %s failed", name)
                values = {"trigger_failed": f"{type(exc).__name__}: {exc}"}
            if values is not None:
                fired[name] = values
        if not any(self._steps - self._last_fired.get(name, -self.min_steps_between) >= self.min_steps_between for name in fired):
            if fired:
                _logger.info("replay log: %s fired but is rate-limited", ", ".join(fired))
            return {}
        for name in fired:
            self._last_fired[name] = self._steps
        return fired

    def _enqueue(self, world: World, snapshot: tuple, before_bytes: bytes, outcome: StepOutcome, fired: dict, seconds: float) -> None:
        seed, elapsed_before, steps_before = snapshot
        now = time.time()
        sidecar = {
            "version": RECORD_VERSION,
            "triggers": fired,
            "step_years": outcome.years,
            "seed": seed,
            "elapsed_years_before": elapsed_before,
            "steps_taken_before": steps_before,
            "node_density": world.node_density,
            "controls": _controls(world),
            "git": self._git,
            "step_seconds": seconds,
            "recorded_at": datetime.fromtimestamp(now, timezone.utc).isoformat(),
            "health_before": asdict(outcome.before),
            "health_after": None if outcome.after is None else asdict(outcome.after),
            # A step that raised left the world half-stepped: nothing to match a replay against.
            "fingerprint_after": None if outcome.error is not None else fingerprint(world),
            "annotations": outcome.annotations,
        }
        name = _record_name(seed, elapsed_before, fired, now)
        message = f"Replay log: {', '.join(sorted(fired))} -> {name}"
        _logger.warning("replay log: %s at %.2f Myr; writing %s", ", ".join(sorted(fired)), elapsed_before / 1e6, self.directory / name)
        with self._mutex:
            self._recent.append({"elapsed_years": elapsed_before, "message": message})
        self._queue.put(_PendingRecord(name, before_bytes, sidecar))
        self._ensure_writer()

    def _ensure_writer(self) -> None:
        with self._mutex:
            if self._writer is None or not self._writer.is_alive():
                self._writer = threading.Thread(target=self._write_loop, name="replay-log-writer", daemon=True)
                self._writer.start()

    def _write_loop(self) -> None:
        while True:
            pending = self._queue.get()
            try:
                self._write(pending)
                self._prune()
            except Exception:  # noqa: BLE001 - e.g. a full disk: the run goes on
                _logger.exception("replay log: couldn't write %s", pending.name)
            finally:
                self._queue.task_done()

    def flush(self) -> None:
        """Wait for every queued record to be written."""
        self._queue.join()

    def _write(self, pending: _PendingRecord) -> None:
        """Written into a temporary directory and renamed into place, so a reader never sees a
        record without its world."""
        temp = Path(tempfile.mkdtemp(prefix=_TEMP_PREFIX, dir=self.directory))
        try:
            for file_name, data in ((BEFORE_FILE, pending.before_bytes), (RECORD_FILE, json.dumps(pending.sidecar, indent=1, default=str).encode())):
                with open(temp / file_name, "wb") as handle:
                    handle.write(data)
                    handle.flush()
                    os.fsync(handle.fileno())
            os.replace(temp, self.directory / pending.name)
        except BaseException:
            shutil.rmtree(temp, ignore_errors=True)
            raise

    def _record_dirs(self) -> list[Path]:
        """Oldest first: names start with a UTC timestamp."""
        return sorted(p for p in self.directory.iterdir() if p.is_dir() and not p.name.startswith("."))

    def _prune(self) -> None:
        """Remove the oldest records until the rest fit in `max_bytes`; the newest always stays."""
        dirs = self._record_dirs()
        sizes = [sum(f.stat().st_size for f in d.iterdir() if f.is_file()) for d in dirs]
        total = sum(sizes)
        for directory, size in zip(dirs[:-1], sizes[:-1]):
            if total <= self.max_bytes:
                break
            shutil.rmtree(directory, ignore_errors=True)
            total -= size
            _logger.info("replay log: removed %s to stay under %.1f GB", directory.name, self.max_bytes / 1e9)

    def records(self) -> list[dict]:
        """A summary of each record on disk, newest first."""
        result = []
        for directory in reversed(self._record_dirs()):
            try:
                sidecar = json.loads((directory / RECORD_FILE).read_text())
            except (OSError, ValueError):
                continue
            result.append({
                "name": directory.name,
                "path": str(directory),
                "triggers": sorted(sidecar.get("triggers", {})),
                "seed": sidecar.get("seed"),
                "elapsed_years_before": sidecar.get("elapsed_years_before"),
                "step_years": sidecar.get("step_years"),
                "recorded_at": sidecar.get("recorded_at"),
            })
        return result

    def recent_events(self) -> list[dict]:
        """This process's records as `{elapsed_years, message}`, the shape of world events."""
        with self._mutex:
            return list(self._recent)


# --- Replay ------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ReplayResult:
    world: World
    sidecar: dict
    # World attributes whose fingerprint differs from the record's.
    mismatched: list[str]
    plate_count_matches: bool
    craton_ledger_matches: bool
    error: BaseException | None = None

    @property
    def matches(self) -> bool:
        return not self.mismatched and self.plate_count_matches and self.craton_ledger_matches


def load_record(record_dir: str | os.PathLike[str]) -> tuple[World, dict]:
    record_dir = Path(record_dir)
    sidecar = json.loads((record_dir / RECORD_FILE).read_text())
    if sidecar.get("version") != RECORD_VERSION:
        raise ValueError(f"replay log record version {sidecar.get('version')!r}, expected {RECORD_VERSION}")
    return persistence.load_world_bytes((record_dir / BEFORE_FILE).read_bytes()), sidecar


def replay(record_dir: str | os.PathLike[str]) -> ReplayResult:
    """Load a record's pre-step world, run its step, and compare against what was logged."""
    world, sidecar = load_record(record_dir)
    error = None
    try:
        for _ in step_world_progress(world, sidecar["step_years"]):
            pass
    except Exception as exc:  # noqa: BLE001 - reproducing a logged crash is a result, not a failure
        error = exc
    expected = sidecar.get("fingerprint_after")
    after = sidecar.get("health_after")
    if expected is None or after is None or error is not None:
        # A logged crash: the replay matches if it crashes too.
        crashed = error is not None and expected is None
        return ReplayResult(world, sidecar, [] if crashed else ["<step outcome>"], crashed, crashed, error)
    actual = fingerprint(world)
    mismatched = sorted(name for name in set(expected) | set(actual) if expected.get(name) != actual.get(name))
    return ReplayResult(
        world,
        sidecar,
        mismatched,
        plate_count_matches=len(world.plates) == after["plate_count"],
        craton_ledger_matches={k: float(v) for k, v in world.craton_ledger.items()} == after["craton_ledger"],
    )
