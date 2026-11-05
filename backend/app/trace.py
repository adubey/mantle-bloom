"""Nested wall-clock spans for the run overlay's phase text and for profiling (#801).

A *span* names a stretch of work and the user-visible group it belongs to. Spans nest by call
stack, but the group is stated at the span, not inferred from its parent: climate, rivers, lakes
and glaciers all run inside `erosion.apply_erosion` yet belong to "Hydrometeorology", and
`docs/simulation-model.md`-style call order would put them under erosion.

Two kinds of span share one mechanism:

* `user=True` spans are the labels the run overlay shows (e.g. "Updating faults"). Entering one
  tells the active recorder's `sink` so the label appears when the work *starts*; the old
  progress labels were yielded after the work they named. Leaving a nested user span re-announces
  the enclosing user span, so the overlay doesn't keep naming finished work.
* Other spans are profiling detail, recorded only. `bin/debug/profile_phases.py` wraps functions
  in these at runtime, so production code carries none of them.

Nothing records unless a `Recorder` is active in the current context (`recording()`), and
inactive cost is one `ContextVar.get()` per span, so the always-on `@traced` decorators on
coarse functions are safe on the step path. A span entered under the same label as its parent
merges into the parent (recursion, or a function reached through two wrappers) instead of
nesting, except that a user-visible span is never merged into a profiling-only one.

`run_traced` runs a progress generator on a worker thread with a recorder whose sink feeds the
same queue as the generator's own messages. The step and generate endpoints use it so a label
from deep inside erosion reaches the client while the work runs; a generator can only yield
between its own statements. The generators themselves (and their tests, replay log and
`_drain_progress`) are unchanged.
"""

from __future__ import annotations

import contextlib
import contextvars
import functools
import queue
import threading
import time
from collections.abc import Callable, Iterator

# Seconds one step takes in each user-visible phase on the young reference world (seed 0,
# 123 K nodes, UI-default densities, diagnostic wind; analysis/issue801/traces.md). They set the
# overlay bar's fraction: a phase's weight is added to the bar when it finishes. A mature world
# spends its time differently (lakes and ice dominate), so the fraction is a guide, not a clock.
# Exclusive of nested user spans: erosion's own work, rivers without their glaciers.
STEP_WEIGHTS: dict[str, float] = {
    "Solving mantle flow and plate motion": 0.53,
    "Detecting collisions": 0.33,
    "Reshaping plates": 0.44,
    "Updating faults": 1.08,
    "Merging and splitting plates": 0.05,
    "Filling seafloor gaps": 0.45,
    "Moving magma": 0.07,
    "Raising mantle plumes": 0.03,
    "Solving wind and ocean currents": 2.7,  # CFD wind only, see step_weights
    "Resolving climate": 0.49,
    "Resolving rivers and lakes": 0.77,
    "Marching glaciers": 0.44,
    "Eroding and depositing sediment": 0.91,
    "Raising volcanoes": 0.05,
    "Forming natural resources": 0.05,
    "Balancing sea level": 0.17,
}
GENERATE_WEIGHTS: dict[str, float] = {
    "Laying out continents and mountain belts": 4.7,
    "Dividing the sphere into plates": 3.2,
    "Building pseudo-history": 5.4,  # synthetic history only, see generate_weights
    "Shaping initial bathymetry": 0.6,
    "Dating the initial rocks": 0.35,
    "Seeding rivers and glaciers": 1.4,
    "Resolving climate": 0.5,
}

# Per-user-span group names, "<bar>/<group>" (analysis/issue801/traces.md).
LAYOUT = "Building world/Laying out the planet"
SEEDING = "Building world/Seeding the surface"
PLATES = "Stepping/Moving plates"
HYDRO = "Stepping/Hydrometeorology (climate, water & ice)"
SURFACE = "Stepping/Reshaping the surface"


# Labels whose spans only run with a given setting, so the weights drop them when it is off.
_PLATE_MOVEMENT_LABELS = frozenset({
    "Solving mantle flow and plate motion", "Detecting collisions", "Reshaping plates", "Updating faults",
    "Merging and splitting plates", "Filling seafloor gaps", "Moving magma", "Raising mantle plumes",
    "Raising volcanoes",
})
_CLIMATE_LABELS = frozenset({
    "Solving wind and ocean currents", "Resolving climate", "Resolving rivers and lakes", "Marching glaciers",
    "Eroding and depositing sediment", "Forming natural resources",
})


def step_weights(wind_model: str, plate_movement: bool = True, climate_biomes: bool = True) -> dict[str, float]:
    """`STEP_WEIGHTS` for a world: only the phases that will run count, or the bar would top
    out early. The CFD solver's label needs `wind_model == "cfd"`, the plate phases
    `simulate_plate_movement`, and the climate, water and erosion phases `simulate_climate_biomes`."""
    skipped: set[str] = set()
    if wind_model != "cfd":
        skipped.add("Solving wind and ocean currents")
    if not plate_movement:
        skipped |= _PLATE_MOVEMENT_LABELS
    if not climate_biomes:
        skipped |= _CLIMATE_LABELS
    return {label: weight for label, weight in STEP_WEIGHTS.items() if label not in skipped}


def generate_weights(synthetic_relief: bool, cratons_layout: bool = True, seeds_rivers: bool = True) -> dict[str, float]:
    """`GENERATE_WEIGHTS` for a generation. The layout phases only run for a random world laid out
    cratons first (not a sketch or premade world), pseudo-history only with synthetic relief on top
    of that, and river seeding only when `cratons_first` is set."""
    skipped: set[str] = set()
    if not cratons_layout:
        skipped |= {"Laying out continents and mountain belts", "Dividing the sphere into plates"}
    if not (synthetic_relief and cratons_layout):
        skipped.add("Building pseudo-history")
    if not seeds_rivers:
        skipped.add("Seeding rivers and glaciers")
    return {label: weight for label, weight in GENERATE_WEIGHTS.items() if label not in skipped}


class _Frame:
    __slots__ = ("label", "group", "user", "start", "child_s")

    def __init__(self, label: str, group: str | None, user: bool, start: float) -> None:
        self.label = label
        self.group = group
        self.user = user
        self.start = start
        self.child_s = 0.0


class Recorder:
    """Accumulates `path -> [inclusive_s, exclusive_s, calls]` for the spans run under it.

    `sink(label, group)` is called when a user span starts (or is resumed after a nested one
    ends) and its label differs from the last one announced; `on_finish(label)` when a user
    span ends. Exclusive time is own time minus time inside child spans, so the rows of one
    parent sum to its inclusive time."""

    def __init__(
        self,
        sink: Callable[[str, str | None], None] | None = None,
        on_finish: Callable[[str], None] | None = None,
    ) -> None:
        self.sink = sink
        self.on_finish = on_finish
        self.stack: list[_Frame] = []
        self.stats: dict[tuple[str, ...], list[float]] = {}
        self.groups: dict[str, str | None] = {}
        self.fired: set[str] = set()
        self._announced: str | None = None

    def reset(self) -> None:
        self.stats = {}
        self._announced = None

    def snapshot(self) -> dict[tuple[str, ...], list[float]]:
        return {path: list(row) for path, row in self.stats.items()}

    def _announce(self, frame: _Frame) -> None:
        if self.sink is not None and frame.label != self._announced:
            self._announced = frame.label
            self.sink(frame.label, frame.group)

    def enter(self, label: str, group: str | None, user: bool) -> _Frame | None:
        self.fired.add(label)
        if group is not None:
            self.groups[label] = group
        top = self.stack[-1] if self.stack else None
        if top is not None and top.label == label and (top.user or not user):
            # Re-entry under the same label merges into the open span, unless that span is
            # profiling-only and this one is user-visible: merging would swallow its announcement.
            return None
        frame = _Frame(label, group, user, time.perf_counter())
        self.stack.append(frame)
        if user:
            self._announce(frame)
        return frame

    def exit(self, frame: _Frame) -> None:
        elapsed = time.perf_counter() - frame.start
        path = tuple(f.label for f in self.stack)
        self.stack.pop()
        if self.stack:
            self.stack[-1].child_s += elapsed
        row = self.stats.setdefault(path, [0.0, 0.0, 0])
        row[0] += elapsed
        row[1] += elapsed - frame.child_s
        row[2] += 1
        if frame.user:
            if self.on_finish is not None:
                self.on_finish(frame.label)
            parent = next((f for f in reversed(self.stack) if f.user), None)
            if parent is not None:
                self._announce(parent)


_active: contextvars.ContextVar[Recorder | None] = contextvars.ContextVar("trace_recorder", default=None)


@contextlib.contextmanager
def recording(recorder: Recorder | None = None) -> Iterator[Recorder]:
    """Make `recorder` (a new one by default) the active recorder for the enclosed code."""
    recorder = recorder if recorder is not None else Recorder()
    token = _active.set(recorder)
    try:
        yield recorder
    finally:
        _active.reset(token)


@contextlib.contextmanager
def span(label: str, group: str | None = None, user: bool = False) -> Iterator[None]:
    """Time the enclosed block under `label`. A no-op unless a recorder is active."""
    recorder = _active.get()
    if recorder is None:
        yield
        return
    frame = recorder.enter(label, group, user)
    try:
        yield
    finally:
        if frame is not None:
            recorder.exit(frame)


def wrap(fn: Callable, label: str, group: str | None = None, user: bool = False) -> Callable:
    """`fn` timed under `label`; the building block of `traced` and of the profiler's hooks."""

    @functools.wraps(fn)
    def traced_fn(*args, **kwargs):
        recorder = _active.get()
        if recorder is None:
            return fn(*args, **kwargs)
        frame = recorder.enter(label, group, user)
        try:
            return fn(*args, **kwargs)
        finally:
            if frame is not None:
                recorder.exit(frame)

    return traced_fn


def traced(label: str, group: str | None = None, user: bool = False) -> Callable[[Callable], Callable]:
    """Decorator form of `span`: `@traced("Updating faults", PLATES, user=True)`."""

    def decorate(fn: Callable) -> Callable:
        return wrap(fn, label, group, user)

    return decorate


class PhaseTracker:
    """The overlay bar's fraction: the share of `weights` already finished, never decreasing
    and held below 1 so the bar only completes when the run does. Feed it `Recorder.on_finish`."""

    def __init__(self, weights: dict[str, float]) -> None:
        self.weights = weights
        self.total = sum(weights.values()) or 1.0
        self.finished: set[str] = set()
        self.fraction = 0.0

    def advance(self, fraction: float) -> float:
        """Raise the fraction to a progress line's own, so span fractions never fall behind it."""
        self.fraction = max(self.fraction, min(fraction, 0.99))
        return self.fraction

    def finish(self, label: str) -> None:
        if label in self.weights:
            self.finished.add(label)
            share = sum(self.weights[name] for name in self.finished) / self.total
            self.fraction = max(self.fraction, min(share, 0.99))


_DONE = object()


def run_traced(progress: Iterator, weights: dict[str, float]) -> Iterator[tuple]:
    """Iterate `progress` on a worker thread, merging span events into its messages.

    Yields what `progress` yields, plus a `("phase", fraction, label, group)` message for each
    user-span change. Once the first of those has been sent, the generator's own
    `("progress", fraction, phase)` messages are dropped: their labels arrive after the work they
    name, and their fractions count units that no longer match the span labels. Until then they
    pass through, so a run with no spans (a scripted debug world) keeps its old progress. The fraction is the share of `weights` already
    finished, never decreasing, and below 1 until the stream ends. `("heartbeat", f)` and
    `("done", ...)` messages pass through, a heartbeat carrying the latest fraction. An
    exception in the worker is re-raised here. Closing the iterator before it ends (or an error on
    the consumer's side) blocks until the worker has finished, so a caller that releases a lock
    after closing it never races the work."""
    events: queue.Queue = queue.Queue()
    tracker = PhaseTracker(weights)
    spans_seen = False

    def work() -> None:
        recorder = Recorder(
            sink=lambda label, group: events.put(("phase", label, group)),
            on_finish=lambda label: events.put(("finished", label)),
        )
        try:
            with recording(recorder):
                for message in progress:
                    events.put(("message", message))
        except BaseException as exc:  # noqa: BLE001 -- re-raised in the consumer
            events.put(("error", exc))
        finally:
            events.put(_DONE)

    thread = threading.Thread(target=work, name="traced-progress", daemon=True)
    thread.start()
    try:
        while True:
            event = events.get()
            if event is _DONE:
                break
            kind = event[0]
            if kind == "error":
                raise event[1]
            if kind == "finished":
                tracker.finish(event[1])
            elif kind == "phase":
                spans_seen = True
                yield ("phase", tracker.fraction, event[1], event[2])
            else:
                message = event[1]
                if message[0] == "progress":
                    if spans_seen:
                        continue
                    # Before any span has spoken: pass the generator's line through, but never
                    # below it later (the first phase message would otherwise drop the bar to 0).
                    # The tracker keeps the span-driven cap below 1; the generator's own line,
                    # including a final 1.0, goes out unchanged.
                    message = ("progress", max(message[1], tracker.advance(message[1])), *message[2:])
                elif message[0] == "heartbeat":
                    message = ("heartbeat", tracker.fraction)
                yield message
    finally:
        # Closed early or failed on the consumer's side: the work must not outlive whoever holds
        # the session lock and compute gate, so wait for it before they are released.
        thread.join()
