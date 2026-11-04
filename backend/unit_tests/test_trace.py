import re
import time
from pathlib import Path

import pytest

from app import trace


def test_span_is_a_no_op_without_a_recorder():
    calls = []

    @trace.traced("Updating faults", trace.PLATES, user=True)
    def work():
        calls.append(1)
        return 7

    with trace.span("anything"):
        assert work() == 7
    assert calls == [1]


def test_nested_spans_record_exclusive_time_that_sums_to_the_parent():
    with trace.recording() as rec:
        with trace.span("outer"):
            time.sleep(0.01)
            with trace.span("inner"):
                time.sleep(0.01)
    outer = rec.stats[("outer",)]
    inner = rec.stats[("outer", "inner")]
    assert inner[2] == 1 and outer[2] == 1
    assert outer[0] == pytest.approx(outer[1] + inner[0])
    assert inner[1] == pytest.approx(inner[0])


def test_same_label_re_entry_merges_into_the_open_span():
    @trace.traced("Computing climate")
    def recurse(depth):
        if depth:
            recurse(depth - 1)

    with trace.recording() as rec:
        recurse(3)
    assert list(rec.stats) == [("Computing climate",)]
    assert rec.stats[("Computing climate",)][2] == 1


def test_group_is_stated_at_the_span_not_inferred_from_the_parent():
    with trace.recording() as rec:
        with trace.span("Eroding and depositing sediment", trace.SURFACE, user=True):
            with trace.span("Resolving climate", trace.HYDRO, user=True):
                pass
    assert rec.groups["Resolving climate"] == trace.HYDRO
    assert rec.groups["Eroding and depositing sediment"] == trace.SURFACE


def test_sink_announces_on_start_and_re_announces_the_parent_after_a_nested_span():
    announced = []
    finished = []
    with trace.recording(trace.Recorder(sink=lambda label, group: announced.append(label), on_finish=finished.append)):
        with trace.span("erosion", user=True):
            with trace.span("detail"):  # profiling-only: never announced
                pass
            with trace.span("rivers", user=True):
                with trace.span("rivers", user=True):  # re-entry merges: no second announcement
                    pass
            with trace.span("rivers", user=True):  # a new span after the parent was re-announced
                pass
    assert announced == ["erosion", "rivers", "erosion", "rivers", "erosion"]
    assert finished == ["rivers", "rivers", "erosion"]


def _steps():
    @trace.traced("Updating faults", trace.PLATES, user=True)
    def faults():
        yield_value = 1
        return yield_value

    @trace.traced("Resolving climate", trace.HYDRO, user=True)
    def climate():
        return 2

    faults()
    yield ("progress", 0.5, "Moving plates")  # the generator's own label: dropped once spans spoke
    yield ("heartbeat", 0.9)
    climate()
    yield ("done", "world")


def test_run_traced_passes_generator_progress_through_until_a_span_has_spoken():
    def no_spans():
        yield ("progress", 0.5, "Creating continents")
        yield ("done", "world")

    assert list(trace.run_traced(no_spans(), {})) == [("progress", 0.5, "Creating continents"), ("done", "world")]


def test_run_traced_replaces_progress_with_phase_messages_and_weighted_fractions():
    weights = {"Updating faults": 3.0, "Resolving climate": 1.0}
    messages = list(trace.run_traced(_steps(), weights))
    kinds = [m[0] for m in messages]
    assert kinds == ["phase", "heartbeat", "phase", "done"]
    assert messages[0] == ("phase", 0.0, "Updating faults", trace.PLATES)
    # Updating faults (3 of 4) had finished by the time the heartbeat went out and the next
    # phase began, and the fraction never reaches 1 before the stream ends.
    assert messages[1] == ("heartbeat", 0.75)
    assert messages[2] == ("phase", 0.75, "Resolving climate", trace.HYDRO)
    assert messages[3] == ("done", "world")


def test_run_traced_re_raises_a_worker_exception_in_the_consumer():
    def boom():
        yield ("progress", 0.1, "x")
        raise RuntimeError("step failed")

    with pytest.raises(RuntimeError, match="step failed"):
        list(trace.run_traced(boom(), {}))


def test_weights_depend_on_the_wind_model_and_on_synthetic_relief():
    assert "Solving wind and ocean currents" in trace.step_weights("cfd")
    assert "Solving wind and ocean currents" not in trace.step_weights("diagnostic")
    assert "Building pseudo-history" in trace.generate_weights(True)
    assert "Building pseudo-history" not in trace.generate_weights(False)


def test_every_user_visible_span_label_has_learn_links_in_the_frontend():
    # The run overlay shows a label's learn-more links by exact name (phaseInfo.ts), so a span
    # label added without an entry would silently show no links.
    app_dir = Path(trace.__file__).parent
    labels = set()
    for source in app_dir.rglob("*.py"):
        labels.update(re.findall(r'@trace\.traced\("([^"]+)",[^\n]*user=True', source.read_text()))
    assert labels, "found no user-visible spans to check"
    phase_info = (app_dir.parents[1] / "frontend" / "src" / "phaseInfo.ts").read_text()
    missing = sorted(label for label in labels if f'"{label}": [' not in phase_info)
    assert not missing, f"no learn links in frontend/src/phaseInfo.ts for: {missing}"


def test_every_weighted_label_is_a_real_user_visible_span():
    app_dir = Path(trace.__file__).parent
    labels = set()
    for source in app_dir.rglob("*.py"):
        labels.update(re.findall(r'@trace\.traced\("([^"]+)",[^\n]*user=True', source.read_text()))
    stale = sorted((set(trace.STEP_WEIGHTS) | set(trace.GENERATE_WEIGHTS)) - labels)
    assert not stale, f"weights for labels that no span emits: {stale}"
