import threading

import pytest
from app.checkpoint import RunStatus
from app.run_control import WorldRunner


class _Harness:
    """A runner over a counter instead of a world. `saved` records each persisted status."""

    def __init__(self, fail_at=None, has_world=True):
        self.steps = 0
        self.saved = []
        self.fail_at = fail_at
        self.has_world = has_world
        self.stepped = threading.Event()
        self.gate = None  # set to an Event to hold the next step until it's set
        self.runner = WorldRunner(self.step_once, self.saved.append, interval_s=0.001)

    def step_once(self, step_years, should_step):
        if self.gate is not None:
            self.gate.wait(5)
        if not self.has_world:
            return False
        if should_step():
            if self.steps == self.fail_at:
                raise RuntimeError("boom")
            self.steps += 1
            self.stepped.set()
        return True

    def wait_for_steps(self, n):
        for _ in range(500):
            if self.steps >= n:
                return
            self.stepped.wait(0.01)
            self.stepped.clear()
        raise AssertionError(f"only {self.steps} steps")


@pytest.fixture
def harness():
    h = _Harness()
    yield h
    h.runner.close()


def _wait_stopped(runner):
    for _ in range(500):
        if not runner.status().running:
            return
        threading.Event().wait(0.01)
    raise AssertionError("still running")


def test_running_steps_until_stopped_and_saves_each_status(harness):
    harness.runner.set_status(RunStatus(running=True, step_years=10.0))
    harness.wait_for_steps(3)
    harness.runner.set_status(RunStatus(running=False, step_years=10.0))
    stopped_at = harness.steps
    threading.Event().wait(0.05)

    assert harness.steps <= stopped_at + 1
    assert harness.saved == [RunStatus(True, 10.0), RunStatus(False, 10.0)]


def test_stop_while_waiting_for_the_world_skips_the_step(harness):
    harness.gate = threading.Event()
    harness.runner.set_status(RunStatus(running=True))
    harness.runner.set_status(RunStatus(running=False))
    harness.gate.set()
    threading.Event().wait(0.05)
    assert harness.steps == 0


def test_a_failing_step_stops_the_run_and_saves_it_stopped():
    h = _Harness(fail_at=2)
    h.runner.set_status(RunStatus(running=True, step_years=5.0))
    _wait_stopped(h.runner)
    h.runner.close()

    assert h.steps == 2
    assert "boom" in h.runner.error()
    assert h.saved[-1] == RunStatus(running=False, step_years=5.0)


def test_no_world_stops_the_run():
    h = _Harness(has_world=False)
    h.runner.set_status(RunStatus(running=True))
    _wait_stopped(h.runner)
    h.runner.close()
    assert "no world" in h.runner.error()


def test_a_restored_status_runs_without_being_saved_again(harness):
    harness.runner.set_status(RunStatus(running=True), persist=False)
    harness.wait_for_steps(1)
    assert harness.saved == []


def test_shutdown_ends_the_loop_but_keeps_the_saved_status_running(harness):
    harness.runner.set_status(RunStatus(running=True))
    harness.wait_for_steps(1)
    harness.runner.close()
    after_close = harness.steps
    threading.Event().wait(0.05)

    assert harness.steps == after_close
    assert harness.runner.status().running
    assert harness.saved == [RunStatus(running=True)]


def test_a_failed_status_save_keeps_the_status_in_memory():
    def failing_save(status):
        raise OSError("volume gone")

    runner = WorldRunner(lambda years, should_step: True, failing_save, interval_s=0.001)
    runner.set_status(RunStatus(running=True))
    assert runner.status().running
    runner.close()
