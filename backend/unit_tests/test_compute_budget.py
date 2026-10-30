"""Tests for bin/compute_budget.py's run/wait decision and worker count (#704), on made-up
machine states: the real sysctl and ps readings depend on whatever else is running."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "compute_budget.py"
spec = importlib.util.spec_from_file_location("compute_budget", SCRIPT)
assert spec is not None and spec.loader is not None
budget = importlib.util.module_from_spec(spec)
sys.modules["compute_budget"] = budget  # dataclasses looks the module up by name
spec.loader.exec_module(budget)


def machine(available_gb=12.0, total_gb=16.0, pressure=1, cpus=10, load=1.0):
    return budget.Machine(total_gb, available_gb, pressure, cpus, load)


def sim(rss_gb=6.0):
    return budget.Run(1, "sim", rss_gb, "python bin/debug/replay.py")


def test_a_run_still_climbing_reserves_the_rest_of_its_peak():
    # Available memory alone would admit this run; a replay that has only just loaded will
    # take most of it.
    assert budget.decide(machine(available_gb=9.0), [sim(rss_gb=6.0)], 6.0)[0]
    assert not budget.decide(machine(available_gb=9.0), [sim(rss_gb=1.0)], 6.0)[0]


def test_critical_pressure_blocks_a_second_run_however_much_looks_free():
    ok, reason = budget.decide(machine(pressure=4), [sim()], 1.0)
    assert not ok
    assert "pressure" in reason


def test_critical_pressure_blocks_the_first_run_too():
    # With no agent runs detected, the pressure comes from the user's own apps.
    ok, reason = budget.decide(machine(pressure=4), [], 1.0)
    assert not ok
    assert "pressure" in reason


def test_warning_pressure_does_not_block_a_run_that_fits():
    assert budget.decide(machine(pressure=2), [], 1.0)[0]


def test_a_run_too_big_for_the_machine_goes_ahead_alone():
    # On the 8 GB machine a mature replay never fits the margin, and waiting wouldn't help.
    ok, reason = budget.decide(machine(available_gb=6.0, total_gb=8.0), [], 6.0)
    assert ok
    assert "alone" in reason


def test_workers_are_limited_by_memory_or_busy_cores_whichever_is_tighter():
    assert budget.worker_count(machine(available_gb=7.0), [], "stress") == 2
    assert budget.worker_count(machine(available_gb=14.0, load=7.0), [], "unit") == 2


@pytest.mark.parametrize(("total_gb", "ceiling"), [(16.0, 4), (8.0, 2)])
def test_installed_ram_caps_workers_on_an_idle_machine(total_gb, ceiling):
    idle = machine(available_gb=100.0, total_gb=total_gb, cpus=64, load=0.0)
    assert budget.worker_count(idle, [], "unit") == ceiling


def test_there_is_always_at_least_one_worker():
    assert budget.worker_count(machine(available_gb=1.0, load=20.0), [], "stress") == 1


@pytest.mark.parametrize(("action", "code", "out"), [("check", 3, "WAIT"), ("workers", 0, "1")])
def test_an_unreadable_machine_never_fails_a_test_run(monkeypatch, capsys, action, code, out):
    def denied():
        raise PermissionError("ps")

    monkeypatch.setattr(budget.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(budget, "read_machine", machine)
    monkeypatch.setattr(budget, "read_runs", denied)
    monkeypatch.setattr(sys, "argv", ["compute_budget.py", action, "--run", "unit"])
    assert budget.main() == code
    assert capsys.readouterr().out.startswith(out)


@pytest.mark.parametrize(
    ("command", "kind"),
    [
        ("/x/Python bin/debug/attribute_hm_ratchet.py save --myr 1", "sim"),
        ("/x/Python -u -c import sys;exec(eval(sys.stdin.readline()))", "stress"),
    ],
)
def test_parse_ps_finds_runs(command, kind):
    parent = "1 0 90000 /x/Python -m pytest stress_tests/ -q -n 2"
    runs = budget.parse_ps(f"{parent}\n2 1 900000 {command}\n")
    assert [(run.pid, run.kind) for run in runs] == [(2, kind)]


def test_parse_ps_ignores_the_users_app_however_large():
    assert budget.parse_ps("3 1 9000000 /x/Python -m app.desktop\n") == []


def test_parse_ps_ignores_a_shell_whose_command_mentions_a_replay():
    assert budget.parse_ps("5 1 2000 /bin/zsh -c python bin/debug/replay.py save\n") == []


def test_parse_ps_counts_an_unrecognised_large_python_process_as_a_simulation():
    assert [r.kind for r in budget.parse_ps("4 1 2000000 /x/Python /tmp/probe.py\n")] == ["sim"]
