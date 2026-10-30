#!/usr/bin/env python3
"""Decide whether a compute-heavy run can start now, and with how many test workers (#704).

Several agents share one Apple Silicon machine, so each agent runs one heavy job at a time
and checks this script first. See .agents/skills/compute-budget/SKILL.md for the procedure.

    bin/compute_budget.py check --run sim [--gb 3.5]   # replay or probe
    bin/compute_budget.py check --run unit             # unit or affected tests
    bin/compute_budget.py check --run stress
    bin/compute_budget.py workers --run unit           # print only pytest's -n value

`check` prints RUN or WAIT and exits 0 or 3. `workers` never waits: bin/unit_test.sh and
friends use it for their default -n, and it prints at least 1, even if the machine can't be
read (a sandboxed shell may deny ps; `check` then says WAIT). Off macOS (CI) there's
nothing to share, so `check` always says RUN and `workers` prints `auto`.

Free memory comes from kern.memorystatus_level, the percentage `memory_pressure` reports,
which unlike vm_stat's free pages counts reclaimable file cache. It doesn't know about a run
that has started but not reached its peak, so each detected run also reserves the rest of
its expected peak. That reservation is what stops agents that check a few seconds apart from
all starting on the same memory, once the first of them is visible in `ps`.
"""

from __future__ import annotations

import argparse
import math
import platform
import re
import subprocess
import sys
from dataclasses import dataclass

GIB = 1024**3
KIB_PER_GB = 1024**2

# Expected peak per process. `sim` is the 90th percentile of the peak RSS values under
# analysis/ (#704), since mature worlds reach it and a run only reveals its size late. The
# test figures sit above each suite's mean worker peak (1.9 GB unit at -n 3, 2.2 GB stress at
# -n 2) and below its worst (2.7, 3.7): a few large tests set those, and rarely coincide.
PEAK_GB = {"sim": 6.0, "unit": 2.5, "stress": 2.5}

# Most test workers a machine gets however idle it is: 4 on 16 GB, 2 on 8 GB.
GB_PER_WORKER_CEILING = 4.0

# Below this RSS a python process with no recognised command line isn't worth reserving for.
UNKNOWN_HEAVY_GB = 1.0

# kern.memorystatus_vm_pressure_level: 1 normal, 2 warning, 4 critical. Only critical blocks a
# run: this machine sits at warning from the user's own apps, with no agent at fault (#704).
CRITICAL_PRESSURE = 4
XDIST_WORKER = "exec(eval(sys.stdin.readline()))"
# Matches the executable only: a shell wrapper's text often mentions python and bin/debug/.
PYTHON_EXE = re.compile(r"\S*python[\d.]*(\s|$)", re.IGNORECASE)
APP_SERVER = re.compile(r"app\.desktop|uvicorn|app\.main")


@dataclass
class Machine:
    total_gb: float
    available_gb: float
    pressure: int
    cpus: int
    load: float


@dataclass
class Run:
    pid: int
    kind: str
    rss_gb: float
    command: str

    @property
    def reserve_gb(self) -> float:
        return max(0.0, PEAK_GB[self.kind] - self.rss_gb)


def read_machine() -> Machine:
    names = [
        "hw.memsize",
        "kern.memorystatus_level",
        "kern.memorystatus_vm_pressure_level",
        "hw.ncpu",
        "vm.loadavg",
    ]
    out = subprocess.run(
        ["sysctl", "-n", *names], check=True, capture_output=True, text=True
    ).stdout.splitlines()
    total_gb = int(out[0]) / GIB
    return Machine(
        total_gb=total_gb,
        available_gb=total_gb * int(out[1]) / 100,
        pressure=int(out[2]),
        cpus=int(out[3]),
        load=float(out[4].strip("{} ").split()[0]),
    )


def classify(command: str, rss_gb: float, parent_command: str) -> str | None:
    """The kind of heavy run a process is, or None if it isn't one. A process that looks
    like a simulation gets the `sim` peak; the user's own app is left out, because it's
    always there and its memory is already in the machine's available figure."""
    if not PYTHON_EXE.match(command) or APP_SERVER.search(command):
        return None
    if XDIST_WORKER in command:
        return "stress" if "stress_tests" in parent_command else "unit"
    if "bin/debug/" in command or rss_gb >= UNKNOWN_HEAVY_GB:
        return "sim"
    return None


def parse_ps(text: str) -> list[Run]:
    rows = []
    for line in text.splitlines():
        fields = line.split(None, 3)
        if len(fields) == 4:
            pid, ppid, rss_kib, command = fields
            rows.append((int(pid), int(ppid), int(rss_kib) / KIB_PER_GB, command))
    commands = {pid: command for pid, _, _, command in rows}
    runs = []
    for pid, ppid, rss_gb, command in rows:
        kind = classify(command, rss_gb, commands.get(ppid, ""))
        if kind is not None:
            runs.append(Run(pid, kind, rss_gb, command))
    return runs


def read_runs() -> list[Run]:
    out = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss=,command="],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return parse_ps(out)


def free_gb(machine: Machine, runs: list[Run]) -> float:
    """Memory left for a new run: what's available, less what running jobs will still claim
    and a margin for the user's apps and the compressor, which degrade badly when squeezed."""
    margin_gb = max(1.0, machine.total_gb / 8)
    return machine.available_gb - sum(r.reserve_gb for r in runs) - margin_gb


def worker_count(machine: Machine, runs: list[Run], kind: str) -> int:
    """Workers that fit in free memory and on idle cores, keeping one core for the user, up
    to the machine's RAM ceiling. Load, not core count, is the CPU limit: a replay holds a
    core for its whole run."""
    by_ram = math.floor(machine.total_gb / GB_PER_WORKER_CEILING)
    by_memory = math.floor(free_gb(machine, runs) / PEAK_GB[kind])
    by_cpu = machine.cpus - 1 - math.ceil(machine.load)
    return max(1, min(by_ram, by_memory, by_cpu))


def decide(machine: Machine, runs: list[Run], need_gb: float) -> tuple[bool, str]:
    # Pressure comes first: with no agent runs, it's the user's apps that are short of memory.
    if machine.pressure >= CRITICAL_PRESSURE:
        return False, "memory pressure is critical"
    free = free_gb(machine, runs)
    if not runs:
        if free < need_gb:
            return True, f"nothing else is running, but {need_gb:.1f} GB may swap: run it alone"
        return True, "nothing else is running"
    if free < need_gb:
        return False, f"{need_gb:.1f} GB needed, {free:.1f} GB free"
    return True, f"{need_gb:.1f} GB needed, {free:.1f} GB free"


def report(machine: Machine, runs: list[Run]) -> None:
    print(
        f"machine: {machine.total_gb:.0f} GB, {machine.available_gb:.1f} GB available, "
        f"pressure {machine.pressure}, {machine.cpus} cores, load {machine.load:.1f}"
    )
    for run in runs:
        print(
            f"running: {run.kind} pid {run.pid}, {run.rss_gb:.1f} GB, "
            f"reserving {run.reserve_gb:.1f} GB: {run.command[-90:]}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("action", choices=["check", "workers"])
    parser.add_argument("--run", choices=sorted(PEAK_GB), required=True)
    parser.add_argument("--gb", type=float, help="expected peak, if better known than the default")
    args = parser.parse_args()

    if platform.system() != "Darwin":
        print("RUN" if args.action == "check" else "auto")
        return 0

    # The budget is advisory, so a sandbox that denies ps mustn't fail a test run (#704).
    try:
        machine, runs = read_machine(), read_runs()
    except (OSError, subprocess.CalledProcessError, ValueError, IndexError) as error:
        print(f"compute_budget: can't read this machine's state: {error}", file=sys.stderr)
        if args.action == "workers":
            print(1)
            return 0
        print("WAIT: can't read this machine's state")
        return 3

    workers = worker_count(machine, runs, args.run)
    if args.action == "workers":
        print(workers)
        return 0

    is_sim = args.run == "sim"
    need_gb = args.gb or PEAK_GB[args.run]
    ok, reason = decide(machine, runs, need_gb)
    report(machine, runs)
    print(f"{'RUN' if ok else 'WAIT'}: {reason}" + ("" if is_sim else f"; -n {workers}"))
    return 0 if ok else 3


if __name__ == "__main__":
    sys.exit(main())
