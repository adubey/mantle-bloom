---
name: compute-budget
description: Decide whether a Mantle-Bloom replay, probe, sweep, stress run or unit test run can start now or must wait for other agents, and choose its test parallelism. Use before starting any compute-heavy or RAM-heavy run.
---

# Compute budget

Several agents share one Apple Silicon machine. **Each agent runs at most one heavy run at a
time**, and starts it only when the machine has room for it. Machine-wide, runs may overlap.

A heavy run is a replay, a probe, a sweep, a stress test run, or a unit or affected test run.
An A/B comparison is two runs, so run its arms one after the other. Lint, `quality_check.sh`,
the frontend build and a single test file without xdist (`pytest -n 0 <file>`) don't count.

## 1. Estimate the peak

| Run | Expected peak | Check with |
|---|---|---|
| Replay or probe, young world (about 40 Myr) | 2–3 GB | `--run sim --gb 3` |
| Replay or probe, mature world (200+ Myr) or CFD wind | 6–7 GB | `--run sim --gb 7` |
| Replay or probe, age unknown | 6 GB | `--run sim` |
| Unit or affected test worker | 0.8–2.7 GB each | `--run unit` |
| Stress test worker | 0.7–3.7 GB each | `--run stress` |

For a replay, use the `peak_rss` the same save reached in `analysis/` when you can find it,
since that beats any default. A run shows its real size late, so don't judge it by its RSS
in the first few minutes.

Replay figures come from the 454 peak-RSS values under `analysis/`: median 3.6 GB, 90th
percentile 5.8 GB, maximum 7.2 GB. Test figures are per-worker peaks from one full run of
each suite on the 16 GB Mac mini (#704): 2,057 unit tests at `-n 3` in 14 minutes and 58
stress tests at `-n 2` in 10 minutes, both alongside two replays.

## 2. Check

```bash
python3 bin/compute_budget.py check --run sim --gb <peak>   # replay or probe
python3 bin/compute_budget.py check --run unit              # unit or affected tests
python3 bin/compute_budget.py check --run stress
```

The script lists what it found running and ends with `RUN` (exit 0) or `WAIT` (exit 3), plus
a worker count for tests. See its docstring for how it estimates free memory. Check that none
of the listed runs is yours: the script can't tell agents apart, and a second run of your own
is over budget even if it says `RUN`.

## 3. Run or wait

- **`RUN`:** start immediately, because the answer goes stale as other agents start. The
  test scripts choose their own worker count. Pass `-n N` only to go lower, or when you've
  decided on a value from the `check` output:

  ```bash
  ./bin/affected_test.sh --base origin/main -n 2
  ./bin/stress_test.sh -n 2
  ```

- **`WAIT`:** run this loop in the background, then start your run when it exits. Do other
  work in the meantime.

  ```bash
  until python3 bin/compute_budget.py check --run sim --gb 6; do sleep 120; done
  ```

Two agents that check close together can both get `RUN`, because `check` can't see a run until
it's started. `bin/debug/` scripts and test workers are listed as soon as they start, but a
script anywhere else is listed only once it passes 1 GB, which can take minutes. An agent the
user prompted directly doesn't need to allow for this, since those rarely start together.
Automatically spawned agents do: whatever spawns a batch of them should start their heavy
runs one at a time, starting each only once `check` lists the one before.

The user can override a `WAIT`: if they tell you to run anyway, start the run.

If `check` says `WAIT` for critical memory pressure with no runs listed, the user's own apps
are short of memory. Keep waiting, and tell the user if it lasts more than 30 minutes. If it says it
can't read the machine, it's probably in a sandbox that denies `ps`, so run it outside the
sandbox rather than looping: the wait would never end.

Never stop or pause another agent's run to make room. If you're the one holding things up,
and your work can use less, shrink it: fewer workers, a shorter horizon or a younger save.

## 4. Size the machine

The check adapts to the machine it's on: the worker ceiling is installed RAM divided by 4, and
free memory and load can lower it. These are what to expect when planning:

| Machine | Test workers, at most | Heavy runs at once |
|---|---|---|
| Mac mini, 16 GB, 10 cores | 4, fewer when busy | 2 mature replays, or 1 plus a small test run |
| MacBook Neo, 8 GB, 6 cores | 2, usually 1 | 1 mature replay, alone; or 2 young replays |

A mature replay nearly fills the Neo, so `check` lets it start only with nothing else
running. Prefer younger saves there, or run mature replays on the mini.

On a loaded machine macOS compresses and swaps before anything fails, so a run slows down
rather than crashes. If your own run stops making progress and `sysctl vm.swapusage` keeps
climbing, say so instead of starting more work.

## 5. Report

Record peak RSS when a report asks for it, and say what else was running. Runtime from a run
that overlapped others isn't comparable: when runtime is a measure in a comparison, wait until
the `check` output lists no other runs before starting each arm.
