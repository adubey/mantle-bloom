# Issue #461: belt marks (RFC 0009)

Replay A/B of the belt marks against `origin/main`. Each run generates a default world (today's
plates-first layout, density 1) for seeds 7, 11 and 2024 and steps it to 1 Gyr in 5 Myr steps,
with a row every 100 Myr. Two runs:

- **main**: `origin/main` at `2c6f668` (#471 merged).
- **belts**: this branch at that commit plus the belt marks, after the PR #474 review fixes
  (no belt on a cell holding craton; a cell retyped oceanic loses its mark, and belt area counts
  continental cells only). Every measure but belt area is identical at every checkpoint to the
  run before those fixes, which counted oceanic-retyped cells as belt.

```bash
backend/.venv/bin/python bin/debug/measure_belt_marks.py run --seeds 7,11,2024 \
  --years 1e9 --every 100e6 --out belts.jsonl
# The baseline: a copy of the script inside the origin/main worktree, which puts its own
# backend first on sys.path (each row records `model_file`).
backend/.venv/bin/python bin/debug/measure_belt_marks.py summarize main.jsonl belts.jsonl
```

The compact comparison is [`compare-belts.json`](compare-belts.json) (every 100 Myr
checkpoint, both runs); the full JSONL is not checked in.

## At 1 Gyr

| Measure | Seed 7 main / belts | Seed 11 main / belts | Seed 2024 main / belts |
|---|---:|---:|---:|
| Formed volume (cumulative) | 331 / 251 M km³ (−24%) | 482 / 412 M km³ (−15%) | 905 / 654 M km³ (−28%) |
| Craton area | 90.2 / 91.1 M km² | 71.7 / 77.4 M km² | 69.5 / 70.3 M km² |
| Genuine continental area | 100.9 / 99.0 M km² | 79.5 / 91.2 M km² | 94.0 / 104.2 M km² |
| Belt area | — / 19.2 M km² (19% of continental) | — / 26.9 M km² (30%) | — / 43.4 M km² (42%) |
| Peak belt area | 67.4 M km² at 200 Myr (35%) | 43.6 M km² at 300 Myr (36%) | 53.7 M km² at 600 Myr (47%) |
| Formation on marked columns | — / 0 | — / 0 | — / 0 |
| Craton ledger balance | ≤ 3e-6 / ≤ 1e-5 km³ | ≤ 1e-6 / ≤ 1e-5 km³ | ≤ 2e-6 / ≤ 1e-6 km³ |
| Continental-material balance | −5.6 / −3.9 M km³ | −9.6 / −13.5 M km³ | −1.7 / −1.8 M km³ |

**The invariant holds.** The audit (`FormationAudit` in the script wraps `cratons.update`)
found no craton formed on a column marked as belt, on any seed at any step.

**Formation falls, craton area doesn't.** Welded and rifted crust no longer cratonises, so 15–28%
less volume forms over 1 Gyr. The live craton area at 1 Gyr ends level with `main` on seeds 7
and 2024, and 8% above it on seed 11. Formed volume counts every column that cratonises, including
columns that form and are later lost; it measures formation activity, not the standing
inventory. Which formation the marks remove (short-lived craton near active margins, or lasting
craton) is not separated here.

**Nothing else moves until formation does.** No craton forms on either run in the first
100 Myr (`CRATON_FORMATION_MYR` is 200 Myr of quiet), and at 100 Myr the two runs match on
every measure but belt area. By 200 Myr formation has started and craton area differs by up to
1.2%, while continental area is still identical. After that the runs diverge as different columns cratonise and
cratons resist differently: continental area at 1 Gyr is 2% lower on seed 7 and 11–15% higher
on seeds 11 and 2024. The gap grows with time on two seeds and changes sign across seeds.
Run-to-run noise is about ±2% at 43 Myr (#315) and is not measured at 1 Gyr. This is evidence
of trajectory divergence, not a demonstrated effect on continental area.

**Belts are large on one seed.** At 1 Gyr, belts cover 19%, 30% and 42% of genuine continental
area, against the 10–30% that RFC 0009 targets for the generated layout at t = 0; they peak at
35–47%. Most new marks come from crust
landing on upper plates at collision fronts. On seed 7 over the first 20 Myr, the new belt area
by writer was: underthrust 10.1 M km², upper-plate thrust-up 8.0 M km², failed rifts 7.8 M km²,
suture accretion 7.5 M km², merge stacking 1.2 M km², split cuts 0.65 M km². RFC 0009 expects
G1's risk to be belts that are lost too fast; here they persist. Belt area peaks at 200–600 Myr
and falls by a fifth to over a half by 1 Gyr, as continents shrink.

**Ledgers.** The craton ledger closes to round-off in both runs. The continental-material
balance error is of the same size on `main` (−1.7 to −9.6 M km³ at 1 Gyr), so it predates this
change.

## Step time

The issue asks for step time to change by less than 2% at density 4. A whole-step A/B on this
machine was not usable: the same benchmark (seed 7, density 4, 2 Myr steps, 2 warm-up and 8
timed steps) gave a median of either ~16 s or ~3.4 s per step, depending on when it ran, on both
trees. The first pair alternated on one venv gave branch 16.06 s and main 16.02 s; the second
pair gave branch 15.79 s and main 3.37 s.

The code the change adds is instead timed inside the same steps (`BeltTimer`: `mark_belt`,
`seed_belts`, `belt_mask`, `core_count` and `majority_history`, wherever they are called from):
it takes **0.09%** of step time at density 4 (two repeats, 0.093% each).

```bash
backend/.venv/bin/python bin/debug/measure_belt_marks.py bench --seed 7 --density 4 --repeats 2 \
  --out bench.jsonl
```

## Runtime

Each 1 Gyr seed took 40–97 minutes of wall time. Other sessions loaded the machine, so these
times are not comparable between runs.
