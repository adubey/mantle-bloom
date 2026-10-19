# RFC 0010 diagnosis gate 1: event-stream reconciliation

## Result so far

A clean replay of the recorded seeds was run on `origin/main` commit
`fcaaf36e95463c078e6441fbee0125da4e275424`, with cratons-first layout, craton formation
off, 2 Myr steps to 1 Gyr, a 1,000 km corridor, cost floor 0.05, and
`PYTHONHASHSEED=0`. The full event output remains in `/tmp/issue487-gate1-current-main-hash0.json`.
The compact results and signature hashes are in [`gate1-summary.json`](gate1-summary.json).

The probe taps `cratons.record_rift_cut`, evaluates the shadow min-cut for that same event,
and then calls the original production recorder. Thus each output row pairs a production
cut event with its shadow result. For seeds 7, 11 and 2024, the row counts (144, 111 and
125) exactly equal the final production ledger cut counts. The output also retains candidate
counts, admissibility decisions, line clearance, band exposure, and residual ratios.

Two independent 100 Myr runs with the same fixed hash seed produced identical event rows.
A clean 100 Myr run at the archived model commit `68f8bf3` also produced the same rows as
current `origin/main` under the same settings. This is a repeatable short-run baseline.

## Why the archived full runs are not combined

The archived production-probe and margin-instrumented series report the same git commit and
nominal seed/layout/step/duration settings, but their event signatures differ substantially
from each other and from the clean run. Signature identity is `(kind, elapsed Myr, plate ID,
node count)`. Only 3–6 events per seed match either archived series in the new 1 Gyr run.
The production series has 125/137/130 events for seeds 7/11/2024; the margin series has
130/116/115; the clean run has 144/111/125. Since these are different trajectories, their
counts and denominators must remain separate.

The archived JSONs omit working-tree fingerprints and runtime/library versions. Their shared
commit hash is not enough to determine whether uncommitted model or diagnostic code, or a
runtime difference, caused the divergence. The current clean short run is repeatable, but the
specific cause of the old full-run mismatch is not recoverable from the saved metadata.

## Remaining gate-1 work

The first full run's compact per-event output records the residual ratio and pass/fail state
for size, pole separation, residual and connectivity gates, but drops the numeric minimum-half
node count and pole-separation margins. Those values are needed to compare the full detailed
margin series on a matched event stream. Add them to the probe output and rerun the baseline
before closing diagnosis gate 1. Until then, keep gate 1 open and do not use the archived
series as a calibration denominator.

The probe has now been extended to retain those four gate margins, assign a sequential event
index, and record the Python/library versions, hash seed, working-tree status and hashes for
the diagnostic and production cut sources. A 10 Myr smoke run produced four tap rows and four
production-ledger cuts, with numeric size, pole, residual and connectivity margins on each
row. The full three-seed replay with this enriched output is still required before the gate
can close.

Reproduction command for the current full run:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
  bin/debug/probe_craton_mosaic.py \
  --layout cratons-first --seeds 7,11,2024 --years 1e9 --step 2e6 \
  --every 1e8 --out /tmp/issue487-gate1-current-main-hash0.json
```
