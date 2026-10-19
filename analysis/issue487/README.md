# RFC 0010 diagnosis gate 1: event-stream reconciliation

## Result

The gate-one probe now records a paired shadow-calculator result for every production
rift-cut event, including the numeric size, pole-separation, residual and connectivity
margins. A clean three-seed replay completed on branch commit
`0118752d0ffe75060464204d5b94f81a60fe9a07` with cratons-first layout, craton formation
off, 2 Myr steps through 1 Gyr, a 1,000 km corridor, cost floor 0.05, and
`PYTHONHASHSEED=0`. Python and library versions, clean worktree status, and SHA256 hashes
of the probe, shadow calculator and production source are recorded in
[`gate1-summary.json`](gate1-summary.json).

| Seed | Tap rows | Production `all_cuts` | Flow-through candidates | Admissible | Line-clear admissible |
|---:|---:|---:|---:|---:|---:|
| 7 | 144 | 144 | 44 | 17 | 6 |
| 11 | 111 | 111 | 42 | 31 | 12 |
| 2024 | 125 | 125 | 44 | 28 | 5 |
| **Pooled** | **380** | **380** | **130** | **76** | **23** |

Each seed's sequential tap rows exactly match its final production ledger cut count. The
pooled through-craton rates are `L=0.177` and `B=0.837`. The gate margins are summarized
over those 130 flow-through candidates in `gate1-summary.json`; individual failure counts
can overlap when an event fails more than one gate.

Two independent 100 Myr runs with the same fixed hash seed produced identical event rows.
A clean 100 Myr run at archived model commit `68f8bf3` also produced the same rows as
`origin/main` under the same settings. These are repeatable short-run baselines. The new
full run records exact source and runtime provenance for a longer baseline.

## Archived stream discrepancy

The archived production-probe and margin-instrumented series report the same git commit and
nominal seed/layout/step/duration settings, but their event signatures differ substantially
from each other and from the clean run. Signature identity is `(kind, elapsed Myr, plate ID,
node count)`. Only 3–6 events per seed match either archived series in the new full run.
The production series has 125/137/130 events for seeds 7/11/2024; the margin series has
130/116/115; the clean run has 144/111/125. Keep these trajectories and denominators
separate.

The archived JSONs omit working-tree fingerprints and runtime/library versions. Their
shared commit hash is not enough to determine whether uncommitted model or diagnostic code,
or a runtime difference, caused the divergence. The archived mismatch cannot be attributed
to a specific cause from the saved metadata. The current matched baseline is reproducible,
but the historical stream discrepancy remains unresolved; diagnosis gate 1 should remain
open until the old runs are either recovered with provenance or explicitly excluded from
calibration.

## Reproduction

Run from the repository root after `./bin/setup.sh`:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
  bin/debug/probe_craton_mosaic.py \
  --layout cratons-first --seeds 7,11,2024 --years 1e9 --step 2e6 \
  --every 1e8 --out /tmp/issue487-gate1-current-main-instrumented.json
```

The full 446 KB output is intentionally left in `/tmp`; the compact statistics,
fingerprints, and archived-series comparison are checked in as `gate1-summary.json`.
