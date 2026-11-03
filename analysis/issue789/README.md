# Issue #789: RFC 0015 D1-D2 implementation replay

## Preregistration

Written before any replay run in this report (replay script and this preregistration
committed first). Branch `muse/issue-789-d1-d2-march` ports the G0/G1 gate prototypes to
production: physical accumulation/surface/basal melt in m/yr of ice inside `rates()`
(D1, G2 knobs DDF_snow 4, DDF_ice 8, T_rs 1.0 C) and the 250 yr settle-and-hold march (D2).

### Tree and saves

Branch at the prereg commit, plus a task-specific detached `origin/main` baseline worktree
for the step-time comparison. Ten 0.25 Myr steps from each diagnosis save (draw 0, each
world's own noise):

| Save | Ice ages |
|---|---|
| Seed 997271774 at 87.4 Myr | off |
| Seed 649658380 at 6.4 Myr | off |
| Seed 374464662 at 132.6 Myr | 8 Myr period |

Held-out saves (F0D, S200.8, S100.2) are not on this machine and are never used, per the
issue constraints -- so the 21.5 Myr F0D replay in the acceptance list cannot be run from
this workstream and is flagged for the owner (see below).

### What the replay measures, per step

`bin/debug/replay_789_march.py` records wall time, the land ledger's residual against its
external inputs (rain + meltwater), the ice ledger's residual against its accumulation,
the march log entry (substeps, marched/held years, settlement, volume path -- null on
the origin/main arm, which has no march), and closing ice volume and sea level. The
baseline arm runs the same script with the baseline worktree's venv and
`MANTLE_BLOOM_BACKEND=<baseline-worktree>/backend`.

### Expected direction and size

- **Ledger closure (acceptance, hard):** every step closes the land ledger to 0.1% of its
  external inputs. Expected PASS: the march sums whole-substep ledgers by construction
  and the unit-level contracts already close to round-off.
- **March settles (acceptance):** settles on the diagnosis replays, with the MedD ice-age
  surcharge noted if the ice-age save settles less often. Expected: settles on most
  steps (G1 measured 119/120 at the 250 yr cap on these saves).
- **Step time (acceptance, owner decision):** median step time within +15% of origin/main
  on the same machine, runs one at a time. Expected: OVER budget (a settled march costs
  ~13-16 ice solves against one; ice was ~9% of a step), reported with the numbers for
  the owner's budget decision. No variability claim: physical units move ice volumes by
  orders of magnitude, so before/after ice volumes are not compared.
- **F0D:** cannot run (save absent, held-out). The issue's F0D acceptance item stays open
  pending the owner's call.

### Controls

Same machine, arms one after the other per `compute-budget`. Step-time arms run with no
other heavy runs listed. Save SHAs recorded at run time.
