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

## Results

Branch commit `d2ea1455`, baseline `origin/main` at `9e2ae610`. Compact numbers in
`replay_ab.json` (checked in); full per-step JSONs stay in scratch
(`~/mantle-bloom-captures/issue-789/`). The replay arms ran before the branch
merged `origin/main` (`cba9875b`); the merge only brings enclosed-ocean hysteresis
(#790) and tests, no ice/march changes, so the acceptance conclusions stand as
recorded.

### Step time (acceptance: median within +15% of origin/main)

| Save | Branch med | Baseline med | Delta |
|---|---|---|---|
| seed 997271774 @ 87.4 Myr | 22 s | 7 s | +213% |
| seed 649658380 @ 6.4 Myr | 15 s | 6 s | +137% |
| seed 374464662 @ 132.6 Myr | 13 s | 8 s | +76% |

OVER budget on all three saves, as preregistered. A settled march costs ~13-32 ice
solves against one (up to 520 substeps on one 997 step: 383 s); ice was ~9% of a step,
so the march alone explains +76-213%. Owner budget decision pending. Peak RSS is
unchanged (~1 GB both arms).

### Ledger closure and settlement (acceptance: land <= 0.1% every step; march settles)

| Save | Branch land max | Baseline land max | Branch ice max | Settled |
|---|---|---|---|---|
| 997271774 | 1.0e-13 | 3.8e-12 | 3.8e-03 | 10/10 |
| 649658380 | 3.9e-13 | 1.6e-12 | 1.9e-02 | 10/10 |
| 374464662 | 3.2e-13 | 6.7e-14 | 2.1e-02 | 10/10 |

Land closes to round-off on both arms: PASS with ~9 orders of margin. Branch ice
residual is 0.4-2% of accumulation on growing worlds (the hold extrapolation net,
kept as residual by design, not as ice). Baseline ice closes to 1e-13 (scaled units,
near-frozen). Per the RFC acceptance the binding metric is the land ledger; ice
residuals are reported, not pass/fail here. No before/after ice-volume comparison
(physical units move volumes by design: e.g. 997 closes at 12.8M km3 branch vs
5.8M km3 baseline).

### Full-suite triage: two failures, both branch-only, both 1 Myr steps

`test_erosion.py::test_deposition_fills_channels_back_in` and
`test_faults.py::test_fault_relief_stamps_a_fault_reason_code_on_nearby_crust` fail
on the branch (isolated and in-suite) and pass on baseline (isolated). Both step
1 Myr worlds; everything else in the suite passes.

**Erosion: fixed in this workstream.** The test's premise ("a fresh world has no ice,
so wherever no river erodes nothing carves") breaks under D1: the 1 Myr step grows
real ice (depth sum 4.7e5 -> 2.9e6, 2423 nodes over 1 m, max 3.5 km) whose scour
deepens `channel_depth` (the mismatch node carries 301 m of ice). 146/375 quiet
nodes have ice; the 229 ice-free quiet nodes satisfy the exact fill relation, so the
fluvial-deposition contract the test pins is intact. Fix: mask `glacier_depth > 0`
nodes out of `quiet` (glacial carving has its own test). Now passes.

**Fault: branch-only lake-flush crash, root-caused to a lake cap-cycle limit under
physical melt throughflow; fixed lake-side in #794 (PR #800), now green.** At the
trip step (seed 7, 9th 1 Myr step) the lake solver gets melt 2.49e13 m3/yr on the
branch vs 4.4e8 on baseline (57,000x = the removed `MODEL_UNIT_SCALE`); dry runoff
is unchanged (2.65e13 vs 2.80e13). The raise was `Lake flush left an unconsumed
cycle owner at node 769`: a capped 23-node/6-node sibling pair sloshed the
throughflow, `held.volume += overflow` lapping it ~200 times to a phantom 2.25e17
m3, while the LCA parent-outlet owner (which had room, completed handoff) never
absorbed it. The melt itself is legitimate equilibrium throughflow (surface PDD +
basal + bare-ground snow; guard-shed overflow/capped are zero on the probed
analog), so this was the legacy lake routing meeting D1's regime, not a D1
accounting bug. The lake-side fix (reuse a parent-outlet owner with room when
re-routing held cycle surplus) landed via #800; the fault test passes on the
branch with it merged in (verified post-merge).

### Open items

- Step-time budget: OVER (+76-213%); owner decision with the numbers above.
- Fault test: green since #794/#800 (was the only red test in the suite).
- F0D replay cannot run (held-out, absent locally and in B2); acceptance item open
  pending the owner's call.
