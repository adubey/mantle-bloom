# Full-duration losses after lake transfers

Issue [#768](https://github.com/adubey/mantle-bloom/issues/768) addresses principle 3:
transfers settled after the timed pass must participate in the receiving lake's implicit
water balance. Principle 2 requires measuring the cost of revisiting affected bodies.

## Preregistration

Recorded before any replay. Baseline: detached `origin/main`, `dea1f6b5` (PR #761 merged).
Branch: the commit containing this preregistration, on `codex/issue-768-lake-transfer-losses`.

The implementation restores a body's previously booked slice losses before re-solving its
current storage and new transfers over the same slice duration. The new loss replaces the
old charge, including losses charged to children before a merger. External supply enters only
the timed pass. Every slice settles transfers and topology before the next slice begins.

**Equation check.** A three-cell fixture starts split, with two 0 m floors, a 4 m shelf,
unit cell areas and a 1 m joining saddle. Supply is 30 m³ over one year and combined loss
is 1 m/yr. The final fully flooded lake must have 27 m³, at level 31/3 m. On unmodified
`dea1f6b5` it holds 29 m³ at 11 m, while its ledger still closes. Evaporation, sublimation,
and their combination must each match the independent answer within 1e-8 m and m³.
Changing which sibling runs first must preserve that answer. An unchanged receiver getting
a late transfer must replace its loss rather than charge another full-duration loss.

**Windows.** Each arm runs 12 steps of 100 kyr (1.2 Myr), ice-age cycle off. Anchors and
checkpoints are those of #757; no alternative checkpoint will be substituted.

| Window | Seed | Opening age (Myr) | SHA-256 | Anchor |
|---|---:|---:|---|---|
| E | 910211954 | 224.0 | `4045b915a2d3d9e534a84089bba7277bbd488690e51d70f4e928f58c8cf674b1` | `analysis/rfc0016/target-anchor-window-e.json` |
| S | 317739385 | 201.6 | `ad8d34d60ac89db97efdfd722bc3db0287db22d8917288871502fe9a7747f94e` | `analysis/issue757/anchor-seed317739385-201.6.json` |

S is absent from the configured save directory and the B2 listing at registration; its
location has been requested. E was retrieved from B2 and its SHA-256 matches.

**Primary metric and decision.** `--audit-losses` captures the solved water column before
the ice split. For every wet final body, compare its trace's booked evaporation and
sublimation against the prescribed rates times final wet area times the full step duration.
The sum of absolute differences, divided by combined expected loss (minimum 1 m³), must be
at most 1e-9 on the branch. Report E's anchor at every step, especially the baseline's refill
steps 2 and 10. Dry bodies are excluded: a supply-limited body can consume water while ending
with zero wet area. The hand-calculated fixture remains the independent equation-level check.

**Accounting gates.** Lake and rain-only land normalized residuals must be <= 1e-9 every
step. Liquid and ice sub-ledgers must close. All lake merge, split, spill, capped-cycle and
silt-accounting tests must pass. A bounded-settling failure is a failure, not a usable result.

**Predictions, no magnitude threshold.** E's artificial refills should shrink because they
now face full-duration losses. The target may still swing through lagged climate and rain
partitioning; neither zero swings nor identical endpoints from different opening states is
required. Report whole-world liquid, ice and combined volume, lag-1 autocorrelation of volume
changes, mean absolute step change from step 3, sea level, largest lakes, runtime and peak RSS.
Do not infer stability from autocorrelation alone: 12 steps contain few independent cycles.
These are deterministic comparisons on two seeds, not an estimate of planetary variability.

**Controls.** The baseline should reproduce #757's final E branch (model commit `e5c697d2`):
step 2 target volume about 27.87 M km³, step 10 about 28.32 M km³, and final whole-world
liquid about 14.30 M km³. Compare exact values from its compact output where available;
use a 0.01% tolerance against unrounded data. Record the SHA and commit for every run.

**Runtime.** Follow `compute-budget` before each arm and run arms sequentially. Prefer no
other heavy run on the machine; disclose overlap. E previously peaked at 4.5 GB, S at 7.1 GB;
reserve 5 and 8 GB respectively. An increase above 25% in comparable median step runtime
triggers investigation before delivery; it is a diagnostic threshold, not proof of a regression.

```bash
backend/.venv/bin/python bin/debug/replay_lake_water_swing.py <checkpoint.mbworld> \
  --anchors <anchor.json> --myr 1.2 --step-years 100000 --checkpoint-myr 0.1 \
  --audit-losses --backend <arm-tree>/backend --out <captures>/issue768-<window>-<arm>.json
```

Full outputs remain outside the repository; only compact results and findings are committed.
