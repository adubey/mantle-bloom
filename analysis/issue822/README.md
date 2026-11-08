# Issue #822: cheaper end-of-step glacier march, initial option comparison

## Preregistration

Written before the comparison runs below (one timing-calibration run of the `d2` arm on
the generated seed-5 world came first and is reported as such). Branch
`claude/issue-822-cheaper-march` implements three of the issue's four options as
`ice_sheet.MarchCost` knobs, selected by preset; the CI-split option is out of scope.
The production default (`d2`) is the fixed-250-yr D2 march unchanged.

| Preset | Option | What it changes |
|---|---|---|
| `d2` | — | Reference: the production march. |
| `cap-hold` | 1 | At most 200 substeps, then the settled-hold fill (window-mean flows, ice frozen). |
| `cap-endpoint` | 1 | At most 200 substeps, then one backward-Euler substep over the whole remainder. |
| `adaptive` | 2 | Substep doubles after each pair whose pair-mean change is < 1%, up to 4 kyr; back to 250 yr above 5%. |
| `cheap-substep` | 3 | Skip the sea-level re-solve while the debit moved < 1 cm of sea level; exact local inside-distance search. |
| `cheap-substep-limited` | 3 | As above, plus the outside-distance search capped at 1000 km (fallback when a member is cut off). |

A per-node linear-trend cap fill was prototyped first and dropped before this
comparison. On the flow-line unit fixture, which keeps growing and never settles, it
stacks ice in place instead of advancing the margin, and the 5 km depth cap then breaks
the ice ledger (78–95% residual).

### Worlds and steps

Every arm of a world opens from the same bytes (an exact save; generated worlds are
generated once, then saved), so accuracy differences are deterministic. Metrics compare
step 0 of each arm against step 0 of `d2`, from the same opening state. Later steps open
from states that have already diverged; for those only cost and ledger closure are compared.

| World | Steps | Why |
|---|---|---|
| `generate_world(5, node_density=0.25)` | 1 x 0.25 Myr, and separately 1 x 1 Myr | The issue's craton probe world: ~1000 substeps, never settles. |
| `generate_world(30, num_plates=8)` | 1 x 0.25 Myr | The issue's monster test world (19m41s for 1 Myr). |
| Save seed 997271774 @ 87.4 Myr | 4 x 0.25 Myr | #789 diagnosis save, mature world. |
| Save seed 649658380 @ 6.4 Myr | 3 x 0.25 Myr | #789 diagnosis save, mostly settles quickly. |

### Metrics and expected direction

- **Cost (primary):** march wall time (the march alone, timed around `_march_glacier`)
  and whole-step wall time, as a ratio to `d2`.
- **Accuracy at step 0, against `d2`:** closing ice volume (relative difference), node-wise
  glacier depth L1 difference over `d2`'s total, sea level after the step (m).
- **Ledgers:** ice-ledger closure |residual| / accumulation, and land-ledger closure
  (the #789 criterion is land <= 0.1% every step).

Expected:

- `cap-hold`: march cost ~0.2x on never-settling steps (200 of ~1000 substeps), unchanged
  where the march settles before 200. Closing ice within ~1% where the sheet oscillates
  near balance. Ice closure degrades by the window drift times the held years, which
  could reach several percent.
- `cap-endpoint`: about the same cost as `cap-hold`, and the ice ledger closes to
  round-off. The closing ice may jump toward the backward-Euler equilibrium: larger
  differences from `d2` than `cap-hold` on still-moving sheets.
- `adaptive`: march cost 0.1–0.5x on never-settling near-balance steps, since pair changes
  of ~5e-4 sit far below the 1% growth threshold. Closing ice within ~1% of `d2` unless
  the alternation grows at long substeps (the #741 10 kyr drops; the 4 kyr cap is
  meant to stay clear of them). Ice closure at round-off.
- `cheap-substep`: march cost 0.75–0.9x (sea-level re-solves ~17% of a march, plus the
  inside-distance dijkstra). Closing ice within 0.1% of `d2`, and sea level within a few cm.
- `cheap-substep-limited`: a few percent cheaper again, accuracy as `cheap-substep`.

### Decision threshold ("promising")

An option is promising for follow-up if, on the never-settling generated worlds, its march
costs <= 0.5x `d2`, its step-0 closing ice is within 1% of `d2`, its node-wise L1 depth
difference is <= 5%, its ice closure is <= 1e-3, and its land closure stays <= 0.1% on
every step of every world. Option 3 cannot reach 0.5x by design; it is judged on
whether its saving is free (accuracy within 0.1%), as a multiplier for whichever of
options 1 and 2 goes forward.

### Controls and variability

One process per arm, arms run one after another and alone (compute-budget), on one 8 GB
Apple Silicon machine. Accuracy metrics are deterministic for a given start state and
code. Wall time varies with machine load; a ratio within ~10% of 1 counts as no change.
The `d2` preset runs exactly the D2 code path, and is checked bit-identical to the
pre-change march by the existing unit tests.

## Results

Branch commit `45bcd1f6` (the code is `6dc9a483` plus a docs commit). All 30 arms ran one
at a time on an 8 GB M-series machine with nothing else heavy running; peak RSS was
0.4–1.35 GB. Compact numbers are in `compare.json`. Per-arm JSONs, depth fields and the
generated start saves are kept out of the repo in the issue's capture directory, with
`SHA256SUMS`. The timing-calibration `d2` run on gen5 (28.2 s step, 1000 substeps) was
reproduced exactly by the batch's `d2` arm (28.6 s, same closing ice).

Commands (placeholders for local paths):

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_822_march_cost.py \
    --generate '{"seed": 5, "node_density": 0.25}' --write-save <scratch>/gen5-d025.mbworld
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_822_march_cost.py \
    <save> --preset <preset> --steps <N> --myr <M> --out <scratch>/<world>-<preset>.json
backend/.venv/bin/python bin/debug/summarize_822_march_cost.py <scratch> \
    --worlds gen5-250k gen5-1m gen30-250k s649 s997 --json analysis/issue822/compare.json
```

| Start state | SHA-256 |
|---|---|
| gen5 (`generate_world(5, node_density=0.25)`, exact save) | `4f933c89...` |
| gen30 (`generate_world(30, num_plates=8)`, exact save) | `7edae919...` |
| seed 649658380 @ 6.4 Myr | `ebbfb5da...` |
| seed 997271774 @ 87.4 Myr | `08ad9461...` (as in `analysis/issue789/replay_ab.json`) |

### Never-settling generated worlds (step 0, same opening state as `d2`)

March x = march wall time / `d2`'s. Ice and sea level are differences from `d2`; L1 is
node-wise |depth − depth_d2| summed, over `d2`'s total.

| World | Preset | March x | Ice vs d2 | Depth L1 | Sea level | Ice closure |
|---|---|---|---|---|---|---|
| gen5, 0.25 Myr (d2: 1000 substeps, 27 s) | cap-hold | 0.20 | −0.01% | 0.01% | +0.00 m | 5.6e-4 |
| | cap-endpoint | 0.20 | **−54.7%** | 53% | +9.70 m | 1.6e-16 |
| | adaptive | 0.17 | +3.1% | 4.4% | −0.55 m | 7.5e-17 |
| | cheap-substep | 0.91 | 0.00% | 0.00% | 0.00 m | 3.8e-15 |
| gen5, 1 Myr (d2: 4000 substeps, 109 s) | cap-hold | **0.05** | 0.00% | 0.00% | 0.00 m | 1.2e-2 |
| | cap-endpoint | 0.05 | **+176%** | 180% | −32.2 m | 8.4e-17 |
| | adaptive | 0.09 | +4.3% | 4.0% | −0.77 m | 1.8e-15 |
| | cheap-substep | 0.92 | 0.00% | 0.00% | 0.00 m | 9.3e-14 |
| gen30, 0.25 Myr (d2: 1000 substeps, 131 s) | cap-hold | 0.19 | +0.66% | 7.8% | −1.06 m | **3.6e-2** |
| | cap-endpoint | 0.20 | **+49.7%** | 63% | −89.5 m | 2.5e-16 |
| | adaptive | 0.15 | +4.9% | 15% | −8.92 m | 2.7e-16 |
| | cheap-substep | 1.00 | +0.21% | 5.3% | −0.36 m | 5.5e-16 |

`cheap-substep-limited` matched `cheap-substep` to the digit on every world (time within
1%), so the 1000 km outside-distance limit buys nothing here and is not listed.

### Diagnosis saves (multi-step; land closure <= 3.5e-13 on every step of every arm)

The cap presets are identical to `d2` until their first capped step, so that step opens
from the same state and is compared directly. `adaptive` and `cheap-substep` diverge
from step 0, so only their step 0 and total cost are like-for-like.

| Save | Preset | March x (all steps) | Compared step | Ice vs d2 | Depth L1 | Sea level | Worst ice closure |
|---|---|---|---|---|---|---|---|
| 649658380 @ 6.4 Myr, 3 x 0.25 Myr (d2: 30/1000/18 substeps, 427 s) | cap-hold | 0.24 | 1 (capped) | 0.00% | 0.00% | 0.00 m | 7.1e-3 |
| | cap-endpoint | 0.24 | 1 (capped) | +0.47% | 10.7% | −0.03 m | 4.7e-3 |
| | adaptive | 0.19 | 0 | +0.45% | 2.5% | −0.04 m | 1.1e-4 |
| | cheap-substep | 0.84 | 0 | 0.00% | 0.00% | 0.00 m | 4.4e-3 |
| 997271774 @ 87.4 Myr, 4 x 0.25 Myr (d2: 34/14/70/298 substeps, all settle, 203 s) | cap-hold | 0.77 | 3 (capped) | +0.82% | 3.2% | −0.30 m | 1.2e-2 |
| | cap-endpoint | 0.77 | 3 (capped) | +8.4% | 14.9% | −3.12 m | 1.7e-3 |
| | adaptive | 0.64 | 0 | +1.0% | 2.8% | −0.42 m | 8.9e-5 |
| | cheap-substep | 0.99 | 0 | 0.00% | 0.00% | 0.00 m | 1.7e-3 |

`d2`'s own worst ice closure on these saves is 4.4e-3 (649) and 1.7e-3 (997), from its
settled hold. That is already over the preregistered 1e-3 ice threshold, so on the saves
ice closure is read relative to `d2`.

### Why the march never settles: a period-5 limit cycle

(Round 2 found the 0.25 Myr step's cycle is period 30, and 649 step 1's is period 20:
the period varies by step. See round 2 below.)

`cap-hold` at 1 Myr on gen5 closes on exactly `d2`'s ice (node-wise max difference
8e-7 m) after 200 substeps instead of 4000. A follow-up run recording every substep's
closing volume on that step shows why. From somewhere between substep 100 and 200, the
march sits on an exact period-5 cycle (rtol 1e-9) with 0.4% volume amplitude. Pair
means of an odd-period cycle never coincide, so the pair-change test straddles its 1e-3
tolerance forever: the issue's "last_change = 5.65e-4, never settles". `cap-hold` and
`d2` both land on the same point of the cycle (substeps 200 and 4000 differ by 3800, a
multiple of 5). The same coincidence shows on the 649 save's capped step. The march
reached its balance long before the cap; only the settle test missed it.

### Verdict against the preregistered threshold

| Option | Preset | Cost (never-settling) | Accuracy vs d2 | Ledger | Verdict |
|---|---|---|---|---|---|
| 1 | cap-hold | 0.05–0.20x | 0.0–0.8% ice; L1 0–7.8% | ice closure 0.06–3.6% (window drift over the held years) | **Promising, with a caveat**: accurate wherever the sheet is on its cycle; the ledger fails 1e-3 when the cap cuts a still-drifting sheet (gen30, 997 step 3) |
| 1 | cap-endpoint | 0.05–0.20x | −55% to +176% ice | closes (round-off) | **Rejected**: one long backward-Euler solve jumps branches (#741's 10 kyr drops, larger) |
| 2 | adaptive | 0.09–0.17x | +0.45% to +4.9% ice, always high; L1 2.5–15% | closes (round-off) | **Not yet**: fastest and best ledger, but a systematic +3–5% ice bias on never-settling worlds misses the 1% bar; it also costs more than `d2` on steps that settle fast (997 steps 0–1) |
| 3 | cheap-substep | 0.84–1.00x | 0.00% on four worlds; +0.21% / L1 5.3% on gen30 | unchanged | **Small, not free**: 0–16% saving; the 1 cm sea-level skip perturbs the never-settling gen30 trajectory. The exact local dijkstra is free |
| 3 | cheap-substep-limited | as cheap-substep | as cheap-substep | unchanged | **No gain** over cheap-substep |

### Follow-ups suggested by these numbers (not implemented here)

- **Settle on recurrence, not on pair change.** The never-settling marches are on an exact
  short cycle, so a test for state recurrence (period p <= ~10) followed by a hold over
  one whole period's mean flows would settle them exactly, with a closed ledger. That is
  the issue comment's "oscillation-vs-drift detection". It would make `cap-hold`'s cap a
  rarely-hit safety net rather than the main saving.
- **Adaptive with a smaller growth cap** (e.g. 1 kyr) to test whether the +3–5% bias is
  the long-substep branch effect, and with growth only while unsettled for >N pairs so
  fast-settling steps keep D2's cost.
- **Option 3:** keep `local_dijkstra` (exact). Drop `dijkstra_limit_m`. The sea-level
  skip is worth re-measuring only once the substep count is fixed, since its saving
  (<= 16%) is small next to options 1 and 2.

## Round 2: settle tests that see the cycle

### Preregistration (round 2)

Written before the round-2 runs. One diagnostic came first: a `d2` run of gen5's
0.25 Myr step recording every substep's depth field. Its depths repeat **exactly**
(L1 difference 0) every **30** substeps from about substep 141 on. Shorter lags only
nearly repeat: lag 6 differs by 1.8e-4 and lag 5 by 1.8e-2. The period-5 cycle reported
above was in the 1 Myr step's volume, so the cycle period changes from step to step.
A 10-substep window mean leaves a 2e-3 change on this step; a 30-substep one leaves 0.

New presets (code in `ice_sheet.MarchCost`, `hydrology._SettleWatch`). Both new tests
run beside D2's pair test, not instead of it:

| Preset | What it changes |
|---|---|
| `cap-hold-recur` | `cap-hold` plus a recurrence settle test: the march settles when its depth field matches the one p substeps back (L1 < 1e-6 of the ice, p <= 64) three substeps running, then holds on whole periods' flows. |
| `trend-hold` | The owner's proposal, hold only: every 30 substeps, the L1 change between consecutive 30-substep window-mean depth fields. The change still to come is estimated as D r/(1-r) while the ratio r of successive changes is < 1, else D times the windows left. The march holds on the last window's flows once that estimate is < 1e-3. |
| `trend-adapt` | The owner's proposal: as `trend-hold`, and while the change shrinks (r < 1) and the march hasn't held, the substep doubles, up to 4 kyr; r > 1.5 drops it back to 250 yr. |
| `trend-adapt` + `trend_max_substep_years=1000` | Tests whether `adaptive`'s ice bias comes from the long substeps. |
| `adaptive` + `adaptive_max_substep_years=1000` | Same test on round 1's adaptive rule. |

`d2`, `cap-hold` and `adaptive` are not rerun: their code paths are unchanged. A `d2`
control rerun on gen5's 0.25 Myr step checks that, and must match round 1's closing
depths bit for bit. Worlds, steps, metrics and the "promising" threshold are as in
round 1.

Expected:

- `cap-hold-recur`: on gen5 (both step lengths) it settles by recurrence before the cap,
  at about 0.15–0.2x march cost. Its closing ice differs from `d2` by at most the
  cycle's swing (0.6% on the 0.25 Myr step), because the march stops at a different
  point of the cycle. Ice closure drops to <= 1e-4, since a whole cycle's flows net to
  zero. Where no exact cycle forms within 200 substeps (possibly gen30 and the 997 save's
  step 3), it falls back to the cap and matches `cap-hold`.
- `trend-hold`: about the cost of `cap-hold-recur` wherever the cycle period divides 30
  (5 and 30 both do). Elsewhere the windowed change plateaus above 1e-3 and the march
  runs to the end of the step, like `d2`. Ledger as `cap-hold-recur`.
- `trend-adapt`: cheaper than `trend-hold` if the change shrinks early, before the
  cycle. Growth may push it onto a different cycle or bias it like `adaptive`. Expect
  its ice error between `trend-hold`'s and `adaptive`'s, and the ledger to close.
- The 1 kyr ceilings: if the long substeps cause the bias, both 1 kyr arms land closer
  to `d2` than their 4 kyr versions, at somewhat higher cost.

### Results (round 2)

The `d2` control rerun on gen5's 0.25 Myr step matched round 1's closing depths
exactly (max difference 0.0 m), so round 1's `d2`, `cap-hold` and `adaptive` arms are
reused unchanged. Full per-world tables: `compare_round2.json` (made by
`summarize_822_march_cost.py --presets ...`). Land closure is <= 1.3e-12 on every step
of every arm.

Each long march is compared with `d2` from the **same opening state**. That is step 0
on the generated worlds. On the saves it is the step with the long march (649 step 1,
997 step 3); every arm except the two `adaptive` ones opens that step from `d2`'s exact
depths, because earlier steps all settled the same way. The `adaptive` arms open
diverged and are marked \*. Cost is the march wall time over the whole run, as a
multiple of `d2`'s.

Cost, closing ice vs `d2`, node-wise L1 depth difference, and ice closure:

| Arm | gen5 0.25 Myr | gen5 1 Myr | 649 step 1 | 997 step 3 | gen30 0.25 Myr |
|---|---|---|---|---|---|
| `cap-hold` | 0.20x, −0.01%, 0.0%, 5.6e-4 | 0.05x, 0.00%, 0.0%, 1.2e-2 | 0.24x, 0.00%, 0.0%, 7.1e-3 | 0.77x, +0.82%, 3.2%, 1.2e-2 | 0.19x, +0.66%, 7.8%, 3.6e-2 |
| `adaptive` | 0.17x, +3.08%, 4.4%, 0 | 0.09x, +4.29%, 4.0%, 0 | 0.19x, +1.12%\*, 13.5%\*, 0 | 0.64x, +0.45%\*, 6.1%\*, 0 | 0.15x, +4.91%, 14.6%, 0 |
| `adaptive` 1 kyr | 0.29x, −0.87%, 2.5%, 0 | 0.02x, +0.59%, 2.2%, 9.9e-5 | 0.31x, +2.08%\*, 7.2%\*, 0 | 0.56x, +0.32%\*, 4.0%\*, 6.9e-5 | 0.30x, +7.12%, 10.8%, 0 |
| `cap-hold-recur` | 0.17x, +0.44%, 1.9%, 4.5e-7 (period 30) | 0.04x, −0.32%, 1.9%, 5.2e-6 (period 5) | 0.16x, −0.15%, 5.8%, 1.0e-6 (period 20) | = `cap-hold` (no cycle) | = `cap-hold` (no cycle) |
| `trend-hold` | 0.14x, +0.46%, 1.9%, 1.4e-5 | 0.04x, 0.00%, 0.0%, 3.1e-5 | 0.13x, −0.01%, 0.1%, 1.6e-4 | 1.00x, = `d2` (never holds) | 1.00x, = `d2` (never holds) |
| `trend-adapt` | 0.21x, −0.87%, 1.8%, 2.7e-6 | 0.04x, +0.03%, 0.5%, 5.3e-4 | = `trend-hold` | 0.61x, +0.23%, 2.4%, 7.6e-4 | 0.22x, −0.15%, 14.5%, 0 |
| `trend-adapt` 1 kyr | 0.24x, −0.81%, 2.5%, 0 | = `trend-adapt` | = `trend-hold` | = `trend-adapt` | 0.21x, +7.87%, 11.0%, 2.3e-3 |
| `d2` at 125 yr (control, `--substep-years 125`) | — | — | — | — | 0.26x, **−7.34%, 10.0%**, 1.8e-2 (pair-settled) |

`d2`'s own ice closure on the saves' long steps is 8.6e-4 (997 step 3); the
"0" entries are below 1e-15.

What the numbers say:

- **Exact cycles are common, and their periods vary:** 30 and 5 on gen5, 20 on 649
  step 1. gen30 and 997 step 3 show no exact cycle within 200 substeps. On those two the
  sheet is still moving: `cap-hold`'s ledger residual there comes from that drift.
- **Recurrence** settles every cycle before the cap, with the ledger closed to
  ~1e-6, as predicted. The ice differs from `d2` by up to 0.44% because the march
  stops at a different point of the cycle than `d2`'s last substep. That is the cycle's
  swing, not an error, and the hold can't remove it. With no cycle it is `cap-hold`
  exactly, ledger failure included.
- **`trend-hold`** is the most accurate early stop on these worlds. Where it holds,
  ice is within 0.46% (L1 <= 1.9%) and closure <= 1.6e-4. It is also the cheapest on
  cycles (0.04–0.14x), because a 30-substep window cancels periods 5 and 30. Period 20
  doesn't divide 30, yet the windowed change still fell below tolerance on 649. On a
  drifting sheet the windowed change stays at ~1.5e-2 and it never holds: no saving,
  but no harm either.
- **`trend-adapt`** saves where nothing else does without breaking the ledger: 0.61x
  on 997 step 3 (L1 2.4%) and 0.22x on gen30. On cycles it either holds as
  `trend-hold` does or grows first and lands a bit further off (−0.87% on gen5 0.25 Myr).
- **The 1 kyr ceiling hypothesis is not supported.** `adaptive` at 1 kyr is closer on
  gen5 (+3.1% → −0.9%, +4.3% → +0.6%) but further on gen30 (+4.9% → +7.1%) and 649;
  `trend-adapt` at 1 kyr is no better. The error isn't a simple function of substep length.
- **gen30 can't rank accuracy.** Halving `d2`'s own substep moves gen30's closing ice
  by −7.3% (L1 10%), as far as any arm. On a sheet that is still changing, 250 yr is not
  a converged reference. Every arm's 7–15% L1 there is within `d2`'s own substep
  sensitivity.

### Verdict (round 2)

| Arm | Cycling sheets (gen5 ×2, 649) | Drifting sheets (997 step 3, gen30) | Verdict |
|---|---|---|---|
| `cap-hold-recur` | 0.04–0.17x, ice ≤ 0.44%, ledger ~1e-6 | Falls back to `cap-hold` (ledger 1.2–3.6%) | **Promising**: strictly better than `cap-hold`; the cap is still the weak part |
| `trend-hold` | 0.04–0.14x, ice ≤ 0.46%, L1 ≤ 1.9%, ledger ≤ 1.6e-4 | No saving (1.00x), exact | **Promising, safe**: passes the threshold on every cycling sheet, never worse than `d2` |
| `trend-adapt` | 0.04–0.21x, ice ≤ 0.87% | 0.22–0.61x, ledger closed, L1 2.4–14.5% | **Promising with a caveat**: the only arm that saves on drifting sheets with a closed ledger; its gen30 field error can't be judged against an unconverged `d2` |
| `adaptive` / `adaptive` 1 kyr | 0.02–0.31x, ice +0.6% to +4.3% | 0.15–0.64x | **Not yet**: the bias moves with the ceiling but doesn't go away |

### Follow-ups suggested by round 2 (not implemented here)

- **Combine:** the trend test's hold plus its substep growth (`trend-adapt`), with
  recurrence as a second settle test. That covers cycles of any period cheaply and
  gives drifting sheets a ledger-closed saving, with no cap and no ledger failure.
- **Reference for drifting sheets:** a 125 yr / 62.5 yr `d2` convergence study on
  gen30 and 997 step 3, so drifting-sheet accuracy can be judged at all.
