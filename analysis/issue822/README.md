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
