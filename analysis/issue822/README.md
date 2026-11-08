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

Pending.
