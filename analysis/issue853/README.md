# Issue #853: why the lake solver's one-step result depends on slice length

Status: exploration complete. No model change was made under this issue, as its
deliverable requires. Probe design was preregistered on the issue before any run.

Driver: [`bin/debug/measure_lake_slice_convergence.py`](../../bin/debug/measure_lake_slice_convergence.py).
Compact data: [`lake-slice-convergence.json`](lake-slice-convergence.json). Raw outputs are in
`~/mantle-bloom-captures/issue853/` with the run logs.

## Summary in plain terms

The question was whether the lake solver's answer after one step depends on how that step is
sliced internally, whether that dependence survives once the slices are short enough, and
whether the RFC 0016 windows D and E "merged inland sea" state is an artefact of the slice
length. The governing principle is principle 3 in
[`docs/design-principles.md`](../../docs/design-principles.md): a long step must represent fast
processes in aggregate, and one long step should be compared with several shorter steps
against a stated tolerance. The metric is the anchor lake's volume after one step against
slice length, at a 1% convergence tolerance set before running.

Three findings:

- **The slice dependence is real on current code, much smaller than on the pinned source
  #811 measured, and it converges.** From window D's checkpoint, one 5 kyr step gives an
  anchor volume 13% low when solved in a single slice; that error falls below 1% at slices
  of about 50 years. At 25 kyr and 100 kyr totals the single-slice error is about 4% and
  about 1%, and slices of 1,000 years or less are converged to within 0.4% or better.
- **The mechanism is the loss terms, not merging or spilling.** Shortening the slices changes
  evaporation and sublimation, and overflow moves only to compensate. The solver charges each
  body its evaporation and sublimation at the water level it reaches at the end of the slice,
  for the whole slice. A lake that is filling therefore pays the full slice's losses at its
  final, largest area even for the part of the slice when it was smaller. Longer slices
  systematically overcharge a filling lake. Merge and split counts barely move with slice
  length, so the discrete merge/spill events that the preregistration predicted would carry
  the effect do not.
- **The merged inland-sea state does not survive.** On current code the window D and E anchor
  lakes do not reach the merged state at any slice length tested, converged or not, and the
  refreshed-climate trajectories at the converged 1,000 year slice never enter it either.
  This agrees with #811, which found the registered D/E events no longer occur on current
  code at all. The merged state on the pinned source was entered only at 2.5 and 5 kyr step
  sizes and was already absent at 1 kyr steps there.

**Recommendation: a bounded implementation issue, not an RFC and not "no change".** Give the
production lake solve an internal slice cap (1,000 years is converged in these windows; the
ice march already marches at 250 years and the river reservoirs at their own substeps), with
a partition test in the style of the existing lake partition gate. The science and the solver
structure do not change; this is a numerical-integration fix. Suggested follow-up issue text
is in the last section.

## Setup and validity

- Source: this branch at `origin/main` commit `5912a7d8`. The #811 numbers were measured on
  the pinned source `0bda513d`; absolute volumes differ between sources, so every comparison
  here is slice against slice on one source, one checkpoint.
- Checkpoints: window D, seed 910211954 at 225.1 Myr
  (`seed910211954-225.1myr.bloom`, SHA-256 `c5b1f5d8…`), and window E at 224.0 Myr
  (`seed910211954-224myr.bloom`, SHA-256 `95dde684…`), from the private captures under
  `rfc0016/window-rescan/checkpoints/`. Anchors are the registered files
  `analysis/rfc0016/target-anchor-window-{d,e}.json`.
- The anchor is a location, not a lake identity. Following the RFC 0016 driver
  (`replay_rfc0016_river_water.py`), the reported body is the wet component at the anchor, or
  the nearest wet component when the anchor cell is dry, tracked across snapshots by largest
  node overlap. On current code the anchor cell is dry after the first step in both windows,
  so the reported volumes are for that nearest tracked body; catchment areas confirm the same
  body is followed in every arm (0.82M km² for D, 1.47M km² for E).
- "Merged inland-sea state" uses #811's definition: the anchor body's upstream catchment
  exceeds 15M km² (the merged state measured 19–25M km²; the small-basin state 3–8M km²).
- Fixed-climate arm: one `step_world` of the stated total, with the lake solver forced by
  monkeypatch in the diagnostic script to `substeps = total / slice`. Climate, routing and
  all other processes see the same single step. Climate-refreshed arm: `step_world` repeated
  in steps of the slice length, so climate and routing are recomputed every slice.
- Validity: the lake ledger's normalized residual stays at or below 5.1e-11 of external
  supply in every run, far inside the preregistered 1e-6 limit, so no arm is invalid on
  accounting grounds. Ice-age cycle off, diagnostic wind, peak RSS 1.0–7.7 GB per run.

## Fixed climate: anchor volume against slice length

Window D, one step from the checkpoint (opening anchor body 20,113 km³ over 0.51M km²;
the tracked nearest body after the step is what is tabulated).

| Total | Slice (substeps) | Anchor volume (km³) | Change from previous row |
|---|---|---:|---:|
| 5 kyr | 5,000 yr (1) | 348,607 | — |
| 5 kyr | 1,000 yr (5) | 366,435 | +5.1% |
| 5 kyr | 200 yr (25) | 389,495 | +6.3% |
| 5 kyr | 50 yr (100) | 393,311 | +1.0% |
| 5 kyr | 10 yr (500) | 394,293 | +0.25% |
| 25 kyr | 25,000 yr (1) | 392,527 | — |
| 25 kyr | 5,000 yr (5) | 408,890 | +4.2% |
| 25 kyr | 1,000 yr (25) | 410,447 | +0.38% |
| 25 kyr | 200 yr (125) | 410,649 | +0.05% |
| 100 kyr | 100,000 yr (1) | 410,814 | — |
| 100 kyr | 10,000 yr (10) | 415,348 | +1.1% |
| 100 kyr | 1,000 yr (100) | 415,348 | 0.00% |
| 100 kyr | 200 yr (500) | 415,348 | 0.00% |

Window E, one 25 kyr step: 999,145 km³ in one slice, 996,545 km³ at 1,000 year slices
(−0.26%), 996,540 km³ at 200 year slices. Window E's anchor body is already converged in a
single slice at this total; its world lake volume still moves 0.8% (below).

Reading: a converged slice exists on current code. It is about 50 years for a 5 kyr total,
about 1,000 years for 25 kyr and 100 kyr totals, under the preregistered 1% tolerance. The
single-slice error that production runs carry today is therefore about 13% on a 5 kyr step's
anchor lake, about 4% at 25 kyr, and about 1% at the production 100 kyr step. These are far
below the factor-of-two-to-seven slice effects #811 measured on the pinned source, where the
same 5 kyr step went from 0.62M km³ in one slice to 0.29M km³ at 200 year slices — and in the
opposite direction (longer slices filled further there; here they store less). The solver has
changed since (sublimation moved into the solve in #757, enclosed basins reclassified in
#758, river substeps in #867), so the pinned numbers describe that source, not this one.

World lake volume moves less than the anchor body, as expected for one marginal basin among
many: at 5 kyr, 3,609,905 km³ in one slice against 3,769,677 km³ converged (+4.4%); at
100 kyr, 3,455,711 against 3,471,588 km³ (+0.46%).

## Mechanism: end-of-slice losses

The solver's per-body update (`_integrate_body_volume` in
[`backend/app/hydroclimate/lakes.py`](../../backend/app/hydroclimate/lakes.py)) solves
implicitly: storage at the end of the slice, plus evaporation and sublimation evaluated at
that end-of-slice water level multiplied by the whole slice length, equals the water
available. Evaporation and sublimation both scale with the wet area at that level, so a
filling lake is charged for the entire slice as if it had been at its final size throughout.
That is a backward-Euler loss term: stable, and first-order accurate in the slice length.
Two related chargings compound it: when bodies merge or receive a transfer mid-slice, the
flush passes re-solve the combined body over the same slice with its losses restored and
recomputed, so the merged body pays losses on its combined area for the full slice however
late in the slice the merger happened (this behaviour is deliberate and tested, #768).

The ledger decomposition at window D, 5 kyr total, one slice against 10 year slices, shows
exactly this signature (km³; external supply identical at 103,061,858):

| Term | One slice | 10 yr slices | Difference |
|---|---:|---:|---:|
| Evaporation | 33,783,859 | 33,483,016 | −300,843 |
| Sublimation | 11,366,329 | 11,138,884 | −227,445 |
| Overflow out | 55,825,108 | 56,193,301 | +368,193 |
| Storage change | 2,086,563 | 2,246,657 | +160,094 |

Shorter slices reduce the two loss terms by 528,288 km³ combined; 368,193 km³ of that saving
leaves as extra overflow from bodies at their caps, and 160,094 km³ stays as storage. Every
term moves monotonically as the slices shorten, and the evaporation error shrinks roughly
fivefold between 1,000 and 200 year slices, consistent with an integration error that
collapses once slice boundaries resolve the filling. Merge events move only from 88 to 94
and splits from 42 to 47 across the same ladder, and the anchor body's catchment never
changes, so discrete merge/spill timing is not the carrier. The preregistered prediction
(overflow and merge counts would carry the effect, losses would move little) is therefore
falsified in its mechanism and confirmed in its accounting: the effect is real, ledger-clean,
and carried by the loss terms.

A synthetic check supports the attribution. In the whole-cell fixtures of
`backend/unit_tests/test_lakes.py` (single basin, and the two-basin merge fixture), stored
volume is identical at 1 to 256 substeps because those bodies fill to their caps, where
volume is pinned; their evaporation still moves by about 3.5% between 1 and 256 substeps,
compensated by overflow. Slice dependence needs a body that is filling below its cap under
level-dependent losses — the situation of the real anchor basins — not merely a saddle or a
spill.

## Climate refreshed: slice length separated from refresh interval

Window D, 25 kyr total, `step_world` repeated at the stated step length (the lake solver
runs one slice per step, as in production; climate, routing, ice and plates are all
recomputed every step).

| Step length | Steps | Anchor volume at 25 kyr (km³) | World lake volume (km³) | Merged steps |
|---|---:|---:|---:|---:|
| 5,000 yr | 5 | 606,839 | 5,035,899 | 0 of 5 |
| 1,000 yr | 25 | 605,159 | 5,615,763 | 0 of 25 |
| 200 yr | 125 | 604,656 | 5,751,589 | 0 of 125 |

Two comparisons separate the effects the issue asked to separate:

- **Climate refresh dominates the one-step slice effect.** At the same 25 kyr total and the
  same 5,000 year lake slice, refreshing the climate every step gives an anchor volume of
  606,839 km³ against 408,890 km³ with the climate held fixed: a 48% difference from the
  refresh (and from every other process being re-stepped), against the 4% single-slice error
  in the fixed arm at that total. The one-step number a gate reads is therefore mostly a
  statement about the step schedule as a whole, and only secondarily about the lake slice.
- **Within the refreshed arm, the anchor volume is already converged at 5,000 year steps**
  (606,839 km³ at 5,000 year steps against 605,159 at 1,000 and 604,656 at 200, a spread of
  0.4%). The fixed arm shows why this does not contradict the one-step result: the lake-slice
  error at a 5 kyr step is about 13% of the anchor body in isolation, but in a refreshed
  trajectory the body is repeatedly re-supplied and sits near the same state each step, so
  the per-step biases do not accumulate into the endpoint the way the one-step measurement
  suggests. World lake volume does move with step length in this arm (5.04M km³ at 5 kyr
  steps against 5.62M at 1 kyr and 5.75M at 200 yr steps, +14%), but that is the whole
  model's step-length sensitivity — plates, ice and routing included — not the lake slice,
  which the fixed arm isolates at +0.46% for a 100 kyr step.

## Does the merged state survive at converged slices?

No. Refreshed-climate trajectories at the converged 1,000 year slice, 100 kyr from each
window's checkpoint, recording every second step (50 records per window):

| Window | Merged records | Largest catchment seen | Largest anchor volume | Anchor volume at 100 kyr |
|---|---:|---:|---:|---:|
| D | 0 of 50 | 1.27M km² | 606,796 km³ | 604,335 km³ |
| E | 0 of 50 | 1.46M km² | 1,003,906 km³ | 614,453 km³ |

The merged state needs a catchment above 15M km²; neither window's tracked body ever
exceeds 1.5M km², a factor of ten below the threshold, at the converged slice. Nor is this
an artefact of choosing the converged slice: in the fixed-climate arm no slice length at
any total (5, 25 or 100 kyr, slices from the whole step down to 10 years) puts either
window's body in the merged state, and in the refreshed 25 kyr ladder no step length from
200 years to 5 kyr does either (0 of 155 steps). On the pinned source, #811 found the
merged state entered only at 2.5 and 5 kyr steps, in the first 20–30 kyr, and already
absent at 1 kyr steps over 0.4 Myr (0 of 400 steps in both windows) — the horizon here,
100 kyr, covers that entry window three times over. The registered D/E merged readings
therefore do not survive at converged slices on either source; on current code the basin
that produced them does not form the merged state at any slice length tested, which is
consistent with #811's finding that the registered events no longer occur on current code
at all (they were traced to the lake-ice sublimation charging fixed by #757 and the basin
reclassification in #758, a separate effect from this one).

## Recommendation

**File one standalone implementation issue; no RFC, and not "no change".**

The mechanism is confirmed, the remedy's shape is already established elsewhere in the
model, and no scientific choice is open: the ice sheet marches in 250 year substeps
(RFC 0015 D2) and the river reservoirs integrate on their own substeps (#867) precisely so
that a long tectonic step does not change their answers. The lake solve is the remaining
solver that takes the tectonic step in one slice. Substepping it changes no equation, no
state and no parameter; it only evaluates the existing loss terms on shorter slices, which
is a verification fix under principle 3, not a model change needing an RFC. It is also too
large to fold into this exploration, which does not authorize a model change.

Proposed issue title: **"Substep the lake solve so one-step lake volumes do not depend on
slice length"**.

Proposed scope and acceptance for that issue:

- In `compute_hydrology`, pass a `substeps` count to `solve_lakes_prognostic` so the lake
  slice is at most 1,000 years (converged to 0.4% or better at 25–100 kyr totals in windows
  D and E here; at most 50 years would be needed to converge a 5 kyr step to 1%, so the
  acceptance run should re-measure the cap against the same tolerance rather than assume
  1,000 years covers every step length).
- Acceptance: from the window D and E checkpoints, one production 100 kyr step's anchor
  volume matches the converged-slice value within 1%; a solver partition test in the style
  of the existing lake tests shows stored volume converging as `substeps` grows on a
  filling, evaporating basin below its cap; ledger residuals stay at current levels.
- Report the step-time cost at density 4 (forcing 100 substeps at a 100 kyr step added
  about 23% to that step's wall time in this probe's fixed arm: 255 s against 313 s, with
  the lake solve only part of that step; the implementation issue should measure the solve
  alone and consider charging losses on a time-weighted level instead if the cost is
  material).
- RFC 0016's gate reads (windows D and E) should be re-read after the change, since #854
  requires gates to state their step size; on current code the expected shift is about 1%
  on the anchor volume at 100 kyr steps, not the factor-of-several that the pinned source
  suggested, and the D/E merged-state conclusion does not change (previous section).

Until that issue lands, any one-step lake number quoted from a single slice should state
its slice length, per #854; at the production 100 kyr step the single-slice bias measured
here is about −1% on the anchor body and −0.5% on world lake volume.

## Reproduction

```bash
# fixed climate, one step, lake slice forced (window D)
backend/.venv/bin/python bin/debug/measure_lake_slice_convergence.py \
  ~/mantle-bloom-captures/rfc0016/window-rescan/checkpoints/seed910211954-225.1myr.bloom \
  --anchors analysis/rfc0016/target-anchor-window-d.json \
  --mode fixed --total-years 5000 --slice-years 50 --myr 0.005 \
  --out ~/mantle-bloom-captures/issue853/fixed-5000-slice50.json

# climate refreshed, same total in 200 yr steps
backend/.venv/bin/python bin/debug/measure_lake_slice_convergence.py \
  ~/mantle-bloom-captures/rfc0016/window-rescan/checkpoints/seed910211954-225.1myr.bloom \
  --anchors analysis/rfc0016/target-anchor-window-d.json \
  --mode refreshed --total-years 25000 --slice-years 200 --myr 0.025 \
  --out ~/mantle-bloom-captures/issue853/refreshed-25000-slice200.json
```
