# RFC 0014 follow-up: large-lake flicker mechanisms

## Scope and provenance

This investigation extends the final RFC 0014 acceptance analysis in [#673](https://github.com/adubey/mantle-bloom/issues/673) and tracks the large-lake follow-up in [#675](https://github.com/adubey/mantle-bloom/issues/675). It does not change model behavior or the accepted flicker threshold.

The registered seed-997271774 save (SHA-256 `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7`) was replayed for 0.8 Myr with ice-age control `never`. The source files for lake solving, hydrology, erosion and world stepping match the #673 model-source checksums. The replay worktree was at `f22bda1785a8792453e4513e95d2655d41267eb5`; its measured model base is `6ebf8cd19ee0098822ff9eae6685a22e20f9db32`. The local save is not available for download. The compact results and eight geographic anchor vectors are in this directory.

Three dynamic schedules were compared: 100 kyr (8 steps), 25 kyr (32 steps), and 10 kyr (80 steps). At each step, `trace_lake_anchor.py` finds the nearest grid cell to each anchor and records all connected wet components containing those cells. Anchors are locations, not persistent lake IDs. The maximum nearest-cell distance was 2.9 km in the longest replay.

## Results

### Timestep changes the connected basin at the target location

The eight anchors identify the shallow 59.62 km³ component at 131.2 Myr in the registered 100 kyr replay. The same locations give very different component sizes under the three dynamic step schedules:

| Elapsed age | 100 kyr | 25 kyr | 10 kyr |
|---:|---:|---:|---:|
| 130.55 Myr | next checkpoint at 130.6: 20,090 km³ / 931,713 km² | 19,870 km³ / 931,713 km² | 67,807 km³ / 1,310,953 km² |
| 131.2 Myr | 59.62 km³ / 31,119 km² | 77,620 km³ / 1,446,251 km² | 84,985 km³ / 1,489,253 km² |
| 131.3 Myr | 157.59 km³ / 65,796 km² | 79,663 km³ / 1,462,044 km² | 86,642 km³ / 1,512,189 km² |

All size entries are volume / flooded area. In the 10 kyr run, the anchor component grows from 19,903 km³ / 931,713 km² at 130.54 Myr to 67,807 km³ / 1,310,953 km² at 130.55 Myr. In the 25 kyr run, connectivity oscillates during the transition: 36,262 km³ / 1,122,339 km² at 130.6 Myr, 67,465 km³ / 1,314,945 km² at 130.65, 45,921 km³ / 1,233,657 km² at 130.7, and 67,088 km³ / 1,314,945 km² at 130.725.

The 100 kyr trajectory also alternates wet and dry at the anchor location: it maps to no wet component at 130.7 Myr, to 49.1 km³ / 23,381 km² at 130.8, 39.7 km³ / 23,381 km² at 130.9, no wet component at 131.0, 147.8 km³ / 61,995 km² at 131.1, 59.62 km³ / 31,119 km² at 131.2, and 157.59 km³ / 65,796 km² at 131.3.

The standard #588 large-body reversal detector counted zero qualifying reversals in the 10 kyr run. That detector asks for a sign reversal in consecutive volume changes above a threshold; it does not flag the sustained connectivity change seen by geographic matching. This is not a pass of the RFC 0014 flicker criterion.

At 131.2 and 131.3 Myr, the 10 and 25 kyr schedules agree within about 10% in component volume and 3% in area, while the 100 kyr trajectory remains a much smaller, intermittently connected basin. The transition timing is not monotonic in step size: the 10 kyr schedule merges earlier than 25 kyr, and the 100 kyr trajectory follows another path. Because every world process updates at each chosen step, this is evidence of dynamic trajectory and topology sensitivity, not a lake-solver-only or single-mechanism test.

All global lake ledgers remain closed. Maximum normalized residuals were `1.93e-13` for 10 kyr and `5.58e-13` for 25 kyr, below the 0.1% requirement.

### Evaporation scale

The current temperature-only evaporation proxy averaged 1.445 m/year on the shallow seed-997271774 footprint immediately before the 100 kyr rebound (10th–90th percentile 1.440–1.451 m/year). Its 59.62 km³ storage over 31,119 km² corresponds to about 1.33 years of no-inflow storage at that rate. For the deep seed-373936961 component immediately before its step-5 drop, the proxy averaged 0.817 m/year; its 136,570 km³ over 185,010 km² corresponds to about 904 years of no-inflow storage.

These are diagnostic turnover ratios, not drought forecasts: catchment supply continues in the replay. The saves do not carry wind, humidity, radiation, or deep-water temperature, so this does not test wind-enhanced evaporation or vertical thermal moderation. The rates do show why shallow storage is sensitive to the interaction of supply, evaporation and event timing.

### Sediment scale and placement

At the current modeled deposition rate, the prior step's silt on the shallow lake's eight-cell footprint is about 56 km³, roughly 57% of its 97.97 km³ rebound. Deposition is applied after that step's water solve; the next solve still expands the connected area and raises storage. On the deep event's 44-cell footprint, prior-step deposition is about 416 km³, roughly 0.3% of the 132,450 km³ component drop. These are equal-area footprint-volume screens, not sediment-only counterfactuals.

For the shallow case, the solver's one-hop spill target before expansion received about 0.04 m of silt in the prior step. The reported target after expansion was a different dry cell. The one-hop target is not a surveyed geomorphic sill, so these figures do not establish connector aggradation or sediment causality.

### Subsurface storage screens

Earlier capacity-only estimates remain upper-bound scale checks, not observed or dynamically exchanged water:

- Shallow-basin pore capacity is about 18.6 km³ for a 1 m sediment layer at 0.3 porosity, or 248 km³ for 10 m at 0.4 porosity. The smaller case is about 19% of the rebound; the larger exceeds it but is not all empty, connected or exchangeable.
- A 60 km nearshore band, 3 m head rise and specific yield 0.2 gives about 32 km³ of shallow-basin bank storage, about one-third of the rebound. Actual capacity and return times depend on lakebed and aquifer geology.

No porewater or shoreline-aquifer storage was added to the model in this investigation. The scratch per-body transfer tap was discarded because hierarchy flush passes recount gross transfers; the global ledger cannot identify event-basin ancestry.

## Interpretation and next test

G6 passes on the registered event windows: each basin budget closes and more than 99% of
component storage change is assigned to measured fluxes or explicit connectivity reassignment.
The owner revised G7 to characterize partition sensitivity rather than require invariance; the
#687 fixed-forcing results are therefore a sensitivity measurement, not by themselves a solver
failure. G8 can proceed after its diagnosis windows and parameter ranges are preregistered, with
the intended production schedule held fixed within each paired comparison. The RFC's existing
flicker threshold remains open and unchanged.

## G6 event replay progress

The two registered #538 saves were replayed from their exact checksums on base commit
`2cefff3b3bb073c0d791e8444ec1fac3db66db15`, with `ice_age_period_years = 0` and the registered
100 kyr schedule. The basin-accounting diagnostics and source-file SHA-256s are recorded in
[`g6-basin-budgets.json`](g6-basin-budgets.json). The saved-world replay also records the
whole-world ledger residual at each event step.

| Seed / event | Geographic trace | Whole-world normalized residual |
|---|---|---:|
| 373936961, step 5 (244.9 Myr) | At step 4, the selected component matches the registered 136,570 km³ / 185,010 km² target. At step 5, eight spatially dispersed anchors map to three wet components: 26,974 km³ (4 anchors), 4,120 km³ (2 anchors), and 4.75 km³ (1 anchor); one anchor is dry. | `3.82e-13` |
| 997271774, steps 7–8 (131.2–131.3 Myr) | The registered eight anchors map to 59.62 km³ / 31,119 km² at step 7 and 157.59 km³ / 65,796 km² at step 8. | `6.27e-15` / `7.09e-16` |

This reproduces the event geometries and confirms whole-world closure. The first basin trace
was invalid because it treated each plate's local surface-node key as globally unique. The
resulting ID collisions collapsed cells from different plates and made tracked storage zero.
The corrected trace pairs each local key with its plate ID, checks uniqueness, and evaluates
old and new storage on the same current basin boundary used for source and spill attribution.
Whole-world source/loss totals agree with the ledger to roundoff. The corrected event budgets
are:

| Seed / event | Component storage change | Basin external supply | Basin residual / supply | Component attribution |
|---|---:|---:|---:|---:|
| 373936961, step 5 | `-105,533 km³` | `29.50 million km³` | `8.4e-17` | `100.0%` |
| 997271774, step 7 | `-88.149 km³` | `106.92 million km³` | `8.6e-18` | `100.0%` |
| 997271774, step 8 | `+97.968 km³` | `67.92 million km³` | `3.8e-18` | `100.0%` |

The per-basin residuals are `-2.47 m³`, `-0.92 m³`, and `+0.26 m³`, respectively. Explicit
connectivity reassignment is `-31.439 km³` for the deep drop, `-34.254 km³` for the shallow
drop, and `+8.175 km³` for the rebound. The destination-aware boundary trace records lake
receiver IDs 269, 741 and 801 in the deep event, IDs 15 and 16 in the shallow drop, and IDs 10
and 13 in the rebound. It classifies the deep event's receiver-less export to stable node
`[360287971766698055, 14]` as a non-lake sink; this deferred hierarchy-cycle export is not
assigned to a receiving lake, and the node is not ocean. None of the three event traces records
an ocean-target transfer. The per-edge records, including source/receiver IDs and direction,
are in the JSON; paired in/out entries are boundary ledger records, while
`net_boundary_transfer_m3` is the quantity used in the closed basin balance. Thus the basin
budget closes and every observed boundary has an explicit lake receiver or terminal sink node,
but this trace does not assign that terminal export to an ocean or receiving lake.
All three event steps meet the G6 limits (residual at most 0.1% of supply and at least 99% of
component storage change assigned). This closes the registered G6 accounting cases; it does not
rank physical mechanisms or itself authorize G8.

## G7 fixed-forcing partition control

The exact #538 save checksums were replayed at the registered 100 kyr schedule to capture the
pre-event lake state and the corresponding one-interval forcing: before deep event step 5 and
before shallow rebound step 8. From each checkpoint, the same terrain, ocean mask, hierarchy
and outlets, initial liquid/ice storage, node-wise lake supply and loss rates, and ice climate
were held fixed for a 100 kyr interval. The interval was integrated as 10 × 10 kyr, 4 × 25 kyr,
and 1 × 100 kyr. Each arm therefore received the same integrated prescribed source input.
These are direct lake-solver controls; no world processes advanced between partition steps.

**The former G7 invariance rule would be exceeded.** Global liquid-plus-ice storage differs
from the 100 kyr arm by 4.963% (10 kyr)
and 8.846% (25 kyr) in the deep case, and by 3.586% and 0.439% in the shallow case. In the
deep case, the anchor component storage differs by 0.334% for both finer partitions, but the
anchor connectivity differs: 10/25 kyr leave four anchors wet across three components, while
100 kyr leaves seven anchors wet across three components. In the shallow case, all eight
anchors remain connected and their 157.586 km³ endpoint storage is identical, even while total
global lake storage differs. The global component count is 43 versus 42 in the deep case and
15 in every shallow arm. Merge/split events occur during the first partition interval in each
arm; finer partitions record more within-interval hierarchy transitions. All normalized
whole-solve liquid-plus-ice residuals are below `7.9e-13`, comfortably within the 0.1% limit.

This isolates a solver/timestep/topology sensitivity under identical prescribed forcing. It
does not identify a physical cause. Under the owner's revised G7 gate, report this sensitivity
without treating it alone as solver failure. The focused trace is tracked in
[#687](https://github.com/adubey/mantle-bloom/issues/687); G8 paired comparisons should use the
intended production schedule and identical schedules within each pair.

Reproduce each case from the repository root (replace `SAVE` with the registered local save):

```bash
NUMBA_CACHE_DIR=/tmp/rfc0014-g7-numba-cache backend/.venv/bin/python \
  bin/debug/gate_lake_partition.py deep SAVE --out /tmp/rfc0014-g7-deep.json
NUMBA_CACHE_DIR=/tmp/rfc0014-g7-numba-cache backend/.venv/bin/python \
  bin/debug/gate_lake_partition.py shallow SAVE --out /tmp/rfc0014-g7-shallow.json
```

### Model-assigned drainage catchments

The replay diagnostic `bin/debug/trace_lake_drainage_phases.py` maps the cells in the governing
lake-hierarchy bodies that overlap the fixed-anchor lake components. Orange cells are dry land
assigned to those bodies; blue cells are connected water. This is the model's catchment
assignment, not an independently derived topographic watershed. In the deep event, member area
is 973,081 km² before and 982,039 km² after the drop, while the connected lake drops from
136,570 to 31,099 km³ and splits into three anchor-containing components. The shallow event's
assigned member area is 2,684,704 km² before the drop, 909,978 km² at the low phase, and
1,901,787 km² on rebound. Relative to the preceding phase, 451 member cells (1,774,726 km²)
leave the assignment at the low phase; 249 cells (991,809 km²) rejoin on rebound. The connected
lake volumes are 147.77, 59.62, and 157.59 km³, respectively. These large membership changes
show why a basin budget must follow cell-level source and transfer ancestry across hierarchy
changes; catchment area change by itself is not water flux.

The superseded trace had collapsed 132,377 pre-drop cells into 59,219 distinct IDs. The corrected
trace uses plate-qualified IDs and rejects duplicates. It also keeps hierarchy membership change
separate: old storage is recomputed over the same current-basin cells used to sum current sources
and spill crossings, while storage change under the old and new memberships is reported as a
separate diagnostic. This prevents water on cells removed from the current hierarchy from being
mistaken for physical outflow.

The phase summary, including the seed checksum and membership transitions, is in
`analysis/issue675/shallow-drainage-phases.json`. The corresponding deep and shallow figures
were rendered from the same replay outputs and are available in the task's visualization
artifacts. Reproduce the shallow map with:

```bash
backend/.venv/bin/python bin/debug/trace_lake_drainage_phases.py <registered-seed-997-save> \
  --anchor-step 6 --anchors analysis/issue675/shallow-lake-anchor.json \
  --phase pre-drop=6 --phase low=7 --phase rebound=8 --myr 0.8 \
  --step-years 100000 --checkpoint-myr 0.1 \
  --out <scratch>/shallow-catchments.npz --figure-out <scratch>/shallow-catchments.png
```

The replay trace helper is `bin/debug/trace_lake_event_anchor.py`. It opts into transient
per-node source/loss arrays, follows wet drainage-body memberships and fixed geographic
anchors, uses plate-qualified node identities, separates anchor-component storage from
tracked-basin storage, evaluates storage on the same current-basin cells as sources and spill
crossings, and reports old/new hierarchy membership stock separately. Traces are diagnostic
only; replay JSON and saves remain outside the repository. The deep event trace can be
reproduced with:

```bash
backend/.venv/bin/python bin/debug/trace_lake_event_anchor.py <registered-seed-373-save> \
  --anchor-step 4 --reference-volume-km3 136570 --myr 0.5 \
  --step-years 100000 --checkpoint-myr 0.1 --out <scratch>/seed373-g6.json
```

For the shallow event, run the same helper with the seed-997 save,
`--anchor-step 6 --anchors analysis/issue675/shallow-lake-anchor.json --myr 0.8`, and the same
100 kyr schedule.

Earlier endpoint and body-overlap attribution experiments used the colliding IDs and are
superseded. The accepted diagnostic boundary is based on stable, plate-qualified endpoint node
IDs; body-overlap classification remains unsuitable because nested hierarchy flush passes can
revisit the same water. Detailed corrected step results are in
`analysis/issue675/g6-basin-budgets.json`.

## Reproduction

Use the registered seed-997271774 `.mbworld` save, which is local and not available for download. The geographic trace script accepts the save and the checked-in anchor coordinates:

```bash
backend/.venv/bin/python bin/debug/trace_lake_anchor.py <save> \
  --anchors analysis/issue675/shallow-lake-anchor.json \
  --backend <backend-at-6ebf8cd-or-source-identical> \
  --myr 0.8 --step-years 10000 --checkpoint-myr 0.01 \
  --out /tmp/seed997271774-anchor-10k.json
```

Repeat with `--step-years 25000 --checkpoint-myr 0.025` and `--step-years 100000 --checkpoint-myr 0.1`. The script sets the ice-age control to `never`, writes each step to the output JSON, and reports geographic anchor-to-component matches and normalized ledger residuals. The model-source SHA-256s used for the registered runs are recorded in `timestep-and-mechanism-screen.json`.
