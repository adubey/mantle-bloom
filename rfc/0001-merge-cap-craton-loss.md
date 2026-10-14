# RFC 0001: Plate merges must not destroy cratons

- **Status:** Draft
- **Date:** 2026-10-13 (re-baselined 2026-10-14 on the year-based split pace, #403)
- **Related:** #274 (cratons), #340 (craton partition rules), #276 (suture accretion),
  #290 (staged suture placement), #358 (split cooldown in years), #403 (re-baseline); RFCs [0002](0002-craton-trench-attrition.md),
  [0004](0004-proto-cratons.md), [0008](0008-craton-collision-jamming.md) (fusion at jam
  onset)

## Problem

When two continental plates merge, their overlapping columns stack. Crust stacked past the
suture cap (`SUTURE_ACCRETION_MAX_HC_M`, 84 km) is overflow, and the merge tries to spread
it over the merged plate. Cratonic crust in that overflow which finds no room is clipped,
and the craton audit around the merge books the clip as `delaminated_m3`.

In the run that prompted this RFC, plate merges destroyed cratons in a few large events
rather than through steady wear. The run is seed 456862899, about 19.2 Gyr, saved locally
as `mantle-bloom-seed456862899-19191800000y.mbworld`; that file is on the author's machine
and **not available for download**. The numbers below come from the save's
`World.stats_history` via `bin/debug/ledger_history.py`, so anyone with a save can make the
same table.

The save predates #358. It ran post-#340 code (exact commit unknown) and counted the split
cooldown in steps: 2 Myr during its 0.1 Myr steps, then **200 Myr** at 10 Myr steps for the
rest of the run. It split and merged far less often than the model does now. Its figures are
kept as the record of how the problem was found; the current ones follow the table.

| When | Craton booked to `delaminated_m3` | Share of live craton |
|---|---|---|
| 602 My | 42 M km³ | 2% |
| **3,052 My** | **364 M km³** | **28%** |
| 6,022 My | 26 M km³ | 5% |
| 6,922 My | 34 M km³ | 15% |
| whole run | 518 M km³ | 19% of all craton destruction |

The single 10 Myr step ending at 3,051.8 My took more than a quarter of the world's
cratons. The same step booked 339 M km³ of continental material to
`collision_subducted_m3`, the account the merge uses for continental material left above
its column (`merge_split._fuse_plates`). The suture path also books its unplaceable
remainder there, so the save can't separate the two; the gate will. Genuine continental
area also fell from 61 to 45 M km².

**On the current model** (origin/main 2562f63, 5 Myr steps; `analysis/issue403/`), merges are
frequent and the merge cap is the largest craton sink while continents are large:

| Run | Merges per Gyr | `delaminated_m3`, % of start craton | Steps booking > 1% of live craton | Largest step, % of live craton |
|---|---|---|---|---|
| Seed 7, 1 Gyr | 55 | 54% | 20 | 17% (365 My) |
| Seed 11, 1 Gyr | 27 | 18% | 7 | 13% (565 My) |
| Seed 2024, 1 Gyr | 36 | 31% | 5 | 14% (715 My) |
| Seed 456862899, 3 Gyr | 12, 1, 2 | 28% (all before 420 My) | 6 | 14% (125 My) |

A step can hold more than one merge, so the step figures bound a single merge's loss from
above. RFC 0008's gate (`analysis/rfc0008/`) found that craton lost in fusions equals each
tuning seed's `delaminated_m3`, and that the cratons the fused pairs stack sum to 53–90% of
the starting craton. The "deep craton-on-craton overlap" row of the gate table below already
applies on every tuning seed.

Two things are wrong:

1. **The booking misleads.** `delaminated_m3` reads as root foundering, which the model
   otherwise forbids for cratons (#340). Here it is the merge cap.
2. **The physics is wrong.** 364 M km³ is about 10 M km² of 35 km craton. Merging welded
   two cratons by *stacking one on top of the other* over that area, then removed what the
   column couldn't hold. Real cratons don't overlap. They jam: the craton brake
   (`torque.CRATON_COLLISION_FRICTION_GAIN`) exists for exactly that reason. A merge that
   finds craton on craton means the overlap built up before the merge.

**What the current placement already does.** `quad_merge._place_stacking_overflow` passes
the overflow to `quad_tectonics._place_suture_crust`. That function's far-field stage walks
the whole adjacency graph and then calls `fill(eligible)`, which takes every remaining
eligible cell on the plate, connected or not. Placement is therefore already global. Volume
goes unplaced only when *every eligible cell* is at the cap.

A merged plate with ~50 M km² of continental cells at ~35 km has ~2.4 billion km³ of room,
seven times the overflow. If the 339 M km³ was mostly the merge's clip, continental
material went unplaced along with the craton, which would mean the Hc overflow itself found
no room, not just the craton riding on it. The candidate explanations:

- **H1: eligible capacity really was exhausted.** The numbers above make this unlikely
  unless the eligible set was small, which is H2.
- **H2: `eligible` is narrower than intended.** It is the merged plate's continental cells,
  resolved from the post-merge `crust_type_code`, or all cells when none are
  continental. A vote or an `INHERIT` resolution could shrink it.
- **H3: placement succeeded, and a later writer in the same audit window clipped the
  craton.**
- **H4: the craton isn't carried with the Hc it rides on.** `craton_moving` and the
  `_carry_material` shares leave craton above some column.

The run is a stress test of the model's long-term behaviour, not a calibration interval.
Its numbers show *that* a merge can remove a quarter of the cratons in one step. They
don't say how much craton a merge should cost; geologically, the answer is close to none.

## Proposed solutions

**What happens in nature.**
- Collisions between continents that both carry cratons (India–Asia is the modern case)
  thicken crust to about 70 km. The thickened crust then escapes laterally along
  strike-slip faults [[Tapponnier et al., 1982]](#ref-tapponnier1982).
- The dense non-cratonic lower crust and mantle root delaminate and sink
  [[Kay & Kay, 1993]](#ref-kay1993).
- Cratonic lithosphere is buoyant and stiff [[Jordan, 1978]](#ref-jordan1978)
  [[Sleep, 2005]](#ref-sleep2005). Thickening by lateral compression is one proposed way
  cratonic roots are *built*, not destroyed [[Lee et al., 2011]](#ref-lee2011)
  [[Wang et al., 2018]](#ref-wang2018).

Nothing in that picture removes tens of percent of a craton in one event.

1. **Give the merge's overflow the suture's outlets.** `_place_stacking_overflow` calls
   `_place_suture_crust` with `strike=None` and `root_capacity=None`, which skips two of
   its five stages: tectonic escape and root delamination.
   - Passing the merged plate's root capacity lets non-cratonic roots shed to make room,
     which matches [[Kay & Kay, 1993]](#ref-kay1993).
   - Passing a strike (the overlap's long axis) enables lateral escape, which matches
     [[Tapponnier et al., 1982]](#ref-tapponnier1982).

   Both add room near the front, so the far-field stage has less to place.
2. **Bound the placement itself.** Today's global `fill(eligible)` can drop crust anywhere
   on the plate: it conserves volume, but it isn't physically plausible. Restrict the
   *primary* call's `eligible` to continental cells within `MERGE_PLACEMENT_MAX_KM` of the
   stacked area and connected to it through continental cells. A second pass can't help:
   the first call already searched everything eligible.

   This makes placement more honest, not more generous, and on its own it may *increase*
   what is clipped. It belongs with solution 3.
3. **Stop craton-on-craton overlap from building.** Measure the overlapping craton volume
   per plate pair. When the overlap is craton on craton, the brake has failed upstream.
   Either the forced merge (`FORCED_MERGE_OVERLAP_FRACTION`, 30% for 30 Myr) should fire
   *before* the overlap is that deep for craton pairs, or retreat of the non-cratonic crust
   between them should take it.

   **Moved to [RFC 0008](0008-craton-collision-jamming.md).** It fuses a craton-jammed pair
   at a hazard rate from jam onset, and stops a stalled jam from resetting its collision
   time, so craton pairs fuse before deep overlap builds. This RFC keeps the measurement
   (the gate records the craton-on-craton overlap before each merge). 0008's acceptance
   criteria include a 75% drop in that overlap at fusion.
4. **Book honestly.** Merge clips go to a new craton account `merge_capped_m3`, kept apart
   from `delaminated_m3`, so the ledger says what happened.

**Recommendation:** a diagnosis gate first. No fix is chosen until the 3,052 My event, or
another merge with the same signature, has been reproduced and its cause established (see
Design). Solution 4 is safe whatever the cause and can land with the diagnosis. The rest
follow from the gate:

| Gate outcome | Fix |
|---|---|
| H1: eligible capacity exhausted | 3, with 1 to add room near the front |
| H2: `eligible` narrower than intended | Fix the mask, then re-measure |
| H3: a later writer clipped placed craton | Fix that writer |
| H4: craton not carried with its Hc | Fix the carry |
| Deep craton-on-craton overlap existed before the merge, whatever else the gate finds | 3 as the primary fix, since the stacking itself is unphysical |

Solution 2 is a validity fix that goes with whichever row applies. It's measured with
solution 3 in place, since alone it would raise the booked loss. Where the table selects 3,
it means RFC 0008's fusion.

**Ordering with RFC 0008.** 0008 makes craton-bearing pairs fuse more often and sooner.
Every extra fusion goes through this RFC's merge placement, so solutions 1, 2 and 4 must land
before 0008's fusion, or with it. With 0008 in place, merges carry less stacked craton, but
there are more of them.

## Design

**Units.** Craton and Hc are thicknesses in metres per cell; volumes are m³ (reported in
M km³). The cap is `SUTURE_ACCRETION_MAX_HC_M` = 84,000 m, and room per cell is
`(cap − Hc) × area`. No per-step or per-node factors are added.

**Diagnosis gate (before any fix but 4).**

*Capturing the merge.* `bin/debug/probe_craton_clips.py` replays a world between two times
and captures the first fusion that leaves craton outside its column. Today it pickles
`(keep, absorb, other_points)`. But `_fuse_plates` calls
`keep.merge_with(absorb, other_points, others)`, and `quad_merge._on_other_plates` uses
the third plates' actual cell footprints to reject cells. The capture must also hold deep
copies of `others`; without them, replaying the tuple can admit cells the real merge
rejected.

*Reproducing the schedule.* Seed plus generation settings is not enough to reproduce a
particular merge: topology and the per-step random streams depend on the number and
length of steps, and on the run's controls.
- **Step schedule.** The save's `stats_history` timestamps are spaced 0.1 Myr apart 218
  times (to 21.8 My), then 10 Myr apart. That suggests one sample per step, which the
  replay must confirm. The event is in the step ending at 3,051.8 My.
- **Controls.** The replay uses the tuning multipliers and simulation flags stored in the
  save, its node density, and its climate setting, since erosion moves crust.
- **Code version.** The run used post-#340 code (its ledger has `audit_clipped_m3` from
  t = 0), but the exact commit is unknown.

If the replay doesn't reproduce a merge at 3,051.8 My, the gate uses the first merge in
the replay with the same signature: craton booked by a merge above 1% of live craton,
together with continental material clipped above its column. The run had three more
(602, 6,022 and 6,922 My).

Since #358 the signature is common on runs anyone can regenerate (table in Problem), so the
gate can start there instead of the save: seed 7 at 5 Myr steps books 17% of its live
craton in the step ending at 365 My, after 365 Myr of stepping rather than 3 Gyr.

*Measuring.* Replay the captured merge alone, printing:
- the overflow volume, the eligible cell count and area, and the room under the cap
  among eligible cells;
- per-stage placement volumes, from `_place_suture_crust`'s `stages`;
- where craton ends up above a column after placement and after the audit;
- for the two plates in the 30 Myr before the merge, their overlap area, the craton
  volume in the overlap, and whether the craton brake was active.

The first leg is about 3 Gyr of stepping. Its wall time must be measured; replays are
RAM-heavy, so run them one at a time. The gate passes when the findings, and the fix the
table selects, are written up in `analysis/` and this RFC is updated with them.

**Merge placement (solutions 1 and 2):**
- Root capacity comes from `orogeny.plate_delamination_capacity_m3(keep, continental,
  years)`, the call the suture path makes. It sheds only non-cratonic roots, restite first
  (#340). The suture placement caps the delaminated share at
  `SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION` (10%) of the volume. `years` is the step
  length, since a merge's roots have had one step to founder like any other belt. The shed
  volume books to `delaminated_lower_crust_m3` with the receiving cells' own provenance,
  as in the suture path.
- Strike is the principal axis of the stacked cells' positions: the eigenvector with the
  largest eigenvalue of their 3×3 covariance, projected onto the tangent plane at their
  centroid.
- The bounded `eligible` (solution 2) is computed once, before the single call: continental
  cells within `MERGE_PLACEMENT_MAX_KM` (great-circle distance) of any stacked cell and in
  the same continental connected component.
- Physical checks on every placement, besides ledger closure:
  - no cell above the 84 km cap;
  - no crust or craton placed on a cell outside the bounded set;
  - placed thickness non-increasing with hop distance from the stacked area (beyond the
    belt's water-fill stage);
  - the merged plate's craton forms no more connected components than the two plates'
    cratons had before the merge.

**Runtime.** Merges happen tens of times per Gyr, not every step.
- The bounded mask is one great-circle distance test and one connected-components pass
  over the merged plate: O(cells), about 10⁵ at density 4.
- Strike is a 3×3 eigenproblem.
- Root capacity is the vector computation the suture path already does.

Together, well under a second per merge, and nothing per ordinary step.

**Save compatibility.** No new per-node field. `merge_capped_m3` is a new key in
`World.craton_ledger`, and `cratons.ensure_ledger` already backfills missing keys at 0.0 on
load. Old saves keep their historical merge losses in `delaminated_m3`; they can't be split
after the fact, and the docs will say so.

**Knobs.** One new knob:
- `MERGE_PLACEMENT_MAX_KM`, starting at 1,000 km. That is the scale over which the Tibetan
  plateau and the blocks extruded from it spread the India–Asia collision's crust
  [[Tapponnier et al., 1982]](#ref-tapponnier1982). It is a calibration assumption, with
  sensitivity runs at 500 and 2,000 km.

Two existing constants gain a second use in the merge, and neither is retuned here:
- `SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION`: how much of a placement may come from
  root shedding.
- `SUTURE_ESCAPE_MAX_ANGLE_DEG`: the escape cone.

If solution 3 is needed, it is built and tuned in RFC 0008 (`CRATON_LOCK_MYR`).

## Detailed design (provisional touchpoints)

- `backend/app/quad_merge.py`
  - `_place_stacking_overflow`: pass `root_capacity` and `strike` to
    `quad_tectonics._place_suture_crust`, with the bounded `eligible` mask; return the shed
    volume per cell for booking.
  - The physical checks in Design, as a debug-mode assertion after each placement.
  - `_transfer` (which calls `_place_stacking_overflow`) and `merge_with` in
    `sparse_quad_patch.py`: thread the step length through, so root capacity can be
    computed. Book the shed roots through `quad_tectonics._book_shed_roots` (the helper
    the suture path uses: `delaminated_lower_crust_m3` and orogeny's
    `restite_delaminated_m3`).
- `backend/app/merge_split.py`: `_merge_plates_with_craton_audit` settles to
  `merge_capped_m3` instead of `delaminated_m3`. Record the craton-on-craton overlap volume
  in `update_overlap_tracking`, as the measurement RFC 0008 builds on.
- `backend/app/cratons.py`: add `merge_capped_m3` to `CRATON_SINK_ACCOUNTS`, with a
  comment citing this RFC.
- `backend/app/stats.py` and the frontend ledger display, if it lists craton accounts by
  name: show the new account.
- `bin/debug/probe_craton_clips.py`
  - Capture `others` (deep copies) with the merge.
  - Accept a step schedule (the 0.1 Myr then 10 Myr pattern, or one read from a save's
    `stats_history`) and the save's controls.
  - Print the per-stage placement breakdown for a captured merge.
- Tests:
  - `backend/unit_tests/test_quad_merge.py`:
    - Two plates whose overlap stacks 2 × 40 km of craton over a patch, with room within
      `MERGE_PLACEMENT_MAX_KM`: after the merge, live craton plus `collision_reworked_m3`
      equals the craton before, and `merge_capped_m3` is 0.
    - The same with room only beyond `MERGE_PLACEMENT_MAX_KM`, or on a disconnected
      continental fragment: nothing lands there, and exactly the excess is booked to
      `merge_capped_m3`.
    - A third plate whose footprint covers candidate cells: those cells are never
      activated, in a direct merge and in a replay of a captured tuple that includes
      `others`.
  - `backend/unit_tests/test_cratons.py`: the ledger balances with the new account, and a
    save without the key loads.
- Docs: `docs/simulation-model.md`, the Cratons section and the merge section of "Merge and
  split".

## Assumptions and uncertainty

- **The cause is unknown.** H1–H4 are hypotheses, and the gate may find another. The fix
  table is written so that each outcome leads somewhere.
- **The event may not replay exactly.** The code version and the per-step random streams
  are only partly recoverable. The gate allows a substitute merge with the same signature.
- **One event, one seed, at the old split pace.** The save's events were measured with a
  200 Myr split cooldown. On the current model the mechanism recurs on every seed measured
  (Problem), so it is no longer a single-seed observation.
- **Escape and root shedding are borrowed.** Both come from the suture path, calibrated for
  ordinary collision belts. Whether they suit a craton-on-craton merge is untested.
  `MERGE_PLACEMENT_MAX_KM` comes from one modern collision.
- **The model stacks; the Earth doesn't.** Real cratons don't overlap, so no published
  figure says where a stacked craton's overflow "should" go. The physical checks bound the
  answer; they don't validate it.

## Acceptance criteria

- [ ] **Gate:** the 3,052 My event (or a substitute with the same signature) is reproduced
      with the full merge inputs, its cause established, and the fix the decision table
      selects is written up in `analysis/` with the per-stage numbers, before any fix
      other than solution 4 is implemented.
- [ ] On a regenerated seed 456862899 to 3 Gyr and 500 Myr runs of seeds 7, 11 and 2024
      (`bin/debug/compare_cratons.py`, 5 Myr steps), merges book less than 1% of the
      initial craton volume in total, and no single merge books more than 0.5% of the live
      craton. On origin/main 2562f63, `delaminated_m3` books 43%, 7.3% and 0.9% of the
      start craton by 500 Myr on the tuning seeds and 28% on seed 456862899. Over 1 Gyr
      (3 Gyr for 456862899), 6–22 steps per run book more than 0.5% of the live craton.
- [ ] Every merge placement passes the physical checks in Design (cap, bounded set,
      monotone thickness, craton component count), checked by an assertion in debug runs.
- [ ] The craton ledger balance error stays within 1e-9 of the live volume
      (`cratons.balance_error_m3`), and the continental-material ledger still closes
      (`continental_ledger.assert_closed`).
- [ ] Merge counts and plate counts on those seeds change by less than 10%, measured
      before RFC 0008's fusion lands. 0008 raises merge counts by design and has its own
      bounds. The baseline (origin/main 2562f63): 37, 16 and 25 merges and 55, 53 and 50
      plates at 500 Myr on seeds 7, 11 and 2024; 15 merges and 58 plates at 3 Gyr on seed
      456862899.
- [ ] The unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.
- [ ] Saves from before the change load, with `merge_capped_m3` = 0.

## Falsification criteria

- **The merge isn't the cause.** If the reproduced event shows the craton was removed by a
  writer other than the merge's placement (it was outside its column before the merge, or
  was clipped on cells the merge never touched), the RFC's premise is wrong. Re-scope it
  to that writer.
- **The loss moves instead of stopping.** If, with the fix, craton destroyed per Gyr on
  the acceptance runs falls by less than half of what merges booked before, the losses have
  moved to another sink and the merge was a symptom. The likely upstream cause is craton
  overlap (solution 3, now RFC 0008's fusion); prioritise that.
- **No plausible room.** If the bounded placement can't hold the overflow in more than half
  of the craton-carrying merges, stacking is the problem, not placement. Solution 3
  becomes the only fix that matters.

## References

- <a id="ref-jordan1978"></a>Jordan, T. H. (1978). Composition and development of the
  continental tectosphere. *Nature*, 274, 544–548.
  [doi:10.1038/274544a0](https://doi.org/10.1038/274544a0)
- <a id="ref-kay1993"></a>Kay, R. W., & Kay, S. M. (1993). Delamination and delamination
  magmatism. *Tectonophysics*, 219, 177–189.
  [doi:10.1016/0040-1951(93)90295-U](https://doi.org/10.1016/0040-1951(93)90295-U)
- <a id="ref-lee2011"></a>Lee, C.-T. A., Luffi, P., & Chin, E. J. (2011). Building and
  destroying continental mantle. *Annual Review of Earth and Planetary Sciences*, 39,
  59–90. [doi:10.1146/annurev-earth-040610-133505](https://doi.org/10.1146/annurev-earth-040610-133505)
- <a id="ref-sleep2005"></a>Sleep, N. H. (2005). Evolution of the continental
  lithosphere. *Annual Review of Earth and Planetary Sciences*, 33, 369–393.
  [doi:10.1146/annurev.earth.33.092203.122643](https://doi.org/10.1146/annurev.earth.33.092203.122643)
- <a id="ref-tapponnier1982"></a>Tapponnier, P., Peltzer, G., Le Dain, A. Y., Armijo, R.,
  & Cobbold, P. (1982). Propagating extrusion tectonics in Asia: New insights from simple
  experiments with plasticine. *Geology*, 10, 611–616.
  [doi:10.1130/0091-7613(1982)10<611:PETIAN>2.0.CO;2](https://doi.org/10.1130/0091-7613(1982)10%3C611:PETIAN%3E2.0.CO;2)
- <a id="ref-wang2018"></a>Wang, H., van Hunen, J., & Pearson, D. G. (2018). Making
  Archean cratonic roots by lateral compression: A two-stage thickening and stabilization
  model. *Tectonophysics*, 746, 562–571.
  [doi:10.1016/j.tecto.2016.12.001](https://doi.org/10.1016/j.tecto.2016.12.001)
