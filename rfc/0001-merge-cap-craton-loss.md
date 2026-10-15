# RFC 0001: Plate merges must not destroy cratons

- **Status:** Accepted (2026-10-15)
- **Date:** 2026-10-13 (evidence rebased 2026-10-14 for #403)
- **Related:** #274 (cratons), #340 (craton partition rules), #276 (suture accretion),
  #290 (staged suture placement); RFCs [0002](0002-craton-trench-attrition.md),
  [0009](0009-craton-mosaic.md) (cratons; supersedes 0004), [0008](0008-craton-collision-jamming.md) (fusion at jam
  onset); #403 (post-#358 baseline audit); #443 (tracking), #440 (diagnosis gate),
  #441 (solution 4)

## Problem

When two continental plates merge, their overlapping columns stack. Crust stacked past the
suture cap (`SUTURE_ACCRETION_MAX_HC_M`, 84 km) is overflow, and the merge tries to spread
it over the merged plate. Cratonic crust in that overflow which finds no room is clipped,
and the craton audit around the merge books the clip as `delaminated_m3`.

This violates design principles 6 (account for material), 7 (preserve cratonic crust and
name its loss path), and 9 (attribute losses before choosing a fix). The gate must expose
the writer and destination; the acceptance runs must bound both individual merge losses
and the cumulative merge sink.

**Historical save evidence (pre-#358; simulation commit unknown).** In the run that
prompted this RFC, plate merges destroyed cratons in a few large events rather than through
steady wear. The run is seed 456862899, about 19.2 Gyr, saved locally
as `mantle-bloom-seed456862899-19191800000y.mbworld`; that file is on the author's machine
and **not available for download**. The numbers below come from the save's
`World.stats_history` via `bin/debug/ledger_history.py`, so anyone with a save can make the
same table. These observations were first recorded in commit `d795a09`; that is the
documentation provenance, not a recovered simulation baseline. The save is post-#340
(its ledger has `audit_clipped_m3` from t = 0) and pre-#358. Its 20-step split cooldown
was **200 Myr during the mostly 10 Myr schedule**, rather than the post-#358 2 Myr
cooldown. All values in the table and the event/room estimates below describe that
historical run, not current `origin/main`.

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

Two things are wrong:

1. **The booking misleads.** `delaminated_m3` reads as root foundering, which the model
   otherwise forbids for cratons (#340). Here it is the merge cap.
2. **The physics is wrong.** 364 M km³ is about 10 M km² of 35 km craton. Merging welded
   two cratons by *stacking one on top of the other* over that area, then removed what the
   column couldn't hold. Real cratons don't overlap. They jam: the craton brake
   (`torque.CRATON_COLLISION_FRICTION_GAIN`) exists for exactly that reason. A merge that
   finds craton on craton means the overlap built up before the merge.

**Current code check (`2562f63`, post-#358).** `quad_merge._place_stacking_overflow` passes
the overflow to `quad_tectonics._place_suture_crust`. That function's far-field stage walks
the whole adjacency graph and then calls `fill(eligible)`, which takes every remaining
eligible cell on the plate, connected or not. Placement is therefore already global. Volume
goes unplaced only when *every eligible cell* is at the cap. The call still passes
`strike=None` and `root_capacity=None`, and returns unplaced material and craton to the
overflowing cells for the merge and phase audit to book. The 2 Myr cooldown changes how
often this path runs; it does not identify H1–H4.

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

**Earlier post-#358 attribution (#395).** The RFC 0008 contact probe on model baseline
`3b9d9c1` plus diagnostic instrumentation, at 5 Myr steps and the 2 Myr cooldown,
attributed fusion losses of 54%, 18% and 31% of starting craton on seeds 7, 11 and 2024
to 1 Gyr. Those losses matched each seed's `delaminated_m3` account; see
[`analysis/rfc0008/README.md`](../analysis/rfc0008/README.md) and its `gate.json`.
This confirms that fusion still carries a large sink after #358. It does not measure
eligible room or per-stage placement, so H1–H4 and the placement gate remain open.

**Fresh post-#358 baseline (#403, `2562f63`).** Generated worlds use the default
controls and node density 1, climate enabled, 5 Myr steps, and the 2 Myr cooldown.
`bin/debug/rebaseline_rfcs.py` wraps `compare_cratons.run_job` with a read-only tap on
each `_merge_plates_with_craton_audit` call. The method and controls are recorded in
[`analysis/issue403/README.md`](../analysis/issue403/README.md); compact measurements
are in [`tuning.json`](../analysis/issue403/tuning.json).

| Seed, interval | Merge loss, M km³ | Share of starting craton | Largest merge, share of pre-merge live craton | Merges / splits |
|---|---|---|---|---|
| 7, 0–1 Gyr | 2,216.52 | 54.09% | 16.66% at 365 Myr | 55 / 115 |
| 11, 0–1 Gyr | 583.77 | 18.42% | 12.54% at 565 Myr | 27 / 109 |
| 2024, 0–1 Gyr | 532.30 | 30.88% | 14.08% at 715 Myr | 36 / 104 |

Merges losing more than 1 km³ number 27, 7 and 9, respectively; their summed losses
equal `delaminated_m3` on every seed. The largest merge on each seed lost 387.96,
280.87 and 157.48 M km³, respectively. Plate counts rose from 20 to 51, 9 to 57 and
11 to 60. These are present substitute events for the historical discovery, but the
tap has not captured their full placement inputs. The recommendation remains
diagnosis first, then the H1–H4 decision table; the per-stage placement gate is not
passed by these measurements. On regenerated seed 456862899, the first Gyr has 12
merges and six losses above 1 km³; the largest removes 181.00 M km³ (13.77% of live
craton) at 125 Myr. No merge-attributed loss occurs in the second or third Gyr. This
adds a substitute capture candidate but still lacks the placement inputs required by
H1–H4; the full window table is in [the #403 evidence](../analysis/issue403/README.md).

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
2. **Bound the placement itself.** The baseline's global `fill(eligible)` (`2562f63`) can drop crust anywhere
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

**Recommendation:** the placement diagnosis gate still comes first. The post-#358 probe
confirms fusion attribution, but no placement fix is chosen until a merge with the
historical event's signature has been captured with full inputs and its cause established
(see Design). Reproducing the historical 3,052 My event is optional because its exact
simulation commit is unknown. Solution 4 is safe whatever the cause and can land with the
diagnosis. The rest
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
- **Code version and cooldown.** The historical run used post-#340, pre-#358 code, but
  its exact simulation commit is unknown (evidence recorded in `d795a09`). Reproducing
  that history requires the old 20-step cooldown: 2 Myr in the initial 0.1 Myr steps,
  then 200 Myr in the 10 Myr steps. A fresh run on the post-#358 2 Myr cooldown is a
  substitute-event search, not a replay of that save's topology history.

If the replay doesn't reproduce a merge at 3,051.8 My, the gate uses the first merge in
the replay with the same signature: craton booked by a merge above 1% of live craton,
together with continental material clipped above its column. The run had three more
(602, 6,022 and 6,922 My).

*Measuring.* Replay the captured merge alone, printing:
- the overflow volume, the eligible cell count and area, and the room under the cap
  among eligible cells;
- per-stage placement volumes, from `_place_suture_crust`'s `stages`;
- where craton ends up above a column after placement and after the audit;
- for the two plates in the 30 Myr before the merge, their overlap area, the craton
  volume in the overlap, and whether the craton brake was active.

The historical event needs about 3 Gyr of stepping; the fresh seed-7 substitute is at
365 Myr. Measure wall time for whichever capture is selected. Replays are RAM-heavy,
so run them one at a time. The gate passes when the findings, and the fix the
table selects, are written up in `analysis/` and this RFC is updated with them.

**Gate outcome (#440, 2026-10-15; [`analysis/rfc0001/`](../analysis/rfc0001/README.md)).**
Seeds 7, 11, 2024 and 456862899 ran to 1 Gyr on `788e654`, which is after RFC 0008's #399–#401,
with the #403 controls. Merges cost 21.0%, 6.0%, 2.7% and 12.2% of starting craton, against
54%, 18% and 31% on `2562f63`. The run captured the 15 merges with the signature (1.0–12.1%
of live craton each, 1,155 of 1,275 M km³ of merge loss), and every replay reproduced its
booked loss exactly.
- **H1 holds** for 97.5% of the loss (14 events). The overflow (207–798 M km³) exceeded all
  the room under the cap on the eligible cells (69–286 M km³). The colliding continents
  already averaged 41–84 km of Hc, not ~35 km, so the room estimate in Problem doesn't hold.
- **H4 holds** for the other 2.5% (seed 11, 195 Myr). It is also present in 10 of 15 events.
  `quad_merge._transfer`'s `extensive` restores each field by its own ratio, which leaves
  new cells with up to 63 M km³ more craton than crust. `_carry_material` then places craton
  in proportion to Hc gain, whatever that gain is.
- **H2 is minor:** thick crust outside the mask has 4–129 M km³ of room per event, never
  enough for the unplaced volume. **H3 does not hold:** no excess is found away from the front
  or the receivers.
- **Deep craton-on-craton overlap** preceded 13 of 15 merges (up to 290 M km³). Every pair
  was jammed throughout the previous 30 Myr.
- **A source outside H1–H4:** 70% of the overflow arises on the merge's *new* cells. The
  absorbed plate's crust on third plates (7–44% of its area, 4,153 M km³ of Hc across the
  events) is never activated (#344), so `extensive` packs it into the new territory past the
  cap.

Selected by the table: **3 (RFC 0008's fusion at jam onset, #402) as the primary fix, and
fixing the carry (H4).** The "no plausible room" falsification criterion is met for solution
2: even unbounded placement couldn't hold the overflow in 14 of 15 merges. Solution 1 adds
at most 10% of the overflow through root shedding. Both are deferred until fusion and the
carry fix are measured. The third-plate compression needs a new solution, which is an
amendment to this RFC.

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

**Runtime estimate (not measured).** This work runs on merges, not every step; the
post-#358 replay records how often that path executes.
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
- **One event, one seed.** The 602, 6,022 and 6,922 My events suggest the mechanism
  recurs, but how often it fires across seeds is unmeasured; the acceptance runs measure
  it.
- **Escape and root shedding are borrowed.** Both come from the suture path, calibrated for
  ordinary collision belts. Whether they suit a craton-on-craton merge is untested.
  `MERGE_PLACEMENT_MAX_KM` comes from one modern collision.
- **The model stacks; the Earth doesn't.** Real cratons don't overlap, so no published
  figure says where a stacked craton's overflow "should" go. The physical checks bound the
  answer; they don't validate it.

## Acceptance criteria

- [x] **Gate:** a post-#358 merge with the historical signature is reproduced
      with the full merge inputs, its cause established, and the fix the decision table
      selects is written up in `analysis/` with the per-stage numbers, before any fix
      other than solution 4 is implemented. #395's fusion attribution is partial
      evidence; it does not satisfy this placement gate. *Done in #440: 15 merges captured
      and reproduced; H1 and H4, with prior craton-on-craton overlap; see the gate outcome
      in Design.*
- [ ] On a 3.5 Gyr replay of seed 456862899 and 500 Myr runs of seeds 7, 11 and 2024
      (`bin/debug/compare_cratons.py`), merges book less than 1% of the initial craton
      volume in total, and no single merge books more than 0.5% of the live craton.
      Use the #403 controls and 5 Myr steps. Its 3 Gyr stress baseline must be extended
      to 3.5 Gyr on pinned `2562f63` before claiming the longer criterion.
- [ ] Every merge placement passes the physical checks in Design (cap, bounded set,
      monotone thickness, craton component count), checked by an assertion in debug runs.
- [ ] The craton ledger balance error stays within 1e-9 of the live volume
      (`cratons.balance_error_m3`), and the continental-material ledger still closes
      (`continental_ledger.assert_closed`).
- [ ] Merge counts and plate counts on those seeds change by less than 10%, measured
      against pinned `2562f63` at matched elapsed times, before RFC 0008's fusion lands.
      0008 raises merge counts by design and has its own bounds.
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
