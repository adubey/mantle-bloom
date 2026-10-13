# RFC 0001: Plate merges must not destroy cratons

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #274 (cratons), #340 (craton partition rules), #276 (suture accretion),
  #290 (staged suture placement); RFCs [0002](0002-craton-trench-attrition.md),
  [0004](0004-proto-cratons.md)

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

| When | Craton booked to `delaminated_m3` | Share of live craton |
|---|---|---|
| 602 My | 42 M km³ | 2% |
| **3,052 My** | **364 M km³** | **28%** |
| 6,022 My | 26 M km³ | 5% |
| 6,922 My | 34 M km³ | 15% |
| whole run | 518 M km³ | 19% of all craton destruction |

The single 10 Myr sample at 3,052 My took more than a quarter of the world's cratons. The
same event booked 339 M km³ of continental material to `collision_subducted_m3`, and
genuine continental area fell from 61 to 45 M km².

Two things are wrong:

1. **The booking misleads.** `delaminated_m3` reads as root foundering, which the model
   otherwise forbids for cratons (#340). Here it is the merge cap.
2. **The physics is wrong.** 364 M km³ is about 10 M km² of 35 km craton. Merging welded
   two cratons by *stacking one on top of the other* over that area, then removed what the
   column couldn't hold. Real cratons don't overlap. They jam: the craton brake
   (`torque.CRATON_COLLISION_FRICTION_GAIN`) exists for exactly that reason. A merge that
   finds craton on craton means the overlap built up before the merge.

We don't yet know why the overflow couldn't be placed. A merged plate with ~50 M km² of
continental cells at ~35 km has ~2.4 billion km³ of room under the cap, seven times the
overflow. The candidate explanations:

- **H1:** the overflow lands only within the hop reach of the stacked area, and the cells
  there were already near the cap.
- **H2:** the `eligible` mask passed to `_place_suture_crust` is narrower than intended.
- **H3:** the placement succeeded but left craton on cells the audit clips for another
  reason (a crust-type vote, or a new cell held to the cap).

## Proposed solutions

**What happens in nature.** Collisions between continents that both carry cratons
(India–Asia is the modern case) thicken crust to about 70 km. The thickened crust then
escapes laterally along strike-slip faults [[Tapponnier et al., 1982]](#ref-tapponnier1982).
The dense non-cratonic lower crust and mantle root delaminate and sink
[[Kay & Kay, 1993]](#ref-kay1993). Cratonic lithosphere is buoyant and stiff
[[Jordan, 1978]](#ref-jordan1978) [[Sleep, 2005]](#ref-sleep2005), and collisional stacking
is one of the ways cratonic mantle is *assembled* rather than destroyed
[[Lee et al., 2011]](#ref-lee2011) [[Capitanio et al., 2020]](#ref-capitanio2020). Nothing
in that picture removes tens of percent of a craton in one event.

1. **Give the merge's overflow the suture's outlets.** `quad_merge._place_stacking_overflow`
   calls `_place_suture_crust` with `strike=None` and `root_capacity=None`, which skips
   two of its five stages: tectonic escape and root delamination. Passing the merged
   plate's root capacity lets non-cratonic roots shed to make room, which matches
   [[Kay & Kay, 1993]](#ref-kay1993). Passing a strike (the overlap's long axis) enables
   lateral escape, which matches [[Tapponnier et al., 1982]](#ref-tapponnier1982).
2. **Place craton overflow over the whole merged plate before clipping any.** If the
   staged placement still leaves craton unplaced, run a second pass over *all* continental
   cells of the merged plate with room, nearest first, whatever their hop distance.
   Clipping should only happen when the merged plate's continental room is truly
   exhausted.
3. **Stop craton-on-craton overlap from building.** Measure the overlapping craton volume
   per plate pair. When the overlap is craton on craton, the brake has failed upstream.
   Either the forced merge (`FORCED_MERGE_OVERLAP_FRACTION`, 30% for 30 Myr) should fire
   *before* the overlap is that deep for craton pairs, or retreat of the non-cratonic crust
   between them should take it.
4. **Book honestly.** Merge clips go to a new craton account `merge_capped_m3`, kept apart
   from `delaminated_m3`, so the ledger says what happened.

**Recommendation:** first reproduce the event to choose between H1–H3. Then do 1, 2 and 4
together, since they are local to the merge and conserve craton by construction. Do 3 only
if the instrumentation shows deep craton-on-craton overlaps are common.

## Design

**Units.** Craton and Hc are thicknesses in metres per cell; volumes are m³ (reported in
M km³). The cap is `SUTURE_ACCRETION_MAX_HC_M` = 84,000 m, and room per cell is
`(cap − Hc) × area`. No per-step or per-node factors are added.

**Reproduction.** `bin/debug/probe_craton_clips.py` already replays a world between two
times and captures the inputs of the first fusion that leaves craton outside its column, as
`(keep, absorb, other_points)`.
- Regenerate seed 456862899 with the generation settings stored in the save, and step to
  3,040 My. Write a snapshot there with `--snapshot`.
- Probe to 3,060 My with `--capture`.
- Replay the captured merge alone, printing per-stage placement volumes, the eligible cell
  count and area, and the room under the cap within reach of the front.

The first leg is about 3 Gyr of stepping. Its wall time must be measured; replays are
RAM-heavy, so run them one at a time.

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
- The fallback pass (2) reuses `_place_suture_crust`'s far-field stage, with `front` set to
  the cells still holding unplaced craton and `eligible` set to every continental cell.

**Runtime.** Merges happen tens of times per Gyr, not every step. The added cost is one
extra `_place_suture_crust` call per merge: a sparse BFS over the merged plate's adjacency,
O(cells × hops), about 10⁵ cells at density 4. That is under a second per merge and nothing
per ordinary step.

**Save compatibility.** No new per-node field. `merge_capped_m3` is a new key in
`World.craton_ledger`, and `cratons.ensure_ledger` already backfills missing keys at 0.0 on
load. Old saves keep their historical merge losses in `delaminated_m3`; they can't be split
after the fact, and the docs will say so.

**Knobs.** None new. Two existing constants gain a second use in the merge, and neither is
retuned here:
- `SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION`: how much of a placement may come from
  root shedding.
- `SUTURE_ESCAPE_MAX_ANGLE_DEG`: the escape cone.

If solution 3 is needed, its threshold (the overlap's craton volume per pair, in m³) is
tuned on the replays below: low enough that no merge finds craton on craton above it,
without changing merge counts by more than 10%.

## Detailed design

- `backend/app/quad_merge.py`
  - `_place_stacking_overflow`: pass `root_capacity` and `strike` to
    `quad_tectonics._place_suture_crust`; add the fallback far-field pass for craton still
    unplaced; return the shed volume per cell for booking.
  - `_transfer` (which calls `_place_stacking_overflow`) and `merge_with` in
    `sparse_quad_patch.py`: thread the step length through, so root capacity can be
    computed. Book the shed roots through `quad_tectonics._book_shed_roots` (the helper
    the suture path uses: `delaminated_lower_crust_m3` and orogeny's
    `restite_delaminated_m3`).
- `backend/app/merge_split.py`: `_merge_plates_with_craton_audit` settles to
  `merge_capped_m3` instead of `delaminated_m3`. For solution 3, record the craton-on-craton
  overlap volume in `update_overlap_tracking`.
- `backend/app/cratons.py`: add `merge_capped_m3` to `CRATON_SINK_ACCOUNTS`, with a
  comment citing this RFC.
- `backend/app/stats.py` and the frontend ledger display, if it lists craton accounts by
  name: show the new account.
- `bin/debug/probe_craton_clips.py`: print the per-stage placement breakdown for a captured
  merge.
- Tests:
  - `backend/unit_tests/test_quad_merge.py`: two plates whose overlap stacks 2 × 40 km of
    craton over a patch, with room elsewhere on the merged plate. After the merge, live
    craton plus `collision_reworked_m3` equals the craton before, and `merge_capped_m3` is
    0. A second case with no room anywhere books exactly the excess to `merge_capped_m3`.
  - `backend/unit_tests/test_cratons.py`: the ledger balances with the new account, and a
    save without the key loads.
- Docs: `docs/simulation-model.md`, the Cratons section and the merge section of "Merge and
  split".

## Acceptance criteria

- [ ] The 3,052 My event is reproduced and its cause (H1, H2, H3 or another) is written up
      in `analysis/` with the per-stage numbers.
- [ ] On a 3.5 Gyr replay of seed 456862899 and 500 Myr runs of seeds 7, 11 and 2024
      (`bin/debug/compare_cratons.py`), merges book less than 1% of the initial craton
      volume in total, and no single merge books more than 0.5% of the live craton.
- [ ] The craton ledger balance error stays within 1e-9 of the live volume
      (`cratons.balance_error_m3`), and the continental-material ledger still closes
      (`continental_ledger.assert_closed`).
- [ ] Merge counts and plate counts on those seeds change by less than 10%.
- [ ] The unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.
- [ ] Saves from before the change load, with `merge_capped_m3` = 0.

## References

- <a id="ref-capitanio2020"></a>Capitanio, F. A., Nebel, O., & Cawood, P. A. (2020).
  Thermochemical lithosphere differentiation and the origin of cratonic mantle. *Nature*,
  588, 89–94. [doi:10.1038/s41586-020-2976-3](https://doi.org/10.1038/s41586-020-2976-3)
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
