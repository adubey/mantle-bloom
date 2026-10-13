# RFC 0004: Proto-cratons: restarting craton formation

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #274 (craton formation), #314 (interior shortening); RFCs
  [0001](0001-merge-cap-craton-loss.md), [0003](0003-rifting-around-cratons.md),
  [0006](0006-intra-oceanic-arcs.md)

## Problem

Today a craton forms in one way only (`cratons.update`). A column must:
- be genuine continental crust, with Hc in [`CRATON_HOST_MIN_HC_M`,
  `CRATON_HOST_MAX_HC_M`] = [28 km, 50 km];
- sit at least `CRATON_FORMATION_MARGIN_KM` = 500 km inside its plate's edge and any
  non-genuine cell;
- stay that way for `CRATON_FORMATION_MYR` = 200 Myr.

Its whole Hc then becomes cratonic. In practice only the deep interiors of large
continents ever qualify, and those are mostly craton already.

In the run behind RFCs 0001–0007, **craton formation stopped for good at about 500 Myr**.
The run is seed 456862899, 19.2 Gyr, saved locally as
`mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine and
**not available for download**. Numbers come from `bin/debug/ledger_history.py` on that
save.

- `formed_m3` reached 357 M km³ by 500 Myr (201 M km³ of it in one sample at 202 My, the
  first time the 200 Myr clock completed). It then stayed at exactly that value for the
  remaining 18.7 Gyr.
- By 500 Myr cratons made up 73% of genuine continental area, and 80–88% for most of the
  rest of the run. The non-cratonic crust that remained was all near margins.

With no source, every craton sink in RFCs 0001–0003 is a one-way drain.

Separately, "quiet" today means only the Hc window plus distance from margins. Interior
shortening (#314), intraplate faults and failed rifts can all deform a column that passes
both tests, so a column can count as quiet while it is being deformed.

The project owner proposed a *proto-craton* stage: crust that can mature into craton, with
continental collisions as a main route.

The run is a stress test of the model's long-term behaviour, not a calibration interval:
it shows that formation shuts off. It doesn't give a formation rate to fit.

## Proposed solutions

**What happens in nature.**
- Cratons are continental cores underlain by thick (about 200–250 km), cold and chemically
  depleted mantle keels [[Pearson et al., 2021]](#ref-pearson2021). Depletion makes the
  keel buoyant enough to offset its cold, dense temperature, which is why it lasts
  [[Jordan, 1978]](#ref-jordan1978).
- Those keels formed mostly in the Archean, by high-degree melting of a hotter mantle
  [[Herzberg & Rudnick, 2012]](#ref-herzberg2012). Lithospheric mantle formed later is
  progressively less depleted [[Griffin et al., 2009]](#ref-griffin2009).
- How keels became thick and stable is debated, and two modelled mechanisms bracket it:
  - **In extension.** [[Capitanio et al., 2020]](#ref-capitanio2020) model large-scale
    decompression melting of a hot early mantle beneath *stretching* lithosphere. Melt
    extraction leaves depleted mantle at depth. Being dehydrated and stiffer, it resists
    deformation, forcing strain to migrate and the depleted mantle to cool into the
    lithosphere during long-lived, diffuse extension.
  - **In compression.** [[Wang et al., 2018]](#ref-wang2018) start from *previously
    depleted* mantle and model tectonic shortening followed by gravitational
    self-thickening. The root stabilises because its compositional buoyancy, rapid cooling
    after shortening, and secular cooling of the mantle prevent it from collapsing.
  - Both require depletion. Neither makes a craton out of thickening alone.
    [[Lee et al., 2011]](#ref-lee2011) review these and other ways cratonic mantle is built
    and destroyed.
- Continents also grow outward: accretionary and collisional belts weld onto older cores
  and become part of the stable shield once they stop deforming
  [[Cawood et al., 2009]](#ref-cawood2009). Laurentia's Proterozoic belts around its
  Archean nuclei are the classic case.

So the proposal is sound, with three caveats:
- A freshly thickened orogen is hot and weak. It collapses, delaminates or rifts before it
  cools, so stabilisation follows the *end* of collision, not the collision itself.
- After the Archean, new cratons should be rarer and less depleted
  [[Griffin et al., 2009]](#ref-griffin2009).
- Collision makes a craton only from mantle that is already depleted. The model doesn't
  track depletion, and its proxies would be wrong in sign (see 2).

1. **Welded growth.** A quiet genuine-continental column within `CRATON_WELD_REACH_KM` of
   existing craton keeps a formation clock running even when it is closer than 500 km to
   the plate edge. It needs only `CRATON_WELD_MARGIN_KM` of interior. Cratons then grow at
   their edges by absorbing the belts welded onto them
   [[Cawood et al., 2009]](#ref-cawood2009).
2. **Collisional keel (path C): deferred to a follow-up RFC.** A collision far from any craton
   should be able to seed a new one, but only by shortening already-depleted mantle
   [[Wang et al., 2018]](#ref-wang2018). That needs a model of the lithospheric mantle's
   depletion, which the model lacks, and simple proxies get it wrong:
   - Hm thickness alone isn't a keel. Hm thickens under any convergent shortening
     (`rheology.apply_convergent_deformation`) and doesn't relax back, so a
     thickness-only rule would cratonise ordinary collision belts.
   - Counting juvenile crust as melt extracted from the column's own root is wrong in
     *reservoir and sign*. Arc magma comes largely from the mantle wedge, and plume and rift
     magma from upwelling mantle below the existing lid. Melt passing through a root can
     refertilise, hydrate and weaken it rather than deplete it
     [[Lee et al., 2011]](#ref-lee2011) [[Foley, 2008]](#ref-foley2008).

   The follow-up RFC must model, with transport and removal rules:
   - melt extracted from the resident lithospheric mantle (depletes it);
   - melt delivered from elsewhere and passing through (refertilises it);
   - depleted residue incorporated into the lithosphere from below, as in the extension
     mechanism.

   This RFC delivers paths A and B and the quiescence signal. It leaves path C off and
   undesigned.
3. **Name the stage.** A *proto-craton* is a non-cratonic column whose formation clock is
   running on path A or B. It gets no extra mechanical resistance: young belts are weak
   (see RFC 0003). It is reported in stats and shown in the existing `craton` debug view,
   so formation can be watched.
4. **A quiescence signal that sees every deforming writer.** A column's clock runs only
   while no tectonic process has disturbed it (Design).
5. **Secular decline (optional).** Scale the cratonic share given at formation by
   `CRATON_FORMATION_DEPLETION`, to mimic less-depleted post-Archean keels. The model
   seeds its starting cratons 2.5 Gyr old, so the simulated era is post-Archean.

**Recommendation:** 1, 3 and 4 now. Path 1 alone restarts formation on the margins of
existing cratons, and 4 stops any path from cratonising a column that is being deformed.
2 waits for its own RFC. Leave 5 at 1.0 (off) until a run shows cratons ratcheting up to
most of the continental area.

## Design

**Quantities and units.**
- Hc, Hm and craton thicknesses in metres.
- Distances in km, converted to graph hops with `cratons._hops` (62.5 km per hop at
  density 4, 125 km at density 1).
- The formation clock `stable_continental_myr` in Myr; it already exists and is a CLOCK
  field.

**Formation rule.** `quiet` currently means `host & (Hc ≤ CRATON_HOST_MAX_HC_M) &
(interior ≥ hops(500 km))`. It becomes
`host & (Hc ≤ CRATON_HOST_MAX_HC_M) & undisturbed & (A | B)`:

- **A (interior):** `interior ≥ hops(CRATON_FORMATION_MARGIN_KM)`, as today.
- **B (welded):** `interior ≥ hops(CRATON_WELD_MARGIN_KM)` and
  `hops_to_craton ≤ hops(CRATON_WELD_REACH_KM)`. Here `hops_to_craton` is
  `quad_tectonics.hop_distance(plate, craton > 0, max_hops)`.

The interior test still excludes plate edges and non-genuine cells, so an active margin
never qualifies. The 200 Myr clock and the Hc window are unchanged: a belt must relax back
to ≤ 50 km of crust and stay quiet for 200 Myr, so stabilisation follows the end of the
orogeny. On formation the craton is `CRATON_FORMATION_DEPLETION × Hc`, 1.0 by default.

**Quiescence (4).** A step-scoped Hc delta from `deform_columns` isn't enough. In
`world.step_world_progress`, faults (`faults.update_faults`) and topology changes, which
include failed rifts and splits, run after the deform loop and before `cratons.update`.
And transform deformation writes elevation, not Hc. So the signal is a persistent per-node
field, set where deformation happens:
- **The field.** `tectonic_disturbance_myr` is a COUNTDOWN surface field. COUNTDOWN
  remaps take the maximum when cells combine in `quad_merge`, so a merged cell is
  disturbed if any part of it was. The refinement and coarsening remaps in
  `sparse_quad_patch` must be checked to do the same, and a test covers it. Cells grown at an edge get their source's value through
  `insert_cells`' field inheritance.
- **Marking.** Each writer that deforms a column calls one helper,
  `cratons.mark_disturbed(plate, mask)`. It sets the field to
  `PROTO_CRATON_DISTURBANCE_HOLD_MYR`, which is at least one maximum step, so the mark
  survives to the step's `cratons.update`. The writers:
  - `lithosphere_plate.deform_columns`: convergent, divergent and arc thickening or
    thinning with |ΔHc| above `PROTO_CRATON_MAX_STRAIN_M_PER_MYR × years`, and every
    transform-band cell with pressure-ridge uplift;
  - `faults.update_faults`: every cell a fault's relief moved;
  - `sparse_quad_patch.apply_failed_rift`: the whole thinning band;
  - `merge_split.maybe_split_plate`: cells within the failed-rift band of a successful
    split's cut;
  - `quad_tectonics`: suture-accretion receivers (`_accrete_onto_survivors`), stretch
    donors (`_allocate_stretch`) and ceiling-overflow receivers;
  - `crust_transfer`: upper-plate receivers;
  - `quad_merge`: stacked cells and overflow receivers;
  - `orogeny.relax_orogens`: cells whose Hc it changes by more than the strain threshold.
- **Clock.** In `cratons.update`, a column with `tectonic_disturbance_myr > 0` is not
  quiet, so its clock resets to 0, as leaving the window does now. The countdown then
  decreases by the step length.
- **Audit for missed writers.** A writer-agnostic check catches writers missing from the
  list. In debug runs, `world.step_world_progress` snapshots each cell's Hc and elevation
  before the deform loop, and compares them, matched by cell key, after the topology and
  gap-fill phase. That is after every tectonic writer and before erosion and volcanism,
  so no surface change needs subtracting. Cells whose keys changed through remap are
  skipped. Any cell with a large change and no disturbance mark is counted and reported.
  The acceptance runs require that count to be 0.

**Runtime.** Path B adds one bounded BFS per plate per step, out to
`hops(CRATON_WELD_REACH_KM)` = 8 hops at density 4. That is the same kind and size of work
as the interior BFS `cratons.update` already does. Marking is a masked write in each
writer, and the countdown one vector update per plate. The debug-only audit keeps one
Hc/elevation snapshot per plate. Expected cost: under 1% of step time, which acceptance
requires measuring.

**Save compatibility.** `tectonic_disturbance_myr` is a new COUNTDOWN surface field,
backfilled at 0 for old saves (undisturbed). Old saves' clocks keep running, as they do
today. Proto-craton status is derived each step. The new stats keys are additive.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `CRATON_WELD_REACH_KM` | 500 km | One belt width: orogens welded onto shields are a few hundred km wide [[Cawood et al., 2009]](#ref-cawood2009) | Calibration assumption; sweep 250–1,000 km |
| `CRATON_WELD_MARGIN_KM` | 250 km | Design choice: half the interior margin keeps active margins and their arcs out | Sweep 125–500 km |
| `PROTO_CRATON_MAX_STRAIN_M_PER_MYR` | 50 m/Myr of tectonic Hc change | Design choice: about 1% of a column's crust per 10 Myr, well below active orogens (hundreds of m/Myr) | Sweep 20–200 |
| `PROTO_CRATON_DISTURBANCE_HOLD_MYR` | 10 Myr | Numerical: at least the longest step the UI offers, so a mark always reaches `cratons.update` | Fixed |
| `CRATON_FORMATION_DEPLETION` | 1.0 | Off; see solution 5 | Lower only if the ratchet criterion fails |
| `CRATON_FORMATION_MYR` | 200 Myr (unchanged) | #274 | Unchanged |

Tune with `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 2 Gyr (formation needs
time to show), plus a regenerated seed 456862899 to 3 Gyr; run the replays one at a time.
Report `formed_m3` per Gyr by path (A, B), craton share of genuine continental area,
proto-craton area, and how often each writer marks cells.

## Detailed design (provisional touchpoints)

- `backend/app/cratons.py`
  - New constants as above, with comments citing this RFC.
  - `update`: compute paths A and B, require `tectonic_disturbance_myr == 0`, combine them
    into `quiet`, count the countdown down, and apply `CRATON_FORMATION_DEPLETION` at
    formation.
  - New `mark_disturbed(plate, mask)`.
  - Book `formed_m3` as now, plus per-path diagnostic counters kept on the ledger as
    `formed_interior_m3` and `formed_welded_m3`. These are sub-accounts of `formed_m3`,
    outside the balance.
  - New `proto_craton_mask(plate)`, shared by `update`, `diagnostics` and rendering.
  - `diagnostics`: add `proto_craton_area_m2`, the per-path formed volumes, and the
    debug audit's unmarked-change count.
  - Update the module docstring's Formation paragraph.
- `backend/app/surface_fields.py`: register `tectonic_disturbance_myr` (new, COUNTDOWN).
- Writers calling `mark_disturbed`: `lithosphere_plate.deform_columns`,
  `faults.update_faults`, `sparse_quad_patch.apply_failed_rift`,
  `merge_split.maybe_split_plate`, `quad_tectonics` (`_accrete_onto_survivors`,
  `_allocate_stretch`, ceiling-overflow placement), `crust_transfer`, `quad_merge`, and
  `orogeny.relax_orogens`.
- `backend/app/world.py`: in debug runs, snapshot Hc and elevation by cell key before the
  deform loop, and run the audit after topology and gap fill.
- `backend/app/stats.py`: `proto_craton_area_km2` in the stats row.
- `backend/app/render_image.py`: the `craton` debug view (`_render_craton_view`) already
  colours quiet interiors by their clock's progress. It takes its mask from
  `proto_craton_mask`, so welded proto-cratons show too, and disturbed cells show as not
  forming. `frontend/src/legendData.ts`'s `craton` legend relabels "Quiet interior near
  cratonising" as "Proto-craton near cratonising".
- `bin/debug/compare_cratons.py` and `bin/debug/ledger_history.py`: report per-path
  formation.
- Tests (`backend/unit_tests/test_cratons.py` unless noted):
  - a quiet belt 300 km from the edge, next to a craton, forms after 200 Myr (path B), and
    the same belt with no craton nearby does not;
  - **near-complete clock, failed rift:** a non-cratonic column at 35 km with its clock at
    199.9 Myr is thinned 10% by `apply_failed_rift` in the same step, stays inside the
    28–50 km window, and does *not* cratonise; its clock resets;
  - the same for a fault's relief, a transform pressure ridge (elevation only), a merge's
    stacking and an orogenic-flow change;
  - a disturbed cell merged with undisturbed ones stays disturbed (COUNTDOWN max);
  - an active margin never forms;
  - each writer in the list marks the cells it deforms (a contract test per writer, in
    that writer's test module);
  - ledger balance is unchanged by the sub-accounts.
- Docs: `docs/simulation-model.md`, the Cratons section ("Formation"), and the render-view
  description if the craton view changes.

## Assumptions and uncertainty

- **Formation may not be the bottleneck.** Formation stopped, but losses (RFCs 0001–0003)
  dominate. Formation is complementary, not a cure.
- **Welded growth is an analogy.** Laurentia-style growth is observed over Gyr, but the
  weld reach, margin and 200 Myr clock are model choices, not measurements. Path B also
  ignores depletion: it assumes belts welded to a craton inherit its stability, which is
  only partly true geologically.
- **The writer list may be incomplete.** The audit is there to catch what the list
  misses; the acceptance criteria require it to report nothing.
- **Secular change is ignored.** Post-Archean keels are less depleted
  [[Griffin et al., 2009]](#ref-griffin2009), so the simulated era should form fewer and
  weaker cratons. Knob 5 is available but off.
- **The ratchet band is a design target.** 30–65% of continental area as craton is the
  range the #274 seeding was calibrated to (shields plus platforms), not a measured Earth
  value for every epoch.

## Acceptance criteria

- [ ] On seeds 7, 11 and 2024 to 2 Gyr, and seed 456862899 to 3 Gyr, `formed_m3` increases
      in every 500 Myr window in which genuine continental area is at least 25% of its
      initial value. Today it stops after 500 Myr.
- [ ] Path B contributes in at least two of the four runs.
- [ ] No column forms craton while being deformed: the unit tests above, plus the debug
      audit reporting zero unmarked tectonic changes on the acceptance runs.
- [ ] No ratchet: craton share of genuine continental area stays between 30% and 65% at
      every 500 Myr checkpoint (it starts at 35%).
- [ ] Step time changes by less than 2% at density 4.
- [ ] Unit tests above pass; ledgers balance; old saves load; affected tests, ruff and
      `bin/quality_check.sh` are green.

## Falsification criteria

- **Formation doesn't matter.** If restarting formation (paths A and B) changes craton
  volume at 3 Gyr on seed 456862899 by less than 10%, formation isn't what limits craton
  survival. Deprioritise this RFC behind 0001–0003.
- **Welded growth is a ratchet.** If path B pushes the craton share above 65% on most
  seeds at any reasonable weld reach (250–1,000 km), welding everything near a craton is
  too generous. Restrict B to belts with a collisional history (an
  `ELEV_CHANGE_COLLISION` stamp or a suture within the reach), or drop it.
- **Quiescence starves formation.** If, with the disturbance signal, almost no column
  stays undisturbed for 200 Myr (path A's formation drops by more than 90% against
  `origin/main`), the model deforms continental interiors far more than the clock
  assumes. Investigate the dominant marking writer before tuning the threshold.

## References

- <a id="ref-capitanio2020"></a>Capitanio, F. A., Nebel, O., & Cawood, P. A. (2020).
  Thermochemical lithosphere differentiation and the origin of cratonic mantle. *Nature*,
  588, 89–94. [doi:10.1038/s41586-020-2976-3](https://doi.org/10.1038/s41586-020-2976-3)
- <a id="ref-cawood2009"></a>Cawood, P. A., Kröner, A., Collins, W. J., Kusky, T. M.,
  Mooney, W. D., & Windley, B. F. (2009). Accretionary orogens through Earth history.
  *Geological Society, London, Special Publications*, 318, 1–36.
  [doi:10.1144/SP318.1](https://doi.org/10.1144/SP318.1)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-griffin2009"></a>Griffin, W. L., O'Reilly, S. Y., Afonso, J. C., & Begg,
  G. C. (2009). The composition and evolution of lithospheric mantle: A re-evaluation and
  its tectonic implications. *Journal of Petrology*, 50, 1185–1204.
  [doi:10.1093/petrology/egn033](https://doi.org/10.1093/petrology/egn033)
- <a id="ref-herzberg2012"></a>Herzberg, C., & Rudnick, R. (2012). Formation of cratonic
  lithosphere: An integrated thermal and petrological model. *Lithos*, 149, 4–15.
  [doi:10.1016/j.lithos.2012.01.010](https://doi.org/10.1016/j.lithos.2012.01.010)
- <a id="ref-jordan1978"></a>Jordan, T. H. (1978). Composition and development of the
  continental tectosphere. *Nature*, 274, 544–548.
  [doi:10.1038/274544a0](https://doi.org/10.1038/274544a0)
- <a id="ref-lee2011"></a>Lee, C.-T. A., Luffi, P., & Chin, E. J. (2011). Building and
  destroying continental mantle. *Annual Review of Earth and Planetary Sciences*, 39,
  59–90. [doi:10.1146/annurev-earth-040610-133505](https://doi.org/10.1146/annurev-earth-040610-133505)
- <a id="ref-pearson2021"></a>Pearson, D. G., Scott, J. M., Liu, J., Schaeffer, A., Wang,
  L. H., van Hunen, J., Szilas, K., Chacko, T., & Kelemen, P. B. (2021). Deep continental
  roots and cratons. *Nature*, 596, 199–210.
  [doi:10.1038/s41586-021-03600-5](https://doi.org/10.1038/s41586-021-03600-5)
- <a id="ref-wang2018"></a>Wang, H., van Hunen, J., & Pearson, D. G. (2018). Making
  Archean cratonic roots by lateral compression: A two-stage thickening and stabilization
  model. *Tectonophysics*, 746, 562–571.
  [doi:10.1016/j.tecto.2016.12.001](https://doi.org/10.1016/j.tecto.2016.12.001)
