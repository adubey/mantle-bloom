# RFC 0009: Cratons as a mosaic: discrete cores and lasting weak belts

- **Status:** Accepted (2026-10-15)
- **Date:** 2026-10-15
- **Related:** #427 (tracking), #426 (gate G1), #274 (cratons), #389 (RFC 0003 tracking), #409 (formation and ratchet
  diagnosis), #411; PR #418 (RFC 0003 solution 1); RFCs
  [0003](0003-rifting-around-cratons.md) and [0004](0004-proto-cratons.md) (both
  superseded by this RFC and [0010](0010-rifts-follow-weak-lithosphere.md)),
  [0007](0007-mantle-plumes-and-plateaus.md), [0008](0008-craton-collision-jamming.md)

## Problem

The model defines a craton by distance from its plate's edge, so every large continent's
interior is a single craton with no internal structure.

**How cratons are made today** (`cratons.py`, origin/main `1227c36`):
- **Seeding** (`seed_initial_cratons`): at generation, every genuine continental column with
  Hc in [28 km, 50 km] that lies at least `CRATON_SEED_MARGIN_KM` = 875 km inside its plate
  and any non-genuine cell becomes craton, dated 2.5 Gyr old.
- **Formation** (`update`): any column meeting the same Hc window at least
  `CRATON_FORMATION_MARGIN_KM` = 500 km inside, and staying so for `CRATON_FORMATION_MYR` =
  200 Myr, becomes craton.

Nothing else decides where a craton is. Three consequences follow.

1. **One craton per continent.** Seeding fills each continent's deep interior with one
   contiguous craton. At `2562f63`, cratons start at 37.6–50.7% of genuine continental area
   on seeds 7, 11 and 2024 (#403). On seed 7 at origin/main, craton covers 104 M km², 20% of
   the planet's surface, at t = 0 (16% and 9% on seeds 11 and 2024). Its largest plate
   (92 M km²) carries 54 M km² (`bin/debug/prototype_mincut_split.py --coverage-only`,
   recorded in [analysis/rfc0010](../analysis/rfc0010/README.md)).
2. **Welds and scars heal into craton.** A collision suture or a failed-rift band is only
   excluded from formation while its Hc is outside the window or it is near an edge. Once
   it has been quiet inside the window for 200 Myr, formation converts it like any other
   interior. The model keeps no record of where continents were joined or stretched.
3. **Continents become all craton.** Cratons take 85–92% of genuine continental area by
   500 Myr and 94–97% by 1 Gyr (#403). #409 showed this is survival bias, not formation:
   at least 91% of the non-cratonic host lost was destroyed or thinned away, not converted.
   The remainder sits near margins.

**Why it matters: rifting.** A plate splits along a cut placed by the mantle flow. RFC
0003 tried to steer that cut off cratons, first by tilting the great circle (PR #418), then,
in an offline probe, by letting it bend around strong lithosphere along a graph minimum cut
(see [RFC 0010](0010-rifts-follow-weak-lithosphere.md) and
[analysis/rfc0010](../analysis/rfc0010/README.md)). Both mostly fail for the same reason:
a tilted great circle found an admissible craton-free cut for 12% of craton-crossing cuts,
and a bent minimum cut for 10% at a 1,000 km corridor and at most 12% at any corridor up
to 2,500 km. Wide corridors find craton-free routes only by slicing a coastal sliver off
the plate. On craton-crossing cuts the median plate is 34% craton, and craton fills a
median 44% of the cut's band. There is no route through a continent that avoids its craton, because the
craton *is* the continent's interior. Steering can't help until interiors have seams.

**Why it matters: formation.** RFC 0004 found formation shutting off by 335–470 Myr, and
#409 traced 75–86% of each run's formation to one generation pulse: seeding uses an 875 km
margin and formation a 500 km one, so the band between matures together at 200 Myr. The
interior-distance rule is the cause in both cases.

What this evidence does *not* show: whether a different craton layout would survive the
model's sinks (RFCs 0001, 0002, 0005, 0008), or whether flow-placed cuts would then miss
cratons. The diagnosis gate below measures the second.

## Proposed solutions

**What happens in nature.**
- A craton is a block of continental crust over a thick (about 200–250 km), cold,
  chemically depleted mantle keel [[Pearson et al., 2021]](#ref-pearson2021)
  [[Jordan, 1978]](#ref-jordan1978). Most keels formed in the Archean, by high-degree
  melting of a hotter mantle [[Herzberg & Rudnick, 2012]](#ref-herzberg2012); later
  lithospheric mantle is progressively less depleted [[Griffin et al., 2009]](#ref-griffin2009).
- Cratons are discrete. The late Archean record survives as roughly 35 separate cratonic
  pieces [[Bleeker, 2003]](#ref-bleeker2003). Continents are assembled from such pieces:
  Laurentia, the North American craton, was built in the Early Proterozoic by joining older
  provinces along orogenic belts [[Hoffman, 1988]](#ref-hoffman1988). Accretionary and
  collisional belts weld onto older cores [[Cawood et al., 2009]](#ref-cawood2009).
- The architecture persists at depth. Africa's lithosphere maps into domains of ancient,
  highly depleted mantle and zones of younger or strongly modified mantle
  [[Begg et al., 2009]](#ref-begg2009).
- Belts stay weak. Rifts run parallel to old orogenic belts, which is explained by a
  mechanical anisotropy frozen into the lithospheric mantle when continents were assembled
  [[Vauchez et al., 1997]](#ref-vauchez1997) [[Tommasi & Vauchez, 2001]](#ref-tommasi2001),
  or by those belts' thinner lithosphere [[Audet & Bürgmann, 2011]](#ref-audet2011). In
  numerical models, structures frozen into the mantle lithosphere by earlier suturing act as
  latent plate boundaries that later forcing reactivates
  [[Heron et al., 2016]](#ref-heron2016). The Atlantic reopened along the
  Caledonian–Appalachian suture [[Wilson, 1966]](#ref-wilson1966).

So a continent should be a mosaic: several cratonic cores, separated and surrounded by
belts that stay weaker than the cores for a long time.

1. **Seed a mosaic.** At generation, partition each continent's eligible interior into
   discrete cratonic cores separated by belts, instead of one interior-distance blob.
   - Deterministic per seed; changes only generation.
   - Also removes #409's generation pulse, if seeding and formation use the same margin.
2. **Mark belts, and keep them.** A new per-node field records where crust was welded or
   stretched: the inter-core belts at generation, then collision sutures and rift bands as
   they happen. A marked column never becomes craton. RFC 0010 uses the same field as a
   low-strength path for rifts.
3. **Welded growth (from RFC 0004, path B).** Quiet columns next to a craton keep a
   formation clock even near the plate edge, so cores grow at their rims. Belt marks still
   exclude columns, so cores grow without closing their seams.
4. **Quiescence signal (from RFC 0004).** A formation clock runs only while nothing deforms
   the column. Sound, but a large change touching every tectonic writer; it doesn't
   affect where cratons are, only when they form.
5. **Craton identity.** A per-node core id, so two cores that collide stay two cores and
   their join is a seam by construction. More state; belt marks (2) give most of the same
   effect.
6. **Depletion-based formation (RFC 0004's deferred path C).** Form cratons only from
   depleted mantle. Still needs a depletion model the code doesn't have; stays deferred.

**Recommendation:** 1 and 2 together, gated by the diagnosis below. Keep 3 and 4 as later
options once the mosaic survives; leave 5 and 6 out.

## Design

**Diagnosis gate G1: does a mosaic give rifts a way around?** A replay-only probe,
`bin/debug/probe_craton_mosaic.py` (new), reseeds generation as a mosaic (solution 1, with
the starting knob values below) and sets `CRATON_FORMATION_MYR` to infinity in-process, so
the seeded layout is the only craton for the run. No saved state changes. On seeds 7, 11
and 2024, 5 Myr steps to 1 Gyr, default density and controls, one at a time, it runs the
min-cut probe of RFC 0010 on every cut. Over the flow cuts whose band crosses craton, at a
1,000 km corridor and with all of RFC 0010's admissibility gates (size, pole separation,
residual, connectivity), it reports two branch measures:
- **L, line clearance:** the share of those cuts with an admissible minimum cut whose line
  touches no craton (no craton node at either end of a cut link). This says whether a
  route between the cores exists.
- **B, band exposure:** craton area in the operational band (`FAILED_RIFT_BAND_MULT` = 2.5
  spacings either side, the band a failed rift thins and `record_rift_cut` books), summed
  over those cuts, using the minimum cut where it is admissible and the flow cut where it
  isn't, divided by the same sum for the flow cuts. This says how much craton a rift would
  still act on.

Both are needed. A 300 km belt is narrower than the 625 km band, so a cut can run along a
belt with a clear line while its band still thins craton on both sides; line clearance
alone doesn't show that the craton a rift acts on falls.

It also reports belt area, as a share of genuine continental area, at each 100 Myr
checkpoint, and the same L and B on today's layout as a control
([analysis/rfc0010](../analysis/rfc0010/README.md)).

Outcomes, fixed before G1 runs:
- **L ≥ 50% and B ≤ 0.50:** build 1 and 2 here, and RFC 0010's minimum cut.
- **L ≥ 20% and B ≤ 0.75**, but not the branch above: build 1 and 2; RFC 0010 builds its
  minimum cut together with its resistance fallback, sized to the remaining loss.
- **Otherwise:** the falsification case "the flow puts regime boundaries through cratons"
  holds even for a mosaic (Falsification). Don't build RFC 0010's cut; take the question
  to the mantle model.

B bounds what failed rifts and the booked exposure see. Stretch donation and divergent
thinning along a new boundary act after the split and can reach beyond the band (#388), and
the probe doesn't follow them; RFC 0010's acceptance runs measure that loss directly.

Report also how fast belts disappear. If belt area falls below a third of its starting
value by 500 Myr, the mosaic won't outlast the model's sinks, and this RFC waits for the
sink repairs (RFCs 0001, 0002, 0005).

**Mosaic seeding (1).** For each plate, take the columns seeding uses today, but with the
formation margin, so the two agree and the 200 Myr pulse disappears:
`host & (Hc ≤ CRATON_HOST_MAX_HC_M) & (interior ≥ hops(CRATON_FORMATION_MARGIN_KM))`.
- Place core nuclei on these columns by Poisson-disc sampling with minimum separation
  `CRATON_CORE_SPACING_KM`, from the world's generation RNG.
- Assign each eligible column to its nearest nucleus by graph hops
  (`quad_tectonics.hop_distance` from each nucleus; a multi-source BFS), giving a discrete
  Voronoi partition of the interior.
- Columns within `CRATON_BELT_WIDTH_KM / 2` of a partition boundary are belts. The rest
  become craton as today (the whole Hc, dated `CRATON_SEED_AGE_YEARS`).
- An interior too small for two nuclei gets one core, so small continents keep a single
  craton, as now.

**Belt marks (2).** `belt_formed_years` (new surface field) holds the time a column last
became a belt, with a "never" sentinel. Writers:
- seeding: belts get `−CRATON_BELT_SEED_AGE_YEARS` (Proterozoic welds, younger than their
  cores);
- collision welds: suture-accretion receivers (`quad_tectonics._accrete_onto_survivors`),
  upper-plate receivers in `crust_transfer`, and stacked cells in `quad_merge`;
- rifting: a failed rift's thinned band (`sparse_quad_patch.apply_failed_rift`) and the band
  along a successful split's cut (`merge_split.maybe_split_plate`).

`cratons.update` excludes marked columns from `quiet`. `CRATON_BELT_HEAL_MYR` (default
infinite) could later let old belts heal; the Assumptions section says why it starts off.

The remap rule must keep a mark when any contributing cell has one and keep its latest
date. `craton_formed_years` is a HISTORY field that keeps the *oldest* contributor, so the
implementation checks whether HISTORY can be parameterised or a max-date class is needed,
and tests merge, refine and coarsen remaps.

**Units.** Distances in km, converted to hops with `cratons._hops`; times in years, like
`craton_formed_years`. The mark changes no thickness; it gates formation only.

**Runtime.** Seeding adds a multi-source BFS per plate at generation. Marking is a masked
write in each writer. `update` adds one mask. Expected well under 1% of step time;
acceptance measures it.

**Save compatibility.** `belt_formed_years` is a new surface field, backfilled with the
"never" sentinel. Old saves keep their existing one-blob cratons: seeding runs only at
generation (and on the first load of a pre-craton save, as now). Their belts appear only as
new collisions and rifts mark them.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `CRATON_CORE_SPACING_KM` | 1,500 km | Design choice: about 35 Archean cratons survive [[Bleeker, 2003]](#ref-bleeker2003), so cores of order 1,000–2,000 km | G1 and acceptance runs; sweep 1,000–2,500 km |
| `CRATON_BELT_WIDTH_KM` | 300 km | Orogens welded onto shields are a few hundred km wide [[Cawood et al., 2009]](#ref-cawood2009) | Sweep 200–600 km; at least two node spacings at density 1 (250 km) |
| `CRATON_BELT_SEED_AGE_YEARS` | 1.8 Gyr | Design choice: Proterozoic assembly, as for Laurentia [[Hoffman, 1988]](#ref-hoffman1988) | Fixed; provenance only |
| `CRATON_BELT_HEAL_MYR` | infinite | Design choice (Assumptions) | Unchanged in this RFC |
| `CRATON_SEED_MARGIN_KM` | retired | Replaced by `CRATON_FORMATION_MARGIN_KM`, removing the generation pulse (#409) | — |

## Detailed design (provisional touchpoints)

- `backend/app/cratons.py`
  - `seed_initial_cratons`: mosaic seeding; new constants above, with comments citing this
    RFC.
  - `update`: exclude columns with a belt mark from `quiet`.
  - New `mark_belt(world, plate, mask)`, shared by the writers.
  - `diagnostics`: `belt_area_m2` and the number of cores (connected cratonic components
    per plate).
  - Module docstring: Formation and a new Belts paragraph.
- `backend/app/surface_fields.py`: `belt_formed_years` (new), with its remap class.
- Writers calling `mark_belt`: `quad_tectonics._accrete_onto_survivors`,
  `crust_transfer.transfer_column`, `quad_merge._transfer`,
  `sparse_quad_patch.apply_failed_rift`, `merge_split.maybe_split_plate`.
- `backend/app/stats.py`: `belt_area_km2` in the stats row.
- `backend/app/render_image.py`: show belts in the `craton` debug view, with
  `frontend/src/legendData.ts`'s `craton` legend.
- `bin/debug/probe_craton_mosaic.py` (new, G1) and `bin/debug/compare_cratons.py` (belt
  area per checkpoint).
- Tests (`backend/unit_tests/test_cratons.py` unless noted):
  - mosaic seeding on a synthetic continent gives several cores separated by belts of the
    configured width, and the same cores for the same seed;
  - a small interior gets a single core;
  - a marked column in the Hc window, 1,000 km inside, quiet for 300 Myr, doesn't form;
  - each writer marks the cells it welds or stretches (a contract test in each writer's
    test module);
  - a merged or coarsened cell keeps a contributor's mark and its latest date;
  - old saves load with no belts and unchanged cratons (`test_persistence.py`).
- Docs: `docs/simulation-model.md`, Cratons (seeding, formation, belts).

## Assumptions and uncertainty

- **The gate uses no formation.** G1 freezes the seeded layout to isolate the geometry. The
  implementation keeps formation (excluding belts), which can fill non-belt interior
  between cores; the acceptance runs measure whether seams survive that.
- **Belts may not survive the sinks.** #409 found 87–96% of non-cratonic host gone by
  500 Myr. Interior belts are farther from margins than most host, but they are host, and
  the same sinks act on them. G1 reports belt survival; if it is short, the sink RFCs come
  first.
- **Belts never heal here.** Geologically, welded belts join the stable shield
  [[Cawood et al., 2009]](#ref-cawood2009), yet remain weaker and reactivate
  [[Heron et al., 2016]](#ref-heron2016) [[Vauchez et al., 1997]](#ref-vauchez1997).
  The model has one craton state, so a belt is either craton or not. Keeping belts
  non-cratonic keeps them weak for rifting, at the cost of leaving them exposed to
  erosion and collision sinks that cratons resist.
- **Core spacing and belt width are design choices.** Bleeker's count constrains how many
  Archean pieces survive, not the size of the model's cores; the sweeps report how the
  results depend on both.
- **The model still has no depletion.** Mosaic seeding places cores geometrically, not
  where a depleted keel would be.

## Acceptance criteria

- [ ] **Gate G1** is run and reported per seed in `analysis/`, with its pinned commit, and
      this RFC is updated with the branch it selects.
- [ ] At t = 0 on seeds 7, 11 and 2024, every interior large enough for two nuclei has at
      least two cores, and belts are 10–30% of the seeded interior area at the starting
      knobs.
- [ ] No generation pulse: formation in the step that first reaches 200 Myr is under 10%
      of the matched baseline's on each seed.
- [ ] Belt area at 1 Gyr is at least a third of its t = 0 value on at least two of the
      three seeds.
- [ ] No marked column ever forms craton (unit tests, and a debug-run audit counting
      formation on marked columns, which must be 0).
- [ ] Step time changes by less than 2% at density 4; generation time by less than 5%.
- [ ] Ledgers balance; old saves load; unit tests above pass; affected tests, ruff and
      `bin/quality_check.sh` are green.

## Falsification criteria

- **The flow puts regime boundaries through cratons.** If G1 falls in neither building
  branch on a mosaic (L < 20%, or B > 0.75), the cut's geometry was never the obstacle: the flow's regime boundaries sit under cratons. Don't
  build RFC 0010's minimum cut; investigate why the mantle flow localises beneath cratons.
- **The mosaic doesn't last.** If belts lose more than two thirds of their area by 500 Myr
  on every seed even without formation, the mosaic only delays the one-blob state. Wait for
  the sink RFCs (0001, 0002, 0005), then rerun G1.
- **Belts become the new sink.** If marked belts lose crust so fast that continental area
  at 1 Gyr falls by more than 10% against the matched baseline, keeping belts
  non-cratonic is too costly; give belts partial craton resistance instead.

## References

- <a id="ref-audet2011"></a>Audet, P., & Bürgmann, R. (2011). Dominant role of tectonic
  inheritance in supercontinent cycles. *Nature Geoscience*, 4, 184–187.
  [doi:10.1038/ngeo1080](https://doi.org/10.1038/ngeo1080)
- <a id="ref-begg2009"></a>Begg, G. C., Griffin, W. L., Natapov, L. M., O'Reilly, S. Y.,
  Grand, S. P., O'Neill, C. J., Hronsky, J. M. A., Djomani, Y. P., Swain, C. J., Deen, T.,
  & Bowden, P. (2009). The lithospheric architecture of Africa: Seismic tomography, mantle
  petrology, and tectonic evolution. *Geosphere*, 5(1), 23–50.
  [doi:10.1130/GES00179.1](https://doi.org/10.1130/GES00179.1)
- <a id="ref-bleeker2003"></a>Bleeker, W. (2003). The late Archean record: a puzzle in ca.
  35 pieces. *Lithos*, 71, 99–134.
  [doi:10.1016/j.lithos.2003.07.003](https://doi.org/10.1016/j.lithos.2003.07.003)
- <a id="ref-cawood2009"></a>Cawood, P. A., Kröner, A., Collins, W. J., Kusky, T. M.,
  Mooney, W. D., & Windley, B. F. (2009). Accretionary orogens through Earth history.
  *Geological Society, London, Special Publications*, 318, 1–36.
  [doi:10.1144/SP318.1](https://doi.org/10.1144/SP318.1)
- <a id="ref-griffin2009"></a>Griffin, W. L., O'Reilly, S. Y., Afonso, J. C., & Begg,
  G. C. (2009). The composition and evolution of lithospheric mantle: A re-evaluation and
  its tectonic implications. *Journal of Petrology*, 50, 1185–1204.
  [doi:10.1093/petrology/egn033](https://doi.org/10.1093/petrology/egn033)
- <a id="ref-herzberg2012"></a>Herzberg, C., & Rudnick, R. (2012). Formation of cratonic
  lithosphere: An integrated thermal and petrological model. *Lithos*, 149, 4–15.
  [doi:10.1016/j.lithos.2012.01.010](https://doi.org/10.1016/j.lithos.2012.01.010)
- <a id="ref-heron2016"></a>Heron, P. J., Pysklywec, R. N., & Stephenson, R. (2016).
  Lasting mantle scars lead to perennial plate tectonics. *Nature Communications*, 7,
  11834. [doi:10.1038/ncomms11834](https://doi.org/10.1038/ncomms11834)
- <a id="ref-hoffman1988"></a>Hoffman, P. F. (1988). United Plates of America, the birth
  of a craton: Early Proterozoic assembly and growth of Laurentia. *Annual Review of Earth
  and Planetary Sciences*, 16, 543–603.
  [doi:10.1146/annurev.ea.16.050188.002551](https://doi.org/10.1146/annurev.ea.16.050188.002551)
- <a id="ref-jordan1978"></a>Jordan, T. H. (1978). Composition and development of the
  continental tectosphere. *Nature*, 274, 544–548.
  [doi:10.1038/274544a0](https://doi.org/10.1038/274544a0)
- <a id="ref-pearson2021"></a>Pearson, D. G., Scott, J. M., Liu, J., Schaeffer, A., Wang,
  L. H., van Hunen, J., Szilas, K., Chacko, T., & Kelemen, P. B. (2021). Deep continental
  roots and cratons. *Nature*, 596, 199–210.
  [doi:10.1038/s41586-021-03600-5](https://doi.org/10.1038/s41586-021-03600-5)
- <a id="ref-tommasi2001"></a>Tommasi, A., & Vauchez, A. (2001). Continental rifting
  parallel to ancient collisional belts: An effect of the mechanical anisotropy of the
  lithospheric mantle. *Earth and Planetary Science Letters*, 185, 199–210.
  [doi:10.1016/S0012-821X(00)00350-2](https://doi.org/10.1016/S0012-821X(00)00350-2)
- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then
  re-open? *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)
