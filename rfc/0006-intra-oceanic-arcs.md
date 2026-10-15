# RFC 0006: Intra-oceanic arcs: new continent from the ocean

- **Status:** Draft
- **Date:** 2026-10-13; evidence refreshed 2026-10-14 (#403)
- **Related:** #120 (arc magmatism vs. land decline), #311 (terrane keels), #318
  (collision polarity fronts), #321 (terrane docking),
  [#406](https://github.com/adubey/mantle-bloom/issues/406) (source-attribution diagnosis); RFCs
  [0009](0009-craton-mosaic.md) (cratons; supersedes 0004), [0005](0005-continental-recycling-balance.md),
  [0007](0007-mantle-plumes-and-plateaus.md) (shared docking and buoyant-crust tracers)

## Problem

The model has no intra-oceanic arc source:
- Arc magmatism needs a *continental* plate facing a converging oceanic one
  (`lithosphere_plate.continental_arc_band` returns nothing for an oceanic plate).
- Volcanic eruptions and arc cell growth book juvenile continental material on
  continental-coded hosts.
- Where two oceanic plates converge, both sides' contested oceanic-coded floor nodes
  subduct (`boundary_context`: an oceanic plate's `shrinkable` is all of `contested`).
  There is no ocean–ocean overriding plate and no arc.

**Current-code qualification (#403).** Continental sources are not strictly restricted to
pre-existing continental codes. `_erupt_melted_nodes` resets a melted column to continental
crust whenever its *pre-melt elevation* is above 0 m, and `deform_columns` and
`quad_tectonics._open_rift` book the resulting juvenile material. An exposed oceanic-coded
column can therefore become continental through an existing elevation rule. Continental
arc production depends on continental margin, but a continent-free state has not been
shown to be absorbing. The source-attribution gate below must distinguish that existing
route from the proposed oceanic-arc maturation route.

**Historical evidence, pre-#358.** The run behind RFCs 0001–0007 showed an end state.
The run is seed 456862899, 19.2 Gyr, saved locally as
`mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine and
**not available for download**. Numbers come from
`bin/debug/ledger_history.py` on that save. They were first recorded in RFC publication
commit `d795a09` (2026-10-13); the generating commit, density and controls remain unknown.
#403 records 0.1 Myr steps through 21.8 Myr (218 steps), then 10 Myr steps: the old 20-step
cooldown was 2 Myr initially and 200 Myr for most of the run. This save predates #358
(`7afecfd`, merged at `1b4f686`), which changed the cooldown to a physical 2 Myr. The
published values below are preserved without recomputing the private save. It cannot
serve as the current-main baseline or be reproduced from its seed and schedule alone.

- Juvenile additions fell with the inventory, from 2–4 km³/yr while continents were
  Earth-sized to 0.4 km³/yr over the last 7 Gyr.
- By then continental material was down to about 130 M km³, 1.6% of its start, and land
  covered about 2% of the planet.

The historical run is a stress test of long-term behaviour, not a calibration interval.
It motivates checking recovery from a continent-free state. The current code trace above
identifies the missing oceanic-arc mechanism; the old end state neither rules out the
elevation-driven route nor measures disappearance or recovery on current main.

**Current-main baseline, post-#358 (`2562f63`, #403).**
[The fresh evidence](../analysis/issue403/README.md) regenerates default quad worlds at
density 1, with 5 Myr steps, 100 Myr checkpoints, a physical 2 Myr split cooldown and no
control overrides; [tuning.json](../analysis/issue403/tuning.json) records the inputs and
ledger checkpoints. The current first-Gyr results are:

| Seed | Continental inventory (M km³) | Inventory change | Genuine continental area (M km²) | Land fraction | Juvenile additions (km³/yr) |
|---|---|---|---|---|---|
| 7 | 11,955 → 1,677 | −86.0% | 260.22 → 23.46 | 45.32% → 5.52% | 6.85 |
| 11 | 6,917 → 1,866 | −73.0% | 161.57 → 26.91 | 27.70% → 6.95% | 6.05 |
| 2024 | 6,522 → 1,642 | −74.8% | 123.36 → 23.55 | 17.59% → 5.84% | 4.22 |

Land fraction is area-weighted from the climate grid; the raw node land fraction also in
the JSON uses another definition. The historical run's 47% first-Gyr inventory loss is
preserved in RFC 0005, but these different experiments do not measure #358's causal effect.

Juvenile additions at 0–100 Myr are 10.70, 10.36 and 6.76 km³/yr on seeds 7, 11 and
2024; by 500–1,000 Myr they are 4.96, 4.22 and 3.14 as the inventories shrink. These are
aggregate ledger additions, not measured arc production. The refresh confirms rapid
loss on these current baselines, but neither a continent-free state nor its recovery is tested. No oceanic-arc
maturation or docking has been implemented or measured, and the source-attribution gate
in [#406](https://github.com/adubey/mantle-bloom/issues/406) remains open.

The existing `stats.plate_overlap_area_fraction` sums each plate's accounting area on
nodes whose `cover_count > 0`, then divides by sphere area (`stats._overlap_once_totals`).
It includes each overlapping plate's represented area, rather than deduplicating the
overlap footprint, and is distinct from raw node overlap. The matched baseline is:

| Seed | Plate count, 0 → 1 Gyr | Overlap area fraction, 0 → 1 Gyr | Overlap fraction, mean of 11 checkpoints |
|---|---|---|---|
| 7 | 20 → 51 | 0.010900 → 0.189386 | 0.142557 |
| 11 | 9 → 57 | 0.007869 → 0.114516 | 0.147621 |
| 2024 | 11 → 60 | 0.009024 → 0.125223 | 0.142188 |

The final column is a sample mean at 0, 100, …, 1,000 Myr, not a continuous-time integral.
The full checkpoint series in the JSON supports comparisons at matching physical times.

A related scaling issue: the continental arc band is `BOUNDARY_FORCE_REACH_MULTIPLIER` (3)
cells wide, 375 km at density 1 and 188 km at density 4, while
`ARC_MAGMATIC_HC_RATE_M_PER_MYR` is a per-cell thickening rate. Arc productivity *per km of
margin* therefore halves from density 1 to density 4. This RFC sets the new arcs' rate per
km of margin, so it doesn't inherit that.

## Proposed solutions

**What happens in nature.**
- Continental crust is made mainly at subduction zones [[Rudnick, 1995]](#ref-rudnick1995)
  [[Hawkesworth et al., 2010]](#ref-hawkesworth2010)
  [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015), and much of it at *intra-oceanic* arcs:
  ocean–ocean subduction zones such as Izu-Bonin-Mariana, Tonga-Kermadec and the
  Aleutians [[Stern, 2010]](#ref-stern2010ioa).
- When two oceanic plates converge, the older, colder, denser one tends to subduct
  [[Stern, 2002]](#ref-stern2002). Age is a clue, not a rule: where a subduction zone
  starts, and its history afterwards, decide polarity in practice. Once established, a
  subduction zone keeps its polarity until something forces a flip, such as a plateau or
  arc jamming the trench [[Stern, 2002]](#ref-stern2002) [[Kerr & Mahoney, 2007]](#ref-kerr2007).
- Intra-oceanic arc crust grows at about 30–95 km³/Myr per km of arc
  [[Dimalanta et al., 2002]](#ref-dimalanta2002) [[Jicha & Jagoutz, 2015]](#ref-jicha2015).
- It builds a middle crust with continental-like seismic velocities (6.0–6.3 km/s) beneath
  Izu-Bonin [[Suyehiro et al., 1996]](#ref-suyehiro1996). The Mariana arc's crust is about
  20 km thick [[Takahashi et al., 2007]](#ref-takahashi2007), against 6–7 km for normal
  oceanic crust [[Kerr & Mahoney, 2007]](#ref-kerr2007).
- **Arc crust is not continental crust.** Bulk arc crust is more mafic than the
  continents. Turning it into continental crust is thought to need the dense, mafic lower
  arc crust to founder into the mantle, leaving the intermediate-to-felsic part
  [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015). Only part of the arc's gross production
  is preserved as continental crust.
- Oceanic arcs don't stay stable as arcs past about 36 km of crust
  [[Clift & Vannucchi, 2004]](#ref-clift2004).
- Arcs, along with oceanic plateaus and other buoyant features, eventually collide with
  each other and with continents and weld on as terranes
  [[Tetreault & Buiter, 2014]](#ref-tetreault2014). The Central Asian Orogenic Belt grew
  largely this way, from arcs and accretionary complexes
  [[Şengör et al., 1993]](#ref-sengor1993).
- Mass balances say this is *required*: erosive continental arcs are net sinks, so
  continental volume holds steady only if oceanic arcs are accreted to the continents
  [[Clift & Vannucchi, 2004]](#ref-clift2004) [[Clift et al., 2009b]](#ref-clift2009b).
  Accreted oceanic plateaus add about 1.1 km³/yr more [[Clift et al., 2009a]](#ref-clift2009a).

1. **Intra-oceanic arcs that mature in place.**
   - Ocean–ocean convergence gets a persistent polarity, kept in the same front records
     that already hold continental collision polarity.
   - The overriding plate grows an arc band that thickens its crust at a physical rate per
     km of margin. This is gross arc production, tracked as arc crust.
   - **Maturation, wherever the arc is.** Once an arc column is thick enough, the
     intermediate-to-felsic share of its arc crust becomes continental material and the
     mafic remainder founders. This happens in place, on an oceanic plate, with no
     continent needed, so a continent-free world can make its first continental crust.
   - A column that is mostly continental material is recoded continental. It is then a
     continental terrane on an oceanic plate, which the model already carries and docks at
     collision (#253, #305, #321), keeping its keel (#311). Two matured arcs that collide
     join as one terrane, so continents can grow from arcs alone.
   - Immature arc crust, still mafic, can resist subduction and dock through the docking
     rule shared with RFC 0007. It matures later on its new host.
2. **Accreted oceanic plateaus.** Needs mantle plumes, which the model doesn't have (no
   arc or hotspot spawning; `elevation_lines.py`). Covered separately by
   [RFC 0007](0007-mantle-plumes-and-plateaus.md), whose plateaus dock through the same
   terrane path [[Tetreault & Buiter, 2014]](#ref-tetreault2014).
3. **A background juvenile source.** Add juvenile crust at a fixed global rate. Rejected:
   it isn't attached to any process, so it can't respond to the tectonics.

**Recommendation:** 1. It is the process the literature says carries continental growth
when continents are small, and it reuses machinery the model already has: polarity front
records, terranes, docking, the juvenile ledger. The #403 refresh retains that mechanism
recommendation because the source is still absent in current code. It does not establish
a need to tune oceanic-arc production to offset the private pre-#358 run's loss rate;
RFC 0005's unpassed sink diagnosis and this RFC's source-attribution gate remain separate
requirements. Measure the existing elevation-driven route before using a stripped world
as evidence that oceanic arcs are the missing route back.
[#406](https://github.com/adubey/mantle-bloom/issues/406) tracks that diagnosis; the arc
design and knobs remain proposed and unchanged. No oceanic-arc output has been measured
in this refresh.

## Design

**Thickness convention: excess crust.** The model's reference oceanic crust is
`REFERENCE_HC_OCEANIC_M` = 11 km; normal oceanic crust on Earth is 6–7 km
[[Kerr & Mahoney, 2007]](#ref-kerr2007). The literature gives *total* crustal thicknesses,
so this RFC and RFC 0007 work in **excess crust**, E = Hc − 11 km on oceanic-coded cells.
Literature totals are converted by subtracting 6.5 km:
- the Mariana arc's ~20 km is E ≈ 13.5 km, a model Hc of ≈ 24.5 km;
- the 36 km arc stability limit is E ≈ 29.5 km, a model Hc of ≈ 40.5 km.

**Polarity: persistent front records.** `collision_polarity` keeps `CollisionFront`
records at the refreshed baseline for continental contacts only. Each record holds a
pair's frozen lower plate, the evidence that decided it, a lineage through topology changes (`end_topology`,
`note_lineage`), and an expiry (`FRONT_EXPIRY_YEARS`, 5 Myr). Decisions are made in the
step's prepass and frozen for the step, so both plates read the same answer.

This RFC adds ocean–ocean fronts to the same records:
- *Detection:* converging contacts where both sides are oceanic-coded, matched to existing
  records as continental fronts are.
- *Memory:* a separate expiry, `IOA_FRONT_EXPIRY_MYR`, longer than the continental 5 Myr,
  so a subduction zone keeps its polarity through a convergence lull. The consumption
  evidence store (`collision_evidence`) expires after 10 Myr; the record is what persists.
- *Decision, at creation only, in order:*
  1. inherited from a parent or same-pair record (lineage), as continental fronts do;
  2. consumption evidence (`collision_evidence`, #318), if asymmetric beyond
     `IOA_POLARITY_MIN_EVIDENCE_SHARE`;
  3. crust age, if both contact cells have one and they differ by more than
     `IOA_POLARITY_MIN_AGE_GAP_MYR`: the older side subducts;
  4. buoyancy: the side with more excess crust (an arc or plateau) at the contact
     overrides;
  5. the existing fallback tiers (motion, size, plate id), which always decide.
- Every ocean–ocean front gets a polarity, and keeps it until the record expires or a
  topology change splits it; there's no symmetric outcome. Stats report each tier's
  share.

*Crust age.* The model has no seafloor-age field: `divergent_age_myr` counts only
continuous divergence, and oceanic Hm relaxes to a fixed reference (#311). This RFC adds
`crust_formed_years`:
- the `elapsed_years` at which the cell's crust was made;
- HISTORY-remapped (oldest wins), with a sentinel for unknown;
- stamped at spreading-edge insertion (`quad_tectonics`, `ELEV_CHANGE_NEW_CRUST`), rift
  melt-through (`_erupt_melted_nodes`) and gap fill (`gaps.py`).

Old saves load with the sentinel everywhere, so tier 3 is skipped until floor is replaced,
and tiers 4–5 decide meanwhile. No backfill is needed: a uniform assumed age couldn't break
a tie anyway.

**Arc growth (gross production).**
- *The band:* the overriding plate's nodes within `reach_rad` of the subducting
  neighbour, with the continental band's intensity taper (1 at the contact, 0.3 at the
  band edge). The overriding side's contested nodes there stop retreating.
- *Rate per km of margin.* For each front segment with measured contact length L_seg (km,
  the summed edge length of its contact cells), the thickening at each band cell is
  `R × intensity_i`, where
  `R_m_per_Myr = 1000 × P × L_seg / Σ_band (A_i × intensity_i)`
  with P = `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` × convergence factor, and A_i in km². By
  construction the segment gains exactly P × L_seg km³ per Myr, at any density and band
  shape.
- *Convergence factor:* `0.4 + 0.6 × v / 5 cm/yr`, capped at 3×, as for the continental
  arc, since convergence controls long-term arc melt production
  [[Clift & Vannucchi, 2004]](#ref-clift2004).
- Hc and a new `arc_crust_m` grow together, and RFC 0007's age moment
  `buoyant_crust_time_m_yr` grows by the added thickness × the current time; Hm is
  unchanged. Hc is capped at 11 km + `IOA_MAX_EXCESS_M`.
- *Example:* at P = 60 km³/Myr/km on a straight segment at density 4 (band ≈ 188 km, mean
  intensity ≈ 0.65), R ≈ 1000 × 60 / (188 × 0.65) ≈ 490 m/Myr at the contact. A column
  reaches the Mariana-like E ≈ 13.5 km in about 28 Myr at the contact and in 30–90 Myr
  across the band. That's consistent with the Izu-Bonin-Mariana system, which is about
  50 Myr old [[Stern, 2010]](#ref-stern2010ioa).

**Maturation (preservation as continent).** When an oceanic-coded column's E reaches
`IOA_MATURE_EXCESS_M`, with f = `IOA_PRESERVED_FRACTION`:
- f × `arc_crust_m` becomes `continental_material_m`, booked once to
  `juvenile_additions_m3` (sub-account `juvenile_oceanic_arc_m3`);
- the mafic remainder, (1 − f) × `arc_crust_m`, founders: it is removed from Hc, and the
  column subsides isostatically, as `cratons.attrite` does after removing crust. It books
  to the buoyant-crust ledger's `foundered_m3`. This is explicit, not left to
  `orogeny.delaminable_root_m`, which keys on crust below 50 km or on restite and wouldn't
  find it;
- `arc_crust_m` is cleared to 0, so the provenance is used once;
- after that, each step's new arc crust on the column matures the same way at once.

The column is recoded `CRUST_TYPE_CONTINENTAL` once `continental_material_m ≥
IOA_RECODE_SHARE × Hc`. From then on, the whole column takes the continental crust
density, although its lower part is still oceanic basement. That's an approximation: it
makes recoded arcs slightly too buoyant. A two-density column would remove it, but is out
of scope.

*Example:* at f = 0.5, a column at E = 13.5 km gains 6.75 km of continental material and
loses 6.75 km of foundered crust (Hc 24.5 → 17.75 km, share 0.38). With 8.5 km more arc
growth, the share passes 0.5 and the column is recoded.

**Docking of immature arc crust.** Oceanic-coded patches with excess crust can resist
subduction and dock as terranes through `terranes.docking_probability`, shared with RFC
0007, which specifies its form and knobs. It uses E and the patch's volume-weighted crust
age, and it is probabilistic and partial. Docked immature arc crust keeps its `arc_crust_m`
on the new host and matures there under the rule above. Matured, continental-coded arc
cells need no docking rule: the existing terrane path already carries continental-coded
cells on oceanic plates.

**Tracer lifecycle.** `arc_crust_m` is a share of Hc, like the plateau crust of RFC 0007,
and goes through one shared lifecycle (RFC 0007, Design):
- it is held to E by the phase audits;
- it scales down with thinning, and is cleared by melt-through resets;
- it is carried by every transfer, including `crust_transfer.ConsumedColumn`, terrane
  accretion and merge overflow;
- it is booked in a shared buoyant-crust ledger whose reservoirs are disjoint.

**Global scale check.** Gross production: at 60 km³/Myr/km, 10,000 km of active
intra-oceanic arc makes 0.6 km³/yr of arc crust. Preservation is smaller: at f = 0.5, at
most 0.3 km³/yr of continental material, less whatever subducts before maturing. That is a
different quantity from the historical pre-#358 run's 2–4 km³/yr of juvenile additions,
which were booked at production; the ledger did not attribute all of them to arc
magmatism. The acceptance criteria keep production, docking and preservation apart.

**Runtime.**
- The band and thickening reuse this step's boundary search (`inputs.dist_to_neighbor`,
  `neighbor_node_index`); the per-segment sums are O(band nodes).
- Ocean–ocean fronts add records to the prepass, which already matches continental fronts.
  The count grows with ocean–ocean contacts, tens per world; the prepass cost is measured
  in acceptance.
- Maturation is a masked update.
- The two new fields add 16 bytes per node to saves, about 2 MB at density 4 (~130k
  nodes).

**Save compatibility.**
- `crust_formed_years` (sentinel default) and `arc_crust_m` (0 default) are registered in
  `surface_fields.py`, so the registry backfills old saves (see
  `docs/save-compatibility.md`).
- `CollisionFront` gains a `kind` field with default `"continental"`, so pickled records
  load unchanged.
- `test_persistence.py` covers both.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` | 60 | Mid-range of 30–95 for western Pacific arcs [[Dimalanta et al., 2002]](#ref-dimalanta2002) | Stay within 30–95, with sensitivity at both ends; not tuned to close a budget |
| `IOA_MATURE_EXCESS_M` | 13,500 m | Mariana-like crust (~20 km total) [[Takahashi et al., 2007]](#ref-takahashi2007), converted to excess | Sensitivity 9,000–18,000 m |
| `IOA_PRESERVED_FRACTION` | 0.5 | Calibration assumption: the mafic lower arc crust must founder for arc crust to become continental [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015), but no single fraction is established | Sensitivity at 0.3 and 0.7 |
| `IOA_RECODE_SHARE` | 0.5 | Design choice: recode when continental material is the majority of the column | Sensitivity 0.35–0.65 |
| `IOA_MAX_EXCESS_M` | 29,500 m | Oceanic arc stability limit of 36 km total [[Clift & Vannucchi, 2004]](#ref-clift2004), converted to excess | Fixed |
| `IOA_FRONT_EXPIRY_MYR` | 30 Myr | Design choice: a convergence lull shorter than this keeps the polarity | Sweep 10–60 |
| `IOA_POLARITY_MIN_EVIDENCE_SHARE` | 0.7 | Design choice: 70% of a pair's recent consumption evidence on one side | Sweep 0.6–0.9 |
| `IOA_POLARITY_MIN_AGE_GAP_MYR` | 10 Myr | Design choice: below this, the age contrast is a weak clue [[Stern, 2002]](#ref-stern2002) | Sweep 5–30; watch overlap metrics in `plate_diagnostics` |
| Docking rule | — | RFC 0007 | RFC 0007 |

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time. Record the generating commit, density, controls and full step
schedule; a regenerated seed 456862899 is a new run, not the historical save. Report:
- arc crust made (gross), matured, foundered, docked and subducted, kept separate;
- each polarity tier's share, and polarity flips per front;
- arc length by type, and recoded continental area on oceanic plates;
- plate overlap and plate count, since stopping retreat on the overriding side could let
  overlaps grow.

## Detailed design (provisional touchpoints)

**Checked against current main, `2562f63` (2026-10-14).**
`lithosphere_plate.continental_arc_band` still requires a continental owning plate;
`boundary_context` sets oceanic-plate `shrinkable = contested` before applying the existing
continental-front masks. `collision_polarity._detect_fronts` still requires continental
codes on both contact sides, so #395's preserved `run_eligible` mask is not an ocean–ocean
polarity implementation. The 11 km reference, three-spacing arc band, and 450 m/Myr
continental-arc rate are unchanged. `crust_formed_years`, `arc_crust_m`, oceanic-front
state, the oceanic rate helper and `terranes.py` below are all new. The existing
`ConsumedColumn`/transfer path remains the place to carry their volumes. The existing
`_erupt_melted_nodes` elevation rule must be included in the source-attribution diagnosis;
its result is continental-coded before the juvenile booking, even when its input was
oceanic-coded.

- `backend/app/surface_fields.py`: register `crust_formed_years` (HISTORY, sentinel) and
  `arc_crust_m` (EXTENSIVE), with comments citing this RFC. The shared
  `buoyant_crust_time_m_yr` (EXTENSIVE) is registered by whichever of this RFC and RFC
  0007 lands first.
- `backend/app/collision_polarity.py`
  - `CollisionFront.kind` ("continental" | "oceanic").
  - Ocean–ocean front detection in the prepass.
  - Decision tiers 1–5 for oceanic fronts, and `IOA_FRONT_EXPIRY_MYR` for their expiry.
  - Lineage through `end_topology`, as for continental fronts.
- `backend/app/quad_tectonics.py`
  - Stamp `crust_formed_years` on inserted ridge cells and on melted-through rift cells.
  - `_retreat`: oceanic patches that `terranes.docking_probability` selects join the
    `terrane` donors.
- `backend/app/gaps.py`: stamp `crust_formed_years` on gap-filled crust.
- `backend/app/lithosphere_plate.py`
  - New `oceanic_arc_band(plate, inputs, closing_rate, frame)`, next to
    `continental_arc_band`, reading the frozen oceanic front records.
  - `boundary_context`: the overriding plate's contested nodes leave `shrinkable`; fill
    `arc_band` / `arc_intensity` for oceanic plates.
  - `deform_columns`: per-segment arc thickening of Hc and `arc_crust_m`, then maturation
    and recoding.
- `backend/app/rheology.py`: new `oceanic_arc_rate_m_per_myr(productivity,
  contact_length_km, band_areas_km2, intensity)` and new constants, with citations.
- `backend/app/terranes.py` (new, shared with RFC 0007): docking, the Hc-share tracer
  lifecycle, and the buoyant-crust ledger.
- `backend/app/continental_ledger.py`: a diagnostic sub-account `juvenile_oceanic_arc_m3`
  of `juvenile_additions_m3`, outside the balance.
- `backend/app/stats.py`: intra-oceanic arc length; the arc-crust fates; polarity tier
  shares; recoded area on oceanic plates.
- `backend/app/render_image.py`: `elevReason` already has `ELEV_CHANGE_SUBDUCTION_ARC`.
  Stamp it on oceanic arc uplift too, so arcs show without a new view.
- Tests:
  - `test_collision_polarity.py`:
    - an ocean–ocean front keeps its polarity through a 15 Myr lull with no convergence;
    - it carries its polarity through a split and a merge (lineage);
    - both plates read the same frozen decision within a step;
    - each tier decides when the ones above it are silent, and the fallback always
      decides.
  - `test_rheology.py`: on sampled band geometries at densities 1 and 4 (straight and
    curved segments), the integrated volume added per measured km of contact is within 5%
    of P. This checks the volume integral, not just the rate helper.
  - `test_lithosphere_plate.py` (new, mirroring the module):
    - maturation books exactly f × `arc_crust_m` as juvenile once, removes
      (1 − f) × `arc_crust_m` from Hc with the isostatic subsidence, and clears the tracer;
    - recoding happens at the share threshold;
    - nothing continental is booked at production.
  - `test_persistence.py`: an old save loads with the sentinel, zero arc crust, and
    continental-kind front records.
- Docs: `docs/simulation-model.md` (Plate motion boundary classification, Volcanism, and
  Known simplifications for the column-density approximation), and
  `docs/save-compatibility.md` for the new fields.

## Assumptions and uncertainty

- **Polarity is simplified.** The tiers stand in for subduction initiation, which the model
  doesn't simulate: no initiation at transforms or passive margins. Flips happen only
  through record expiry or jamming (RFC 0007).
- **Maturation is a parameterisation.** One threshold and one preserved fraction stand in
  for differentiation and lower-crust foundering. The literature agrees much of the mafic
  lower arc crust is lost [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015), not on a number.
  The sensitivity runs bound the effect.
- **Recoded columns use one density.** A recoded arc column takes continental density for
  its oceanic basement too, slightly overstating its buoyancy.
- **Productivity is modern and regional.** 30–95 km³/Myr/km is from western Pacific arcs
  [[Dimalanta et al., 2002]](#ref-dimalanta2002).
- **The thickness conversion is approximate.** Subtracting 6.5 km maps literature totals
  onto the model's 11 km reference, but the model's crust isn't calibrated layer by layer.
- **Docking is shared with RFC 0007** and inherits its uncertainty.
- **Overriding plates no longer retreat.** That changes ocean–ocean consumption, which
  could let overlaps grow. Acceptance watches it.

## Acceptance criteria

- [ ] **Source-attribution gate ([#406](https://github.com/adubey/mantle-bloom/issues/406),
      raised by #403):** run the stripped world below on the matched,
      pinned current-main baseline first. Record continental recoding and juvenile volume
      from `_erupt_melted_nodes` (both column and `_open_rift` callers), with prior effective
      crust type and elevation. Establish which existing routes can create continent before
      attributing recovery to intra-oceanic arcs. The aggregate juvenile account cannot
      pass this gate by itself.
- [ ] On seeds 7, 11 and 2024 to 1 Gyr, gross arc production, matured, foundered, docked
      and subducted arc crust are each reported, and the buoyant-crust ledger closes.
      Matured continental material from intra-oceanic arcs is 0.05–0.5 km³/yr at the
      starting knobs, and at least one matured arc terrane docks onto another terrane or a
      continent per run.
- [ ] **Stripped world.** A test world with every plate's `crust_type` set to oceanic,
      every cell oceanic-coded, and all continental material and craton removed at t = 0.
      With no continental plates, the continental-arc machinery can't fire. Over 500 Myr,
      arc-derived recoded continental area appears and grows. Report
      `juvenile_oceanic_arc_m3` separately from the existing elevation-driven reset
      source measured by the gate; do not require every juvenile addition to be from
      oceanic arcs.
- [ ] Combined with RFC 0005, the continental inventory on seed 456862899 changes by less
      than 15% per Gyr over the first 3 Gyr, and never falls below 50% of its initial value
      over 5 Gyr.
- [ ] Ocean–ocean polarity flips less than once per front per 50 Myr, and the fallback
      tier decides less than 30% of fronts.
- [ ] `stats.plate_overlap_area_fraction` and plate count change by less than 15%
      against the matched, pinned current-main run at the same physical checkpoints
      (same density, controls and step schedule). Use the current area-fraction and
      count baselines above, not raw node overlap or the historical pre-#358 save.
- [ ] Integrated arc volume per km of contact is within 5% of P at densities 1 and 4 (unit
      test).
- [ ] Step time changes by less than 3%; save size by less than 3%.
- [ ] Ledgers close; old saves load; affected tests, ruff and `bin/quality_check.sh` are
      green.

## Falsification criteria

- **Arcs don't mature or survive.** If under 10% of gross arc production ever matures or
  docks over 1 Gyr, at every sensitivity case, arcs in the model are consumed too fast.
  Then this mechanism can't supply continental crust whatever its production, and the
  RFC should be revised or dropped.
- **Polarity doesn't settle.** If ocean–ocean polarity flips more than once per front per
  50 Myr even with the front records, the front matching or lineage is wrong, and the
  polarity design needs rework.
- **The stripped world doesn't recover.** If, with intra-oceanic arcs, the stripped test
  world gains no arc-derived continental area in 500 Myr, the proposed maturation route
  fails for this model. If the baseline already recovers through elevation-driven resets,
  revise the claim that arcs are required for recovery; assess their geological role
  separately from that premise.

## References

- <a id="ref-clift2004"></a>Clift, P., & Vannucchi, P. (2004). Controls on tectonic
  accretion versus erosion in subduction zones: Implications for the origin and recycling
  of the continental crust. *Reviews of Geophysics*, 42, RG2001.
  [doi:10.1029/2003RG000127](https://doi.org/10.1029/2003RG000127)
- <a id="ref-clift2009a"></a>Clift, P. D., Vannucchi, P., & Morgan, J. P. (2009a). Crustal
  redistribution, crust–mantle recycling and Phanerozoic evolution of the continental
  crust. *Earth-Science Reviews*, 97, 80–104.
  [doi:10.1016/j.earscirev.2009.10.003](https://doi.org/10.1016/j.earscirev.2009.10.003)
- <a id="ref-clift2009b"></a>Clift, P. D., Schouten, H., & Vannucchi, P. (2009b).
  Arc-continent collisions, sediment recycling and the maintenance of the continental
  crust. *Geological Society, London, Special Publications*, 318, 75–103.
  [doi:10.1144/SP318.3](https://doi.org/10.1144/SP318.3)
- <a id="ref-dimalanta2002"></a>Dimalanta, C., Taira, A., Yumul, G. P., Tokuyama, H., &
  Mochizuki, K. (2002). New rates of western Pacific island arc magmatism from seismic and
  gravity data. *Earth and Planetary Science Letters*, 202, 105–115.
  [doi:10.1016/S0012-821X(02)00761-6](https://doi.org/10.1016/S0012-821X(02)00761-6)
- <a id="ref-hawkesworth2010"></a>Hawkesworth, C. J., Dhuime, B., Pietranik, A. B.,
  Cawood, P. A., Kemp, A. I. S., & Storey, C. D. (2010). The generation and evolution of
  the continental crust. *Journal of the Geological Society*, 167, 229–248.
  [doi:10.1144/0016-76492009-072](https://doi.org/10.1144/0016-76492009-072)
- <a id="ref-jagoutz2015"></a>Jagoutz, O., & Kelemen, P. B. (2015). Role of arc processes
  in the formation of continental crust. *Annual Review of Earth and Planetary Sciences*,
  43, 363–404.
  [doi:10.1146/annurev-earth-040809-152345](https://doi.org/10.1146/annurev-earth-040809-152345)
- <a id="ref-jicha2015"></a>Jicha, B. R., & Jagoutz, O. (2015). Magma production rates for
  intraoceanic arcs. *Elements*, 11, 105–111.
  [doi:10.2113/gselements.11.2.105](https://doi.org/10.2113/gselements.11.2.105)
- <a id="ref-kerr2007"></a>Kerr, A. C., & Mahoney, J. J. (2007). Oceanic plateaus:
  Problematic plumes, potential paradigms. *Chemical Geology*, 241, 332–353.
  [doi:10.1016/j.chemgeo.2007.01.019](https://doi.org/10.1016/j.chemgeo.2007.01.019)
- <a id="ref-rudnick1995"></a>Rudnick, R. L. (1995). Making continental crust. *Nature*,
  378, 571–578. [doi:10.1038/378571a0](https://doi.org/10.1038/378571a0)
- <a id="ref-sengor1993"></a>Şengör, A. M. C., Natal'in, B. A., & Burtman, V. S. (1993).
  Evolution of the Altaid tectonic collage and Palaeozoic crustal growth in Eurasia.
  *Nature*, 364, 299–307. [doi:10.1038/364299a0](https://doi.org/10.1038/364299a0)
- <a id="ref-stern2002"></a>Stern, R. J. (2002). Subduction zones. *Reviews of
  Geophysics*, 40(4), 1012. [doi:10.1029/2001RG000108](https://doi.org/10.1029/2001RG000108)
- <a id="ref-stern2010ioa"></a>Stern, R. J. (2010). The anatomy and ontogeny of modern
  intra-oceanic arc systems. *Geological Society, London, Special Publications*, 338,
  7–34. [doi:10.1144/SP338.2](https://doi.org/10.1144/SP338.2)
- <a id="ref-suyehiro1996"></a>Suyehiro, K., Takahashi, N., Ariie, Y., Yokoi, Y., Hino,
  R., Shinohara, M., Kanazawa, T., Hirata, N., Tokuyama, H., & Taira, A. (1996).
  Continental crust, crustal underplating, and low-Q upper mantle beneath an oceanic
  island arc. *Science*, 272, 390–392.
  [doi:10.1126/science.272.5260.390](https://doi.org/10.1126/science.272.5260.390)
- <a id="ref-takahashi2007"></a>Takahashi, N., Kodaira, S., Klemperer, S. L., Tatsumi, Y.,
  Kaneda, Y., & Suyehiro, K. (2007). Crustal structure and evolution of the Mariana
  intra-oceanic island arc. *Geology*, 35, 203–206.
  [doi:10.1130/G23212A.1](https://doi.org/10.1130/G23212A.1)
- <a id="ref-tetreault2014"></a>Tetreault, J. L., & Buiter, S. J. H. (2014). Future
  accreted terranes: A compilation of island arcs, oceanic plateaus, submarine ridges,
  seamounts, and continental fragments. *Solid Earth*, 5, 1243–1275.
  [doi:10.5194/se-5-1243-2014](https://doi.org/10.5194/se-5-1243-2014)
