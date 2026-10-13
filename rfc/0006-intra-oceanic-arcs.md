# RFC 0006: Intra-oceanic arcs: new continent from the ocean

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #120 (arc magmatism vs. land decline), #311 (terrane keels), #318
  (consumption history), #321 (terrane docking); RFCs
  [0004](0004-proto-cratons.md), [0005](0005-continental-recycling-balance.md)

## Problem

In the model, new continental crust can only form where continental crust already is:
- Arc magmatism needs a *continental* plate facing a converging oceanic one
  (`lithosphere_plate.continental_arc_band` returns nothing for an oceanic plate).
- Rift melts, volcanic eruptions and arc cell growth book juvenile continental material
  only on continental hosts.
- Where two oceanic plates converge, both sides' contested nodes simply subduct
  (`boundary_context`: an oceanic plate's `shrinkable` is all of `contested`). There is no
  overriding plate and no arc.

Production therefore scales with the length of continental margin. As continents shrink,
production shrinks with them, and a world that loses its continents can't make new ones.

The run behind RFCs 0001–0006 shows the end state. The run is seed 456862899, 19.2 Gyr,
saved locally as `mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the
author's machine and **not available for download**. Numbers come from
`bin/debug/ledger_history.py` on that save.

- Juvenile additions fell with the inventory, from 2–4 km³/yr while continents were
  Earth-sized to 0.4 km³/yr over the last 7 Gyr.
- By then continental material was down to about 130 M km³, 1.6% of its start, and land
  covered about 2% of the planet.

The run is a stress test of the model's long-term behaviour, not a calibration interval:
it shows that the model has no route back from a continent-free state. It does not say how
fast continents should regrow.

A related scaling issue: the continental arc band is `BOUNDARY_FORCE_REACH_MULTIPLIER` (3)
cells wide, 375 km at density 1 and 188 km at density 4, while `ARC_MAGMATIC_HC_RATE_M_PER_MYR`
is a per-cell thickening rate. Arc productivity *per km of margin* therefore halves from
density 1 to density 4. This RFC sets the new arcs' rate per km of margin, so it doesn't
inherit that.

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
  Izu-Bonin [[Suyehiro et al., 1996]](#ref-suyehiro1996). The Mariana arc has grown about
  20 km of crust [[Takahashi et al., 2007]](#ref-takahashi2007).
- **Arc crust is not continental crust.** Bulk arc crust is more mafic than the
  continents. Turning it into continental crust is thought to need the dense, mafic lower
  arc crust to founder into the mantle, leaving the intermediate-to-felsic part
  [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015). Only part of the arc's gross production
  is preserved as continental crust.
- Oceanic arcs don't stay stable as arcs past about 36 km of crust
  [[Clift & Vannucchi, 2004]](#ref-clift2004).
- Arcs, along with oceanic plateaus and other buoyant features, eventually collide with
  and weld onto continents as terranes [[Tetreault & Buiter, 2014]](#ref-tetreault2014).
  The Central Asian Orogenic Belt grew largely this way
  [[Şengör et al., 1993]](#ref-sengor1993).
- Mass balances say this is *required*: erosive continental arcs are net sinks, so
  continental volume holds steady only if oceanic arcs are accreted to the continents
  [[Clift & Vannucchi, 2004]](#ref-clift2004) [[Clift et al., 2009b]](#ref-clift2009b).
  Accreted oceanic plateaus add about 1.1 km³/yr more [[Clift et al., 2009a]](#ref-clift2009a).

1. **Intra-oceanic arcs.**
   - Give ocean–ocean convergence a persistent polarity. It is set when subduction
     starts, from boundary history where there is any and from crust age where there
     isn't, and it is kept afterwards.
   - The overriding plate grows an arc band that thickens its crust at a physical rate per
     km of margin. This is gross arc production, tracked as arc crust on the oceanic
     plate, not as continental material.
   - Thick arc crust resists subduction and docks as a terrane, through the same docking
     rule as RFC 0007's plateaus.
   - When it docks onto a continent, an explicit preserved fraction of the arc crust
     becomes continental material. The rest is the mafic root that, in nature, founders.
     That docked part then joins the continent like other terranes (#253, #305, #321) and
     keeps its keel (#311).
2. **Accreted oceanic plateaus.** Needs mantle plumes, which the model doesn't have (no
   arc or hotspot spawning; `elevation_lines.py`). Covered separately by
   [RFC 0007](0007-mantle-plumes-and-plateaus.md), whose plateaus dock through the same
   terrane path [[Tetreault & Buiter, 2014]](#ref-tetreault2014).
3. **A background juvenile source.** Add juvenile crust at a fixed global rate. Rejected:
   it isn't attached to any process, so it can't respond to the tectonics.

**Recommendation:** 1. It is the process the literature says carries continental growth
when continents are small, it uses machinery the model already has (terranes, docking,
the juvenile ledger), and it gives RFC 0004's collisional cratons more collisions to work
with.

## Design

**Polarity.** Decided per ocean–ocean boundary segment (a plate pair) when convergence
starts, then kept.
1. **History first.** `world.collision_evidence` (#318) already records, per plate pair,
   which plate's floor went down (`ROLE_LOWER`) and which overrode it (`ROLE_UPPER`), with
   a lookback window. If the pair has consumption evidence that is asymmetric beyond
   `IOA_POLARITY_MIN_EVIDENCE_SHARE`, that polarity holds. Once arcs exist, this makes
   polarity persistent: the established upper plate keeps overriding.
2. **Age at initiation.** With no or symmetric evidence, compare crust ages at the
   contact. The model has no seafloor-age field (`divergent_age_myr` counts only
   continuous divergence, and oceanic Hm relaxes to a fixed reference, #311), so this RFC
   adds `crust_formed_years`:
   - the `elapsed_years` at which the cell's crust was made, HISTORY-remapped (oldest
     wins), with a sentinel for unknown;
   - stamped at spreading-edge insertion (`quad_tectonics`, `ELEV_CHANGE_NEW_CRUST`),
     rift melt-through (`_erupt_melted_nodes`) and gap fill (`gaps.py`);
   - the older side by more than `IOA_POLARITY_MIN_AGE_GAP_MYR` subducts.
3. **Otherwise symmetric.** If ages are unknown or too close, today's symmetric
   consumption is kept and no arc forms.

**How long the symmetric fallback lasts on old saves.** Cells keep the unknown sentinel
until they are replaced by new crust. Most ocean floor is recycled within ~200 Myr, but
old floor far from trenches can survive much longer, and the sentinel doesn't spread on
its own. Two things bound it:
- step 1 needs no ages, so a boundary that has consumed asymmetrically keeps its
  polarity;
- the share of ocean–ocean contact length with an unknown age on either side is
  reported in stats, and an acceptance criterion caps it 300 Myr after load.

If that cap fails, old saves backfill unknown oceanic ages once, on first load, with
`elapsed_years − OLD_FLOOR_ASSUMED_AGE_MYR`. That makes every old cell the same age, so
polarity there comes from history alone.

**Arc growth (gross production).**
- The arc band is the overriding plate's nodes within `reach_rad` of the converging
  subducting neighbour, with the same intensity taper as the continental band (1 at the
  contact, 0.3 at the band edge, mean ≈ 0.65). The overriding side's contested nodes there
  stop retreating.
- Productivity is set per km of margin: `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` = P. It is
  converted per step to a thickening at the contact, `R = P / (W × 0.65)` in m/Myr, where
  `W` is the band width in km at the world's density, so P holds at any density.
- P scales with convergence the way the continental arc does (`0.4 + 0.6 × v / 5 cm/yr`,
  capped at 3×), since convergence controls long-term arc melt production
  [[Clift & Vannucchi, 2004]](#ref-clift2004).
- Hc and a new `arc_crust_m` (EXTENSIVE, m) grow together; Hc is capped at
  `IOA_MAX_HC_M`, and Hm is unchanged. The cell stays oceanic-coded: arc crust is not
  continental crust yet.

At P = 60 km³/Myr/km and density 4 (W ≈ 188 km), R ≈ 490 m/Myr at the contact, and a column
goes from 11 km to 20 km in about 20 Myr there and in 30–60 Myr across the band. That is
consistent with the Mariana and Izu-Bonin arcs, which are about 50 Myr old
[[Stern, 2010]](#ref-stern2010ioa) with about 20 km of crust
[[Takahashi et al., 2007]](#ref-takahashi2007).

**Docking.** Oceanic cells whose arc (or plateau) crust makes them buoyant resist
subduction and dock as terranes. One shared docking rule serves this RFC and RFC 0007
(`terranes.docking_probability`, a new module). Following
[[Cloos, 1993]](#ref-cloos1993) as summarised by
[[Kerr & Mahoney, 2007]](#ref-kerr2007), it is probabilistic, not a threshold:
- docking probability rises with Hc and falls with the time since the crust was built;
- a cell that doesn't dock subducts, as thick plateaus and arcs sometimes do;
- docking cells go through the existing terrane path in `quad_tectonics._retreat`.

The rule's form and knobs are specified in RFC 0007.

**Preservation as continent.** When docked arc crust lands on a continental upper plate:
- `IOA_PRESERVED_FRACTION` × `arc_crust_m` becomes `continental_material_m`, booked to
  `juvenile_additions_m3` (sub-account `juvenile_oceanic_arc_m3`), and the cell is
  recoded continental;
- the remaining `(1 − f) × arc_crust_m` is the mafic root that, in nature, founders
  [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015). It stays in Hc as non-continental crust
  (continental material < Hc is allowed), where the existing root-delamination path
  (`orogeny`) can shed it like any other dense root.

Gross production (arc thickening), docking, and preserved continental material are
reported separately.

**Global scale check.** Gross production: at 60 km³/Myr/km, 10,000 km of active
intra-oceanic arc makes 0.6 km³/yr of arc crust. Net preservation is smaller: only arc
crust that docks becomes continental, and only `IOA_PRESERVED_FRACTION` of it, so at a
docked share d and f = 0.5 the continental gain is 0.3·d km³/yr. Those are different
quantities from the 2–4 km³/yr of continental-arc juvenile additions the ledger books
today, which are booked at production, not at preservation. Comparisons in the acceptance
criteria keep the three apart.

**Runtime.**
- The band and thickening reuse this step's boundary search (`inputs.dist_to_neighbor`,
  `neighbor_node_index`); the polarity test is one field read at the neighbour's node.
- Both are O(boundary nodes) per oceanic plate, under 1 ms each.
- The evidence lookup per boundary segment reuses the store #318 maintains.
- The two new fields add 16 bytes per node to saves, about 2 MB at density 4 (~130k
  nodes).

**Save compatibility.** `crust_formed_years` (sentinel default) and `arc_crust_m` (0
default) are registered in `surface_fields.py`, so the registry backfills old saves (see
`docs/save-compatibility.md`); `test_persistence.py` covers it. Unknown ages fall back as
described under Polarity.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` | 60 | Mid-range of 30–95 for western Pacific arcs [[Dimalanta et al., 2002]](#ref-dimalanta2002) | Stay within 30–95, with sensitivity at both ends; not tuned to close a budget |
| `IOA_PRESERVED_FRACTION` | 0.5 | Calibration assumption: the mafic lower arc crust must founder for arc crust to become continental [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015), but no single fraction is established | Sensitivity at 0.3 and 0.7 |
| `IOA_MAX_HC_M` | 36,000 m | Oceanic arc stability limit [[Clift & Vannucchi, 2004]](#ref-clift2004) | Fixed |
| `IOA_POLARITY_MIN_EVIDENCE_SHARE` | 0.7 | Design choice: 70% of a pair's recent consumption evidence on one side | Sweep 0.6–0.9 |
| `IOA_POLARITY_MIN_AGE_GAP_MYR` | 10 Myr | Design choice: below this, the age contrast is a weak clue [[Stern, 2002]](#ref-stern2002) | Sweep 5–30; watch overlap metrics in `plate_diagnostics` |
| `OLD_FLOOR_ASSUMED_AGE_MYR` (only if needed) | 60 Myr | Design choice for the one-off backfill | — |
| Docking rule | — | RFC 0007 | RFC 0007 |

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time. Report:
- arc crust made (gross), arc crust docked, and arc crust preserved as continental, kept
  separate;
- how each ocean–ocean polarity was decided (history, age, or fallback);
- arc length by type;
- terrane area on oceanic plates, and terrane area docked per Gyr;
- plate overlap and plate count, since stopping retreat on the overriding side could let
  overlaps grow.

## Detailed design (provisional touchpoints)

- `backend/app/surface_fields.py`: register `crust_formed_years` (HISTORY, sentinel) and
  `arc_crust_m` (EXTENSIVE), with comments citing this RFC.
- `backend/app/quad_tectonics.py`
  - Stamp `crust_formed_years` on inserted ridge cells and on melted-through rift cells.
  - `_retreat`: oceanic cells that `terranes.docking_probability` selects join the
    `terrane` donors.
  - Where docked crust lands on a continental upper plate (`_accrete_onto_survivors` /
    `crust_transfer`): convert `IOA_PRESERVED_FRACTION × arc_crust_m` to continental
    material and book `juvenile_oceanic_arc_m3`.
- `backend/app/gaps.py`: stamp `crust_formed_years` on gap-filled crust.
- `backend/app/collision_polarity.py`: a helper `ocean_pair_polarity(world, a, b)`,
  returning the established upper plate from consumption evidence, or None.
- `backend/app/lithosphere_plate.py`
  - New `oceanic_arc_band(plate, inputs, closing_rate, polarity)`, next to
    `continental_arc_band`.
  - `boundary_context`: decide polarity per pair (history, then age, then symmetric); the
    overriding plate's contested nodes leave `shrinkable`; fill `arc_band` /
    `arc_intensity` for oceanic plates.
  - `deform_columns`: apply intra-oceanic thickening to Hc and `arc_crust_m`. No
    continental booking at production.
- `backend/app/terranes.py` (new, shared with RFC 0007): `docking_probability(hc, age)`
  and the deterministic draw.
- `backend/app/rheology.py`: `apply_oceanic_arc_thickening(hc, closing_rate, years_myr,
  intensity, band_width_km)` with the per-km conversion; constants with citations.
- `backend/app/continental_ledger.py`: a diagnostic sub-account `juvenile_oceanic_arc_m3`
  of `juvenile_additions_m3`, outside the balance.
- `backend/app/stats.py`: intra-oceanic arc length; arc crust made, docked and preserved;
  how polarity was decided; the share of ocean–ocean contact with unknown age.
- `backend/app/render_image.py`: `elevReason` already has `ELEV_CHANGE_SUBDUCTION_ARC`.
  Stamp it on oceanic arc uplift too, so arcs show without a new view.
- Tests:
  - `test_lithosphere_plate.py` (new, mirroring the module):
    - established consumption evidence decides polarity even against the age contrast;
    - with no evidence, the older side subducts;
    - unknown ages with no evidence keep symmetric retreat.
  - `test_rheology.py`: per-km productivity is density-independent; the cap holds at
    36 km.
  - `test_quad_tectonics.py`: docked arc crust on a continental upper plate books exactly
    `IOA_PRESERVED_FRACTION × arc_crust_m` as juvenile, and nothing is booked at
    production.
  - `test_persistence.py`: an old save loads with the sentinel and zero arc crust.
- Docs: `docs/simulation-model.md` (Plate motion boundary classification, Volcanism, and
  Known simplifications), and `docs/save-compatibility.md` for the new fields.

## Assumptions and uncertainty

- **Polarity is simplified.** History and age stand in for subduction initiation, which
  the model doesn't simulate: no initiation at transforms or passive margins. Flips
  happen only through jamming (RFC 0007).
- **The preserved fraction is a calibration assumption.** The literature agrees much of
  the mafic lower arc crust is lost [[Jagoutz & Kelemen, 2015]](#ref-jagoutz2015), not on
  a number. The sensitivity runs bound its effect.
- **Productivity is modern and regional.** 30–95 km³/Myr/km is from western Pacific
  arcs [[Dimalanta et al., 2002]](#ref-dimalanta2002).
- **Docking is shared with RFC 0007** and inherits its uncertainty.
- **Overriding plates no longer retreat.** That changes ocean–ocean consumption, which
  could let overlaps grow. Acceptance watches it.

## Acceptance criteria

- [ ] On seeds 7, 11 and 2024 to 1 Gyr, gross arc production, docked arc crust and
      preserved continental material are each reported. Preserved continental material
      from intra-oceanic arcs is 0.1–1.0 km³/yr at the starting knobs, and at least one
      arc terrane docks onto a continent per run.
- [ ] 300 Myr after loading an old save, less than 10% of ocean–ocean contact length
      depends on the symmetric fallback. Otherwise the one-off age backfill is enabled.
- [ ] Combined with RFC 0005, the continental inventory on seed 456862899 changes by less
      than 15% per Gyr over the first 3 Gyr, and never falls below 50% of its initial value
      over 5 Gyr.
- [ ] In a world stripped of continents (a test world with all continental crust removed
      at t = 0), preserved continental material stays above zero and genuine continental
      area grows over 500 Myr.
- [ ] Plate overlap area fraction and plate count change by less than 15% against
      `origin/main`.
- [ ] Intra-oceanic arc productivity per km of margin is within 10% between density 1 and
      density 4 (unit test).
- [ ] Step time changes by less than 3%; save size by less than 3%.
- [ ] Ledgers close; old saves load; affected tests, ruff and `bin/quality_check.sh` are
      green.

## Falsification criteria

- **Arcs don't dock.** If almost no arc crust docks (under 10% of gross production over
  1 Gyr, at every docking sensitivity case), arcs in the model are consumed before they
  meet a continent. Then this mechanism can't supply continental crust, whatever its
  production, and the RFC should be revised or dropped.
- **Polarity doesn't settle.** If ocean–ocean polarity flips often (more than once per
  boundary per 50 Myr) once history decides it, history isn't a stable proxy. The
  polarity design is wrong.
- **The stripped world doesn't recover.** If, with intra-oceanic arcs, the stripped test
  world still gains no continental area in 500 Myr, the premise that arcs are the missing
  route back is wrong for this model. Look to RFC 0007 or to continent seeding instead.

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
- <a id="ref-cloos1993"></a>Cloos, M. (1993). Lithospheric buoyancy and collisional
  orogenesis: Subduction of oceanic plateaus, continental margins, island arcs, spreading
  ridges, and seamounts. *Geological Society of America Bulletin*, 105, 715–737.
  [doi:10.1130/0016-7606(1993)105<0715:LBACOS>2.3.CO;2](https://doi.org/10.1130/0016-7606(1993)105%3C0715:LBACOS%3E2.3.CO;2)
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
