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
- When two oceanic plates converge, the older, colder, denser one subducts
  [[Stern, 2002]](#ref-stern2002).
- Intra-oceanic arc crust grows at about 30–95 km³/Myr per km of arc
  [[Dimalanta et al., 2002]](#ref-dimalanta2002) [[Jicha & Jagoutz, 2015]](#ref-jicha2015).
- It builds a middle crust with continental-like seismic velocities (6.0–6.3 km/s) beneath
  Izu-Bonin [[Suyehiro et al., 1996]](#ref-suyehiro1996). The Mariana arc has grown about
  20 km of crust [[Takahashi et al., 2007]](#ref-takahashi2007).
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
   - Give ocean–ocean convergence a polarity: the older floor subducts, the younger side
     overrides.
   - The overriding plate grows an arc band that thickens its crust at a physical rate per
     km of margin.
   - Arc crust past `IOA_CONTINENTAL_HC_M` is recoded continental, and its arc-built part
     is booked as juvenile continental material.
   - It is then a continental terrane on an oceanic plate. The model already carries such
     terranes and docks them at collision (#253, #305, #321), and it keeps their keel
     (#311).
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

**Polarity.** The model has no seafloor-age field: `divergent_age_myr` counts only
continuous divergence, and oceanic Hm relaxes to a fixed reference (#311), so neither
stands in for age.
- Add a per-node `crust_formed_years`, the `elapsed_years` at which the cell's crust was
  made. It is HISTORY-remapped (oldest wins), with a sentinel for unknown.
- It is stamped where oceanic crust is made: new cells at spreading edges
  (`quad_tectonics` cell insertion, `ELEV_CHANGE_NEW_CRUST`), rift melt-through resets
  (`_erupt_melted_nodes`), and gap fill (`gaps.py`).
- At an ocean–ocean contact, the side whose contact cell is older by more than
  `IOA_POLARITY_MIN_AGE_GAP_MYR` subducts. The other side overrides: its contested nodes
  stop retreating and form the arc band.
- If either age is unknown, or the gap is smaller, today's symmetric consumption is kept
  and no arc forms. Old saves are fully dated within one ocean-floor lifetime, about 200
  Myr, as new floor replaces old.

**Arc growth.**
- The arc band is the overriding plate's nodes within `reach_rad` of the converging older
  neighbour, with the same intensity taper as the continental band (1 at the contact, 0.3
  at the band edge, mean ≈ 0.65).
- Productivity is set per km of margin: `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` = P. It is
  converted per step to a thickening at the contact, `R = P / (W × 0.65)` in m/Myr, where
  `W` is the band width in km at the world's density, so P holds at any density.
- P scales with convergence the way the continental arc does (`0.4 + 0.6 × v / 5 cm/yr`,
  capped at 3×), since convergence controls long-term arc melt production
  [[Clift & Vannucchi, 2004]](#ref-clift2004).
- Hc is capped at `IOA_MAX_HC_M`; Hm is unchanged.

**Becoming continent.**
- When Hc ≥ `IOA_CONTINENTAL_HC_M`, the cell's `crust_type_code` becomes
  `CRUST_TYPE_CONTINENTAL`.
- `continental_material_m` gains `max(Hc − REFERENCE_HC_OCEANIC_M, 0)`, the arc-built part,
  booked to `juvenile_additions_m3`. The oceanic basement under it isn't counted as
  continental.
- After conversion, further arc growth books as juvenile in the normal way.

At P = 60 km³/Myr/km and density 4 (W ≈ 188 km), R ≈ 490 m/Myr at the contact, and a
column goes from 11 km to 20 km in about 20 Myr there and in 30–60 Myr across the band. That
is consistent with the Mariana and Izu-Bonin arcs, which are about 50 Myr old
[[Stern, 2010]](#ref-stern2010ioa) with about 20 km of crust
[[Takahashi et al., 2007]](#ref-takahashi2007).

**Global scale check.** At 60 km³/Myr/km, 10,000 km of active intra-oceanic arc makes
0.6 km³/yr of crust. Only the arc-built part of converted cells counts as continental, so
the expected juvenile contribution is ~0.3–0.6 km³/yr on an Earth-like world. That is
meaningful next to the 2–4 km³/yr of continental arcs, and it does not run away.

**Runtime.**
- The band and thickening reuse this step's boundary search (`inputs.dist_to_neighbor`,
  `neighbor_node_index`); the polarity test is one field read at the neighbour's node.
- Both are O(boundary nodes) per oceanic plate, under 1 ms each.
- The new field adds 8 bytes per node to saves, about 1 MB at density 4 (~130k nodes).

**Save compatibility.** `crust_formed_years` is registered in `surface_fields.py` with its
sentinel as the default, so the registry backfills old saves (see
`docs/save-compatibility.md`); `test_persistence.py` covers it. Unknown ages fall back to
today's behaviour, as above.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `IOA_PRODUCTIVITY_KM3_PER_MYR_PER_KM` | 60 | Mid-range of 30–95 [[Dimalanta et al., 2002]](#ref-dimalanta2002) | Sweep 30–95. Choose so total juvenile additions match RFC 0005's sinks within ±15% when continents are Earth-sized |
| `IOA_CONTINENTAL_HC_M` | 20,000 m | Mariana/Izu-Bonin crust with felsic middle crust [[Suyehiro et al., 1996]](#ref-suyehiro1996) [[Takahashi et al., 2007]](#ref-takahashi2007) | Fixed |
| `IOA_MAX_HC_M` | 36,000 m | Oceanic arc stability limit [[Clift & Vannucchi, 2004]](#ref-clift2004) | Fixed |
| `IOA_POLARITY_MIN_AGE_GAP_MYR` | 10 Myr | Below this the density contrast is small and polarity is ambiguous [[Stern, 2002]](#ref-stern2002) | Sweep 5–30; watch the overlap metrics in `plate_diagnostics` |

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time. Report:
- juvenile additions split by source (continental arc, intra-oceanic arc, rift, volcanic);
- arc length by type;
- terrane area on oceanic plates, and terrane area docked per Gyr;
- plate overlap and plate count, since stopping retreat on the overriding side could let
  overlaps grow.

## Detailed design

- `backend/app/surface_fields.py`: register `crust_formed_years` (HISTORY, sentinel),
  with a comment citing this RFC.
- `backend/app/quad_tectonics.py`: stamp `crust_formed_years` on inserted ridge cells and
  on melted-through rift cells.
- `backend/app/gaps.py`: stamp it on gap-filled crust.
- `backend/app/lithosphere_plate.py`
  - New `oceanic_arc_band(plate, inputs, closing_rate, ages)`, next to
    `continental_arc_band`.
  - `boundary_context`: an oceanic overriding plate's contested nodes facing a younger-side
    contact are removed from `shrinkable`, and `arc_band` / `arc_intensity` are filled for
    oceanic plates.
  - `deform_columns`: apply intra-oceanic thickening, convert crust type past
    `IOA_CONTINENTAL_HC_M`, and book juvenile through
    `continental_ledger.add_material_thickness`.
- `backend/app/rheology.py`: `apply_oceanic_arc_thickening(hc, closing_rate, years_myr,
  intensity, band_width_km)` with the per-km conversion; constants with citations.
- `backend/app/continental_ledger.py`: a diagnostic sub-account `juvenile_oceanic_arc_m3`
  of `juvenile_additions_m3`, outside the balance.
- `backend/app/stats.py`: intra-oceanic arc length, juvenile by source, and terrane area on
  oceanic plates.
- `backend/app/render_image.py`: `elevReason` already has `ELEV_CHANGE_SUBDUCTION_ARC`.
  Stamp it on oceanic arc uplift too, so arcs show without a new view.
- Tests:
  - `test_lithosphere_plate.py` (new, mirroring the module): polarity picks the older
    side; the younger side doesn't retreat; an unknown age keeps symmetric retreat.
  - `test_rheology.py`: per-km productivity is density-independent.
  - `test_lithosphere_plate.py`: conversion at 20 km books exactly `Hc − 11 km` as
    juvenile; the cap holds at 36 km.
  - `test_persistence.py`: an old save loads with the sentinel.
- Docs: `docs/simulation-model.md` (Plate motion boundary classification, Volcanism, and
  Known simplifications for plateaus), and `docs/save-compatibility.md` for the new field.

## Acceptance criteria

- [ ] On seeds 7, 11 and 2024 to 1 Gyr, intra-oceanic arcs make 0.2–1.0 km³/yr of juvenile
      continental material, and at least one arc terrane docks onto a continent per run.
- [ ] Combined with RFC 0005, the continental inventory on seed 456862899 changes by less
      than 15% per Gyr over the first 3 Gyr, and never falls below 50% of its initial value
      over 5 Gyr.
- [ ] In a world stripped of continents (a test world with all continental crust removed
      at t = 0), juvenile additions stay above zero and genuine continental area grows over
      500 Myr.
- [ ] Plate overlap area fraction and plate count change by less than 15% against
      `origin/main`.
- [ ] Intra-oceanic arc productivity per km of margin is within 10% between density 1 and
      density 4 (unit test).
- [ ] Step time changes by less than 3%; save size by less than 3%.
- [ ] Ledgers close; old saves load; affected tests, ruff and `bin/quality_check.sh` are
      green.

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
