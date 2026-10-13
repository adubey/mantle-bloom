# RFC 0007: Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #205 (magma transport); RFCs [0003](0003-rifting-around-cratons.md)
  (plume weakening of cratons), [0004](0004-proto-cratons.md),
  [0005](0005-continental-recycling-balance.md), [0006](0006-intra-oceanic-arcs.md)
  (terranes, `crust_formed_years`)

## Problem

The model has no mantle plumes. Its mantle (`mantle.py`) is a fixed set of convection
centres that push or pull plates sideways; nothing rises through it. Volcanoes appear only
where rifts thin the crust (`lithosphere_plate._erupt_melted_nodes` and
`_ignite_early_rift_volcanoes`); `elevation_lines.py` notes that "there is no arc/hotspot
spawning yet". Four Earth processes are therefore missing:

1. **Oceanic plateaus.** Earth's largest volcanic provinces, such as Ontong Java, form when
   a plume head reaches oceanic lithosphere. They are thick enough to resist subduction,
   dock onto continents as terranes, and add about 1.1 km³/yr of crust to the continents
   [[Clift et al., 2009]](#ref-clift2009). The model has no such source, so RFC 0006's
   intra-oceanic arcs are its only route to new continent from the ocean.
2. **Continental flood basalts**, such as the Deccan and Siberian Traps: large, short-lived
   additions of crust to existing continents.
3. **Hotspot tracks**, such as the Hawaiian–Emperor chain: age-progressive volcanic chains
   that record plate motion over a fixed source.
4. **Plume-weakened lithosphere.** Cratons rift only after heating or melt weakens them
   [[Foley, 2008]](#ref-foley2008) [[Koptev et al., 2015]](#ref-koptev2015). RFC 0003
   leaves this out because there is nothing in the model to do the weakening.

The run behind RFCs 0001–0006 shows the cost of item 1. The run is seed 456862899, 19.2
Gyr, saved locally as `mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on
the author's machine and **not available for download**. Once its continents were gone,
juvenile continental additions fell to 0.4 km³/yr
(`bin/debug/ledger_history.py`). On Earth, plateau accretion alone would add about three
times that.

## Proposed solutions

**What happens in nature.**
- *Plumes.* Plumes rise from deep thermal boundary layers [[Morgan, 1971]](#ref-morgan1971).
  A starting plume has a large head followed by a narrow tail. The head flattens beneath
  the lithosphere and erupts flood basalt over an equant area 2,000–2,500 km across in a
  short time [[Campbell & Griffiths, 1990]](#ref-campbell1990). The tail then feeds a
  hotspot track for tens of Myr [[Richards et al., 1989]](#ref-richards1989).
- *How many.* Of about 49 catalogued hotspots, only 7 meet the criteria for a deep
  (core–mantle boundary) origin; the rest are shallower [[Courtillot et al., 2003]](#ref-courtillot2003).
- *Where.* Reconstructed large igneous provinces (LIPs) of the last 200–300 Myr erupted
  mostly above the edges of the two large low-velocity regions at the base of the mantle
  [[Burke & Torsvik, 2004]](#ref-burke2004). Plume sources are fixed in the deep mantle
  and long-lived.
- *How often and how big.* LIPs have occurred about once every 20 Myr since the Archean,
  or about once every 10 Myr once the oceanic record since 250 Ma is counted
  [[Bryan et al., 2010]](#ref-bryan2010). Each covers more than 0.1 M km² and typically
  erupts at least 1 M km³ [[Bryan et al., 2010]](#ref-bryan2010) in pulses of about 1–5
  Myr [[Bryan & Ernst, 2008]](#ref-bryan2008). Oceanic plateaus are the largest, up to
  2 M km² [[Kerr & Mahoney, 2007]](#ref-kerr2007). Ontong Java holds 44–57 M km³ of crust
  [[Gladczenko et al., 1997]](#ref-gladczenko1997).
- *Plateau structure and fate.*
  - Plateau crust is far thicker than the 6–7 km of normal oceanic crust: at least 30 km
    over much of Ontong Java, and 8–15 km for the Caribbean plateau. Plateaus stand 2–3 km
    above the surrounding abyssal floor [[Kerr & Mahoney, 2007]](#ref-kerr2007)
    [[Coffin & Eldholm, 1994]](#ref-coffin1994).
  - Their resistance to subduction depends on both thickness and age: a plateau that meets
    a trench soon after forming is more buoyant than an older one of the same thickness
    [[Cloos, 1993]](#ref-cloos1993) [[Kerr & Mahoney, 2007]](#ref-kerr2007).
  - The Caribbean plateau (8–15 km) resisted subduction because it reached a trench only a
    few Myr after forming. Ontong Java (≥30 km) resisted it 100 Myr after forming; only its
    lower-crustal edge has been subducted [[Mann & Taira, 2004]](#ref-mann2004).
  - Plateaus that resist dock onto the upper plate as terranes
    [[Kerr & Mahoney, 2007]](#ref-kerr2007) [[Tetreault & Buiter, 2014]](#ref-tetreault2014),
    and plume episodes have been tied to pulses of continental growth
    [[Stein & Hofmann, 1994]](#ref-stein1994).
- *Hotspot tracks.* The Hawaiian–Emperor chain's igneous volume flux peaks at about
  8 m³/s (≈ 0.25 km³/yr) and has varied by about a factor of two over its history
  [[Van Ark & Lin, 2004]](#ref-vanark2004). Hawaii is the strongest hotspot on Earth.
- *Plumes and rifting.* Plume heads beneath continents produce volcanic rifted margins and
  breakup [[White & McKenzie, 1989]](#ref-white1989). A plume beneath a craton deflects
  rifting to the craton's edges [[Koptev et al., 2015]](#ref-koptev2015), and reheating and
  melt weaken cratonic lithosphere [[Foley, 2008]](#ref-foley2008).

1. **A plume population in the mantle frame (recommended).** Plumes are discrete,
   long-lived objects fixed in the world (mantle) frame. Each has a head stage that
   emplaces a LIP, then a tail stage that builds a hotspot track as plates move over it.
   - On an oceanic host, the head builds a plateau that can dock as a terrane.
   - On a continental host, it adds flood basalt and weakens the lithosphere.
   - Tail-only plumes, the shallower majority of
     [[Courtillot et al., 2003]](#ref-courtillot2003), provide the remaining hotspots.
2. **Statistical LIP source only.** Add plateau crust at random places at the observed LIP
   rate, with no plume objects, tails or weakening. Simpler, but it gives none of items 3
   or 4 and no link between LIPs and the mantle flow. Rejected.
3. **Plumes from the convection model.** Turn `mantle.py`'s upwelling centres themselves
   into plumes. They are broad flow features (falloff 0.6–1.4 rad) with no head or tail,
   and there are only ~8 of them. Rejected as the plume model, but used below to place
   plumes, by analogy with [[Burke & Torsvik, 2004]](#ref-burke2004).

**Recommendation:** 1.

## Design

**State.** Plumes live in a world-level list `World.mantle_plumes`. Each `MantlePlume` is a
small dataclass:

| Field | Units | Meaning |
|---|---|---|
| `plume_id` | int | Stable id for events and the UI |
| `position` | unit vector, world frame | Fixed, like `mantle.ConvectionCenter.position`; plates move over it |
| `born_years` | yr | When the plume started |
| `head_volume_m3` | m³ | Head volume still to emplace (0 for a tail-only plume) |
| `head_end_years` | yr | When the head stage ends |
| `head_radius_km` | km | Head footprint radius |
| `tail_flux_m3_per_yr` | m³/yr | Tail volume flux |
| `tail_end_years` | yr | When the plume dies |

Three per-node surface fields are added:
- `plateau_crust_m` (EXTENSIVE, m): the plume-built share of Hc on oceanic hosts. Remaps
  conserve its volume, like `craton_crust_m`.
- `plateau_formed_years` (HISTORY, yr): when that crust was emplaced, for the age-dependent
  docking rule.
- `plume_heat_myr` (COUNTDOWN, Myr): time left on plume weakening.

**Spawning.** Each step draws events from a Poisson process, with a deterministic generator
keyed on `(seed, elapsed_years, "plume")` as other step RNGs are:
- head plumes at `PLUME_HEAD_RATE_PER_MYR`;
- tail-only plumes at `PLUME_TAIL_ONLY_RATE_PER_MYR`.

Positions are drawn from a density proportional to
`Σ_upwelling strength × exp(−(d/falloff)²)` over the upwelling convection centres, plus a
uniform floor `PLUME_UNIFORM_SHARE`. Plumes thus favour the mantle's upwellings, as LIPs
favour the margins of deep low-velocity regions [[Burke & Torsvik, 2004]](#ref-burke2004),
without being tied to them.
- Head volume is drawn log-uniformly in [`PLUME_HEAD_VOLUME_MIN_KM3`,
  `PLUME_HEAD_VOLUME_MAX_KM3`].
- Tail flux is drawn log-uniformly in [`PLUME_TAIL_FLUX_MIN_KM3_PER_YR`,
  `PLUME_TAIL_FLUX_MAX_KM3_PER_YR`].
- Every head plume also gets a tail, after its head stage.

**Head stage.** Over `PLUME_HEAD_DURATION_MYR`, the head emplaces
`head_volume × years / duration` per step. It is spread over the cells under the head
footprint with weights `w = exp(−(d / head_radius_km)²)` in great-circle distance d,
normalised by `Σ w·A`, so the volume placed is exact.
- **Oceanic host cells.** Hc and `plateau_crust_m` grow by the cell's share. Elevation
  follows by isostasy (`lithosphere.isostatic_elevation`, as `cratons.attrite` uses), which
  gives the 2–3 km plateau relief. `plateau_formed_years` is stamped. The crust stays
  oceanic-coded, since plateaus are mafic, not continental.
- **Continental host cells.** Hc grows by the share and is booked as juvenile continental
  material through `continental_ledger.add_material_thickness` (sub-account
  `juvenile_flood_basalt_m3`). Mobile cover under the lava is sealed, as
  `volcanism.apply_volcanic_activity` seals it.
- **All covered cells.** `plume_heat_myr` is set to `PLUME_HEAT_MYR`, and covered cells
  near the centre are made volcanoes (`is_volcano`, `volcano_active_years_remaining`), so
  the existing eruption lifecycle draws the surface.
- A plate whose continental cells take a head gets `PLUME_HEAD_STRESS` added to its
  `internal_stress`, the breakup pressure `merge_split.maybe_split_plate` reads. This
  follows plume-driven breakup [[White & McKenzie, 1989]](#ref-white1989).
- One event is logged per head: plate, volume, and oceanic or continental.

**Tail stage.** Until `tail_end_years`, the tail delivers `tail_flux × years` to the cells
within `PLUME_TAIL_RADIUS_KM`, with the same Gaussian weights.
- On oceanic hosts it goes to Hc and `plateau_crust_m` (seamounts and aseismic ridges); on
  continental hosts it is booked as juvenile.
- The cell nearest the plume is made a volcano, so the track is a chain of volcanoes left
  behind as the plate moves: a hotspot track with an age progression, for free.

**Weakening (for RFC 0003).** While `plume_heat_myr > 0`, craton strength in the rift
paths (`lithosphere_plate.deform_columns`'s divergent thinning,
`quad_tectonics._allocate_stretch`, failed rifts) is multiplied by
`1 − PLUME_CRATON_WEAKENING`. In RFC 0003's split-cut score, heated cells count as weak,
so cuts are drawn toward plume-heated lithosphere and cratons rift where a plume has
weakened them, after [[Foley, 2008]](#ref-foley2008) and
[[Koptev et al., 2015]](#ref-koptev2015). Erosion, attrition and the collision brake are
unchanged.

**Plateau docking.** In `quad_tectonics._retreat`, an oceanic cell being consumed counts as
a terrane donor rather than subducting when:

`plateau_crust_m ≥ PLATEAU_DOCK_HC_M`, or
(`plateau_crust_m ≥ PLATEAU_DOCK_YOUNG_HC_M` and age < `PLATEAU_YOUNG_MYR`).

Here age = `elapsed_years − plateau_formed_years`. This is the thickness-and-age rule of
[[Cloos, 1993]](#ref-cloos1993) and [[Kerr & Mahoney, 2007]](#ref-kerr2007). Donor cells go
through the existing continental-terrane path (#253, #305, #321) to the upper plate. When
the receiving plate is continental, the docked plateau crust is recoded continental and
booked as juvenile (`juvenile_plateau_m3`). Otherwise it stays plateau crust on the new
host and can dock again later.

A plateau's thick root can also feed RFC 0004's keel path: Ontong Java is underlain by a
mantle root over 300 km thick [[Kerr & Mahoney, 2007]](#ref-kerr2007). No extra rule is
needed, because Hm is shared.

**Global scale check.** At one head per 10 Myr with volumes log-uniform in [1, 50] M km³
(mean ≈ 12.5 M km³), heads deliver about 1.25 km³/yr on long-run average. That matches the
~1.1 km³/yr of plateau accretion in [[Clift et al., 2009]](#ref-clift2009), assuming most
plateau crust eventually docks or is subducted. With ~25 active tails at a mean of ~0.06
km³/yr, tails deliver about 1.5 km³/yr, mostly to the ocean floor.

**Runtime.** Each active plume tests the cells of each plate it might touch: a dot product
against its position, after a bounding-cap test per plate. With ≤ 50 plumes and ~1.3 × 10⁵
nodes at density 4, that is ~6 × 10⁶ multiply-adds per step, a few ms. Head emplacement
happens a few steps per 10 Myr.

**Save compatibility.**
- `World.mantle_plumes` is a `default_factory` list, backfilled in
  `persistence._backfill_added_fields`. An old save starts with no plumes and spawns them
  from its next step.
- The three surface fields are registered in `surface_fields.py` with defaults (0, sentinel,
  0), which the registry backfills.
- `MantlePlume` gets defaults on every field, so later additions stay loadable (see
  `docs/save-compatibility.md`).

**Live knob.** `plume_multiplier` joins `World.TUNING_MULTIPLIER_FIELDS` (default 1.0) and
scales spawn rates. At 0, no plume ever spawns, which reproduces `origin/main`. Debug and
sketch worlds can turn plumes off that way.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `PLUME_HEAD_RATE_PER_MYR` | 0.1 | One LIP per ~10 Myr including oceanic LIPs [[Bryan et al., 2010]](#ref-bryan2010) | Fixed; 0.05 (one per 20 Myr) is the sensitivity case |
| `PLUME_HEAD_VOLUME_MIN_KM3`, `_MAX_KM3` | 1 × 10⁶, 5 × 10⁷ | ≥ 1 M km³ typical [[Bryan et al., 2010]](#ref-bryan2010); Ontong Java 44–57 M km³ [[Gladczenko et al., 1997]](#ref-gladczenko1997) | Fixed |
| `PLUME_HEAD_RADIUS_KM` | 1,100 km | Heads 2,000–2,500 km across [[Campbell & Griffiths, 1990]](#ref-campbell1990) | Fixed |
| `PLUME_HEAD_DURATION_MYR` | 3 Myr | Pulses of ~1–5 Myr [[Bryan & Ernst, 2008]](#ref-bryan2008) | Fixed; must exceed one step |
| `PLUME_TAIL_FLUX_MIN/MAX_KM3_PER_YR` | 0.01, 0.25 | Hawaii, the strongest hotspot, peaks near 0.25 [[Van Ark & Lin, 2004]](#ref-vanark2004) | Fixed |
| `PLUME_TAIL_LIFETIME_MYR` | 80 Myr | Length of the Hawaiian–Emperor record [[Van Ark & Lin, 2004]](#ref-vanark2004) | Sweep 40–120 |
| `PLUME_TAIL_RADIUS_KM` | 150 km | A volcanic edifice and its apron; at least two cells at density 1 | Fixed |
| `PLUME_TAIL_ONLY_RATE_PER_MYR` | 0.3 | Sets ~25 active hotspots with the lifetime above; Earth has ~49, 7 deep [[Courtillot et al., 2003]](#ref-courtillot2003) | Tune so 20–50 hotspots are active |
| `PLUME_UNIFORM_SHARE` | 0.3 | Not every plume sits on an upwelling [[Burke & Torsvik, 2004]](#ref-burke2004) | Fixed |
| `PLATEAU_DOCK_HC_M` | 25,000 m | Between Caribbean (8–15 km) and Ontong Java (≥ 30 km), which resisted at any age [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Sweep 20–30 km; target the docked share below |
| `PLATEAU_DOCK_YOUNG_HC_M`, `PLATEAU_YOUNG_MYR` | 10,000 m, 10 Myr | The young Caribbean case [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Fixed |
| `PLUME_HEAT_MYR` | 20 Myr | The order of thermal re-equilibration of heated lithosphere; no direct constraint | Tune with RFC 0003: craton rifting should cluster on plume-heated cells |
| `PLUME_CRATON_WEAKENING` | 0.75 | Cancels `CRATON_RIFT_RESISTANCE` at full heat | Sweep 0.5–1.0 |
| `PLUME_HEAD_STRESS` | in `internal_stress` units | Breakup often follows a continental plume head [[White & McKenzie, 1989]](#ref-white1989); no rate is given | Tune to a design target: about half of continental plume heads followed by a split of their plate within 20 Myr |

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time. Report:
- LIP count and volumes per Gyr, with the oceanic/continental split;
- active hotspot count;
- plateau crust docked versus subducted;
- juvenile continental by source (flood basalt, docked plateau, tail);
- craton rifted volume inside versus outside plume-heated cells;
- step time.

## Detailed design

- `backend/app/plumes.py` (new):
  - `MantlePlume`;
  - `spawn(world, years)`;
  - `apply(world, years)` (head, tail and heat, per plate);
  - `docking_mask(plate, elapsed_years)`;
  - `diagnostics(world)`;
  - constants with citations.
- `backend/app/world.py`: the `mantle_plumes` field and `plume_multiplier` knob; call
  `plumes.spawn` then `plumes.apply` in `step_world_progress` after topology changes and
  before `volcanism.apply_volcanic_activity`, so new volcanoes roll the same step.
- `backend/app/persistence.py`: backfill `mantle_plumes`.
- `backend/app/surface_fields.py`: register `plateau_crust_m`, `plateau_formed_years` and
  `plume_heat_myr`.
- `backend/app/quad_tectonics.py`
  - `_retreat`: include `plumes.docking_mask` in the `terrane` donors.
  - At the point docked crust lands on a continental upper plate (`_accrete_onto_survivors`
    and `crust_transfer`): recode it and book `juvenile_plateau_m3`.
  - `_allocate_stretch`: plume weakening on craton strength.
- `backend/app/lithosphere_plate.py`, `deform_columns`: plume weakening on
  `CRATON_RIFT_RESISTANCE`.
- `backend/app/sparse_quad_patch.py`, `apply_failed_rift`: the same weakening (with RFC
  0003's resistance).
- `backend/app/merge_split.py`: RFC 0003's strength proxy counts heated cells as weak; add
  the head stress to `internal_stress`.
- `backend/app/continental_ledger.py`: sub-accounts `juvenile_flood_basalt_m3`,
  `juvenile_plateau_m3` and `juvenile_hotspot_m3` under `juvenile_additions_m3`.
- `backend/app/elevation_lines.py`: correct the "no hotspot spawning yet" comment.
- `backend/app/render_image.py`: a `plumes` debug view (plume positions, head footprints,
  tail age, plateau crust thickness), added to `DEBUG_VIEWS`; `frontend/src/api.ts` view
  union and `legendData.ts` legend.
- `backend/app/stats.py`: active plumes, LIP events, plateau area and volume.
- `backend/app/main.py`: `GET /world/plumes` listing live plumes for the Plate Inspector or
  a future panel. Optional; the debug view alone meets the criteria.
- Tests:
  - `backend/unit_tests/test_plumes.py` (new): spawn counts match the rate over many draws,
    deterministically; head volume placed equals the scheduled volume; a plate moving over
    a fixed tail leaves volcanoes in a line parallel to its motion; the docking rule
    (thick, young-thin, old-thin); weakening halves effective craton resistance at 0.5;
    `plume_multiplier = 0` changes nothing.
  - `test_persistence.py`: an old save loads with no plumes and default fields.
  - `test_continental_ledger.py`: the juvenile sub-accounts sum correctly.
- Docs:
  - `docs/simulation-model.md`: a new "Mantle plumes (`plumes.py`)" section; update the
    Volcanism, Cratons (Resistance) and Known simplifications sections.
  - `docs/architecture.md` (module map), `docs/api-reference.md` (view, endpoint, stats
    keys) and `docs/save-compatibility.md` (new fields).

## Acceptance criteria

- [ ] Over 1 Gyr on each of seeds 7, 11 and 2024: 70–130 plume heads, a volume
      distribution matching the configured range, and between 20 and 50 active hotspots at
      every 100 Myr checkpoint.
- [ ] Hotspot tracks: on a pinned-motion test world, volcanoes created by one tail lie
      within 20° of the plate's motion direction at the plume and get younger toward it.
- [ ] Plateaus reach 2–3 km above the surrounding abyssal floor in the median head, and
      plateaus thicker than `PLATEAU_DOCK_HC_M` dock rather than subduct (unit test plus a
      replay count).
- [ ] Long-run juvenile continental material from docked plateaus and flood basalts is
      0.5–1.5 km³/yr on the tuning seeds, against ~1.1 km³/yr on Earth
      [[Clift et al., 2009]](#ref-clift2009).
- [ ] With RFC 0003 in place, at least half of the craton volume lost to rifting is lost in
      plume-heated cells.
- [ ] Step time changes by less than 3% at density 4; save size by less than 3%.
- [ ] `plume_multiplier = 0` gives step results identical to `origin/main` on a 10-step
      replay.
- [ ] Ledgers close; old saves load; affected tests, ruff and `bin/quality_check.sh` are
      green.

## References

- <a id="ref-bryan2008"></a>Bryan, S. E., & Ernst, R. E. (2008). Revised definition of
  Large Igneous Provinces (LIPs). *Earth-Science Reviews*, 86, 175–202.
  [doi:10.1016/j.earscirev.2007.08.008](https://doi.org/10.1016/j.earscirev.2007.08.008)
- <a id="ref-bryan2010"></a>Bryan, S. E., Ukstins Peate, I., Peate, D. W., Self, S.,
  Jerram, D. A., Mawby, M. R., Marsh, J. S., & Miller, J. A. (2010). The largest volcanic
  eruptions on Earth. *Earth-Science Reviews*, 102, 207–229.
  [doi:10.1016/j.earscirev.2010.07.001](https://doi.org/10.1016/j.earscirev.2010.07.001)
- <a id="ref-burke2004"></a>Burke, K., & Torsvik, T. H. (2004). Derivation of Large
  Igneous Provinces of the past 200 million years from long-term heterogeneities in the
  deep mantle. *Earth and Planetary Science Letters*, 227, 531–538.
  [doi:10.1016/j.epsl.2004.09.015](https://doi.org/10.1016/j.epsl.2004.09.015)
- <a id="ref-campbell1990"></a>Campbell, I. H., & Griffiths, R. W. (1990). Implications of
  mantle plume structure for the evolution of flood basalts. *Earth and Planetary Science
  Letters*, 99, 79–93.
  [doi:10.1016/0012-821X(90)90072-6](https://doi.org/10.1016/0012-821X(90)90072-6)
- <a id="ref-clift2009"></a>Clift, P. D., Vannucchi, P., & Morgan, J. P. (2009). Crustal
  redistribution, crust–mantle recycling and Phanerozoic evolution of the continental
  crust. *Earth-Science Reviews*, 97, 80–104.
  [doi:10.1016/j.earscirev.2009.10.003](https://doi.org/10.1016/j.earscirev.2009.10.003)
- <a id="ref-cloos1993"></a>Cloos, M. (1993). Lithospheric buoyancy and collisional
  orogenesis: Subduction of oceanic plateaus, continental margins, island arcs, spreading
  ridges, and seamounts. *Geological Society of America Bulletin*, 105, 715–737.
  [doi:10.1130/0016-7606(1993)105<0715:LBACOS>2.3.CO;2](https://doi.org/10.1130/0016-7606(1993)105%3C0715:LBACOS%3E2.3.CO;2)
- <a id="ref-coffin1994"></a>Coffin, M. F., & Eldholm, O. (1994). Large igneous provinces:
  Crustal structure, dimensions, and external consequences. *Reviews of Geophysics*, 32,
  1–36. [doi:10.1029/93RG02508](https://doi.org/10.1029/93RG02508)
- <a id="ref-courtillot2003"></a>Courtillot, V., Davaille, A., Besse, J., & Stock, J.
  (2003). Three distinct types of hotspots in the Earth's mantle. *Earth and Planetary
  Science Letters*, 205, 295–308.
  [doi:10.1016/S0012-821X(02)01048-8](https://doi.org/10.1016/S0012-821X(02)01048-8)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-gladczenko1997"></a>Gladczenko, T. P., Coffin, M. F., & Eldholm, O. (1997).
  Crustal structure of the Ontong Java Plateau: Modeling of new gravity and existing
  seismic data. *Journal of Geophysical Research: Solid Earth*, 102, 22711–22729.
  [doi:10.1029/97JB01636](https://doi.org/10.1029/97JB01636)
- <a id="ref-kerr2007"></a>Kerr, A. C., & Mahoney, J. J. (2007). Oceanic plateaus:
  Problematic plumes, potential paradigms. *Chemical Geology*, 241, 332–353.
  [doi:10.1016/j.chemgeo.2007.01.019](https://doi.org/10.1016/j.chemgeo.2007.01.019)
- <a id="ref-koptev2015"></a>Koptev, A., Calais, E., Burov, E., Leroy, S., & Gerya, T.
  (2015). Dual continental rift systems generated by plume–lithosphere interaction.
  *Nature Geoscience*, 8, 388–392. [doi:10.1038/ngeo2401](https://doi.org/10.1038/ngeo2401)
- <a id="ref-mann2004"></a>Mann, P., & Taira, A. (2004). Global tectonic significance of
  the Solomon Islands and Ontong Java Plateau convergent zone. *Tectonophysics*, 389,
  137–190. [doi:10.1016/j.tecto.2003.10.024](https://doi.org/10.1016/j.tecto.2003.10.024)
- <a id="ref-morgan1971"></a>Morgan, W. J. (1971). Convection plumes in the lower mantle.
  *Nature*, 230, 42–43. [doi:10.1038/230042a0](https://doi.org/10.1038/230042a0)
- <a id="ref-richards1989"></a>Richards, M. A., Duncan, R. A., & Courtillot, V. E. (1989).
  Flood basalts and hot-spot tracks: Plume heads and tails. *Science*, 246, 103–107.
  [doi:10.1126/science.246.4926.103](https://doi.org/10.1126/science.246.4926.103)
- <a id="ref-stein1994"></a>Stein, M., & Hofmann, A. W. (1994). Mantle plumes and episodic
  crustal growth. *Nature*, 372, 63–68.
  [doi:10.1038/372063a0](https://doi.org/10.1038/372063a0)
- <a id="ref-tetreault2014"></a>Tetreault, J. L., & Buiter, S. J. H. (2014). Future
  accreted terranes: A compilation of island arcs, oceanic plateaus, submarine ridges,
  seamounts, and continental fragments. *Solid Earth*, 5, 1243–1275.
  [doi:10.5194/se-5-1243-2014](https://doi.org/10.5194/se-5-1243-2014)
- <a id="ref-vanark2004"></a>Van Ark, E., & Lin, J. (2004). Time variation in igneous
  volume flux of the Hawaii-Emperor hot spot seamount chain. *Journal of Geophysical
  Research: Solid Earth*, 109, B11401.
  [doi:10.1029/2003JB002949](https://doi.org/10.1029/2003JB002949)
- <a id="ref-white1989"></a>White, R., & McKenzie, D. (1989). Magmatism at rift zones: The
  generation of volcanic continental margins and flood basalts. *Journal of Geophysical
  Research: Solid Earth*, 94, 7685–7729.
  [doi:10.1029/JB094iB06p07685](https://doi.org/10.1029/JB094iB06p07685)
