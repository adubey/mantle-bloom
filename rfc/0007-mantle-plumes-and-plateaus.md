# RFC 0007: Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks

- **Status:** Draft
- **Date:** 2026-10-13; evidence refreshed 2026-10-14 (#403)
- **Related:** #205 (magma transport),
  [#407](https://github.com/adubey/mantle-bloom/issues/407) (source attribution and breakup baseline); RFCs [0010](0010-rifts-follow-weak-lithosphere.md)
  (plume weakening of cratons and conditional plume-assisted local rifting; supersedes 0003), [0009](0009-craton-mosaic.md) (cratons;
  supersedes 0004, keeps its depletion follow-up deferred),
  [0005](0005-continental-recycling-balance.md), [0006](0006-intra-oceanic-arcs.md)
  (shared docking, buoyant-crust tracers and ledger),
  [0012](0012-physics-first-simulation.md) (force-response interface for plume loads)

## Problem

The model has no mantle plumes. Its mantle (`mantle.py`) is a fixed set of convection
centres that push or pull plates sideways; nothing rises through it. Volcanoes appear only
where rifts thin the crust (`lithosphere_plate._erupt_melted_nodes` and
`_ignite_early_rift_volcanoes`); `elevation_lines.py` notes that "there is no arc/hotspot
spawning yet". Four Earth processes are therefore missing:

1. **Oceanic plateaus.** Earth's largest volcanic provinces, such as Ontong Java, form when
   a plume head reaches oceanic lithosphere. They are thick enough to resist subduction,
   dock onto continents as terranes, and add about 1.1 km³/yr of crust to the continents
   [[Clift et al., 2009]](#ref-clift2009). The model has no plume-plateau source. RFC 0006 proposes intra-oceanic arcs as another
   oceanic source; the existing elevation-driven reset below is a separate model route.
2. **Continental flood basalts**, such as the Deccan and Siberian Traps: large, short-lived
   additions of crust to existing continents.
3. **Hotspot tracks**, such as the Hawaiian–Emperor chain: age-progressive volcanic chains
   that record plate motion over a fixed source.
4. **Plume-weakened lithosphere.** Heating and melt can weaken lithosphere and influence some
   rifts [[Foley, 2008]](#ref-foley2008) [[Koptev et al., 2015]](#ref-koptev2015). RFC 0010
   proposes this as a conditional local contribution, but its production rift mechanics and plume
   coupling have not yet passed Gate 3. This RFC supplies the plume field; it does not make a plume
   necessary for every rift or let a plume event force a plate-wide split.

**Historical evidence, pre-#358.** The run behind RFCs 0001–0007 motivated item 1. The
run is seed 456862899, 19.2 Gyr, saved locally as
`mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine and
**not available for download**. Once its continents had almost disappeared,
juvenile continental additions fell to 0.4 km³/yr
(`bin/debug/ledger_history.py`). This claim was first recorded in RFC publication
commit `d795a09` (2026-10-13); the generating commit, density and controls remain unknown.
#403 records 0.1 Myr steps through 21.8 Myr (218 steps), then 10 Myr steps. The old 20-step
cooldown therefore lasted 2 Myr initially and 200 Myr for most of the run. The save
predates #358 (`7afecfd`, merged at `1b4f686`), which made the cooldown a physical 2 Myr.
The published late flux is preserved without recomputing the private save and is not
a current-main baseline. Earth's estimated plateau accretion is about three times that
historical flux, but the quantities have different process scopes and do not calibrate plume production.

The historical run is a stress test of long-term behaviour, not a calibration interval.
The absence of plume emplacement is a current code finding. A broader claim that there
is no route to continent without existing continental crust needs qualification:
`lithosphere_plate._erupt_melted_nodes` recodes any column melting through above 0 m as
continental, including an exposed oceanic-coded one, and its callers book juvenile
material. This existing route also matters if a proposed plume plateau later melts
through. Neither the old end state nor the current-main refresh below measures that
route or says how much plumes should add.

**Current-main baseline, post-#358 (`2562f63`, #403).**
[The fresh evidence](../analysis/issue403/README.md) regenerates default quad worlds at
density 1, using 5 Myr steps, 100 Myr checkpoints, a physical 2 Myr split cooldown and no
control overrides. [Tuning results](../analysis/issue403/tuning.json) record the inputs
and ledger checkpoints; [RFC 0005](0005-continental-recycling-balance.md) reports the
continental fluxes. Seeds 7, 11 and 2024 lose 86.0%, 73.0% and 74.8% of continental
inventory in 1 Gyr, while juvenile additions average 6.85, 6.05 and 4.22 km³/yr. On seed 7
they fall from 10.70 at 0–100 Myr to 4.96 at 500–1,000 Myr. The old 0.4 km³/yr end-state rate
covers a different late interval under the pre-#358 schedule; neither it nor the new
early-window rates calibrate a plume source. The current code contains no plume
population, plateau tracer or buoyant-crust ledger, so this refresh supplies no plume production, docking or preservation measurement.
The missing plume mechanism remains a design motivation; the magnitude of preserved
plume growth stays subject to the isolated-head and source-attribution gates, sensitivity
runs and falsification criteria.

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
  [[Bryan et al., 2010]](#ref-bryan2010). Both figures are averages over an incomplete
  record: older and oceanic LIPs are lost to erosion and subduction. LIPs also cluster in
  time rather than arriving at a steady rate, and crustal growth may itself be episodic
  [[Stein & Hofmann, 1994]](#ref-stein1994). They constrain the order of magnitude, not a
  stationary arrival process. Each covers more than 0.1 M km² and typically
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

**Recommendation:** 1 remains the proposed mechanism after #403's current-code check.
The isolated-head scale gate still comes before tuning, and the refreshed continental
budget does not select plume rate, plateau preservation or head stress to close a loss
measured in the pre-#358 save. Those calibration assumptions retain their sensitivity
runs and falsification criteria. Before interpreting plume-driven continental growth,
the source-attribution gate below must distinguish direct plume emplacement and plateau
docking from the existing elevation-driven reset route.
[#407](https://github.com/adubey/mantle-bloom/issues/407) tracks that diagnosis alongside
RFC 0006's [#406](https://github.com/adubey/mantle-bloom/issues/406). The plume and docking
design stay unchanged, and this refresh measures no plume output.

## Design

**Conventions shared with RFC 0006.**
- *Excess crust.* On oceanic-coded cells, E = Hc − `REFERENCE_HC_OCEANIC_M` (11 km). The
  literature gives total crustal thicknesses over a normal oceanic crust of 6–7 km
  [[Kerr & Mahoney, 2007]](#ref-kerr2007), so totals are converted to excess by
  subtracting 6.5 km.
- *Buoyant crust.* "Buoyant crust" is the plume-built `plateau_crust_m` plus RFC 0006's
  `arc_crust_m`. Both are tracked shares of Hc with one shared lifecycle and ledger
  (below).

**State.**
- `World.mantle_plumes`: a list of `MantlePlume` dataclasses.

  | Field | Units | Meaning |
  |---|---|---|
  | `plume_id` | int | The event index k (below); stable id for events and the UI |
  | `position` | unit vector, world frame | Fixed, like `mantle.ConvectionCenter.position`; plates move over it |
  | `born_years` | yr | The event's arrival time |
  | `head_volume_m3`, `head_remaining_m3` | m³ | Head volume, and the part not yet emplaced (0 for a tail-only plume) |
  | `head_end_years` | yr | `born_years + PLUME_HEAD_DURATION_MYR` |
  | `tail_flux_m3_per_yr`, `tail_end_years` | m³/yr, yr | Tail volume flux, and when the plume dies |

- `World.plume_stream`: the next event index and its arrival time, per stream (head,
  tail-only).
- `World.buoyant_crust_ledger`: one ledger shared with RFC 0006.
- Three new per-node surface fields:
  - `plateau_crust_m` (EXTENSIVE, m): the plume-built share of Hc.
  - `buoyant_crust_time_m_yr` (EXTENSIVE, m·yr): Σ thickness × emplacement time over the
    cell's buoyant crust. The cell's volume-weighted emplacement time is this divided by
    `plateau_crust_m + arc_crust_m`. It is conserved by remaps, unlike a HISTORY date,
    which keeps only the oldest contributor.
  - `plume_heat_myr` (COUNTDOWN, Myr): time left on plume weakening.

**An event stream independent of step length.** Each event's draws are keyed by its
*index*, not by the step it falls in, so ten 1 Myr steps and one 10 Myr step produce the
same plumes:
- *Arrivals.* Arrivals of stream s form a continuous-time process at a ceiling rate
  λ_max,s. Event k arrives at t_k = t_(k−1) + g_k, where g_k is exponentially distributed
  with mean 1/λ_max,s, drawn from a generator keyed on `(seed, "plume", s, k)`.
- *Thinning.* Event k is kept if `u_k < λ_s(t_k) / λ_max,s`, with u_k from the same
  generator. Here λ_s(t) = base rate × `plume_multiplier` × modulation(t). The modulation
  is 1 for the stationary default. For the Cox variant it is a deterministic function of
  t, built from a seed-keyed log-normal noise series with correlation time
  `PLUME_RATE_CORRELATION_MYR` and clipped to at most `PLUME_RATE_MODULATION_MAX`, since
  thinning needs a finite ceiling. λ_max,s is the base rate × the knob's upper bound ×
  that maximum.
- *Attributes.* Position, volume and tail flux come from the same generator.
- *Each step [t₀, t₁]* consumes every event with t_k ≤ t₁, born at t_k, and stores the
  next index and arrival time in `World.plume_stream`.

Kept events at `plume_multiplier` m are a subset of those at any larger multiplier, so
changing the knob doesn't reshuffle the rest of the run.

*Placement.* Positions are drawn from a density proportional to
`Σ_upwelling strength × exp(−(d/falloff)²)` over the upwelling convection centres, plus a
uniform floor `PLUME_UNIFORM_SHARE`. Plumes thus favour the mantle's upwellings, as LIPs
favour the margins of deep low-velocity regions [[Burke & Torsvik, 2004]](#ref-burke2004),
without being tied to them.

*Sizes.* Head volume is drawn log-uniformly in [`PLUME_HEAD_VOLUME_MIN_KM3`,
`PLUME_HEAD_VOLUME_MAX_KM3`], and tail flux log-uniformly in
[`PLUME_TAIL_FLUX_MIN_KM3_PER_YR`, `PLUME_TAIL_FLUX_MAX_KM3_PER_YR`]. Every head plume gets
a tail after its head stage.

The rate and the volume distribution are **calibration assumptions** anchored to the LIP
record's order of magnitude, not values it establishes (see the knob table). Three
sensitivity variants are run:
- **Rate:** 0.05, 0.1 and 0.2 heads per Myr.
- **Clustering:** the Cox variant, against the stationary default.
- **Volumes:** log-uniform (default), log-normal with the same median and range, and a
  power law, to show how much the size distribution's tail matters.

**Timestep integration.** Every plume stage is an interval: the head is
[`born_years`, `head_end_years`], the tail [`head_end_years`, `tail_end_years`], and heat
counts down over [emplacement, + `PLUME_HEAT_MYR`]. For a step covering [t₀, t₁], a stage
delivers its rate times the *overlap* of the two intervals, never `rate × years`:

`delivered = min(remaining, V × |[t₀, t₁] ∩ [t_start, t_end]| / (t_end − t_start))`

where `remaining` is the stage's undelivered volume, which then decreases by `delivered`.
- A 10 Myr step containing a whole 3 Myr head delivers exactly V. A head that starts late
  in one step and ends in the next splits V between them.
- Total delivery can never exceed V, whatever the step lengths.
- Plumes whose tail has ended are dropped at the end of the step.

With the index-keyed event stream, the same seed gives the same events and the same
delivered volumes under any step partition, except where plate positions differ between
partitions (which changes *where* the volume lands, not how much).

**Head footprint: relief consistent with volume.** Spreading a median head over the 2,000–
2,500 km flood-basalt area [[Campbell & Griffiths, 1990]](#ref-campbell1990) would be
inconsistent with plateau relief.
- *The problem.* The median log-uniform head (V ≈ 7.1 M km³) over a Gaussian of radius
  1,100 km (∫ w dA ≈ 2.4 M km²) adds only ~2.9 km of crust at the centre. The model's
  water-loaded oceanic isostasy turns that into ~0.46 km of relief, against the 2–3 km
  plateaus stand above the abyssal floor [[Kerr & Mahoney, 2007]](#ref-kerr2007).
- *The fix.* The footprint depends on where the head lands, judged at the host cell under
  the plume centre at `born_years`. The oceanic/continental branch follows that host
  cell, so two step partitions whose plates have drifted apart can branch the same
  scheduled event differently; the timestep acceptance below measures that branch
  stability on a pinned-motion test.
- **Oceanic: a plateau.** The area is set by volume:
  `A = clip(V / PLATEAU_MEAN_EXCESS_M, PLATEAU_AREA_MIN_KM2, PLATEAU_AREA_MAX_KM2)`.
  Weights are flat-topped, `w = exp(−(d/r)⁴)`, with r chosen so the kernel's own integral
  ∫ w dA equals A, i.e. `r = sqrt(A / (π Γ(3/2)))` ≈ `sqrt(A / 2.78)`. The peak excess is
  then exactly V/A, and the volume placed is exact after normalising over the cells.
  - In the model's isostasy, each km of excess oceanic crust raises the surface by
    ~155 m. `PLATEAU_MEAN_EXCESS_M` = 15 km therefore gives ~2.3 km of relief, inside the
    2–3 km range.
  - With the area bounds (0.1 to 2 M km², [[Bryan et al., 2010]](#ref-bryan2010)
    [[Kerr & Mahoney, 2007]](#ref-kerr2007)), a 1 M km³ head covers 0.1 M km² at ~10 km
    excess (~1.6 km relief). The median 7.1 M km³ head covers 0.47 M km² (radius ~410 km)
    at 15 km peak (~2.3 km). A 50 M km³ head is capped at 2 M km² and averages 25 km excess
    (~3.9 km); that overshoots Ontong Java's relief, whose crust is ≥ 30 km total
    (≥ 23.5 km excess).
- **Continental: flood basalt** over the Gaussian of radius `PLUME_HEAD_RADIUS_KM`
  (1,100 km), matching the 2,000–2,500 km eruptive area. Crust added on a continent
  includes lava and underplated intrusions; the median head adds ~2.9 km at the centre.
- Cells under a footprint get their share on their own host type: oceanic-coded cells
  gain Hc and `plateau_crust_m`, continental-coded cells gain Hc and juvenile continental
  material.
- **Heat** (weakening) always uses the broad head radius, since the plume head heats a
  wide region whatever erupts.

The isolated-head scale check and the source-attribution gate are both required before
any tuning replay.

**Head and tail emplacement.**
- *Oceanic host cells:* Hc and `plateau_crust_m` grow by the cell's share, and
  `buoyant_crust_time_m_yr` by share × emplacement time. Elevation follows by isostasy
  (`lithosphere.isostatic_elevation`, as `cratons.attrite` uses).
- *Continental host cells:* Hc grows by the share and is booked as juvenile continental
  material through `continental_ledger.add_material_thickness` (sub-account
  `juvenile_flood_basalt_m3` for heads, `juvenile_hotspot_m3` for tails). Mobile cover
  under the lava is sealed, as `volcanism.apply_volcanic_activity` seals it.
- *Tails* deliver to the cells within `PLUME_TAIL_RADIUS_KM` with the same flat-topped
  weights: a volcanic edifice on whatever cell is over the plume, left behind as the plate
  moves, which is a hotspot track.
- *Stamps.* Every cell that gains crust is stamped `ELEV_CHANGE_VOLCANO`.
- **One volume owner.** `plumes.py` owns all plume volume. It does **not** set
  `is_volcano` on plume cells. `volcanism.apply_volcanic_activity` adds Hc on every
  eruption (`lithosphere.back_elevation_gain_fields`) and books juvenile on continental
  hosts, so marking plume cells as volcanoes would add an unbudgeted second magma source,
  and on oceanic hosts it would also miss `plateau_crust_m`. Plume relief and tracks come
  from the plume's own budget alone. (Mineral deposits, which `volcanism.py` grows per
  eruption, are not grown on plume cells; out of scope.)
- *Breakup stress.* A plate whose continental cells take a head gets `PLUME_HEAD_STRESS`
  added to its `internal_stress`, the breakup pressure `merge_split.maybe_split_plate`
  reads, after plume-driven breakup [[White & McKenzie, 1989]](#ref-white1989).
- *Events.* One event is logged per head: plate, volume, and oceanic or continental.

**Weakening (for RFC 0010).** While `plume_heat_myr > 0`, craton strength in the rift
paths (`lithosphere_plate.deform_columns`'s divergent thinning,
`quad_tectonics._allocate_stretch`, failed rifts) is multiplied by
`1 − PLUME_CRATON_WEAKENING`, with the product clamped to [0, 1] as `cratons.rift_resistance`
already clamps it (#385). Heat integrates by interval overlap like head and tail volume:
`plume_heat_myr` decreases by the step's elapsed years, and when a step extends beyond
the remaining heat, that step's weakening is scaled by `remaining_heat / step_years`
instead of applying full weakening for the whole step. In RFC 0010's cut strength, heated cells count as weak, so cuts are drawn toward
plume-heated lithosphere and cratons rift where a plume has weakened them, after
[[Foley, 2008]](#ref-foley2008) and [[Koptev et al., 2015]](#ref-koptev2015). Erosion,
attrition and the collision brake are unchanged.

**Docking: probabilistic and partial.** Thick, young buoyant crust is *more likely* to
resist subduction [[Cloos, 1993]](#ref-cloos1993) [[Kerr & Mahoney, 2007]](#ref-kerr2007),
but it isn't guaranteed. And even a plateau that resists loses part of itself: Ontong
Java's lower-crustal edge has gone down the Solomon trench [[Mann & Taira, 2004]](#ref-mann2004).
- **Patches.** At each step's retreat, consumable oceanic-coded cells with buoyant crust
  are grouped into connected patches. Each patch has an area-weighted mean excess E (the
  whole column's excess, not the tracer, so it matches how the literature measures
  thickness) and a volume-weighted buoyant-crust age, from `buoyant_crust_time_m_yr`.
- **Threshold with age.** In excess terms, the thickness needed to resist rises as the
  patch ages and cools:
  `E₅₀(age) = E_young + (E_old − E_young) × (1 − exp(−age / τ_cool))`.
  - *E_young* = 5 km: the young Caribbean plateau resisted with 8–15 km total crust, which
    is 1.5–8.5 km excess.
  - *E_old* = 15 km: old plateaus with Caribbean thickness would subduct, while Ontong Java
    (≥ 30 km total, ≥ 23.5 km excess) resisted 100 Myr after forming
    [[Kerr & Mahoney, 2007]](#ref-kerr2007).
- **Probability.** `P = σ((E − E₅₀(age)) / ΔE)`, with σ the logistic function, so P is in
  [0, 1]. For example:
  - an Ontong-Java-like patch (E = 23.5 km, 100 Myr old) gets P ≈ 0.95;
  - an old Caribbean-like one (E = 5 km) gets P ≈ 0.04;
  - an active arc, whose crust is young, docks readily.
- **The draw.** One deterministic draw per patch, keyed on (seed, plate, the patch's
  lowest cell key, volume-weighted buoyant-crust emplacement time quantized to 1 Myr),
  decides once per encounter; the result is kept until the patch leaves the trench.
  Keying on emplacement time rather than encounter start time keeps the draw stable
  when different step partitions first detect the same patch on different steps.
- **Dock.** A docking patch goes through the existing terrane path in `_retreat` (#253,
  #305, #321), losing `DOCK_LOWER_CRUST_LOSS` of its buoyant crust down the trench.
- **Subduct.** A patch that doesn't dock subducts whole.
- **Shared.** This lives in a new shared module, `terranes.py`, used by RFC 0006 too.

**Preservation as continent.** Buoyant crust becomes continental material in one of two
ways:
- arcs mature in place (RFC 0006);
- plateau crust converts when docked onto continental-coded receiving cells, which
  include RFC 0006's matured arcs, so the first continent doesn't have to predate the
  plateau.

On conversion, `PLATEAU_PRESERVED_FRACTION` of the plateau crust becomes continental
material, booked once as juvenile (`juvenile_plateau_m3`). The rest founders and is
removed from Hc, as for arcs, and the tracer is cleared. Plateau crust docked onto another
oceanic-coded host stays live buoyant crust there.

Plateau crust is basaltic, but mass balances count accreted plateaus as continental growth
[[Clift et al., 2009]](#ref-clift2009). Hence a default near 1, with a lower sensitivity
case.

**Tracer lifecycle and a disjoint ledger.** Registering a field as EXTENSIVE only
conserves it through refinement and coarsening remaps. Every other Hc writer has to update
the buoyant-crust tracers (`plateau_crust_m`, `arc_crust_m`, and the age moment, which
scales with them) explicitly, or they drift above Hc or vanish without a booking.
`terranes.py` holds one lifecycle for them:

| Writer | Effect on the tracers | Ledger account |
|---|---|---|
| Plume emplacement, arc growth | + the added Hc | `emplaced_m3` (by source) |
| Order-preserving thinning: divergent deformation, stretch donors, failed rifts | Scaled by the Hc ratio, like `cratons.thin_with_column` | `thinned_m3` |
| Rift melt-through reset | Cleared | `reset_m3` |
| Transfers: `crust_transfer.ConsumedColumn`, `_accrete_onto_survivors`, merge overflow (`quad_merge._place_stacking_overflow` → `quad_tectonics._carry_material`) | Carried with the Hc they ride on; `ConsumedColumn` gains the tracer volumes and age moment | none: they stay live |
| Consumption at a trench, including docking's lower-crust loss | Removed | `subducted_m3` |
| Conversion: arc maturation, plateau docking onto continental-coded cells | Cleared after converting | `converted_m3`, `foundered_m3` |
| Everything else (erosion, column caps, crust-type changes) | A phase audit, like `cratons.PhaseAudit`, holds the tracers to `0 ≤ plateau + arc ≤ E` on oceanic-coded cells and to 0 on continental-coded ones | `clipped_m3`, with the phase recorded (erosion's share is "eroded") |

- *Balance:* `emplaced = live + thinned + reset + subducted + converted + foundered +
  clipped`, plus a signed `unattributed_m3` that should stay ~0.
- `docked_m3` counts docking *events* and sits outside the balance, because docked crust
  stays live and can dock again on a later host.
- The ledger lives in `World.buoyant_crust_ledger`, keyed by source (`arc`,
  `plume_head`, `plume_tail`).

**Roots.** Ontong Java is underlain by a mantle root over 300 km thick
[[Kerr & Mahoney, 2007]](#ref-kerr2007). Plume residue incorporated into the lithosphere
is one of the reservoirs the deferred depletion RFC (from RFC 0004, now RFC 0009) must model. This RFC adds no
depletion tracer.

**Global scale check.** These are gross figures; preservation is smaller and is measured,
not assumed.
- At one head per 10 Myr with volumes log-uniform in [1, 50] M km³ (mean ≈ 12.5 M km³),
  heads emplace about 1.25 km³/yr on long-run average.
- With ~32 active tails (24 tail-only plus ~8 head tails at the starting rates) at a mean of ~0.075 km³/yr, tails emplace about 2.4 km³/yr, mostly
  on the ocean floor.
- [[Clift et al., 2009]](#ref-clift2009)'s ~1.1 km³/yr is a *preserved* (accreted) figure
  from a Cenozoic budget. It is the right comparison for `converted_m3`, not for
  emplacement.

**Runtime.**
- Each active plume tests the cells of each plate it might touch: a dot product against
  its position, after a bounding-cap test per plate. With ≤ 50 plumes and ~1.3 × 10⁵ nodes
  at density 4, that is ~6 × 10⁶ multiply-adds per step, a few ms.
- Patch grouping is a connected-components pass over consumable buoyant cells.
- The phase audits are O(cells) per phase.
- Budget: under 3% of step time.

**Save compatibility.**
- `World.mantle_plumes`, `World.plume_stream` and `World.buoyant_crust_ledger` are
  `default_factory` fields, backfilled in `persistence._backfill_added_fields`. An old save
  starts with no plumes and a fresh stream from its current time.
- The three surface fields are registered in `surface_fields.py` with defaults (0, 0, 0),
  which the registry backfills.
- `MantlePlume` gets defaults on every field, so later additions stay loadable (see
  `docs/save-compatibility.md`).

**Live knob.** `plume_multiplier` joins `world.TUNING_MULTIPLIER_FIELDS` (default 1.0,
bounded by the stream's ceiling) and scales acceptance. At 0, no plume ever spawns, which
reproduces a matched, pinned current-main run without plume code. Debug and sketch
worlds can turn plumes off that way.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `PLUME_HEAD_RATE_PER_MYR` | 0.1 | **Calibration assumption** anchored to the LIP record's average of one per ~10–20 Myr [[Bryan et al., 2010]](#ref-bryan2010), which is incomplete and not stationary | Sensitivity at 0.05 and 0.2; Cox variant |
| `PLUME_RATE_CORRELATION_MYR`, `PLUME_RATE_MODULATION_MAX` (Cox variant) | 100 Myr, 3 | Design choices | Sensitivity case only |
| `PLUME_HEAD_VOLUME_MIN_KM3`, `_MAX_KM3`, log-uniform | 1 × 10⁶, 5 × 10⁷ | **Calibration assumption**: the range spans ≥ 1 M km³ typical [[Bryan et al., 2010]](#ref-bryan2010) to Ontong Java's 44–57 M km³ [[Gladczenko et al., 1997]](#ref-gladczenko1997); the distribution's shape is not established | Log-normal and power-law sensitivity cases |
| `PLUME_HEAD_DURATION_MYR` | 3 Myr | Pulses of ~1–5 Myr [[Bryan & Ernst, 2008]](#ref-bryan2008) | Fixed. Any step length works, by interval-overlap integration |
| `PLUME_HEAD_RADIUS_KM` (continental heads and heat) | 1,100 km | Heads 2,000–2,500 km across [[Campbell & Griffiths, 1990]](#ref-campbell1990) | Fixed |
| `PLATEAU_MEAN_EXCESS_M` | 15,000 m | Derived: the model's isostasy turns it into ~2.3 km of relief, inside the 2–3 km plateaus stand above the floor [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Re-derive if the isostasy changes |
| `PLATEAU_AREA_MIN_KM2`, `_MAX_KM2` | 1 × 10⁵, 2 × 10⁶ | LIPs cover > 0.1 M km² [[Bryan et al., 2010]](#ref-bryan2010); plateaus up to 2 M km² [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Fixed |
| `PLUME_TAIL_FLUX_MIN/MAX_KM3_PER_YR` | 0.01, 0.25 | Hawaii, the strongest hotspot, peaks near 0.25 [[Van Ark & Lin, 2004]](#ref-vanark2004) | Fixed |
| `PLUME_TAIL_LIFETIME_MYR` | 80 Myr | Length of the Hawaiian–Emperor record [[Van Ark & Lin, 2004]](#ref-vanark2004) | Sweep 40–120 |
| `PLUME_TAIL_RADIUS_KM` | 150 km | A volcanic edifice and its apron; at least two cells at density 1 | Fixed |
| `PLUME_TAIL_ONLY_RATE_PER_MYR` | 0.3 | Sets ~24 active tail-only hotspots with the lifetime above (~32 with head tails); Earth has ~49, 7 deep [[Courtillot et al., 2003]](#ref-courtillot2003) | Tune so 20–50 hotspots are active |
| `PLUME_UNIFORM_SHARE` | 0.3 | Not every plume sits on an upwelling [[Burke & Torsvik, 2004]](#ref-burke2004) | Fixed |
| E_young, `DOCK_E_YOUNG_M` | 5,000 m excess | The young Caribbean plateau, 8–15 km total [[Kerr & Mahoney, 2007]](#ref-kerr2007), converted | Calibration assumption; sensitivity 2–8 km |
| E_old, `DOCK_E_OLD_M` | 15,000 m excess | Between old Caribbean-like crust and Ontong Java's ≥ 23.5 km excess [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Calibration assumption; sensitivity 12–20 km |
| τ_cool, `DOCK_COOLING_MYR` | 20 Myr | Design choice: "a few Myr" is young and 100 Myr old in [[Kerr & Mahoney, 2007]](#ref-kerr2007) | Sensitivity 10–40 Myr |
| ΔE, `DOCK_WIDTH_M` | 3,000 m | Design choice: the spread of outcomes for similar plateaus | Sensitivity 1,500–6,000 m |
| `DOCK_LOWER_CRUST_LOSS` | 0.2 | Design choice: Ontong Java lost lower crust at its edge [[Mann & Taira, 2004]](#ref-mann2004); no fraction is given | Sensitivity 0–0.4 |
| `PLATEAU_PRESERVED_FRACTION` | 0.9 | Calibration assumption: budgets count accreted plateaus as continental growth [[Clift et al., 2009]](#ref-clift2009), though they are mafic | Sensitivity at 0.5 |
| `PLUME_HEAT_MYR` | 20 Myr | The order of thermal re-equilibration of heated lithosphere; no direct constraint | Provisional; tune for plume-assisted sensitivity only after RFC 0010 Gate 3 propagation/linkage checks (see cross-RFC interface) |
| `PLUME_CRATON_WEAKENING` | 0.75 | Cancels `CRATON_RIFT_RESISTANCE` at full heat | Provisional 0.5–1.0 sweep; heated-as-weak cannot replace local damage/growth (see cross-RFC interface) |
| `PLUME_HEAD_STRESS` | 1.5 (`internal_stress` units; `STRESS_CERTAIN_RIFT` = 3.0) | Breakup often follows a continental plume head [[White & McKenzie, 1989]](#ref-white1989); no rate is given | Provisional design target: about half of continental plume heads followed by a split of their plate within 20 Myr; report matched plume-off/background incidence before stress tuning. Sensitivity at 0.75 and 3.0, plus 0 to disable |

**Breakup calibration after #358.** The #403 tuning worlds record 115, 109 and 104 splits
per Gyr on seeds 7, 11 and 2024, respectively. These whole-world counts do not measure the
conditional 20 Myr split incidence at a head's host. A split within 20 Myr of a head can
also occur under the current 2 Myr cooldown without plume stress. Before the half-of-heads target supports
stress tuning, report it beside matched plume-off/background split incidence at the same
scheduled head times and host plates, with the same step schedule and controls. The half
fraction remains a provisional design choice, and this refresh sets no new stress value.
[#407](https://github.com/adubey/mantle-bloom/issues/407) includes that baseline comparison.

**RFC 0010 target class and current evidence.** RFC 0010's first Earth-facing comparison is the
plume-associated continental-rift subset. It does not claim that every Earth rift is plume-driven;
the plume-off member of a paired test is a causal control under the same initial state and tectonic
forcing, not a model of the full non-plume rift population. The 6–8 successful continental
split-systems per 200 Myr reported by the RFC 0010 agent remain a whole-population plausibility
check, not a plume-event rate or tuning target.

RFC 0010 Gate 3 remains open ([tracking issue #487](https://github.com/adubey/mantle-bloom/issues/487)).
Its current local-damage and dike diagnostics have not demonstrated
causal tip propagation, segment linkage and a viable production split. The dike-screen figures
reported by the RFC 0010 agent (arrest bracketed between 60 and 60.5 km, event volume about
0.814 km³ against a 1–2 km³ comparison range) come from uncommitted Gate 3 diagnostics and are
not yet in-repo evidence: they identify mechanics still to solve, do not calibrate plume
weakening, and do not show that a plume caused a rift. Keep the half-of-heads
split incidence as a provisional design choice, and do not tune plume coupling to it until local
rift growth and the matched plume-off/on comparison are supported.

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time, after the isolated-head scale check passes. Record the generating
commit, density, controls and full step schedule; a regenerated seed 456862899 is a new
run, not the historical save. Report:
- LIP count and volumes per Gyr, with the oceanic/continental split;
- active hotspot count;
- every buoyant-crust ledger account, by source;
- juvenile continental by source (flood basalt, converted plateau, tail);
- craton rifted volume inside versus outside plume-heated cells;
- step time.

**Cross-RFC interface (for RFC 0010 and RFC 0012).** Other RFCs consume plume state; they do not
reimplement plume emplacement or ownership. RFC 0007 owns the plume event stream, emplacement,
the `plume_heat_myr` lifecycle, and the buoyant-crust ledger. RFC 0010 owns local cut-strength
consumption, damage calibration, propagation and final boundary formation. RFC 0012 owns the
motion-solve force interface. RFC 0010 may use plume event time, position, host class and heat as
read-only context for attribution and paired tests; plume presence alone must not create a cut or
plate split. Tectonic extension and RFC 0010's local failure/growth conditions remain necessary.

RFC 0010's Gate 2 diagnostic pairs the extensional membrane resultant with integrated capacity
converted through the same effective thickness (`sigma_ext = N_ext / H_eff`; the current state
audit uses `Hc + Hm` as shell thickness). RFC 0007's dimensionless `PLUME_HEAD_STRESS` in
`Plate.internal_stress` is only a provisional split-eligibility signal: it is not a local stress,
force or substitute for those fields. The Gate 3 plume contribution has not yet been tested against
that local failure and growth path. Therefore the RFC 0010 rows below describe available state and
future sensitivity bounds, not calibrated values. Do not tune RFC 0010's plume-coupling knobs until
the local mechanics pass their propagation/linkage checks and the matched plume-off/on comparison
is run. Any eventual tuning must keep this RFC's scale, ledger and source-attribution gates green
and report the values used.

| Consumed by | Knob or field | What the other RFC may change | Limits |
|---|---|---|---|
| RFC 0010 | `PLUME_HEAT_MYR` | Heat duration for a later plume-assisted sensitivity | Keep the 5–50 Myr range provisional; do not tune against Gate 3 until local propagation/linkage and the paired plume comparison are supported |
| RFC 0010 | `PLUME_CRATON_WEAKENING` | Strength reduction on heated cells | Keep the 0.5–1.0 sweep provisional; the product stays clamped to [0, 1], and heated-as-weak cannot replace local damage/growth |
| RFC 0010 | Plume-event schedule (read-only) | Attribute and pair local rift cases by event index, `born_years`, position, head volume, host class and heat | Event semantics are owned here; use as context, not as an automatic cut/split trigger. Branch agreement target is in the timestep acceptance |
| RFC 0010 | `PLUME_HEAD_STRESS` | Background candidate/split-eligibility signal only, if retained | Provisional `internal_stress` units; not local loading. It cannot create a cut or bypass a connected viable boundary |
| RFC 0012 | Plateau excess thickness | Derive a column ΔGPE from plume-built `plateau_crust_m`, total layer thicknesses and densities, relative to a declared compensation/reference level; project it onto the boundary | Read-only plume/material fields; the resulting load is in N/m. Thickness by itself is not a force |
| RFC 0012 | `PLUME_HEAD_STRESS` | Nothing yet | The `internal_stress` bump is not a force response and stays out of the motion solve until re-expressed as a bounded (N/m) response or retired |

**RFC 0012 boundary learned from #547 and #535.** Issue #547 is narrowed to removing the
merge path's angular-momentum rate assignment; its matched merge-only replay is the check
for that scope. The earlier combined merge-and-split experiment changed boundary
classification on the first solve after a split, and stored rates are read by boundary
classification and topology gates. Split-rate assignment therefore remains outside #547
and is deferred to RFC 0010's split work; this is not evidence that every topology rate
assignment can be removed together. Keep plume events and `PLUME_HEAD_STRESS` on the
RFC 0010 eligibility/trigger side of the interface, not as motion loads. RFC 0012's G3
transferability screen also used prescribed synthetic rift-strength profiles: it did not
derive that strength from `World` state or validate a production plume-load law. Treat the
column ΔGPE coupling above as a candidate physical interface until its World-derived inputs
and plume response are separately checked.

## Detailed design (provisional touchpoints)

**Checked against current main, `2562f63` (2026-10-14).** `mantle.ConvectionCenter` is
still tangential forcing without a plume source; `world.step_world_progress` has the
existing topology, erosion and eruption phases, and no plume pass. `plumes.py`,
`terranes.py`, the three world fields, three surface fields and the knob below are new.
The eruption writer remains `volcanism.apply_volcanic_activity`, so its independent Hc
and juvenile bookings still require the one-owner contract. Merge overflow is handled
by `quad_merge._place_stacking_overflow`, which calls `quad_tectonics._carry_material`;
those existing writers, not a new merge path, must carry the tracers. #358's 2 Myr split
cooldown is now part of the breakup-stress tuning baseline. The existing
`_erupt_melted_nodes` elevation-based recoding and juvenile booking in both column and
stretch-rift callers are additional source-attribution touchpoints; they are not plume
emplacement or plateau docking.

- `backend/app/plumes.py` (new):
  - `MantlePlume`;
  - the index-keyed event stream, with thinning and the Cox modulation;
  - `apply(world, years)` (head, tail and heat, per plate, by interval overlap);
  - the footprint kernels;
  - `diagnostics(world)`;
  - constants with citations.
- `backend/app/terranes.py` (new, shared with RFC 0006): patch grouping,
  `docking_probability(excess, age)`, the per-encounter draw, partial loss, conversion,
  the tracer lifecycle helpers (`thin_tracers`, `clear_tracers`, `carry_tracers`), the
  phase audit for tracers, and the buoyant-crust ledger.
- `backend/app/world.py`: the `mantle_plumes`, `plume_stream` and `buoyant_crust_ledger`
  fields and the `plume_multiplier` knob; call `plumes.apply(world, years)` in
  `step_world_progress` after topology, gap-fill and magma-transport have settled
  (inside `if world.simulate_plate_movement`, after the magma-transport block) and
  before `gather_node_positions` / `_advance_fluid_dynamics`, so fresh plume crust
  exists before climate, erosion and volcanism read it in the same step. Bracket the
  existing phases with the tracer audit, as the craton audit is.
- `backend/app/persistence.py`: backfill the three world fields.
- `backend/app/surface_fields.py`: register `plateau_crust_m`, `buoyant_crust_time_m_yr`
  and `plume_heat_myr`.
- `backend/app/quad_tectonics.py`
  - `_retreat`: patches that `terranes` docks join the `terrane` donors, minus their
    lower-crust loss; the rest subduct, booked.
  - `_accrete_onto_survivors` and the conversion on continental-coded receivers.
  - `_allocate_stretch`: thin tracers with donors; plume weakening on craton strength.
- `backend/app/crust_transfer.py`: `ConsumedColumn` gains the tracer volumes and age
  moment, and the placement carries them.
- `backend/app/quad_merge.py`: carry the tracers with merge overflow.
- `backend/app/lithosphere_plate.py`, `deform_columns`: thin tracers under divergence;
  clear them on melt-through; apply plume weakening to `CRATON_RIFT_RESISTANCE`.
- `backend/app/sparse_quad_patch.py`, `apply_failed_rift`: thin tracers; apply the same
  weakening, with `cratons.rift_resistance`'s clamp.
- `backend/app/merge_split.py`: RFC 0010's cut strength counts heated cells as weak; add
  the head stress to `internal_stress`.
- `backend/app/hydroclimate/volcanism.py`: unchanged. Its untouched eruption path is part
  of the volume-ownership contract, covered by a test.
- `backend/app/continental_ledger.py`: sub-accounts `juvenile_flood_basalt_m3`,
  `juvenile_plateau_m3` and `juvenile_hotspot_m3` under `juvenile_additions_m3`.
- `backend/app/elevation_lines.py`: correct the "no hotspot spawning yet" comment.
- `backend/app/render_image.py`: a `plumes` debug view (plume positions, head footprints,
  tail age, plateau crust thickness), added to `DEBUG_VIEWS`; `frontend/src/api.ts` view
  union and `legendData.ts` legend.
- `backend/app/stats.py`: active plumes, LIP events, buoyant-crust ledger and areas.
- `backend/app/main.py`: `GET /world/plumes` listing live plumes, optional; the debug
  view alone meets the criteria.
- Tests:
  - `backend/unit_tests/test_plumes.py` (new):
    - **isolated-head scale check:** a median head on a reference oceanic plate yields
      2.0–2.6 km of peak relief and the configured area; a minimum and a maximum head stay
      within 1.2–4.5 km;
    - the event stream yields identical events (indices, times, positions, volumes) for
      ten 1 Myr steps and one 10 Myr step, and for any random partition;
    - over random step partitions, including steps longer than the head, total head and
      tail delivery equals the scheduled volume and never exceeds it;
    - kept events at multiplier 0.5 are a subset of those at 1.0;
    - a plate moving over a fixed tail leaves a line of crust maxima parallel to its
      motion, getting younger toward the plume;
    - **one volume owner:** with plumes and the eruption lifecycle both on, no plume cell
      becomes `is_volcano`, and the step's total Hc change equals plume delivery plus
      non-plume eruptions, exactly;
    - weakening halves effective craton resistance at 0.5, clamped, and a step
      extending beyond the remaining heat scales that step's weakening by the
      heat/step overlap;
    - `plume_multiplier = 0` changes nothing.
  - `backend/unit_tests/test_terranes.py` (new):
    - docking probability is in [0, 1], rises with excess and falls with age, and matches
      the example values above;
    - a patch's draw is stable across steps within one encounter;
    - a docked patch loses exactly its lower-crust share;
    - **tracer lifecycle:** each writer in the table, applied to a synthetic plate
      (erosion, divergent thinning, a failed rift, melt-through, a column cap, a
      `crust_transfer` consumption, a merge overflow, two successive dockings on oceanic
      hosts, then conversion on a continental one), keeps `0 ≤ tracers ≤ E` and books
      exactly what it removes, and the buoyant-crust ledger closes at every stage;
    - `docked_m3` counts both dockings while the balance stays exact.
  - `test_persistence.py`: an old save loads with no plumes, an empty ledger and default
    fields.
  - `test_continental_ledger.py`: the juvenile sub-accounts sum correctly.
- Docs:
  - `docs/simulation-model.md`: a new "Mantle plumes (`plumes.py`)" section; update the
    Volcanism, Cratons (Resistance) and Known simplifications sections.
  - `docs/architecture.md` (module map, including `terranes.py`), `docs/api-reference.md`
    (view, endpoint, stats keys) and `docs/save-compatibility.md` (new fields).

## Assumptions and uncertainty

- **The arrival process is assumed.** The LIP record is incomplete and clustered; the
  stationary default is a simplification, and the rate and volume distribution are
  calibration assumptions. The sensitivity variants show how much the results depend on
  them.
- **Fixed plumes.** Real hotspots drift slowly relative to each other
  [[Courtillot et al., 2003]](#ref-courtillot2003). A fixed mantle frame is a
  simplification, matching the model's fixed convection centres.
- **Placement is an analogy.** Biasing plumes toward upwellings stands in for the
  deep-mantle structure [[Burke & Torsvik, 2004]](#ref-burke2004) that the model doesn't
  have.
- **The plateau footprint is derived from the isostasy, not observed.** Tying area to
  volume through one mean excess thickness reproduces typical relief, but real plateaus
  vary in thickness (Caribbean versus Ontong Java). Very large heads overshoot.
- **Docking is uncertain.** The age-dependent threshold, its width and the lower-crust
  loss are shaped by a few well-studied plateaus, and the conversion from total to excess
  thickness is approximate. The probabilistic form admits that.
- **Plateau crust is mafic.** Counting most of it as continental once docked follows mass
  balances, not composition. `PLATEAU_PRESERVED_FRACTION` exposes the choice.
- **Weakening has no direct constraint.** `PLUME_HEAT_MYR` and the strength reduction are
  tuned to an outcome (craton rifting clustering at plumes), not measured.

## Acceptance criteria

- [ ] **Source-attribution gate ([#407](https://github.com/adubey/mantle-bloom/issues/407),
      before tuning):** record the current-main
      elevation-driven reset route on the stripped baseline shared with RFC 0006 and on
      an exposed oceanic-column example, covering both the column pass
      (`lithosphere_plate.deform_columns`) and the stretch-rift donor and rifted-cell
      paths (`quad_tectonics`). With plume code, distinguish direct plume
      juvenile volume, docked plateau conversion, and any juvenile volume the existing
      `_erupt_melted_nodes` path makes after plume uplift (donor and rifted-cell
      bookings separately). Report the latter separately;
      aggregate juvenile growth does not establish plume emplacement or docking.
- [ ] **Isolated-head scale check (before tuning):** the unit test above passes. A median
      head on reference oceanic lithosphere gives 2.0–2.6 km of peak relief.
- [ ] Over 1 Gyr on each of seeds 7, 11 and 2024: plume-head counts consistent with the
      configured rate (Poisson 95% interval), a volume distribution matching the
      configured one, and 20–50 active hotspots at every 100 Myr checkpoint.
- [ ] Timestep robustness:
  - [ ] the event stream is identical across step partitions (unit test);
  - [ ] for a shared predetermined schedule, total delivered volume is identical within
        1e-9 relative for 1 Myr and 10 Myr steps;
  - [ ] on a pinned-motion test the oceanic/continental branch agrees on over 95% of
        events between 1 Myr and 10 Myr partitions;
  - [ ] on full replays, totals agree within 2% (positions differ).
- [ ] Hotspot tracks: on a pinned-motion test world, one tail's crust maxima lie within 20°
      of the plate's motion direction at the plume and get younger toward it.
- [ ] The buoyant-crust ledger closes on every acceptance run, with `unattributed_m3`
      below 1e-6 of emplaced volume, and the tracer invariants hold at every step.
- [ ] Gross emplacement, docking events and conversion are reported separately. Long-run
      `converted_m3` from plateaus is 0.3–1.5 km³/yr on the tuning seeds at the starting
      knobs. That is a sanity range around ~1.1 km³/yr of Cenozoic plateau accretion
      [[Clift et al., 2009]](#ref-clift2009), not a fit.
- [ ] The sensitivity runs (rate, clustering, volume distribution, docking knobs) are
      reported, with the spread of converted continental material each produces.
- [ ] Any RFC 0010/0012 tuning of the cross-RFC coupling knobs is reported with the
      values used, and this RFC's scale, ledger and source-attribution gates are
      re-run green at those values. RFC 0010 tuning also requires its local propagation/linkage
      checks and matched plume-off/on comparison to be supported.
- [ ] Side effects stay bounded: report plate counts and mean overlap fraction at every
      100 Myr checkpoint beside the matched plume-off baseline; flag over 20% drift
      for review before tuning further.
- [ ] With RFC 0010 in place, report craton volume lost to rifting inside and outside plume-heated
      cells, split between plume-associated and other rift episodes, with a matched plume-off/on
      comparison for the plume-associated subset. Do not use a single whole-world 50% localization
      threshold: non-plume rifts are outside the initial plume-assisted target class.
- [ ] Before tuning `PLUME_HEAD_STRESS`, report the fraction of continental heads followed
      by a split within 20 Myr beside matched plume-off/background incidence under the
      physical 2 Myr cooldown. The half-of-heads target remains a provisional design
      choice, not evidence that a later split was caused by a plume.
- [ ] Step time changes by less than 3% at density 4; save size by less than 3%.
- [ ] `plume_multiplier = 0` gives step results identical to a matched, pinned
      current-main run without plume code on a 10-step replay; record both commits,
      density, controls and step schedule.
- [ ] Ledgers close; old saves load; affected tests, ruff and `bin/quality_check.sh` are
      green.

## Falsification criteria

- **Plateaus don't survive.** If, across the docking sensitivity range, under 10% of
  emplaced oceanic plateau crust ever docks, plateaus in the model meet trenches too old
  or too thin. They can't be a crustal source here; drop the docking part and keep plumes
  only for tracks and weakening.
- **The source is all in the knobs.** If converted continental material varies by more
  than a factor of five across the rate and volume sensitivity cases, the model's crust
  budget would be set by unconstrained assumptions. Don't merge plateau preservation until
  a better constraint exists.
- **Weakening doesn't localise.** If, with RFC 0010, craton rifting doesn't concentrate in
  plume-heated cells at any `PLUME_HEAT_MYR` from 5 to 50 Myr, the weakening mechanism
  isn't what drives craton rifting in the model. Remove it.

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
