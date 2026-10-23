# RFC 0015: Frozen water in physical units: glacier mass balance, freezing lakes and rivers, ocean ice and polar snow

- **Status:** Draft
- **Date:** 2026-10-23
- **Related:** #368 (evidence; PR #536), #353 and #359 (ice-sheet model, decision B), #355 (water units), #503 and #504 (lake freezing at `WATER_TO_ICE`), #297 (balanced lakes), RFC [0014](0014-lake-water-volume-and-freezing.md) (lake water volume, seasonal lake ice), RFC [0012](0012-physics-first-simulation.md) (physics-first simulation); design principles 3, 5, 6 and 8

## Problem

Water changes phase in five places in Mantle Bloom: snow falling on land, lakes freezing,
rivers freezing, ice reaching the ocean, and ice and snow drawn on the maps. Each place uses
its own units and rules, and they don't agree.

- **Liquid water is in physical units.** Runoff and lakes balance in m³/yr (#297, #355).
- **Glacier mass balance is in model units.** Snowfall builds ice at
  `hydrology.GLACIER_ACCUMULATION_RATE` = 0.02 m per Myr per mm/yr (decision B in #353).
  That is 54,526× slower than the physical exchange of 1000/917 mm of ice per mm of water.
  Geothermal basal melt (`ice_sheet.GEOTHERMAL_BASAL_MELT_M_PER_MYR`) and the exchange for
  freezing rivers and lakes (`hydrology.WATER_TO_ICE`) are defined through the same factor.
- **Stored ice is in physical units.** `glacier_depth` is metres of ice: eustasy debits it at
  0.917 m of water per metre (`eustasy._ICE_WATER_EQUIVALENT`), and lithosphere loads it.

So ice is real once it exists, but it forms and melts tens of thousands of times too slowly,
and the conversions at its borders with water, the ocean and the maps are inconsistent.

**Principles violated:**
- #5, consistent units: one process, two unit systems.
- #6, account for material through every process: snowfall water vanishes, and clipped
  river ice changes reservoir unbooked.
- #3, fast processes in aggregate over long steps: annual-mean switches for precipitation
  phase and freezing.
- #8, each process follows its own physics: a floating sea-ice cap runs the glacier mass
  balance, up to the 5 km glacier depth guard.

### Evidence

All replay numbers come from #368 ([`analysis/issue368/`](../analysis/issue368/README.md)),
on `origin/main` `30766af`, with two saves:
- seed 997271774 at 87.4 Myr (SHA-256 `08ad9461…ed762b`);
- seed 649658380 at 6.4 Myr (`ebbfb5da…73f7b3`).

Both are local and **not available for download**; the numbers are quoted so the argument
stands without them. Each arm replayed 10 Myr in 0.25 Myr steps. `kphys` scales snowfall
accumulation, basal melt and `WATER_TO_ICE` by 54,526 and surface melt by the same factor.
`phys` uses a physical degree-day melt (139,408×) and physical abrasion. Neither save has
ice-age forcing.

The map and sea-level numbers come from a new read-only probe,
[`bin/debug/measure_frozen_water.py`](../bin/debug/measure_frozen_water.py), on the same
two saves at their saved state ([`analysis/rfc0015/`](../analysis/rfc0015/README.md)).

1. **At physical rates the ice doesn't settle, and the result depends on the step.**
   - In the current model, ice volume is nearly frozen in place: 5.7 → 6.0 M km³ over 10 Myr
     on seed 997271774. It is still the save's inherited ice, far from its own balance.
   - Scaling accumulation and melt by 100 brings it to a new steady level (about 18–19 M km³)
     within 4 Myr.
   - At physical rates, volume swings 10–15× between 1 Myr checkpoints (`kphys`: 748, 115,
     177, 1,016, 349 … M km³), and sea level swings by up to 3 km. On seed 649658380 the swing
     repeats about every 3 Myr.
   - `kphys` ends at 125 M km³ with 0.25 Myr steps and 927 M km³ with 1 Myr steps.
   - A physical sheet responds in about 10 kyr (a few km of ice at ~300 m/kyr of
     accumulation), much shorter than the step. Sea level, the bed's isostatic deflection
     (`ice_load_deflection_m`) and frozen river ice reach the ice solve only a step late. The
     swings look like the overshoot of those lagged couplings, but #368 didn't isolate which
     coupling drives them (gate G0).
2. **Erosion's per-step limits do the abrasion law's work.**
   - Raw abrasion with today's coefficient rose about 1,500× at `kphys`, while applied
     erosion rose 9×. The neighbor-drop and Hc-headroom limits clipped 99% of it.
   - With a physical coefficient (`phys`), 93–96% is clipped.
   - With today's calibration divided by 54,526 (`phys-neutral`), glacial erosion falls
     100–200×, because sliding grew only 120–260×, not 54,526×.
3. **Freezing rivers fill the depth guard.** At physical `WATER_TO_ICE`, river water settling
   at a frozen node becomes kilometres of ice in one step. The guard clips that back to the
   ocean at 17–57% of meltwater's volume (`IceLedger.frozen_river_capped_m3`). Ice shed from
   sheets stays below 0.01%.
4. **Snowfall water vanishes today.**
   - On a node whose annual mean is below 0 °C, `compute_hydrology` removes all
     precipitation from physical runoff (`liquid_precip`), but the ice gains only 1/54,526 of
     it.
   - In the current model that is about 7,100 km³/yr on seed 997271774 (~225,000 m³/s, a
     quarter of discharge into the ocean) and about 1,200 km³/yr on seed 649658380. The water
     leaves no reservoir; it is simply not routed.
5. **Precipitation phase is an annual-mean switch.** On land below 0 °C annual mean, all
   precipitation is snow. Under the model's own sinusoidal annual cycle, 28% (seed 997271774)
   and 44% (seed 649658380) of that precipitation falls in months above 0 °C.
6. **Lakes and rivers freeze at the glacier's exchange.**
   - A standing lake below 0 °C annual mean freezes into glacier ice at `WATER_TO_ICE`
     (`frozen_from_lake`, since #504); a settling river does the same a step later
     (`frozen_river_ice`).
   - Until #504, lakes froze 1:1 and fed the ice about 35× faster than all snowfall (#503).
     The 1/54,526 factor is a patch at the same broken boundary.
   - RFC 0014 (D2) removes both conversions for liquid water, and leaves the hand-off of
     frozen-to-bed lake ice to the glacier "until glacier mass balance is physical (#368)".
7. **Floating ice lowers sea level.**
   - `eustasy.trapped_water_volume_m3` debits every metre of `glacier_depth`, including the
     ice-age sea-ice cap over the ocean and marine ice grounded below sea level.
   - Floating ice already displaces its own weight of water and shouldn't move sea level. Only
     the part of marine ice above flotation should
     [[Goelzer et al., 2020]](#ref-goelzer2020).
   - On both saves today, floating and marine ice is under 0.1% of ice volume (0.4 and 71 km³
     water equivalent floating; 3,551 and 2,803 km³ marine). The ice-age cap, though, runs the
     glacier mass balance per node up to the 5 km guard (`ice_sheet.step_ice`), while real
     multiyear sea ice is about 3 m thick [[Maykut & Untersteiner, 1971]](#ref-maykut1971).
8. **The maps show the bed, not the ice.**
   - The Elevation view colours each cell by **bed** elevation and paints white where ice is
     at least `GLACIER_VISIBLE_DEPTH_M` (10 m) (`render_image.elevation_colors`, the glacier
     overlay).
   - Climate's air temperature is the bed's (the ice solve lapses it itself:
     `ice_sheet.SurfaceClimate`), and the Biome view classifies Köppen zones from that
     temperature. Land ice only affects ocean cells (`biomes._smooth_biome_classify`'s
     `has_sea_ice`).
   - On seed 649658380, 33% of land under visible ice is classified EF (Ice Cap) and 60% ET
     (Tundra). For 37% of it, the warmest month is above 0 °C at the bed but would be below
     0 °C at the ice surface (mean ice 853 m, 5.5 °C colder).
   - On seed 997271774, 93% of iced land is EF, but 39% of EF land carries no ice.

### What the evidence does not show

- Which lagged coupling drives the instability at physical rates. #368 changed every rate
  together and recorded only 1 Myr checkpoints.
- Whether the maps' mismatch survives a physical mass balance. Ice in the current model
  barely responds to climate, so EF land without ice may simply be ice that hasn't grown yet.
- The size of the eustatic error at physical rates or in an ice age. #368 didn't split
  floating, marine and land ice.
- Anything beyond two saves at one node density, without ice-age forcing.

## Proposed solutions

### What happens in nature

- **Snow builds ice by mass, and ice melts by energy.** A water-equivalent metre of snowfall
  becomes 1000/917 m of glacier ice once compacted. Temperature-index melt models relate melt
  to positive degree-days. At glacier equilibrium lines, degree-day factors for melting snow
  average 3.5 ± 1.4 and 4.6 ± 1.4 mm d⁻¹ K⁻¹ for low and high accumulation estimates on 66
  glaciers, derived with a sinusoidal annual temperature cycle like the model's own
  [[Braithwaite, 2008]](#ref-braithwaite2008).
- **Precipitation phase depends on the temperature at the time it falls.** Across the
  Northern Hemisphere, the air temperature at which rain and snow fall equally often averages
  1.0 °C, ranging from −0.4 to 2.4 °C at 95% of stations
  [[Jennings et al., 2018]](#ref-jennings2018). A year with a sub-zero mean but thawing months
  gets part of its precipitation as rain.
- **Ice sheets respond in millennia, and ice-age models couple them asynchronously.**
  Asynchronous schemes run ice and atmospheric forcing together for short synchronous periods,
  then run the ice alone with forcing extrapolated from them
  [[Pollard et al., 1990]](#ref-pollard1990). The lesson for a model whose step is far longer
  than the ice's response: the ice must see the feedbacks it changes within the step, not a
  step later.
- **Glacial erosion scales with sliding.** The dimensionless erosion rate of Variegated
  Glacier's surge, metres eroded per metre slid, was of order 1.0 × 10⁻⁴
  [[Humphrey & Raymond, 1994]](#ref-humphrey1994).
  - A global compilation of 38 glaciers gives the same mean coefficient, K_G = 1 × 10⁻⁴ for
    erosion in m/yr, with a sub-linear exponent of 0.65 (0.69 including Meserve Glacier)
    [[Cook et al., 2020]](#ref-cook2020).
  - One Alpine glacier gives an exponent of 2 [[Herman et al., 2015]](#ref-herman2015).
- **Rivers don't freeze into glaciers.**
  - River ice covers many cold-region rivers for part of the year; its hydrologic effects are
    ice jams, floods and low winter flows [[Prowse & Beltaos, 2002]](#ref-prowse2002).
  - Where winter water is expelled onto frozen ground or ice, it builds aufeis (icings) in
    layers [[Ensom et al., 2020]](#ref-ensom2020).
  - Both are seasonal stores. Liquid water reaching a cold sink pools there (RFC 0014 D2b).
- **Floating ice is thin and doesn't move sea level.**
  - Sea-ice growth is set by conduction through the ice [[Leppäranta, 1993]](#ref-leppaesranta1993).
  - A thermodynamic model of central-Arctic sea ice gives an equilibrium average thickness of
    288 cm, in close agreement with observations
    [[Maykut & Untersteiner, 1971]](#ref-maykut1971).
  - For sea level, ice counts only by its volume above flotation, and that convention needs
    care when the bed or sea level move for other reasons
    [[Goelzer et al., 2020]](#ref-goelzer2020).
  - Scale check: at the Last Glacial Maximum, grounded ice was about 52 × 10⁶ km³ greater than
    today, and sea level about 134 m lower [[Lambeck et al., 2014]](#ref-lambeck2014).
- **Ice cap is a climate class and a surface state.** Köppen's polar group has a warmest
  month below 10 °C: ET (tundra) above 0 °C, EF (frost) at or below 0 °C
  [[Peel et al., 2007]](#ref-peel2007). An ice sheet's surface is EF partly because the ice
  lifts it kilometres into colder air, while the bed beneath may be warmer.

### Alternatives

**A. Change the constants.** Set accumulation, melt, basal melt and `WATER_TO_ICE` to
physical values.
- *For:* one line each.
- *Against:* #368 shows the result swings 10–15× and depends 4–7× on step size, the erosion
  caps clip 93–99%, and the depth guard clips river ice by the cubic kilometre. Rejected.

**B. Keep model units and fix the borders.** Book snowfall water consistently at 1/54,526,
rescale meltwater, and leave the ice slow.
- *For:* small and stable.
- *Against:* ice doesn't reach its balance in hundreds of Myr, so its extent records the
  initial state, not the climate. The model's other water is physical, so every new coupling
  needs a conversion. Principle #5 stays violated.

**C. Physical mass balance, solved to balance within the step where the ice responds fast,
with the fast feedbacks inside the solve; physical borders with lakes, rivers, the ocean
and the maps (recommended).**
- Accumulation, melt and basal melt in physical units, with precipitation phase from the
  seasonal cycle.
- Each ice body's state solved to step's end with sea level, bed deflection and ice-surface
  climate updated inside the solve. Ice bodies slower than the step are integrated in
  substeps.
- One owner for each frozen volume, with conversions by mass at every border.
- Maps that show ice as a surface.

**D. Full sub-stepping of the hydroclimate.** Run climate, hydrology and ice at a short
sub-step inside every tectonic step.
- *For:* closest to the physics.
- *Against:* the hydroclimate costs seconds per step today; ×25–100 sub-steps is out of
  budget (principle #2). Kept as the fallback per ice body, not globally.

**Recommendation: C.** It follows RFC 0014's pattern for lakes: an equilibrium endpoint
where the response time is short against the step, with integration for bodies where it
isn't, and an auditable ledger. G0 tests which couplings must move inside the solve before
the design is fixed.

## Design

### D1. Physical mass balance

All rates are metres of ice per year; `ice_sheet.SurfaceClimate.rates` integrates them over
the step as today.

- **Accumulation** = snowfall in water equivalent × ρ_w / ρ_i, with ρ_w / ρ_i = 1000/917.
  This replaces `GLACIER_ACCUMULATION_RATE`'s model unit and removes `MODEL_UNIT_SCALE`.
- **Snow fraction from the seasonal cycle.**
  - Precipitation is assumed spread evenly over the year (a stated proxy; see Assumptions).
    Its snow fraction is the share of the year below the rain–snow threshold `T_s`:
    `f_snow = 1 − arccos(−(T − T_s) / A) / π`, clamped to [0, 1].
  - `T` and `A` are the surface mean and half-amplitude: the ice surface on ice, the ground
    elsewhere.
  - This is the closed form of the cycle `hydrology.seasonal_positive_degrees` already uses.
  - It replaces the annual-mean switch (`is_frozen` → `liquid_precip = 0`). The rain fraction
    `1 − f_snow` becomes runoff everywhere, ice included.
- **Snow on bare ground.** The share of snow that the year's positive degree-days melt on
  ice-free ground returns as meltwater (`_melt_bare_snow`, as today, in physical units). The
  remainder starts or feeds ice.
- **Surface melt** = DDF × positive degree-days × ρ_w / ρ_i.
  - Positive degree-days = 365.25 × `seasonal_positive_degrees` on the lapse-corrected surface.
  - `DDF_snow` = 4 mm w.e. d⁻¹ K⁻¹ (Braithwaite's range), applied to the year's snow.
  - `DDF_ice` is applied to bare ice once that snow is gone. Its value must be sourced and
    verified before implementation (G2). #368 used 7 mm w.e. d⁻¹ K⁻¹ unverified.
- **Geothermal basal melt** = G / (ρ_i L_f) on a thawed bed: about 6.2 mm/yr of ice at
  G = 0.06 W m⁻², with `MODEL_UNIT_SCALE` dropped.
- **Bounds:** every rate ≥ 0; `f_snow` ∈ [0, 1]; melt in a step can't exceed the opening ice
  plus that step's accumulation (the solve's existing volume contract).
- **Meltwater** is routed by mass: ice volume × ρ_i / ρ_w (the TODO at
  `hydrology.WATER_TO_ICE`). Calved ice leaves for the ocean unrouted, as today.

### D2. Solving the ice to the end of the step

**The premise.** A body's volume response time is `τ = H / |b|`, with `H` its mean thickness
and `b` its mean absolute mass-balance rate. Using #368's numbers, a 3 km sheet under
0.3 m/yr of accumulation gives about 10 kyr, short against a 0.1–1 Myr step. Where τ ≪ Δt,
the ice at the end of the step is the balanced sheet for that step's climate, as with
balanced lakes in RFC 0014.

**What must be inside the solve.** The ice changes sea level, the bed it rests on, and its
own surface climate within the step. Each of these is a feedback that must be evaluated
with the ice it responds to, not last step's ice:
- **Ice-surface climate:** already inside. `rates(thickness)` lapses temperature by the ice
  thickness on every root-find evaluation.
- **Sea level:** today lagged one step. The solve would iterate a 1-D fixed point: ice
  volume above flotation → `eustasy` sea level → flotation and calving → ice volume. This is
  monotonic and cheap (one bisection per iteration).
- **Isostatic deflection:** today lagged one step (`ice_load_deflection_m`). This design assumes
  the mantle's relaxation under an ice load is also short against the step, so the
  deflection is near equilibrium by step's end (unsourced; G0 tests whether the lag
  matters). It would be solved with the ice load (`lithosphere.grounded_ice_load_kg_m2`) in the
  same fixed point.
- **Frozen river and lake ice:** removed as a glacier source (D4, D5); this lag disappears.

G0 identifies which lag drives #368's instability before this list is fixed.

**Bodies slower than the step.** For bodies with τ > 0.1 Δt (thin, dry or very large),
integrate `dV/dt` in substeps of at most τ / 4 on the step's frozen climate. This mirrors
RFC 0014's per-body fallback. G1 reports their share and cost.

**Timestep robustness:** rates are annual, and the endpoint state doesn't depend on Δt for
fast bodies. Substepped bodies integrate along their path, not at their end state.

### D3. Glacial erosion

- **Law:** `E = K_G · u_b^l`, with `u_b` the physical sliding speed (m/yr) and E in m/yr.
  Starting values: K_G = 1 × 10⁻⁴, l = 1
  [[Humphrey & Raymond, 1994]](#ref-humphrey1994) [[Cook et al., 2020]](#ref-cook2020). This
  replaces `GLACIER_ABRASION_COEFFICIENT` = 430 per model-unit velocity.
- **Integrated along the step.** A sliding rate times a long step overshoots the neighbor-drop
  limit, which then does the law's work (Evidence item 2). Glacial erosion deepens the bed
  only until the ice's surface slope, and so its sliding, changes. Within a step, applied
  glacial erosion is bounded by the relief the ice can remove, and that bound should be the
  physical one (valley-floor grade under the ice surface), not a clip. G3 measures how much
  applied erosion still depends on step size and on the clip.
- **Calibration.** The current coefficient was set to match the old depth × slope law
  (#353). Under physical sliding, `K_G` is literature-constrained, and the global glacial
  erosion volume is reported against `origin/main` as a sanity check, not a target.

### D4. Lakes that freeze (with RFC 0014)

RFC 0014 owns lake water: its volume, seasonal ice cover and the bounded lake ice cap
`lake_ice_m` (D2). This RFC owns the border with the glacier.
- **No lake water becomes glacier ice by exchange factor.** `frozen_from_lake` at
  `WATER_TO_ICE` is removed, as 0014 D2(b) proposes.
- **A cap frozen to the bed hands over by mass.** When a lake cell's `lake_ice_m` (water
  equivalent) reaches its water depth, the cell has no liquid left. Its ice then joins
  `glacier_depth` as `lake_ice_m × ρ_w / ρ_i`, and `lake_ice_m` is zeroed in the same step.
  One owner at a time.
- **Glacier meltwater into lakes** is booked by mass (D1), so lake budgets close with it.
- **Ice over a lake.** A glacier advancing over a lake cell displaces its water, which spills
  through the lake solve. Ice doesn't float on lakes in this design (out of scope; see
  Assumptions).

### D5. Rivers that freeze

- **No river water becomes glacier ice.** `frozen_river_ice` and `WATER_TO_ICE` are retired.
  Liquid water reaching a cold sink pools there (RFC 0014 D2b), and the lake's cap holds it
  as ice by D4's rules. This removes the depth guard's 17–57% of meltwater in #368's
  physical runs.
- **River ice cover** is seasonal [[Prowse & Beltaos, 2002]](#ref-prowse2002). Over a step
  of a year or more it changes when water arrives, not how much, so it isn't represented.
- **Aufeis** [[Ensom et al., 2020]](#ref-ensom2020) is a seasonal store that melts out each
  summer where the year has positive degree-days. Where a cell has none, the water already
  pools and freezes through D4. No separate store is added.
- **Snowfall leak:** with D1, precipitation leaves runoff only as snow, and snow becomes ice
  or meltwater by mass. Diagnostic: liquid water lost outside the ledger = 0.

### D6. The ocean

- **Sea level counts ice above flotation.** `trapped_water_volume_m3` debits only grounded
  ice's volume above flotation, `max(0, H − (ρ_w / ρ_i) · max(0, sea level − bed))`, in water
  equivalent. Floating ice is excluded. The ocean area lost to grounded marine ice is
  already handled by the connected-ocean mask.
  - The generalised approach of [[Goelzer et al., 2020]](#ref-goelzer2020) handles bed and
    sea-level changes not caused by the ice. It would need the bed change booked separately;
    adopt it only if G4 shows the simple form drifts.
- **Sea ice is a thin, seasonal cover, not a glacier body.**
  - Each ocean cell carries a derived sea-ice fraction of the year,
    `f_ice = 1 − arccos(−(T_sst − T_f) / A_sst) / π`, with `T_f` = −1.8 °C (seawater
    freezing). Its thickness approaches a conduction-limited equilibrium of a few metres
    [[Leppäranta, 1993]](#ref-leppaesranta1993)
    [[Maykut & Untersteiner, 1971]](#ref-maykut1971), capped at `SEA_ICE_MAX_THICKNESS_M`
    (new, 3 m).
  - The ice-age sea-ice cap (`allow_ocean_ice`, `SEA_ICE_FORMATION_TEMP_C`,
    `ICE_AGE_SEA_ICE_MIN_INTENSITY`) no longer runs the glacier mass balance or accumulates
    up to the 5 km guard.
  - Sea ice is floating, so it doesn't change sea level and isn't persisted. It feeds the
    Biome view's "Polar Sea Ice" and the climate's ocean surface.
- **Ice shelves are out of scope.** Ice that floats calves, as today. Calved ice returns to
  the ocean by mass.

### D7. Polar snow and ice on the maps

- **Elevation.** The Elevation view and the inspector show the **surface** elevation,
  bed + ice, on ice cells: its colour, hillshade and value. The bed stays available in the
  inspector. A plastic ice sheet's smooth dome then reads as a dome, not white paint over
  buried relief.
- **Temperature.** The Temperature view and the climate fields used by maps read at the
  surface, lapse-corrected by the ice: the same correction the ice solve applies, from one
  shared function. The mass balance keeps using its own per-evaluation thickness, so nothing
  is lapsed twice.
- **Biome.** A land cell with perennial ice (≥ `GLACIER_VISIBLE_DEPTH_M`) is shown as Ice Cap
  whatever its Köppen class. Its climate class at the surface stays in the inspector. EF
  without ice is still possible (cold and dry) and is reported, not overridden.
- **Seasonal snow** isn't drawn: it melts within the year by D1's degree-day term. A
  snow-cover-duration layer could later use `f_snow` and the melt season (out of scope).

### Diagnosis gates

**Seed sets** (recorded on the tracking issue before G0, disjoint):
- a **diagnosis set** of ≥ 3 saves, including seed 997271774 at 87.4 Myr and seed 649658380 at
  6.4 Myr, for G0–G4;
- a **held-out acceptance set** of ≥ 3 other saves, at least one with ice ages on, never used
  for gates or knob choices.

- **G0: which lag drives the instability?** At D1's physical rates, on the diagnosis set,
  over 10 Myr at 0.25 Myr and 1 Myr steps, with checkpoints every step:
  - (a) today's lags;
  - (b) sea level frozen at its opening value;
  - (c) bed deflection frozen;
  - (d) frozen river and lake ice disabled;
  - (e) all three.

  Measure the step-to-step swing in ice volume and the 0.25 vs 1 Myr gap.
  *Outcome:* each knockout that removes the swing and closes the gap moves that coupling into
  D2's solve. If (e) still swings, the instability is in the ice solve itself: D2 falls back
  to substepping all bodies (alternative D, per body), and the RFC is revised.
- **G1: where does the endpoint premise hold?** On the diagnosis set, per body: `τ = H / |b|`
  at step's end, and for bodies above 10⁴ km² a direct check that substepping the step's
  frozen climate reaches the endpoint volume within 1%.
  *Pass:* ≥ 95% of ice volume has τ ≤ 0.1 Δt at 0.25 Myr steps. Bodies that fail are
  substepped; report their share and cost.
- **G2: melt factors and snow phase.** Source and verify `DDF_ice`. Run sensitivity of ice
  area and volume to `DDF_snow` ∈ {3, 4, 5}, `DDF_ice` across its sourced range and
  `T_s` ∈ {0, 1, 2} °C on the diagnosis set. *Outcome:* fix the starting values. A knob whose
  range moves ice volume by more than 2× is flagged as weakly constrained in the RFC.
- **G3: erosion and its limits.** With D3's law at l ∈ {0.65, 1, 2}, report applied glacial
  erosion, its applied/raw ratio and its 0.25 vs 1 Myr gap.
  *Outcome:* if the neighbor-drop clip still removes more than half of raw abrasion, D3's
  within-step bound is designed explicitly before implementation.
- **G4: sea level.** On the diagnosis set and one ice-age save, report floating, grounded
  marine and land ice under D1–D2 and the eustatic difference between the current debit and
  D6's. *Outcome:* adopt D6's simple form if the generalised correction changes sea level by
  less than 1 m.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| ρ_w / ρ_i | 1000/917 | Densities | None |
| `T_s` (rain–snow threshold) | 1.0 °C | [[Jennings et al., 2018]](#ref-jennings2018) | G2 sensitivity 0–2 °C |
| `DDF_snow` | 4 mm w.e. d⁻¹ K⁻¹ | [[Braithwaite, 2008]](#ref-braithwaite2008) | G2 sensitivity 3–5 |
| `DDF_ice` | To be sourced in G2 | Literature, unverified | G2 |
| G (geothermal flux) | 0.06 W m⁻² | Existing `ICE_GEOTHERMAL_FLUX_W_M2` | Out of scope |
| `K_G`, `l` | 1 × 10⁻⁴, 1 | [[Cook et al., 2020]](#ref-cook2020), [[Humphrey & Raymond, 1994]](#ref-humphrey1994) | G3 sensitivity l ∈ {0.65, 1, 2} |
| `T_f` (seawater freezing) | −1.8 °C | Physical | None |
| `SEA_ICE_MAX_THICKNESS_M` | 3 m | [[Maykut & Untersteiner, 1971]](#ref-maykut1971) | None |
| Endpoint threshold τ / Δt | 0.1 | Design choice | Owner may revise after G1 |
| `GLACIER_MAX_DEPTH_M` | 5,000 m, unchanged | Existing guard | Report how often it binds |

### Cost and saves

- **Runtime.** The ice solve is about 9% of step time today (#353). D2's fixed point adds a
  few outer iterations, each a body solve plus a sea-level bisection and a load update.
  Substepped bodies multiply their own cost.
  - Target: ≤ +15% median step time at density 4. G1 reports the substepped share.
  - D6 and D7 are closed-form per cell.
- **Saves.** No new persisted state. `glacier_depth` has always been metres of ice, so old
  saves' ice loads as it is and re-balances under physical rates on its first steps (#368
  shows it moves within 1–4 Myr). Sea ice is derived, not stored. `lake_ice_m` is RFC 0014's.
  Old saves load (`unit_tests/test_persistence.py`).

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/hydrology.py`:
  - `GLACIER_ACCUMULATION_RATE`, `GLACIER_MELT_RATE_M_PER_MYR`,
    `GLACIER_MELT_REFERENCE_DEGREES_C` and `GLACIER_MELT_MAX_FACTOR` are replaced by
    physical constants; `WATER_TO_ICE` is removed;
  - `compute_hydrology`: the snow fraction replaces `liquid_precip = where(is_frozen, 0, …)`;
    `frozen_from_lake` and the frozen-river-ice block go (shared with RFC 0014);
    `meltwater_m3_per_yr` converts by mass;
  - `allow_ocean_ice`, `SEA_ICE_FORMATION_TEMP_C` and `ICE_AGE_SEA_ICE_MIN_INTENSITY` are
    replaced by a derived sea-ice fraction (new, e.g. `sea_ice_fraction`);
  - `seasonal_positive_degrees`: reused, plus a sibling for the share of the year below a
    threshold (new, shared with RFC 0014's `f_free`).
- `backend/app/hydroclimate/ice_sheet.py`:
  - `SurfaceClimate.rates` in physical units; `MODEL_UNIT_SCALE` and the model-unit
    `GEOTHERMAL_BASAL_MELT_M_PER_MYR` removed;
  - `step_ice`: the sea-ice cap block is removed; an outer fixed point with sea level and
    deflection (new); per-body substepping for slow bodies (new);
  - `IceLedger`: `frozen_river_ice_m3`, `frozen_river_capped_m3` and `frozen_lake_m3` retired;
    a lake-cap hand-off term (new).
- `backend/app/hydroclimate/erosion.py`: `GLACIER_ABRASION_COEFFICIENT`,
  `GLACIER_ABRASION_EXPONENT` and the glacier term's within-step bound.
- `backend/app/hydroclimate/eustasy.py`: `trapped_water_volume_m3` by height above flotation;
  `_ICE_WATER_EQUIVALENT` kept for the conversion.
- `backend/app/lithosphere.py`: `grounded_ice_load_kg_m2` called inside D2's fixed point.
- `backend/app/hydroclimate/biomes.py`: `_smooth_biome_classify` takes land ice for the Ice
  Cap override and the derived sea-ice fraction for `has_sea_ice`.
- `backend/app/render_image.py`: `elevation_colors` and hillshade from surface elevation on
  ice cells; the glacier overlay in `_render_combined_view`; `_biome_fields` carries surface
  elevation and lapse-corrected temperature.
- `backend/app/main.py`: the inspector shows surface and bed elevation.
- `frontend/src/legendData.ts`: legend text for ice and sea ice.
- **Tests:** `test_ice_sheet.py` (physical rates, the fixed point's convergence, substepped
  bodies match the endpoint at steady state, ledger closure), `test_hydrology.py` (snow
  fraction bounds and continuity, no liquid water deleted, no river ice), `test_eustasy.py`
  (floating ice leaves sea level unchanged; marine ice counts above flotation),
  `test_biomes.py` (ice override), `test_render_image.py`, `test_persistence.py`, and a
  `backend/stress_tests/` replay comparing 0.25 and 1 Myr steps.
- **Docs:** `docs/simulation-model/hydrology.md` (glaciers, freezing, sea ice),
  `docs/simulation-model/climate.md` (surface temperature),
  `docs/simulation-model/rendering-and-inspectors.md` (surface elevation, Ice Cap).
- **Diagnostics:** `bin/debug/measure_ice_sheet.py` (#368's scaling flags; G0's knockouts
  and per-step checkpoints, new) and `bin/debug/measure_frozen_water.py` (this RFC).

## Assumptions and uncertainty

- **Endpoint ice.** The premise that most ice volume balances within a step rests on #368's
  rough τ, not a measurement (G1). Large, dry polar sheets have the longest τ.
- **The coupling list** (sea level, deflection, river ice) is inferred, not shown (G0).
  Climate itself doesn't see the ice: there is no ice albedo and no ice-surface orography in
  climate. Those are real feedbacks this RFC leaves out; with them, the coupling problem
  would grow.
- **Precipitation spread evenly through the year** stands in for the seasonality of
  precipitation that `biomes.classify_koppen` already synthesizes. Using that split for
  `f_snow` is a refinement G2 can test.
- **Degree-day melt** is a temperature-index proxy for the energy balance. `DDF_ice` isn't yet
  verified. Shortwave and albedo effects are folded into the factors.
- **Abrasion only.** Quarrying and subglacial water are not represented, and the coefficient
  comes from a few dozen glaciers with erosion rates spanning five orders of magnitude
  [[Cook et al., 2020]](#ref-cook2020).
- **No ice shelves, no ice over lakes.** Floating ice calves.
- **Sea ice** is a derived cover from SST, without a dynamic thickness or salt budget.
- **Two saves, one density, no ice ages** for all evidence so far. The held-out set includes
  an ice-age save.

## Acceptance criteria

- [ ] **G0–G4 run and recorded** in this RFC, with D2's coupling list fixed by G0.
- [ ] **Step robustness:** on each held-out save, ice volume and area at 10 Myr agree within
  10% between 0.25 Myr and 1 Myr steps (proposed).
- [ ] **No unforced swings:** with the ice-age control at "never", no step-to-step change in
  global ice volume exceeds 25% after the first 2 Myr (proposed).
- [ ] **Ledger closes:** glacier ledger residual ≤ 0.1% of the step's accumulation on every
  step of a 10-step replay of each held-out save.
- [ ] **No water deleted:** precipitation leaves runoff only as snow, which is booked as ice
  or meltwater; liquid water lost outside the ledgers is 0.
- [ ] **Units:** `WATER_TO_ICE` and `MODEL_UNIT_SCALE` are gone; every ice–water conversion
  goes through ρ_w / ρ_i in one named place.
- [ ] **Sea level:** a unit test shows floating ice leaves sea level unchanged; G4's
  difference is recorded.
- [ ] **Sea ice:** no ocean cell carries more than `SEA_ICE_MAX_THICKNESS_M` of floating ice,
  and none carries `glacier_depth`.
- [ ] **Maps:** ≥ 95% of land under visible ice shows as Ice Cap in the Biome view, and the
  Elevation view's value on ice cells is bed + ice.
- [ ] **Scale check** (reported, not pass/fail): global ice volume and its sea-level effect
  against the LGM's ~52 × 10⁶ km³ of extra grounded ice and ~134 m
  [[Lambeck et al., 2014]](#ref-lambeck2014), with the world's ocean area stated.
- [ ] **Runtime:** median step time at density 4 within +15% of `origin/main`.
- [ ] **Saves:** old saves load.
- [ ] **Docs** updated as listed above.

## Falsification criteria

- **G0 (e) still swings:** with sea level, deflection and river ice all frozen, physical
  rates still give step-to-step swings above 25% or a 0.25 vs 1 Myr gap above 10%. The
  instability isn't from lagged coupling. *Revise:* substep every body (alternative D per
  body), or find the solver's own fault, before D2.
- **G1 fails broadly:** more than 20% of ice volume has τ > 0.1 Δt. The endpoint premise is
  wrong for this model's ice. *Revise:* make substepping the default and the endpoint an
  optimisation.
- **Physical rates can't hold a plausible ice volume:** after G0's fix, global grounded ice
  on a held-out save with ordinary climate exceeds the ocean's ability to supply it (sea
  level falls more than 1 km below its start in 10 Myr without ice-age forcing). Either the
  climate's polar temperatures and precipitation are too cold or wet for physical ice, or a
  sink is missing. *Revise:* raise it with climate (#358) before shipping physical rates;
  alternative B stays the interim.
- **The maps' mismatch is not an artefact of slow ice:** after D1–D2, more than 20% of EF land
  still carries no ice on a held-out save with positive net snowfall there. *Revise:* the
  snow-phase or melt law, not the maps.

## References

- <a id="ref-braithwaite2008"></a>Braithwaite, R. J. (2008). Temperature and precipitation climate at the equilibrium-line altitude of glaciers expressed by the degree-day factor for melting snow. *Journal of Glaciology*, 54, 437–444. [https://doi.org/10.3189/002214308785836968](https://doi.org/10.3189/002214308785836968)
- <a id="ref-cook2020"></a>Cook, S. J., Swift, D. A., Kirkbride, M. P., Knight, P. G., & Waller, R. I. (2020). The empirical basis for modelling glacial erosion rates. *Nature Communications*, 11, 759. [https://doi.org/10.1038/s41467-020-14583-8](https://doi.org/10.1038/s41467-020-14583-8)
- <a id="ref-ensom2020"></a>Ensom, T., Makarieva, O., Morse, P., Kane, D., Alekseev, V., & Marsh, P. (2020). The distribution and dynamics of aufeis in permafrost regions. *Permafrost and Periglacial Processes*, 31(3), 383–395. [https://doi.org/10.1002/ppp.2051](https://doi.org/10.1002/ppp.2051)
- <a id="ref-goelzer2020"></a>Goelzer, H., Coulon, V., Pattyn, F., de Boer, B., & van de Wal, R. (2020). Brief communication: On calculating the sea-level contribution in marine ice-sheet models. *The Cryosphere*, 14, 833–840. [https://doi.org/10.5194/tc-14-833-2020](https://doi.org/10.5194/tc-14-833-2020)
- <a id="ref-herman2015"></a>Herman, F., Beyssac, O., Brughelli, M., Lane, S. N., Leprince, S., Adatte, T., Lin, J. Y., Avouac, J.-P., & Cox, S. C. (2015). Erosion by an Alpine glacier. *Science*, 350, 193–195. [https://doi.org/10.1126/science.aab2386](https://doi.org/10.1126/science.aab2386)
- <a id="ref-humphrey1994"></a>Humphrey, N. F., & Raymond, C. F. (1994). Hydrology, erosion and sediment production in a surging glacier: Variegated Glacier, Alaska, 1982–83. *Journal of Glaciology*, 40, 539–552. [https://doi.org/10.3189/S0022143000012429](https://doi.org/10.3189/S0022143000012429)
- <a id="ref-jennings2018"></a>Jennings, K. S., Winchell, T. S., Livneh, B., & Molotch, N. P. (2018). Spatial variation of the rain–snow temperature threshold across the Northern Hemisphere. *Nature Communications*, 9, 1148. [https://doi.org/10.1038/s41467-018-03629-7](https://doi.org/10.1038/s41467-018-03629-7)
- <a id="ref-lambeck2014"></a>Lambeck, K., Rouby, H., Purcell, A., Sun, Y., & Sambridge, M. (2014). Sea level and global ice volumes from the Last Glacial Maximum to the Holocene. *Proceedings of the National Academy of Sciences*, 111, 15296–15303. [https://doi.org/10.1073/pnas.1411762111](https://doi.org/10.1073/pnas.1411762111)
- <a id="ref-leppaesranta1993"></a>Leppäranta, M. (1993). A review of analytical models of sea-ice growth. *Atmosphere-Ocean*, 31(1), 123–138. [https://doi.org/10.1080/07055900.1993.9649465](https://doi.org/10.1080/07055900.1993.9649465)
- <a id="ref-maykut1971"></a>Maykut, G. A., & Untersteiner, N. (1971). Some results from a time-dependent thermodynamic model of sea ice. *Journal of Geophysical Research*, 76(6), 1550–1575. [https://doi.org/10.1029/JC076i006p01550](https://doi.org/10.1029/JC076i006p01550)
- <a id="ref-peel2007"></a>Peel, M. C., Finlayson, B. L., & McMahon, T. A. (2007). Updated world map of the Köppen-Geiger climate classification. *Hydrology and Earth System Sciences*, 11, 1633–1644. [https://doi.org/10.5194/hess-11-1633-2007](https://doi.org/10.5194/hess-11-1633-2007)
- <a id="ref-pollard1990"></a>Pollard, D., Muszynski, I., Schneider, S. H., & Thompson, S. L. (1990). Asynchronous coupling of ice-sheet and atmospheric forcing models. *Annals of Glaciology*, 14, 247–251. [https://doi.org/10.3189/S0260305500008685](https://doi.org/10.3189/S0260305500008685)
- <a id="ref-prowse2002"></a>Prowse, T. D., & Beltaos, S. (2002). Climatic control of river-ice hydrology: a review. *Hydrological Processes*, 16, 805–822. [https://doi.org/10.1002/hyp.369](https://doi.org/10.1002/hyp.369)
