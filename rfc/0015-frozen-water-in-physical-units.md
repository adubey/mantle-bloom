# RFC 0015: Frozen water in physical units: glacier mass balance, ice flow and meltwater, freezing lakes and rivers, ocean ice and polar snow

- **Status:** Accepted (2026-10-23)
- **Date:** 2026-10-23
- **Related:** #544 (tracking); #368 (evidence), #353 and #359 (ice-sheet model, decision B), #355 (water units), #503 and #504 (lake freezing at `WATER_TO_ICE`), #297 (balanced lakes), RFC [0014](abandoned/0014-lake-water-volume-and-freezing.md) (lake water volume and seasonal lake ice; its D2(b) and D2(c) are prerequisites here), RFC [0012](0012-physics-first-simulation.md) (physics-first simulation); superseded on acceptance: #197, #203, #358 item 4, #313 Stage E (see "Open issues this RFC supersedes"); design principles 3, 5, 6 and 8

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

So ice is real once it exists, but it forms and melts tens of thousands of times too slowly.
Its meltwater is correspondingly small and mostly never reaches a river, and the
conversions at its borders with water, the ocean and the maps are inconsistent.

**Principles violated:**
- #5, consistent units: one process, two unit systems.
- #6, account for material through every process: snowfall water vanishes, and clipped
  river ice changes reservoir unbooked.
- #3, fast processes in aggregate over long steps: annual-mean switches for precipitation
  phase, freezing and meltwater routing.
- #8, each process follows its own physics: a floating sea-ice cap runs the glacier mass
  balance, up to the 5 km glacier depth guard, and meltwater follows the ground under the
  ice instead of the ice's own drainage.

### Evidence

All replay numbers come from #368 ([`analysis/issue368/`](../analysis/issue368/README.md)),
on `origin/main` `30766af`, with two saves:
- seed 997271774 at 87.4 Myr (SHA-256 `08ad9461…ed762b`);
- seed 649658380 at 6.4 Myr (`ebbfb5da…73f7b3`).

Both are local and **not available for download**; the numbers are quoted so the argument
stands without them. Each arm replayed 10 Myr in 0.25 Myr steps:
- `k1` is the current model;
- `kphys` scales snowfall accumulation, basal melt and `WATER_TO_ICE` by 54,526, and surface
  melt by the same factor;
- `phys` uses a degree-day surface melt (139,408×) and physical abrasion.

Neither save has ice-age forcing.

The map, polar-ocean and sea-level numbers come from a read-only probe,
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
   - The per-body solve is already implicit: `ice_sheet._solve_plastic_body` is a
     backward-Euler volume update, with the ice-surface lapse inside the mass balance. So the
     swings need a mechanism that defeats an implicit solve: an input taken from last step's
     ice or sea level, or the solver's choice between roots (D2, G0).
   - **Response time differs between the saves.** Mean ice volume over mean accumulation in
     `phys` is ≈ 10 kyr on seed 997271774 (295 M km³ ÷ 29,966 km³/yr), but ≈ 110 kyr on seed
     649658380 (43.3 M km³ ÷ 399 km³/yr), 0.43 of a 0.25 Myr step. On that save, surface melt
     (541 km³/yr) exceeds snowfall accumulation, and frozen river ice (157 km³/yr kept, 8,986
     km³/yr clipped at the depth guard) is a large share of what feeds the ice.
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
   - On a node whose annual mean is below 0 °C, or that carries ≥ 10 m of last step's ice,
     `compute_hydrology` removes all precipitation from physical runoff (`liquid_precip`), but
     the ice gains only 1/54,526 of it.
   - In the current model that is about 7,100 km³/yr on seed 997271774 (~225,000 m³/s, a
     quarter of discharge into the ocean) and about 1,200 km³/yr on seed 649658380. The water
     leaves no reservoir; it is simply not routed.
5. **Precipitation phase is an annual-mean switch.** On land below 0 °C annual mean, all
   precipitation is snow. Under the model's own sinusoidal annual cycle, 28% (seed 997271774)
   and 44% (seed 649658380) of that precipitation falls in months above 0 °C.
6. **Meltwater rarely reaches a river.**
   - The ledger's meltwater has two sources: melt from the ice itself (surface and basal
     melt, overflow, and ice shed at the depth guard), and snow melting where it fell on
     ice-free ground (`_melt_bare_snow`).
   - In the current model, all of it together is 0.38 and 0.33 km³/yr on the two saves,
     about 12 and 10 m³/s worldwide: 54,526× too small to make a river.
   - At physical rates (`phys`) it is 59,652 and 15,829 km³/yr. Snow melting on bare ground
     is 65% and 97% of that (38,691 and 15,280 km³/yr); melt from the ice is 20,961 and
     549 km³/yr.
   - Discharge into the ocean rose by only about 11,300 and 2,700 km³/yr, roughly a fifth of
     all meltwater. That is a net figure: ice area also changed, and with it the rain on
     frozen nodes.
   - The routing explains the gap. `_compute_flow_direction` sets `flow_target = −1` on every
     land node whose annual mean is below 0 °C (`water_frozen`). Melt from the ice happens
     mostly on such nodes, since the model's ice margins sit near −12.5 °C annual mean
     (`hydrology.py`'s melt-threshold comment), and bare-ground snowmelt on such nodes is
     trapped the same way. That water settles where it melts and refreezes as river ice (item 3), so meltwater rivers
     from ice margins seldom form.
   - Melt from the ice also follows the ground under the ice, not the ice's own drainage.
7. **Lakes and rivers freeze at the glacier's exchange.**
   - A standing lake below 0 °C annual mean freezes into glacier ice at `WATER_TO_ICE`
     (`frozen_from_lake`, since #504); a settling river does the same a step later
     (`frozen_river_ice`).
   - Until #504, lakes froze 1:1 and fed the ice about 35× faster than all snowfall (#503).
     The 1/54,526 factor is a patch at the same broken boundary.
   - RFC 0014 (D2) removes both conversions for liquid water, and leaves the hand-off of
     frozen-to-bed lake ice to the glacier "until glacier mass balance is physical (#368)".
8. **Marine and floating ice are mis-booked in sea level.**
   - `eustasy.trapped_water_volume_m3` debits every metre of `glacier_depth` at 0.917 m of
     water, including the ice-age sea-ice cap over the ocean and marine ice grounded below
     sea level.
   - The connected-ocean mask (`eustasy._ocean_connected_mask`) is elevation-only, so it keeps
     grounded marine cells as ocean capacity. The full debit then counts the capacity under
     the ice as empty ocean and the ice that fills it as trapped water (D6 derives the
     correct debit). Floating ice should leave sea level unchanged, and grounded marine ice
     should count only above flotation [[Goelzer et al., 2020]](#ref-goelzer2020).
   - On both saves today, floating and marine ice is under 0.1% of ice volume (0.4 and 71 km³
     water equivalent floating; 3,551 and 2,803 km³ marine). The ice-age cap, though, runs the
     glacier mass balance per node up to the 5 km guard (`ice_sheet.step_ice`), while
     multiyear sea ice in the central Arctic averages about 3 m
     [[Maykut & Untersteiner, 1971]](#ref-maykut1971).
9. **The ocean temperature can't show perennial sea ice.**
   - The ocean's baseline is `climate.compute_ocean_temperature_baseline`,
     `WATER_TEMP_MIN_C + WATER_TEMP_RANGE_C × insolation` (−2.0 °C + 30 °C × insolation),
     and over the ocean the air temperature is the ocean's in both the cold-start and
     diagnostic climates (`compute_air_temperature`, `compute_air_temperature_diagnostic`).
   - On both saves, the ocean poleward of 60° (14% of the ocean) never falls below +1.2 °C;
     none of it reaches the −2.0 °C floor or seawater's freezing point (−1.8 °C), in the
     annual mean or in the warmest month. The air over it averages 8.4 °C (seed 997271774)
     and 8.7 °C (seed 649658380).
   - With the maritime seasonal half-amplitude (`biomes.seasonal_temp_amplitude` at
     continentality 0: 9–15 °C at 60–89°), the coldest month of the air over that ocean falls
     below −1.8 °C on about half of it (52% and 49%), and its warmest month on almost none
     (under 0.01%). In the diagnostic climate that air is the ocean temperature itself. Even
     with the full ice-age cooling of `ICE_AGE_MAX_COOLING_C` = 8 °C, the pole's warmest month
     stays above freezing. So the ocean temperature can mark seasonal ice at most, never
     perennial ice.
   - Today's Polar Sea Ice class (`biomes.classify_pelagic`, `sst < −1.0`) is therefore empty
     on both saves outside cells carrying `glacier_depth`.
   - Physically this is expected: the water under sea ice stays at its freezing point all
     year, so its temperature can't tell whether the ice persists. The surface above the ice
     can.
10. **The maps show the bed, not the ice, and ice over water is drawn only by accident.**
    - The Elevation view colours each cell by **bed** elevation and paints white where ice is
      at least `GLACIER_VISIBLE_DEPTH_M` (10 m) (`render_image.elevation_colors` and the
      glacier overlay).
    - Climate's air temperature is the bed's (the ice solve lapses it itself:
      `ice_sheet.SurfaceClimate`), and the Biome view classifies Köppen zones from that
      temperature. Land ice affects only ocean cells (`biomes._smooth_biome_classify`'s
      `has_sea_ice = is_ocean & (glacier_depth > 0)`), which also draws grounded marine
      glacier ice as Polar Sea Ice.
    - On seed 649658380, 33% of land under visible ice is classified EF (Ice Cap) and 60% ET
      (Tundra). For 37% of it, the warmest month is above 0 °C at the bed but would be below
      0 °C at the ice surface (mean ice 853 m, 5.5 °C colder).
    - On seed 997271774, 93% of iced land is EF, but 39% of EF land carries no ice.
    - Ice over water is drawn only where it is `glacier_depth`: the ice-age sea-ice cap and
      lakes frozen into glacier ice. Once RFC 0014 moves lake ice into `lake_ice_m`, and D6
      removes the cap, nothing would draw a perennially frozen lake or sea.

### What the evidence does not show

- Which mechanism drives the instability at physical rates. #368 changed every rate together
  and recorded only 1 Myr checkpoints.
- How much meltwater reaches the ocean once routing stops dead-ending it: item 6's fifth is a
  net figure from runs that still had the dead end.
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
- **Precipitation phase depends on the temperature when it falls.** Across the Northern
  Hemisphere, the air temperature at which rain and snow fall equally often averages 1.0 °C,
  ranging from −0.4 to 2.4 °C at 95% of stations [[Jennings et al., 2018]](#ref-jennings2018).
  A year with a sub-zero mean but thawing months gets part of its precipitation as rain.
- **Ice flows, and its bed melts where it is warm.** Glacier flow carries snow from the
  accumulation area to the ablation area. On a thawed bed, the heat that melts ice is the
  Earth's heat plus friction at the sliding base, less the heat conducted up into the colder
  ice. A perfectly plastic profile, which the model already uses (`ice_sheet.py`, after Nye),
  is the equilibrium shape of that flow.
- **Meltwater in temperate ice drains to its margin.** In temperate ice, at the melting
  point, water percolates through a network of passages. Its direction is set mainly by the
  pressure in the ice, governed primarily by the slope of the ice surface and secondarily by
  the bed. The passages become arborescent: a river-like system on the surface, passages
  through the ice, and tunnels along the bed carrying water and sediment. This explains
  ice-marginal meltwater streams [[Shreve, 1972]](#ref-shreve1972). Cold ice has no such
  network, so melt on cold ice stays on the surface until it reaches the margin or warmer
  ice.
- **Coupling schemes depend on which side is fast.** Ice-age models couple ice sheets and
  atmosphere asynchronously: they run the slow ice alone, with forcing extrapolated from
  short synchronous periods of the fast atmosphere [[Pollard et al., 1990]](#ref-pollard1990).
  Mantle Bloom is in the opposite regime: its ice can respond faster than its step. That
  the ice must then see the feedbacks it changes within the step is this RFC's own
  inference, tested by G0.
- **Glacial erosion scales with sliding.** The dimensionless erosion rate of Variegated
  Glacier's surge, metres eroded per metre slid, was of order 1.0 × 10⁻⁴; the same abstract
  calls that erosion extreme next to non-surging glaciers
  [[Humphrey & Raymond, 1994]](#ref-humphrey1994).
  - A global compilation of 38 glaciers fits `E = K_G · u^l` with l = 0.65 (0.69 including
    Meserve Glacier) and K_G = 1 × 10⁻⁴ for E in m/yr, and derives the relationship for use
    in erosion models [[Cook et al., 2020]](#ref-cook2020).
  - One Alpine glacier gives an exponent of 2 [[Herman et al., 2015]](#ref-herman2015). Its
    coefficient isn't in the abstract and isn't verified here.
  - Since K_G's units are m^(1−l) yr^(l−1), each exponent goes with its own fitted
    coefficient.
- **Rivers don't freeze into glaciers.**
  - River ice covers many cold-region rivers for part of the year; its hydrologic effects are
    ice jams, floods and low winter flows [[Prowse & Beltaos, 2002]](#ref-prowse2002).
  - Where winter water is expelled onto frozen ground or ice, it builds aufeis (icings) in
    layers [[Ensom et al., 2020]](#ref-ensom2020).
  - Both are seasonal stores. Liquid water reaching a cold sink pools there (RFC 0014 D2b).
- **Sea ice is thin, set by the air above it, and doesn't move sea level.**
  - Sea-ice growth is classically modelled analytically as conduction through the ice
    [[Leppäranta, 1993]](#ref-leppaesranta1993). Its simplest form, Stefan's law, grows ice as
    the square root of accumulated freezing degree-days.
  - A thermodynamic model of central-Arctic multiyear ice gives an equilibrium average
    thickness of 288 cm, in close agreement with observations. It needs the ocean to supply
    1–2 kcal cm⁻² yr⁻¹ (about 1.3–2.7 W m⁻²) to the ice; an additional 4 kcal cm⁻² yr⁻¹
    (about 5.3 W m⁻²) would make the ice vanish [[Maykut & Untersteiner, 1971]](#ref-maykut1971).
    Metres of floating ice are negligible as a load or a reservoir next to an ice sheet.
  - For sea level, marine ice counts only by its volume above flotation. That convention needs
    care when the bed or sea level move for other reasons
    [[Goelzer et al., 2020]](#ref-goelzer2020).
  - Scale check: at the Last Glacial Maximum, grounded ice was about 52 × 10⁶ km³ greater than
    today, and sea level about 134 m lower [[Lambeck et al., 2014]](#ref-lambeck2014).
- **The bed sinks under ice.** How fast the mantle relaxes depends on its viscosity; in
  low-viscosity regions the relaxation time is years to decades
  [[Whitehouse, 2018]](#ref-whitehouse2018). The model takes ~10 kyr in general
  (`lithosphere.ice_load_deflection`), a value this RFC hasn't sourced.
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
  initial state, not the climate. Meltwater stays too small to make rivers. The model's other
  water is physical, so every new coupling needs a conversion. Principle #5 stays violated.

**C. Physical mass balance in the existing implicit solve, with its remaining lagged inputs
and root choice fixed as G0 finds; meltwater routed through the ice to its margin; sea ice
from the air above the ocean; physical borders with lakes, rivers, the ocean and the maps
(recommended).**

**D. Full sub-stepping of the hydroclimate.** Run climate, hydrology and ice at a short
sub-step inside every tectonic step.
- *For:* closest to the physics.
- *Against:* the hydroclimate costs seconds per step today; ×25–100 sub-steps is out of
  budget (principle #2). Kept as the per-body fallback, not globally.

**E. A shallow-ice flow solve.** Replace the plastic profile with a time-stepped shallow-ice
approximation.
- *For:* transient flow, thinning and surges.
- *Against:* an explicit solve needs steps of years at this resolution, and the plastic
  profile already gives the equilibrium shape that flow reaches within a step. Rejected
  while G1 shows most ice balances within a step.

**Recommendation: C.** It follows RFC 0014's pattern for lakes: an equilibrium endpoint
where the response time is short against the step, integration for bodies where it isn't,
and an auditable ledger. G0 identifies what makes the existing implicit solve swing before
D2's fixes are chosen.

## Design

### Prerequisites and ownership

RFC 0014 D2(b) (liquid water reaches its sink whatever the temperature, and pools there)
and D2(c) (the lake ice cap `lake_ice_m`) are **prerequisites** of D4, D5 and D8. Without
them, removing frozen river ice deletes the water that settles at cold nodes, which breaks
this RFC's "no water deleted" criterion. D1–D3 should also wait for 0014 D2(b): at physical
rates, `WATER_TO_ICE` turns that settled water into kilometres of river ice (Evidence item 3).

If D1–D3 must ship first, the interim rule is: keep today's river and lake freezing, with
`WATER_TO_ICE` held at its model-unit value as a listed exception, and book what the depth
guard clips as ocean water.

Each removal has one owner:

| Item | Removed or changed by |
|---|---|
| `_compute_flow_direction`'s freeze override (`water_frozen`) for liquid routing | RFC 0014 D2(b) |
| `water_frozen`, `lake_frozen`, `freezing_lake`, `frozen_from_lake` in `compute_hydrology` | RFC 0014 D2(b) |
| The frozen-river-ice block for liquid water in `compute_hydrology` | RFC 0014 D2(b) |
| Writing `IceLedger.frozen_river_ice_m3`, `frozen_river_capped_m3`, `frozen_lake_m3` | RFC 0014 stops writing them |
| Deleting those three fields | RFC 0015 D5, after RFC 0014 |
| `WATER_TO_ICE` itself, once 0014 removes its uses | RFC 0015 D5 |
| `is_frozen`'s use as the snow switch, `liquid_precip`, `rates()`'s `snows`, and the breach gate's `WaterBalance` source | RFC 0015 D1 |
| RFC 0014 D2(b)'s "Unchanged: … snow onto a lake's ice cap joins the glacier accumulation" | Replaced by RFC 0015 D4; a note in RFC 0014 records it |
| The ice-age sea-ice cap (`allow_ocean_ice` block in `step_ice`) | RFC 0015 D6 |

### D1. Physical mass balance

**Units.** The physical rates below are metres of ice per year. `SurfaceClimate.rates` keeps
its interface in m/Myr, which its callers expect (`dt_myr` in the body solve, `/ 1e6` in
`_diagnose_flow`, `_melt_bare_snow`, the ledgers), so the m/yr rates are multiplied by 10⁶ in
one place, inside `rates()`. Every conversion to m/yr of ice:

| Quantity | Its unit | To m/yr of ice |
|---|---|---|
| Precipitation `P` | mm/yr water equivalent | × 10⁻³ × ρ_w / ρ_i |
| Degree-day factor × degree-days | mm w.e. d⁻¹ K⁻¹ × K·d/yr | × 10⁻³ × ρ_w / ρ_i |
| A heat flux `q` | W m⁻² | × `W_M2_TO_M_ICE_PER_YR` = Y / (ρ_i L_f) ≈ 0.103 |
| Frictional heat `τ_b u_b` | Pa × m/yr = J m⁻² yr⁻¹ | ÷ (ρ_i L_f), i.e. `W_M2_TO_M_ICE_PER_YR × τ_b u_b / Y` |
| Sliding speed `u_b` | m/yr (`IceStep.sliding_speed_m_per_yr`) | — |

Here Y = 3.15576 × 10⁷ s per year, ρ_w = 1000 and ρ_i = 917 kg m⁻³, and L_f = 334 kJ kg⁻¹.
`W_M2_TO_M_ICE_PER_YR` is one named constant (new), used by D1 and D6. A unit test checks
each row.

- **Accumulation** = snowfall × 10⁻³ × ρ_w / ρ_i. This replaces
  `GLACIER_ACCUMULATION_RATE`'s model unit and removes `MODEL_UNIT_SCALE`.
- **Snow fraction from the seasonal cycle.**
  - Precipitation is assumed spread evenly over the year (a stated proxy; see Assumptions).
    Its snow fraction is the share of the year below the rain–snow threshold `T_rs`:
    `f_snow(T, A) = ½ + arcsin(clip((T_rs − T) / A, −1, 1)) / π`.
  - Where `A` → 0 (a zero-tilt world), `f_snow` falls back to the existing continuous
    `ICE_SNOW_RAMP_C` ramp around `T_rs`, so the root find never meets a jump.
  - On ice, `T` is the ice surface at the trial thickness; elsewhere, the ground.
- **One partition of precipitation per cell.** Each cell's precipitation `P` is split once,
  after the ice solve, into disjoint parts that sum to `P`. With `w` the lake's wet fraction
  (RFC 0014 D1) and `c` the closing ice cover from the body solve (`_plastic_body`'s cover):
  - **lake:** `w · P` goes to the lake, as liquid supply or onto its cap (D4);
  - **ice:** snow `(1 − w) c f_snow(ice) P` is the body's accumulation, and rain
    `(1 − w) c (1 − f_snow(ice)) P` is routed as water on the ice (D8). `rates()` returns both
    at the same trial thickness, so the root sees the same split that is booked;
  - **bare ground:** snow `(1 − w)(1 − c) f_snow(ground) P`, and rain
    `(1 − w)(1 − c)(1 − f_snow(ground)) P`, which is runoff. As today, where the bare part of a
    cell inside a body's domain has a positive balance and so feeds the body (`_plastic_body`'s
    `feeds`, `_Shape`'s `bare` term), its snow is the body's accumulation; elsewhere it melts
    where it fell within the step (`_melt_bare_snow`).

  `SurfaceClimate` receives `(1 − w) P`, so the lake's share isn't also counted as glacier
  accumulation or bare-ground melt. The liquid source routed by hydrology is the sum of the
  rain shares and melt, computed after `step_ice` returns. A unit test checks that the parts
  sum to `P` on every cell, and that the ledger's accumulation and `snow_melted_m3` match them.
- **The breach gate reads a stated lag.** `breaching.breach_depressions` runs before
  `step_ice` and needs a liquid supply. Its `WaterBalance` reads the same partition at the
  **opening** state: `w_open · P + (1 − w_open) [(1 − c_open)(1 − f_snow(ground)) + c_open
  (1 − f_snow(ice at H_open))] · P`, with `w_open` RFC 0014's last-step wet fraction and
  `c_open`, `H_open` the opening ice cover and thickness. The lake's whole share counts, as in
  D4, since snow on a lake is lake water. This is a one-step lag on the cover only, in place
  of today's `is_frozen` switch.
- **Surface melt** = DDF × positive degree-days × 10⁻³ × ρ_w / ρ_i.
  - Positive degree-days = 365.25 × `seasonal_positive_degrees` on the lapse-corrected surface.
  - `DDF_snow` = 4 mm w.e. d⁻¹ K⁻¹ (Braithwaite's range), applied to the year's snow.
  - `DDF_ice` is applied to bare ice once that snow is gone. Its value must be sourced and
    verified before implementation (G2). #368's runs used 7 mm w.e. d⁻¹ K⁻¹, not verified.
- **Basal melt** on the bed is an energy balance, in m/yr of ice:
  `m_b = W_M2_TO_M_ICE_PER_YR × max(0, Q_b)`, with the net basal heat in W m⁻²
  `Q_b(H) = G + Q_f − k · max(0, T_pmp(H) − T_surf(H)) / H`.
  - `G` = 0.06 W m⁻² geothermal flux. Alone, it melts about 6.2 mm/yr of ice.
  - `Q_f = τ_b u_b / Y` is frictional heating at the sliding base, in W m⁻². `τ_b` is the
    profile's yield stress (`ice_sheet.yield_stress_pa`, 30–100 kPa) and `u_b` the sliding
    speed in m/yr. At 100 kPa and 100 m/yr it is 0.32 W m⁻², about 33 mm/yr of ice.
  - `k · max(0, T_pmp − T_surf) / H` is the heat conducted up into the ice, with
    `k = ICE_THERMAL_CONDUCTIVITY_W_M_K` = 2.1 W m⁻¹ K⁻¹, `T_pmp` the pressure-melting point
    (`−ICE_PRESSURE_MELT_C_PER_M · H`) and `T_surf` the ice-surface temperature. Under 1.5 km
    of ice and a −20 °C surface it is about 0.026 W m⁻², which cuts geothermal melt by about
    40%. Ice can't be warmer than its melting point, so a surface above `T_pmp` conducts
    nothing to the bed: its heat goes into surface melt, which the degree-day term already
    counts. The clamp also keeps the term bounded as H → 0.
  - **Friction is lagged as a whole.** `Q_f` is last step's frictional heat, `τ_b u_b / Y`
    with both `τ_b` and `u_b` from last step's state (the hydrology cache's `IceStep`), and 0
    on nodes that didn't exist then. Taking `τ_b` from this step would make the test below
    circular, since `τ_b` depends on the thermal state and the thermal state on `τ_b`. The lag
    is G0 arm (j)'s subject.
  - **One thermal state.** With `Q_f` fixed for the step, the bed is thawed where
    `Q_b(H) ≥ 0`, and frozen elsewhere. Each use reads it at its own thickness:
    - melt, inside `rates()`, at the trial thickness, so melt is zero exactly at the boundary
      and continuous across it, and the `ICE_THAW_RAMP_C` ramp goes;
    - yield stress (`yield_stress_pa`'s softness) and body membership, at the opening
      thickness, as today (G0 arm (f) tests re-evaluating them at the solved thickness);
    - the diagnostic `warm_based` and sliding, at the closing thickness.

    Without friction this is today's test: `_bed_excess_c ≥ 0` is the condition
    `G ≥ k (T_pmp − T_surf) / H`.
  - **Inside the root.** `m_b` is a per-cell rate in `rates()`, so the solve's balance `B(d)`
    includes it.
- **Bounds:** every rate ≥ 0; `f_snow` ∈ [0, 1]; melt in a step can't exceed the opening ice
  plus that step's accumulation (the solve's existing volume contract).
- **Meltwater** is converted by mass: ice volume × ρ_i / ρ_w (the TODO at
  `hydrology.WATER_TO_ICE`), and routed by D8. Calved ice leaves for the ocean unrouted, as
  today.

### D2. Ice flow and the end-of-step solve

**What the solve does.** Each ice body is solved implicitly. `_solve_plastic_body` finds the
margin position `d` with `V(d) = V_open + Δt · B(d)`, where `B` is evaluated on the trial
plastic surface (`_shape` calls `SurfaceClimate.rates` with the mean thickness). The plastic
profile is the equilibrium shape of ice flow. The balance flux routed down that surface gives
each cell's depth-averaged velocity and, on a thawed bed, a sliding speed
(`ICE_SLIDING_FRACTION` = 0.7 of it; `_diagnose_flow`). Till rides the ice-surface flow
(`hydrology._ice_surface_flow_target`). Under D1 these velocities become physical m/yr.

**What is lagged or ambiguous** (each is a G0 arm):
- the snow switch from last step's ice (fixed by D1);
- body membership and yield stress, from last step's ice (`reference_balance` and the
  thawed-bed test at the opening thickness);
- the root: `B(d)` can change sign more than once, and the bracket walk keeps "the root
  nearest last step's extent". When the opening volume stops mattering (fast ice), either
  branch can win, and the winner can flip between steps;
- sea level, read by flotation, calving, eligibility and climate at the step's opening
  value;
- the mass-balance surface, which uses last step's isostatic deflection.
  `erosion._apply_ice_load` already brings the bed into balance with the closing ice in the
  same step;
- frictional heat, from last step's yield stress and sliding speed (D1);
- frozen river and lake ice, injected a step late (removed by D4, D5).

**Fixes, each adopted only if G0 shows it is needed:**
- **Root choice.** Every step enumerates every sign change of the body's residual on
  `[d_min, d_max]` and logs, per body, the number of stable roots (where `B_V` decreases with
  `V`) and which one was kept. If G0 shows root choice causes the swings, the rule becomes
  the fixed rule that passed (largest stable root, or smallest), and every case with more than
  one stable root is logged.
- **Membership and yield stress** are re-evaluated at the solved thickness and re-solved, at
  most 5 times; a body that hasn't converged keeps its last iterate and is logged.
- **Sea level** becomes an outer unknown. Bracket `g(h) = h − SL(V_trapped(h)) = 0` with
  Brent's method, to 1 m tolerance and at most 8 evaluations. `V_trapped(h)` is D6's
  above-flotation debit after the ice solve at `h`.
  - At each trial `h`, `is_ocean` is recomputed from the connected-ocean mask, so ice can
    ground on newly exposed shelf within the step: that is the main positive feedback. The
    flotation, eligibility (`eligible`) and membership (`positive`) masks follow it.
  - Only bodies that are marine at some `h` in the bracket are re-solved at each trial; the
    rest keep their solution. G1 counts bodies that change marine status within the bracket.
    Their climate's relief term stays at the step's value.
  - Beds: flotation and the thickness solve use the unloaded bed, elevation − deflection, as
    `is_grounded` does today. Eustasy's hypsometry and D6's water depth use the deflected
    elevation it already integrates.
  - The feedback is positive (lower sea level grounds more ice), so convergence comes from
    the bracket, not from fixed-point iteration. A bracket that can't be found (no sign
    change within ±500 m) is logged and the step falls back to the opening sea level.
  - **Discontinuities.** As `h` crosses a sill or a node's threshold, the connected-ocean,
    flotation, eligibility and membership masks change, and `g(h)` can jump, so a sign change
    need not hold a zero. The solve therefore checks water, not bracket width. Its residual is
    the water-volume imbalance `r(h) = C(h) − (V_total − V_trapped(h))`, expressed as metres
    of sea level over the ocean area.
    - Converged: `|r|` ≤ 1 m at the final `h`.
    - Otherwise the final bracket spans an active-set transition: the masks differ between
      its two ends. The solve evaluates `r` on both sides, keeps the side with the smaller
      `|r|`, and logs the step, the masks that changed, and both residuals.
    - G1 reports the distribution of `|r|`, the share of transition steps, and the largest
      transition residual, including steps where body membership changed.
- **Deflection** for the mass-balance surface is updated inside the same outer loop.

**Bodies slower than the step.** A body's response time is `τ = 1 / (−∂B_V/∂V)` at its root,
measured by finite difference: re-solve with `V_open` perturbed by ±5% and ±10%, with the
root tolerance tightened tenfold for the measurement, and report both. (For a fast body,
`∂V_close/∂V_open = τ / (τ + Δt)`, about 0.09 at τ = 0.1 Δt, so smaller perturbations fall
within today's `ICE_ROOT_TOLERANCE_M`.) `H/|b|` is reported alongside.
- A body with τ > 0.1 Δt is integrated in substeps of at most τ/4.
- A body on an unstable branch (`∂B_V/∂V ≥ 0`) has no decay time; it is integrated in
  substeps of a quarter of its growth time `1 / (∂B_V/∂V)`, with a floor of Δt/100.
- During a body's substeps, the other bodies stay at their endpoints, and sea level and
  deflection are updated each substep if G0 kept their coupling.

**Timestep robustness:** rates are annual, and the endpoint doesn't depend on Δt for fast
bodies. Substepped bodies integrate along their path.

### D3. Glacial erosion

- **Law:** `E = K_G · u_b^l`, with `u_b` the physical sliding speed (m/yr) and E in m/yr.
  Each exponent takes its own fitted coefficient:
  - **start:** K_G = 1 × 10⁻⁴, l = 0.65 [[Cook et al., 2020]](#ref-cook2020), the
    38-glacier fit derived for erosion models;
  - **alternative:** K_G = 1 × 10⁻⁴, l = 1 [[Humphrey & Raymond, 1994]](#ref-humphrey1994),
    a surging glacier's rate;
  - **l = 2** [[Herman et al., 2015]](#ref-herman2015) only with its coefficient verified from
    the paper; otherwise this option is dropped.

  This replaces `GLACIER_ABRASION_COEFFICIENT` = 430 per model-unit velocity.
- **Integrated along the step.** A sliding rate times a long step overshoots the neighbor-drop
  limit, which then does the law's work (Evidence item 2). Glacial erosion deepens the bed
  only until the ice's surface slope, and so its sliding, changes. Within a step, applied
  glacial erosion is bounded by the relief the ice can remove, and that bound should be the
  physical one (valley-floor grade under the ice surface), not a clip. G3 measures how much
  applied erosion still depends on step size and on the clip.
- **Calibration.** The current coefficient was set to match the old depth × slope law
  (#353). Under physical sliding, `K_G` is literature-constrained, and the global glacial
  erosion volume is reported against `origin/main` as a sanity check, not a target.
- **Flattening.** The glacial flattening blur (`GLACIER_FLATTEN_RATE_PER_MYR`, scaled by ice
  depth up to `GLACIER_EROSION_MAX_FACTOR`) reads ice thickness, not a mass-balance rate, so
  the unit change doesn't rescale it. It has no literature value, and physical rates change
  where ice sits. G3 reports its volume beside abrasion's and its share of
  `glacier_erosion_multiplier`'s effect.
- **The tuning knob.** `World.glacier_erosion_multiplier` stays a dimensionless user
  multiplier on both terms, 1.0 by default. #197 and #203 measured its effect under the
  model-unit law; G3 re-measures it under this law.

### D4. Lakes that freeze, and lakes under ice (with RFC 0014)

RFC 0014 owns lake water: its volume, seasonal ice cover and the bounded lake ice cap
`lake_ice_m` (D2). This RFC owns the border with the glacier.
- **No lake water becomes glacier ice by exchange factor.** `frozen_from_lake` at
  `WATER_TO_ICE` goes (RFC 0014 D2(b)).
- **Snow on a lake belongs to the lake.** Snow on a cell's wet fraction `w` (D1's partition)
  joins the lake: onto `lake_ice_m` in water equivalent while the surface is frozen,
  otherwise as lake supply. Snow is a mass input to the cap's budget, separate from its
  conductive growth law. Only the dry part `(1 − w)` feeds glacier accumulation. A perennially
  frozen lake is then fed, not emptied, by its own snowfall, as Lake Untersee is (RFC 0014's
  example).
- **At the cap's bound.** RFC 0014 bounds `lake_ice_m` by the cell's water depth. Once the cap
  reaches it, the cell is frozen to its bed and holds no liquid, so further snow on it is
  glacier accumulation on that cell.
- **A cap hands over to the glacier by mass.** Its ice joins `glacier_depth` as
  `lake_ice_m × ρ_w / ρ_i`, and `lake_ice_m` is zeroed in the same step, in either of two
  cases:
  - **glacier ice covers the cell** (visible `glacier_depth`, from an advancing body or from
    snow above the cap's bound): the cap hands over whatever the cell's balance, so the cell
    never carries both reservoirs;
  - **a bare cap frozen to its bed:** it hands over only when `lake_ice_m` reaches the water
    depth **and** the cell's glacier mass balance at that thickness is ≥ 0, so the ice would
    persist as glacier.

  One owner at a time; the hand-off unit test covers both cases.
- **Every ice-covered cell is ice, not lake basin, to the lake solve.** On every cell carrying
  visible glacier ice (≥ `GLACIER_VISIBLE_DEPTH_M`), handed-off or not, the lake hierarchy and
  hypsometry (`build_lake_hierarchy`, `CellHypsometry`) read the bed plus the ice, and liquid
  routing on it follows D8.
  - **Which ice:** the hierarchy, hypsometry, spill routing and `_compute_flow_direction` all
    run before `step_ice` in `compute_hydrology`, and are built once per step. They read the
    **opening** ice: a one-step lag, like the breach gate's. D8's φ routing exists only after
    the solve and reads the **closing** ice; at the margin it hands water to the ground's
    `flow_target`, built on the opening ice.
  - **The lag's measure and trigger.** G1 measures the liquid `lake_depth` under closing
    glacier ice, as a share of global lake volume, at every step. If its maximum over the
    steps of any diagnosis save exceeds 1%, the hierarchy is rebuilt on the closing ice after
    the solve instead, at the cost G1 reports.
  - **The one exception to one owner per cell.** For the step in which ice covers it, a cell
    holds both closing `glacier_depth` and liquid `lake_depth`. For that step lithosphere
    loads the ice (`grounded_ice_load_kg_m2`, as today), eustasy counts the lake water as
    lake water and the ice by D6, and the maps draw the ice. The next step the water spills.
  - So a ground depression under an established ice sheet is not a lake basin, and holds no
    water under the ice (consistent with D8's omission of subglacial lakes).
  - When glacier ice advances over a lake, the `lake_depth` it now covers is booked as lake
    water for the step in which it was covered, and spills through the lake solve on the next
    step.
  - A cold closed basin with steady inflow reaches a steady state instead of converting its
    inflow into ice every step.

  Unit tests check the last two.
- **Glacier meltwater into lakes** is booked by mass (D1, D8), so lake budgets close with it.
- Ice doesn't float on lakes in this design (see Assumptions).

### D5. Rivers that freeze

- **No river water becomes glacier ice.** `frozen_river_ice` and `WATER_TO_ICE` are retired
  (after RFC 0014 D2(b); see Prerequisites). Liquid water reaching a cold sink pools there,
  and the lake's cap holds it as ice by D4's rules. This removes the depth guard's 17–57% of
  meltwater in #368's physical runs.
- **River ice cover** is seasonal [[Prowse & Beltaos, 2002]](#ref-prowse2002). Over a step
  of a year or more it changes when water arrives, not how much, so it isn't represented.
- **Aufeis** [[Ensom et al., 2020]](#ref-ensom2020) is a seasonal store that melts out each
  summer where the year has positive degree-days. Where a cell has none, the water already
  pools and freezes through D4. No separate store is added.
- **Snowfall leak:** with D1, precipitation leaves runoff only as snow, and snow becomes ice
  or meltwater by mass. Diagnostic: liquid water lost outside the ledgers = 0.

### D6. The ocean

- **Sea level counts ice above flotation.** In water equivalent, the debit per marine ice
  cell, per unit area, is `max(0, (ρ_i / ρ_w) · H − D)`, with `D = max(0, h − z)` the water
  depth over the deflected elevation `z` that eustasy's hypsometry integrates. Land ice
  (`D` = 0) debits `0.917 H`, as today.
  - **Derivation.** Eustasy solves `C(h) = V_total − V_trapped`, where `C(h)` is the
    *geometric* capacity below `h` summed over the connected-ocean mask
    (`eustasy._ocean_connected_mask`). The mask is elevation-only, so it includes every
    ice-covered marine cell, with capacity `D`. The physical constraint is on *liquid* water:
    the liquid ocean is the total water minus the water held as ice (`0.917 H`, by mass) and
    in lakes, and it fills only the capacity the ice doesn't occupy. Per marine ice cell:
    - **grounded** (`0.917 H > D`): the ice fills the whole capacity `D` and holds no liquid,
      so the liquid ocean is `C_other(h)`, the capacity of every other cell. The physical
      equation is `C_other(h) = V_total − 0.917 H − lakes`. The model's left side is
      `C(h) = C_other(h) + D`, so the debit that makes the two agree is
      `V_trapped = 0.917 H − D` (plus lakes): the ice's volume above flotation.
    - **floating** (`0.917 H ≤ D`): the ice displaces liquid equal to its own mass, `0.917 H`,
      and liquid fills the remaining `D − 0.917 H`. The liquid ocean is
      `C_other(h) + D − 0.917 H = V_total − 0.917 H − lakes`, so `C(h) = V_total − lakes`: a
      debit of 0.
    - Both cases are `max(0, 0.917 H − D)`. Today's debit of `0.917 H` on grounded marine cells
      is too large by `D`.
  - A unit test checks both cases, including ice grounded exactly at flotation, where the
    debit is 0 and sea level doesn't move.
  - If a later change excluded ice cells from the mask, this debit would have to change with
    it; a comment at both places records that.
  - The generalised approach of [[Goelzer et al., 2020]](#ref-goelzer2020) handles bed and
    sea-level changes not caused by the ice. Adopt it only if G4 shows the simple form drifts.
- **Sea ice is set by the air above it, not by the water below.** The ocean temperature can't
  mark perennial ice (Evidence item 9). Each ocean cell instead gets a derived sea-ice state
  from a degree-day balance:
  - **Surface temperature over ice, `T_si`.** Where ice covers the sea, the ocean's heat no
    longer reaches the air, and the surface behaves like land at sea level. `T_si` is the
    climate's land-temperature construct at sea level for the cell's insolation
    (`climate.compute_land_temperature` at elevation 0), with the ice-age cooling the climate
    already applies, **plus** the cell's ocean-current anomaly, `ocean_temperature_c −
    ocean_baseline_c`, the heat the model's currents carry there. A warm current then holds
    the ice edge poleward and a cold one pushes it equatorward, so the state isn't a function
    of latitude alone. Its seasonal half-amplitude `A` is the maritime one
    (`biomes.seasonal_temp_amplitude` at continentality 0). Both are stated proxies, checked
    by G5.
  - **Freezing and thawing degree-days** per year, in closed form from the same sinusoid:
    `FDD` below `T_f` (`SEAWATER_FREEZING_C` = −1.8 °C, new), and `TDD` above 0 °C
    (`seasonal_positive_degrees` at threshold 0), each × 365.25.
  - **Annual growth** by Stefan's law, `h_g = sqrt(2 k FDD · 86,400 / (ρ_i L_f))`, with `k`
    the ice conductivity. At 4,000 °C·days this gives about 2.2 m.
  - **Annual loss** = surface melt `DDF_ice × TDD × ρ_w / ρ_i` (D1's factor) plus bottom melt
    from the ocean heat flux, `W_M2_TO_M_ICE_PER_YR × F_w` m/yr. `F_w` starts at 2 W m⁻²,
    inside the 1.3–2.7 W m⁻² Maykut & Untersteiner's Arctic ice needs
    [[Maykut & Untersteiner, 1971]](#ref-maykut1971), so bottom melt is about 0.21 m/yr.
  - **State:** *perennial* where `FDD > 0` and `h_g` exceeds the annual loss, so ice survives
    the summer and multiyear ice builds; *seasonal* where `FDD > 0` but it doesn't;
    *ice-free* where `FDD` = 0. The ice-covered share of the year, `f_ice`, is 1 for
    perennial cells and the share of the year below `T_f` otherwise.
  - Ice-age cooling lowers `T_si`, so sea ice grows in an ice age without a separate cap.
  - **Scale check.** At Earth tilt, without the current anomaly, this balance gives perennial
    ice poleward of about 66° (65.75° on a 0.25° grid), seasonal ice equatorward of that to
    about 52°, and open water at 50°. Poleward of 70°, ice grows 1.9–2.9 m in a winter from
    open water. That is Stefan's upper bound for
    first-year ice without snow cover; multiyear ice equilibrates near 288 cm
    [[Maykut & Untersteiner, 1971]](#ref-maykut1971), the same magnitude but a different
    quantity. The 65.75° edge is for `DDF_ice` = 4 mm w.e. d⁻¹ K⁻¹. At 7 it is 67.4°, and at
    G2's starting value of 8 it is 67.9° (G5; `measure_frozen_water.py`'s 5° latitude rows
    could not separate them).
- **What reads sea ice:** the maps (D7: perennial sea ice drawn as ice; the Polar Sea Ice
  class where `f_ice` ≥ 0.5) and the inspector. On perennial sea-ice cells the Temperature
  view and the inspector show `T_si` as the surface temperature, not the air over the open
  ocean. Climate doesn't read sea ice: there is no ice albedo, and the ocean temperature isn't
  changed (see Assumptions). Sea ice has no persisted
  thickness: metres of floating ice are negligible as a load or reservoir
  [[Maykut & Untersteiner, 1971]](#ref-maykut1971), and floating ice doesn't move sea level.
- **One freezing point.** `SEAWATER_FREEZING_C` replaces `SEA_ICE_FORMATION_TEMP_C` (−1.0 °C)
  and `biomes.classify_pelagic`'s `sst < −1.0` cutoff. The Polar Sea Ice class then follows
  `f_ice`, not the ocean temperature; G5 reports its area against today's, which is empty
  on both evidence saves.
- **Grounded marine glacier ice is glacier ice.** Ice on an `is_ocean` cell that is grounded
  (`ice_sheet.is_grounded`) is part of an ice body and is drawn and counted as glacier ice;
  only floating ice is excluded, and floating ice calves.
- **The ice-age sea-ice cap goes.** The `allow_ocean_ice` block in `step_ice`,
  `ICE_AGE_SEA_ICE_MIN_INTENSITY` and its accumulation up to the 5 km guard are removed.
- **Ice shelves are out of scope.** Ice that floats calves, as today. Calved ice returns to
  the ocean by mass.

### D7. Polar snow and ice on the maps

- **Year-round ice is drawn wherever it is, land or water.** In the Elevation, Combined and
  Biome views, four kinds of perennial ice are drawn as ice:
  - glacier ice of at least `GLACIER_VISIBLE_DEPTH_M` on land;
  - grounded marine glacier ice on ocean cells (D6), the same way as land ice;
  - perennial lake ice: a lake cell whose ice-free fraction (RFC 0014's `f_free`) is 0 and
    whose `lake_ice_m` > 0;
  - perennial sea ice (D6).

  Ice over water (lake and sea ice) stays flat-shaded, as the Combined view already does for
  water ("water reads flat"). Seasonal sea ice is not drawn as ice; it shows through the
  Polar Sea Ice class.
- **Elevation.** The Elevation view and the inspector show the **surface** elevation,
  bed + ice, on glacier ice, land or marine: its colour, hillshade and value. The bed stays
  available in the inspector. A plastic ice sheet's smooth dome then reads as a dome, not
  white paint over buried relief.
- **Temperature.** The Temperature view and the climate fields used by maps read at the
  surface, lapse-corrected by the ice: the same correction the ice solve applies, from one
  shared function. The mass balance keeps using its own per-evaluation thickness, so nothing
  is lapsed twice.
- **Biome.** Cells drawn as glacier ice or perennial lake ice are shown as Ice Cap, and cells
  with perennial sea ice as Polar Sea Ice, whatever their Köppen class. Their surface climate
  class stays in the inspector. EF land without ice is still possible (cold and dry) and is
  reported, not overridden.
- **Seasonal snow** isn't drawn: it melts within the year by D1's degree-day term. A
  snow-cover-duration layer could later use `f_snow` and the melt season (out of scope).

### D8. Meltwater rivers

- **No dead ends at cold nodes.** With RFC 0014 D2(b), liquid routing no longer stops at land
  nodes below 0 °C annual mean. That is what keeps both kinds of meltwater out of rivers
  today (Evidence item 6). Bare-ground snowmelt then follows the ground's own routing like
  any runoff.
- **Water on and under the ice follows the ice's drainage.** Melt from the ice and rain on ice
  (D1) are routed according to the ice's thermal state:
  - **warm-based ice** (thawed bed, `ice_sheet.warm_based`): surface melt, rain and basal melt
    follow steepest descent of the hydraulic potential `φ = ρ_w g z_b + ρ_i g H`, with `z_b`
    the bed and `H` the ice thickness [[Shreve, 1972]](#ref-shreve1972);
  - **cold-based ice:** surface melt and rain follow the ice surface
    (`hydrology._ice_surface_flow_target`, which already routes till), until they reach the
    margin or warm-based ice, where they join the φ routing. Basal melt is zero on a frozen
    bed.

  Either way, the water reaches the ice margin and leaves on the ground's own `flow_target`
  as an outlet stream.
- **No subglacial storage.** A closed low in φ under the ice (where a real ice sheet could
  hold a subglacial lake) stores nothing. Each ice body's φ is depression-filled (priority
  flood, as `build_lake_hierarchy` does for the ground), and water is routed on the filled
  surface, so it passes through a pit within the step and leaves by the pit's lowest spill
  point in φ. There is no new state: nothing persists between steps, nothing enters
  eustasy, and the meltwater ledger has no storage term under the ice. G6 reports how many
  pits were filled and how much water passed through them, to show what storage is left out.
- **By mass and annual.** Meltwater enters routing as water, ice volume × ρ_i / ρ_w. It is
  an annual mean; summer peak flows are not represented, and river visibility keeps its
  existing annual-discharge threshold.
- **Ledger, by source.** The meltwater ledger keeps two sources apart: melt from the ice (and
  rain on ice), routed here, and snow melting on bare ground, routed as runoff. Each source's
  water is accounted for as discharge to the ocean, lake supply, transit loss, or ice again
  only through RFC 0014's lake cap, and each closes.

### Diagnosis gates

**Seed sets** (recorded on the tracking issue before G0, disjoint):
- a **diagnosis set** of ≥ 3 saves for G0–G6, including:
  - seed 997271774 at 87.4 Myr;
  - seed 649658380 at 6.4 Myr;
  - one save replayed with ice ages on, with a period of at least 8 of the longest step
    tested (8 Myr when 1 Myr steps are run), so no arm aliases the forcing. The Ice Age
    Frequency slider stops at 1 Myr (`ICE_AGE_MAX_YEARS` in `ControlsModal.tsx`), so the
    period is set through `POST /world/controls`, which accepts any non-negative
    `ice_age_period_years`, or by the replay script setting `World.ice_age_period_years`.
  - a report-only run of the same save at a period users can pick, 1 Myr, with 0.25 Myr steps,
    so the selectable range is seen at physical rates. It isn't pass/fail: 4 steps per cycle
    sample the forcing coarsely (see Assumptions).
- a **held-out acceptance set** of ≥ 3 other saves, at least one replayed with ice ages on
  under the same period rule, never used for gates or knob choices.

**Chosen (2026-10-23, #548).** All six are density 4, at a tilt of 23.5°, with ice ages off
as saved. Checksums and the survey are in [`analysis/issue548/`](../analysis/issue548/README.md).
- Diagnosis: seed 997271774 at 87.4 Myr, seed 649658380 at 6.4 Myr, and seed 374464662 at
  132.6 Myr, replayed with ice ages at an 8 Myr period (report-only arm: 1 Myr).
- Held out: seed 447276134 at 116.4 Myr (replayed with ice ages at an 8 Myr period), seed
  628977350 at 226.7 Myr, and seed 295569319 at 180.1 Myr.

Every gate run includes an `origin/main` arm on the same save and steps, as the baseline.

- **G0: what makes the implicit solve swing?** At D1's physical rates, on the diagnosis set,
  at 0.25 Myr and 1 Myr steps, with checkpoints every step: over 10 Myr on saves without ice
  ages, and over at least two forcing periods plus 2 Myr of spin-up (≥ 18 Myr at an 8 Myr
  period) on the ice-age save. Every arm runs with D1's lagged frictional heat unless it says
  otherwise.

  **Arms:**
  - (a) today's inputs and root rule;
  - (b) sea level held at the run's opening value for all 10 Myr (feedback removed);
    (b′) sea level coupled within the step by D2's bracket;
  - (c) the mass-balance surface's deflection held at the run's opening value;
    (c′) deflection updated within the step by D2's loop;
  - (d) frozen river and lake ice disabled;
  - (e) D1's `f_snow` in place of today's `is_frozen` switch;
  - (f) membership and yield stress re-evaluated at the solved thickness;
  - (j) the passing combined arm, with frictional heat from the current step's state by an
    outer iteration on the solved shape, in place of last step's (run only once a combined
    arm passes);
  - (g) all of (b′), (c′), (d), (e) and (f);
  - (g−x) for each x in {(b′), (c′), (d), (e), (f)}: (g) without x;
  - (h1) today's inputs, always keeping the largest stable root; (h2) the smallest;
  - (i1), (i2): (g) with (h1) or (h2);
  - (i1−x), (i2−x) for each input fix x: (i1) or (i2) without x, run only when that arm
    passes.

  Every arm enumerates all roots on `[d_min, d_max]` each step and logs, per body, the number
  of stable roots and the root kept. A body with no stable root under (h1) or (h2) takes D2's
  unstable-branch substeps, and is logged.

  **Metrics,** after the first 2 Myr: the *swing*, the largest step-to-step change in global
  ice volume, |ΔV| / V; and the *step gap*, |V(0.25 Myr steps) − V(1 Myr steps)| / V at the
  run's end. On saves without ice ages, an arm *passes* when its swing is below 25% and its
  step gap below 10%. On the ice-age save, whose swing is forced, an arm passes when its step
  gap is below 10% and its ice volume follows the forcing:
  - the second cycle's range of V exceeds the first's by less than 25% (no growing
    oscillation), measured on the 0.25 Myr run;
  - V's lag behind the intensity curve, at its maximum, is at most one step of the longest
    arm (1 Myr).

  An arm passes the gate only if it passes on every diagnosis save. The forcing-following
  condition is diagnosis-only; acceptance carries the step gap and the period rule.

  **Outcome:**
  - An input fix is adopted if its own arm passes where (a) fails, **or** if a combined arm
    ((g), (i1) or (i2)) passes and that arm without the fix fails. The leave-one-out runs from
    whichever combined arm passes. The coupled arms (b′) and (c′) decide the sea-level and
    deflection fixes; (b) and (c) only show whether those feedbacks matter at all.
  - The root rule changes only if (h1) or (h2) passes where (a) fails, or (i1) or (i2) passes
    where (g) fails. It becomes the rule that passed. A correlation between root counts and
    volume jumps in (a) is reported but doesn't decide it.
  - The outer iteration on the current step's frictional heat is adopted if (j) reduces the
    passing combined arm's swing or step gap by more than 5 percentage points, or if that arm
    fails without it on any save where (j) passes; otherwise the lagged heat stays.
  - If (g), (i1) and (i2) all fail, the instability is in the solve itself: D2 substeps all
    bodies (alternative D, per body), and the RFC is revised.

  **Gate result (2026-10-24, #554; prototype on `claude/issue-554-g0-ice-swing` at `ea9a187`):
  the instability is in the solve itself.** Details in
  [`analysis/issue554/`](../analysis/issue554/README.md).
  - Since #550, `origin/main` has no frozen river or lake ice, so arm (d) is `origin/main` and
    every arm includes it.
  - Every single-input arm fails on all three diagnosis saves: (a), (b), (b′), (c), (c′), (e),
    (f), (h1) and (h2). So do four diagnostic arms, each removing one feedback:
    - climate noise frozen;
    - yield stress fixed;
    - no frictional heat;
    - no erosion.
  - The (f) and (c′) in-step iterations never converge.
  - (g), (i1), (i2), the leave-one-outs and (j) weren't run. (g) nests the non-converging
    loops, and the owner chose to test alternative D instead.
  - **Alternative D, run as a march.** Each body is marched in physical-time substeps (adaptive,
    250 yr–10 kyr), with every fast input recomputed between substeps, until it settles or the
    step ends.
    - Seed 649658380: passes (swing 9.4%, gap 0.3%).
    - Seed 997271774: passes with the climate's per-step noise frozen (14.9%, 2.9%); with the
      noise, one 27% drop occurs at a step boundary.
    - The ledger closes to ≤ 4 × 10⁻⁵ of accumulation where the march runs the whole step.
  - **The ice-age save runs away under the march at both step sizes** (gap 9.2%) to
    600–800 M km³, with sea level near −2,800 m, about 15× the LGM excess. That is for G4 and
    the ice-age forcing.
  - The diagnosed sliding speeds aren't physical: the fastest thin margin nodes slide ~7 km/yr.
  - D2's revision (a march, possibly with shallow-ice flow) and the climate's per-step noise
    go to the owner and `rfc-authoring`.
- **G1: where does the endpoint premise hold?** On the diagnosis set, per body: τ from
  `−∂B_V/∂V` at the root (D2's perturbations), `H/|b|` alongside, and for bodies above
  10⁴ km² a direct check that substepping reaches the endpoint volume within 1%. Also report:
  the outer iteration count and per-iteration cost of any loop G0 kept; the sea-level water
  residual `|r|`, the share of active-set transition steps and their largest residual (D2);
  the bodies that change marine status within the sea-level bracket; the cost of the root
  enumeration; the membership iterations and non-convergences; frictional melt's share of
  total melt; and the liquid `lake_depth` under closing glacier ice during D4's lag, as a
  share of global lake volume, with its maximum per save.
  *Pass:* ≥ 95% of ice volume has τ ≤ 0.1 Δt at 0.25 Myr steps. Bodies that fail are
  substepped; report their share and cost. Seed 649658380's `H/|b|` of ≈ 110 kyr suggests it
  may fail there (see Falsification), though frozen river ice, which D5 removes, feeds much
  of that ice.
- **G2: melt factors and snow phase.** Source and verify `DDF_ice`. Run sensitivity of ice
  area and volume to `DDF_snow` ∈ {3, 4, 5}, `DDF_ice` across its sourced range and
  `T_rs` ∈ {0, 1, 2} °C on the diagnosis set. *Outcome:* fix the starting values. A knob whose
  range moves ice volume by more than 2× is flagged as weakly constrained in the RFC.

  **`DDF_ice` sourced (2026-10-23, #548).** Start 8 mm w.e. d⁻¹ K⁻¹, range 5.5–8.2: the
  multi-season fits on bare ice in [[Hock, 2003]](#ref-hock2003), Table 1, and 7.5–8.2 observed
  in West Greenland [[Braithwaite, 1995]](#ref-braithwaite1995). Single-season values reach 20,
  which G5 runs as an extreme. The sensitivity runs wait for the D1 prototype.
- **G3: erosion and its limits.** For each (K_G, l) pair in D3, report applied glacial
  erosion, its applied/raw ratio and its 0.25 vs 1 Myr gap. Report flattening's volume beside
  abrasion's, and the effect of `glacier_erosion_multiplier` at 0.25× and 4× on ice-cap area
  and land volume. *Outcome:* if the neighbor-drop clip still removes more than half of raw
  abrasion, D3's within-step bound is designed explicitly before implementation.
- **G4: sea level.** On the diagnosis set, including its ice-age save, report floating,
  grounded marine and land ice under D1–D2, and the eustatic difference between today's
  debit and D6's. *Outcome:* adopt D6's simple form if the generalised correction changes sea
  level by less than 1 m.
- **G5: sea ice.** On the diagnosis set, with and without ice ages, report the perennial and
  seasonal sea-ice areas from D6's balance, by hemisphere and latitude band, the spread of the
  perennial ice edge's latitude by longitude, and the Polar Sea Ice class area before and
  after. Run `F_w` ∈ {1.3, 2, 2.7} W m⁻² and `DDF_ice` across G2's range. Also compute the
  zonal bands, without the current anomaly, at axial tilts of 15°, 23.5° and 35° from the pure
  functions.
  - **Reference.** An observed climatology of perennial and seasonal sea-ice extent by
    latitude (for example NSIDC's), sourced and verified in this gate. Earth-tilt saves are
    compared with it as a sanity check, since the model's continents and currents differ from
    Earth's.
  - **Edge statistic.** In each hemisphere, the *equivalent-area edge*: the latitude that
    encloses the same area as the perennial ice, poleward of it. A hemisphere with no ocean
    poleward of 60° has no edge, and is reported and excluded.
  - **Expected margin.** Without the current anomaly the zonal edge at Earth tilt is 65.75°
    (pure functions, 0.25° grid). The anomaly can move it by at most about ±4°
    (`advect_ocean_temperature` samples the 30 °C-range ocean baseline 6° upstream, against
    `T_si`'s 55 °C-range land construct). So a late-summer observed edge near 75° would put the
    10° test near its bound before any run; the gate records that margin when the reference is
    sourced.
  - *Outcome:* keep the proxies if, without ice ages on an Earth-tilt save, the equivalent-area
    edge falls within 10° of latitude of the observed one in each hemisphere that has one, and
    the current anomaly moves the local edge by at least 2° somewhere along it (so the ocean's
    heat transport is felt). Otherwise revise the proxies before implementation.

  **Gate result (2026-10-23, #548, commit `402f009`): pass; keep the proxies.** The reference
  is ERA5's 1990–2019 sea-ice concentration (≥ 15% every day of the year is perennial), with
  and without coastal cells filled. The two versions bracket NSIDC's 1990–2019 extents. The
  edge is measured against ocean area poleward of it. Observed perennial edges are
  75.3–76.7°N and 69.5–70.5°S.
  - Without ice ages, at `F_w` = 2 and `DDF_ice` = 8, the model edge is 68.2–68.3°N and
    68.4–68.6°S on the three diagnosis saves: gaps of at most 8.5° and 2.1°. The anomaly
    (±2.7 °C poleward of 50°) moves the local edge by ≥ 2° in each hemisphere of each save.
  - `F_w` across 1.3–2.7 moves the edge by ≤ 1°, and `DDF_ice` across 5.5–20 by about 4°.
    Only `DDF_ice` = 5.5 takes one save's north past 10° (10.4°), against the far end of the
    observed range.
  - Three findings remain open before D6's implementation issue (owner decisions, through
    `rfc-authoring`):
    1. The north is biased 6–10° equatorward, which the knobs can't fix.
    2. Seasonal ice reaches about 50.5° (observed 60–71°), so `f_ice` ≥ 0.5 gives 15–34 M km²
       per hemisphere against 9–15 M km² observed ice-covered for half the year.
    3. At the full 8 °C ice-age cooling, perennial ice covers 63–84% of the ocean poleward of
       50°. The falsification criterion is read as applying without ice ages.
  - Zonal bands, perennial / freezing at `DDF_ice` = 8: 64.4° / 54.9° at 15° tilt, 67.9° /
    51.1° at 23.5°, and 81.9° / 45.9° at 35°.
- **G6: meltwater rivers.** Against the `origin/main` arm on the diagnosis set, report for
  each meltwater source the water reaching the ocean, lakes and transit loss; the number of
  filled φ pits and the water passing through them; the number of visible river networks
  whose source lies within one cell of an ice margin; their largest discharges; and the
  surface melt routed by φ on warm-based ice thicker than 1 km, where a thawed bed is least
  like temperate ice (see Assumptions). *Outcome:* report; the acceptance criteria below set
  the bar.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| ρ_w / ρ_i | 1000/917 | Densities | None |
| `T_rs` (rain–snow threshold) | 1.0 °C | [[Jennings et al., 2018]](#ref-jennings2018) | G2 sensitivity 0–2 °C |
| `DDF_snow` | 4 mm w.e. d⁻¹ K⁻¹ | [[Braithwaite, 2008]](#ref-braithwaite2008) | G2 sensitivity 3–5 |
| `DDF_ice` | 8 mm w.e. d⁻¹ K⁻¹ | [[Hock, 2003]](#ref-hock2003); [[Braithwaite, 1995]](#ref-braithwaite1995) (5.5–8.2) | G2 sensitivity |
| G (geothermal flux) | 0.06 W m⁻² | Existing `ICE_GEOTHERMAL_FLUX_W_M2` | Out of scope |
| k (ice conductivity) | 2.1 W m⁻¹ K⁻¹ | Existing `ICE_THERMAL_CONDUCTIVITY_W_M_K` | None |
| τ_b (frictional melt) | The profile's yield stress, 30–100 kPa | Existing `ICE_YIELD_STRESS_*_PA` | None |
| (K_G, l) | (1 × 10⁻⁴, 0.65) | [[Cook et al., 2020]](#ref-cook2020) | G3 pairs: (1 × 10⁻⁴, 1) [[Humphrey & Raymond, 1994]](#ref-humphrey1994); l = 2 only with a verified coefficient |
| `SEAWATER_FREEZING_C` | −1.8 °C | Physical | None |
| `F_w` (ocean heat flux under sea ice) | 2 W m⁻² | [[Maykut & Untersteiner, 1971]](#ref-maykut1971) (1.3–2.7) | G5 sensitivity |
| `T_si`, `A` over sea ice | Land construct at sea level; maritime amplitude | Design choice (proxies) | G5 |
| Polar Sea Ice class threshold | `f_ice` ≥ 0.5 | Design choice | G5 |
| Endpoint threshold τ / Δt | 0.1 | Design choice | Owner may revise after G1 |
| τ perturbations | ±5%, ±10%, root tolerance ÷ 10 | Design choice | G1 |
| Unstable-branch substep | ¼ growth time, floor Δt/100 | Design choice | G1 |
| Membership iterations | ≤ 5 | Design choice | G1 |
| Sea-level bracket | 1 m tolerance, ≤ 8 evaluations, ±500 m | Design choice | G1 cost report |
| `GLACIER_MAX_DEPTH_M` | 5,000 m, unchanged | Existing guard | Report how often it binds |

### Cost and saves

- **Runtime.** The ice solve is about 9% of step time today (#353).
  - A full `step_ice` per outer iteration would cost ~9% each, so the sea-level loop
    re-solves only bodies that are marine somewhere in the bracket (D2), and runs at all only
    if G0 adopts it.
  - The eustasy bisection (`eustasy._bisect_sea_level`, up to 80 iterations, each rebuilding
    the connected-ocean mask) runs once per outer evaluation.
  - Root enumeration scans each body's `[d_min, d_max]` instead of walking out from last
    step's extent.
  - G1 reports the actual counts and costs. Substepped bodies multiply their own cost.
  - D6's sea ice, D7 and D8's routing are closed-form or one pass per cell, plus a priority
    flood of φ per ice body.
  - If G1 finds more than 1% of lake volume under new ice during D4's one-step lag, the lake
    hierarchy is rebuilt after the solve, adding a second `build_lake_hierarchy` per step;
    G1 reports its cost.
  - Target: ≤ +15% median step time at density 4. If G1 shows the adopted fixes can't meet
    it, the owner decides between the budget and the fix.
- **Saves.** No new persisted state. `glacier_depth` has always been metres of ice, so old
  saves' ice loads as it is and re-balances under physical rates on its first steps (#368
  shows it moves within 1–4 Myr). Frictional melt reads the hydrology cache's sliding speed
  and is 0 where it's missing. Sea ice is derived, not stored. `lake_ice_m` is RFC 0014's. Old
  saves load (`unit_tests/test_persistence.py`).

### Open issues this RFC supersedes

All 55 open issues were searched (2026-10-23) for glacier, ice-sheet, ice-cap, frozen,
freezing, sea-ice, snow, `WATER_TO_ICE`, abrasion, calving and eustasy terms, and each match
was read. Each issue below is mapped to the part of this RFC that replaces it, so a partial
withdrawal can be handled part by part.

| Issue | What it asks | Covered by | On acceptance |
|---|---|---|---|
| #368 | Physical glacier mass-balance scaling and cap behaviour | Whole RFC; its analysis is the Evidence | Already closed (PR #536). Comment there with a link to this RFC. |
| #197 | A dedicated sweep of `glacier_erosion_multiplier` | D3, G3 | Close as superseded: the sweep measures the model-unit abrasion law that D3 replaces. G3 re-measures the multiplier. |
| #203 | Confirm or taper `glacier_erosion_multiplier`'s ice-cap and land-volume lever | D3 (abrasion and flattening), G3 | Close as superseded: it cites the depth × slope law #353 already replaced, and the remaining lever is re-measured by G3. |
| #358 item 4 | Glacier mass balance in physical units | D1, D2 | Keep #358 open (its other items stand). Mark item 4 as superseded by this RFC in its body. |
| #313 Stage E | Glacier mass-balance controls in the calibration campaign | D1–D3, Knobs | Keep #313 open. Replace Stage E's constant list with this RFC's knobs table, and point the stage here. |

**Related but not superseded** (each keeps its own scope; this RFC changes their inputs):
- #176, #289: tectonics that raise land and so grow ice caps. Their ice-cap metrics move under
  D1, so their baselines need re-measuring after implementation.
- #310, #493: craton work that notes seed 580014012's repeated 5 km glacier cap.
  Interpretation changes once the cap rarely binds.
- #505, #519, #520: lake water, owned by RFC 0014. #519's item on rivers freezing at the
  annual mean is RFC 0014's D2(b); this RFC owns only the glacier border (D4, D5, D8).

**(a) When this RFC is accepted** (the owner merges the RFC PR), on the tracking issue:
1. For #197 and #203, comment with the section and gate that covers each, set its project
   status to `Deprecated`, and close it with `gh issue close <N> --reason "not planned"`
   (the project's procedure for work that won't be pursued).
2. Edit #358 and #313 as the table says, with a link to this RFC.
3. Comment on #368 with a link to this RFC.
4. List every closed or edited issue on the tracking issue, with the RFC section each maps
   to.

**(b) If this RFC is withdrawn**, or the part an issue maps to is dropped or falsified:
1. Review each issue the table maps to that part. Read its last state and check whether its
   premise still holds on current `origin/main` (for #197 and #203, whether
   `glacier_erosion_multiplier` still acts on the same law).
2. Reopen it if the premise holds, with a comment linking the withdrawal or the gate result
   and saying what changed since it was closed. Set its project status back to `Backlog`.
   Leave it closed, with a comment saying why, if its premise no longer holds.
3. Revert the #358 and #313 edits, or point them at whatever supersedes this RFC.
4. Record each decision in this RFC's Outcome section (`complete-rfc`).

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/hydrology.py`:
  - `GLACIER_ACCUMULATION_RATE`, `GLACIER_MELT_RATE_M_PER_MYR`,
    `GLACIER_MELT_REFERENCE_DEGREES_C` and `GLACIER_MELT_MAX_FACTOR` are replaced by
    physical constants; `WATER_TO_ICE` is removed after RFC 0014 removes its uses;
  - `compute_hydrology`: D1's three-part precipitation partition after `step_ice`, replacing
    `liquid_precip = where(is_frozen, 0, …)`; the breach gate's `WaterBalance` source at the
    opening cover; `meltwater_m3_per_yr` converts by mass;
  - meltwater routing on ice by thermal state: hydraulic potential over a depression-filled φ
    on warm-based ice, `_ice_surface_flow_target` on cold-based ice (new; no persisted state);
  - the lake solve's inputs read bed + ice on every visibly ice-covered cell (with RFC 0014);
  - `allow_ocean_ice`, `SEA_ICE_FORMATION_TEMP_C` and `ICE_AGE_SEA_ICE_MIN_INTENSITY` are
    replaced by a derived sea-ice state (new, e.g. `sea_ice_state`) and `SEAWATER_FREEZING_C`
    (new);
  - `seasonal_positive_degrees`: reused, plus siblings for the share of the year below a
    threshold and freezing degree-days (new, shared with RFC 0014's `f_free`).
- `backend/app/hydroclimate/ice_sheet.py`:
  - `SurfaceClimate.rates` in physical units, returning snow and rain on ice and the basal
    energy balance; `W_M2_TO_M_ICE_PER_YR` (new); `MODEL_UNIT_SCALE`, the model-unit
    `GEOTHERMAL_BASAL_MELT_M_PER_MYR` and `ICE_THAW_RAMP_C` removed;
  - `warm_based` and `yield_stress_pa` read one thermal state from the full basal balance,
    with last step's frictional heat, each at its own thickness (D1);
  - `SurfaceClimate` takes `(1 − w) P`, so the lake's share isn't counted twice;
  - `_solve_plastic_body`: root enumeration and logging; the root rule as G0 decides;
  - `step_ice`: the sea-ice cap block removed; membership, sea level and deflection updates as
    G0 decides; per-body substepping for slow and unstable bodies (new); last step's frictional
    heat (yield stress × sliding speed) as an input;
  - `IceLedger`: `frozen_river_ice_m3`, `frozen_river_capped_m3` and `frozen_lake_m3` deleted;
    lake-cap hand-off and frictional-melt terms (new).
- `backend/app/hydroclimate/lakes.py`: hypsometry and hierarchy read bed + opening glacier ice
  on every visibly ice-covered cell (with RFC 0014's solve), or are rebuilt on the closing ice
  if G1 requires it; snow on the wet fraction joins the lake; snow above the cap's bound is
  glacier accumulation; the cap hands over whenever glacier ice covers the cell.
- `backend/app/hydroclimate/breaching.py`: `WaterBalance`'s liquid source from D1.
- `backend/app/hydroclimate/erosion.py`: `GLACIER_ABRASION_COEFFICIENT`,
  `GLACIER_ABRASION_EXPONENT` and the glacier term's within-step bound.
- `backend/app/hydroclimate/eustasy.py`: `trapped_water_volume_m3` by height above flotation;
  `_ICE_WATER_EQUIVALENT` kept for the conversion; a comment tying it to
  `_ocean_connected_mask`.
- `backend/app/hydroclimate/climate.py`: `compute_land_temperature` reused at sea level, and
  the ocean-current anomaly (`ocean_temperature_c` minus the baseline) exposed, for `T_si`
  (no change to the ocean temperature).
- `backend/app/hydroclimate/biomes.py`: `classify_pelagic` and `_smooth_biome_classify` take
  the sea-ice state in place of `has_sea_ice`; the Ice Cap override for drawn glacier and lake
  ice, land or marine.
- `backend/app/render_image.py`: `elevation_colors` and hillshade from surface elevation on
  glacier ice; perennial ice drawn on land, lake and sea cells in the Elevation, Combined and
  Biome views; `_biome_fields` carries surface elevation, lapse-corrected temperature,
  `lake_ice_m` and the sea-ice state.
- `backend/app/main.py`: the inspector shows surface and bed elevation, sea-ice state, and
  `T_si` as the surface temperature on perennial sea ice.
- `frontend/src/legendData.ts`: legend text for ice, lake ice and sea ice.
- `rfc/abandoned/0014-lake-water-volume-and-freezing.md`: a note on D2(b)'s snow bullet and the
  `IceLedger` fields, pointing here.
- **Tests:**
  - `test_ice_sheet.py`: physical rates and every row of D1's units table, including
    frictional heat from Pa × m/yr; root enumeration and
    the fixed root rules; basal melt continuous at the thaw boundary, including a cell with
    non-zero sliding; basal melt bounded as H → 0 under a surface above 0 °C; frictional melt
    inside the root; substepped bodies match the endpoint at steady state; ledger closure;
  - `test_hydrology.py`: the three-part partition sums to P per cell; snow-fraction bounds,
    continuity and the A → 0 fallback; no liquid water deleted; melt from warm-based ice
    reaches the margin by φ, and from cold-based ice by the ice surface, and leaves as a river;
    the sea-ice state's perennial, seasonal and ice-free cases;
  - `test_lakes.py`: snow on a frozen lake feeds it; snow above the cap's bound becomes glacier
    accumulation; the cap hands over when glacier ice covers the cell and, bare, only with a
    non-negative balance; a cold closed basin with steady inflow reaches a steady state with no
    per-step hand-off; ice advancing over a lake leaves no `lake_depth` under it after the
    next step, booked as spill;
  - `test_eustasy.py`: floating ice, and ice grounded exactly at flotation, leave sea level
    unchanged; marine ice counts above flotation;
  - `test_biomes.py`, `test_render_image.py`: perennial ice drawn on land, marine, lake and sea
    cells in each view; the surface-elevation colour on glacier ice;
  - `test_persistence.py`; a `backend/stress_tests/` replay comparing 0.25 and 1 Myr steps.
- **Docs:** `docs/simulation-model/hydrology.md` (glaciers, ice flow, meltwater, freezing, sea
  ice), `docs/simulation-model/climate.md` (surface temperature),
  `docs/simulation-model/rendering-and-inspectors.md` (surface elevation, perennial ice).
- **Diagnostics:** `bin/debug/measure_ice_sheet.py` (#368's scaling flags; G0's arms, root
  enumeration and per-step checkpoints, new) and `bin/debug/measure_frozen_water.py`.

## Assumptions and uncertainty

- **Endpoint ice.** The premise that most ice volume balances within a step holds on one save
  and is doubtful on the other (`H/|b|` ≈ 10 and 110 kyr). G1 measures the response time from
  the solve itself.
- **The swing's cause** is unknown; G0's arms cover the inputs that lag and the root choice,
  but not every possibility.
- **Climate doesn't see the ice:** no ice albedo, no ice-surface orography in the climate, no
  sea-ice albedo, and no change to the ocean temperature under sea ice. Those are real
  feedbacks this RFC leaves out; with them, the coupling problem would grow.
- **Precipitation spread evenly through the year** stands in for the seasonality of
  precipitation that `biomes.classify_koppen` already synthesizes. Using that split for
  `f_snow` is a refinement G2 can test.
- **Degree-day melt** is a temperature-index proxy for the energy balance. `DDF_ice` isn't yet
  verified. Shortwave and albedo effects are folded into the factors.
- **Frictional melt** uses last step's frictional heat (its yield stress and sliding speed). The basal energy balance assumes steady
  conduction through the ice, without advection or strain heating.
- **Meltwater routing** follows the hydraulic potential on warm-based ice and the ice surface
  on cold-based ice, at the annual mean. A thawed bed stands in for temperate ice: much of a
  large warm-based sheet is cold above its bed, and its surface water reaches the bed only
  through moulins or fractures. The water ends at the margin either way; G6 reports how much
  surface melt is routed by φ on warm-based ice thicker than 1 km. Seasonal drainage evolution, moulins, supraglacial
  lakes and outburst floods aren't represented, and subglacial lakes are omitted: water passes
  through φ pits within the step.
- **Abrasion only.** Quarrying and subglacial water pressure are not represented, and the
  coefficient comes from a few dozen glaciers with erosion rates spanning five orders of
  magnitude [[Cook et al., 2020]](#ref-cook2020).
- **No ice shelves, no floating lake ice from glaciers.** Floating glacier ice calves.
- **Sea ice** comes from a land-temperature construct standing in for the air over ice, shifted
  by the ocean-current anomaly, a borrowed maritime amplitude, Stefan growth without snow cover
  (an upper bound on first-year growth), and a single ocean heat flux. It has no persisted
  thickness, salt budget, drift or climate feedback. The Temperature view shows `T_si` over
  perennial sea ice, but the climate's own air temperature there stays the ocean's. The
  current anomaly can move the ice edge by at most about ±4° of latitude, far less than
  Earth's strongest currents do.
- **Ice-age forcing is sampled once per step.** The climate reads the ice-age intensity at the
  step's time, so a period shorter than a few steps is aliased. Periods users can pick
  (0.1–1 Myr) are sampled 1–10 times per cycle at 0.1–1 Myr steps.
- **Isostatic relaxation** is taken as short against the step; the model's ~10 kyr isn't
  sourced here.
- **Two saves, one density, no ice ages** for all evidence so far. The diagnosis and
  held-out sets each include an ice-age replay.

## Acceptance criteria

The thresholds below were confirmed by the project owner on 2026-10-23.

- [ ] **Prerequisites:** RFC 0014 D2(b) and D2(c) landed before D4, D5 and D8; or D1–D3 shipped
  first under the interim rule, stated in the PR.
- [ ] **Superseded issues settled** at acceptance as planned in "Open issues this RFC
  supersedes" (a), and listed on the tracking issue.
- [ ] **G0–G6 run and recorded** in this RFC, with D2's fixes chosen by G0.
- [ ] **Step robustness:** on each held-out save, ice volume and area at 10 Myr agree within
  10% between 0.25 Myr and 1 Myr steps. Ice-age saves use a period of at least 8
  of the longest step.
- [ ] **No unforced swings:** with the ice-age control at "never", no step-to-step change in
  global ice volume exceeds 25% after the first 2 Myr.
- [ ] **Ledgers close:** glacier ledger residual ≤ 0.1% of the step's accumulation, and both
  meltwater sources' ledgers (D8) close to round-off, on every step of a 10-step replay of
  each held-out save.
- [ ] **No water deleted:** the three-part partition of precipitation sums to P on every cell;
  liquid water lost outside the ledgers is 0.
- [ ] **Units:** `WATER_TO_ICE` and `MODEL_UNIT_SCALE` are gone; every ice–water conversion
  goes through ρ_w / ρ_i in one named place.
- [ ] **Meltwater rivers:** on each held-out save, at least half of the melt from the ice
  (surface, basal and rain on ice) reaches the ocean or a lake as discharge, and the number of
  visible river networks sourced at ice margins is reported against `origin/main`
  (G6 informs the report). Bare-ground snowmelt is reported separately.
- [ ] **Sea level:** unit tests show floating ice, and ice grounded exactly at flotation,
  leave sea level unchanged; G4's difference is recorded. If D2's sea-level loop is adopted,
  on each held-out save the water residual `|r|` is ≤ 1 m of sea level on at least 99% of
  steps, and ≤ 10 m on every active-set transition step.
- [ ] **Sea ice:** no ocean cell carries floating `glacier_depth` (every ice-bearing ocean
  cell is grounded by `is_grounded`); G5's extents recorded and its outcome met.
- [ ] **Maps:** unit tests show perennial ice drawn on land, grounded marine, lake and sea
  cells in the Elevation, Combined and Biome views, and glacier ice's Elevation value equal to
  bed + ice. As a model check (reported, not pass/fail): the share of drawn land ice whose
  lapse-corrected warmest month is above 0 °C.
- [ ] **Scale check** (reported, not pass/fail): global ice volume and its sea-level effect
  against the LGM's ~52 × 10⁶ km³ of extra grounded ice and ~134 m
  [[Lambeck et al., 2014]](#ref-lambeck2014), with the world's ocean area and the velocity
  distribution of the ice (interior and margins) stated.
- [ ] **Runtime:** median step time at density 4 within +15% of `origin/main`.
- [ ] **Saves:** old saves load.
- [ ] **Docs** updated as listed above.

## Falsification criteria

- **G0's combined arms still swing:** with every input fix applied and either fixed root rule
  ((g), (i1) and (i2)), physical rates still give step-to-step swings above 25% or a 0.25 vs
  1 Myr gap above 10%. The instability isn't from lagged inputs or root choice. *Revise:*
  substep every body (alternative D per body), or find the solver's own fault, before D2.
- **G1 fails broadly:** more than 20% of ice volume has τ > 0.1 Δt. The endpoint premise is
  wrong for this model's ice. *Revise:* make substepping the default and the endpoint an
  optimisation. Seed 649658380 is already close to this on `H/|b|`.
- **Physical rates can't hold a plausible ice volume:** after G0's fixes, with ice ages off
  and each save's own polar climate, grounded ice exceeds the same save's `origin/main` arm by
  more than 3× the LGM's excess (52 × 10⁶ km³, about 134 m of sea level) within 10 Myr.
  Either the climate's polar temperatures and precipitation are too cold or wet for physical
  ice, or a sink is missing. *Revise:* raise it with climate (#358) before shipping physical
  rates; alternative B stays the interim. (#368's `phys` arm, ~800 m below `k1`'s mean sea
  level on seed 997271774, would meet this today; the criterion applies after G0's fixes.)
- **Meltwater still doesn't reach rivers:** after D5 and D8, less than a fifth of the melt from
  the ice reaches the ocean or a lake on a held-out save. The routing premise is wrong.
  *Revise:* D8.
- **Sea ice can't be placed from the air above it:** at every G5 setting, perennial sea ice
  covers no ocean without ice ages on an Earth-tilt save, or covers most of the ocean poleward
  of 50°. *Revise:* D6's `T_si` proxy, or give sea ice an energy balance with the ocean.
- **The maps' mismatch is not an artefact of slow ice:** after D1–D2, more than 20% of EF land
  still carries no ice on a held-out save with positive net snowfall there. *Revise:* the
  snow-phase or melt law, not the maps.

Whenever a criterion here leads to withdrawing the RFC or dropping a design part, the issues
mapped to that part are reviewed and possibly reopened by plan (b) in "Open issues this RFC
supersedes".

## References

- <a id="ref-braithwaite1995"></a>Braithwaite, R. J. (1995). Positive degree-day factors for ablation on the Greenland ice sheet studied by energy-balance modelling. *Journal of Glaciology*, 41(137), 153–160. [https://doi.org/10.3189/S0022143000017846](https://doi.org/10.3189/S0022143000017846)
- <a id="ref-braithwaite2008"></a>Braithwaite, R. J. (2008). Temperature and precipitation climate at the equilibrium-line altitude of glaciers expressed by the degree-day factor for melting snow. *Journal of Glaciology*, 54, 437–444. [https://doi.org/10.3189/002214308785836968](https://doi.org/10.3189/002214308785836968)
- <a id="ref-cook2020"></a>Cook, S. J., Swift, D. A., Kirkbride, M. P., Knight, P. G., & Waller, R. I. (2020). The empirical basis for modelling glacial erosion rates. *Nature Communications*, 11, 759. [https://doi.org/10.1038/s41467-020-14583-8](https://doi.org/10.1038/s41467-020-14583-8)
- <a id="ref-ensom2020"></a>Ensom, T., Makarieva, O., Morse, P., Kane, D., Alekseev, V., & Marsh, P. (2020). The distribution and dynamics of aufeis in permafrost regions. *Permafrost and Periglacial Processes*, 31(3), 383–395. [https://doi.org/10.1002/ppp.2051](https://doi.org/10.1002/ppp.2051)
- <a id="ref-goelzer2020"></a>Goelzer, H., Coulon, V., Pattyn, F., de Boer, B., & van de Wal, R. (2020). Brief communication: On calculating the sea-level contribution in marine ice-sheet models. *The Cryosphere*, 14, 833–840. [https://doi.org/10.5194/tc-14-833-2020](https://doi.org/10.5194/tc-14-833-2020)
- <a id="ref-herman2015"></a>Herman, F., Beyssac, O., Brughelli, M., Lane, S. N., Leprince, S., Adatte, T., Lin, J. Y., Avouac, J.-P., & Cox, S. C. (2015). Erosion by an Alpine glacier. *Science*, 350, 193–195. [https://doi.org/10.1126/science.aab2386](https://doi.org/10.1126/science.aab2386)
- <a id="ref-hock2003"></a>Hock, R. (2003). Temperature index melt modelling in mountain areas. *Journal of Hydrology*, 282, 104–115. [https://doi.org/10.1016/S0022-1694(03)00257-9](https://doi.org/10.1016/S0022-1694(03)00257-9)
- <a id="ref-humphrey1994"></a>Humphrey, N. F., & Raymond, C. F. (1994). Hydrology, erosion and sediment production in a surging glacier: Variegated Glacier, Alaska, 1982–83. *Journal of Glaciology*, 40, 539–552. [https://doi.org/10.3189/S0022143000012429](https://doi.org/10.3189/S0022143000012429)
- <a id="ref-jennings2018"></a>Jennings, K. S., Winchell, T. S., Livneh, B., & Molotch, N. P. (2018). Spatial variation of the rain–snow temperature threshold across the Northern Hemisphere. *Nature Communications*, 9, 1148. [https://doi.org/10.1038/s41467-018-03629-7](https://doi.org/10.1038/s41467-018-03629-7)
- <a id="ref-lambeck2014"></a>Lambeck, K., Rouby, H., Purcell, A., Sun, Y., & Sambridge, M. (2014). Sea level and global ice volumes from the Last Glacial Maximum to the Holocene. *Proceedings of the National Academy of Sciences*, 111, 15296–15303. [https://doi.org/10.1073/pnas.1411762111](https://doi.org/10.1073/pnas.1411762111)
- <a id="ref-leppaesranta1993"></a>Leppäranta, M. (1993). A review of analytical models of sea-ice growth. *Atmosphere-Ocean*, 31(1), 123–138. [https://doi.org/10.1080/07055900.1993.9649465](https://doi.org/10.1080/07055900.1993.9649465)
- <a id="ref-maykut1971"></a>Maykut, G. A., & Untersteiner, N. (1971). Some results from a time-dependent thermodynamic model of sea ice. *Journal of Geophysical Research*, 76(6), 1550–1575. [https://doi.org/10.1029/JC076i006p01550](https://doi.org/10.1029/JC076i006p01550)
- <a id="ref-peel2007"></a>Peel, M. C., Finlayson, B. L., & McMahon, T. A. (2007). Updated world map of the Köppen-Geiger climate classification. *Hydrology and Earth System Sciences*, 11, 1633–1644. [https://doi.org/10.5194/hess-11-1633-2007](https://doi.org/10.5194/hess-11-1633-2007)
- <a id="ref-pollard1990"></a>Pollard, D., Muszynski, I., Schneider, S. H., & Thompson, S. L. (1990). Asynchronous coupling of ice-sheet and atmospheric forcing models. *Annals of Glaciology*, 14, 247–251. [https://doi.org/10.3189/S0260305500008685](https://doi.org/10.3189/S0260305500008685)
- <a id="ref-prowse2002"></a>Prowse, T. D., & Beltaos, S. (2002). Climatic control of river-ice hydrology: a review. *Hydrological Processes*, 16, 805–822. [https://doi.org/10.1002/hyp.369](https://doi.org/10.1002/hyp.369)
- <a id="ref-shreve1972"></a>Shreve, R. L. (1972). Movement of water in glaciers. *Journal of Glaciology*, 11, 205–214. [https://doi.org/10.3189/S002214300002219X](https://doi.org/10.3189/S002214300002219X)
- <a id="ref-whitehouse2018"></a>Whitehouse, P. L. (2018). Glacial isostatic adjustment modelling: historical perspectives, recent advances, and future directions. *Earth Surface Dynamics*, 6, 401–429. [https://doi.org/10.5194/esurf-6-401-2018](https://doi.org/10.5194/esurf-6-401-2018)
