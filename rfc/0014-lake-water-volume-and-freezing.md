# RFC 0014: Lakes conserve their water: in-step fill-spill-merge, seasonal freezing and water-surface temperature

- **Status:** Accepted (2026-10-23)
- **Date:** 2026-10-22 (revised 2026-10-23 after review; follow-up evidence added 2026-10-28)
- **Related:** #538 (tracking), #539 (D2b), #540 (gate G3), #541 (seed sets, gate G5), #615 (D2c design screen), #632 (D2c implementation); #519 (exploration and evidence), #520 (perched-lake false merges, river lines across lakes, floor-drain flips), #521 (split reset), #505 (lake merges lag routing), #503 and #297 (breaching, balanced lakes), #355 (water units), #363 (lake silt), #353 and #368 (glacier ice and its units), #358 items 10–11 (water-surface temperature, one evaporation flux), #118 (wind speed), #675 (large-lake flicker follow-up)

## Problem

Lakes in Mantle Bloom are not conserved bodies of water. Each step `lakes.resolve_lakes` places every lake at its *balanced level*, where evaporation over its flooded area equals its inflow, capped at its rim. That premise is often reasonable: many closed lakes adjust within years to centuries, short next to a 10 kyr–10 Myr step (G1 tests where it holds). But the machinery around it creates and deletes water, and it decides merging from last step's water rather than this step's. On a mature world this produces seas that appear from nowhere, vanish the next step, and hold far more water than their budget supports.

### Evidence

All numbers are from one save, `seed308161711-121.5Myr` (128,340 nodes, 121.5 Myr, 100 kyr steps, `origin/main` `59d2ef3`). The save is local and **not available for download**. The probes, knockouts and their preregistrations are on #519; the numbers below are quoted so the argument stands without the save. The evidence is **single-save**: no other seed, age or resolution has been tested (see diagnosis gate G0).

1. **The split reset creates water.** When a lake resolved as one merged body balances below its saddle, `lakes._resolve` resets *both* children to the saddle level without a balance check (`reset_level = new_level if tier_capped else node.min_depth`).
   - Resets added **4.7–5.5 M km³ of water per step**, about the world's entire visible lake volume.
   - A 4,028-member basin (sea 26) went from 1,231 km³ to 161,953 km³ in one step and collapsed the next.
   - Lake 13 (16.6 M km², 4.4 M km³) sat exactly on its saddle for 8 of 10 steps, evaporating 1.75–2.3× its inflow.
   - Knockout K2b, which re-resolves each child at its own balance (run with K1 below), kept sea 26 at ≤ 5 wet nodes instead of 1,713 and cut the step-to-step swing in global lake area from 2.4× to 1.4×.
2. **The frozen rule dries whole lakes, and its water goes unbooked.** `_water_balance` returns the dry level if *any* member of the lake's catchment, wet or dry, has an annual mean below 0 °C.
   - About 98% of the water created by resets came from resets this rule triggered on a merged parent.
   - Knockout K1, which judges freezing from wet cells only, removed lake 13's saddle pin (4.4 → 0.5–1.3 M km³) and cut global lake area by 62%.
   - K1 left a conservation hole: a lake dried because one wet cell froze leaves its warm sink with inflow and no lake, and that water (up to 181,000 m³/s) is booked nowhere.
3. **Merging is judged from last step's water, and perched lakes fool it.** `lakes._prev_level` returns the highest level of any wet member. A perched lake upstream can therefore make a parent look already merged; one cell 0.44 m above a 602.83 m saddle did so for sea 26. The same test chooses governing rims in `lakes.governing_bodies`, so it also redirects sinks' spill targets: river 307's 1,669 km sink-to-rim jump. Merging on last step's water is also #505's lag against routing, which already unites basins within the step (`hydrology._spill_within_step`).
4. **Rivers freeze on the annual mean and their water is deleted.** A river reaching a land sink whose annual mean is below 0 °C becomes river ice at `hydrology.WATER_TO_ICE`, which is 1/50,000 of its volume in the glacier model's units. On the evidence save that was 2.5–48% of land-to-ocean discharge per step, mostly one river of ~150,000 m³/s at −0.9 °C.
   - Under K2b, the loss fell ~75% and that sink dropped out of the top list. That is *consistent with* item 1's false flood having fed it, but no probe traced the water; G0 checks it.
   - The rule also contradicts the model's own principle #3, which uses exactly this case (a sub-zero annual mean doesn't mean no summer melt) as its example.
5. **Lakes read the wrong temperature.** Lake cells use the node's air temperature, read at the **bed** elevation:
   - a cell whose bed is below sea level reads the *ocean* temperature (`erosion.py`: `temperature = where(elevation <= sea_level, ocean_temperature_c, air_temperature_c)`), even inside an enclosed lake;
   - the lapse rate cools by the bed elevation clipped at 0, not by the water surface, so lake 13's shallower cells read about 3 °C too warm.

   Climate's own grid already uses a connectivity-aware ocean mask (`climate._sample_elevation_and_crust`); the problem is only the erosion-side temperature pick.
6. **Rain on a lake is counted as land runoff.** `compute_hydrology` computes Budyko runoff on every cell, wet lake cells included (`runoff_m3 = breaching.runoff_mm(liquid_precip, evaporation_mm) / 1000 * areas`). Rain on a lake surface therefore first loses land evapotranspiration, and the lake balance then evaporates over the same area again. Water crossing a lake also takes the per-node in-transit loss. Unmeasured; G0 measures it. For a lake of lake 13's size (~16.6 M km²) the missing direct precipitation could be thousands of km³/yr against an inflow of ~12,000 km³/yr.
7. **A depth cap stands in for the balance.** Once item 2 is fixed, lake 13 rests at the sea-tier depth cap (floor + `SEA_MAX_DEPTH_M`, 3,500 m) on its largest steps, with inflow (3,100–10,300 km³/yr) above its evaporation at the cap (2,200–2,800 km³/yr). The cap was added to stop a runaway basin when nothing else bounded it (`lakes.py`, "Lake-vs-sea tiers"). Whether it is still needed once the balance is right is open.

### What the evidence does not show

- Whether these mechanisms are as large on other seeds, ages and densities (G0).
- How far #520 part 1 (connected-water merge test) and #521 (no split reset), both standalone fixes, already remove items 1 and 3 by themselves.
- How large item 6 is anywhere.
- Whether evaporation's temperature-only rate matters to lake size. It runs 15–35% low for continental deserts per #519's literature check, but that belongs to #358 and #118, not here.
- Whether the two large-lake flickers in the #538 held-out saves are caused by missing groundwater storage, wind-enhanced evaporation, vertical heat storage, sediment placement, or integration/topology timing. The #675 follow-up measured the existing evaporation and sediment scales and compared dynamic 10, 25 and 100 kyr schedules; it did not add a physical mechanism or change the acceptance threshold (see [`analysis/issue675/README.md`](../analysis/issue675/README.md)).

### Follow-up diagnosis from #675

The #538 held-out saves show two very different geometries at the detector's reversal:

| Seed and event | Pre-event connected component | Geometry and temperature screen | Event |
|---|---|---|---|
| 373936961, step 5 at 244.8 Myr | 136,570 km³ liquid over 185,010 km²; mean depth 738.2 m, maximum 1,706.7 m; 156.6 km³ lake ice (water equivalent); area-weighted modeled surface-temperature proxy −1.86 °C (p10/p50/p90 −22.53/9.52/13.22 °C) | Widest cell-centre transect about 313 km, with 7 samples at density 4 | Component fell to 4,120 km³ liquid, a 132,450 km³ decrease |
| 997271774, step 7–8 at 131.1–131.3 Myr | 147.8 km³ over 61,995 km²; mean depth 2.38 m, maximum 4.10 m; modeled surface-temperature proxy 19.05 °C (p10/p50/p90 18.91/19.05/19.18 °C) | Widest cell-centre transect about 174 km, with 4 samples at density 4 | Fell to 59.62 km³ in step 7 (−88.15 km³), then rebounded to 157.59 km³ in step 8 |

These transects show how sparsely the model grid samples the cross-lake geometry; they do not resolve a smooth bathymetric profile. The saved-world model has a water-surface-temperature proxy, but no deep-water temperature state. Whole-world ledger closure in #673 rules out an unbooked global lake-water loss for these replays; it does not identify whether the component's water moved to another basin, crossed an outlet, or left the geographically matched component through a connectivity change. The local saves are checksum-registered but are not available for download; the quoted measurements and reproduction details are in [`analysis/issue675/README.md`](../analysis/issue675/README.md).

The shallow event is strongly schedule-dependent in a **dynamic** replay. At 131.2–131.3 Myr the same eight geographic anchors lie in a 1.45–1.51 million km² component under 10 and 25 kyr schedules, while the registered 100 kyr trajectory has only 31,119–65,796 km² around them and alternates wet and dry. The finer schedules also contain abrupt earlier expansions: the 10 kyr run grows by 47,904 km³ in one step at 130.54–130.55 Myr, and the 25 kyr run changes connectivity repeatedly around 130.6–130.725 Myr. Thus shorter steps do not simply smooth the event. The standard large-body detector finds zero reversals in the 10 kyr run, but misses this sustained geographic connectivity change. Ledger residuals are ≤5.58×10⁻¹³ in the 10 and 25 kyr runs, within the 0.1% bound. Since climate, erosion, routing and topology all evolve at each schedule, these runs establish path and event-timing sensitivity, not timestep convergence of the lake solver or a single-process cause.

The measured mechanism scales rank testable hypotheses; none is a causal result:

- The existing temperature-only evaporation proxy is about 1.445 m/year over the shallow footprint and 0.817 m/year over the deep footprint. Dividing the pre-event liquid volume by area and this rate gives no-inflow turnover screens of about 1.33 and 904 years, respectively. These are not drought predictions because catchment inflow continues. Wind, humidity, radiation and vertical water-column temperature are absent from the saved-world forcing. Lake Mead's measured 2.074 and 1.881 m/year annual evaporation was correlated monthly with wind speed, vapor-pressure difference and atmospheric stability; stored water-column heat delayed peak evaporation and sustained fall and nighttime losses [[Moreo & Swancar, 2013]](#ref-moreo2013). Global estimates likewise include a lake heat-storage term [[Zhao et al., 2022]](#ref-zhao2022). These sources motivate measuring the missing terms but do not calibrate these synthetic lakes: Lake Mead is a site-specific warm-arid reservoir, and Lake Winnipeg is a cold-climate, smaller analog.
- Prior-step deposition over equal-area event footprints is about 56 km³ on the shallow case (57% of its 97.97 km³ rebound) and 416 km³ on the deep case (0.3% of its 132,450 km³ decrease). The shallow one-hop spill target accumulated about 0.041 m of silt in that step, and the target after expansion was a different cell. This does not measure sill aggradation or sediment-only causality. A rate of 0.02 mm/year would accumulate 2 m in 100 kyr if maintained at the relevant connector; the test must measure local deposition against the actual sill/floor clearance, rather than assume that a basin-wide rate raises the outlet.
- Capacity-only shallow-basin screens give 18.6 km³ for a 1 m sediment layer at 0.3 porosity, 248 km³ for 10 m at 0.4 porosity, and about 32 km³ for a 60 km nearshore band with 3 m head change and specific yield 0.2. These scenarios span a useful test range, but are not estimates of empty, connected storage or of the volume exchanged during the event. No porewater or shoreline-aquifer exchange was implemented. A future reservoir must store only the water actually infiltrated and return it with equal-and-opposite ledger transfers.
- Lake Winnipeg is a cold-climate, smaller analog than either modeled lake. It can inform a climate/geometry comparison, but its evaporation and stratification parameters should not be transferred directly to these warm, much larger synthetic basins.

**Interpretation.** The dynamic replays establish schedule-dependent topology, but they do not isolate a cause. The mechanism tests therefore need the explicit follow-up gates G6–G8 below: close event attribution first, isolate schedule/topology under prescribed forcing next, then screen one physical process at a time. The two #538 saves are used only for read-only event attribution; candidate ranges are selected on disjoint diagnosis cases. Dynamic 10/25/100 kyr replay results remain diagnostic, not causal tests.

The existing **Flicker** threshold below is unchanged and remains open: #673 failed it on the two held-out saves, and #675 neither passes nor revises it. Only the project owner may change the threshold.

**Principles violated:** #6 (account for material through every process: water is created by resets and deleted at frozen sinks and frozen lakes, and rain on lakes is mis-booked); #3 (model fast processes in aggregate over long steps: annual-mean freezing and last-step merging); #5 (consistent units: water becomes ice at 1/50,000).

## Proposed solutions

### What happens in nature

- **Closed lakes adjust toward a balanced budget, at a rate set by their shape.** A closed lake's area adjusts until evaporation from its surface matches inflow, and it responds to climate change as a low-pass filter with an equilibrium response time that depends on the basin's morphology and climate. The publisher's abstract reports estimates of 1.5–350 years for its closed-lake examples and 2.9 days–2 years for open ones [[Mason et al., 1994]](#ref-mason1994). (Those figures were checked by the independent review on PR #523; the abstract wasn't accessible from the author's environment.)
  - These are adjustment times, not the turnover time `V / inflow`, and a steep-sided lake can turn over quickly yet adjust slowly (review on PR #523; G1).
  - On the evidence save, turnover times were 0.1–140 years for sea 26's basins and about 350–430 years for lake 13 at its inflated size. Its adjustment times weren't measured.
- **Water fills, spills and merges.** Runoff fills depressions, which overflow into their neighbours; when both a depression and its neighbour fill, they merge [[Barnes et al., 2021]](#ref-barnes2021). The *depression hierarchy* records how sub-depressions merge into meta-depressions at their spill points [[Barnes et al., 2020]](#ref-barnes2020). Mantle Bloom's `lakes.build_lake_hierarchy` is this structure on a k-NN graph. Fill–Spill–Merge is mass-conserving and keeps water held within a subtree separate from the excess passed on, but it has no evaporation; adding it is this RFC's own extension.
- **Lake ice is seasonal and set mostly by air temperature.** Observations at 513 Northern Hemisphere lakes identify air temperature, lake depth, elevation and shoreline complexity as key controls on ice cover, and an estimated 14,800 lakes have intermittent winter ice today [[Sharma et al., 2019]](#ref-sharma2019). The Global Lake and River Ice Phenology Database holds freeze and breakup dates for 865 lakes and rivers [[Benson et al., 2000]](#ref-benson2000). Ice cover is a fraction of the year, not an annual-mean switch.
- **Ice cover doesn't stop a lake's water budget.** Lake Untersee in Antarctica has been continuously ice-covered for at least 500 years. Yet it holds a 169 m deep basin, receives its only liquid water from the adjoining glacier, and loses water only by sublimation of its ice cover, with about 60 cm/yr of summer ablation [[Weisleitner et al., 2019]](#ref-weisleitner2019). Ice cover, liquid storage and phase change therefore need separate rules.
- **Ice grows by conduction, so its thickness is bounded.** Ice growth is classically modelled as one-dimensional vertical heat transfer, which gives analytic first-order approximations of ice thickness [[Leppäranta, 1993]](#ref-leppaesranta1993).
- **Lake surface temperature can be modelled from air temperature.** Air2Water relates the surface-layer temperature of lakes to air temperature alone, with a lumped heat exchange with the atmosphere and the deeper water [[Piccolroaz et al., 2013]](#ref-piccolroaz2013).
- **Evaporation depends on energy, the vapour-pressure deficit and heat storage.** A Penman-based global estimate, with a lake heat-storage term, puts long-term lake evaporation at 1500 ± 150 km³/yr, excluding the Caspian Sea [[Zhao et al., 2022]](#ref-zhao2022). The physical rate is out of scope here (#358, #118), but the lake solve takes evaporation from a function that a Penman rate can later replace.

### Alternatives

**A. Patch the current resolver** (#520 part 1, #521 and a wet-cells-only frozen rule).
- *For:* small, and the knockouts show it removes the largest artifacts.
- *Against:* merges still follow last step's water (#505), and freezing stays all-or-nothing on the annual mean, with frozen water unbooked. Rain on lakes stays mis-booked, and nothing audits conservation.
- These patches are worth landing now; #520 and #521 do. They are not the end state.

**B. In-step fill-spill-merge with evaporation, conservative step budgets, seasonal ice cover with a bounded ice cap, and a water-surface temperature (recommended).**
- Resolve all lakes each step from *this* step's water in one pass over the hierarchy, with no continuity pins and no "already merged" test.
- Close every body's water budget over the step.
- Model ice cover, liquid storage and phase change separately.
- Read temperature at the water surface.

**C. Prognostic lake storage.** Carry each lake's volume and integrate `dV/dt = inflow − evaporation − overflow` through the step with substeps.
- *For:* handles lakes whose adjustment time approaches the step.
- *Against:* more state and cost.
- **Kept as the per-body fallback** for bodies G1 or the allocation bounds (D1) flag.

**D. A full lake thermal model** (FLake-style mixed layer and ice).
- Too data-hungry for this model's inputs (no physical humidity, wind or radiation yet; #358, #118) and too costly at 10⁴–10⁵ lakes per world.

**Recommendation: B, with C applied per body where equilibrium fails.** B keeps the equilibrium premise where G1 supports it, and removes every mechanism in the Problem by construction rather than patch by patch. It makes routing and lakes one solve, and closes the water budget with an auditable ledger. #520 and #521 are interim fixes that B supersedes.

## Design

### D1. In-step balanced fill-spill-merge with a conservative step budget

**Sources.** Each cell's supply in m³/yr is split by its wet fraction `w` (the hypsometric wet fraction of last step's lake surface, a one-step lag like D3's):
- **wet part:** liquid precipitation `P × w × area`. The lake's evaporation is the only evaporation there, and there's no in-transit loss on wet cells;
- **dry part:** Budyko runoff `R(P, PET) × (1 − w) × area`, as today;
- **plus** glacier meltwater, as today.

Runoff is routed to sinks by `route_downstream` with the existing transit losses on dry cells. Each leaf's settled water is its **external supply**.

**The endpoint solve.** For each body `d` (a leaf, or a parent once both its children are full), its external supply `S_d` is:
- the external supply of all its member leaves, plus
- spills that enter `d` from **outside** its subtree.

Transfers between `d`'s descendants cancel, so no volume is counted twice (review on PR #523).
- If `d`'s evaporation filled to its rim, `E_d(rim)`, is at least `S_d`, `d` sits at the balanced level `h_d` where `E_d(h_d) = S_d`, found by bisection (`lakes.balanced_level`), and spills nothing.
- Otherwise `d` is full: it evaporates `E_d(rim)` and spills `O*_d = S_d − E_d(rim)` through its outlet.
- **Where a spill goes:** to the sibling subtree's leaf that owns `outlet_target_idx` in the *hierarchy's own* catchment assignment (pure steepest descent, `build_lake_hierarchy` phase 1), as Barnes et al.'s geolinks do; from a root, to the ocean. The river carrying it is still routed along `flow_target` for discharge and erosion, but only the hierarchy decides which basin receives the water. Any path whose `flow_target` reaches a different sink is counted as a diagnostic, and must not change water amounts.
- **Merging:** a parent floods as one body only when both children are full. Its balance then replaces its descendants'. The children's spills into each other are internal and excluded from `S`.
- **No merge decision depends on last step's state.** `_prev_level`, the `already_merged` test, the merge pin, the split reset and its tier-capped variant are removed. This removes #505's lag and the perched-lake false merge. The solve becomes the single source of routing's spill decisions (`should_spill`, `filled_elevation`, `spill_target`), replacing `governing_bodies`' use of last step's water.

  Three explicit lags remain: the wet fraction for sources above, D3's water-surface elevation, and `_classify_tier`'s sea-tier test, which reads last step's wet extent (an input to climate's humidity ceiling and to rendering).
- **Processing order:** leaves before parents. A spill into a sibling subtree raises that subtree's supply, and it is re-solved. Since `E_d(h)` is non-decreasing in `h` (`balanced_level` relies on this; `SHALLOW_EVAPORATION_BOOST < 1` keeps it so), adding water never lowers a level, and the endpoint should be unique and independent of order. G2 tests this against an **independent** reference: a time-stepped integration (alternative C) run to steady state, and nested and asymmetric hand cases checked against total source supply. A fixed-point solver of the same rules would reproduce the same rule errors.

**The conservative step budget** (review on PR #523). The endpoint gives each body's equilibrium level and volume `V_new`, but the fluxes booked for the step must close the budget:

`V_new − V_old = (S̄ − Ē − Ō) × Δt`, with `V_old = Σ lake_depth × area` from the persisted state.

Here `S̄`, `Ē` and `Ō` are step-averaged supply, evaporation (including D2's sublimation) and overflow, and `S̄` includes upstream bodies' `Ō`. Write `R = (V_new − V_old) / Δt` and let `E*`, `O*` be the endpoint rates.

- **Growing (`R > 0`):** filling water comes first from what would otherwise overflow: `Ō = max(O* − R, 0)`, then `Ē = S̄ − Ō − R`.
- **Shrinking (`R < 0`):** released water `−R` leaves by the outlet if the body lost its rim or drains at its floor (`Ō = O* − R`, a real flood that adds to downstream discharge and erosion). Otherwise it evaporates (`Ē = E* − R`).
- **Bounds:** `Ē ≥ 0`, `Ō ≥ 0`, and `Ē ≤ max(E(h_old), E(h_new))`, since evaporation over the step can't exceed the larger surface's rate.

  A body that violates a bound couldn't reach equilibrium within the step and is integrated with alternative C instead. The flagged share of volume is a G1 metric.
- **Downstream:** each spill passed on is `Ō`, so receivers' budgets close too. The endpoint solve uses the endpoint rates `O*`; the step budget uses `Ō`. Each body's budget closes exactly, and transfers match sender to receiver.

**Ledger** (new, per step, m³/yr; every reservoir disjoint):
- **sources:** precipitation on wet area, routed runoff and meltwater into basins;
- **sinks and reservoirs:**
  - evaporation and sublimation;
  - overflow to the ocean;
  - liquid storage change;
  - ice-cap storage change (D2), in water equivalent;
- **internal transfers**, outside the balance: spills between basins.

The residual should be zero to round-off; it's an acceptance criterion.

**Boundary:** this ledger covers the represented surface-water, lake and ice reservoirs. It
does not include subsurface storage or groundwater exchange, so ledger closure is not a claim
of whole-land water-cycle closure.

**Silt** (#363) is deposited once per body, from its step-averaged supply `S̄`, including the spill-through of a spilling lake (as today), over its final wet area, with `SILT_ACCUMULATION_COEFFICIENT`. A body whose resolution is superseded by its parent's deposits nothing, since its silt is part of the parent's.

**Bounds:** levels stay within [dry level, rim]; every booked flux is ≥ 0; ice-cap thickness stays within [0, water depth].

**Timestep:** rates are annual and the step budget is exact for any `Δt`. Equilibrium validity is gated by G1, and per-body fallback to C handles bodies that fail it.

### D2. Ice cover, liquid water under ice, and a bounded ice cap

Following the review on PR #523, ice cover, liquid storage and routing, and phase change are three separate rules.

**(a) Ice cover and open-water evaporation.**
- Each wet cell's ice-free fraction of the year is `f_free = φ / π`, with `φ = arccos(−(T_w − T_ice) / A_w)`. This uses the same sinusoidal annual cycle the glacier model already uses (`hydrology.seasonal_positive_degrees`), with `T_w` and `A_w` the water-surface mean and amplitude (D3). It is 1 when the whole year is above `T_ice` and 0 when the whole year is below.
- **Evaporation:** `e_open(T_open) × f_free` (today `e_open` is `breaching.open_water_evaporation_mm`, later Penman, #358). `T_open` is the mean water temperature over the ice-free part of the year, in closed form from the sinusoid.
- **Sublimation:** `s_ice × (1 − f_free)` from the ice-covered part.
- **Calibration shift.** `300 + 60·T` was fitted to *annual-mean* temperatures that already include ice seasons, so this changes the calibration. A first check on Lake Superior (rough T ≈ 4 °C, A ≈ 13 °C) gives f_free ≈ 0.60, T_open ≈ 10.6 °C and e ≈ 560 mm/yr, against ~557 observed and 540 from today's formula. An acceptance check covers #519's whole reference-lake table.

**(b) Liquid storage and routing continue under ice.**
- Ice cover doesn't stop inflow, storage or spill. A lake with `f_free = 0` everywhere still holds liquid water and balances its supply against sublimation and overflow, as Lake Untersee does [[Weisleitner et al., 2019]](#ref-weisleitner2019).
- Rivers carry their water to the sink they reach whatever the cell temperature. Water arriving at a cold sink forms or feeds a lake there, so liquid water is never deleted.
- **This replaces:**
  - the annual-mean `water_frozen` / `lake_frozen` tests and `frozen_from_lake` in `hydrology.compute_hydrology`;
  - the any-member `is_frozen[node.members].any()` in `lakes._resolve`;
  - the frozen-river-ice conversion of liquid water at `WATER_TO_ICE`.
- **Unchanged:** precipitation still falls as snow on frozen or ice-capped nodes (principle #8); snow onto a lake's ice cap joins the glacier accumulation as today.
  - *Note (RFC 0015, if accepted):* RFC 0015 D1 replaces the frozen-node snow switch with a seasonal snow fraction, and its D4 gives snow on a lake's wet fraction to the lake (onto `lake_ice_m` while frozen), not to glacier accumulation. This RFC stops writing `IceLedger.frozen_river_ice_m3`, `frozen_river_capped_m3` and `frozen_lake_m3`; RFC 0015 D5 deletes them. See [RFC 0015](0015-frozen-water-in-physical-units.md#d4-lakes-that-freeze-and-lakes-under-ice-with-rfc-0014).

**(c) A bounded ice cap for the non-permafrost lake scope.**
- The owner-directed rule assumes modeled lakes are non-permafrost and caps physical ice at 2 m. This encodes the owner's premise that surface ice on non-permafrost lakes rarely exceeds 2 m as a deterministic model ceiling; it is not a universal observed maximum and does not predict permafrost.
- Each wet cell's stored cap remains water equivalent `H_w`; convert physical ice thickness as `h_i = H_w ρ_w/ρ_i`. The water-equivalent bound is `H_w,max = min(2 m × ρ_i/ρ_w, local water depth)`, or at most 1.834 m water equivalent for 2 m of ice at ρ_i = 917 and ρ_w = 1000 kg m⁻³. Equivalently, `h_i,max = H_w,max × ρ_w/ρ_i`. Retain the existing seasonal cover and below-cap growth/loss treatment, but clamp at this bound.
- At depths ≤2 m, ice may reach the bed; report this as shallow bedfast area separately. Bottomfast ice is observed in shallow tundra lakes [[Hirose et al., 2008]](#ref-hirose2008). At depths >2 m, bedfast area must remain ≤1% of deep-water lake area. The threshold does not apply to shallow water. Mantle Bloom has no permafrost field or ground-temperature evolution, so no permafrost classification is implied.
- **Bounds and ledger:** `0 ≤ h_i ≤ min(2 m, water depth × ρ_w/ρ_i)` and `0 ≤ H_w ≤ water depth`. Freezing transfers liquid water into the cap; existing loss transfers remain as currently modeled. Cap water equivalent must stay within the lake-water budget. At the bed limit, additional growth stops and inflow remains available to the liquid lake budget and routing.
- **Design boundary:** no snow, sediment, water-column, radiation or atmospheric energy states are added for this rule. Detailed models such as CLIMo, CSLM and CLM4-LISSS represent those processes [[Brown & Duguay, 2012]](#ref-brownduguay2012) [[MacKay, 2019]](#ref-mackay2019) [[Subin et al., 2012]](#ref-subin2012), and a shallow-lake CSLM evaluation shows that those added processes still leave measurable ice and layer-partition errors [[Kallel et al., 2024]](#ref-kallel2024). They remain useful background in the issue #615 comparison, but are out of scope for the owner-directed cap. The hard ceiling is the selected empirical simplification; perennial ice such as Lake Untersee is outside this scope [[Steel et al., 2015]](#ref-steel2015).
- **G5 gate:** the #615 postprocessing sensitivity is not a full simulation result. The bounded offline prototype in #541 demonstrates the physical cap, water-equivalent bound and cap-volume plateau on the diagnosis saves, with shallow and deep bedfast shares reported. This clears the bounded diagnosis gate; #632 still must validate the production implementation and combined liquid/ice ledger. Keep the three held-out saves reserved for acceptance.
- **Glaciers:** the cap stays apart from `glacier_depth` until glacier mass balance is physical (#368). After that, a cap frozen to the bed can be handed to the glacier model by mass.

### D3. Water-surface temperature

- **Elevation.** Lake cells take their temperature at the **water-surface elevation**: last step's lake level, a one-step lag that climate already tolerates.
  - Enclosed lake cells use the air temperature even where the bed is below sea level. Only connected ocean cells use the ocean temperature, through the hydrology cache's `is_ocean`, as climate already does.
  - The lapse rate cools by the water-surface elevation.
- **Mean and amplitude.** `T_w` starts as the air temperature at the water surface. The seasonal amplitude can be damped with mean lake depth, `A_w = A × max(1 − D / D_ref, a_min)`, representing heat storage, in the spirit of Air2Water's exchange with deeper water [[Piccolroaz et al., 2013]](#ref-piccolroaz2013).
- **Gate.** G3 fits this before it's used; if the fit doesn't beat `A_w = A`, damping is dropped.
- **Climate.** This is the water temperature #358 item 10 asks for. Climate's own lake moisture source is unchanged here; one shared evaporation flux is #358 item 11.

### D4. The depth caps

After D1–D3, measure how often `LAKE_MAX_DEPTH_M` and `SEA_MAX_DEPTH_M` bind, on which lakes, and with what supply surplus (G4). The options were to keep both caps as model limits, retain them only as numerical guards with a warning (like `GLACIER_MAX_DEPTH_M`), or remove them.

**Owner decision after G4:** retain both only as numerical guards with a warning; they are not validated physical lake-depth limits. When either tier cap routes overflow, the solver emits one aggregate warning for the solve with affected-body counts and gross routed overflow by tier. Water a cap holds back remains booked through the lake hierarchy and ledger.

### Implementation order

At acceptance (2026-10-23), the owner scoped implementation first to what [RFC 0015](0015-frozen-water-in-physical-units.md) needs. RFC 0015 needs D2(b) and D2(c) as prerequisites, and D2(a)'s `f_free` and seasonal helpers.

- **Phase 1:** D2(b), then D2(a) after G3, then D2(c) after G5, with the seed sets chosen before G5.
- **Phase 2**, deferred until the owner reopens it: G0, G1, G2 and G4, with D1, D3 and D4.

G0's outcome decides only D1's scope, so it runs before D1 rather than before D2. Until D3 lands, D2 reads today's lake temperature. The tracking issue (#538) holds the checklist.

### Diagnosis gates

**Seed sets** (chosen and recorded on the tracking issue before G0, disjoint):
- a **diagnosis set** of ≥ 3 mature saves for G0, G1, G4 and G5;
- a **held-out acceptance set** of ≥ 3 other mature saves, never used for gates, knob choices or threshold revisions.

Neither set includes the development save, seed 308161711. All runs use the ice-age control at "never".

- **G0: baseline after #520 and #521.** Once those land, on the diagnosis set at 100 kyr steps, measure:
  - reset-created water;
  - false merges;
  - whole-lake freezing;
  - frozen-sink water loss, and whether its largest sinks are fed by spills from basins that flood and empty (Problem item 4);
  - the size of the rain-on-lakes mis-booking (item 6);
  - flicker reversals (defined in the acceptance criteria);
  - the balance residual of non-spilling lakes.

  *Outcome:* if the interim fixes already meet D1's acceptance thresholds and item 6 is below 1% of lake supply on every save, D1 shrinks to the ledger, the precipitation fix and the #505 unification, and the RFC is revised. Otherwise D1 proceeds as designed.
- **G1: where is equilibrium valid?** On the diagnosis set at 10 kyr and 100 kyr steps, measure per body:
  - the **adjustment time** `τ_relax = (dV/dh) / (dF/dh)` at the endpoint level, where `F` is total water loss (evaporation, sublimation and the overflow's sensitivity to level). A flat `F(h)` gives `τ_relax = ∞`;
  - the **fill time** `τ_fill = V_new / S̄` for bodies forming or growing from dry;
  - a **direct transient check**: integrate alternative C on frozen terrain and forcing for every body above 10⁴ km², plus a random sample of smaller ones, and compare the end-of-step volume with the endpoint `V_new`.

  *Pass:* ≥ 99% of lake volume has `max(τ_relax, τ_fill) ≤ 0.1 × Δt` at 10 kyr, **and** the transient-to-endpoint volume difference is ≤ 1% for ≥ 99% of lake volume. Bodies that fail use C (alternative C becomes per-body, not global). Report the flagged share and its runtime cost.
- **G2: is the endpoint solve right?** On random hierarchies (unit-test scale) and on the diagnosis saves' real hierarchies:
  - compare the endpoint solve with alternative C run to steady state (agreement ≤ 0.1% of volume);
  - check that total held plus evaporated plus spilled water equals total source supply on nested and asymmetric hand cases.

  If either fails, find the cause before implementation.
- **G3: ice cover.**
  - **(a)** From the Benson et al. database, take ≥ 50 lakes with long records. Compare observed ice duration with `1 − f_free` from D2 and D3, using each lake's observed mean and amplitude of air temperature (ERA5 monthly 2 m temperature climatology at the lake). Run it with and without depth damping. Adopt damping only if it lowers the median absolute error; target ≤ 30 days.
  - **(b)** At the same lakes' latitudes and continentality, report the model's synthesized amplitude (`biomes.seasonal_temp_amplitude`) against the observed one, since (a) validates the relation and not the model's input.

  **Result (issue #540; report:** [`analysis/issue540/README.md`](../analysis/issue540/README.md)**):** The 1990–2018 sample contains 54 lakes with at least 20 observed years. The undamped proxy's median absolute error is 10.61 days, meeting the ≤30-day target. Depth damping is not adopted: the full-sample fit improves by only 0.10 day and its leave-one-lake-out median error is worse (11.28 days vs. 10.61 days undamped). The model's synthesized seasonal amplitude has 2.08 °C median absolute error and +0.75 °C median bias. The evaporation calibration check passes across all nine #519 reference lakes; its largest increase in relative absolute error over the published observation-range endpoints is 3.91%, below the 10% allowance. G3 passes. The lake sample, source provenance, per-lake results and limitations are in the linked report and its compact JSON. G3 validates cover duration only, not D2(b) or D2(c).
- **G4: depth caps.** After D1–D3 on the diagnosis set, report how often each cap binds, on which lakes and with what surplus. *Outcome:* retain both as numerical guards with warnings when they route overflow; see the result below.

**G4 result and D4 decision (issue #553; report:** [`analysis/issue553/README.md`](../analysis/issue553/README.md)**):** On the three checksum-matched diagnosis saves, each replayed for ten 100 kyr steps on the merged C-default solver, the lake cap bound for 34, 2 and 13 body-steps (seeds 910211954, 635056485 and 317739385); the sea cap bound for 12, 0 and 1 body-steps. The largest gross per-body overflow surpluses were `4.5591e20` and `2.2244e20 m³` across the lake and sea tiers on seed 910211954, `5.6545e17 m³` for the lake tier on seed 635056485, and `1.0229e20` plus `2.8775e17 m³` for the lake and sea tiers on seed 317739385. See the report for binding-step frequencies, body areas, sink identities, and the surplus accounting. Both caps bind repeatedly on two saves; seed 635056485 has no sea-cap events. G4 is complete. **D4:** retain both ceilings only as numerical guards and warn when they route overflow; the thresholds are not validated physical limits. Overflow stays booked through the lake hierarchy and ledger.
- **G5: ice-cap budget.** On the diagnosis set over 10 Myr, report cap thickness, cap volume over time, global sublimation flux, and bedfast area split into lake depth ≤2 m and >2 m. Bedfast area divided by total lake area deeper than 2 m must be ≤1%; shallow bedfast area is reported separately. *Outcome:* adopt or revise the empirical cap rule after checking the full simulation and ledger.

**Result (issue #541; report:** [`analysis/issue541/README.md`](../analysis/issue541/README.md)**):** The three #538 diagnosis saves reach a volume plateau by 0.1 Myr under fixed climate and lake depth. Under the prior unbounded starting law, 5.09% and 18.45% of wet area freeze to the bed on seeds 910211954 and 317739385; seed 635056485 is at 0.54%. The prior law fails the revised depth-aware diagnosis because deep-water bedfast area is 2.344% and 3.777% of total wet area on the two failing seeds. Its baseline physical-ice thickness p50s were 4.50, 9.50 and 1.71 m respectively. This prior result is historical and does not test the 2 m cap.

**Updated G5 result (issue #541; report:** [`analysis/issue541/README.md`](../analysis/issue541/README.md)**):** The bounded offline prototype reran all nine sublimation/conduction settings on the three checksum-matched diagnosis saves over 10 Myr. Under the 1×/1× setting, bedfast area in water deeper than 2 m is 0% for all saves; shallow bedfast area is 2.570%, 0.537% and 11.231% of wet area (seeds 910211954, 635056485 and 317739385). All settings keep physical ice at or below 2 m and cap volume plateaus by 0.1 Myr. **The revised depth-aware G5 diagnostic passes.** This remains a fixed-climate, fixed-depth prototype, not production-code or dynamic ledger validation; the three held-out saves remain reserved.

- **G6: attribute the two #675 events before causal interpretation.** Use the exact #538 held-out save checksums and the registered 100 kyr schedule only to trace the seed-373936961 step-5 event and seed-997271774 step-7/8 event; do not tune parameters or change criteria from these runs. Starting before each event and continuing through the drop/reversal, match geographic basin lineages and their descendants. For each step, book liquid and ice-equivalent storage, external precipitation/runoff/melt supply, evaporation and sublimation, overflow to named receiving basins or ocean, and any connectivity-only reassignment. Count each internal transfer once; do not use gross hierarchy-flush totals. **Pass:** each event-step basin budget has normalized residual ≤0.1% of its external lake supply, and at least 99% of the measured component storage change is attributed to those terms or an explicitly measured connectivity reassignment. **If it fails:** improve lineage/flux instrumentation and stop; do not rank physical mechanisms from a component-volume difference. **If it passes:** use the signed terms to select which G8 screens are relevant; a closed budget alone does not establish physical plausibility.

  **G6 result (issue #675; report:** [`analysis/issue675/README.md`](../analysis/issue675/README.md)**):** The corrected plate-qualified lineage trace closes all three registered event steps. Normalized basin residuals are `8.37e-17`, `8.59e-18`, and `3.83e-18`; measured component storage change is assigned at >99% in each case after explicit connectivity reassignment. The deep step-5 loss is −105,533 km³, the shallow step-7 loss is −88.149 km³, and the shallow step-8 rebound is +97.968 km³. Destination-aware boundary records identify receiver lake IDs for all routed lake transfers and one deep-event export to a non-lake sink; none of these event traces records an ocean-target transfer. The corrected trace records boundary transfers without double-counting internal hierarchy flushes. **G6 passes.** This resolves event water accounting for these windows; it does not identify a physical mechanism or validate the flicker threshold.
- **G7: separate timestep/topology sensitivity from process response.** From a common event checkpoint, hold terrain, hierarchy/outlets, non-lake climate, and a prescribed basin-supply time series fixed. Compare the lake solve at 10, 25 and 100 kyr partitions, splitting the same piecewise-constant forcing intervals so every schedule receives the same integrated inputs. Record end-of-interval liquid-plus-ice-equivalent storage, anchor connectivity, component count, event timing, and ledger residual. **Pass:** matched-time storage agrees within 1% and anchor connectivity/event classification agrees across partitions, with residual ≤0.1% of supply. **If it fails:** treat the event as solver/timestep/topology-sensitive, fix or characterize that path before attributing it to a physical mechanism. **If it passes:** carry the fixed schedule into G8. This controlled result is separate from the earlier coupled-world dynamic replays.
- **G8: run one-factor mechanism screens on disjoint diagnosis cases.** Before running, record at least three mature diagnosis event windows that contain comparable large-lake reversals, the source commit and forcing, and parameter ranges; exclude the #538 held-out saves from range selection and tuning. For each case, use paired runs from the same checkpoint, the G7 schedule and identical non-target inputs. Run porewater storage and shoreline-aquifer exchange as separate tests with equal-and-opposite reservoir transfers; test sediment against the measured connector/sill clearance with a sediment-off control; test wind-enhanced evaporation and vertical heat storage separately and only when explicit meteorological and mixed-layer inputs exist. Record event ΔV, flooded area/connectivity, event timing, each process flux/store, ledger residual, and the range source. The capacity figures in the #675 report are illustrative screen points, not these physical parameter ranges.

  **Proposed materiality rule (design choice; owner confirmation required before it is used):** a mechanism is a material contributor when, within a source-supported range, it changes `|ΔV|` by ≥10% relative to its same-checkpoint paired baseline, or changes anchor connectivity/event classification for at least one 10 kyr interval, in at least two of the three diagnosis windows. A material result advances that mechanism to a separately scoped implementation/design decision and a fresh held-out validation; it does not prove that mechanism is the unique cause. **Negative result:** if every source-supported value stays below that response threshold in all diagnosis windows, deprioritize the mechanism as a primary cause. **Blocked result:** if physical ranges or required inputs are unavailable, defer the screen and report the missing evidence; do not call it falsified. If the event budget fails closure, return to G6 instead of interpreting the perturbation.

  Select a new held-out set, disjoint from diagnosis and the already-used #538 saves, before any future acceptance replay. Do not tune ranges, metrics, or the existing Flicker threshold on held-out results. Lake Winnipeg remains a cold, smaller comparison case and is not a coefficient source for these large synthetic basins.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| Maximum physical lake-ice thickness | 2 m, with water-equivalent cap also limited by local water depth | Owner-directed empirical assumption for the non-permafrost lake scope | No fit; verify the bound and report shallow/deep bedfast fractions on the G5 diagnosis set |
| Ice/water density ratio | 917/1000 | Physical densities; converts the physical cap to the water-equivalent ledger | None |
| Existing seasonal cover and below-cap growth/loss terms | Retain their current values | Existing D2(c) treatment | Not retuned by this design decision; G5 reports resulting cap volume and global sublimation |
| `SHALLOW_EVAPORATION_BOOST` | 0.3 (unchanged) | Existing (#297) | Out of scope |
| `LAKE_MAX_DEPTH_M`, `SEA_MAX_DEPTH_M` | Unchanged; numerical guards with warnings on binding | Existing thresholds, not physically calibrated | G4, owner decision |
| G1 thresholds | 0.1 × Δt, 1% | Design choice | Owner may revise |

### Cost and saves

- **Runtime.**
  - The endpoint solve visits each hierarchy node a small number of times, with one bisection per wet body: the same order as today's `_resolve`. It replaces `_resolve` and the fixed-point rounds of `_spill_within_step`. Report the number of spill rounds, since each re-route calls `route_downstream`.
  - The ice-free fraction, degree-days and cap update are closed-form per wet cell.
  - Alternative C runs only on flagged bodies; G1 reports their share and cost.
  - Target: ≤ +10% median step time at density 4.
- **Save compatibility.**
  - `lake_depth` is already persisted and gives `V_old`.
  - The ice-cap thickness `H` is new per-node state: a new `SURFACE_FIELDS` entry, default 0, extensive remap, with old saves loading at 0 (`docs/save-compatibility.md`, `unit_tests/test_persistence.py`).
  - Lagged inputs read the persisted `lake_depth`.

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/lakes.py`:
  - `_resolve`, `_water_balance`, `_prev_level`, `governing_bodies`, `compute_spill_routing` and `resolve_lakes` are replaced or rewritten around a new in-step solve and step budget (new, e.g. `solve_lakes`);
  - kept: `balanced_level`, `lake_evaporation_m3_per_yr` (gains `f_free` and sublimation), `CellHypsometry`, `build_lake_hierarchy` (its catchment assignment decides spill receivers), `_mark_floor_drains` and `_classify_tier` (lagged).
- `backend/app/hydroclimate/hydrology.py`:
  - `compute_hydrology`: the source split by wet fraction, no transit loss on wet cells, the order of routing, solve and spill routing; removes `water_frozen`, `lake_frozen`, `freezing_lake`, `frozen_from_lake` and the frozen-river-ice block for liquid water;
  - `_spill_within_step` and `SpillBasins`: folded into the solve;
  - `_compute_flow_direction`: freeze override removed for liquid routing;
  - `seasonal_positive_degrees`: reused, plus new siblings for the ice-free fraction, the open-season mean and freezing degree-days;
  - `WATER_TO_ICE` is no longer used for liquid water.
- `backend/app/hydroclimate/erosion.py`: the `temperature` assembly (~L1926: enclosed lake cells and water-surface lapse), and silt from the solve's `S̄`.
- `backend/app/hydroclimate/ice_sheet.py`: `IceLedger` gains the ice-cap reservoir and sublimation terms; `frozen_river_ice_m3` is retired for liquid water.
- `backend/app/surface_fields.py`: new `lake_ice_m` (ice-cap thickness, water equivalent).
- `backend/app/hydroclimate/breaching.py`: `open_water_evaporation_mm` is unchanged and called with `T_open`. `_overflows`' rim-level demand switches to the solve's evaporation and sublimation, so the gate and lakes agree.
- `backend/app/main.py`: `_lake_basin_summary` and `_smallest_lake_containing` read `is_spilling` and levels from the new solve.
- **Tests:**
  - `backend/unit_tests/test_lakes.py`:
    - merge, split, tier and frozen tests rewritten for the solve;
    - new: nested and asymmetric supply cases, a growing and a shrinking single basin with no plug term, perched lakes, spills delivered through the hierarchy when a channel crosses catchments, silt totals and the ice-cap bounds;
  - `test_hydrology.py`: the wet-fraction source split, rivers into cold sinks;
  - `test_persistence.py`: the new field;
  - `backend/stress_tests/`: multi-step replays.
- **Docs:** `docs/simulation-model/hydrology.md` (balanced lakes, lake tree, tiers, frozen lakes, rivers, the ledger), `docs/simulation-model/climate.md` (lake temperature), `docs/save-compatibility.md`.
- **Diagnostics (new):** a `bin/debug/` replay script for G0, G1, G4 and G5 metrics and the lake ledger, built from #519's probes (in the shared handoff directory), with a usage docstring.

## Assumptions and uncertainty

- **Equilibrium within a step** holds only where adjustment times are short. Turnover times on one save (0.1–430 years) don't establish adjustment times, and a steep-sided lake can turn over fast and adjust slowly. G1 measures this directly and C covers the rest.
- **A unique endpoint** under fill-spill-merge with evaporation is argued from `E(h)` being non-decreasing, not proven. G2 tests it against an independent integration.
- **Ice cover:**
  - A sinusoidal annual cycle with the model's synthesized amplitude stands in for real seasons, and air temperature stands in for water temperature, apart from optional depth damping.
  - Sharma et al. also name elevation and shoreline complexity as controls; neither is represented.
  - G3 measures the error of the relation and of the input separately.
- **The 2 m physical-ice ceiling is an owner-supplied empirical assumption** for modeled non-permafrost lakes. It does not resolve seasonal lake thermodynamics or represent permafrost. Treat the hard ceiling as a simplification of the owner's “rarely thicker than 2 m” premise; revisit it if the target scope includes permafrost lakes or observations contradict the bound.
- **Liquid water under ice** is assumed to keep its storage and spill; a frozen-to-bed cell is the only exception.
- **Evaporation is still the temperature-only rate.** It runs 15–35% low for continental deserts and ~70% low at high altitude (#519). Lake sizes will move again when #358 and #118 give a Penman rate; that interface is preserved.
- **Subsurface water storage is outside the implemented design.** Water can infiltrate soil, sediment and porous or fractured rock, where it may be stored and later return to surface waters; rates and storage depend on porosity, permeability and material type [[U.S. Geological Survey, n.d.]](#ref-usgs-infiltration). The model has no soil-moisture or groundwater reservoir, infiltration flux, or delayed groundwater return. The #675 values (1 m/0.3 and 10 m/0.4 pore-capacity cases; 60 km/3 m/0.2 shoreline case) are analyst-selected geometric scenarios, not literature-supported natural bounds or measurements of empty/exchangeable water. Before G8, source a plausible range and sensitivity plan or keep the test explicitly exploratory; any implemented transfers need disjoint storage ownership and equal-and-opposite ledger entries.
- **Evidence scope.** The original RFC's mechanism magnitudes came from one 121.5 Myr world. The #675 scale screens add two registered held-out events and dynamic schedules for one of those seeds, but do not establish causality, prevalence across worlds, or broad validation. The #538 saves are reserved from parameter selection; select separate diagnosis cases and a fresh held-out set for G8 and later validation.

## Acceptance criteria

- [x] **G0–G5 run and recorded** in this RFC, with the owner's D4 decision. G5's prior unbounded candidate failed; its bounded offline diagnosis passed with a cap-volume plateau, bounded water-equivalent cap and 0% deep-water bedfast area across the diagnosis saves. The production implementation and combined-ledger validation now pass on three diagnosis saves and all coded bounds pass on three reserved held-out saves (two had informative lake ice; the third had almost no lake water). See [#632 analysis](../analysis/issue632/README.md). Report shallow bedfast area separately. The #615 postprocessed sensitivity does not satisfy the production diagnosis criterion.
- [x] **G6 event attribution:** for the registered deep step-5 event and shallow steps 7–8, basin residuals are below 0.1% of supply and more than 99% of component storage change is assigned to fluxes or explicit connectivity reassignment. See [#675 analysis](../analysis/issue675/README.md). G6 passes for these registered windows; this does not establish a physical cause.
- [ ] **G7 fixed-forcing partition control:** the same prescribed event forcing across 10/25/100 kyr partitions yields ≤1% matched-time storage difference and the same anchor connectivity/event classification, with ≤0.1% ledger residual; otherwise characterize/fix the solver/timestep path before G8.
- [ ] **G8 mechanism screens:** ≥3 preregistered diagnosis event windows, disjoint from #538 and future held-out saves; source-supported ranges; one-factor paired runs; ledger closure. The proposed materiality rule is ≥10% change in `|ΔV|` relative to the paired baseline or a one-10-kyr-interval connectivity/event-class change in ≥2/3 windows. This is a design choice requiring owner confirmation before use. Record supported, deprioritized, and blocked mechanisms separately.
- [ ] **Unit tests:**
  - [ ] the endpoint solve matches alternative C at steady state (≤ 0.1% of volume) on random hierarchies;
  - [ ] nested and asymmetric supply cases hold, spill, and evaporate exactly their total source supply;
  - [ ] growing and shrinking single basins close their budget with no plug term;
  - [ ] a perched lake above a saddle doesn't merge its parent;
  - [ ] a merged lake whose balance falls below its saddle splits with each child at its own balance;
  - [ ] a spill crossing catchments via a channel is delivered to the hierarchy's receiving leaf;
  - [ ] the ice-free fraction is 0 below and 1 above the threshold range and continuous between;
  - [ ] a lake with `f_free = 0` keeps liquid water and balances its supply by sublimation and overflow;
  - [ ] physical ice thickness stays within 2 m and the water-equivalent cap stays within available water depth, closes the liquid/ice budget, and cap volume approaches a plateau under constant forcing (unit tests and fixed-forcing prototype; production replay verifies bounds and reports the evolving trajectory);
  - [ ] basal melt and surface melt return water to liquid; sublimation is the only D2(c) atmospheric mass loss;
  - [ ] water at a cold sink is never deleted;
  - [ ] an enclosed lake cell with its bed below sea level uses the air temperature at its water surface;
  - [ ] silt equals coefficient × supply × years.
- [ ] **Ledger closes:** |residual| ≤ 0.1% of lake supply on every step of a 10-step replay on each held-out save.
- [ ] **No water deleted:** liquid water lost outside the ledger is 0 on held-out saves.
- [ ] **No pinned lakes:** no non-spilling, uncapped body (not flagged to C) whose endpoint evaporation differs from its supply by more than 5% (proposed, owner to confirm).
- [ ] **Flicker:** per step, the fraction of bodies above 10⁴ km² whose volume change reverses sign on consecutive steps with |ΔV| > 50% of V is ≤ 5% (proposed; ice-age control "never"). Global volume max/min is reported but isn't pass/fail. **Still unresolved:** #673 failed this threshold on the two informative held-out saves; #675 found strong 10/25/100 kyr dynamic topology dependence, but did not revise or satisfy the gate.
- [ ] **Step robustness of the lake solve:** on frozen terrain and climate from each held-out save, solving with Δt = 10 kyr and 100 kyr gives the same endpoint volumes for unflagged bodies (≤ 0.1%) and volumes within 5% for flagged ones. Whole-model 100 × 10 kyr vs 10 × 100 kyr comparisons are reported, not pass/fail, because tectonics differs.
- [ ] **Resolution:** global lake area agrees within 20% between node densities 4 and 6 on one held-out seed (proposed).
- [ ] **Evaporation calibration:** D2's `e_open(T_open) × f_free` on #519's reference-lake table, using ERA5 1990–2019 means and amplitudes sampled at HydroLAKES v1 polygon centroids, doesn't worsen any lake's fit by more than 10% against today's formula. G3 reports this check passing across all nine lakes; the largest increase over published observation-range endpoints is 3.91%.
- [ ] **Side effects:** breached-pit counts and floor-drain toggles reported against `origin/main`; the #503 and #218 tests pass.
- [ ] **Runtime:** median step time at density 4 within +10% of `origin/main`.
- [ ] **Saves:** old saves load (`test_persistence.py`); `lake_ice_m` defaults to 0.
- [ ] **Docs** updated as listed above.

## Falsification criteria

- **G1 fails broadly:** more than 20% of lake volume is flagged to alternative C at 10 kyr steps. The equilibrium premise is wrong for this model's lakes. *Revise:* make C the default solve and the endpoint solve an optimisation for flagged-fast bodies.
- **G2 fails:** the endpoint solve disagrees with steady-state integration beyond 0.1% of volume, or depends on processing order. *Revise:* fix the rules (order, transfer accounting) or move to C.
- **G0 shows the interim fixes suffice:** #520 and #521 already meet the flicker, no-pinned-lakes and no-water-deleted criteria, and item 6 is under 1%, on every diagnosis save. *Revise:* withdraw D1's solve; keep the precipitation fix, the ledger, D2 and D3.
- **The 2 m cap implementation exceeds its physical bound, fails ledger closure, does not reach a cap-volume plateau, or leaves >1% of deep-water area bedfast:** *Revise:* fix the cap/ledger implementation or reconsider the owner-directed empirical bound. Shallow bedfast area is reported separately and does not trigger the deep-water threshold.
- **G3 gives no skill:** `1 − f_free` misses observed ice duration by more than 60 days median, even with observed inputs and damping. *Revise:* the sinusoid-and-air-temperature proxy is too crude; consider an ice-duration regression instead.
- **G6 does not close:** event attribution residual exceeds 0.1% of external lake supply or more than 1% of component storage change remains unassigned. *Consequence:* the lineage/accounting hypothesis is unresolved; improve instrumentation and do not interpret mechanism tests.
- **G7 remains partition-sensitive:** matched-time storage differs by more than 1% or the anchor connectivity/event class changes across partitions under the same prescribed forcing. *Consequence:* treat timestep/topology handling as an active numerical hypothesis and defer physical attribution until characterized.
- **G8 finds no material response:** all source-supported settings change `|ΔV|` by <10% and preserve anchor connectivity/event class across all diagnosis cases. *Consequence:* deprioritize that mechanism as the primary flicker explanation, without claiming the process is absent. If physical ranges or meteorological/mixed-layer inputs are unavailable, mark the candidate blocked/deferred rather than falsified.

## References

- <a id="ref-barnes2020"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2020). Computing water flow through complex landscapes – Part 2: Finding hierarchies in depressions and morphological segmentations. *Earth Surface Dynamics*, 8, 431–445. [https://doi.org/10.5194/esurf-8-431-2020](https://doi.org/10.5194/esurf-8-431-2020)
- <a id="ref-barnes2021"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2021). Computing water flow through complex landscapes – Part 3: Fill–Spill–Merge: flow routing in depression hierarchies. *Earth Surface Dynamics*, 9, 105–121. [https://doi.org/10.5194/esurf-9-105-2021](https://doi.org/10.5194/esurf-9-105-2021)
- <a id="ref-benson2000"></a>Benson, B., Magnuson, J., & Sharma, S. (2000). *Global Lake and River Ice Phenology Database, Version 1*. NSIDC (G01377). [https://nsidc.org/data/G01377/versions/1](https://nsidc.org/data/G01377/versions/1)
- <a id="ref-brownduguay2012"></a>Brown, L. C., & Duguay, C. R. (2012). Simulation of surface temperature and ice cover of large northern lakes with 1-D models: a comparison with MODIS satellite data and in situ measurements. *Tellus A*, 64. [https://doi.org/10.3402/tellusa.v64i0.17614](https://doi.org/10.3402/tellusa.v64i0.17614)
- <a id="ref-kallel2024"></a>Kallel, H., Nadeau, D. F., Bouchard, B., Thiboult, A., MacKay, M. D., & Anctil, F. (2024). Modeling seasonal ice and its impact on the thermal regime of a shallow boreal lake using the Canadian Small Lake Model. *Frontiers in Environmental Science*, 12. [https://doi.org/10.3389/fenvs.2024.1371108](https://doi.org/10.3389/fenvs.2024.1371108)
- <a id="ref-leppaesranta1993"></a>Leppäranta, M. (1993). A review of analytical models of sea-ice growth. *Atmosphere-Ocean*, 31(1), 123–138. [https://doi.org/10.1080/07055900.1993.9649465](https://doi.org/10.1080/07055900.1993.9649465)
- <a id="ref-mackay2019"></a>MacKay, M. D. (2019). Incorporating wind sheltering and sediment heat flux into 1-D models of small boreal lakes: a case study with the Canadian Small Lake Model V2.0. *Geoscientific Model Development*, 12, 3045–3054. [https://doi.org/10.5194/gmd-12-3045-2019](https://doi.org/10.5194/gmd-12-3045-2019)
- <a id="ref-hirose2008"></a>Hirose, T., Kapfer, M., Bennett, J., Cott, P., Manson, G., & Solomon, S. (2008). Bottomfast Ice Mapping and the Measurement of Ice Thickness on Tundra Lakes Using C-Band Synthetic Aperture Radar Remote Sensing. *Journal of the American Water Resources Association*, 44(2), 285–292. [https://doi.org/10.1111/j.1752-1688.2007.00161.x](https://doi.org/10.1111/j.1752-1688.2007.00161.x)
- <a id="ref-mason1994"></a>Mason, I. M., Guzkowska, M. A. J., Rapley, C. G., & Street-Perrott, F. A. (1994). The response of lake levels and areas to climatic change. *Climatic Change*, 27, 161–197. [https://doi.org/10.1007/BF01093590](https://doi.org/10.1007/BF01093590)
- <a id="ref-piccolroaz2013"></a>Piccolroaz, S., Toffolon, M., & Majone, B. (2013). A simple lumped model to convert air temperature into surface water temperature in lakes. *Hydrology and Earth System Sciences*, 17, 3323–3338. [https://doi.org/10.5194/hess-17-3323-2013](https://doi.org/10.5194/hess-17-3323-2013)
- <a id="ref-sharma2019"></a>Sharma, S., Blagrave, K., Magnuson, J. J., et al. (2019). Widespread loss of lake ice around the Northern Hemisphere in a warming world. *Nature Climate Change*, 9, 227–231. [https://doi.org/10.1038/s41558-018-0393-5](https://doi.org/10.1038/s41558-018-0393-5)
- <a id="ref-steel2015"></a>Steel, M., McKay, C. P., & Andersen, D. T. (2015). Modeling circulation and seasonal fluctuations in perennially ice-covered and ice-walled Lake Untersee, Antarctica. *Limnology and Oceanography*, 60(4), 1231–1245. [https://doi.org/10.1002/lno.10086](https://doi.org/10.1002/lno.10086)
- <a id="ref-subin2012"></a>Subin, Z. M., Riley, W. J., & Mironov, D. (2012). An improved lake model for climate simulations: Model structure, evaluation, and sensitivity analyses in CESM1. *Journal of Advances in Modeling Earth Systems*, 4, M02001. [https://doi.org/10.1029/2011MS000072](https://doi.org/10.1029/2011MS000072)
- <a id="ref-usgs-infiltration"></a>U.S. Geological Survey. (n.d.). *Infiltration and the Water Cycle*. [https://www.usgs.gov/water-science-school/science/infiltration-and-water-cycle](https://www.usgs.gov/water-science-school/science/infiltration-and-water-cycle)
- <a id="ref-weisleitner2019"></a>Weisleitner, K., Perras, A., Moissl-Eichinger, C., Andersen, D. T., & Sattler, B. (2019). Source environments of the microbiome in perennially ice-covered Lake Untersee, Antarctica. *Frontiers in Microbiology*, 10, 1019. [https://doi.org/10.3389/fmicb.2019.01019](https://doi.org/10.3389/fmicb.2019.01019)
- <a id="ref-zhao2022"></a>Zhao, G., Li, Y., Zhou, L., & Gao, H. (2022). Evaporative water loss of 1.42 million global lakes. *Nature Communications*, 13, 3686. [https://doi.org/10.1038/s41467-022-31125-6](https://doi.org/10.1038/s41467-022-31125-6)
- <a id="ref-moreo2013"></a>Moreo, M. T., & Swancar, A. (2013). Evaporation from Lake Mead, Nevada and Arizona, March 2010 through February 2012. *U.S. Geological Survey Scientific Investigations Report 2013–5229*. [https://doi.org/10.3133/sir20135229](https://doi.org/10.3133/sir20135229)
