# RFC 0014: Lakes conserve their water: in-step fill-spill-merge, seasonal freezing and water-surface temperature

- **Status:** Accepted (2026-10-23)
- **Date:** 2026-10-22 (revised 2026-10-23 after review)
- **Related:** #538 (tracking), #539 (D2b), #540 (gate G3), #541 (seed sets, gate G5); #519 (exploration and evidence), #520 (perched-lake false merges, river lines across lakes, floor-drain flips), #521 (split reset), #505 (lake merges lag routing), #503 and #297 (breaching, balanced lakes), #355 (water units), #363 (lake silt), #353 and #368 (glacier ice and its units), #358 items 10–11 (water-surface temperature, one evaporation flux), #118 (wind speed)

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

**(c) A bounded, energy-limited ice cap.**
- Each wet cell carries an ice-cap thickness `H_w` (water equivalent), persisted between steps. For heat transfer, use physical ice thickness `h_i = H_w ρ_w / ρ_i`; all water-ledger changes remain in water-equivalent volume.
- **Winter growth:** replace the bare-ice, one-sided Stefan approximation with the seasonal interface heat budget. The upward conductive flux is `q_c = (T_f − T_s) / (h_i/k_i + h_s/k_s)`, where `T_s` is the ice/snow surface temperature and `h_s/k_s` is snow's thermal resistance. The Stefan balance is `ρ_i L_f dh_i/dt = q_c − q_w`, where positive `q_w` is heat delivered from the lake water to the ice base. Snow insulation and basal heat therefore reduce the no-snow, zero-basal-flux Stefan growth, which is an upper bound rather than a central estimate [[Leppäranta, 1993]](#ref-leppaesranta1993).
- **Summer and seasonal loss:** keep surface melt and sublimation as distinct fluxes. Melt returns liquid water to the lake; only sublimation leaves the water system. The observed ~0.6 m summer ablation at Lake Untersee is compensated by winter basal refreezing and is not an annual sublimation rate [[Weisleitner et al., 2019]](#ref-weisleitner2019). Do not use it as `s_ice`. Ice-covered lake observations show that snow and water-to-ice heat flux vary substantially with season and lake conditions, so neither can be represented as an unqualified universal constant [[Huang et al., 2022]](#ref-huang2022).
- **Reduced forcing and gate:** the current climate has no incident shortwave radiation, snow-depth, humidity, wind, or under-ice water-temperature state from which to calculate `T_s`, `h_s`, `q_w`, or sublimation mechanistically. RFC 0015 D4 currently sends snow on a frozen lake directly into `lake_ice_m`; a snow-insulation term therefore requires either a distinct lake-snow reservoir and a coordinated D4 change, or an observed-data-screened reduced proxy. Keep snow resistance, lake-to-ice heat flux and sublimation as separate uncertain terms in the G5 prototype. Do not pick final values by fitting cap thickness alone. If the reduced inputs cannot predict seasonal thickness and bedfast occurrence on observed lakes, defer D2(c) or scope the missing thermal state before implementation.
- **Bounds:** integrate the annual forcing over its seasonal intervals (including partial-step overlaps); clip `H_w` to `[0, water depth]`. If `q_w > q_c`, basal melt returns ice mass to liquid rather than forcing negative thickness. At the bed limit, additional freezing stops and liquid inflow remains in the lake budget and routing.
- **Ledger:** freezing transfers water from liquid storage to the cap; basal or surface melt returns it. Sublimation alone is an external water loss. The cap is a disjoint, bounded reservoir booked by water-equivalent volume, never at `WATER_TO_ICE`.
- **Bedfast ice:** report the fraction and depth distribution of wet cells at the bed limit separately from cap thickness. Bottomfast ice is observed in shallow tundra lakes [[Hirose et al., 2008]](#ref-hirose2008), so the fraction is a diagnostic to compare with observations, not evidence by itself that water was lost. Preserve the current >1% G5 falsification threshold for this rerun unless the owner explicitly revises it.
- **Glaciers:** the cap stays apart from `glacier_depth` until glacier mass balance is physical (#368). After that, a cap frozen to the bed can be handed to the glacier model by mass.

### D3. Water-surface temperature

- **Elevation.** Lake cells take their temperature at the **water-surface elevation**: last step's lake level, a one-step lag that climate already tolerates.
  - Enclosed lake cells use the air temperature even where the bed is below sea level. Only connected ocean cells use the ocean temperature, through the hydrology cache's `is_ocean`, as climate already does.
  - The lapse rate cools by the water-surface elevation.
- **Mean and amplitude.** `T_w` starts as the air temperature at the water surface. The seasonal amplitude can be damped with mean lake depth, `A_w = A × max(1 − D / D_ref, a_min)`, representing heat storage, in the spirit of Air2Water's exchange with deeper water [[Piccolroaz et al., 2013]](#ref-piccolroaz2013).
- **Gate.** G3 fits this before it's used; if the fit doesn't beat `A_w = A`, damping is dropped.
- **Climate.** This is the water temperature #358 item 10 asks for. Climate's own lake moisture source is unchanged here; one shared evaporation flux is #358 item 11.

### D4. The depth caps

After D1–D3, measure how often `LAKE_MAX_DEPTH_M` and `SEA_MAX_DEPTH_M` bind, on which lakes, and with what supply surplus (G4). Options:
- keep both caps;
- keep them only as numerical guards with a warning, like `GLACIER_MAX_DEPTH_M`;
- remove them, leaving lakes bounded by their rim and their balance.

The choice is the owner's after G4. Water a cap holds back is booked as overflow at the cap, never dropped.

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
- **G4: depth caps.** After D1–D3 on the diagnosis set, report how often each cap binds, on which lakes and with what surplus. *Outcome:* the owner chooses among the D4 options.
- **G5: ice-cap budget.** On the diagnosis set over 10 Myr, run D2(c) with its starting knobs and report:
  - the distribution of cap thickness;
  - the share of cells frozen to the bed;
  - cap volume over time, which should approach a plateau;
  - the global sublimation flux.

  Compare cap thicknesses with published perennial and seasonal lake-ice thicknesses, sourced and verified in G5. *Outcome:* adopt the knobs, recalibrate them on the diagnosis set, or revise D2(c) before implementation.

  **Result (issue #541; report:** [`analysis/issue541/README.md`](../analysis/issue541/README.md)**):** The three #538 diagnosis saves reach a volume plateau by 0.1 Myr under fixed climate and lake depth. At baseline, 5.09% and 18.45% of wet area freeze to the bed on seeds 910211954 and 317739385, exceeding the 1% falsification threshold; seed 635056485 is at 0.54%. Across the requested sublimation (×0.25/×1/×4) and conduction (×0.5/×1/×2) sensitivities, seed 317739385 remains above threshold (minimum 1.63%). The baseline physical-ice thickness p50s are 4.50, 9.50 and 1.71 m respectively; thickness and bed-freeze shares vary widely across settings. **G5 fails for the starting D2(c) law. Revise its bed-freeze and loss treatment before implementation.** The baseline thicknesses are compared with seasonal lake ice (Hirose et al., 2008 report a 1.25 m regional mean seasonal maximum [[Hirose et al., 2008]](#ref-hirose2008)) and perennial lake ice (Lake Untersee's measured range is 2.2–3.8 m, summer mean 3.0 m [[Steel et al., 2015]](#ref-steel2015)). The 0.6 m/yr Untersee reference is summer ablation compensated by winter refreezing, not a measured annual sublimation flux [[Weisleitner et al., 2019]](#ref-weisleitner2019); this first-order prototype omits that basal refreezing and holds climate and depth fixed.

  **Revision after G5:** the diagnosis did not falsify conductive growth itself; it falsified the tested combination of bare-ice growth, fixed positive-degree-day melt, a uniform `s_ice` borrowed from summer ablation, and zero basal heat exchange. The next G5 candidate uses a snow-resisted, two-sided seasonal heat budget; records surface melt, basal melt/refreezing and sublimation as separate transfers; and reports bedfast area as a state. The starting `s_ice = 0.6 m/yr` and the no-snow Stefan coefficient are retired as adopted defaults. The new candidate must first show that its reduced forcing terms can be constrained against independent seasonal ice observations; otherwise D2(c) remains deferred. The existing diagnosis set may screen sensitivities, but the three held-out saves remain untouched.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `T_ice` (water freezing threshold) | 0 °C | Physical | None |
| `D_ref`, `a_min` (amplitude damping with depth) | Off (`A_w = A`) | Design choice; Air2Water's concept [[Piccolroaz et al., 2013]](#ref-piccolroaz2013) | G3 on observed lakes, never on world seeds |
| Snow thermal resistance `h_s/k_s` and surface temperature `T_s` | Not set; requires a separate lake-snow reservoir coordinated with RFC 0015 D4, or a reduced proxy validated against observed seasonal data | Snowless Stefan growth is an upper bound; current model has no snow-depth or surface-energy forcing | New G5 design screen against independent observed lakes before tuning on diagnosis saves |
| Basal lake-to-ice heat flux `q_w` | Not set; test a literature-bounded seasonal range rather than a global constant | Observed heat flux varies strongly by season, snow, insolation and lake stratification [[Huang et al., 2022]](#ref-huang2022) | New G5 design screen; select only a range supported by available model inputs and observations |
| `s_ice` (true sublimation from ice cover) | No adopted value; `0.6 m/yr` is removed as a default | Untersee's ~0.6 m summer ablation is compensated by winter basal refreezing, not a measured annual sublimation rate [[Weisleitner et al., 2019]](#ref-weisleitner2019) | Independently bound from observations; then sensitivity-test on the diagnosis set |
| Ice water-equivalent ratio | 1000/917 | Densities of water and ice | None |
| `SHALLOW_EVAPORATION_BOOST` | 0.3 (unchanged) | Existing (#297) | Out of scope |
| `LAKE_MAX_DEPTH_M`, `SEA_MAX_DEPTH_M` | Unchanged | Existing | G4, owner decision |
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
- **The ice cap's reduced heat budget remains underdetermined.** The no-snow Stefan law overestimates winter growth, while basal heat, snow insulation, solar energy, wind and humidity vary among lakes and are not all represented in the current model. The revised design separates those fluxes, but D2(c) cannot proceed until independent observations constrain a reduced parameterization; otherwise it needs a scoped lake thermal state or deferral.
- **Liquid water under ice** is assumed to keep its storage and spill; a frozen-to-bed cell is the only exception.
- **Evaporation is still the temperature-only rate.** It runs 15–35% low for continental deserts and ~70% low at high altitude (#519). Lake sizes will move again when #358 and #118 give a Penman rate; that interface is preserved.
- **Single-save evidence.** All mechanism sizes come from one 121.5 Myr world. G0 is required before implementation.

## Acceptance criteria

- [ ] **G0–G5 run and recorded** in this RFC, with the owner's D4 decision. G5's first candidate failed; the revised candidate must pass the independent forcing/observations screen before its diagnosis-set sensitivity rerun.
- [ ] **Unit tests:**
  - [ ] the endpoint solve matches alternative C at steady state (≤ 0.1% of volume) on random hierarchies;
  - [ ] nested and asymmetric supply cases hold, spill, and evaporate exactly their total source supply;
  - [ ] growing and shrinking single basins close their budget with no plug term;
  - [ ] a perched lake above a saddle doesn't merge its parent;
  - [ ] a merged lake whose balance falls below its saddle splits with each child at its own balance;
  - [ ] a spill crossing catchments via a channel is delivered to the hierarchy's receiving leaf;
  - [ ] the ice-free fraction is 0 below and 1 above the threshold range and continuous between;
  - [ ] a lake with `f_free = 0` keeps liquid water and balances its supply by sublimation and overflow;
  - [ ] ice-cap thickness stays within [0, water depth], closes the liquid/ice budget, and approaches a plateau under constant forcing;
  - [ ] basal melt and surface melt return water to liquid; sublimation is the only D2(c) atmospheric mass loss;
  - [ ] water at a cold sink is never deleted;
  - [ ] an enclosed lake cell with its bed below sea level uses the air temperature at its water surface;
  - [ ] silt equals coefficient × supply × years.
- [ ] **Ledger closes:** |residual| ≤ 0.1% of lake supply on every step of a 10-step replay on each held-out save.
- [ ] **No water deleted:** liquid water lost outside the ledger is 0 on held-out saves.
- [ ] **No pinned lakes:** no non-spilling, uncapped body (not flagged to C) whose endpoint evaporation differs from its supply by more than 5% (proposed, owner to confirm).
- [ ] **Flicker:** per step, the fraction of bodies above 10⁴ km² whose volume change reverses sign on consecutive steps with |ΔV| > 50% of V is ≤ 5% (proposed; ice-age control "never"). Global volume max/min is reported but isn't pass/fail.
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
- **The revised ice cap doesn't plateau or still freezes >1% of wet area to the bed:** *Revise:* D2(c)'s reduced heat budget again, or defer the cap if the missing thermal forcing cannot be constrained. The bedfast share must also be compared with bottomfast-lake observations before treating it as a physical failure; the >1% threshold remains the gate until the owner decides otherwise.
- **G3 gives no skill:** `1 − f_free` misses observed ice duration by more than 60 days median, even with observed inputs and damping. *Revise:* the sinusoid-and-air-temperature proxy is too crude; consider an ice-duration regression instead.

## References

- <a id="ref-barnes2020"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2020). Computing water flow through complex landscapes – Part 2: Finding hierarchies in depressions and morphological segmentations. *Earth Surface Dynamics*, 8, 431–445. [https://doi.org/10.5194/esurf-8-431-2020](https://doi.org/10.5194/esurf-8-431-2020)
- <a id="ref-barnes2021"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2021). Computing water flow through complex landscapes – Part 3: Fill–Spill–Merge: flow routing in depression hierarchies. *Earth Surface Dynamics*, 9, 105–121. [https://doi.org/10.5194/esurf-9-105-2021](https://doi.org/10.5194/esurf-9-105-2021)
- <a id="ref-benson2000"></a>Benson, B., Magnuson, J., & Sharma, S. (2000). *Global Lake and River Ice Phenology Database, Version 1*. NSIDC (G01377). [https://nsidc.org/data/G01377/versions/1](https://nsidc.org/data/G01377/versions/1)
- <a id="ref-leppaesranta1993"></a>Leppäranta, M. (1993). A review of analytical models of sea-ice growth. *Atmosphere-Ocean*, 31(1), 123–138. [https://doi.org/10.1080/07055900.1993.9649465](https://doi.org/10.1080/07055900.1993.9649465)
- <a id="ref-hirose2008"></a>Hirose, T., Kapfer, M., Bennett, J., Cott, P., Manson, G., & Solomon, S. (2008). Bottomfast Ice Mapping and the Measurement of Ice Thickness on Tundra Lakes Using C-Band Synthetic Aperture Radar Remote Sensing. *Journal of the American Water Resources Association*, 44(2), 285–292. [https://doi.org/10.1111/j.1752-1688.2007.00161.x](https://doi.org/10.1111/j.1752-1688.2007.00161.x)
- <a id="ref-huang2022"></a>Huang, W., Zhao, W., Zhang, C., Leppäranta, M., Li, Z., Li, R., & Lin, Z. (2022). Sunlight penetration dominates the thermal regime and energetics of a shallow ice-covered lake in arid climate. *The Cryosphere*, 16, 1793–1806. [https://doi.org/10.5194/tc-16-1793-2022](https://doi.org/10.5194/tc-16-1793-2022)
- <a id="ref-mason1994"></a>Mason, I. M., Guzkowska, M. A. J., Rapley, C. G., & Street-Perrott, F. A. (1994). The response of lake levels and areas to climatic change. *Climatic Change*, 27, 161–197. [https://doi.org/10.1007/BF01093590](https://doi.org/10.1007/BF01093590)
- <a id="ref-piccolroaz2013"></a>Piccolroaz, S., Toffolon, M., & Majone, B. (2013). A simple lumped model to convert air temperature into surface water temperature in lakes. *Hydrology and Earth System Sciences*, 17, 3323–3338. [https://doi.org/10.5194/hess-17-3323-2013](https://doi.org/10.5194/hess-17-3323-2013)
- <a id="ref-sharma2019"></a>Sharma, S., Blagrave, K., Magnuson, J. J., et al. (2019). Widespread loss of lake ice around the Northern Hemisphere in a warming world. *Nature Climate Change*, 9, 227–231. [https://doi.org/10.1038/s41558-018-0393-5](https://doi.org/10.1038/s41558-018-0393-5)
- <a id="ref-steel2015"></a>Steel, M., McKay, C. P., & Andersen, D. T. (2015). Modeling circulation and seasonal fluctuations in perennially ice-covered and ice-walled Lake Untersee, Antarctica. *Limnology and Oceanography*, 60(4), 1231–1245. [https://doi.org/10.1002/lno.10086](https://doi.org/10.1002/lno.10086)
- <a id="ref-weisleitner2019"></a>Weisleitner, K., Perras, A., Moissl-Eichinger, C., Andersen, D. T., & Sattler, B. (2019). Source environments of the microbiome in perennially ice-covered Lake Untersee, Antarctica. *Frontiers in Microbiology*, 10, 1019. [https://doi.org/10.3389/fmicb.2019.01019](https://doi.org/10.3389/fmicb.2019.01019)
- <a id="ref-zhao2022"></a>Zhao, G., Li, Y., Zhou, L., & Gao, H. (2022). Evaporative water loss of 1.42 million global lakes. *Nature Communications*, 13, 3686. [https://doi.org/10.1038/s41467-022-31125-6](https://doi.org/10.1038/s41467-022-31125-6)
