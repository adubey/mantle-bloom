# RFC 0014: Lakes conserve their water: in-step fill-spill-merge, seasonal freezing and water-surface temperature

- **Status:** Draft
- **Date:** 2026-10-22
- **Related:** #519 (exploration and evidence), #520 (perched-lake false merges, river lines across lakes, floor-drain flips), #521 (split reset), #505 (lake merges lag routing), #503 and #297 (breaching, balanced lakes), #355 (water units), #353 and #368 (glacier ice and its units), #358 items 10–11 (water-surface temperature, one evaporation flux), #118 (wind speed)

## Problem

Lakes in Mantle Bloom are not conserved bodies of water. Each step `lakes.resolve_lakes` places every lake at its *balanced level*, where evaporation over its flooded area equals its inflow, capped at its rim. That premise is sound: real closed lakes approach balance quickly compared with a 10 kyr–10 Myr step. But the machinery around it creates and deletes water, and it decides merging from last step's water rather than this step's. On a mature world this produces seas that appear from nowhere, vanish the next step, and hold far more water than their budget supports.

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
4. **Rivers freeze on the annual mean and their water is deleted.** A river reaching a land sink whose annual mean is below 0 °C becomes river ice at `hydrology.WATER_TO_ICE`, which is 1/50,000 of its volume in the glacier model's units. On the evidence save that was 2.5–48% of land-to-ocean discharge per step, mostly one river of ~150,000 m³/s at −0.9 °C. Most of that river was itself fed by item 1's false flood spilling north; under K2b the loss fell ~75%. The rule also contradicts the model's own principle #3, which uses exactly this case (a sub-zero annual mean doesn't mean no summer melt) as its example.
5. **Lakes read the wrong temperature.** Lake cells use the node's air temperature, and the climate reads it at the **bed** elevation:
   - a cell whose bed is below sea level reads the *ocean* temperature (`erosion.py`: `temperature = where(elevation <= sea_level, ocean_temperature_c, air_temperature_c)`), even inside an enclosed lake;
   - the lapse rate cools by the bed elevation clipped at 0, not by the water surface, so lake 13's shallower cells read about 3 °C too warm.
6. **A depth cap stands in for the balance.** Once item 2 is fixed, lake 13 rests at the sea-tier depth cap (floor + `SEA_MAX_DEPTH_M`, 3,500 m) on its largest steps, with inflow (3,100–10,300 km³/yr) above its evaporation at the cap (2,200–2,800 km³/yr). The cap was added to stop a runaway basin when nothing else bounded it (`lakes.py`, "Lake-vs-sea tiers"). Whether it is still needed once the balance is right is open.

### What the evidence does not show

- Whether these mechanisms are as large on other seeds, ages and densities (G0).
- How far #520 part 1 (connected-water merge test) and #521 (no split reset), both standalone fixes, already remove items 1 and 3 by themselves.
- Whether evaporation's temperature-only rate matters to lake size. It runs 15–35% low for continental deserts per #519's literature check, but that belongs to #358 and #118, not here.

**Principles violated:** #6 (account for material through every process: water is created by resets and deleted at frozen sinks and frozen lakes); #3 (model fast processes in aggregate over long steps: annual-mean freezing and last-step merging); #5 (consistent units: water becomes ice at 1/50,000).

## Proposed solutions

### What happens in nature

- **Closed lakes balance their budget.** A closed lake's area adjusts until evaporation from its surface matches inflow, and its level responds to a change in climate as a low-pass filter. The response time runs from years to a few centuries for closed lakes, and far shorter for lakes with outlets [[Mason et al., 1994]](#ref-mason1994) (*unverified: the abstract isn't publicly accessible; G1 measures the model's own response times instead*). On the evidence save, residence times (volume / inflow) were 0.1–140 years for sea 26's basins and about 350–430 years for lake 13 at its inflated, pinned size; that's at most 4% of the shortest 10 kyr step the UI offers.
- **Water fills, spills and merges.** Runoff fills depressions, which overflow into their neighbours; when both a depression and its neighbour fill, they merge [[Barnes et al., 2021]](#ref-barnes2021). The *depression hierarchy* records how sub-depressions merge into meta-depressions at their spill points [[Barnes et al., 2020]](#ref-barnes2020). Mantle Bloom's `lakes.build_lake_hierarchy` is this structure on a k-NN graph. Fill–Spill–Merge is mass-conserving but has no evaporation; adding it is this RFC's own extension.
- **Lake ice is seasonal and set mostly by air temperature.** Observations at 513 Northern Hemisphere lakes identify air temperature, lake depth, elevation and shoreline complexity as key controls on ice cover, and an estimated 14,800 lakes have intermittent winter ice today [[Sharma et al., 2019]](#ref-sharma2019). The Global Lake and River Ice Phenology Database holds freeze and breakup dates for 865 lakes and rivers [[Benson et al., 2000]](#ref-benson2000). Ice cover is a fraction of the year, not an annual-mean switch.
- **Lake surface temperature can be modelled from air temperature.** Air2Water relates the surface-layer temperature of lakes to air temperature alone, with a lumped heat exchange with the atmosphere and the deeper water [[Piccolroaz et al., 2013]](#ref-piccolroaz2013).
- **Evaporation depends on energy, the vapour-pressure deficit and heat storage.** A Penman-based global estimate, with a lake heat-storage term, puts long-term lake evaporation at 1500 ± 150 km³/yr, excluding the Caspian Sea [[Zhao et al., 2022]](#ref-zhao2022). The physical rate is out of scope here (#358, #118), but the lake solve should take evaporation from a function that a Penman rate can later replace.

### Alternatives

**A. Patch the current resolver** (#520 part 1, #521 and a wet-cells-only frozen rule).
- *For:* small, and the knockouts show it removes the largest artifacts.
- *Against:* merges still follow last step's water (#505), so routing and lakes still disagree for a step. Freezing stays all-or-nothing on the annual mean, frozen water stays unbooked, and nothing audits conservation.
- These patches are worth landing now; #520 and #521 do. They are not the end state.

**B. In-step fill-spill-merge with evaporation, seasonal per-cell freezing and a water-surface temperature (recommended).**
- Resolve all lakes each step from *this* step's routed water in one conservative pass over the hierarchy.
- No continuity pins, no "already merged" test.
- Freezing per cell from an ice-free season.
- Frozen water booked by water-equivalent volume.
- Temperature at the water surface.

**C. Prognostic lake storage.** Carry each lake's volume and integrate `dV/dt = inflow − evaporation − overflow` through the step with substeps.
- *For:* handles lakes whose response time approaches the step.
- *Against:* more state and cost, and unnecessary if response times are far shorter than steps.
- Kept as the fallback if G1 fails.

**D. A full lake thermal model** (FLake-style mixed layer and ice).
- Too data-hungry for this model's inputs (no physical humidity, wind or radiation yet; #358, #118) and too costly at 10⁴–10⁵ lakes per world.

**Recommendation: B, with C held as the fallback that G1 decides.** B keeps the equilibrium premise, which the residence times support, and removes every mechanism in the Problem by construction rather than patch by patch. It makes routing and lakes one solve, and gives the water budget an auditable ledger. #520 and #521 are interim fixes that B supersedes.

## Design

### D1. In-step balanced fill-spill-merge

Each step, after runoff and meltwater are routed to sinks (`route_downstream` with the existing transit losses), lakes are solved once over the hierarchy from this step's water:

- **Inputs, all in m³/yr:**
  - `settled[s]`, the water arriving at each leaf's sink;
  - each cell's evaporation rate `e[i]` (m/yr; D2 and D3);
  - each cell's hypsometry (`CellHypsometry`: area and sub-cell relief, as today).
- **The rule for a depression d with inflow `I_d`:**
  - If the evaporation of d filled to its rim, `E_d(rim)`, is at least `I_d`, d holds water at the balanced level `h_d` where `E_d(h_d) = I_d`, found by bisection (`lakes.balanced_level`, unchanged), and spills nothing.
  - Otherwise d is full: it evaporates `E_d(rim)` and spills `O_d = I_d − E_d(rim)` through its outlet, `Lake.outlet_target_idx`.
- **Where the spill goes:**
  - into the sibling's subtree: injected at the outlet target and routed by steepest descent to the receiving leaf's sink, as `_spill_within_step` does today;
  - from a root, to the ocean.
- **Merging:** a parent floods as one body only when *both* children are full in this step's solve. Its inflow is then the sum of its leaves' inflows (including spills received), and it is balanced over all its members, capped at its own rim.
- **No state from last step decides any of this.** `_prev_level`, the `already_merged` test, the merge pin, the split reset and its tier-capped variant are removed. This also removes #505's lag and the perched-lake false merge, and makes the solve the single source of the spill decisions routing uses (`should_spill`, `filled_elevation`, `spill_target`). Those then come from this step's solve, replacing `governing_bodies`' dependence on last step's water.
- **Processing order:** leaves before parents. A spill into a sibling subtree raises that subtree's inflow, so the sibling is re-solved. Since `E_d(h)` is non-decreasing in `h` (`balanced_level` relies on this; `SHALLOW_EVAPORATION_BOOST < 1` keeps it so), adding water never lowers a level and the final state should be unique and independent of order. That claim is a **unit-test obligation**: compare against a brute-force fixed-point iteration on random small hierarchies.

**Storage change.** The solve gives this step's balanced volume `V_new` for each lake. The difference from last step's volume `V_old = Σ lake_depth × area` is real water: a growing lake takes `V_new − V_old` from the step's inflow, and a shrinking one returns it to evaporation and overflow.
- In the ledger it is spread over the step, as `ΔV / years` in m³/yr.
- This is the only use of last step's state. It keeps the budget closed without changing where water ends up.
- Equilibrium is valid while `V_new − V_old` is small next to `inflow × years`; G1 checks that. Where it isn't, alternative C applies to that lake.

**Ledger** (new, per step, m³/yr):
- **sources:** settled runoff and meltwater into basins;
- **sinks and reservoirs:**
  - evaporation;
  - overflow to the ocean;
  - storage change;
  - perennial ice (D2), in water equivalent;
- **internal transfers**, outside the balance: spills between basins.

The residual is a test target (acceptance). It replaces today's absence of any lake-water audit.

**Bounds:** levels stay within [dry level, rim]; spills are ≥ 0; evaporation is ≥ 0 and never exceeds the inflow plus the water released by storage change.

**Timestep:** no step-count dependence. Flows are annual rates and storage change is booked per year. G1 checks that response times stay below 10% of the step at the shortest step, 10 kyr.

### D2. Per-cell seasonal freezing, with frozen water booked

- **Ice-free fraction.** Each lake or river cell's fraction of the year free of ice is `f_free = φ / π`, with `φ = arccos(−(T_w − T_ice) / A)` from the same sinusoidal annual cycle the glacier model already uses (`hydrology.seasonal_positive_degrees`; amplitude `A` from `biomes.seasonal_temp_amplitude`).
  - `T_w` is the water-surface temperature (D3) and `T_ice` the freezing threshold (0 °C).
  - `f_free = 1` when the whole year is above the threshold and 0 when the whole year is below.
- **Evaporation** from a wet cell is `e = e_open(T_open) × f_free`, where `T_open` is the mean water temperature over the ice-free part of the year, in closed form from the same sinusoid. Today `e_open` is `breaching.open_water_evaporation_mm`; a Penman rate can replace it later (#358).
- **No whole-lake freezing rule.** A lake whose cells all have `f_free = 0` evaporates nothing; partly ice-covered lakes simply evaporate less.
- **Perennially frozen water.** Inflow that a fully frozen lake can't hold below its rim becomes **perennial lake ice**. So does water reaching a land sink whose `f_free = 0`. It is booked in a declared ice account in **water-equivalent m³** (×1000/917 for ice volume), not converted at `WATER_TO_ICE`.
  - Until glacier mass balance is physical (#368), this account is kept apart from `glacier_depth`, so it can't inflate glaciers by 50,000×.
  - After #368, it can feed the glacier model by mass.
- **Rivers.** A river cell with `f_free > 0` carries its annual flow; seasonal freezing delays water but doesn't stop it over a step. Only `f_free = 0` stops flow.
- **This replaces:**
  - the annual-mean `water_frozen` / `lake_frozen` tests in `hydrology.compute_hydrology`;
  - `frozen_from_lake`;
  - the any-member `is_frozen[node.members].any()` in `lakes._resolve`;
  - the frozen-river-ice conversion for liquid water.
- **Unchanged:** precipitation still falls as snow on frozen or ice-capped nodes (a separate rule; principle #8).

### D3. Water-surface temperature

- **Elevation.** Lake cells take their temperature at the **water-surface elevation**: the lake level from the D1 solve of the previous step (a one-step lag that climate already tolerates), not the bed.
  - Enclosed lake cells use the land and air temperature path even where the bed is below sea level. Only ocean cells (`is_ocean`, connected to the world ocean) use the ocean temperature.
  - The lapse rate cools by the water-surface elevation.
- **Mean and amplitude.** `T_w` starts as the air temperature at the water surface. The seasonal amplitude is damped with mean lake depth, `A_w = A × max(1 − D / D_ref, a_min)`, representing heat storage, in the spirit of Air2Water's exchange with deeper water [[Piccolroaz et al., 2013]](#ref-piccolroaz2013). `D_ref` and `a_min` are calibration knobs (below).
- **Gate.** G3 fits these to observed ice duration before they are used; if the fit doesn't beat `A_w = A`, damping is dropped.
- **Climate.** This is the water temperature #358 item 10 asks for. Climate's own lake moisture source is unchanged here; one shared evaporation flux is #358 item 11.

### D4. The sea-tier depth cap

After D1–D3, measure how often `LAKE_MAX_DEPTH_M` and `SEA_MAX_DEPTH_M` bind, on which lakes, and with what inflow surplus (G4). Options:
- keep both caps;
- keep them only as numerical guards with a warning, like `GLACIER_MAX_DEPTH_M`;
- remove them, leaving lakes bounded by their rim and their balance.

The choice is the owner's after G4. Any water a cap holds back must be booked (overflow at the cap, or a declared sink), not dropped.

### Diagnosis gates

- **G0: baseline after #520 and #521.** Once those land, measure items 1–4 of the Problem on **three or more new mature saves** (other seeds; not 308161711, which is the development save) at 100 kyr steps:
  - reset-created water;
  - false merges;
  - whole-lake freezing;
  - frozen-sink water loss;
  - flicker (global lake volume max/min over 10 steps);
  - balance residual of non-spilling lakes.

  *Outcomes:* if the interim fixes already meet D1's acceptance thresholds on every save, D1 shrinks to the ledger plus the #505 unification, and the RFC is revised. Otherwise D1 proceeds as designed.
- **G1: is equilibrium valid?** On the same saves at 10 kyr and 100 kyr steps, the distribution of lake response times `τ = V / inflow` and of `|ΔV| / (inflow × years)`.

  *Outcome:* if ≥ 99% of lake volume has `τ ≤ 0.1 ×` step at 10 kyr, D1 stands. Otherwise use alternative C for lakes over the threshold.
- **G2: is the fill-spill-merge solve unique?** On random hierarchies (unit-test scale) and on the G0 saves' real hierarchies, compare the D1 solve with brute-force fixed-point iteration of the same rules.

  *Outcome:* if they agree to ≤ 0.1% of volume, D1 stands. Otherwise, find the order dependence before implementation.
- **G3: ice-free season and water temperature.** Using the Benson et al. database, take ≥ 50 lakes with long records and an independent mean-air-temperature estimate. Compare observed ice duration with `1 − f_free` from D2 and D3, with and without depth damping.

  *Outcome:* adopt damping only if it lowers the median absolute error. Report that error, with a target of ≤ 30 days.
- **G4: depth caps.** After D1–D3 on the G0 saves, report how often each cap binds, on which lakes and with what surplus. *Outcome:* the owner chooses among the D4 options.

### Knobs

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `T_ice` (water freezing threshold) | 0 °C | Physical | None |
| `D_ref`, `a_min` (amplitude damping with depth) | Off (`A_w = A`) | Design choice; Air2Water's concept [[Piccolroaz et al., 2013]](#ref-piccolroaz2013) | G3 on the ice-phenology lakes, never on held-out world seeds |
| Ice water-equivalent ratio | 1000/917 | Densities of water and ice | None |
| `SHALLOW_EVAPORATION_BOOST` | 0.3 (unchanged) | Existing (#297) | Out of scope |
| `LAKE_MAX_DEPTH_M`, `SEA_MAX_DEPTH_M` | Unchanged | Existing | G4, owner decision |
| G1 response-time threshold | 0.1 × step | Design choice | Owner may revise |

**Seeds.** Development: 308161711 (#519). Tuning: none on world seeds; G3 uses observed lakes. Held-out acceptance: the G0 saves, three or more seeds chosen before G0 runs and recorded in the tracking issue.

### Cost and saves

- **Runtime.** The solve visits each hierarchy node a small number of times, with one bisection per wet body: the same order as today's `_resolve`. It replaces `_resolve` and the fixed-point rounds of `_spill_within_step`. The ice-free fraction is closed-form per node. Target: ≤ +10% median step time at density 4 (acceptance).
- **Save compatibility.**
  - `lake_depth` is already persisted and gives `V_old`.
  - The perennial-ice account is new. If it is per node, it's a new `SURFACE_FIELDS` entry, default 0, an extensive remap, with old saves loading at 0 (`docs/save-compatibility.md`, `unit_tests/test_persistence.py`). If it is aggregate, it lives on the ice ledger and needs nothing.
  - The lag in water-surface elevation reads the persisted `lake_depth`.

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/lakes.py`:
  - `_resolve`, `_water_balance`, `_prev_level`, `governing_bodies`, `compute_spill_routing` and `resolve_lakes` are replaced or rewritten around a new in-step solve (new, e.g. `solve_lakes`);
  - `balanced_level`, `lake_evaporation_m3_per_yr` (gains `f_free`), `CellHypsometry`, `build_lake_hierarchy`, `_mark_floor_drains` and `_classify_tier` are kept.
- `backend/app/hydroclimate/hydrology.py`:
  - `compute_hydrology`: the order of routing, solve and spill routing; `water_frozen`, `lake_frozen`, `freezing_lake`, `frozen_from_lake` and the frozen-river-ice block;
  - `_spill_within_step` and `SpillBasins`: folded into the solve;
  - `_compute_flow_direction`: freeze override from `f_free`;
  - `seasonal_positive_degrees`: reused, plus a new sibling for the ice-free fraction and open-season mean;
  - `WATER_TO_ICE` is no longer used for liquid water.
- `backend/app/hydroclimate/erosion.py`: the `temperature` assembly (~L1926; enclosed lake cells and water-surface lapse) and where silt reads the solve's inflow.
- `backend/app/hydroclimate/ice_sheet.py`: the `IceLedger` perennial-lake-ice account (new), and `frozen_river_ice_m3` semantics.
- `backend/app/hydroclimate/breaching.py`: `open_water_evaporation_mm` is unchanged and called with `T_open`. Check `_overflows`' demand against the solve's evaporation so the breach gate and lakes agree.
- `backend/app/main.py`: `_lake_basin_summary` and `_smallest_lake_containing` read `is_spilling` and levels from the new solve.
- **Tests:**
  - `backend/unit_tests/test_lakes.py`: merge, split, tier and frozen tests rewritten for the solve; new uniqueness (G2), conservation and perched-lake tests;
  - `test_hydrology.py`: seasonal freezing, frozen sinks;
  - `test_persistence.py`, if a per-node ice field is added;
  - `backend/stress_tests/` for the multi-step replays.
- **Docs:** `docs/simulation-model/hydrology.md` (balanced lakes, lake tree, tiers, frozen lakes, rivers), `docs/simulation-model/climate.md` (lake temperature), `docs/save-compatibility.md` if needed.
- **Diagnostics (new):** a `bin/debug/` replay script measuring G0, G1 and G4 metrics and the lake ledger, built from #519's probes (in the shared handoff directory), with a usage docstring.

## Assumptions and uncertainty

- **Equilibrium within a step.** Supported by residence times on one save (0.1–430 years). Mason et al.'s range is unverified here. G1 tests it at 10 kyr steps on new saves.
- **A unique fill-spill-merge state** with evaporation added is argued from `E(h)` being non-decreasing, not proven. G2 tests it.
- **A sinusoidal annual cycle** with the existing synthesized amplitude stands in for real seasonal temperature. Air temperature stands in for water temperature, apart from optional depth damping. Sharma et al. name depth, elevation and shoreline complexity as further controls; only depth is represented. G3 measures the error.
- **The water-equivalent ice account** is a holding reservoir until glacier units are physical (#368). Its growth on held-out seeds is a falsification check (below).
- **Evaporation is still the temperature-only rate.** It runs 15–35% low for continental deserts and ~70% low at high altitude (#519). Lake sizes will move again when #358 and #118 give a Penman rate; that interface is preserved.
- **Single-save evidence.** All mechanism sizes come from one 121.5 Myr world. G0 is required before implementation.

## Acceptance criteria

- [ ] **G0–G4 run and recorded** in this RFC, with the owner's D4 decision.
- [ ] **Unit tests:**
  - [ ] the D1 solve equals brute-force fixed-point iteration on random hierarchies (≤ 0.1% of volume);
  - [ ] a perched lake above a saddle doesn't merge its parent;
  - [ ] a merged lake whose balance falls below its saddle splits with each child at its own balance (no water created);
  - [ ] the ice-free fraction is 0 below and 1 above the annual threshold range and continuous between;
  - [ ] a partly frozen lake keeps water and evaporates less;
  - [ ] water at an `f_free = 0` sink is booked as perennial ice, not lost;
  - [ ] an enclosed lake cell with its bed below sea level uses the air path at its water surface.
- [ ] **Lake water ledger closes:** |residual| ≤ 0.1% of lake inflow on every step of a 10-step replay on each held-out save.
- [ ] **No water deleted:** on held-out saves, water reaching frozen sinks or frozen lakes that isn't booked to the ice account is 0.
- [ ] **No pinned lakes:** no non-spilling, uncapped lake whose evaporation differs from its inflow by more than 5% (proposed threshold, owner to confirm).
- [ ] **Stability:** global visible lake volume, max/min over 10 × 100 kyr steps, ≤ 1.5 on each held-out save (proposed; the evidence save measured 2.4 before and 1.4 under K2b).
- [ ] **Step robustness:** global lake volume and area after 1 Myr agree within 15% between 100 × 10 kyr and 10 × 100 kyr steps on each held-out save.
- [ ] **Resolution:** global lake area agrees within 20% between node densities 4 and 6 on one held-out seed (proposed).
- [ ] **Runtime:** median step time at density 4 within +10% of `origin/main`.
- [ ] **Saves:** old saves load (`test_persistence.py`); any new field defaults safely.
- [ ] **Docs** updated as listed above.

## Falsification criteria

- **G1 fails:** more than 1% of lake volume has `τ > 0.1 ×` step at 10 kyr. The equilibrium premise is wrong at short steps. *Revise:* adopt alternative C for those lakes; don't reject the RFC.
- **G2 fails:** the solve depends on processing order beyond 0.1% of volume. A unique balanced state doesn't exist under these rules. *Revise:* fix the order (for example, elevation order with iteration) or move to alternative C.
- **G0 shows the interim fixes suffice:** #520 and #521 already meet the stability, no-pinned-lakes and no-water-deleted criteria on every held-out save. *Revise:* withdraw D1's solve, keep the ledger, D2 and D3.
- **Runaway perennial ice:** the account grows without bound over 10 Myr on a held-out save, or high-latitude closed basins accumulate ice faster than any precipitation could supply. *Revise:* D2's frozen-sink rule; the premise that only `f_free = 0` stops flow may be wrong at coarse resolution.
- **G3 gives no skill:** the D2–D3 ice-free fraction misses observed ice duration by more than 60 days median even after damping. *Revise:* the sinusoid-and-air-temperature proxy is too crude; consider an ice-duration regression instead.

## References

- <a id="ref-barnes2020"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2020). Computing water flow through complex landscapes – Part 2: Finding hierarchies in depressions and morphological segmentations. *Earth Surface Dynamics*, 8, 431–445. [https://doi.org/10.5194/esurf-8-431-2020](https://doi.org/10.5194/esurf-8-431-2020)
- <a id="ref-barnes2021"></a>Barnes, R., Callaghan, K. L., & Wickert, A. D. (2021). Computing water flow through complex landscapes – Part 3: Fill–Spill–Merge: flow routing in depression hierarchies. *Earth Surface Dynamics*, 9, 105–121. [https://doi.org/10.5194/esurf-9-105-2021](https://doi.org/10.5194/esurf-9-105-2021)
- <a id="ref-benson2000"></a>Benson, B., Magnuson, J., & Sharma, S. (2000). *Global Lake and River Ice Phenology Database, Version 1*. NSIDC (G01377). [https://nsidc.org/data/G01377/versions/1](https://nsidc.org/data/G01377/versions/1)
- <a id="ref-mason1994"></a>Mason, I. M., Guzkowska, M. A. J., Rapley, C. G., & Street-Perrott, F. A. (1994). The response of lake levels and areas to climatic change. *Climatic Change*, 27, 161–197. [https://doi.org/10.1007/BF01093590](https://doi.org/10.1007/BF01093590)
- <a id="ref-piccolroaz2013"></a>Piccolroaz, S., Toffolon, M., & Majone, B. (2013). A simple lumped model to convert air temperature into surface water temperature in lakes. *Hydrology and Earth System Sciences*, 17, 3323–3338. [https://doi.org/10.5194/hess-17-3323-2013](https://doi.org/10.5194/hess-17-3323-2013)
- <a id="ref-sharma2019"></a>Sharma, S., Blagrave, K., Magnuson, J. J., et al. (2019). Widespread loss of lake ice around the Northern Hemisphere in a warming world. *Nature Climate Change*, 9, 227–231. [https://doi.org/10.1038/s41558-018-0393-5](https://doi.org/10.1038/s41558-018-0393-5)
- <a id="ref-zhao2022"></a>Zhao, G., Li, Y., Zhou, L., & Gao, H. (2022). Evaporative water loss of 1.42 million global lakes. *Nature Communications*, 13, 3686. [https://doi.org/10.1038/s41467-022-31125-6](https://doi.org/10.1038/s41467-022-31125-6)
