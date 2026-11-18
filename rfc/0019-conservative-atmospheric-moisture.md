# RFC 0019: Conserve atmospheric moisture in diagnostic rain and transport

**Status:** Accepted (revised 2026-11-18 after three-round G3 failure and owner-approved exception; revised 2026-11-17 after owner-approved G3 regional checks and held-out result; revised 2026-11-16 after G3 tuning and owner scope clarification; revised 2026-11-15 after G3 tuning screens; revised 2026-11-14 after G0; owner confirmation 2026-11-14)
**Date:** 2026-11-03
**Revisions:** See [Revision note](#revision-note).
**Related:** [tracking #898](https://github.com/adubey/mantle-bloom/issues/898), [#771](https://github.com/adubey/mantle-bloom/issues/771), RFC 0015, RFC 0016; issues #522, #526, #530, #358

## Problem

The climate estimates humidity with two passes across its grid: one along longitude (east-west) and one along latitude (north-south). Over land, the model carries less moisture farther along the wind path. Windward slopes remove more, and the model tracks that slope-related amount separately. The ordinary loss as moisture crosses land has no recorded destination. After the passes, the model calculates rain from the remaining humidity at each cell and adds rain on slopes and where winds flow together. This later rain is not subtracted from moisture carried downwind; only the slope-related loss reduces the carried amount. The model also discards moisture above a limit, adds a small minimum amount over dry land, and resets ocean humidity based on water temperature without recording evaporation as added water. These steps add or remove water without accounting for it, so the model cannot balance water entering the atmosphere with rain, other losses, and changes in atmospheric storage.

This violates **design principle 6, account for material through every process**: the same atmospheric water can rain repeatedly without replenishment, and other water is silently removed or introduced at caps and clamps. It also strains principle 5, explicit and consistent units: near-surface specific humidity in g/kg is used as a surrogate for a moisture source and annual precipitation through a fitted conversion.

The downstream effect is measured. RFC 0016's [G1b result](../analysis/rfc0016/README.md#g1b-result) reports that coupling plant transpiration to humidity raised land rain by 38–51% and increased water exported to the ocean. The report attributes that multiplication to rain not being debited from humidity. Those are measurements on its pinned source and windows, not a forecast of a conserving model's result. G1b cannot evaluate the intended land-water coupling until this RFC's question is resolved.

The scope is the annual-mean diagnostic water cycle used by erosion, hydrology, biomes, and ice. The model should produce visually recognizable rain belts, wetter windward slopes and drier leeward rain shadows, dry continental interiors, and plausible precipitation on cold land where moisture can reach it. Approximate intensities and geographic positions are acceptable on synthetic worlds. The water budget must still close, and the delivered rain and snow must support plausible downstream behavior. Clouds, individual storms, seasonal variability, and vertical atmospheric structure are represented only through bulk effects.

## Proposed solutions

### A. Conservative two-dimensional column-water budget (recommended for diagnosis)

Represent moisture as the total water in an air column. On each climate update, calculate a steady balance: evaporation adds water, wind moves it between grid cells, and rain removes it. The same amount that leaves one cell enters the next, so transport cannot create or lose water overall. Repeat the calculation until each cell stops gaining or losing water. This is a numerical way to find a steady state, not another simulated day or year.

Use a positive, bounded precipitation sink linked to column water, initially `P = W / τp`, and test a saturation-triggered or orographic enhancement only as a sink that removes available water. Wind convergence is already in the finite-volume flux divergence; do not add it a second time as unbudgeted rain. A source tag (ocean, lakes/rivers, vegetation) is useful for diagnostics, but total water—not each tag's separate transport—is the closure requirement.

This option replaces the two blended one-dimensional sweeps with a low-order, conservative two-dimensional solve. It can represent meridional and zonal transport together, wrap longitude, close fluxes at the poles, and use both diagnostic and CFD winds. Its uncertainties are the effective column-mean transport speed, rainout law, and how annual-mean evaporation is represented.

### B. Source-fed Markov transport along local wind arrows

Represent the climate grid as a directed graph. Each cell's local wind arrow sets transition rates from that cell to its downwind neighbors; a calm cell has no horizontal transitions. Inject evaporation as water mass, and let the same rainout law as option A remove water from each cell. Solve the source-fed stationary balance for the unnormalized water mass in each cell. This is the stationary humidity field; the water removed at each cell per unit time is its rainfall. Normalize only when reporting source shares, not when calculating humidity or rain.

For cell mass `M_i = A_i W_i` in kg, source `S_i = A_i E_i` in kg/s, wind-driven transition rate `q_ij` from cell `i` to `j` in 1/s, and rainout rate `λ_i` in 1/s, solve `0 = S_i + Σ_j(q_ji M_j) - M_i(λ_i + Σ_j q_ij)`. Set `λ_i = 1/τp` initially. Derive the outgoing `q_ij` from the downwind component of the local wind arrow and cell geometry; outgoing rates are nonnegative and sum the transfers out of a cell. Rainfall is `R_i = λ_i M_i` in kg/s. This balance gives absolute water mass; normalizing `M` would discard the source scale and is only appropriate for reporting source shares.

Incoming water adds at converging cells, and outgoing water splits across downwind branches at diverging cells. Closed arrow loops are handled by solving the full stationary balance, not by stopping a sweep at an arbitrary point. With a positive rain sink, loop water can leave through rainfall. A closed class with no sink has no finite unique stationary solution and must be detected and reported as a failed prototype, not silently drained or discarded. Source tags can show where rain came from, but total water must close. This is an arrow-following alternative to option A's grid-based finite-volume solve; its runtime advantage is unproven and must be measured.

### C. Prognostic moisture inside the atmosphere CFD

Add column vapour, evaporation, condensation, and precipitation to `AtmosphereCFDState`, and advance them with atmospheric flow. This makes storage explicit but couples annual-mean surface hydrology to a solver that currently advances one simulated day per tectonic step. It also adds persistent state, save migration, CFL/runtime questions, and a substantial calibration problem. The retired CFD humidity/precipitation path produced near-zero, uncalibrated rainfall and starved erosion and hydrology; adding state alone would repeat that failure mode. Revisit only if the steady diagnostic approaches fail a named criterion and a longer-term atmosphere/hydrology coupling is designed.

### Recommendation

Keep option A as the next G3 prototype. G2 established that options A and B can both conserve water. The G3 screens show that the next design must improve regional water delivery as well as the cost of the solve. Use a coarse finite-volume budget with explicit mixing, then place the rain it supplies on the native climate grid. Use current diagnostic and CFD winds as separate inputs to the same budget. Option C remains a fallback if a steady diagnostic cannot produce useful climates at acceptable cost.

## Design

### Budget and transport

For every climate cell `i`, use area `A_i` in m² and column water `W_i` in kg/m². Let `E_i` and `P_i` be surface-to-air and air-to-surface fluxes in kg/m²/s. Let `F_ij` be signed vapour mass transport across a shared edge in kg/s. Over solve interval `Δt`, audit

`A_i(W_i,new - W_i,old) = Δt[A_i(E_i - P_i) - Σ_j F_ij]`.

The global residual is `Σ A_i W_i,new - Σ A_i W_i,old - Δt Σ A_i(E_i-P_i)`; internal edge fluxes cancel, leaving only documented open boundaries (none on the closed sphere). The steady solution must have `W_i,new = W_i,old` within solver tolerance and global evaporation equal to precipitation. Report both per-cell normalized residuals and the global residual. Do not count solver relaxation as physical storage change.

Use a conservative donor-cell flux on the spherical climate grid. Cell area is calculated from latitude-edge bands and longitude width, rather than treating every latitude row as equal area. Longitude wraps; polar faces have zero area; east/west fluxes use the same value with opposite signs in their two adjacent cells. Clip negative numerical values to zero only with their clipped mass added to the local precipitation/sink ledger; otherwise fail the solve. Preserve the climate-density meaning of `World.climate_density`: cell fluxes and totals are area-weighted and should converge with resolution.

Use `surface_wind` from either wind mode for the first prototype. It is in m/s and targets a 10 m wind; the CFD solver currently reaches only 0.42–0.66 of its target on the reported diagnosis. The moisture-bearing column wind is not observed by this model. Start with an effective transport multiplier of **1.0 (design choice)**. Retain 0.5, 1.0, and 2.0 as sensitivities and report diagnostic and CFD wind separately. The ERA5 vertically integrated vapour transport comparison informs the uncertainty in this proxy; any selected multiplier must earn its place by improving delivery on the tuning worlds. Label column water times near-surface wind as modeled transport.

### Coarse transport and fine rain placement

The initial G3 screen held evaporation and rainout fixed and tested
shared-face mixing with both wind modes at 2°.
Mean flow omits moisture carried by transient circulation
[[Seager et al., 2010]](#ref-seager2010); down-gradient exchange is a
simple test of redistribution, with no guarantee that it reproduces the
direction of that missing transport. Preserve water transferred across
coasts and between latitude bands in the ledger. Keep separate land and
ocean rain budgets wherever a coarse cell includes both.

Spend each coarse rain budget exactly once. Use bounded terrain exposure and
a modest rain-belt cue inspired by the current climate to distribute it over
the native cells. This can improve windward slopes and rain shadows without
introducing an extra local rain sink that intercepts moisture before it can
reach cold land. The broad wet belts, dry interiors and cold-land supply still
depend on the regional solution. Compare rain maps, biomes and those regional
metrics alongside the approved G3 tests. A transport-only multiplier may be
selected separately for each mode from the existing sensitivity range on tuning
worlds if a paired comparison supports it. Keep it separate from the actual
surface wind and evaporation calculation.

The [G3 implementation plan](../analysis/rfc0017/g3-coarse-transport-plan.md)
specifies the grid, fluxes, coastal bookkeeping, weights, checks and timing
protocol. The separate [wind follow-up](../analysis/rfc0017/wind-transport-handoff.md)
describes how to isolate and test possible changes to wind direction, drag and
the CFD surface conversion. Those changes are outside this RFC's transport
prototype.

The current tuning candidate uses a convergence-adjusted carrier with a
limited broad meridional contribution. A local saturation cap removes
excess vapor once. A bounded assignment of that condensate to nearby
cells spreads the wet belts and preserves ocean, land and freezing-land
rain totals. Terrain weighting inside the coarse budgets produced weak
rain shadows. Intercepting vapor as it crosses a ridge made the windward
cells wetter and the lee drier, while reducing cold inland supply on A CFD.
Trace the moisture reaching that region before selecting a ridge rule.

### Rain and the atmospheric column

Start with a first-order sink `P = W / τp`, with `τp = 9.1 days` as a **literature-constrained initial global scale**, from the ratio of 24.6 mm global mean precipitable water to 2.69 mm/day precipitation in Trenberth (1998). This is an Earth global-mean sanity check, not a per-cell law or a target precipitation climatology. Test a 5-day and 18-day sensitivity and compare latitude- and surface-class-dependent depletion times inferred from ERA5 before selecting any spatial variation. The sink must be nonnegative and cannot remove more water than is present in the finite-volume update. A cell's rain plus exported flux cannot exceed its opening column inventory and incoming flux plus evaporation over the solve interval.

The first recommendation has no hard humidity floor or vapor cap. Saturation is represented by routing excess column water to precipitation, with the local excess debited exactly once. If the atmospheric layer cannot support that treatment at the selected grid and timestep, propose a named temporary reservoir or source/sink and show its budget. Keep `humidity_g_kg` as the near-surface interface for existing readers by deriving it from `W` using an explicit effective mixing depth and air density; diagnose the conversion against observed precipitable water and ERA5 near-surface humidity. A maximum near-surface humidity remains a diagnostic bound, not a place to discard water.

Orographic lift may increase the precipitation sink on windward slopes only when it removes moisture from `W`; the corresponding rain is credited to `P`. Moisture-flux convergence and subsidence emerge from the divergence of `F`. Any extra convective efficiency or subsidence parameter must be bounded, have units, and repartition available water rather than create or delete it. Coherent climate noise may perturb transport or the rainout rate, but the field must remain nonnegative and its area-integrated source and sink effects must be included in the ledger.

### Physical surface sources and ownership

- **Ocean evaporation:** replace the humidity reset with a surface evaporation flux in kg/m²/s. Prototype a bulk aerodynamic flux from sea-surface saturation humidity, near-surface air humidity, surface wind speed, air density, and a transfer coefficient. This bulk form is used in ocean latent-heat-flux products [[EUMETSAT SAF CM, 2011]](#ref-eumetsat2011). Bound evaporation to available ocean water. Diagnose the coefficient from ERA5 and test sensitivity; the model's ocean temperature is a diagnostic and must be used at its actual cell value. The G3 tuning-only prototype using legacy near-surface humidity did not resolve spatial precipitation error, so it is not selected for implementation.
- **Lakes and rivers:** convert the current depth-based humidity proxies into explicit annualized mass fluxes. Freeze or dry water sources contribute zero. Use only liquid water the surface owns; lake/river evaporation already removed by `hydrology.py` must not be charged again to the atmosphere. The gate must inventory the current `RIVER_EVAPORATION_*` and `LAKE_EVAPORATION_*` terms before choosing one owner.
- **Vegetation:** express transpiration as a mass flux. RFC 0016's current 0.16 recycling fraction of prior precipitation may be a starting calibration hypothesis, but it is not a measured ET law. Convert the prior rain from mm/yr to kg/m²/s, apply the biome coefficient once, and ensure RFC 0016's land-water store or another explicit surface reservoir owns the debit. When its water source is not available, transpiration is zero rather than an unbooked input.
- **Background vapour:** remove `MIN_BACKGROUND_HUMIDITY_G_KG` as an independent source. If a trace precipitation floor is still required for ice accumulation, make it an explicit transfer from a named surface reservoir or retain no floor and measure the glaciation impact.

[[Oki & Kanae, 2006]](#ref-okikanae2006) report that the ocean supplies about 87% of global evaporation and receives about 78% of global precipitation. Use these as broad Earthlike scale checks alongside land/ocean rain share; do not force local cells or each synthetic world to those fractions. [[van der Ent et al., 2010]](#ref-vanderent2010) estimate 40% of terrestrial precipitation originates from land evaporation in their ERA-Interim accounting. This supports source-tagged recycling diagnostics, not an exact per-cell transpiration fraction.

### Scales, outputs, and compatibility

[[Trenberth, 1998]](#ref-trenberth1998)'s global-mean column water and precipitation give a depletion time near nine days. A deliberately rough transport check using 25 kg/m² of column water and 5 m/s gives about 125 kg/(m·s) of moisture flux. At 1,000 mm/yr of rain (about 3.2×10⁻⁵ kg/m²/s), a one-dimensional flux-to-sink ratio is about 4,000 km. This is only a consistency check: it assumes uniform flow, a single transport direction, and no replenishment, and is not a fitted inland decay distance. It shows why a 380 km half-life should not be retained as the rain law without a moisture budget.

Convert `P` to the existing annualized precipitation field in mm/yr at one named interface. Resample it to terrain nodes with area-aware weights when integrating totals; do not sum nearest-neighbor grid values as if their areas were equal. Preserve the existing `ClimateFields` and `climate_cache` readers where possible. The recommended solve adds no persistent `World` field, so old saves need no migration. If diagnosis selects persistent atmospheric storage instead, it must specify a default for old saves and a persistence test before implementation.

### Knobs and tuning plan

| Knob | Initial value | Basis | Tuning and sensitivity |
|---|---:|---|---|
| Effective column-wind multiplier | 1.0 | Design choice; sparse 2000 ERA5 sample gives global signed regression 1.06, land 1.36, ocean 1.05, with 15° bands about 0.72–2.83 | Retain 0.5/1/2 as calibration sensitivities; evaluate diagnostic and CFD separately because the sample is only 24 days and the pressure-level integration is approximate |
| Mean depletion time `τp` | 9.1 days | [[Trenberth, 1998]](#ref-trenberth1998) global means | Compare 5/9.1/18 days and spatial ERA5 estimates; tune with global water balance and climate metrics |
| Effective mixing depth / humidity conversion | ERA5 screen: global 1.69–2.45 km (10th–90th percentile); land 1.90–2.63 km; ocean 1.65–2.38 km | Ratio of sampled column water to near-surface vapour density; an effective conversion depth, not a diagnosed physical boundary-layer depth | Carry 1.5–3.0 km as a provisional calibration range; keep latitude-class diagnostics (sample medians 1.81–2.89 km) and test sensitivity before implementation |
| Ocean transfer coefficient | Diagnose, no fixed value selected | Bulk-flux formulation; coefficient needs source and ERA5 check | Fit on tuning data; compare ocean evaporation share and wind/temperature bands |
| Lake, river, and plant source coefficients | Existing values converted to flux only for initial comparison | Existing humidity-equivalent proxies; not physical measurements | Compare source volumes to RFC 0016's water account; separate fitting from held-out evaluation |
| Orographic enhancement | Existing relationship as a candidate only | Existing tunable diagnostic | Sensitivity includes zero enhancement; every increment must debit column water |
| Background humidity floor | 0 | Conservation requirement; no unaccounted source | Track snow/glacier response; any nonzero floor requires an explicit owner |

Use windows A, D, and E and the second, colder source world listed in #771 for diagnosis and tuning only, after checking their exact hashes and current-source compatibility. Keep at least two named, unexposed saves as held-out acceptance cases; preregister them and their hashes before any tuning or model replay. Do not use the issue's rough values as observations.

## Diagnosis gates

Run these only after the project owner accepts this Draft RFC. Record scripts, source commit, input hashes, outputs, and decisions in `analysis/rfc0017/README.md`; keep saves and full outputs in private captures.

1. **G0 — Reproduce and inventory.** On the exact #771 source worlds/checkpoints, reproduce current `measure_climate_fields.py` outputs and the RFC 0016 G1b result where its inputs and pinned code permit. Trace every humidity, rain, evaporation, and runoff writer from climate through hydrology, vegetation, ice, and RFC 0016. Establish area-weighted source/sink totals by ocean, lake/river, and vegetation, and verify the reported 38–51% land-rain rise. If it does not reproduce, stop and revise the premise before choosing parameters.
   **G0 result, 2026-11-14:** The pinned G1b rain increases reproduce on the converted checkpoints, and the current-source inventory is recorded in [the G0 report](../analysis/rfc0017/README.md#g0-assessment). The owner accepted the documented limitations on 2026-11-14: converted inputs with matching checked outputs, the source snapshot without its incompatible cached land/ice summary, and unavailable physical atmospheric source/closure totals satisfy G0 for proceeding to G1 and G2. These limits do not establish atmospheric closure.
2. **G1 — Literature and unit screen.** Verify reference equations and units from the papers, obtain ERA5 column water, surface humidity, surface wind, vertically integrated vapour transport and precipitation at a documented period/resolution, and compare them on the same grid. Select the mixing-depth and wind-flux mapping ranges. If ERA5 cannot distinguish them, retain the ranges as calibration uncertainties; do not present them as observed constants.
   **G1 result, 2026-11-15:** The literature and units are checked in [the G1 report](../analysis/rfc0017/README.md#g1-era5-sample-and-calibration-screen). A reproducible WeatherBench2 ERA5 probe sampled twelve 48-hour windows around month starts in 2000 at 1.5° (24 days, eight six-hourly samples per window). It found an area-weighted global column-water/precipitation ratio of 8.10 days, versus 8.28 days over land and 8.05 over ocean. The pressure-integrated column water matches the supplied total-column field within 1.3% in global, land and ocean means. The approximate IVT-to-column-water-times-10m-wind signed regression is 1.06 globally, 1.36 over land, and 1.05 over ocean; latitude-band fits range from 0.72 to 2.83. The effective depth's global 10th–90th percentile is 1.69–2.45 km, with land 1.90–2.63 km and ocean 1.65–2.38 km. Use these only to set broad calibration ranges, not as physical constants or climatology. Keep the 0.5/1/2 wind sweep and use 1.5–3.0 km as a provisional effective-depth sensitivity range; examine latitude-dependent values separately. The pressure integration has 13 levels, surface/top extrapolation, and native-grid dewpoint sampled at coarse-grid centers, so these ranges remain provisional pending broader sampling and tuning-set checks.
3. **G2 — Conservative offline prototype.** Implement only a diagnostic script or isolated prototype first. On analytic source/transport cases, test closed-sphere global closure, local availability, source attribution, orographic removal, zero sources, cap behavior, pole/seam handling, and climate-density convergence. For option B, verify that transitions follow local arrows and test convergence, divergence, calm cells, closed loops with rainout, and closed classes without a sink. Compare option A with option B for runtime and closure. If neither can close to `1e-8` relative to total source in double precision and remain positive, revise the discretization before production code.
   **G2 result, 2026-11-15:** The isolated prototype in [the G2 analysis](../analysis/rfc0017/README.md#g2-preliminary-analytic-prototype-screen) compares an upwind finite-volume face-flux scheme with local wind-arrow transitions. Both close below `1.3e-16` globally and `1.6e-17` per-cell on an 8×16 ocean-source case, remain nonnegative, preserve source-tag linearity, produce land rain only from transported ocean moisture, and handle zero source, calm flow, convergence/divergence, and a closed loop with positive rainout. A localized enhanced sink raises local rain and global closure remains within `5.7e-16`; a one-cell finite-interval check caps rain to opening mass. The prototype also checks longitude wrap, latitude bounds, wind-arrow alignment, sinkless closed-class detection, and uniform-source invariance from 4×8 to 16×32. At 90×180 (climate density 1), both schemes remain positive and close globally within `1.3e-16` and locally within `6.0e-18`. Separate uncontended sparse-solve measurements are 0.02865 s for finite volume and 0.02824 s for arrows, effectively equal at this scale. Peak RSS and density-4 (360×720) scaling remain unmeasured; the one-cell cap check is not a time-integrated transport update. G2's analytic and density-1 prototype gate passes.
   **Owner review, 2026-11-15:** The owner approved moving to G3 after reviewing the G2 results. G3's numerical climate-fidelity and runtime limits are preregistered in [the G3 plan](../analysis/rfc0017/g3-preregistration.md) before any G3 comparison.
4. **G3 — Rainfall patterns, budget, and runtime.** Use the registered tuning windows and, after the candidate and checks are frozen, the registered held-out saves. Compare both wind modes and climate densities 1 and 4. Evaluate wet belts, rain shadows, dry interiors, and cold-land precipitation using the [owner-approved regional checks](../analysis/rfc0017/g3-preregistration.md). Continue reporting legacy precipitation RMS, latitude-band means and Köppen shares as diagnostics, along with land/ocean rain shares, source amounts, transport, column water, closure and climate-update time. The two held-out saves have been evaluated for the frozen `p = 0.25` ridge-rainout candidate; the dry-interior and wind-arrow held-out land/ocean-share limits fail. The density-4 offline candidate-update runtime proxy is below 2×, while production integration cost remains unmeasured. **G3 failed after the owner-authorized three tuning rounds. By explicit owner exception dated 2026-11-18, RFC 0019 may continue to G4 and implementation planning with this known modeling simplification. This exception waives the G3 sequencing block only; it does not mark G3 as passed or select a tuning-only candidate for production.** The inland/coast pattern remains unmet and approximate on synthetic worlds; see the [three-round tuning record](../analysis/rfc0017/g3-three-round-plan.md#three-round-disposition-g3-failed-by-owner-exception).
5. **G4 — Downstream coupling.** After G3 passes or the owner grants an explicit exception, rerun the eligible RFC 0016 coupling screen under the atmospheric version selected for implementation. Report land rain, lake supply/volume, ocean export, frozen precipitation share, and ledger changes against its control. Recheck glacier accumulation, erosion, and biome distributions. This gate informs RFC 0016; it does not itself implement land-water storage or RFC 0015. The 2026-11-18 exception permits continuing to G4 while G3 remains failed; it does not select a tuning-only candidate.

   **G3 tuning screen, 2026-11-15:** The preregistered proxy-to-source hypothesis closes the stationary water budget (maximum relative residual below `3e-15`) and is resolution-convergent (density-1/4 global rain changes by at most 0.10%). The tuning screen indicates it cannot meet the held-out climate-pattern target and it exceeds the density-4 runtime limit. Even after fitting one uniform precipitation scale, wind-arrow precipitation rel-RMS is 0.669–0.695 at density 1 and 0.897–0.911 at density 4; the density-4 climate update takes 2.97–3.65× current-source time at the selected parameters. A complete window-A bracket sweep found no setting below 0.50 best-fit RMS or 2× runtime. These are tuning-screen results, not held-out acceptance measurements. The two acceptance-only saves remain unopened, G3 is incomplete, and G4 is not authorized. See the [G3 tuning report](../analysis/rfc0017/README.md#g3-tuning-screen-initial-source-mapping-does-not-pass). The source mapping and/or transport approach needs a design revision before held-out evaluation.

   **G3 source, rainout, and steady-formulation screens, 2026-11-15:** A shared land-proxy multiplier calibrated over the tuning set to the broad 13% land-origin source scale check did not materially improve density-4 spatial RMS (0.896–0.922). A second screen replaced the ocean humidity target with a bulk humidity-deficit flux, `E = rho_air * C_E * |u| * max(q_sat(SST) - q_air, 0)`, while retaining the scaled land proxies. Across transfer coefficients 0.001/0.003/0.005 and land-proxy multipliers 1/5/10, the best density-4 wind-arrow RMS after a uniform amplitude fit was 0.889, still far above 0.50. At the broad-source-share setting (`C_E=0.003`, land multiplier 5), the ocean source share ranged 81–89% and density-4 climate-update time remained 2.8–3.75× baseline. The screen used the legacy climate humidity as ambient `q_air`; this is not the candidate stationary column humidity and remains a limitation. A follow-up used G1's 24-day class-mean column-water/precipitation ratios as a latitude- and surface-dependent rainout time; best-fit RMS worsened to 0.766–0.805 at density 1 and 0.934–0.946 at density 4, with density-4 runtime 3.08–3.93× baseline. At the owner's direction, a one-day excess-water rain sink was then screened on tuning windows A, D, and E at density 4. With 80% and 60% saturation thresholds, best-fit RMS remained 0.898–0.907 and update time was 17.6–21.2× baseline; closure was exact at displayed precision. This rejects the tested thresholded closure, not every steady formulation. No source, rainout, or thresholded-rain candidate is selected. The tuning evidence still fails the climate-pattern and runtime criteria. Do not open held-out saves. Next constrain a steady humidity/transport field against independent column-water and vertically integrated vapor-transport data; if a steady closure still misses the registered limits, return for a prognostic-atmosphere design. Because atmosphere CFD advances one simulated day per tectonic step, that path needs an explicit multi-day spin-up cadence and runtime estimate. G3 remains incomplete and G4 is not authorized. See the [source-mapping and rainout report](../analysis/rfc0017/README.md#g3-source-mapping-revision-screen) and the [steady-formulation screen](../analysis/rfc0017/README.md#g3-steady-formulation-screen-saturation-triggered-rainout).

   **G3 local-lift screen, 2026-11-15 (#913):** One preregistered bounded sink based on smoothed positive terrain ascent was compared with the uniform sink on tuning windows A and cold at density 1, using finite-volume transport and diagnostic and CFD winds. It improved windward/leeward contrast and modestly reduced best-fit precipitation RMS in all four pairs, with lower Köppen dry-group share. It also reduced ≥60° land rain in all four pairs (9–53% on A and 13–38% on cold). The preregistered cold-land stop condition is met; this local-lift law is not selected for further tuning or held-out evaluation. This result does not reject conservative moisture accounting or every orographic formulation. G3 remains incomplete; acceptance thresholds and runtime limits are unchanged. See the [issue #913 report](../analysis/issue913/README.md) and [lift preregistration](../analysis/issue913/lift-preregistration.md).

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/climate.py`: `compute_humidity`, `compute_precipitation`, `_land_moisture_source`, `_vegetation_transpiration_source`, `_orographic_retained_fraction`, and `compute_moisture_flux_convergence`; likely replace or retire `_humidity_zonal_sweep`, `_humidity_meridional_sweep`, `_retention_factor`, `_zonal_base_retention`, and `_meridional_base_retention`. Add a conservative column-water solver and named unit conversions.
- `backend/app/hydroclimate/atmosphere_cfd.py`: supply its current wind through the same `ClimateFields` interface; do not add moisture state in option A. Diagnostic wind remains through `climate.py`'s existing wind path.
- `backend/app/hydroclimate/erosion.py`, `hydrology.py`, `biomes.py`, RFC 0015 glacier code, and RFC 0016 land-water code: audit precipitation consumers and source ownership; avoid independent water writers in climate.
- Tests: add focused coverage in `backend/unit_tests/test_climate.py` for cell/global closure, available-water limits, source totals, cap conversion, source-free transport, both wind modes, pole/seam faces, resolution changes, and persistence compatibility. Add only a stress test if the production solve cost or coupled long-run behavior requires it.
- `docs/simulation-model/climate.md`, `docs/simulation-model/hydrology.md`, and `docs/architecture.md`: update the climate field contract, physical source units, diagnostic/prognostic boundary, and request flow after implementation.
- `bin/debug/measure_climate_fields.py`: reuse for the existing metrics; add an RFC-specific ledger/ERA5 diagnostic only if the current script cannot provide the required measurements.

## Assumptions and uncertainty

- A steady annual-mean moisture budget is an adequate representation for the simulator's geological timestep and current diagnostic precipitation consumers. It omits transient weather and atmospheric storage changes between climate calls.
- One effective wind transports the whole water column. Variation in wind and humidity with height, and variations within the averaging period, can change the strength and direction of moisture delivery. Their contribution to the tuning-world errors is unmeasured. Retain an explicit calibration range and track improvements to vertical transport in [#944](https://github.com/adubey/mantle-bloom/issues/944).
- A first-order rainout time is a bulk closure law, not cloud microphysics. The nine-day global scale does not set local rainout, convective efficiency, or storm tracks.
- Ocean evaporation, lake/river evaporation, and transpiration are surface fluxes with different drivers. Existing g/kg-equivalent source formulas do not constrain their physical rates without additional evidence.
- A biome-weighted recycling fraction is an acceptable provisional transpiration proxy when a named surface reservoir pays for the water. Vegetation contributed 1.4–1.7% of total source in the A and cold control snapshots, so direct plant supply explains little of the measured rain concentration at those settings. The runs do not quantify coupled feedback. Improved transpiration physics is tracked in [#944](https://github.com/adubey/mantle-bloom/issues/944).
- Land ET can cycle water repeatedly and therefore gross precipitation can exceed oceanic input locally and over a chosen interval. Conservation means each evaporation, atmospheric transport, precipitation, and surface storage transfer has one owner; it does not prohibit recycling.
- ERA5 is a modern-Earth comparison, while synthetic topography, wind, and land cover vary. Use it to constrain scales and patterns, not as a universal geological forcing.
- Changing rain affects erosion, lakes, Köppen classes, and glacier accumulation. Their movement is an impact to measure, not a reason to introduce a hidden rain floor.

## Acceptance criteria

- [x] G0 reproduces its pinned evidence and records the owner's accepted limitations on converted input bytes, the skipped source cache summary, and unavailable physical atmosphere closure.
- [x] G1 documents verified units, source URLs/chunk hashes, and provisional calibration ranges; the sparse ERA5 sample and integration/regridding uncertainty are stated before model calibration.
- [x] G2 prototypes demonstrate nonnegative water, local availability limits, no source-free land rain beyond transported ocean moisture, and per-cell/global water closure within `1e-8` relative to total source for analytic and climate-density-1 grid cases.
- [ ] With a zero land source, total land precipitation is no greater than ocean-origin moisture entering land, within ledger tolerance.
- [ ] A saturation cap and each precipitation enhancement convert their removed moisture to rain exactly once; there is no unaccounted humidity floor or ocean reset.
- [ ] Area-weighted closure holds at climate densities 1 and 4; global rainfall changes by no more than 5% relative and land/ocean rainfall shares by no more than 5 percentage points between densities on registered cases.
- [ ] Meet the owner-approved held-out limits: equatorial ocean rain-belt center-to-flank ratio ≥1.5; at least 60% correctly directed ridge crossings with positive median contrast gain; inland (>1,500 km) mean rain ≤0.6× the near-coast (<500 km) mean; and positive precipitation on exposed cold coast, with no minimum cold-inland rate. Preserve the held-out land/ocean share bound, closure, density convergence, and runtime limits. Report legacy precipitation RMS, latitude-band means, and Köppen shares as diagnostics. See the [G3 preregistration](../analysis/rfc0017/g3-preregistration.md) and [held-out result](../analysis/rfc0017/g3-coarse-transport-plan.md#held-out-pattern-screen-for-frozen-ridge_025-candidate).
- [ ] Effective column water stays physically plausible under the stated density and depth assumptions. Report the upper tail and spatial concentration of water and precipitation. Any additional rain used to limit accumulation debits the atmospheric budget exactly once.
- [ ] Both diagnostic and CFD wind complete the same budget; report step-time change and keep median density-4 climate-update time within the owner-approved 2× current-source budget.
- [ ] RFC 0016's rerun reports the coupled and control changes with closed land and atmosphere ledgers; no treatment is used to tune G0–G3.
- [ ] Old saves load; option A adds no persistent state or save-size cost.

## Falsification criteria

- If the RFC 0016 rain multiplication does not reproduce on its pinned inputs, the specific coupling premise is false or stale; revise the RFC around the defect that does reproduce.
- If an area-conservative column solver cannot maintain positive water and close the finite-volume ledger at the specified tolerance on a closed sphere, reject the proposed discretization and try the streamline alternative before production work.
- If the only parameters that pass the ledger make precipitation patterns or downstream hydrology unusable on held-out worlds, the assumption that a steady diagnostic is sufficient is falsified. Reopen the choice between a more complete prognostic atmosphere and a simpler explicitly calibrated proxy.
- If the diagnostic and CFD winds imply irreconcilable transport mappings, do not share one fitted multiplier; revise the design to mode-specific transport or narrow the supported wind modes.
- If plausible rainout and transport sensitivities produce order-of-magnitude differences in land rain or lake supply, the current literature and diagnostics do not constrain this design enough for implementation; gather new evidence or withdraw the recommendation.

## Revision note

**2026-11-18 — G3 failed after three tuning rounds; owner exception to continue.** The
owner authorized at most three carefully selected tuning rounds, with an explicit failure
and exception if no eligible candidate passed. The best eligible tuning candidate,
`alpha=1, L=1,500 km`, had a worst inland/coast rain ratio of `1.515` on the cold case with
CFD winds, against the `0.6` limit. The more aggressive `alpha=2` setting reduced that ratio
to `1.300` but retained only `77.5%` of baseline inland rain, below the `80%` safeguard, so
it was ineligible. No candidate was frozen or evaluated on new held-outs. The owner explicitly
waived the G3 sequencing block so RFC work may continue to G4 and implementation planning
with this known approximation. G3 remains formally failed; the exception does not claim that
rain shadows, dry interiors, or other regional patterns passed. Results and fingerprints are
in the [three-round tuning record](../analysis/rfc0017/g3-three-round-plan.md#three-round-disposition-g3-failed-by-owner-exception).

**2026-11-17 — G3 owner-approved pattern checks and held-out result.** The owner
replaced the held-out latitude-band and Köppen-share differences (each ≤10
percentage points) and area-weighted precipitation relative RMS (≤0.50) with
measurable rain-belt, rain-shadow, dry-interior and exposed cold-coast checks.
The held-out land/ocean share difference remains limited to 5 percentage
points. The new limits are an equatorial belt ratio ≥1.5, ≥60% correctly
directed ridge crossings with positive median contrast gain, inland/coast rain
ratio ≤0.6, and positive exposed cold-coast precipitation, with no cold-inland
minimum. This change reflects the owner's preference for visible regional
patterns on synthetic worlds over pixelwise agreement with current climate.
Legacy metrics remain diagnostics. On the two held-out
saves, the frozen `p = 0.25` candidate passed rain-belt, shadow, cold-coast,
budget-closure and density-convergence checks. The inland/coast limit failed in
5 of 16 wind/solver/density combinations, all under CFD winds. Wind-arrow
land/ocean shares also exceeded the retained held-out limit at density 4 on
both saves. A density-4 offline candidate-update timing proxy on the four
registered tuning worlds measured 1.22–1.42× current-climate time by mode and
solver, below the 2× limit; production integration cost remains unmeasured.
G3 has not passed; do not tune from the opened held-outs. See the impact table in the
[held-out results](../analysis/rfc0017/g3-coarse-transport-plan.md#held-out-pattern-screen-for-frozen-ridge_025-candidate).

**2026-11-16 — G3 ridge rainout screen.** Intercepting water at a crest
improved paired windward/lee contrasts in both wind modes. The weak arm
still reduced cold inland rain on A CFD by 14%, explained by lower net
inflow into those regions. Turning down the effect in dry air weakened
the shadows and retained the A deficit. The atmospheric budget closed in
all arms. G3 remains open. See the
[ridge comparisons](../analysis/rfc0017/g3-coarse-transport-plan.md#rainout-during-ridge-crossings).

**2026-11-16 — G3 rain belt and terrain screens.** Modulating the local
rainout rate produced little belt contrast. Restoring broad moisture
transport and spreading only saturation overflow produced broader rain
belts with exact water accounting. The fixed-budget fine terrain test
changed median windward/lee differences by less than 1 mm/yr in every
case. A ridge-scale rainout candidate still needs tuning, and the native
density-4 budget, runtime and held-out behavior remain unverified. See the
[tuning report](../analysis/rfc0017/g3-coarse-transport-plan.md#broad-rainout-tuning-screen).

**2026-11-16 — G3 bounded-column and transport screen.** The local saturation
cap closed the budget and removed the vapor excess, but 62–64% of rain came
from immediate overflow and narrow wet stripes remained. Removing the
convergent component of transport reduced overflow to 3–6% with the same cap
and sources. It softened the stripes, improved some cold-land supply, and
left weak rain belts and rain shadows. The next tuning candidate separates
the broad rainout cue and terrain-weighted placement from the corrected
transport. G3 remains open. See the [paired results](../analysis/rfc0017/g3-coarse-transport-plan.md#convergence-adjustment-with-the-local-cap-in-both-modes).

**2026-11-16 — G3 rainfall-pattern scope.** The owner asked RFC 0019 to produce
recognizable rain belts, rain shadows, dry interiors, and cold-land precipitation
at useful approximation for synthetic worlds. The proposed G3 change is to
replace numerical agreement with the current-source rain and biome fields,
including the 0.50 relative RMS limit, with checks of those four regional
patterns. Keep the original comparisons as diagnostics for the size and location
of downstream changes. Retain closure, nonnegativity, availability, resolution,
the 2× density-4 runtime limit, and G4 downstream checks. Freeze measurable
regional checks with the owner before evaluating the held-out saves.

The [latest A and cold tuning screen](../analysis/rfc0017/g3-coarse-transport-plan.md#targeted-convergence-result)
reduced excessive convergence, but regional water delivery moved in both
directions and the effective water column still exceeded its assumed air mass
in one cold CFD snapshot. Stop tuning that setting. The next candidate should
keep the coarse conservative solve, bound accumulation using a budgeted sink,
and test regional delivery before spending time on fine rain placement. Keep
the single-layer wind and transpiration simplifications documented. The more
complete physics is a deferred [follow-up in #944](https://github.com/adubey/mantle-bloom/issues/944).
G3 remains incomplete and G4 follows G3.

**2026-11-16 — G3 coarse transport plan.** After reviewing #913 and the
exploratory comparison in [#923](https://github.com/adubey/mantle-bloom/issues/923),
the owner chose coarse transport, conservative mixing, and bounded rain placement
as the next design to plan. The #923 full-resolution run stopped without a result;
it is a supporting diagnostic and does not measure either wind solver's error.
G0–G2 stand on their original evidence; new operator and allocation checks
extend G2 before the new G3 screen. G3 limits and held-out cases remain as
registered, and G4 depends on G3.

**2026-11-14 — G0 interpretation ([PR #903](https://github.com/adubey/mantle-bloom/pull/903)).** The original gate required the exact #771 source
worlds, current climate-field output and area-weighted physical source/sink totals.
The owner accepted G0 on the documented evidence: the converted checkpoints reproduce
the checked G1b quantities exactly; the source snapshot omits an incompatible cached
land/ice summary; and current atmospheric source/closure volumes cannot be computed
from the humidity proxies. The [G0 report](../analysis/rfc0017/README.md#g0-assessment)
preserves those limits. G1 and G2 may proceed; this does not establish atmospheric
closure or change the later gate thresholds.

**2026-11-15 — G3 preregistration.** After reviewing G2, the owner authorized G3 and
approved its held-out climate-fidelity and runtime limits before any G3 comparison. The
limits, registered tuning/held-out separation, and initial humidity-proxy-to-flux hypothesis
are recorded in the [G3 preregistration](../analysis/rfc0017/g3-preregistration.md).

**2026-11-15 — G3 tuning screen.** The initial mapping passed stationary closure and
density convergence but its tuning results indicate a miss on precipitation-pattern fidelity
and density-4 runtime. The held-out saves remain unopened; G3 has not passed and no
implementation or G4 work should begin until the owner selects a revised design direction.
See the [tuning report](../analysis/rfc0017/README.md#g3-tuning-screen-initial-source-mapping-does-not-pass).

**2026-11-15 — G3 source and rainout revision screens ([PR #903](https://github.com/adubey/mantle-bloom/pull/903)).** At the owner's direction, the tuning-only source mapping was revised to test a pooled land-proxy scale and a bulk ocean humidity-deficit flux. The initial land proxy supplied only 2–4% of total source; scaling it to approach the RFC's broad land-source scale did not materially improve density-4 spatial error. The bulk ocean flux used legacy near-surface humidity as ambient humidity and produced an 81–89% ocean share at the selected coefficient, but density-4 amplitude-adjusted spatial RMS remained at least 0.889 and runtime exceeded budget. A class-dependent rainout sensitivity derived from G1 also worsened the best-fit RMS at both densities. Thus the source mapping and first rainout revision are insufficient; held-out saves stay unopened. The numerical criteria and G3 gates did not change. Results and limitations are in the [source-mapping and rainout report](../analysis/rfc0017/README.md#g3-source-mapping-revision-screen).

**2026-11-15 — G3 steady-formulation screen ([PR #903](https://github.com/adubey/mantle-bloom/pull/903)).** The owner directed the investigation to revisit the steady formulation. A conservative excess-water rain sink at 60% and 80% saturation closed its budget but left density-4 best-fit precipitation RMS at 0.898–0.907 and increased update time to 17.6–21.2× baseline on three tuning windows. This specific closure fails G3; it does not prove that all steady closures fail. The held-out saves remain sealed. The next diagnostic will constrain steady humidity and transport with independent column-water and vertically integrated transport data. A prognostic alternative needs an explicit multi-day spin-up design because the existing atmosphere advances one simulated day per tectonic step. G3 is incomplete and G4 is not authorized. See the [steady-formulation screen](../analysis/rfc0017/README.md#g3-steady-formulation-screen-saturation-triggered-rainout).

**2026-11-15 — G3 local-lift screen (#913).** The preregistered 250 km smoothed ascent sink improved coarse-scale RMS and windward/leeward contrast on both tuning saves at density 1, but decreased ≥60° land rain in diagnostic and CFD modes on both saves. Stop this local-lift formulation under the preregistered cold-land criterion. G3 remains incomplete, held-out saves remain unopened, and no G3 limit changes. Findings and private capture hashes are recorded in the [issue #913 report](../analysis/issue913/README.md).

## References

- <a id="ref-seager2010"></a>
Seager, R., Naik, N., & Vecchi, G. A. (2010). Thermodynamic and dynamic mechanisms for large-scale changes in the hydrological cycle in response to global warming. *Journal of Climate*, 23, 4651–4668. [doi](https://doi.org/10.1175/2010JCLI3655.1)

- <a id="ref-trenberth1998"></a>
Trenberth, K. E. (1998). Atmospheric moisture residence times and cycling: Implications for rainfall rates and climate change. *Climatic Change*, 39, 667–694. [doi](https://doi.org/10.1023/A:1005319109110)

- <a id="ref-okikanae2006"></a>
Oki, T., & Kanae, S. (2006). Global hydrological cycles and world water resources. *Science*, 313(5790), 1068–1072. [doi](https://doi.org/10.1126/science.1128845)

- <a id="ref-vanderent2010"></a>
van der Ent, R. J., Savenije, H. H. G., Schaefli, B., & Steele-Dunne, S. C. (2010). Origin and fate of atmospheric moisture over continents. *Water Resources Research*, 46, W09525. [doi](https://doi.org/10.1029/2010WR009127)

- <a id="ref-eumetsat2011"></a>
EUMETSAT Satellite Application Facility on Climate Monitoring. (2011). *Product User Manual: SSM/I products, release 3.2*, section 3.1.3, “Latent heat flux / Evaporation”. [PDF](https://www.geo.uni-hamburg.de/en/service/icdc/data/atmosphere/docs-atmo/saf-cm-dwd-pum-hoaps3-2-1-1.pdf)
