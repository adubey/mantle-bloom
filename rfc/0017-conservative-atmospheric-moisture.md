# RFC 0017: Conserve atmospheric moisture in diagnostic rain and transport

**Status:** Draft
**Date:** 2026-11-03
**Related:** [#771](https://github.com/adubey/mantle-bloom/issues/771), RFC 0015, RFC 0016; issues #522, #526, #530, #358

## Problem

The climate estimates humidity with two passes across its grid: one along longitude (east-west) and one along latitude (north-south). Over land, the model carries less moisture farther along the wind path. Windward slopes remove more, and the model tracks that slope-related amount separately. The ordinary loss as moisture crosses land has no recorded destination. After the passes, the model calculates rain from the remaining humidity at each cell and adds rain on slopes and where winds flow together. This later rain is not subtracted from moisture carried downwind; only the slope-related loss reduces the carried amount. The model also discards moisture above a limit, adds a small minimum amount over dry land, and resets ocean humidity based on water temperature without recording evaporation as added water. These steps add or remove water without accounting for it, so the model cannot balance water entering the atmosphere with rain, other losses, and changes in atmospheric storage.

This violates **design principle 6, account for material through every process**: the same atmospheric water can rain repeatedly without replenishment, and other water is silently removed or introduced at caps and clamps. It also strains principle 5, explicit and consistent units: near-surface specific humidity in g/kg is used as a surrogate for a moisture source and annual precipitation through a fitted conversion.

The downstream effect is measured. RFC 0016's [G1b result](../analysis/rfc0016/README.md#g1b-result) reports that coupling plant transpiration to humidity raised land rain by 38–51% and increased water exported to the ocean. The report attributes that multiplication to rain not being debited from humidity. Those are measurements on its pinned source and windows, not a forecast of a conserving model's result. G1b cannot evaluate the intended land-water coupling until this RFC's question is resolved.

The scope is the annual-mean diagnostic water cycle used by erosion, hydrology, biomes, and ice. It is not a weather model and does not resolve clouds, convective storms, seasons, or vertical atmospheric structure. The proposal must preserve a useful long-term precipitation field without claiming those omitted processes.

## Proposed solutions

### A. Conservative two-dimensional column-water budget (recommended for diagnosis)

Represent atmospheric water as column-integrated vapour, `W` in kg/m². Represent evaporation and precipitation as mass fluxes, `E` and `P` in kg/m²/s, and horizontal vapour transport as `F` in kg/(m·s). On the spherical climate grid, solve the finite-volume balance

`d(A_i W_i)/dt = A_i(E_i - P_i) - net outward edge flux_i`,

where `A_i` is the spherical area of cell `i`. Use wind to transport moisture, and remove water from the column through explicit precipitation sinks. Each flux crosses a shared cell edge once, so transport cancels exactly in the global sum. The recommended first implementation is a steady climatological solve: the temporary atmospheric inventory reaches a converged field on each climate update, so its net storage change is zero. Keep the inventory and iteration ledger for the solve; it need not become saved `World` state. A future prognostic atmosphere would need a separately justified physical time interval and save-compatibility plan.

Use a positive, bounded precipitation sink linked to column water, initially `P = W / τp`, and test a saturation-triggered or orographic enhancement only as a sink that removes available water. Wind convergence is already in the finite-volume flux divergence; do not add it a second time as unbudgeted rain. A source tag (ocean, lakes/rivers, vegetation) is useful for diagnostics, but total water—not each tag's separate transport—is the closure requirement.

This option replaces the two blended one-dimensional sweeps with a low-order, conservative two-dimensional solve. It can represent meridional and zonal transport together, wrap longitude, close fluxes at the poles, and use both diagnostic and CFD winds. Its uncertainties are the effective column-mean transport speed, rainout law, and how annual-mean evaporation is represented.

### B. Conservative source-tagged streamline sweeps

Retain directional sweeps for low runtime and familiar behavior, but carry source-tagged water mass rather than a humidity-like scalar. Every reduction in carried water becomes precipitation or a named boundary export; local evaporation is an explicit source. The zonal and meridional transports would need a conservative joint rule, because blending two independently computed humidity fields is not itself a mass-conserving transport. This is a useful prototype and fallback if the two-dimensional solve is too costly, but it does not naturally represent cross-flow transport or wind convergence.

### C. Prognostic moisture inside the atmosphere CFD

Add column vapour, evaporation, condensation, and precipitation to `AtmosphereCFDState`, and advance them with atmospheric flow. This makes storage explicit but couples annual-mean surface hydrology to a solver that currently advances one simulated day per tectonic step. It also adds persistent state, save migration, CFL/runtime questions, and a substantial calibration problem. The retired CFD humidity/precipitation path produced near-zero, uncalibrated rainfall and starved erosion and hydrology; adding state alone would repeat that failure mode. Revisit only if the steady diagnostic approaches fail a named criterion and a longer-term atmosphere/hydrology coupling is designed.

### Recommendation

Take option A through diagnosis, with option B as the bounded-cost fallback. Do not restore the retired CFD humidity equation as-is. The key difference is the conserved extensive column quantity and shared-edge flux accounting, with each water source and rain sink in physical units. The diagnosis gates below decide whether this representation retains useful climate structure and is affordable before an implementation is authorized.

## Design

### Budget and transport

For every climate cell `i`, use area `A_i` in m² and column water `W_i` in kg/m². Let `E_i` and `P_i` be surface-to-air and air-to-surface fluxes in kg/m²/s. Let `F_ij` be signed vapour mass transport across a shared edge in kg/s. Over solve interval `Δt`, audit

`A_i(W_i,new - W_i,old) = Δt[A_i(E_i - P_i) - Σ_j F_ij]`.

The global residual is `Σ A_i W_i,new - Σ A_i W_i,old - Δt Σ A_i(E_i-P_i)`; internal edge fluxes cancel, leaving only documented open boundaries (none on the closed sphere). The steady solution must have `W_i,new = W_i,old` within solver tolerance and global evaporation equal to precipitation. Report both per-cell normalized residuals and the global residual. Do not count solver relaxation as physical storage change.

Use a conservative donor-cell flux on the spherical climate grid. Cell area is calculated from latitude-edge bands and longitude width, rather than treating every latitude row as equal area. Longitude wraps; polar faces have zero area; east/west fluxes use the same value with opposite signs in their two adjacent cells. Clip negative numerical values to zero only with their clipped mass added to the local precipitation/sink ledger; otherwise fail the solve. Preserve the climate-density meaning of `World.climate_density`: cell fluxes and totals are area-weighted and should converge with resolution.

Use `surface_wind` from either wind mode for the first prototype. It is in m/s and targets a 10 m wind; the CFD solver currently reaches only 0.42–0.66 of its target on the reported diagnosis. The moisture-bearing column wind is not observed by this model. Start with an effective transport multiplier of **1.0 (design choice)** and a separately named, dimensionless multiplier fitted only on the tuning set by comparing model flux against ERA5 vertically integrated vapour transport resampled to this grid. Sweep 0.5, 1.0, and 2.0 as sensitivities; report diagnostic and CFD wind separately. Do not call the product of column water and 10 m wind a measured atmospheric flux.

### Rain and the atmospheric column

Start with a first-order sink `P = W / τp`, with `τp = 9.1 days` as a **literature-constrained initial global scale**, from the ratio of 24.6 mm global mean precipitable water to 2.69 mm/day precipitation in Trenberth (1998). This is an Earth global-mean sanity check, not a per-cell law or a target precipitation climatology. Test a 5-day and 18-day sensitivity and compare latitude- and surface-class-dependent depletion times inferred from ERA5 before selecting any spatial variation. The sink must be nonnegative and cannot remove more water than is present in the finite-volume update. A cell's rain plus exported flux cannot exceed its opening column inventory and incoming flux plus evaporation over the solve interval.

The first recommendation has no hard humidity floor or vapor cap. Saturation is represented by routing excess column water to precipitation, with the local excess debited exactly once. If the atmospheric layer cannot support that treatment at the selected grid and timestep, propose a named temporary reservoir or source/sink and show its budget. Keep `humidity_g_kg` as the near-surface interface for existing readers by deriving it from `W` using an explicit effective mixing depth and air density; diagnose the conversion against observed precipitable water and ERA5 near-surface humidity. A maximum near-surface humidity remains a diagnostic bound, not a place to discard water.

Orographic lift may increase the precipitation sink on windward slopes only when it removes moisture from `W`; the corresponding rain is credited to `P`. Moisture-flux convergence and subsidence emerge from the divergence of `F`. Any extra convective efficiency or subsidence parameter must be bounded, have units, and repartition available water rather than create or delete it. Coherent climate noise may perturb transport or the rainout rate, but the field must remain nonnegative and its area-integrated source and sink effects must be included in the ledger.

### Physical surface sources and ownership

- **Ocean evaporation:** replace the humidity reset with a surface evaporation flux in kg/m²/s. Prototype a bulk aerodynamic flux from sea-surface saturation humidity, near-surface air humidity, surface wind speed, air density, and a transfer coefficient. Bound evaporation to available ocean water. Diagnose the coefficient from ERA5 and test sensitivity; the model's ocean temperature is a diagnostic and must be used at its actual cell value.
- **Lakes and rivers:** convert the current depth-based humidity proxies into explicit annualized mass fluxes. Freeze or dry water sources contribute zero. Use only liquid water the surface owns; lake/river evaporation already removed by `hydrology.py` must not be charged again to the atmosphere. The gate must inventory the current `RIVER_EVAPORATION_*` and `LAKE_EVAPORATION_*` terms before choosing one owner.
- **Vegetation:** express transpiration as a mass flux. RFC 0016's current 0.16 recycling fraction of prior precipitation may be a starting calibration hypothesis, but it is not a measured ET law. Convert the prior rain from mm/yr to kg/m²/s, apply the biome coefficient once, and ensure RFC 0016's land-water store or another explicit surface reservoir owns the debit. When its water source is not available, transpiration is zero rather than an unbooked input.
- **Background vapour:** remove `MIN_BACKGROUND_HUMIDITY_G_KG` as an independent source. If a trace precipitation floor is still required for ice accumulation, make it an explicit transfer from a named surface reservoir or retain no floor and measure the glaciation impact.

[[Oki & Kanae, 2006]](#ref-okikanae2006) report that the ocean supplies about 86% of global evaporation and receives about 78% of global precipitation. Use these as broad Earthlike scale checks alongside land/ocean rain share; do not force local cells or each synthetic world to those fractions. [[van der Ent et al., 2010]](#ref-vanderent2010) estimate 40% of terrestrial precipitation originates from land evaporation in their ERA-Interim accounting. This supports source-tagged recycling diagnostics, not an exact per-cell transpiration fraction.

### Scales, outputs, and compatibility

[[Trenberth, 1998]](#ref-trenberth1998)'s global-mean column water and precipitation give a depletion time near nine days. A deliberately rough transport check using 25 kg/m² of column water and 5 m/s gives about 125 kg/(m·s) of moisture flux. At 1,000 mm/yr of rain (about 3.2×10⁻⁵ kg/m²/s), a one-dimensional flux-to-sink ratio is about 4,000 km. This is only a consistency check: it assumes uniform flow, a single transport direction, and no replenishment, and is not a fitted inland decay distance. It shows why a 380 km half-life should not be retained as the rain law without a moisture budget.

Convert `P` to the existing annualized precipitation field in mm/yr at one named interface. Resample it to terrain nodes with area-aware weights when integrating totals; do not sum nearest-neighbor grid values as if their areas were equal. Preserve the existing `ClimateFields` and `climate_cache` readers where possible. The recommended solve adds no persistent `World` field, so old saves need no migration. If diagnosis selects persistent atmospheric storage instead, it must specify a default for old saves and a persistence test before implementation.

### Knobs and tuning plan

| Knob | Initial value | Basis | Tuning and sensitivity |
|---|---:|---|---|
| Effective column-wind multiplier | 1.0 | Design choice; 10 m wind is only a proxy | Tune 0.5/1/2 against ERA5 integrated moisture transport on tuning seeds; evaluate diagnostic and CFD separately |
| Mean depletion time `τp` | 9.1 days | [[Trenberth, 1998]](#ref-trenberth1998) global means | Compare 5/9.1/18 days and spatial ERA5 estimates; tune with global water balance and climate metrics |
| Effective mixing depth / humidity conversion | Diagnose, no fixed value selected | Not resolved by the present model | Derive from column water and near-surface humidity using ERA5; report uncertainty and sensitivity before implementation |
| Ocean transfer coefficient | Diagnose, no fixed value selected | Bulk-flux formulation; coefficient needs source and ERA5 check | Fit on tuning data; compare ocean evaporation share and wind/temperature bands |
| Lake, river, and plant source coefficients | Existing values converted to flux only for initial comparison | Existing humidity-equivalent proxies; not physical measurements | Compare source volumes to RFC 0016's water account; separate fitting from held-out evaluation |
| Orographic enhancement | Existing relationship as a candidate only | Existing tunable diagnostic | Sensitivity includes zero enhancement; every increment must debit column water |
| Background humidity floor | 0 | Conservation requirement; no unaccounted source | Track snow/glacier response; any nonzero floor requires an explicit owner |

Use windows A, D, and E and the second, colder source world listed in #771 for diagnosis and tuning only, after checking their exact hashes and current-source compatibility. Keep at least two named, unexposed saves as held-out acceptance cases; preregister them and their hashes before any tuning or model replay. Do not use the issue's rough values as observations.

## Diagnosis gates

Run these only after the project owner accepts this Draft RFC. Record scripts, source commit, input hashes, outputs, and decisions in `analysis/rfc0017/README.md`; keep saves and full outputs in private captures.

1. **G0 — Reproduce and inventory.** On the exact #771 source worlds/checkpoints, reproduce current `measure_climate_fields.py` outputs and the RFC 0016 G1b result where its inputs and pinned code permit. Trace every humidity, rain, evaporation, and runoff writer from climate through hydrology, vegetation, ice, and RFC 0016. Establish area-weighted source/sink totals by ocean, lake/river, and vegetation, and verify the reported 38–51% land-rain rise. If it does not reproduce, stop and revise the premise before choosing parameters.
2. **G1 — Literature and unit screen.** Verify reference equations and units from the papers, obtain ERA5 column water, surface humidity, surface wind, vertically integrated vapour transport and precipitation at a documented period/resolution, and compare them on the same grid. Select the mixing-depth and wind-flux mapping ranges. If ERA5 cannot distinguish them, retain the ranges as calibration uncertainties; do not present them as observed constants.
3. **G2 — Conservative offline prototype.** Implement only a diagnostic script or isolated prototype first. On analytic source/transport cases, test closed-sphere global closure, local availability, source attribution, orographic removal, zero sources, cap behavior, pole/seam handling, and climate-density convergence. Compare option A with option B for runtime and closure. If neither can close to `1e-8` relative to total source in double precision and remain positive, revise the discretization before production code.
4. **G3 — Climate impact and runtime.** After owner review of G2, run climate-only comparisons on the registered windows and held-out saves with both wind modes and climate densities 1 and 4. Report area-weighted ocean/land rain shares, latitude-band means, Köppen shares, precipitation rel-RMS to current source, surface-source amounts, wind transport, per-update closure, and median step-time cost. If held-out acceptance thresholds cannot be met without violating closure or exceeding the runtime budget, return to the RFC for a design choice before implementation.
5. **G4 — Downstream coupling.** Only after G3, rerun the eligible RFC 0016 coupling screen under the accepted atmospheric version. Report land rain, lake supply/volume, ocean export, frozen precipitation share, and ledger changes against its control. Recheck glacier accumulation, erosion, and biome distributions. This gate informs RFC 0016; it does not itself implement land-water storage or RFC 0015.

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/climate.py`: `compute_humidity`, `compute_precipitation`, `_land_moisture_source`, `_vegetation_transpiration_source`, `_orographic_retained_fraction`, and `compute_moisture_flux_convergence`; likely replace or retire `_humidity_zonal_sweep`, `_humidity_meridional_sweep`, `_retention_factor`, `_zonal_base_retention`, and `_meridional_base_retention`. Add a conservative column-water solver and named unit conversions.
- `backend/app/hydroclimate/atmosphere_cfd.py`: supply its current wind through the same `ClimateFields` interface; do not add moisture state in option A. Diagnostic wind remains through `climate.py`'s existing wind path.
- `backend/app/hydroclimate/erosion.py`, `hydrology.py`, `biomes.py`, RFC 0015 glacier code, and RFC 0016 land-water code: audit precipitation consumers and source ownership; avoid independent water writers in climate.
- Tests: add focused coverage in `backend/unit_tests/test_climate.py` for cell/global closure, available-water limits, source totals, cap conversion, source-free transport, both wind modes, pole/seam faces, resolution changes, and persistence compatibility. Add only a stress test if the production solve cost or coupled long-run behavior requires it.
- `docs/simulation-model/climate.md`, `docs/simulation-model/hydrology.md`, and `docs/architecture.md`: update the climate field contract, physical source units, diagnostic/prognostic boundary, and request flow after implementation.
- `bin/debug/measure_climate_fields.py`: reuse for the existing metrics; add an RFC-specific ledger/ERA5 diagnostic only if the current script cannot provide the required measurements.

## Assumptions and uncertainty

- A steady annual-mean moisture budget is an adequate representation for the simulator's geological timestep and current diagnostic precipitation consumers. It omits transient weather and atmospheric storage changes between climate calls.
- Column water multiplied by near-surface wind is not the true vertically integrated transport. The model needs a measured mapping or an explicit calibration range.
- A first-order rainout time is a bulk closure law, not cloud microphysics. The nine-day global scale does not set local rainout, convective efficiency, or storm tracks.
- Ocean evaporation, lake/river evaporation, and transpiration are surface fluxes with different drivers. Existing g/kg-equivalent source formulas do not constrain their physical rates without additional evidence.
- Land ET can cycle water repeatedly and therefore gross precipitation can exceed oceanic input locally and over a chosen interval. Conservation means each evaporation, atmospheric transport, precipitation, and surface storage transfer has one owner; it does not prohibit recycling.
- ERA5 is a modern-Earth comparison, while synthetic topography, wind, and land cover vary. Use it to constrain scales and patterns, not as a universal geological forcing.
- Changing rain affects erosion, lakes, Köppen classes, and glacier accumulation. Their movement is an impact to measure, not a reason to introduce a hidden rain floor.

## Acceptance criteria

- [ ] G0 reproduces its pinned evidence or the RFC is revised to explain why the premise changed.
- [ ] G1 documents verified units, source URLs/hashes, and the selected parameter ranges before model calibration.
- [ ] G2 demonstrates nonnegative water, rain no greater than local available water, no source-free land rain beyond transported ocean moisture, and per-cell/global water closure within `1e-8` relative to total source for analytic and climate-grid cases.
- [ ] With a zero land source, total land precipitation is no greater than ocean-origin moisture entering land, within ledger tolerance.
- [ ] A saturation cap and each precipitation enhancement convert their removed moisture to rain exactly once; there is no unaccounted humidity floor or ocean reset.
- [ ] Area-weighted closure holds at climate densities 1 and 4; global precipitation and land/ocean shares change by no more than 5% between them on the registered cases.
- [ ] On preregistered held-out saves, report land/ocean rain share, latitude bands, Köppen shares, and precipitation rel-RMS. Final thresholds are set by owner review of G0 before G3 and cannot be tuned on held-out saves.
- [ ] Both diagnostic and CFD wind complete the same budget; report step-time change and keep it within the owner-approved budget set before G3.
- [ ] RFC 0016's rerun reports the coupled and control changes with closed land and atmosphere ledgers; no treatment is used to tune G0–G3.
- [ ] Old saves load; option A adds no persistent state or save-size cost.

## Falsification criteria

- If the RFC 0016 rain multiplication does not reproduce on its pinned inputs, the specific coupling premise is false or stale; revise the RFC around the defect that does reproduce.
- If an area-conservative column solver cannot maintain positive water and close the finite-volume ledger at the specified tolerance on a closed sphere, reject the proposed discretization and try the streamline alternative before production work.
- If the only parameters that pass the ledger make precipitation patterns or downstream hydrology unusable on held-out worlds, the assumption that a steady diagnostic is sufficient is falsified. Reopen the choice between a more complete prognostic atmosphere and a simpler explicitly calibrated proxy.
- If the diagnostic and CFD winds imply irreconcilable transport mappings, do not share one fitted multiplier; revise the design to mode-specific transport or narrow the supported wind modes.
- If plausible rainout and transport sensitivities produce order-of-magnitude differences in land rain or lake supply, the current literature and diagnostics do not constrain this design enough for implementation; gather new evidence or withdraw the recommendation.

## References

- <a id="ref-trenberth1998"></a>
Trenberth, K. E. (1998). Atmospheric moisture residence times and cycling: Implications for rainfall rates and climate change. *Climatic Change*, 39, 667–694. [doi](https://doi.org/10.1023/A:1005319109110)

- <a id="ref-okikanae2006"></a>
Oki, T., & Kanae, S. (2006). Global hydrological cycles and world water resources. *Science*, 313(5790), 1068–1072. [doi](https://doi.org/10.1126/science.1128845)

- <a id="ref-vanderent2010"></a>
van der Ent, R. J., Savenije, H. H. G., Schaefli, B., & Steele-Dunne, S. C. (2010). Origin and fate of atmospheric moisture over continents. *Water Resources Research*, 46, W09525. [doi](https://doi.org/10.1029/2010WR009127)
