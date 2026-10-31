# RFC 0016: Effective land-water storage and upstream groundwater feedback

- **Status:** Draft
- **Date:** 2026-10-30
- **Related:** #718 (exploration); #680 (standalone follow-up on hold; formerly RFC 0014 G8); [PR #721](https://github.com/adubey/mantle-bloom/pull/721) (RFC 0014 withdrawal proposal); #675 (flicker characterization); #544 (RFC 0015 tracking); RFC 0015 (accepted); design principles 1, 3, 4, 5, 6, 8 and 9

## Problem

Mantle-Bloom has no persistent land-water reservoir. In `hydrology.compute_hydrology`, rain on previously wet lake area goes directly to lake supply. On dry area, the Fu–Budyko relation immediately partitions precipitation into routed runoff and implied evapotranspiration. There is no infiltration, soil-water inventory, groundwater storage, or return flow from a stored subsurface reservoir. The persisted `soil_depth` is geomorphic regolith thickness, not soil moisture or water capacity.

The climate module already recycles a vegetation moisture proxy into humidity, but estimates it as a biome-weighted share of previous local precipitation. It is not a transpiration flux withdrawn from land storage, and plants are not water-limited. Climate runs before hydrology, so a physically withdrawn transpiration flux can affect diagnostic humidity only on the next step unless the simulation order changes.

Rivers carry physical discharge in m³/s, but `discharge_m3_s` is a through-flow rate. The water currently held in each river cell is not represented in state or the lake ledger. `channel_depth` and `channel_width` describe persistent geomorphic incision and width; they do not describe the active wetted cross-section. `compute_river_speed` returns a stylized unitless value, not a physical velocity that can turn discharge into water volume.

This prevents the model from testing the selected land-water question: does rain first replenish a local subsurface store, with only water beyond storage capacity becoming surface flow, and does that upstream storage and plant use alter water delivered to lakes? The recent large-lake flickers do not prove that missing groundwater caused them. Issue #680 registered RFC 0014 G8 screens for local porewater and shoreline-aquifer exchange; it is now a standalone follow-up on hold after the owner's RFC 0014 withdrawal decision (closeout PR #721). This RFC asks a broader question about a per-node groundwater store and all upstream catchments, including its transpiration feedback. Reuse #680's registered checkpoints and any equivalent completed evidence, but do not repeat its distinct local screens.

For the intended land-water partition question, the current immediate runoff/ET split does not represent storage memory or water-limited transpiration (principle 1). A new store and river inventory must treat terrain nodes as areas and retain water across remapping (principles 4 and 6), use consistent physical units (principle 5), integrate fast hydrology over long tectonic steps (principle 3), and keep groundwater and river behavior tied to their own processes (principle 8). These are scope-specific omissions, not a claim that the model already violates whole-atmosphere conservation: diagnostic precipitation is intentionally an external forcing to the proposed terrestrial budget.

### Evidence

The source trace and decision history are recorded in [#718](https://github.com/adubey/mantle-bloom/issues/718). Current hydrology computes source volumes in m³/yr, routes them through the land graph, and reports discharge in m³/s. Lakes have a prognostic volume solve and ledger; rivers do not have a water-volume inventory. RFC 0015 owns physical snowfall, snowmelt, glacier and meltwater accounting. Until its physical-unit precipitation partition lands, the combined rain-and-snow terrestrial budget cannot close against current glacier rates.

The existing flicker characterization found whole-world lake-ledger closure and event-scale topology sensitivity, but did not identify a groundwater cause ([#675](https://github.com/adubey/mantle-bloom/issues/675)). Issue #680 records three exact A/B/C diagnosis windows and the owner's 10% event-magnitude / fixed-anchor connectivity screen; its G8 work is incomplete and on hold, and its local porewater and shoreline-exchange screens are not a test of all dryland stores upstream of a lake. RFC 0014's withdrawal is proposed in PR #721; this RFC reuses #680's window definitions and owner-selected diagnostic threshold without treating either as an accepted RFC 0014 criterion.

## Proposed solutions

Nature partitions precipitation among interception, evaporation, plant uptake, soil storage, surface runoff, recharge, and groundwater discharge. Storage and discharge depend on soil and aquifer properties, topography, river stage, and geology. Saturation excess can produce surface runoff, while groundwater contributions to rivers depend on hydraulic connection and head. Lakes may gain or lose groundwater depending on local gradients. A single effective store can represent water retention and a threshold overflow, but not these separate hydraulic states or reversible exchanges. Saturation-bucket behavior is used in conceptual catchment storage models, while explicit groundwater stores are added when water-table levels and flow direction matter [[Krasnostein & Oldham, 2004]](#ref-krasnostein2004). Storage-discharge laws are not universally linear, so the selected threshold is a testable abstraction rather than a general groundwater rule [[Botter et al., 2009]](#ref-botter2009).

**A. Keep immediate Fu–Budyko runoff and implied ET.** This is the least costly option and remains the control arm. It cannot answer whether stored water, saturation-priority partitioning, or water-limited plant use changes river/lake persistence.

**B. Add one effective store per terrain node, with saturation excess routed to surface water (recommended, gated).** On dry land, actual evaporation and water-limited transpiration are accounted first; the remaining rain and eligible melt replenish the node's store up to capacity, and only water beyond capacity follows existing river/lake routing. On the experimental arm only, transpiration is derived from available node storage and is passed one step later into diagnostic humidity. A paired gate measures the whole upstream effect on lake flicker before implementation is planned.

*For:* directly tests the user's stated storage priority, water-limited transpiration, and total catchment effect on lake inflow. It fits the node-area substrate and keeps water-table head out of the first model.
*Against:* store capacity and ET demand are poorly constrained on synthetic worlds; without head it cannot represent unsaturated baseflow or groundwater/lake exchange. Its total effect includes transpiration-driven humidity and any resulting diagnostic precipitation feedback.

**C. Add separate root-zone and aquifer stores with head-driven exchange.** This can distinguish plant-accessible water from slower groundwater and represent gaining/losing lakes and rivers.

*For:* supports more groundwater questions and physical exchange directions.
*Against:* Mantle-Bloom lacks porosity, hydraulic conductivity, water-table head, channel stage and aquifer geometry. It would add state and parameters the present model cannot constrain. Deferred unless the effective-store gates show that the missing distinction changes the target outcomes.

**Recommendation: B, only if the diagnosis gates pass.** First establish a river-water inventory estimate and closed whole-land ledger in the control. Then compare the control to the user's groundwater treatment on matched flicker events. A material result supports a scoped implementation design, but does not prove groundwater is the unique cause or establish coefficients. A nonmaterial result rejects this package as a primary explanation for these events, not the existence or broader importance of groundwater.

## Design

### Boundary and stores

Diagnostic climate precipitation is an imposed terrestrial input in each arm. Climate remains diagnostic: experimental transpiration is used as a humidity source one step later, and the resulting precipitation may differ between arms. The paired result therefore measures the total effect of the groundwater package, including its climate moisture-recycling path; it does not close atmospheric or ocean water.

Use RFC 0015's physical precipitation partition. Rain over existing lake area remains direct lake input. Frozen precipitation remains snow/ice under RFC 0015 until it melts; liquid bare-ground snowmelt may then enter the land-water partition. Rain and melt on ice follow RFC 0015's glacier-routing rules. The first groundwater treatment applies to liquid water on dry, non-ice land. It must not silently count snowfall as groundwater recharge or mix current scaled glacier rates with physical liquid-water volumes.

Represent the effective store as a per-node water-equivalent depth (g_i) in metres over node area (A_i), with stored volume (S_i = g_i A_i) in m³ and a nonnegative capacity C_i in m³. This store combines plant-accessible soil water and slower subsurface water as an effective quantity. It is not a water-table elevation or a literal aquifer. Do not infer C_i from `soil_depth` alone.

For node `i`, let `L_i` be liquid rain plus eligible bare-ground snowmelt over the interval (m³), `Eᵖ_i` be separately sourced potential land ET demand (m³), and `T_i` be the transpiration part of actual ET. The experimental arm must integrate, rather than endpoint-clamp, the interval balance in this order:
- compute actual soil evaporation and transpiration from demand and water available from the store plus current liquid input; debit these actual losses exactly once;
- fill the effective store from the remaining liquid input up to `C_i`; route only the remaining saturation excess to surface flow;
- saturation excess, with provenance recorded, enters the existing river/lake routing;
- transpiration is part of actual ET, is debited exactly once from the node's available land-water pool, and its volume is passed to the next step's diagnostic humidity calculation;
- no groundwater discharge occurs before saturation in the primary test. A source-supported leakage/baseflow law below capacity is an alternative for a later RFC revision, not silently added to this treatment.

This states the experiment's priority explicitly: after actual evaporation and transpiration, the rain/melt remainder becomes groundwater until the store saturates; only then does it become surface flow. Transpiration is limited by available water and comes from this combined store, while its humidity contribution is delayed one climate update. The gate must document the within-interval solver and conserve the sequence at all tested timestep partitions. The present `open_water_evaporation_mm` function is an open-water proxy, not a validated land PET law. Do not use it as land demand without a separate source and scale check. Do not retain the current precipitation-based vegetation proxy on top of store-derived transpiration in the treatment; the control retains current behavior.

The user-selected “groundwater first” rule is the treatment hypothesis, not an established universal hydrologic law. Per-node storage may fill unevenly and every terrain node represents finite area. The treatment does not infer subsurface lateral flow from a topographic gradient; once saturated, its water follows the existing surface-flow network to a downstream lake, river reach or ocean. Therefore the experiment tests a distributed saturation-storage/runoff pathway, not a head-resolved groundwater model.

### Per-cell river water estimate

Both control and treatment require a comparable in-channel inventory. For each routed reach, let `V_r,i` be persistent river water (m³), `I_i` its interval-mean inflow (m³/s), and `τ_i = L_i / U_i` its travel time, where `L_i` is the geodesic length of the flow edge (m) and `U_i` is a separately sourced physical reach velocity (m/s). A candidate linear-reservoir law is `dV_r,i/dt = I_i − V_r,i/τ_i`, with reach outflow `Ȯ_i = V_r,i/τ_i`. For constant interval-mean inflow over `Δt` seconds, integrate exactly: `V₁ = V₀ exp(−Δt/τ_i) + I_i τ_i (1 − exp(−Δt/τ_i))`; then book interval outflow as `O_i = V₀ + I_i Δt − V₁`. This makes each reach update conservative; route `O_i` to its downstream reach, lake or ocean. `Q_i L_i / U_i` is only the steady-flow inventory implied by this law, not a sufficient checkpoint-to-checkpoint storage history. Convert annualized model rates to m³/s before the update. Hydraulic geometry relates width, depth and velocity to discharge, but its coefficients vary with channel setting [[Leopold & Maddock, 1953]](#ref-leopoldmaddock1953) [[Singh, 2003]](#ref-singh2003). Manning-based depth estimates require channel width, slope, roughness and uniform-flow assumptions; they can be unreliable at confluences, constrictions and backwater controls [[Neal et al., 2021]](#ref-neal2021).

Do not use `compute_river_speed` as U_i: it is explicitly unitless. Do not treat `channel_depth` as water depth or `channel_width` as active wetted width without an independently justified mapping; they persist as geomorphic fields. Persist one extensive `V_r,i` for each eligible routed reach, including a reach with zero current inflow but remaining stored water; do not use `is_river`, which is a presentation threshold over high-discharge nodes, as the storage mask. Exclude lake and ocean storage already counted in their own reservoirs. Initialize old saves at the registered checkpoint with `V_r,i = Q_i τ_i` where the saved routed flow is positive, zero otherwise, and record this as the starting inventory (not an interval input); thereafter evolve only by the conservative reach law.

The river inventory remains an uncertain, source-supported hydraulic estimate, not a calibrated channel model. `V_r,i = Q_i τ_i` supplies a reproducible checkpoint initialization and steady-flow comparison; only the persistent state updated by the equation above enters the closed budget. G0 must choose and document the physical velocity/geometry method and uncertainty range, demonstrate equal-and-opposite reach transfers and interval closure in the control, and compare at least two plausible parameterizations. Verify the same transfer identity in both paired arms before interpreting G1. If an implementation cannot retain and update `V_r,i` conservatively, no river-inclusive closure or groundwater-effect conclusion is valid; revise the state/routing design before G1.

### Whole-land water ledger

Report all disjoint reservoirs as interval volumes in m³:
- effective groundwater storage `sum_i S_i`;
- persistent river water `sum_i V_r,i`;
- lake liquid volume and lake ice water equivalent;
- snow and glacier ice water equivalent under RFC 0015.

For each interval, close `ΔW_land = P_diag - ET_soil - T - E_river - E_lake - S_snow/ice - Q_ocean`.
All terms are integrated volumes over the interval. `P_diag` is the arm's own diagnostic rain plus snow input. `ET_soil` excludes separately reported transpiration `T`; `E_river`, `E_lake` and `S_snow/ice` are river evaporation, open-water evaporation and sublimation. `Q_ocean` includes river export and calving. Infiltration/recharge, saturation excess, groundwater-to-river transfer, river-to-lake transfer, lake spill within the land boundary, and snow/glacier melt are internal transfers and cancel in the total. If RFC 0015 has not landed, the rain-only diagnostic budget may be reported separately, but the full rain-plus-snow closure criterion is not met.

The river term must include both estimated storage change and routed loss/export. A change in `V_r` is not an unexplained residual or an additional input. Store-internal transfers must be equal and opposite, and no volume may be counted in groundwater, river and lake reservoirs at once. Record source categories (direct rain, bare-ground snowmelt, groundwater saturation overflow, glacier melt) through routing sufficiently to report how much of each target lake's supply came from the upstream groundwater treatment.

### Diagnosis gates

**G0 — conservative river inventory and control ledger.** Before the groundwater treatment:
1. Establish and document the river residence-time source and its uncertainty range. Choose at least two source-supported hydraulic parameterizations for `U_i`, calculate `τ_i`, and implement one persistent extensive river-water state per routed reach with the linear-reservoir update above. A derived `Q_i τ_i` snapshot without this state and update is not eligible for the strict water ledger.
2. At the registered diagnosis checkpoints, initialize control river state from saved `Q_i τ_i` and record the initial inventory. Run the current control hydrology through each interval, recording local water input, per-reach inflow, state change, outflow and routed destination. Every reach must satisfy `V₁ − V₀ = I_i Δt − O_i` within floating-point tolerance; internal transfers must cancel globally.
3. Extend the control ledger to include groundwater (zero), persistent river water, lake liquid/ice and RFC 0015 physical snow/glacier accounts. Reproduce the control on the fixed A/B/C windows with normalized interval residual ≤0.1% of total external water input. Report each river parameterization's effect on its inventory and the target-lake supply.
4. If the river law, source-supported range, initial-state rule or full control ledger fails, stop before any G1 treatment. Do not substitute a discharge snapshot or an unbooked inventory change.

Use the exact #680 A/B/C checkpoint windows below as G0/G1 defaults. Record the saved-world hashes and the pinned production source commit in every run report. A window may be excluded or replaced only for a preregistered input/data-integrity reason documented before running any treatment; identify the replacement checkpoint and hash before treatment. Do not use held-out saves to select the river or ET ranges. #680's distinct local porewater and shoreline-exchange screens remain on hold and must not be rerun. If an equivalent whole-upstream treatment is completed before this RFC's G1 starts, reuse its results rather than duplicate the arm.

**G1a — pre-register land ET and transpiration response.** Before any paired treatment run, select and document: (1) a source-supported potential land-evaporation and potential-transpiration method, with available climate/vegetation inputs and units; (2) how potential demand is split between soil evaporation and transpiration; (3) the water-limitation response using available `S_i` and current liquid input, including bounds and the rule that actual ET cannot exceed available water; and (4) source-backed ranges for demand, plant fraction and storage sensitivity, kept distinct from calibration assumptions. The response must debit each actual flux once, treat transpiration as the store-derived portion, and pass its volume to diagnostic humidity one climate update later. Freeze the equations, parameter ranges, seeds/checkpoints and readouts before G1b. If no defensible demand/split law can be supported with available inputs, G1 is blocked and must be revised before paired comparison.

**G1b — total upstream groundwater effect.** After G0 and G1a pass, compare same-checkpoint paired arms on the fixed A/B/C diagnosis windows, disjoint from all held-out acceptance saves:
- **Window A:** seed `910211954`, 221.6–222.0 Myr, target production-schedule ΔV `+495,624 km³`, checkpoint SHA-256 `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340`.
- **Window B:** seed `910211954`, 222.1–222.5 Myr, target production-schedule ΔV `+131,625 km³`, checkpoint SHA-256 `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544`.
- **Window C:** seed `317739385`, 201.0–201.4 Myr, target production-schedule ΔV `−261,592 km³`, checkpoint SHA-256 `041b7430a6cf8867f2e2c76f18c3328e440add585ee4a2e4b3055e7766c2217e`.
- The registered checkpoint runs use production source commit `f4438e1c85042de789bd464ce5d2f58c8c0717eb` (#680). Verify this source and all three saved worlds are available before preregistration; if not, document the missing input and select/hash replacements before any paired arm.
- **Control:** current hydrology and climate behavior, plus the same persistent river state and ledger as G0.
- **Treatment:** per-node effective storage; the frozen G1a land-ET/transpiration law; storage filled from remaining liquid water before surface saturation excess; store-derived transpiration passed one step later to diagnostic humidity; all other model rules unchanged. Each arm uses the same step schedule and common random inputs where possible. Record per-arm precipitation differences caused by the humidity feedback.
- Report target component (|ΔV|), flooded area, fixed-anchor wet/dry/connectivity sequence, event timing, per-cell river inventory, groundwater inventory, ET/T, river/lake/ocean transfers and ledger residual. Attribute whether upstream store inputs reach the target lake through saturation overflow or change humidity/precipitation.
- On the primary production schedule, report the event magnitude comparison. Separately report the registered 10 kyr connectivity trace; do not treat it as a replacement for the production-schedule comparison or as timestep invariance.
- Use #680's owner-selected interpretation: a material response is ≥10% change in event (|ΔV|), or a fixed-anchor connectivity/event-class mismatch lasting two consecutive 10 kyr outputs, in at least two of the three windows. This is a diagnostic mechanism-screen threshold, not an accepted RFC 0014 criterion; RFC 0014 is being withdrawn in PR #721.

A material G1 result shows that the combined catchment store and plant-water feedback can affect the tested events. It does not separate direct saturation overflow from transpiration-driven climate feedback; report both pathways and propose an ablation only if the total effect is material. A below-threshold result across the source-supported range deprioritizes this package as the primary cause of those events, without rejecting groundwater's broader role. Missing source-supported ranges or missing river hydraulics are blocked results, not negative evidence.

### Timestep, bounds, units and knobs

Use physical units throughout: water depth in m, cell area in m², inventory in m³, discharge in m³/s, annual forcing in m³/yr, and interval duration in years. Clamp fractions to [0, 1], storage to [0, C_i], and all flux volumes to nonnegative values. Bounds must not delete excess water; capacity overflow is explicitly routed.

The tectonic interval is much longer than river residence time and may exceed soil-storage and ET response times. A one-update-per-step bucket could miss fill, saturation and depletion transitions. Integrate the finite-capacity store over the path of its input and outflows, or use a closed-form update for the selected flux law. Compare equivalent integrated input under multiple step partitions before implementation; RFC 0014's withdrawn G7 characterization and #687 show that event connectivity can depend on schedule.

No global capacity, hydraulic roughness, active channel geometry, or land-PET coefficient is selected here. Storage observations for the upper metre/root zone show substantial spatial and definitional variation [[Gupta et al., 2023]](#ref-gupta2023); they are not estimates of combined soil-plus-groundwater capacity. Global recharge/baseflow estimates show that groundwater can matter to rivers, with broad model spread, but are not transferable per-node rates [[Xie et al., 2024]](#ref-xie2024). Derive the treatment range from sources appropriate to each quantity, keep capacity and ET demand separate, and sensitivity-test all design choices on diagnosis cases only. Never tune against held-out flicker events. #675's hypothetical capacity scenarios are not parameter bounds.

### Save compatibility and expected cost

A future persisted groundwater field must be a water-equivalent extensive surface field whose area-integrated volume is conserved through cell refinement, coarsening, splitting, merging and plate reassignment. Old saves must load with zero groundwater storage and a documented initialization rule. River volume is persistent state, not a derived ledger term; remap it conservatively and initialize it from `Q_i τ_i` at the preregistered old-save checkpoints before starting each closed continuation ledger.

The bucket, reach reservoir and ledger add `O(N)` work, one persistent groundwater float per terrain node and one river-water float per routed reach (plus any explicitly justified provenance state). Gate measurements must report median step-time impact at density 4, with a target no more than +10% over the paired control. Reuse the existing flow graph and routing pass where possible; do not add a separate whole-world graph traversal without measuring its cost.

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/hydrology.py`: `_partition_precipitation_sources`, `_route_external_water_sources`, `compute_hydrology`, river discharge and flow-edge geometry. Add persistent reach water and its conservative update, the effective-store update, and the source-specific water ledger after the RFC's diagnosis/design gates pass.
- `backend/app/hydroclimate/climate.py`: `_vegetation_transpiration_source`, `compute_humidity`, `compute_climate`. Replace the current rainfall proxy in the treatment with previous-step store-derived transpiration; preserve the control behavior.
- `backend/app/hydroclimate/lakes.py`: `solve_lakes_prognostic` and its water ledger, to receive groundwater-derived routed supply without double booking.
- `backend/app/hydroclimate/ice_sheet.py` and RFC 0015 D1/D4/D8 interfaces: consume its physical rain/snow/melt partition; do not repeat RFC 0015 glacier-conversion work.
- `backend/app/surface_fields.py`: register groundwater and river-water inventories as persistent extensive fields with conservative remapping semantics.
- `backend/app/plates.py`, `backend/app/hydroclimate/erosion.py`, `backend/app/world.py`: gather/commit the field and keep within-step climate/hydrology state ordering explicit.
- Tests: `backend/unit_tests/test_hydrology.py`, `test_climate.py`, `test_persistence.py`, and surface-field/remapping tests. Add unit checks for volume conversion, nonnegative capped storage, one-step transpiration lag, source attribution, remap conservation and full ledger closure.
- Docs: `docs/architecture.md`, `docs/simulation-model/hydrology.md`, `docs/simulation-model/climate.md`, and `docs/save-compatibility.md`.
- New diagnostic tooling, if needed: `bin/debug/measure_land_water.py` (new), with a usage docstring; compact findings in `analysis/rfc0016/README.md`. Keep saves and full replay outputs outside the repository.

## Assumptions and uncertainty

- The target is land-water partition, persistence, plant water limitation, and the total upstream effect on the registered lake-flicker events. Diagnostic precipitation is not an atmospheric water store.
- One effective store per node is an intentionally unresolved mixture of root-zone soil water and deeper groundwater. It cannot predict water-table elevation, lateral aquifer gradients, losing lakes, or plant access to a separate root zone.
- The initial treatment uses saturation excess as the groundwater-to-surface return. It does not include unsaturated baseflow; adding it requires a source-supported relation and an RFC revision or accepted design choice.
- Potential land ET, storage capacity and timing are not calibrated. `open_water_evaporation_mm` is not a land transpiration law. Earth global statistics are sanity checks only; lithology, soil, relief and climate differ on synthetic worlds.
- River water storage is particularly uncertain because the model lacks active wetted geometry and physical reach velocity. A diagnostic estimate is useful only if its uncertainty is explicit and its changes reconcile with the ledger.
- Since transpiration changes diagnostic humidity and precipitation one step later, a paired treatment may differ in precipitation. The gate estimates the total feedback effect; it does not close the atmosphere or identify every mediator without further ablations.
- RFC 0015 is accepted but its physical snow/ice ledger is not yet implemented. Full precipitation-to-land closure, including snow and glaciers, depends on RFC 0015. Do not mix the current scaled glacier boundary with physical rain.
- G1 reuses #680's registered A/B/C events and owner-selected diagnostic materiality rule; #680's local mechanism screens are distinct. Its whole-upstream treatment is not evidence that direct shoreline exchange caused flicker.

## Acceptance criteria

Diagnosis gates come first.

- [ ] **G0 ledger and river inventory:** establish the conservative reach-state law and source range, then reproduce the control on fixed windows A/B/C with all river, groundwater, lake and snow/ice inventories and transfers reconciled; normalized interval residual is ≤0.1% of total external water input. If not, stop; no groundwater-effect conclusion.
- [ ] **G1a ET/transpiration pre-registration:** before any paired run, freeze a source-supported potential-demand method, its soil-evaporation/transpiration split, water-limitation equations, input units, sensitivity ranges and selected checkpoints. If this response cannot be specified, G1 remains blocked.
- [ ] **G1b upstream groundwater screen:** compare control and treatment on fixed windows A/B/C using their exact hashes and pinned source commit, disjoint from held-out saves. Report the owner-selected 10% (|ΔV|) and two-consecutive-output 10 kyr connectivity indicators, event timing, per-cell river/GW storage, store-derived transpiration, diagnostic precipitation response, source transfers and ledger residual.
- [ ] **Precipitation partition:** with RFC 0015 physical units, rain/snow/melt/ice partitions sum to the per-cell precipitation volume; no water is lost outside a named reservoir, transfer or sink.
- [ ] **Land-ledger closure:** effective storage + river inventory + lake liquid/ice + snow/glacier water equivalent close to ≤0.1% of total external input per interval. Internal transfers are equal and opposite and each water volume has one owner.
- [ ] **Store bounds and remapping:** groundwater stays in [0, C_i]; saturation excess is routed; area-integrated storage is conserved through every surface remap and plate reassignment to within floating-point tolerance.
- [ ] **ET and climate coupling:** actual ET never exceeds available water; transpiration is included once in the land sink and reaches diagnostic humidity exactly one step later. Control and treatment source definitions are recorded.
- [ ] **Timestep sensitivity:** equivalent integrated forcing at the tested partitions gives bounded storage, ET, runoff and lake-response differences with an owner-approved tolerance set before the acceptance run. Report connectivity separately; do not infer convergence from a single schedule.
- [ ] **Held-out validation:** after gates and parameter choices are frozen, evaluate on at least three held-out mature saves disjoint from all diagnosis and tuning. Report lake storage/flicker, land storage, water partition, budget residual and uncertainty range.
- [ ] **Runtime:** median full-step time at density 4 is within +10% of matched control; report storage and ledger memory/save-size impact.
- [ ] **Compatibility:** old saves load with zero groundwater and river state initialization consistent with the RFC; new saves round-trip the water inventory.

## Falsification criteria

- If no source-supported river geometry or physical velocity range can be established, the proposed river inventory and whole-land ledger cannot support the gate. Revise the river representation or withdraw this RFC; do not use the unitless speed proxy to force a result.
- If the ledger cannot close after adding persistent river storage, or its residual changes materially across equivalent timestep partitions, the current reach-storage law is not a defensible reservoir. Stop before interpreting flicker.
- If every source-supported treatment setting stays below both registered materiality outcomes in all diagnosis windows, conclude that this saturation-store plus transpiration package is not a primary cause of the tested flicker. Do not claim groundwater is unimportant to other land-water outcomes.
- If the treatment materially changes flicker but the change comes only from an out-of-range or unsupported capacity, ET or hydraulic assumption, the causal screen is invalid; revise its ranges and preregister before rerunning.
- If the whole-upstream arm duplicates an experiment completed by #680 or its follow-up, reuse that evidence and remove the duplicate gate. A repeated or overlapping run is not new support.
- If the one-store abstraction cannot meet both mass closure and the selected persistence/plant-water outcomes without separately resolving root-zone and aquifer stores, replace it with a multi-store RFC revision; do not label one scalar as a water table.

## References

- <a id="ref-krasnostein2004"></a>
Krasnostein, A. L., & Oldham, C. E. (2004). Predicting wetland water storage. *Water Resources Research*, 40, W10203. [doi](https://doi.org/10.1029/2003WR002899)

- <a id="ref-leopoldmaddock1953"></a>
Leopold, L. B., & Maddock, T. (1953). The hydraulic geometry of stream channels and some physiographic implications. *U.S. Geological Survey Professional Paper 252*. [doi](https://doi.org/10.3133/pp252)

- <a id="ref-singh2003"></a>
Singh, V. P. (2003). Downstream hydraulic geometry relations: 1. Theoretical development. *Water Resources Research*, 39(6), 1146. [doi](https://doi.org/10.1029/2003WR002484)

- <a id="ref-neal2021"></a>
Neal, J., Hawker, L., Savage, J., Durand, M., Bates, P., & Sampson, C. (2021). Estimating river channel bathymetry in large scale flood inundation models. *Water Resources Research*, 57, e2020WR028301. [doi](https://doi.org/10.1029/2020WR028301)

- <a id="ref-gupta2023"></a>
Gupta, S., Lehmann, P., Bickel, S., Bonetti, S., & Or, D. (2023). Global mapping of potential and climatic plant-available soil water. *Journal of Advances in Modeling Earth Systems*, 15, e2022MS003277. [doi](https://doi.org/10.1029/2022MS003277)

- <a id="ref-xie2024"></a>
Xie, J., Liu, X., Jasechko, S., et al. (2024). Majority of global river flow sustained by groundwater. *Nature Geoscience*, 17, 770–777. [doi](https://doi.org/10.1038/s41561-024-01483-5)

- <a id="ref-botter2009"></a>
Botter, G., Porporato, A., Rodriguez-Iturbe, I., & Rinaldo, A. (2009). Nonlinear storage-discharge relations and catchment streamflow regimes. *Water Resources Research*, 45, W10427. [doi](https://doi.org/10.1029/2008WR007658)
