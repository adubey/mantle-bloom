# RFC 0016: Effective land-water storage and upstream groundwater feedback

- **Status:** Draft
- **Date:** 2026-10-30
- **Related:** #718 (exploration); #680 (RFC 0014 G8 lake-flicker screens); #538 (RFC 0014 tracking); #675 (flicker characterization); #544 (RFC 0015 tracking); RFC [0014](0014-lake-water-volume-and-freezing.md); RFC [0015](0015-frozen-water-in-physical-units.md); design principles 1, 3, 4, 5, 6, 8 and 9

## Problem

Mantle-Bloom has no persistent land-water reservoir. In `hydrology.compute_hydrology`, rain on previously wet lake area goes directly to lake supply. On dry area, the Fu–Budyko relation immediately partitions precipitation into routed runoff and implied evapotranspiration. There is no infiltration, soil-water inventory, groundwater storage, or return flow from a stored subsurface reservoir. The persisted `soil_depth` is geomorphic regolith thickness, not soil moisture or water capacity.

The climate module already recycles a vegetation moisture proxy into humidity, but estimates it as a biome-weighted share of previous local precipitation. It is not a transpiration flux withdrawn from land storage, and plants are not water-limited. Climate runs before hydrology, so a physically withdrawn transpiration flux can affect diagnostic humidity only on the next step unless the simulation order changes.

Rivers carry physical discharge in m³/s, but `discharge_m3_s` is a through-flow rate. The water currently held in each river cell is not represented in state or the lake ledger. `channel_depth` and `channel_width` describe persistent geomorphic incision and width; they do not describe the active wetted cross-section. `compute_river_speed` returns a stylized unitless value, not a physical velocity that can turn discharge into water volume.

This prevents the model from testing the selected land-water question: does rain first replenish a local subsurface store, with only water beyond storage capacity becoming surface flow, and does that upstream storage and plant use alter water delivered to lakes? The recent large-lake flickers do not prove that missing groundwater caused them. RFC 0014 G8 (#680) is screening local porewater and shoreline-aquifer exchange, among other event mechanisms. This RFC asks a broader question about a per-node groundwater store and all upstream catchments, including its transpiration feedback; it must reuse rather than repeat any G8 result that already tests the same whole-catchment process.

For the intended land-water partition question, the current immediate runoff/ET split does not represent storage memory or water-limited transpiration (principle 1). A new store and river inventory must treat terrain nodes as areas and retain water across remapping (principles 4 and 6), use consistent physical units (principle 5), integrate fast hydrology over long tectonic steps (principle 3), and keep groundwater and river behavior tied to their own processes (principle 8). These are scope-specific omissions, not a claim that the model already violates whole-atmosphere conservation: diagnostic precipitation is intentionally an external forcing to the proposed terrestrial budget.

### Evidence

The source trace and decision history are recorded in [#718](https://github.com/adubey/mantle-bloom/issues/718). Current hydrology computes source volumes in m³/yr, routes them through the land graph, and reports discharge in m³/s. Lakes have a prognostic volume solve and ledger; rivers do not have a water-volume inventory. RFC 0015 owns physical snowfall, snowmelt, glacier and meltwater accounting. Until its physical-unit precipitation partition lands, the combined rain-and-snow terrestrial budget cannot close against current glacier rates.

The existing flicker characterization found whole-world lake-ledger closure and event-scale topology sensitivity, but did not identify a groundwater cause ([#675](https://github.com/adubey/mantle-bloom/issues/675)). RFC 0014 G8 has registered event windows and a 10% event-magnitude / fixed-anchor connectivity materiality interpretation; its paired mechanism work is still in progress ([#680](https://github.com/adubey/mantle-bloom/issues/680)). Its local porewater and shoreline-exchange screens are not a test of all dryland stores upstream of a lake.

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

Both control and treatment need a comparable in-channel inventory estimate. For each active routed reach, a candidate volume is `V_r,i = Q_i × L_i / U_i`, where `Q_i` is physical discharge through the node (m³/s), `L_i` is the geodesic length of its flow edge (m), and `U_i` is a separately estimated physical reach velocity (m/s). This is discharge times reach travel time; it returns m³. Estimate U_i from a source-supported hydraulic-geometry relation or a Manning-style cross-section using discharge, slope, wetted width and roughness. Hydraulic geometry relates width, depth and velocity to discharge, but its coefficients vary with channel setting [[Leopold & Maddock, 1953]](#ref-leopoldmaddock1953) [[Singh, 2003]](#ref-singh2003). Manning-based depth estimates require channel width, slope, roughness and uniform-flow assumptions; they can be unreliable at confluences, constrictions and backwater controls [[Neal et al., 2021]](#ref-neal2021).

Do not use `compute_river_speed` as U_i: it is explicitly unitless. Do not treat `channel_depth` as water depth or `channel_width` as active wetted width without an independently justified mapping; they persist as geomorphic fields. Include every active flow reach with `Q_i > 0` and a valid downstream edge, not only `is_river`, which is a presentation threshold over high-discharge nodes. Exclude lake and ocean storage already counted in their own reservoirs.

The river estimate is initially a diagnostic inventory, not a new calibrated hydraulic model. The G0 gate must choose and document the physical velocity/geometry estimate and its uncertainty range on both arms. It must also demonstrate how changes in this inventory enter the water ledger. If a diagnostic V_r cannot be made compatible with the reach transfers and interval storage change, the gate is blocked until the RFC is revised to add an explicit conservative river-water state and routing law.

### Whole-land water ledger

Report all disjoint reservoirs as interval volumes in m³:
- effective groundwater storage `sum_i S_i`;
- estimated river water `sum_i V_r,i`;
- lake liquid volume and lake ice water equivalent;
- snow and glacier ice water equivalent under RFC 0015.

For each interval, close `ΔW_land = P_diag - ET_soil - T - E_river - E_lake - S_snow/ice - Q_ocean`.
All terms are integrated volumes over the interval. `P_diag` is the arm's own diagnostic rain plus snow input. `ET_soil` excludes separately reported transpiration `T`; `E_river`, `E_lake` and `S_snow/ice` are river evaporation, open-water evaporation and sublimation. `Q_ocean` includes river export and calving. Infiltration/recharge, saturation excess, groundwater-to-river transfer, river-to-lake transfer, lake spill within the land boundary, and snow/glacier melt are internal transfers and cancel in the total. If RFC 0015 has not landed, the rain-only diagnostic budget may be reported separately, but the full rain-plus-snow closure criterion is not met.

The river term must include both estimated storage change and routed loss/export. A change in `V_r` is not an unexplained residual or an additional input. Store-internal transfers must be equal and opposite, and no volume may be counted in groundwater, river and lake reservoirs at once. Record source categories (direct rain, bare-ground snowmelt, groundwater saturation overflow, glacier melt) through routing sufficiently to report how much of each target lake's supply came from the upstream groundwater treatment.

### Diagnosis gates

**G0 — control-arm river inventory and whole-land ledger.** Before the treatment is interpreted:
1. Apply the per-cell river inventory estimator to the control at the registered diagnosis checkpoints. Compare at least two plausible, independently sourced hydraulic parameterizations and report their volume range and its share of total land water and event-lake supply.
2. Extend the existing lake ledger to include groundwater (zero in control), estimated river inventory, and RFC 0015 physical ice/snow accounts. Report every external source, sink, inventory change and internal transfer.
3. Reproduce the control on at least three registered event windows, with interval residual ≤0.1% of total external water input. Verify node-area weighting and that river storage changes reconcile to reach inflow/outflow.
4. If the river-volume estimate changes flicker or lake supply materially across its defensible range, or the ledger cannot close, stop and revise the estimator/ledger before the groundwater arm.

Reuse #680's registered diagnosis windows and event baselines where suitable; do not use its held-out saves for range selection or tuning. #680's current work screens local lake porewater and shoreline exchange. Do not rerun those tests. If #680 later completes a whole-upstream groundwater treatment under equivalent controls, reuse its results rather than duplicate the arm.

**G1 — total upstream groundwater effect.** After G0 passes and the #680 overlap check is complete, compare same-checkpoint paired arms on at least three diagnosis flicker windows, disjoint from all held-out acceptance saves:
- **Control:** current hydrology and climate behavior, plus the same river inventory estimator and ledger as G0.
- **Treatment:** per-node effective storage; water-limited actual ET/transpiration; storage filled before surface saturation excess; store-derived transpiration passed one step later to diagnostic humidity; all other model rules unchanged. Each arm uses the same step schedule and common random inputs where possible. Record per-arm precipitation differences caused by the humidity feedback.
- Report target component (|ΔV|), flooded area, fixed-anchor wet/dry/connectivity sequence, event timing, per-cell river inventory, groundwater inventory, ET/T, river/lake/ocean transfers and ledger residual. Attribute whether upstream store inputs reach the target lake through saturation overflow or change humidity/precipitation.
- On the primary production schedule, report the event magnitude comparison. Separately report the registered 10 kyr connectivity trace; do not treat it as a replacement for the production-schedule comparison or as timestep invariance.
- Reuse #680's owner-selected interpretation: a material response is ≥10% change in event (|ΔV|), or a fixed-anchor connectivity/event-class mismatch lasting two consecutive 10 kyr outputs, in at least two of the three windows. This is a mechanism-screen threshold, not a new flicker acceptance criterion. Keep RFC 0014's existing acceptance threshold unchanged.

A material G1 result shows that the combined catchment store and plant-water feedback can affect the tested events. It does not separate direct saturation overflow from transpiration-driven climate feedback; report both pathways and propose an ablation only if the total effect is material. A below-threshold result across the source-supported range deprioritizes this package as the primary cause of those events, without rejecting groundwater's broader role. Missing source-supported ranges or missing river hydraulics are blocked results, not negative evidence.

### Timestep, bounds, units and knobs

Use physical units throughout: water depth in m, cell area in m², inventory in m³, discharge in m³/s, annual forcing in m³/yr, and interval duration in years. Clamp fractions to [0, 1], storage to [0, C_i], and all flux volumes to nonnegative values. Bounds must not delete excess water; capacity overflow is explicitly routed.

The tectonic interval is much longer than river residence time and may exceed soil-storage and ET response times. A one-update-per-step bucket could miss fill, saturation and depletion transitions. Integrate the finite-capacity store over the path of its input and outflows, or use a closed-form update for the selected flux law. Compare equivalent integrated input under multiple step partitions before implementation; RFC 0014 G7 already shows that event connectivity can depend on schedule.

No global capacity, hydraulic roughness, active channel geometry, or land-PET coefficient is selected here. Storage observations for the upper metre/root zone show substantial spatial and definitional variation [[Gupta et al., 2023]](#ref-gupta2023); they are not estimates of combined soil-plus-groundwater capacity. Global recharge/baseflow estimates show that groundwater can matter to rivers, with broad model spread, but are not transferable per-node rates [[Xie et al., 2024]](#ref-xie2024). Derive the treatment range from sources appropriate to each quantity, keep capacity and ET demand separate, and sensitivity-test all design choices on diagnosis cases only. Never tune against held-out flicker events. #675's hypothetical capacity scenarios are not parameter bounds.

### Save compatibility and expected cost

A future persisted groundwater field must be a water-equivalent extensive surface field whose area-integrated volume is conserved through cell refinement, coarsening, splitting, merging and plate reassignment. Old saves must load with zero groundwater storage and a documented initialization rule. The river inventory can remain derived if G0 proves that its endpoint estimate and changes can be booked consistently; otherwise it needs explicit persistent state and the same remap/save compatibility treatment.

The bucket and ledger add `O(N)` work and one persistent float per node for storage (plus any explicitly justified provenance state). Gate measurements must report median step-time impact at density 4, with a target no more than +10% over the paired control. Reuse the existing flow graph and routing pass where possible; do not add a separate whole-world graph traversal without measuring its cost.

## Detailed design (provisional touchpoints)

- `backend/app/hydroclimate/hydrology.py`: `_partition_precipitation_sources`, `_route_external_water_sources`, `compute_hydrology`, river discharge and flow-edge geometry. Add the effective-store update and source-specific water ledger only after gates pass.
- `backend/app/hydroclimate/climate.py`: `_vegetation_transpiration_source`, `compute_humidity`, `compute_climate`. Replace the current rainfall proxy in the treatment with previous-step store-derived transpiration; preserve the control behavior.
- `backend/app/hydroclimate/lakes.py`: `solve_lakes_prognostic` and its water ledger, to receive groundwater-derived routed supply without double booking.
- `backend/app/hydroclimate/ice_sheet.py` and RFC 0015 D1/D4/D8 interfaces: consume its physical rain/snow/melt partition; do not repeat RFC 0015 glacier-conversion work.
- `backend/app/surface_fields.py`: register the new persistent water-equivalent field with conservative extensive remapping semantics.
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
- G1 reuses RFC 0014 G8 event-level attribution/materiality results where possible. Its whole-upstream treatment is not evidence that direct shoreline exchange caused flicker.

## Acceptance criteria

Diagnosis gates come first.

- [ ] **G0 ledger and river inventory:** on at least three registered diagnosis windows, the river estimate is physically dimensioned, its source range is documented, all water inventories and transfers reconcile, and normalized interval residual is ≤0.1% of total external water input. If not, stop and revise; no groundwater-effect conclusion.
- [ ] **G1 upstream groundwater screen:** compare control and treatment under the same checkpoints and schedules on at least three diagnosis windows, disjoint from held-out saves. Report the registered 10% (|ΔV|) and two-consecutive-output 10 kyr connectivity materiality results, event timing, per-cell river/GW storage, store-derived transpiration, diagnostic precipitation response, source transfers and ledger residual. Do not change RFC 0014's flicker threshold.
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
- If the ledger cannot close after adding estimated river inventory, or its residual changes materially across equivalent timestep partitions, the current diagnostic river estimator is not a defensible reservoir. Stop before interpreting flicker.
- If every source-supported treatment setting stays below both registered materiality outcomes in all diagnosis windows, conclude that this saturation-store plus transpiration package is not a primary cause of the tested flicker. Do not claim groundwater is unimportant to other land-water outcomes.
- If the treatment materially changes flicker but the change comes only from an out-of-range or unsupported capacity, ET or hydraulic assumption, the causal screen is invalid; revise its ranges and preregister before rerunning.
- If the whole-upstream arm duplicates an experiment completed by RFC 0014 G8 or its follow-up, reuse that evidence and remove the duplicate gate. A repeated or overlapping run is not new support.
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
Gupta, H. V., et al. (2023). Global mapping of potential and climatic plant-available soil water. *Journal of Advances in Modeling Earth Systems*, 15, e2022MS003277. [doi](https://doi.org/10.1029/2022MS003277)

- <a id="ref-xie2024"></a>
Xie, J., et al. (2024). Majority of global river flow sustained by groundwater. *Nature Geoscience*, 17, 501–506. [doi](https://doi.org/10.1038/s41561-024-01483-5)

- <a id="ref-botter2009"></a>
Botter, G., Porporato, A., Rodriguez-Iturbe, I., & Rinaldo, A. (2009). Nonlinear storage-discharge relations and catchment streamflow regimes. *Water Resources Research*, 45, W10427. [doi](https://doi.org/10.1029/2008WR007658)
