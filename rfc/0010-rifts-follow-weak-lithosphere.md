# RFC 0010: Rifts follow weak lithosphere

- **Status:** Accepted (revised 2026-10-30 after fixed-geometry Gate 3b result)
- **Date:** 2026-10-15
- **Last revised:** 2026-10-30; Gate 3b requires an evolving 3-D slab body coupled to deforming plate geometry before slab-pull rifting can count
- **Related:** #699 (Gate 3b tracking), #487 (closed Gate 3a evidence), #481 (superseded tracker), #480 (min-cut diagnostic), #472 / PR #479 (RFC 0009 G3), #389 (RFC 0003 tracking), #378, #385, #388 (RFC 0003 gates), PR #418 (RFC
  0003 solution 1); RFCs [0003](abandoned/0003-rifting-around-cratons.md) (superseded by RFC 0010),
  [0009](0009-craton-mosaic.md) (cratons-first layout; G3 initially selected the min-cut),
  [0007](0007-mantle-plumes-and-plateaus.md) (plume weakening),
  [0012](0012-physics-first-simulation.md) (force-response interfaces)

## Revision note

- **2026-10-30 ([PR #717](https://github.com/adubey/mantle-bloom/pull/717)):** The prior Gate 3b source criterion looked for a broad,
  sustained load, then coupled it to local crack rules. The corrected 1 Gyr diagnostic created broad
  load and damage but no viable daughter plates; it explicitly held plate geometry fixed and used an
  effective slab-sink length instead of an evolving slab body ([Gate 3b result](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6538397727)). This
  showed that fixed-geometry load screens cannot answer whether slab pull bends a slab and stretches
  a plate into a rift. The revision gates slab-pull attribution on a reduced-order 3-D slab body,
  endogenous plate finite strain and geometry, and a controlled pair of outward slab pulls before
  geological-history acceptance runs.

  | Changed item | Previous criterion | Revised criterion | Reason |
  |---|---|---|---|
  | Gate 3b event/source | Broad sustained extension supplied by a low-cost current field, then local crack response | For slab-pull cases, require evolving attached 3-D slab bodies, measured bending, and force-transmitted finite strain in the plate | Fixed-geometry runs can ramp load and accumulate damage but cannot demonstrate that the slab body or plate geometry evolves |
  | Controlled case | At least one broad event on a generated plate | Opposed outward slab tractions on one plate, plus no-slab and detached-slab controls, before the event/source criterion | Separates slab-causal plate extension from an imposed trench-force pattern |
  | Plate-scale histories | Run 200 Myr histories once the load source and crack response are coupled | Run histories only after the controlled slab/plate coupling demonstrates deformation, local crack linkage, and closed force/material accounts | Avoids spending long histories on fixed geometry that cannot produce physical breakup |
  | State and runtime | Persist damage and segment history with O(N+E) expected work | Also persist slab-body and plate-deformation state, with O(N+E+S) expected work for S bounded slab stations | The added mechanics changes geometry and has continuing state across saves and steps |

## Problem

A splitting plate breaks along a line the mantle flow alone decides, so continents break
wherever the flow puts that line, through cratons and strong lithosphere as readily as
through weak belts.

**How a plate splits today** (`merge_split.maybe_split_plate`, origin/main `1227c36`):
1. A plate at least 2 Myr old (`SPLIT_MIN_AGE_YEARS`) whose flow one rigid rotation fits
   poorly has its node velocities clustered into two regimes (`kmeans2`).
2. `_refine_split_cut` places a great circle midway between the two clusters' centroids
   and iterates until the halves stop changing. The great circle is there for a
   computational reason: the raw cluster labels aren't spatially contiguous, and a plane
   through the centre always gives two clean halves.
3. If the halves' poles are far enough apart, the plate either splits along that circle
   (`RIFT_SUCCESS_PROBABILITY` = 0.55) or fails to rift, and `apply_failed_rift` thins a
   band `FAILED_RIFT_BAND_MULT` = 2.5 spacings either side of it.

The cut never looks at the lithosphere. On origin/main `02d35b9`, 25–36% of successful
splits cross craton on seeds 7, 11 and 2024 (PR #418); rifting takes 8.0–9.6% of the
starting craton per Gyr, and stretch donation and divergent thinning along the new
boundaries are most of that loss (#378, #388).

The current split path also has `Plate.internal_stress`, a single dimensionless, per-plate
accumulator based on plate radius and overlap. It relaxes the flow-fit and pole-separation
gates; it is not a local stress field in Pa and does not locate or grow a rift. Likewise,
`torque.BoundaryForceInputs` contains node-wise geometry and material values, but
`torque.predict_plate_motion` integrates their force contributions into a whole-plate torque
and one angular acceleration. Neither is, by itself, a local extensional load suitable for
damage accumulation. The revised design must build or explicitly reject a defensible local
loading field before choosing a constitutive law; it must not relabel either existing signal
as local stress.

**RFC 0003 tried to steer the cut and failed.**
- **Tilting the great circle** (PR #418): of 158 flow cuts that crossed craton, a
  craton-free candidate passing the admissibility gates existed for 19 (12%) within 20°,
  and at most 30 (19%) anywhere in the sweep ranges. A great circle can only move as a
  whole, so missing a craton sweeps a large strip of plate to the other side.
- **Letting the cut bend** (the RFC 0010 strength-graph probe): a minimum cut through a strength
  graph follows the flow as well as the great circle does, but still crosses craton in
  most cases. Of 135 craton-crossing cuts on the same seeds, 13 (10%) get a craton-free
  admissible minimum cut within a 1,000 km corridor, and at most 16 (12%) at any corridor
  up to 2,500 km. Wider corridors find craton-free routes only by slicing a coastal sliver
  off the plate, which the size gate rejects ([analysis/rfc0010](../analysis/rfc0010/README.md)).

The bent cut initially failed for a different reason than the tilted one: craton fills each
continent's interior, leaving no seam to follow. RFC 0009 addressed that layout problem by
placing weak belts around cratonic cores. #480 then tested the strength-graph cut on the new
layout; the static separator still missed its routing target. The remaining design question is
how local damage accumulates and propagates, not how to find a cheaper whole-plate line.

RFC 0009's G3 initially selected a min-cut prototype on pooled L/B thresholds, but this was a
branch-selection gate, not proof that it met the per-seed acceptance target. In #480's 1 Gyr
same-state diagnostic, cut-link crossings fell by only 11.8%, 20.0% and 28.6% for seeds 7, 11
and 2024. The detailed admitted-only diagnostic gave B = 0.690, 0.747 and 0.643 (pooled 0.703),
with 30 of 78 admitted routes line-clear but only 3 with no craton in the wider band. The margin
instrumentation rerun has different event totals from the earlier production-probe aggregate;
those denominators are not combined ([analysis/issue480](../analysis/issue480/README.md)). These
same-state diagnostics do not test every possible edge cost or replace a full matched flow-only
replay. They do show that the current static-cut strategy misses its line-routing target and does
not represent local damage accumulation.

## Proposed solutions

**What happens in nature.** Continental rifts exploit inherited weaknesses, including ancient
orogenic belts and lithospheric fabrics, though the structures that reactivate depend on their
orientation and strength [[Wilson, 1966]](#ref-wilson1966)
[[Vauchez et al., 1997]](#ref-vauchez1997) [[Tommasi & Vauchez, 2001]](#ref-tommasi2001)
[[Audet & Bürgmann, 2011]](#ref-audet2011) [[Heron et al., 2016]](#ref-heron2016). Models
of the East African Rift show how plume–lithosphere interaction can branch rifts around a cratonic
keel [[Koptev et al., 2015]](#ref-koptev2015); cratonic survival and rejuvenation depend on
coupled mantle/lithosphere behavior and heating or metasomatism
[[Lenardic et al., 2003]](#ref-lenardic2003) [[Foley, 2008]](#ref-foley2008). Continental rifts also
depend on interacting forces and evolving lithospheric strength. Rift systems commonly begin
as separate fault segments that propagate, interact and link; progressive thinning changes the
force balance, and some rifts arrest before breakup [[Brune et al., 2023]](#ref-brune2023)
[[Kolawole & Ajala, 2024]](#ref-kolawole2024). The Afar linkage observations show that transfer
between rift tips can involve distributed extension and rift-perpendicular shear; a transform-like
link is not necessarily a single strike-slip fault that slices a plate [[Pagli et al., 2019]](#ref-pagli2019).
Plate motion and rifting are affected by slab, ridge, mantle and gravitational-potential-energy
loads, whose relative importance varies by setting [[Forsyth & Uyeda, 1975]](#ref-forsyth1975)
[[Brune et al., 2023]](#ref-brune2023).

1. **Minimum cut through a strength graph (rejected after #480).** It bends a whole-plate
   separator toward weak links, but minimizes a static path cost; it does not accumulate strain,
   model local failure or propagate rift tips. The 1 Gyr same-state diagnostics reduced cut-link
   craton crossings by only 11.8%, 20.0% and 28.6% for seeds 7, 11 and 2024. In the detailed
   admitted-only diagnostic, B was 0.690, 0.747 and 0.643 (pooled 0.703): a 29.7% reduction in
   craton area in the measured band, not 50%. Only 30 of 78 admitted routes were line-clear and
   only 3 had a craton-free band. These results reject min-cut as the model strategy for #480;
   changing its edge weights would continue to tune a geometric proxy rather than the failure
   process ([analysis/issue480](../analysis/issue480/README.md) and
   [gate-margin summary](../analysis/issue480/mincut-gate-margins-1gyr.json)). The detailed
   margin replay's event totals differ from the earlier production-probe table, so its counts are
   kept as a separate diagnostic and must not be pooled with that table.
2. **Nucleate and propagate (an incomplete local-growth proposal).** A rift begins where extension
   overcomes local strength and advances through neighboring lithosphere; it can stall rather than
   invoke a success coin flip. This proposal does not define persistent damage, a segment-linkage
   rule or the local loading field, so it is not yet an implementable model.
3. **Stress-field trigger (incomplete on its own).** Derive intraplate loading from the plate
   force balance and compare tensile force with lithospheric strength. This supplies a physical
   driver, but an instantaneous threshold still cannot represent accumulated damage, intermittent
   activity or arrested segments.
4. **Raise craton resistance (RFC 0003 solution 3).** This can reduce craton loss, but cannot
   choose a geologically plausible rift path or express segment growth. Keep it as a separate
   material-retention control, not as a substitute for a rift model.
5. **Plume weakening (RFC 0007).** Heating and melt can weaken lithosphere and may influence
   initiation or propagation. Keep it coupled through RFC 0007; a plume is not required for every
   rift.
6. **Local cuts that may link into a rift (new; recommended).** Evaluate short cuts between
   adjacent surface cells under local extensional loading. Repeated loading can damage and
   weaken one segment, advance its tips into neighboring links, or let it arrest. Multiple
   separated cuts can link through distributed extension and transfer/shear deformation. No
   step searches for or applies a complete plate-wide separator. The plate splits only after
   the accumulated local cuts form a connected, viable boundary with two material-retaining
   sides and compatible divergence; until then the plate continues evolving as one plate.
   Transform motion transfers displacement between offset segments; it is not an independent
   whole-plate cut.

**Recommendation (revised after #480 evidence, 2026-10-18).** Treat solution 1 as failed for this
problem. Replace it with solution 6. RFC 0012 supplies a plate-scale torque model. A new
per-node source-force interface and reduced thin-shell diagnostic now transfer those loads into
local membrane stress resultants, following the broad thin-shell precedent in Kong & Bird (1995).
The diagnostic now compares the largest positive principal *increment* of its in-plane membrane
resultant with the column's integrated additional extensional strength from an assumed `K0`
background. The pairing is a normal-faulting criterion only while vertical compression remains
the major principal stress; it does not represent intact-rock tensile rupture or a prestressed,
yielding shell. Extensional lithosphere commonly fails by normal-fault frictional sliding at
shallow depth and ductile creep at depth
[[Brace & Kohlstedt, 1980]](#ref-brace1980) [[Kohlstedt et al., 1995]](#ref-kohlstedt1995).
The standalone column-envelope diagnostic uses published wet-quartzite, plagioclase,
mafic-granulite and wet/dry olivine flow-law parameterizations, with an explicit frictional
branch [[Byerlee, 1978]](#ref-byerlee1978) [[Chen & Gerya, 2016]](#ref-chen2016)
[[Wang et al., 2012]](#ref-wang2012)
[[Hirth & Kohlstedt, 2003]](#ref-hirth2003). Its scenario ranges are hypotheses to test, not a spatially
calibrated Earth strength field; the model lacks the local temperature, composition, fluid and
stress state needed to choose among them.
An independent curved-sphere manufactured solution now validates that stress invariant with
decreasing error through 17 × 33 nodes. Free in-plane displacement is selected as the boundary
condition for the isolated single-plate diagnostic because its inputs provide no displacement
anchor; clamps are sensitivity probes. A shared-node two-plate manufactured benchmark now checks
numerical interface transfer for uniform isotropic tension, but it does not define a geologic
contact law. The matched diagnostic replay validates equilibrium eligibility and compares physical and explicitly
phenomenological capacities on Gate 1 candidates whose source cells form one exact edge-connected
component. Gate 2 is explicitly scoped to this connected-source population: the 56 genuinely
disconnected-source events are counted and reported separately, but excluded from the Gate 2 coverage
denominator because this isolated shell solve cannot represent their disconnected material domains.
This scope does not imply that connected events represent all events, and their capacity scores must
not be generalized to disconnected-source events. A conditional scale audit now brackets a
100 km mobile/rift-like column at 6.28 TN/m and a 150 km thick-craton scenario at 20.86 TN/m,
with wide constitutive sensitivity ranges; only the 100 km value has a regional strength-scale
comparison. These are reference columns used to set scale, not a selected or calibrated spatial
field. Gate 2 is closed as a validation diagnostic, not a production rupture-threshold decision:
the connected-source shell solve meets the stated coverage, equilibrium and boundary-condition checks,
and the capacity alternatives are reported with their assumptions and sensitivity. The replay
has no independent rupture labels and the local material state remains unresolved, so its capacity
scores are conditional scenarios; they cannot calibrate a production damage threshold. Physical
capacity selection and coupled neighbor/interface mechanics remain outside Gate 2 and require separate
evidence before any production damage law uses those quantities.
The existing flow fit and plate-level `internal_stress` may remain candidate-event or
background-eligibility signals, but they do not choose a complete cut or stand in for local loading. Keep solution 4 separately sized against the
craton-loss gate, and take plume weakening through RFC 0007. The connected-boundary rule and
diagnosis gates are proposed here; constitutive values remain unsettled, so no implementation issue
should start until those gates are resolved and the revised design is accepted.

## Design

### Diagnosis gates

1. **Reconcile the diagnostic event stream before calibrating — closed 2026-10-19 by retiring
   the archives from calibration.** Two provenance-recorded paired 1 Gyr runs on the same model
   sources and runtime reproduced every production cut row for seeds 7, 11 and 2024; tap totals
   matched the production ledger (144, 111 and 125). Use those paired runs as the repeatable
   baseline. The older production and margin archives have unrecoverable runtime/worktree
   provenance and diverge from that baseline within their first 3–6 events; explicitly exclude
   them from calibration and keep their denominators separate. The archived discrepancy is not
   explained or erased, but it no longer blocks use of the documented baseline. See
   [`analysis/issue487/README.md`](../analysis/issue487/README.md).
2. **Gate 2 — validate the conditional local-load/capacity diagnostic — closed 2026-10-22.
   Incremental extensional measure, column scenario diagnostic, curved-shell benchmark, curved-mesh
   checks and isolated-plate edge policy
   added; connected-source mesh repairs and a 2026-10-21 skip audit completed. Scope is craton-crossing
   candidates whose source cells form one exact edge-connected component. The 56 disconnected-source
   candidates are counted separately and excluded from the coverage denominator. The connected-source
   coverage requirement passes. Under the scope decision below, Gate 2 validates a conditional
   single-plate diagnostic; it does not select a production capacity field or validate coupled stress.
   A thickness-aware spatial capacity screen now uses per-node Hc/Hm,
   crust type, craton share, belt marker and a bracketed Moho-lag profile on the pinned Gate 1
   replay; the capacity field remains scenario-dependent and underidentified.**
   `torque.gather_per_node_force_field` preserves each force
   source as a per-node vector in N. `thin_shell.solve_force_field` projects the combined loads
   onto the shell tangent plane, removes the area-weighted rigid-rotation component handled by
   plate motion, and solves small-strain triangular membrane equilibrium. It returns local
   incremental depth-integrated stress resultants and compares their maximum positive principal
   increment with an explicit additional normal-fault capacity, all in N/m. This is a diagnostic
   normal-faulting measure: extension failure is usually governed by brittle
   frictional sliding at shallow depth and temperature/rate-dependent creep at depth. A candidate
   column diagnostic calculates effective lithostatic pressure, a Mohr-Coulomb normal-fault differential
   strength and a layer-specific ductile flow stress, then integrates the weaker branch through
   depth. The result is an extensional differential-strength resultant in N/m and is paired with the
   shell increment using an assumed background horizontal stress. The candidate must bracket pore pressure, geotherm,
   rock rheology, strain rate and background stress because the current state does not resolve them;
   details and caveats are recorded in
   [`analysis/issue487/gate2-load-diagnosis.md`](../analysis/issue487/gate2-load-diagnosis.md).
   A scale-calibration audit evaluates the full 1,728-case constitutive grid at 70, 90, 100, 135
   and 150 km total thickness, holding the existing 35 km continental crust fixed. The reference
   uses `K0=0.8`, freshwater hydrostatic pressure `P_f=rho_water*g*z`, friction coefficient
   `mu=0.6` (`phi=atan(mu)`), and the model's 1330 C base temperature. The 100 km column gives
   6.28 TN/m; the factorial median is 4.88 TN/m (p10-p90: 1.49-12.07). This falls within the
   3-8 TN/m integrated strength Stamps et al. report for a regional 100 km East African Rift
   lithosphere, but definitions and spatial averaging differ [[Stamps et al., 2014]](#ref-stamps2014).
   This is an order-of-magnitude check, not a fitted target. Hydrostatic pressure is computed from
   water density and depth; the former fixed ratio `P_f/P_lith=0.42` is now clearly treated as an
   overpressure comparison. Friction is bracketed by `mu=0.6-0.85`, based on Byerlee's rock-friction
   experiments [[Byerlee, 1978]](#ref-byerlee1978). K0, cohesion, temperature, strain rate and
   water/rheology state remain assumptions because the world has no corresponding local fields.
   The K0 audit finds that the force interface supplies external loads and the shell solve supplies
   stress increments, while neither contains tectonic prestress; literature stress observations
   constrain regime/orientation better than a globally transferable magnitude ratio. Thus `K0=0.8`
   remains a reproduction anchor only, not a selected physical value. The independent assumption-basis
   audit records that EARS stress observations support the normal-faulting regime but do not identify
   the model's spatial `K0`; the geodetic strain-rate
   upper scale is about `6.3e-16 s-1`, so `1e-14 s-1` is an extreme localization sensitivity;
   and the wet-olivine law is an endmember with a stated water content, not a mapped mantle field.
   The existing 1250/1330/1400 C thermal boundaries are checked against their implied heat flow
   and remain conductive scenarios, not local geotherm constraints. Hydrostatic pressure is
   supported only as a crustal benchmark, not as a full mantle-lid fluid model. See
   [`analysis/issue487/gate2-assumption-basis.md`](../analysis/issue487/gate2-assumption-basis.md)
   and its generated JSON. No legacy event outcome is used to set these assumptions.
   A fixed-load crossed screen then scored all 22,896 connected-source node loads at
   `6.3e-16 s-1`, crossing conditional `K0≈0.32/0.6/0.8/1.0`, dry/wet olivine and full-column versus
   capped crust-hydrostatic pressure. The added `K0≈0.32` is the cohesionless near-critical limit
   for the same assumed `mu=0.6` friction law, not an independent estimate. The extended finite
   screen changes 57.59% of node classifications, with 3.07% over capacity in every case and 60.66%
   in at least one. The original central cases show lower exceedance (3.07–19.03%); the capped
   pressure profile changes fewer than 1% of classifications at matched K0/rheology there. No
   scenario is assigned a probability. The screen's robust/possible bounds and scenario details are in
   [`analysis/issue487/gate2-assumption-screen.json`](../analysis/issue487/gate2-assumption-screen.json)
   and are not probabilities. These checks show that the result is highly sensitive to an
   unobserved background stress state; Gate 2 still lacks an independently selected production
   capacity field. The current world state and force interface do not contain tectonic prestress, and
   stress literature does not identify one globally transferable `K0`, so `0.8` remains a replay
   anchor only. The K0 audit and expanded result are documented in
   [`analysis/issue487/gate2-assumption-basis.md`](../analysis/issue487/gate2-assumption-basis.md).
   The pressure-profile, Mohr-Coulomb, and dry/wet creep-inversion formulas pass six analytical unit checks. The shell-load
   transfer has separate manufactured-solution checks documented in
   [`analysis/issue487/gate2-load-diagnosis.md`](../analysis/issue487/gate2-load-diagnosis.md).
   A 150 km Tanzania/Congo craton thickness is used in one regional geodynamic model, motivating a
   distinct strong-column scenario, not a strength observation [[Rajaonarison et al., 2021]](#ref-rajaonarison2021).
   A separate rifting model's 22 TN/m isotropic peak is a driving force, not capacity, and is not
   fitted [[Duretz et al., 2025]](#ref-duretz2025). Full results, assumptions and reproducibility are
   in [`analysis/issue487/gate2-capacity-calibration.md`](../analysis/issue487/gate2-capacity-calibration.md).
   The old-algorithm replay candidates are not rupture/no-rupture labels, so their load distribution
   is not used to fit strength. Replay nodes carry local `Hc`, `Hm`, crust type, craton share, belt
   history and Moho thermal lag, which now feed the matched spatial screen. The world still lacks
   depth-resolved temperature, composition, pore pressure, water, strain-rate and fabric state.
   This is conditional scale calibration, not production calibration. On the same fixed 22,896 node
   loads, the hydrostatic reference gives 1,563 nodes (6.83%) above capacity, compared with 1,718
   (7.50%) using fixed `P_f/P_lith=0.42`. Those are scenario screens, not rupture labels.
   Across the finite one-factor envelope, 795 nodes (3.47%) exceed capacity in every tested variant
   and 2,701 (11.80%) in at least one; node rank correlation remains at least 0.9936. These are the
   earlier one-factor results; the expanded 16-case screen above includes the low-K0 limit. The `K0`
   range 0.6-1.0 moves the fraction from 11.23% to 3.47%; wet olivine yields 9.74%. These bounds
   omit interactions and do not identify a production threshold. Full results and field limitations
   are in `analysis/issue487/gate2-capacity-sensitivity.md`,
   `analysis/issue487/gate2-capacity-sensitivity.json` and
   `analysis/issue487/gate2-capacity-state-audit.md`.
   The elastic thickness and
   modulus, boundary conditions, default convex-hull mesh and rigid-mode removal are documented
   diagnostic assumptions, not calibrated production behavior. Explicit mesh connectivity now
   preserves concave outlines. A manufactured uniaxial load on a
   shallow cap recovers its analytic resultant within 0.16% across 3 × 3 to 9 × 9 node grids, with
   less than 0.004% spread in the recovered mean. An irregular convex cap also recovers the
   analytic free-edge load within 0.12%; clamping its west edge shifts the mean resultant by 2.82%
   with a small equilibrium residual. An L-shaped concave mesh recovers a manufactured uniaxial
   free-edge load within 0.16% at 9 × 9 and 17 × 17 resolution, with less than 0.001% spread. On
   that finer mesh, west and south clamps shift the mean resultant by 16.5%; these are sensitivity
   checks and do not establish a preferred physical edge condition. A
   deterministic generated sample (seed 7, 12 plates, node density 0.5) selected its largest
   hole-bearing oceanic plate and largest single-loop continental plate. Their point-set hulls
   included 2.48% and 2.38% outside-territory triangle centroids, respectively. Explicit
   territory-clipped meshes retained two and one boundary loops, and recovered a controlled
   shallow-cap manufactured resultant within 0.09% for both. This tests actual generated node
   distributions and topology after mapping to a shallow cap. The same two meshes were also solved
   on their original broad curved coordinates under projected manufactured edge tractions; residuals
   were below `8.5e-5` of load. Four localized clamp placements shifted mean resultants by up to
   23.5% on the oceanic plate and 8.9% on the continental plate. These curved cases show numerical
   convergence and boundary sensitivity, not production load distributions. A separate
   uniform-tension spherical-shell benchmark decreases mean error from 2.94% to 1.18%, with
   fine-grid p95 absolute error 2.75%. Free displacement is selected for the isolated diagnostic;
   paired clamps shift means by up to 46.6% and full-boundary clamps suppress the mean by about
   100%, demonstrating why arbitrary anchors are not adopted.
   A two-plate shared-node benchmark passes manufactured force/torque balance, interface reaction and
   resultant checks. A corrected six-event sample (one split and one failed-rift per seed) measures
   5.3–97.0% mean absolute source-band N_ext change on a band-scale denominator when the dominant
   nearest neighbor is included. For four events with nonzero mutual contact drag, pair force closure
   residuals are 2.0–21.4% and torque closure residuals are 4.6–21.6%; two events have no mutual
   contact drag and do not test closure. The initial 26.2%/51.2% per-node-relative summaries were
   invalid where isolated stress was zero and are superseded. The tied/no-slip variant omits only
   mutual pair drag, but remains an idealized bonded interface, not a physical finite-slip law.
   See `analysis/issue487/gate2-interface-sensitivity.md`. These are sensitivity checks, not
   physical coupled predictions. The pair-level projection audit finds 3.0–95.1% cut-band stress sensitivity between joint and separate-group rigid-mode removal. The joint projection matches the bonded mesh; the separate-group variant does not represent free slip because it retains bonded displacements. A mechanically consistent interface/load-transfer formulation remains unresolved.

   **Gate 2 scope decision (2026-10-22).** Gate 2 evaluates the single connected source plate using its per-node force field, free in-plane edge condition, incremental shell resultant and column capacity scenario. Coupled neighbor displacement and a finite-slip interface law are not prerequisites for closing this diagnostic gate. The six-event pair replay remains sensitivity evidence only: its projection choices represent different kinematics, and current pair contact forces do not consistently close force and moment. These pair results must not be presented as physical coupled stress predictions. A future coupled-stress claim requires a separately specified and validated interface model. This clarifies the earlier wording that required resolving the interface before treating coupled values as physical; that requirement applies only to future coupled predictions, not Gate 2’s single-plate diagnostic.
   A direct phenomenological capacity candidate,
   `T_cap = T_ref * (1 - lambda * w)`, was swept on a prescribed weak band: increasing `T_ref`
   monotonically reduced failed nodes, and at `2e10 N/m` the weak-band failure fraction exceeded
   the background fraction for both plates. This verifies the candidate rule's response, not its
   physical parameters or geological realism; it remains an unaccepted design alternative. These
   remain synthetic diagnostics. `orogeny.steady_moho_temperature_c` and `moho_temperature_c` provide a coarse
   conductive Moho-temperature proxy, while `rheology.yield_excess` applies a boundary-local
   scalar Mohr-Coulomb proxy. Neither supplies a validated temperature-dependent vertical yield
   envelope or integrated extensional differential capacity. The current maximum-positive-tension
   shell comparison uses a positive principal stress increment paired with additional normal-fault
   capacity only under the stated `K0` and normal-faulting assumptions.
   The first full replay scored 130 identical Gate 1 craton-crossing candidates (39,467
   cut-band node samples) and reported 29 shell solves above the 1% equilibrium-residual screen.
   All nine inspected high-residual cases had disconnected triangle components that touched only
   at vertices. The subsequent exact-source audit classified 56/130 candidates as genuinely
   disconnected under `plate.adjacency()`; these remain excluded as complete events because
   retaining one source island or solving islands independently would invent unsupported load paths
   or boundary conditions. For the other 74 source-connected candidates, the replay now repairs all
   ten previously edge-disconnected meshes: six by clipping a spherical hull with centroid and
   edge-midpoint probes and four by retaining a dominant component containing at least 99% of
   triangles. Repaired meshes must cover at least 95% of both candidate cut-band and craton-band
   nodes. Every accepted solve meets stop code 1 and residual/load ≤1%; the maximum is 0.4514%.
   Production cut counts and ordered candidate keys match Gate 1 exactly for all seeds. The old
   29/130 and 64/130 metrics remain historical pre-repair results, not current equilibrium
   estimates.

   **Coverage population and threshold.** Gate 2 covers craton-crossing candidates whose source
   plate is one connected component under exact cell-edge adjacency. Require valid shell solves for
   at least 90% of these eligible candidates overall and 80% in each seed, event kind (`split` /
   `failed_rift`), world-time bin (0–200, 200–500, 500–1000 Myr), and nonempty event-kind × time
   stratum. A valid solve must also meet the per-mesh 95% cut-band and craton-band sample coverage,
   stop code 1 and residual/load ≤1%. These are project completeness thresholds, not physical
   constants. Empty strata are not applicable. The 74/74 connected-source candidates pass: overall
   and every nonempty stratum have 100% solver coverage. The other 56/130 Gate 1 candidates (43.1%)
   have genuinely disconnected source domains; report them separately by seed, kind and time, but
   exclude them from this Gate 2 denominator because this isolated shell diagnostic has no defensible
   load path for those components. Their frequency is uneven (source-connected coverage of the full
   population was 47.8% for splits versus 67.2% for failed rifts, and 68.0%, 41.2%, and 42.9% across
   successive time bins), so connected-source results are conditional and must not be generalized to
   disconnected-source events or the full population. On 22,896 connected-source cut-band samples,
   the physical median-capacity failure fractions are 0.274% continental, 78.03% young-ridge oceanic
   proxy, and 2.84% mature oceanic; at the phenomenological 1× reference, 2.82% of samples fail.
   These remain screening results: the oceanic-age proxy is coarse, the capacity anchors are
   conditional rather than spatially calibrated, and natural interplate continuity and production
   load distributions remain outside this single-plate solve. Gate 2 closure is scoped to this
       single-plate diagnostic; it does not validate coupled neighbor stresses or imply a finite-slip
       contact law. Do not calibrate the damage law or
   infer local stress from applied traction, rigid rotation, flow clusters or `Plate.internal_stress`.
   See the audits in `analysis/issue487/README.md`, `analysis/issue487/gate2-load-diagnosis.md`,
   `analysis/issue487/gate2-capacity-calibration.md`, and `analysis/issue487/gate2-matched-replay.json`.
3. **Gate 3a — local crack mechanics; Gate 3b — plate-scale rift-potential event and source
   (scope revised 2026-10-30; tracking issue #699).** Gate 3a's local crack, propagation, arrest
   and linkage rules are the response mechanism for this work. Gate 3b now works backward from at
   least one candidate event: first define and detect a broad, sustained extensional loading event
   on a real spherical plate; characterize its footprint, direction, magnitude relative to the
   existing local crack threshold, duration/accumulation, force and work accounts, and timestep
   sensitivity; then find the simplest low-cost mantle-current representation that can supply the
   measured load pattern; finally couple that source to the Gate 3a local response and check that
   cracks grow within the loaded footprint before a viable boundary can split the plate. A
   tens-of-millions-of-years persistence window (initially 10–30 Myr) is a hypothesis to examine,
   not a fixed timer or accepted cutoff. A slab-pull explanation is eligible only when subducted
   lithosphere has an evolving three-dimensional body and its force is transmitted into a plate
   whose in-plane geometry and strain can evolve. The current boundary-derived slab-pull term and
   fixed-geometry shell screens do not meet that condition. Freeze the operational event definition,
   source amplitude bracket, footprint criterion and controlled cases in #699 before new probes.
   The event label is diagnostic bookkeeping, not a new force or damage law. Existing prescribed-axis
   routes and static coarse shell maps are exploratory screens; they neither identify a mantle/magma
   source nor establish event persistence. If the current source representation cannot produce an
   eligible event, revise the source hypothesis before tuning local crack thresholds or starting
   broad geological-history acceptance runs.

### Mechanics (provisional)

- **State and loading.** After diagnosis gate 2, store cumulative irreversible extensional strain
  (`epsilon_rift`, dimensionless) on the state element whose area and ownership are defined by the
  load discretization. Each local cut is an adjacency link between cells; a failed cell or link
  becomes a short rift segment, not a plate partition. The thin-shell diagnostic returns tensile
  membrane stress increment `N_ext = max(0, lambda_max(Delta_N))` and compares it with the column's
  additional normal-fault strength resultant `N_cap`, both in N/m. The diagnostic pairing is
  defined. Gate 2 validates only this conditional diagnostic; production damage use waits for a
  separately justified capacity selection and the later damage-calibration gate. If the
  damage law remains
  in Pa, convert both through the same effective thickness, `sigma_ext = N_ext / H_eff` and
  `Y_t = N_cap / H_eff`; otherwise formulate an integrated-viscosity law with units. The current
  torque balance drives whole-plate motion; its loads do not by themselves define intraplate stress.
- **Damage update.** Where the extensional increment exceeds local yield, integrate the inelastic extension
  rate over the portion of the step in which it is positive:
  `d epsilon_rift / dt = max(0, sigma_t - Y_t) / eta_rift`, with `eta_rift` in Pa·s. Cap
  `epsilon_rift` at `epsilon_fail` and convert it to a bounded damage fraction
  `D = clip(epsilon_rift / epsilon_fail, 0, 1)`. A provisional bounded softening form is
  `Y_t(D) = Y_t(0) * [1 - D * (1 - q_res)]`, where the residual-strength fraction `q_res` is in
  [0, 1]. Reaching `D = 1` marks local failure; this creates/extends a segment but does not split the
  plate unless the failed links form a viable connected boundary. `epsilon_fail`, `eta_rift`,
  `q_res`, any healing time, and whether this softening form is adequate are calibration assumptions
  with no approved numeric values yet. They must be set after diagnosis gate 3, with physical scale
  checks and sensitivity runs; no parameter should be selected to hit the craton-crossing percentage
  alone.
- **Local cuts, nucleation and growth.** Mark an adjacency link failed only when its neighboring
  state has reached the calibrated damage/failure condition and its displacement can be
  accommodated without creating an unowned gap. Nucleate a finite segment when accumulated strain
  and local loading satisfy the calibrated failure rule. Extend a segment tip only into adjacent
  cells where the local tensile load can sustain damage; permit multiple segments and abandonment.
  Link overlapping or interacting segments through distributed extension and transfer/shear zones.
  The graph identifies
  local adjacency and connectedness only; it does not search globally for the least-cost separator.
  Strike-slip/transform response transfers relative displacement between segments and may change tip
  loading, but cannot split the plate by itself.
- **Rift completion.** Create a new plate only when linked local cuts form a connected, viable boundary
  with two material-retaining sides and the required kinematic divergence. Until then, keep one plate
  and apply incremental thinning/opening through existing deformation and ledger paths. Failed or
  abandoned rifts retain their modeled damage/weakness unless a measured healing law removes it.
  Do not apply the current whole-band one-shot thinning plus success coin flip to the same event.
- **Crust type.** Permit a propagating spreading boundary to enter oceanic lithosphere when its local
  loading and damage criteria allow it. Track continental, oceanic and mixed-crust segment length
  separately so continental breakup is not inferred from a coastline crossing.
- **Slab body and deforming plate geometry (required before slab-pull attribution).** The current
  model removes subducted surface cells and derives slab pull from the trench neighborhood; it has
  no persistent body below the surface and does not accumulate in-plane deformation through the
  plate interior. Those terms can move a rigid plate, but cannot demonstrate that a bent slab pulls
  a plate apart into a rift. The candidate low-cost representation is a three-dimensional slab
  ribbon attached to the subducting plate: a small set of material stations along each active trench
  segment and down dip, carrying position, width, thickness, density/temperature anomaly and
  attachment identity. Subducted surface-column material transfers into this body once, so the slab
  is not counted both as live plate material and as a separate pull source. Its path and curvature
  evolve as material enters the trench and sinks. Pull arises from
  integrated slab buoyancy, with bending resistance and mantle drag opposing motion; slab length,
  thickness, dip, curvature and consumed material are measured state, not an effective sink-depth
  parameter. A sparse two-dimensional plate shell embedded in three-dimensional space must accumulate
  in-plane displacement and strain, update cell geometry and thickness, and pass the resulting
  extension into the Gate 3a crack/damage law. Damage can feed back on shell strength, and local
  geometry must change before a new plate is split. Slab attachment, equal-and-opposite force,
  moment and work accounting, slab consumption, and cell-area/material ledgers remain explicit
  through subduction, rifting, split/merge and save/load. Where the body exits the modeled mantle,
  transfer remaining mass and material to a named deep reservoir and account for the lost buoyancy
  and work; never silently drop the body. Slab-to-plate attachment forces are an
  internal action/reaction pair; slab buoyancy and mantle drag are external force/work terms and
  must be budgeted separately. This reduced-order coupling is the target;
  full mantle convection and a full 3-D finite-element solution are not prerequisites. Subduction
  force balances include slab pull, interplate resistance, bending and basal traction
  [[Li & Gurnis, 2023]](#ref-li2023), while 3-D rift models show distributed strain and fault-network
  evolution can control localization and linkage [[Naliboff et al., 2020]](#ref-naliboff2020).
- **Gravity and mantle coupling.** Use loads already represented in the plate torque model—slab pull,
  ridge push and basal drag—as contributors only where a defensible spatial load distribution exists.
  Include lithospheric gravitational-potential-energy gradients when column thickness/topography can
  support a dimensionally consistent force. Mantle-flow vectors may assist or resist motion through
  traction; they are not a prescribed rift line. When an explicit slab body drives a plate, retire or
  subtract the existing boundary-proxy slab-pull contribution for that same slab. Do not count the
  same load twice between rigid plate torque, slab-body force, shell equilibrium and membrane stress.
- **Rift-potential source diagnosis (Gate 3b, revised 2026-10-30).** Work from an observed
  candidate broad extension event back to the smallest current field that can supply its measured
  load, then run that source forward through the existing Gate 3a local crack response. The current
  `mantle.flow_at` field (summed `ConvectionCenter` velocities) and
  `torque.basal_drag_force_vectors_n` (relative-flow basal drag at owned nodes) are candidate
  starting interfaces, not evidence that magma currents or rift-driving regional traction are
  represented. Audit their spatial scales, time behavior and force ownership first. A source model
  may use a small number of broad current lobes or another comparably simple field, but its
  amplitudes and event footprint must be fixed by the event/source diagnostic rather than tuned to
  force a desired crack route. Preserve equal-and-opposite traction and global force/torque/work
  accounts, and ensure each basal load is counted once across rigid plate motion and membrane
  equilibrium. Test persistence over a preregistered duration bracket including 10–30 Myr as a
  hypothesis; the event state accumulates actual load or extension and never acts as an independent
  timer.

**Units, bounds and timestep.** Stress and strength use Pa; strain is dimensionless; viscosity uses
Pa·s; damage stays in [0, 1]; thickness stays nonnegative. Integrate `max(0, sigma_t - Y_t)` along
its within-step trajectory rather than sampling only the step endpoint. Segment nucleation events are
keyed to node/segment identity and continuous time, not step count. Validate with steps longer than
the nucleation, propagation and healing durations, and compare accumulated strain, rift length and
thinning, not only final topology.

**Conservation and ownership.** Store damage as an intensive state with an explicit lifecycle through
split, merge, remap, node activation/deactivation and save/load. Default missing state in old saves
to zero. A failed local cut does not itself duplicate or delete crust: thinning, opening and new crust
continue through one existing material-budget owner per step. Damage is a state variable, not an
extra crust-volume reservoir. Check crust/mantle ledgers and spatial contiguity independently:
conservation alone does not validate a plausible rift.

**Runtime and saves.** The added work is expected to scale as O(N + E + S) per step for N plate nodes,
their E local-neighbor links, and S slab-body stations, plus a connectedness check when segments
change. Slab stations are bounded and owned by active trench segments; avoid all-pairs slab coupling.
No measured runtime is available before a prototype; benchmark at density 4 against the existing 2%
step-time budget and do not extrapolate from min-cut's per-attempt runtime. Old saves load without slab
or deformation history, initialized compatibly. New saves retain slab geometry, plate deformation,
damage and segment state needed to continue the same rift history deterministically.

**Knobs.** No numeric starting values are approved for the new constitutive parameters; selecting them
before validation would be false precision. The design revision must set each value, its physical
source or explicit design-choice status, and the sensitivity plan before implementation.

| Knob | Start | Source | Tuning plan |
|---|---|---|---|
| Extensional load and yield strength | Incremental largest positive principal membrane resultant paired with additional normal-fault capacity from explicit `K0` (gate 2 diagnostic, 2026-10-19) | Candidate integrated brittle–ductile differential-strength envelope; material and background state remain scenario assumptions | Use the matched replay as screening evidence; improve equilibrium and edge validation, then retain scenario brackets and justify any selected threshold |
| Failure strain (`epsilon_fail`) | Pending gate 3 | Calibration assumption | Sensitivity on controlled loads and candidate events |
| Damage viscosity (`eta_rift`, Pa·s) | Pending gate 3 | Calibration assumption | Sweep strain rates and timestep partitions |
| Residual strength (`q_res`, [0, 1]) | Pending gate 3 | Design choice | Sweep weak-belt/craton cases; check localization and runaway |
| Healing time (Myr) | Pending gate 3 | Calibration assumption; may be omitted | Check arrested-segment persistence and reactivation |
| Tip propagation/link thresholds | Pending gate 3 | Design choice; use km and Myr | Tune on synthetic linkage cases, then report named-seed outcomes |
| Slab bending resistance | Pending coupled slab coupon | Calibration assumption; bracket effective flexural strength | Compare slab curvature and transmitted pull across the bracket |
| Slab-mantle drag | Pending coupled slab coupon | Calibration assumption; use the existing mantle drag law only where its geometry applies | Check sinking rate, force closure and plate response |
| Slab station spacing | Pending resolution study | Numerical choice; report km | Refine stations until length, curvature and force change within preregistered tolerances |
| Plate-shell effective modulus/yield | Pending finite-strain shell coupon | Material scenario | Compare extension, strain localization and craton diversion across material brackets |

## Detailed design (provisional touchpoints)

- `backend/app/torque.py`: `gather_per_node_force_field` exposes existing source forces in N
  before torque aggregation. `backend/app/thin_shell.py` supplies a reduced linear membrane
  equilibrium and exposes the largest positive principal stress increment for comparison with an
  integrated additional normal-fault capacity in N/m. `bin/debug/diagnose_gate2_column_envelope.py`
  is a standalone scenario diagnostic; a production column-envelope module remains a proposed
  touchpoint, not implemented. Its spherical benchmark and free-edge single-plate baseline are
  checked by `bin/debug/validate_thin_shell_plate_meshes.py`; coupled interplate boundaries and
  nonlinear load redistribution remain outside this diagnostic.
- `backend/app/world.py` and a proposed slab-body module: own persistent, three-dimensional slab
  material state and advance it with subduction, sinking, bending and mantle resistance. Derive
  trench pull from the body and its attachment; do not add that force on top of
  `torque.slab_pull_torque` for the same consumed plate material.
- `backend/app/thin_shell.py` is a linear displacement/stress diagnostic only. A production shell
  response must remove rigid motion already handled by `torque.py`, integrate finite in-plane strain
  in the co-rotating frame, and update/remap the sparse quadrilateral geometry and column thickness
  without losing cell area or ownership. Reuse the diagnostic for benchmark comparisons, not as a
  substitute for that evolving geometry path.
- `backend/app/merge_split.py`: move rift-state advancement out of the all-at-once split decision;
  retain a split/topology check that runs only after local failed links form a viable connected
  boundary. Existing flow fit may remain a diagnostic or a candidate-event signal, not the rift path.
- `backend/app/faults.py`, `backend/app/quad_tectonics.py` and `backend/app/rheology.py`: reuse or
  extend local fault/extension response so segment damage, slip and thinning share one constitutive
  response and ownership path. Exact function names are to be confirmed after diagnosis gate 2.
- `backend/app/sparse_quad_patch.py` and `backend/app/plates.py`: persist/remap the new intensive
  state with area-aware ownership; missing fields in old saves initialize to zero.
- `backend/app/cratons.py`: retain separate line-crossing and band-exposure measures and add segment
  growth, arrest, linkage, crust-type and damage diagnostics. Do not redefine B as the sole physical
  success measure.
- Tests: timestep-partition invariance; damage bounds and no damage under compression; localized
  weakening; tip advance/arrest; segment linkage via distributed shear; oceanic propagation; split
  only after a connected viable boundary; remap/save compatibility; ledger closure.
- Docs: update `docs/simulation-model/tectonics.md` and the architecture module map when the model
  design is accepted for implementation.

## Assumptions and uncertainty

- **The strength model is an explicit physical hypothesis, not a resolved lithosphere.** A first
  candidate combines frictional normal-fault failure with temperature- and strain-rate-dependent
  creep, then integrates differential strength with depth. It omits fluid state, vertical
  composition, inherited faults, anisotropy, grain-size evolution and 3-D load transfer. Existing
  Moho temperature and boundary Mohr-Coulomb proxies do not fill those gaps. Weak/wet and strong/dry
  brackets may show sensitivity but cannot identify the real material state.
- **Local loading is the key unresolved interface.** Node-wise force inputs are currently reduced
  into a whole-plate torque; neither that torque nor dimensionless `Plate.internal_stress`
  uniquely specifies intraplate stress. The replacement must derive a supported local load or stop
  and revise; it must not disguise the min-cut flow partition as mechanics.
- **Damage proxies brittle/ductile weakening.** A single strain/damage scalar omits temperature,
  fluids, grain-scale rheology and depth-dependent faulting. It is an efficient model proxy, not a
  resolved fault-zone simulation.
- **Segment linkage is not always a discrete transform.** Rift transfer may be distributed or
  extension-dominated; the model must not force a strike-slip “scissor” at every offset.
- **Craton avoidance is probabilistic geology, not an absolute rule.** Strong cratons should resist
  extension, while sufficiently large loads or inherited weaknesses may still rupture them. The
  acceptance target measures a comparative model behavior, not a universal geological law.
- **Continental and oceanic rifting differ.** A boundary may propagate from continental extension
  into oceanic crust, but the two modes need not share identical thresholds.
- **The 1 Gyr replays are stress tests, not calibration intervals.** Tune constitutive assumptions on
  shorter controlled cases, then use long runs to check side effects and budgets.

## Acceptance criteria

1. [x] **Gate 1:** use the repeatable paired runs as the baseline and retire the provenance-
       incomplete archives from calibration (2026-10-19; see `analysis/issue487/README.md`).
   [x] **Gate 2 — validated conditional diagnostic (2026-10-22):** column-envelope scenarios, the
       incremental shell stress measure, curved-sphere benchmark, original curved-mesh checks and
       isolated-plate free-edge policy are documented.
       Scope is craton-crossing events whose source cells form one exact edge-connected component.
       All 74/74 eligible source-connected candidates pass the 95% cut-band and craton-band node
       coverage checks, stop-code and 1% residual screens (maximum 0.4514%); this is 100% overall
       and in every nonempty seed, event-kind and time stratum, above the 90% overall / 80% stratum
       minimum. The 56/130 genuinely disconnected-source events are counted and reported separately,
       outside the Gate 2 denominator. Their unequal distribution means connected-source conclusions
       are conditional and cannot be generalized to all events. Capacity scores over 22,896 accepted
       cut-band samples remain screening evidence; the oceanic-age proxy is coarse. The
       conditional 100 km rift-like and 150 km thick-craton capacity anchors and their factorial
       uncertainty are recorded in `analysis/issue487/gate2-capacity-calibration.md`; they do not
       use old-algorithm cuts as calibration labels. A thickness-aware diagnostic mapping using
       local Hc/Hm and proxy fields with steady-versus-lagged temperature assumptions is complete;
       its output remains a scenario screen because profile and rheology state are unresolved.
       The manufactured shared-node benchmark passes. A corrected six-event sample covers one
       connected-source split and failed-rift per seed; it records 5.3–97.0% band-scale load
       sensitivity and 2.0–21.4% force / 4.6–21.6% torque closure residuals for active mutual
       contact drag. The tied/no-slip comparison is an idealized bonded endmember, not a selected
       finite-slip law; the previous two-event relative-change summaries were invalid where baseline
       node stress was zero. These pair results are excluded from Gate 2 acceptance. A future
       coupled-stress claim requires a consistent interface model, rigid-mode projection and traction
       ownership; the projection audit
       in `analysis/issue487/gate2-rigid-projection-audit.md` does not select one. Capacity and source
       assumptions remain conditional and are reported as uncertainty scenarios; Gate 2 does not
       require selecting a production capacity field, rupture label or damage threshold (see
       `analysis/issue487/README.md`, `analysis/issue487/gate2-load-diagnosis.md`,
       `analysis/issue487/gate2-column-envelope.json`, `analysis/issue487/gate2-matched-replay.json`,
       `analysis/issue487/gate2-spatial-capacity-replay.json`, and
       `analysis/issue487/gate2-capacity-sensitivity.json`,
       `analysis/issue487/gate2-coupled-shell-benchmark.json`,
       `analysis/issue487/gate2-interface-sensitivity.json`, and
       `analysis/issue487/gate2-rigid-projection-audit.md` and
       `analysis/issue487/gate2-rigid-projection-audit.json`).
2. [x] **Gate 3a — local-scale evidence phase closed by owner (2026-10-30).** The local work is
       recorded in `analysis/issue487/README.md`,
       the [Gate 3a acceptance matrix](../analysis/issue487/gate3-acceptance-matrix.md). It includes
       constitutive and conservation coupons; weak-seed versus intact-control initiation evidence
       (3 MPa sensitivity initiates, 14 MPa control does not in 9 tested flaws); local finite-strain,
       layered brittle–ductile and moving-tip growth/arrest diagnostics; and matched finite-source
       stress-transfer evidence. These results are reference evidence, not calibrated natural-rift
       parameters or a predictive plate model. Known limitations are retained: some local outcomes
       remain mesh/assumption sensitive; the recent finite-source host test has no damage outside its
       prescribed seed and does not advance its fixed crack; and its pressure-work comparison omits
       boundary and inelastic terms. Those limitations do not keep Gate 3a open; plate-scale
       consequences are Gate 3b.
3. [ ] **Gate 3b — rift-potential existence and source (tracked by #699):** publish the event
       definition and source/load accounting before probing. A slab-pull case must use evolving 3-D
       slab bodies coupled to in-plane deformation of the plate; static trench forces on a rigid plate
       do not qualify. In a controlled opposed-trench case, make the slab traction vectors point
       apart across the plate interior; characterize both slab bodies' length, thickness, depth, dip
       and curvature over time; show force transfer from each slab into the
       plate, and show the plate's finite strain and geometry evolve. No-slab and detached-slab
       controls must remove or reduce that response. Bending resistance and mantle drag must oppose
       motion with force, moment and work residuals reported. The hinge interaction must be equal and
       opposite on slab and plate, while buoyancy and mantle reaction remain separately accounted.
       Then identify at least one broad,
       sustained extension event on a generated spherical plate and couple it to Gate 3a local crack
       rules to produce a viable rift before splitting. Test a 10–30 Myr duration bracket as a
       hypothesis, alongside timestep partitioning; do not turn this span into an arbitrary timer.
       Record footprint, direction, extension/load relative to local capacity, force/torque/work
       residuals, crack response and slab geometry/load history. Exact operational cutoffs are to be
       preregistered in #699 before the probe; a missing event or failed accounting requires revising
       the source/event definition before broad history runs.
4. [ ] **Gate 3b plate-scale breakup:** on real spherical plate meshes at production resolution and
       timestep, matched whole-planet histories spanning at least 200 Myr demonstrate that rifting
       initiates preferentially in weak lithosphere, remains suppressed or diverted by strong/cratonic
       regions, and forms a connected viable boundary before a plate split. Report failed and
       successful rift systems separately, including duration, length, crust type, segment
       linkage/arrest, daughter-plate sizes, craton intersections and slab geometry/load history.
       No instantaneous plate-wide
       split may stand in for propagation. Repeat decisive comparisons at finer surface resolution
       and with a halved timestep; qualitative outcomes and named split/rift measures must remain
       within preregistered tolerances. Material advection, crust/mantle ledgers and boundary work
       must close, and runtime must be reported against the production budget. Local coupons can
       verify pieces of this closure but do not satisfy this criterion.
5. [ ] In synthetic tests, damage is timestep-partition invariant, bounded, absent under compression
       or sub-yield loading, and advances only through adjacent cells with supported positive loading.
6. [ ] Segments can nucleate independently, propagate, arrest and link; a transform/transfer zone
       transfers shear without making an instantaneous plate-wide split.
7. [ ] A plate splits only after a connected viable rift boundary forms; split frequency and mean plate
       count remain within 10% of the matched baseline unless the RFC is revised with a justified new
       target. Also report successful continental split-system counts per 200 Myr against the broad
       6–8-system Earth-history context in Assumptions and uncertainty; this is a plausibility
       diagnostic, not an additional pass threshold.
8. [ ] On seeds 7, 11 and 2024 over 1 Gyr, the successful split line-craton rate and band-craton
       exposure each fall at least 50% against a matched old-algorithm replay, with continental,
       oceanic and mixed crust reported separately. These are project acceptance targets, not claims
       that natural rifts obey a universal 50% rule.
9. [ ] Rifting-related craton loss is at most 2% of live craton per Gyr on average across seeds 7, 11
       and 2024, and in each of the first three Gyr on seed 456862899, subject to matched replay and
       ledger checks.
10. [ ] Crust/mantle ledgers balance; old saves load; deterministic event identities and timestep tests
       pass; runtime overhead is measured at density 4 and remains within the existing 2% step-time
       budget or the RFC is revised with a measured trade-off. Old saves initialize slab and
       deformation history safely; new saves preserve both for deterministic continuation.
11. [ ] **Gate 3b conditional mantle/magma attribution at plate scale:** preserve a tectonic-only or
       mantle-independent control. Where plume or magma forcing is represented, compare matched
       geological-duration whole-sphere histories with and without that forcing. Attribute changes
       separately in rift initiation, growth, linkage, and successful split-system counts, with
       source volume and force/torque budgets closed. A plume-assisted case may increase rifting, but
       the model must also permit rifts without plumes. Local finite-source coupons and their stress
       transfer measures are supporting verification only. Do not infer plate-scale causal benefit
       from a prescribed crack-tip test or tune forcing to match a single observed endpoint.

## Falsification criteria

- **No defensible local loading field:** the 2026-10-19 trace found current force/torque data
  insufficient to support a spatial extensional load without an added transfer/equilibrium rule.
  Do not calibrate friction, creep or fracture from these data. Extend the force interface and repeat
  the diagnosis, or revise/reject the physical local-mechanics premise in favor of a clearly labeled
  phenomenological rift model.
- **Slab force does not deform the plate:** if slab buoyancy is not transmitted through an attached
  body, bending/drag do not oppose the pull, or plate geometry and finite strain remain unchanged
  under the resulting load, slab-pull rifting is unsupported. Revise the coupling or reject that
  source attribution; do not use a fixed-geometry crack screen as evidence of breakup.
- **Local deformation does not localize:** if weak inherited zones do not attract slip and fracture
  under controlled loading, or deformation spreads diffusely regardless of strength, revise the
  constitutive law; do not tune only the craton-resistance parameter.
- **Segments do not produce viable boundaries:** if propagation/linkage cannot produce a connected
  separator without violating daughter size, fit or material budgets, reject the proposed topology
  transition and revisit the segment model.
- **No improvement against baseline:** if line-craton and band-craton measures miss the 50% target
  on any named seed after the approved parameter sensitivity range, reject or revise the mechanism;
  do not claim acceptance from improved B alone.
- **Timestep dependence or runaway deformation:** if event timing, stress, thinning or rift geometry
  changes materially under equivalent step partitions/refinements, or softening causes uncontrolled
  rupture, the integration/constitutive law is falsified and must be revised; do not hide it with a
  strain or stress clamp.
- **Unphysical magma-assisted growth:** if a dike proxy fails force/moment/work closure, exceeds its
  independent aperture or magma-volume bound, creates fracture without a supplied source, or fails to
  show the expected local stress relief alongside tip stress transfer, revise or reject the
  intrusion coupling; do not replace it with a fitted propagation multiplier.

## References

- <a id="ref-wilson1992"></a>Wilson, C. D. (1992). *Linear Elastic Fracture Mechanics Primer*.
  NASA Technical Memorandum 103591. [NASA Technical Reports Server](https://ntrs.nasa.gov/citations/19920021173).
- <a id="ref-schultz2008"></a>Schultz, R. A., Mège, D., & Diot, H. (2008). Emplacement
  conditions of igneous dikes in Ethiopian Traps. *Journal of Volcanology and Geothermal Research*,
  178(4), 683–692. [doi:10.1016/j.jvolgeores.2008.08.012](https://doi.org/10.1016/j.jvolgeores.2008.08.012)
- <a id="ref-carpenter2022"></a>Carpenter, B. M., Williams, J. N., Fagereng, Å., Wedmore, L. N. J.,
  Biggs, J., Mphepo, F., Mdala, H., Dulanya, Z., & Manda, B. (2022). Comparing intrarift and border
  fault structure in the Malawi Rift: Implications for normal fault growth. *Journal of Structural
  Geology*, 165, 104761. [doi:10.1016/j.jsg.2022.104761](https://doi.org/10.1016/j.jsg.2022.104761)
- <a id="ref-vachon2019"></a>Vachon, R., & Hieronymus, C. F. (2019). Mechanical energy balance and
  apparent fracture toughness for dykes in elastoplastic host rock with large-scale yielding.
  *Geophysical Journal International*, 219(3), 1786–1804.
  [doi:10.1093/gji/ggz383](https://doi.org/10.1093/gji/ggz383)
- <a id="ref-hatton1994"></a>Hatton, C. G., Main, I. G., & Meredith, P. G. (1994). Non-universal
  scaling of fracture length and opening displacement. *Nature*, 367, 160–162.
  [doi:10.1038/367160a0](https://doi.org/10.1038/367160a0)
- <a id="ref-buck2006"></a>Buck, W. R., Einarsson, P., & Brandsdóttir, B. (2006). Tectonic
  stress and magma chamber size as controls on dike propagation: Constraints from the 1975–1984
  Krafla rifting episode. *Journal of Geophysical Research: Solid Earth*, 111, B12404.
  [doi:10.1029/2005JB003879](https://doi.org/10.1029/2005JB003879)
- <a id="ref-sigmundsson2015"></a>Sigmundsson, F., et al. (2015). Segmented lateral dyke growth
  in a rifting event at Bárðarbunga volcanic system, Iceland. *Nature*, 517, 191–195.
  [doi:10.1038/nature14111](https://doi.org/10.1038/nature14111)
- <a id="ref-maccaferri2016"></a>Maccaferri, F., E. Rivalta, & Y. Aoki (2016). On the mechanisms
  governing dike arrest: Insight from the 2000 Miyakejima dike injection. *Earth and Planetary
  Science Letters*, 434, 64–74. [doi:10.1016/j.epsl.2015.11.024](https://doi.org/10.1016/j.epsl.2015.11.024)
- <a id="ref-ueda2005"></a>Ueda, H., E. Fujita, M. Ukawa, E. Yamamoto, M. Irwan, & F. Kimata
  (2005). Magma intrusion and discharge process at the initial stage of the 2000 activity of
  Miyakejima, Central Japan, inferred from tilt and GPS data. *Geophysical Journal International*,
  161(3), 891–906. [doi:10.1111/j.1365-246X.2005.02602.x](https://doi.org/10.1111/j.1365-246X.2005.02602.x)
- <a id="ref-corti2023"></a>Corti, N., et al. (2023). Feeders vs arrested dikes: A case study from
  the Younger Stampar eruption in Iceland. *Journal of Volcanology and Geothermal Research*, 443,
  107914. [doi:10.1016/j.jvolgeores.2023.107914](https://doi.org/10.1016/j.jvolgeores.2023.107914)
- <a id="ref-konstantinou2024"></a>Konstantinou, K. I. (2024). Dike volume derived from seismicity
  as a gauge of fracture toughness and propagation dynamics. *Scientific Reports*, 14, 17593.
  [doi:10.1038/s41598-024-67724-0](https://doi.org/10.1038/s41598-024-67724-0)
- <a id="ref-geshi2020"></a>Geshi, N., J. Browning, & S. Kusumoto (2020). Magmatic overpressures,
  volatile exsolution and potential explosivity of fissure eruptions inferred via dike aspect ratios.
  *Scientific Reports*, 10, 9406. [doi:10.1038/s41598-020-66226-z](https://doi.org/10.1038/s41598-020-66226-z)

- <a id="ref-grandin2009"></a>Grandin, R., et al. (2009). September 2005 Manda Hararo-Dabbahu
  rifting event, Afar (Ethiopia): Constraints provided by geodetic data. *Journal of Geophysical
  Research: Solid Earth*, 114, B08404. [doi:10.1029/2008JB005843](https://doi.org/10.1029/2008JB005843)
- <a id="ref-grandin2010"></a>Grandin, R., A. Socquet, E. Jacques, N. Mazzoni, J.-B. de
  Chabalier, and G. C. P. King (2010). Sequence of rifting in Afar, Manda-Hararo rift, Ethiopia,
  2005-2009: Time-space evolution and interactions between dikes from interferometric synthetic
  aperture radar and static stress change modeling. *Journal of Geophysical Research: Solid Earth*,
  115, B10413. [doi:10.1029/2009JB000815](https://doi.org/10.1029/2009JB000815)
- <a id="ref-grandin2011"></a>Grandin, R., et al. (2011). Seismicity during lateral dike
  propagation: Insights from new data in the recent Manda Hararo-Dabbahu rifting episode (Afar,
  Ethiopia). *Geochemistry, Geophysics, Geosystems*, 12, Q0AB08.
  [doi:10.1029/2010GC003434](https://doi.org/10.1029/2010GC003434)
- <a id="ref-belachew2011"></a>Belachew, M., C. J. Ebinger, D. Coté, D. Keir, J. V. Rowland,
  J. O. S. Hammond, and A. Ayele (2011). Comparison of dike intrusions in an incipient
  seafloor-spreading segment in Afar, Ethiopia: Seismicity perspectives. *Journal of Geophysical
  Research: Solid Earth*, 116, B06405.
  [doi:10.1029/2010JB007908](https://doi.org/10.1029/2010JB007908)

- <a id="ref-bendick2006"></a>Bendick, R., S. McClusky, R. Bilham, L. Asfaw, and S. Klemperer (2006). Distributed Nubia–Somalia relative motion and dike intrusion in the Main Ethiopian Rift. *Geophysical Journal International*, 165(1), 303–310. [doi:10.1111/j.1365-246X.2006.02904.x](https://doi.org/10.1111/j.1365-246X.2006.02904.x)

- <a id="ref-kogan2012"></a>Kogan, L., S. Fisseha, R. Bendick, R. Reilinger, S. McClusky, R. King, and T. Solomon (2012). Lithospheric strength and strain localization in continental extension from observations of the East African Rift. *Journal of Geophysical Research: Solid Earth*, 117, B03402. [doi:10.1029/2011JB008516](https://doi.org/10.1029/2011JB008516)

- <a id="ref-simo1988"></a>Simo, J. C. (1988). A framework for finite strain elastoplasticity
  based on maximum plastic dissipation and the multiplicative decomposition: Part I. Continuum
  formulation. *Computer Methods in Applied Mechanics and Engineering*, 66(2), 199–219.
  [doi:10.1016/0045-7825(88)90104-1](https://doi.org/10.1016/0045-7825(88)90104-1)

- <a id="ref-allken2012"></a>Allken, V., Huismans, R. S., & Thieulot, C. (2012). Factors controlling the mode of rift interaction in brittle–ductile coupled systems: A 3D numerical study. *Geochemistry, Geophysics, Geosystems*, 13, Q05010. [doi:10.1029/2012GC004077](https://doi.org/10.1029/2012GC004077)
- <a id="ref-glerum2018"></a>Glerum, A., Thieulot, C., Fraters, M., Blom, C., & Spakman, W. (2018). Nonlinear viscoplasticity in ASPECT: Benchmarking and applications to subduction. *Solid Earth*, 9, 267–294. [doi:10.5194/se-9-267-2018](https://doi.org/10.5194/se-9-267-2018)

- <a id="ref-brace1980"></a>Brace, W. F., & Kohlstedt, D. L. (1980). Limits on lithospheric
  stress imposed by laboratory experiments. *Journal of Geophysical Research: Solid Earth*,
  85(B11), 6248–6252. [doi:10.1029/JB085iB11p06248](https://doi.org/10.1029/JB085iB11p06248)
- <a id="ref-kohlstedt1995"></a>Kohlstedt, D. L., Evans, B., & Mackwell, S. J. (1995). Strength
  of the lithosphere: Constraints imposed by laboratory experiments. *Journal of Geophysical
  Research: Solid Earth*, 100(B9), 17587–17602.
  [doi:10.1029/95JB01460](https://doi.org/10.1029/95JB01460)
- <a id="ref-hirth2003"></a>Hirth, G., & Kohlstedt, D. (2003). Rheology of the upper mantle and
  the mantle wedge: A view from the experimentalists. In *Inside the Subduction Factory* (pp.
  83–105). American Geophysical Union. [doi:10.1029/138GM06](https://doi.org/10.1029/138GM06)
- <a id="ref-byerlee1978"></a>Byerlee, J. D. (1978). Friction of rocks. *Pure and Applied
  Geophysics*, 116, 615–626. [doi:10.1007/BF00876528](https://doi.org/10.1007/BF00876528)
- <a id="ref-chen2016"></a>Chen, L., & Gerya, T. V. (2016). The role of lateral lithospheric
  strength heterogeneities in orogenic plateau growth: Insights from 3-D thermo-mechanical
  modeling. *Journal of Geophysical Research: Solid Earth*, 121, 3118–3138.
  [doi:10.1002/2016JB012872](https://doi.org/10.1002/2016JB012872)
- <a id="ref-wang2012"></a>Wang, Y. F., Zhang, J. F., Jin, Z. M., & Green, H. W., II (2012).
  Mafic granulite rheology: Implications for a weak continental lower crust. *Earth and Planetary
  Science Letters*, 353–354, 99–107. [doi:10.1016/j.epsl.2012.08.004](https://doi.org/10.1016/j.epsl.2012.08.004)
- <a id="ref-stamps2014"></a>Stamps, D. S., Flesch, L. M., Calais, E., & Ghosh, A. (2014).
  Current kinematics and dynamics of Africa and the East African Rift System. *Journal of
  Geophysical Research: Solid Earth*, 119(6), 5161–5186.
  [doi:10.1002/2013JB010717](https://doi.org/10.1002/2013JB010717)
- <a id="ref-rajaonarison2021"></a>Rajaonarison, T. A., et al. (2021). Role of lithospheric buoyancy
  forces in driving deformation in East Africa from 3D geodynamic modeling. *Geophysical Research
  Letters*, 48, e2020GL090483. [doi:10.1029/2020GL090483](https://doi.org/10.1029/2020GL090483)
- <a id="ref-duretz2025"></a>Duretz, T., et al. (2025). Lithospheric deformation with mechanical
  anisotropy: A numerical model and application to continental rifting. *Geochemistry, Geophysics,
  Geosystems*. [doi:10.1029/2025GC012409](https://doi.org/10.1029/2025GC012409)

- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., et al. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
- <a id="ref-forsyth1975"></a>Forsyth, D. W., & Uyeda, S. (1975). On the relative importance
  of the driving forces of plate motion. *Geophysical Journal of the Royal Astronomical Society*,
  43(1), 163–200. [doi:10.1111/j.1365-246X.1975.tb00631.x](https://doi.org/10.1111/j.1365-246X.1975.tb00631.x)
- <a id="ref-li2023"></a>Li, Y., & Gurnis, M. (2023). A simple force balance model of subduction
  initiation. *Geophysical Journal International*, 232(1), 128–146.
  [doi:10.1093/gji/ggac332](https://doi.org/10.1093/gji/ggac332)
- <a id="ref-naliboff2020"></a>Naliboff, J. B., Glerum, A., Brune, S., Péron-Pinvidic, G., & Wrona, T.
  (2020). Development of 3-D Rift Heterogeneity Through Fault Network Evolution. *Geophysical
  Research Letters*, 47(13), e2019GL086611.
  [doi:10.1029/2019GL086611](https://doi.org/10.1029/2019GL086611)
- <a id="ref-kolawole2024"></a>Kolawole, F., & Ajala, R. (2024). Propagating rifts: the roles of
  crustal damage and ascending mantle fluids. *Solid Earth*, 15, 747–762.
  [doi:10.5194/se-15-747-2024](https://doi.org/10.5194/se-15-747-2024)
- <a id="ref-pagli2019"></a>Pagli, C., Yun, S.-H., Ebinger, C., Keir, D., & Wang, H. (2019).
  Strike-slip tectonics during rift linkage. *Geology*, 47(1), 31–34.
  [doi:10.1130/G45345.1](https://doi.org/10.1130/G45345.1)
- <a id="ref-audet2011"></a>Audet, P., & Bürgmann, R. (2011). Dominant role of tectonic
  inheritance in supercontinent cycles. *Nature Geoscience*, 4, 184–187.
  [doi:10.1038/ngeo1080](https://doi.org/10.1038/ngeo1080)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-heron2016"></a>Heron, P. J., Pysklywec, R. N., & Stephenson, R. (2016).
  Lasting mantle scars lead to perennial plate tectonics. *Nature Communications*, 7,
  11834. [doi:10.1038/ncomms11834](https://doi.org/10.1038/ncomms11834)
- <a id="ref-koptev2015"></a>Koptev, A., Calais, E., Burov, E., Leroy, S., & Gerya, T.
  (2015). Dual continental rift systems generated by plume–lithosphere interaction.
  *Nature Geoscience*, 8, 388–392. [doi:10.1038/ngeo2401](https://doi.org/10.1038/ngeo2401)
- <a id="ref-lenardic2003"></a>Lenardic, A., Moresi, L.-N., & Mühlhaus, H. (2003).
  Longevity and stability of cratonic lithosphere: Insights from numerical simulations of
  coupled mantle convection and continental tectonics. *Journal of Geophysical Research:
  Solid Earth*, 108(B6), 2303.
  [doi:10.1029/2002JB001859](https://doi.org/10.1029/2002JB001859)
- <a id="ref-tommasi2001"></a>Tommasi, A., & Vauchez, A. (2001). Continental rifting
  parallel to ancient collisional belts: An effect of the mechanical anisotropy of the
  lithospheric mantle. *Earth and Planetary Science Letters*, 185, 199–210.
  [doi:10.1016/S0012-821X(00)00350-2](https://doi.org/10.1016/S0012-821X(00)00350-2)
- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then
  re-open? *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)
