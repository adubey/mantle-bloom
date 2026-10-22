# RFC 0010: Rifts follow weak lithosphere

- **Status:** Accepted
- **Date:** 2026-10-15
- **Last revised:** 2026-10-18 after #480 routing diagnostics; revised local-cut design accepted and tracked in #487
- **Related:** #487 (tracking), #481 (superseded tracker), #480 (min-cut diagnostic), #472 / PR #479 (RFC 0009 G3), #389 (RFC 0003 tracking), #378, #385, #388 (RFC 0003 gates), PR #418 (RFC
  0003 solution 1); RFCs [0003](0003-rifting-around-cratons.md) (superseded by RFC 0010),
  [0009](0009-craton-mosaic.md) (cratons-first layout; G3 initially selected the min-cut),
  [0007](0007-mantle-plumes-and-plateaus.md) (plume weakening),
  [0012](0012-physics-first-simulation.md) (force-response interfaces)

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
field. Gate 2 remains open
because the replay provides no independent rupture labels, local material state is unresolved,
and broader load/boundary validation remains outstanding.
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
2. **Gate 2 — establish local loading and strength — incremental extensional measure, column
   scenario diagnostic, curved-shell benchmark, curved-mesh checks and isolated-plate edge policy
   added; connected-source mesh repairs and a 2026-10-21 skip audit completed. Scope is craton-crossing
   candidates whose source cells form one exact edge-connected component. The 56 disconnected-source
   candidates are counted separately and excluded from the coverage denominator. The connected-source
   coverage requirement passes. A thickness-aware spatial capacity screen now uses per-node Hc/Hm,
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
   and 2,701 (11.80%) in at least one; node rank correlation remains at least 0.9936. The `K0`
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
   physical coupled predictions; the pair-level rigid-mode projection remains unaudited.
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
   load distributions remain outside this single-plate solve. Do not calibrate the damage law or
   infer local stress from applied traction, rigid rotation, flow clusters or `Plate.internal_stress`.
   See the audits in `analysis/issue487/README.md`, `analysis/issue487/gate2-load-diagnosis.md`,
   `analysis/issue487/gate2-capacity-calibration.md`, and `analysis/issue487/gate2-matched-replay.json`.
3. **Calibrate accumulation and segment growth.** On synthetic plates and replayed candidate
   events, measure damage/strain accumulation, nucleation, tip propagation, arrest and linkage under
   timestep refinement. Select thresholds only after checking physical scale and sensitivity. If no
   local load can be supported from the model's forces and material state, revise the premise rather
   than add a free-standing damage timer.

### Mechanics (provisional)

- **State and loading.** After diagnosis gate 2, store cumulative irreversible extensional strain
  (`epsilon_rift`, dimensionless) on the state element whose area and ownership are defined by the
  load discretization. Each local cut is an adjacency link between cells; a failed cell or link
  becomes a short rift segment, not a plate partition. The thin-shell diagnostic returns tensile
  membrane stress increment `N_ext = max(0, lambda_max(Delta_N))` and compares it with the column's
  additional normal-fault strength resultant `N_cap`, both in N/m. The diagnostic pairing is
  defined, but production use waits for Gate 2 validation and a justified capacity selection. If the
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
- **Gravity and mantle coupling.** Use loads already represented in the plate torque model—slab pull,
  ridge push and basal drag—as contributors only where a defensible spatial load distribution exists.
  Include lithospheric gravitational-potential-energy gradients when column thickness/topography can
  support a dimensionally consistent force. Mantle-flow vectors may assist or resist motion through
  traction; they are not a prescribed rift line. Do not count the same load twice between plate torque
  and segment stress.

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

**Runtime and saves.** The added work is expected to scale as O(N + E) per step for N plate nodes and
their E local-neighbor links, plus a connectedness check when segments change. No measured runtime is
available before a prototype; benchmark at density 4 against the existing 2% step-time budget and do
not extrapolate from min-cut's per-attempt runtime. Old saves load with zero damage/segment state. New
saves retain all state needed to continue the same rift history deterministically.

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

## Detailed design (provisional touchpoints)

- `backend/app/torque.py`: `gather_per_node_force_field` exposes existing source forces in N
  before torque aggregation. `backend/app/thin_shell.py` supplies a reduced linear membrane
  equilibrium and exposes the largest positive principal stress increment for comparison with an
  integrated additional normal-fault capacity in N/m. `bin/debug/diagnose_gate2_column_envelope.py`
  is a standalone scenario diagnostic; a production column-envelope module remains a proposed
  touchpoint, not implemented. Its spherical benchmark and free-edge single-plate baseline are
  checked by `bin/debug/validate_thin_shell_plate_meshes.py`; coupled interplate boundaries and
  nonlinear load redistribution remain outside this diagnostic.
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
   [ ] **Gate 2:** column-envelope scenarios, the incremental shell stress measure, curved-sphere
       benchmark, original curved-mesh checks and isolated-plate free-edge policy are documented.
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
       node stress was zero. Audit pair-level rigid-mode projection before treating coupled values
       as physical. Justify the remaining source/capacity assumptions before Gate 2 closes (see
       `analysis/issue487/README.md`, `analysis/issue487/gate2-load-diagnosis.md`,
       `analysis/issue487/gate2-column-envelope.json`, `analysis/issue487/gate2-matched-replay.json`,
       `analysis/issue487/gate2-spatial-capacity-replay.json`, and
       `analysis/issue487/gate2-capacity-sensitivity.json`,
       `analysis/issue487/gate2-coupled-shell-benchmark.json`, and
       `analysis/issue487/gate2-interface-sensitivity.json`).
2. [ ] In synthetic tests, damage is timestep-partition invariant, bounded, absent under compression
       or sub-yield loading, and advances only through adjacent cells with supported positive loading.
3. [ ] Segments can nucleate independently, propagate, arrest and link; a transform/transfer zone
       transfers shear without making an instantaneous plate-wide split.
4. [ ] A plate splits only after a connected viable rift boundary forms; split frequency and mean plate
       count remain within 10% of the matched baseline unless the RFC is revised with a justified new
       target.
5. [ ] On seeds 7, 11 and 2024 over 1 Gyr, the successful split line-craton rate and band-craton
       exposure each fall at least 50% against a matched old-algorithm replay, with continental,
       oceanic and mixed crust reported separately. These are project acceptance targets, not claims
       that natural rifts obey a universal 50% rule.
6. [ ] Rifting-related craton loss is at most 2% of live craton per Gyr on average across seeds 7, 11
       and 2024, and in each of the first three Gyr on seed 456862899, subject to matched replay and
       ledger checks.
7. [ ] Crust/mantle ledgers balance; old saves load; deterministic event identities and timestep tests
       pass; runtime overhead is measured at density 4 and remains within the existing 2% step-time
       budget or the RFC is revised with a measured trade-off.

## Falsification criteria

- **No defensible local loading field:** the 2026-10-19 trace found current force/torque data
  insufficient to support a spatial extensional load without an added transfer/equilibrium rule.
  Do not calibrate the damage law from these data. Extend the force interface and repeat the
  diagnosis, or revise/reject the physical local-damage premise in favor of a clearly labeled
  phenomenological rift model.
- **Damage does not localize:** if weak inherited zones do not attract segments under controlled
  loading, or damage spreads diffusely regardless of strength, revise the constitutive law; do not
  tune only the craton-resistance parameter.
- **Segments do not produce viable boundaries:** if propagation/linkage cannot produce a connected
  separator without violating daughter size, fit or material budgets, reject the proposed topology
  transition and revisit the segment model.
- **No improvement against baseline:** if line-craton and band-craton measures miss the 50% target
  on any named seed after the approved parameter sensitivity range, reject or revise the mechanism;
  do not claim acceptance from improved B alone.
- **Timestep dependence or runaway damage:** if event timing, final strain or rift geometry changes
  materially under equivalent step partitions, or damage softening causes uncontrolled rupture, the
  integration/constitutive law is falsified and must be revised.

## References

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
