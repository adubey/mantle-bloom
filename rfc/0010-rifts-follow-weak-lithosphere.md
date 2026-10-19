# RFC 0010: Rifts follow weak lithosphere

- **Status:** Accepted
- **Date:** 2026-10-15
- **Last revised:** 2026-10-18 after #480 routing diagnostics
- **Related:** #481 (tracking), #480 (implementation), #472 / PR #479 (RFC 0009 G3), #389 (RFC 0003 tracking), #378, #385, #388 (RFC 0003 gates), PR #418 (RFC
  0003 solution 1); RFCs [0003](0003-rifting-around-cratons.md) (superseded by RFC 0010),
  [0009](0009-craton-mosaic.md) (cratons-first layout; its G3 gate selects RFC 0010),
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
[[Lenardic et al., 2003]](#ref-lenardic2003) [[Foley, 2008]](#ref-foley2008). Continental rifts
also depend on interacting forces and evolving lithospheric strength. Rift systems commonly begin
as separate fault segments that propagate, interact and link; progressive thinning changes the
force balance, and some rifts arrest before breakup [[Brune et al., 2023]](#ref-brune2023)
[[Kolawole & Ajala, 2024]](#ref-kolawole2024). The Afar linkage observations show that transfer
between rift tips can involve distributed extension and rift-perpendicular shear; a transform-like
link is not necessarily a single strike-slip fault that slices a plate [[Pagli et al., 2018]](#ref-pagli2018).
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
6. **Incremental, damage-driven rift segments (new; recommended).** Use local extensional loading
   and strength to accumulate irreversible strain/damage through time. Damage nucleates discrete
   rift segments, their tips advance into adjacent material when local loading can sustain growth,
   and nearby segments link through extensional and transfer/shear deformation. A plate splits only
   after a damaged, sufficiently connected boundary separates it; before then it may continue,
   migrate, or arrest. Transform motion transfers displacement between offset segments; it is not
   an independent whole-plate cut.

**Recommendation (revised after #480 evidence, 2026-10-18).** Treat solution 1 as failed for this
problem. Replace it with solution 6, developed with RFC 0012's process-specific force-response
interfaces. The nucleate-and-propagate proposal anticipates local growth, but does not define the
state and linkage rules needed for local damage accumulation. Retain existing plate-scale
flow/force information as a possible source of loading, not as a command to draw a complete cut.
Keep solution 4 separately
sized against the craton-loss gate, and take plume weakening through RFC 0007.

## Design

### Diagnosis gates

1. **Reconcile the diagnostic event stream before calibrating.** Re-run the recorded seeds and
   compare event signatures, counts and admissibility margins from the production tap and shadow
   calculator. The detailed margin reruns and the earlier production table have different event
   totals despite nominally matching seed, layout, step and duration. Preserve each output; do not
   combine denominators until the cause is understood. This does not rescue min-cut's observed
   shortfall, but it is necessary for a reproducible baseline for the replacement.
2. **Establish the local load available to damage.** Trace RFC 0012's torque/process interfaces and
   determine whether a spatially varying extensional traction or strain-rate field can be derived
   from current loads without inventing one from rigid plate rotation. Record its direction and
   magnitude in Pa (stress), N/m (vertically integrated force per boundary length), or s⁻¹ (strain
   rate), and compare the integrated tensile load with vertically integrated lithosphere strength.
   If the current force interface only yields whole-plate torque, extend that interface first; do not
   use the min-cut's two flow clusters as a surrogate stress field.
3. **Calibrate accumulation and segment growth.** On synthetic plates and replayed candidate
   events, measure damage/strain accumulation, nucleation, tip propagation, arrest and linkage under
   timestep refinement. Select thresholds only after checking physical scale and sensitivity. If no
   local load can be supported from the model's forces and material state, revise the premise rather
   than add a free-standing damage timer.

### Mechanics (provisional)

- **State and loading.** Add a per-node cumulative irreversible extensional strain
  `epsilon_rift` (dimensionless, initially zero and capped at `epsilon_fail`). Obtain local tensile
  stress `sigma_t` in Pa and local yield strength `Y_t` in Pa from a reduced-order load/strength
  calculation tied to RFC 0012. The current torque balance drives whole-plate motion; its loads do
  not by themselves define intraplate stress, so a spatial loading distribution is an explicit
  prerequisite, not an assumed by-product.
- **Damage update.** Where the tensile load exceeds local yield, integrate the inelastic extension
  rate over the portion of the step in which it is positive:
  `d epsilon_rift / dt = max(0, sigma_t - Y_t) / eta_rift`, with `eta_rift` in Pa·s. Keep
  `epsilon_rift` stays in `[0, epsilon_fail]`; convert it to a bounded damage fraction
  `D = clip(epsilon_rift / epsilon_fail, 0, 1)`. A provisional bounded softening form is
  `Y_t(D) = Y_t(0) * [1 - D * (1 - q_res)]`, where the residual-strength fraction `q_res` is in
  [0, 1]. Reaching `D = 1` marks local failure; this creates/extends a segment but does not split the
  plate unless the failed links form a viable connected boundary. `epsilon_fail`, `eta_rift`,
  `q_res`, any healing time, and whether this softening form is adequate are calibration assumptions
  with no approved numeric values yet. They must be set after diagnosis gate 3, with physical scale
  checks and sensitivity runs; no parameter should be selected to hit the craton-crossing percentage
  alone.
- **Nucleation and growth.** Nucleate a finite segment when accumulated strain and local loading
  satisfy the calibrated failure rule. Extend a segment tip only into adjacent cells where the local
  tensile load can sustain damage; permit multiple segments and abandonment. Link overlapping or
  interacting segments through distributed extension and transfer/shear zones. The graph identifies
  local adjacency and connectedness only; it does not search globally for the least-cost separator.
  Strike-slip/transform response transfers relative displacement between segments and may change tip
  loading, but cannot split the plate by itself.
- **Rift completion.** Create a new plate only when damage/extension forms a connected, viable boundary
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

**Conservation and ownership.** Store damage as an intensive per-node state with an explicit lifecycle
through split, merge, remap, node activation/deactivation and save/load. Default missing state in old
saves to zero. Rift thinning and new crust continue through one existing material-budget owner per
step; damage is a state variable, not an extra crust-volume reservoir. Check crust/mantle ledgers and
spatial contiguity independently: conservation alone does not validate a plausible rift.

**Runtime and saves.** The added work is expected to scale as O(N + E) per step for N plate nodes and
their E local-neighbor links, plus a connectedness check when segments change. No measured runtime is
available before a prototype; benchmark at density 4 against the existing 2% step-time budget and do
not extrapolate from min-cut's per-attempt runtime. Old saves load with zero damage/segment state. New
saves retain all state needed to continue the same rift history deterministically.

**Knobs.** No numeric starting values are approved for the new constitutive parameters; selecting them
before the load diagnosis would be false precision. The diagnosis-gate revision must set each value,
its physical source or explicit design-choice status, and the sensitivity plan before implementation.

## Detailed design (provisional touchpoints)

- `backend/app/torque.py` and RFC 0012's force-response interfaces: provide or aggregate a spatially
  resolved tensile load/strain-rate proxy; this may require a new reduced-order load distributor.
- `backend/app/merge_split.py`: replace the whole-plate min-cut choice and success coin flip with
  segment eligibility, state advancement and a connected-boundary completion check. Existing flow
  fit may remain a diagnostic/kinematic consistency check, not the rift path.
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

- **Local loading is the key unresolved interface.** Rigid plate rotation and whole-plate torque do
  not uniquely specify intraplate stress. The replacement must derive a supported local load or stop
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

1. [ ] **Diagnosis gates first:** reconcile the differing diagnostic event streams and publish a
       repeatable, matched baseline; provide a spatial load/strength derivation with units; set and
       sensitivity-test the constitutive parameters before full replays.
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

- **No defensible local loading field:** if current force/torque data cannot support a spatial
  extensional load without an arbitrary flow-derived line, stop this design and extend the force
  interface or retain a clearly labeled phenomenological rift model.
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

- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., & Stamps, D. S. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
- <a id="ref-forsyth1975"></a>Forsyth, D. W., & Uyeda, S. (1975). On the relative importance
  of the driving forces of plate motion. *Geophysical Journal of the Royal Astronomical Society*,
  43(1), 163–200. [doi:10.1111/j.1365-246X.1975.tb00631.x](https://doi.org/10.1111/j.1365-246X.1975.tb00631.x)
- <a id="ref-kolawole2024"></a>Kolawole, F., & Ajala, R. (2024). Propagating rifts: the roles of
  crustal damage and ascending mantle fluids. *Solid Earth*, 15, 747–762.
  [doi:10.5194/se-15-747-2024](https://doi.org/10.5194/se-15-747-2024)
- <a id="ref-pagli2018"></a>Pagli, C., Yun, S.-H., Ebinger, C., Keir, D., & Wang, H. (2018).
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
