# RFC 0012: Physics-first simulation

- **Status:** Draft
- **Date:** 2026-10-16
- **Related:** RFCs [0009](0009-craton-mosaic.md) (cratons-first layout),
  [0010](0010-rifts-follow-weak-lithosphere.md) (rift routing), and
  [0011](https://github.com/adubey/mantle-bloom/blob/claude/rfc-0011-contact-strength/rfc/0011-contact-strength.md)
  (collision resistance); design principles 1–5 and 8–9

## Problem

Mantle Bloom needs tectonic processes to respond to the physical causes that control their
outcomes, while keeping each step fast enough to run interactively and over long histories.
“Physics-first” does not mean resolving every force on every rock element. It means choosing
state variables, process interfaces, constitutive responses, and numerical approximations
from the physical question the model is meant to answer.

Plate motion already combines rigid-sphere kinematics with force-inspired terms. The torque
engine integrates slab pull, ridge push, basal drag, sunk-slab drag, and continental contact
drag into each plate's angular motion (`torque.py`). Other mechanics are represented in
separate paths: split eligibility uses velocity clustering and geometric cuts
(`merge_split.py`); boundary deformation classifies local convergence, divergence, and
transform motion (`quad_tectonics.py`, `rheology.py`); fault systems update local relief and
slip (`faults.py`). This is a useful reduced model, but a process can change geometry or
material without returning its resistance to the motion solve. Conversely, a force term can
change a plate's rigid rotation without resolving the internal stress and strain it would
cause in the lithosphere.

RFC 0011 makes this gap concrete for blocked convergence: its proposed yield-limited contact
force feeds resistance back into plate motion, and its yield branch allocates the convergence
to deformation. Rifting and transform sliding raise the same coupling question. Rifting
currently thins crust under divergent motion, while split initiation is selected by a
velocity-fit and geometric-cut procedure. Transform faults can accommodate local slip, but
their shear resistance is not an explicit term in the plate torque solve.

These are not arguments for replacing every process with a fine-grid continuum solve. Away
from boundaries, treating plates as rigid is a useful approximation; strong, broad interiors
and localized weak boundaries are central features of plate tectonics, though diffuse and
deforming plate regions show the approximation has limits [[Gordon, 1998]](#ref-gordon1998)
[[Mulyukova et al., 2019]](#ref-mulyukova2019). The design problem is to pass each process
enough information to calculate its physically relevant response, aggregate that response at
the scale needed by the global motion solve, and retain the invariants that matter.

**Design principles at issue:** process-specific physics at clear interfaces (principle 8),
stated model boundaries and scales (principle 1), efficient approximations (principle 2),
timestep robustness (principle 3), explicit units and conservation (principles 5–6), and
verification separated from validation (principle 9).

## Proposed solutions

**What happens in nature.** Plate interiors are often approximated as rigid blocks separated
by deforming boundaries. Their angular motion can be studied through a torque balance of
driving and resisting forces, including forces from slabs, ridges, the mantle, and plate
boundaries [[Forsyth & Uyeda, 1975]](#ref-forsyth1975). Where the lithosphere deforms, its
response depends on strength and rheology; that response can localize near boundaries rather
than requiring every plate interior to deform equally [[Mulyukova et al., 2019]](#ref-mulyukova2019).
These observations support a reduced model with rigid plate interiors and process-specific
boundary mechanics, while leaving open which forces and constitutive laws are adequate for
Mantle Bloom's scales.

1. **Continue with process-specific algorithms and add feedback case by case.** This has low
   up-front cost and fits the existing code. But each new mechanism risks inventing its own
   units, contact definition, averaging rule, and coupling order; the same physical quantity
   may be measured differently by motion, deformation, and topology code.
2. **Resolve continuum mechanics throughout the lithosphere.** A fine spatial model could
   derive deformation and plate boundaries from stress and rheology in one framework. It
   would be a major change in state, numerics, runtime, and calibration, and is not needed
   for every question the current simulator asks.
3. **Use physics-based process interfaces with reduced-order coupling (recommended).** Keep
   rigid plate interiors where appropriate. Represent mechanically active boundaries as
   segments or patches. Each process supplies a physical force-response law and its material
   response; the motion solver uses aggregated loads, while deformation receives the solved
   relative motion and transfers material through its own rules and ledgers. Use sufficient
   statistics and implicit or bounded solves where needed, rather than resolving all rock at
   fine scale.

**Recommendation:** adopt solution 3 for tectonic mechanics, incrementally. Treat RFC 0011 as
a candidate collision pilot, not as a blanket template: first express its contact behavior in
terms of the interface below, then test whether the same interface can describe rift and
transform resistance without forcing those processes into collision-specific rules.

## Design

### Model boundary and state

For the first stage, retain the current useful abstraction: a plate interior has one angular
velocity `omega` and rotates as a rigid body on the sphere. Boundary zones may deform, change
thickness and composition, transfer material, or create/remove plates. This RFC does not
propose a full mantle-convection solve, earthquake-cycle model, or continuum stress field
through every interior cell.

Represent a mechanically active boundary as a set of physical segments. A segment's shared
state should be limited to what its processes need, for example:

- geometry: endpoints or representative position, normal and tangent directions, and segment
  length in metres;
- both sides' material state: upper-crust, lower-crust, and mantle-lithosphere thickness,
  their compact composition and representative temperature, plus process-relevant strength
  or buoyancy derived from those properties;
- kinematics: relative normal closing/opening rate and relative tangential slip rate, in
  metres per year;
- ownership/provenance needed to apply the resulting deformation and book material.

The values are computed from the current world state and boundary geometry. They are not an
extra authoritative copy of per-node fields. Processes can add private state only when it
represents a real memory of the material (for example, inherited weakness), and must define
its lifecycle through remapping, splitting, merging, and old-save loading.

#### Minimal layered lithosphere state

The first physics-first material model should expose three persistent layers per represented
column: **upper crust, lower crust, and lithospheric mantle**. Their thicknesses define their
boundary depths; do not separately store boundary positions when they can be derived from the
thicknesses. Existing total crust and mantle-lithosphere thickness fields must remain
consistent with the layer sum or become derived values if this state is adopted. The
asthenosphere is the underlying mechanical/thermal medium, not a fourth plate layer.

Keep composition compact and specific to each layer rather than using one global rock-type
enum:

- Upper crust: sediment fraction, plus felsic fraction of the non-sedimentary material; mafic
  material is the remainder.
- Lower crust: felsic fraction, with mafic material as the remainder.
- Lithospheric mantle: a bounded depletion index from fertile to strongly melt-depleted.

These are proposed coarse mixture variables, not a claim that sedimentary, felsic, mafic, and
depleted are mutually exclusive rock classes. Fractions are mass fractions and must be
bounded in [0, 1], with transfer and mixing preserving material mass. A typed process boundary
may distinguish upper-crust, lower-crust, and mantle composition variants, but runtime
per-node state remains dense numerical arrays; do not allocate a Python object union for each
cell. Any variant-specific methods belong to the concrete type or shared process protocol,
not to the union alias itself.

Use a coarse thermal state rather than a fine vertical temperature mesh in the first
implementation. The provisional minimum is a crustal and a lithospheric-mantle thermal state;
derive upper- and lower-crust representative temperatures from a stated vertical profile.
Processes sample representative temperatures for strength and melt response. Thermal updates
may run less often than mechanical steps if timestep comparisons show that lag is acceptable.
Temperature is persistent only if the chosen thermal scheme evolves it; otherwise it is
derived from an explicit thermal-age or geotherm proxy. Water content, damage, and metamorphic
progress are not default fields: add them only when a process needs their history and cannot
derive the response from existing composition, temperature, and pressure state.

Melt generation and transport are process fluxes during a step, not automatically a permanent
fourth layer. Persist a melt reservoir only if melt survives between steps. Once magma
solidifies, book its volume and composition to the receiving crustal or mantle layer. An
optional intrusive body/layer is justified only where its depth and geometry affect later
mechanics or transport.

The existing `coal_deposit_m` and `oil_gas_deposit_m` fields remain surface resource
inventories, not buried sediment fractions, pore-fluid pressure, or lithology. They must not
directly weaken a layer. Sedimentary material can have a different baseline rheology, and
pore-fluid pressure can change effective stress and deformation response
[[Lockner & Stanchits, 2002]](#ref-lockner2002). Hydrocarbon abundance alone is not a proxy
for pore pressure or weakness. If sediment mechanics need this effect, add a coarse pressure
or fluid-state variable with its own generation, drainage, transfer, and timestep rules; do
not infer it from the resource-map deposits.

### Process response interface

Each mechanical process consumes the shared segment state and returns a response in physical
units. The initial contract should express:

1. **Loads:** force per unit boundary length (N/m) in normal and tangential directions, or
   distributed traction (Pa) for a basal/body force. Integrating over the segment or represented
   area gives force (N); applying its lever arm gives torque (N·m).
2. **Motion dependence:** how the load varies with relative velocity or plate angular
   velocities, including a linearized resistance when that is a valid approximation. A
   process may instead return a bounded nonlinear response or an admissible range.
3. **Capacity and yield:** the force or stress a segment can sustain, and what deformation
   takes up excess motion after yield. A geometric label such as “blocked” may identify a
   model condition, but it is not itself a constitutive law.
4. **Material response:** the rate or budget of thinning, thickening, subduction, accretion,
   slip, or transfer that follows from the solved motion. This is applied by the responsible
   deformation process and booked to existing or explicitly new ledgers.
5. **Diagnostics:** integrated load, work/dissipation where defined, yield state, and material
   budget, so the interface can be checked without exposing its internal discretization.

The interface is a contract, not a demand that all processes use the same constitutive law.
Subduction, collision, rifting, transform sliding, basal drag, and ridge push retain their own
driving physics and state. Shared segment geometry and units make their responses comparable
and composable.

### Timestep coupling

Coupling across successive timesteps is a first-class option. The step itself is an iteration:
motion changes relative boundary rates; deformation changes thickness, strength and geometry;
those changed fields alter the loads on the next step. A process does not need to be solved
simultaneously with every other process if this lagged feedback is accurate and stable at the
chosen timestep.

The provisional **staggered step** is:

1. Read a consistent pre-step state and construct/reuse boundary geometry and process
   classifications.
2. Ask each active process for its current load or force-response approximation on the
   segments it owns. For a response that depends on the unknown motion, evaluate it from the
   previous/current relative rate or use a bounded linearization.
3. Aggregate those loads and update plate angular motion. Compute all plates from the same
   pre-step snapshot, then apply their rotations together so plate iteration order does not
   determine which neighbor's motion is used.
4. From the resulting relative normal and tangential rates, let boundary processes calculate
   their deformation and material budgets. Apply those updates and close their ledgers.
5. Update topology (rift splits, collision merges, consumed plates) and rebuild only the
   interfaces invalidated by that change. The next timestep sees the changed geometry,
   material state, and resistance.

This is an operator-split approximation: a step may use a force or strength evaluated before
that step's deformation. It is acceptable when the lagged response converges to the intended
large-scale behavior as the step is refined, and when it does not skip a relevant threshold.
The process interface must state which values are sampled at step start, which are updated
after motion, and which feedback arrives on the next step.

Use a **local inner solve, substep, or fully coupled solve only where needed**: for example,
when a stiff drag causes explicit overshoot, when a yield threshold changes within one step,
or when lagged action/reaction produces an order-dependent or timestep-sensitive result. Do
not impose a global nonlinear solve on processes whose slower feedback is adequately captured
between steps. The choice is made by G2/G3 measurements, not by a blanket coupling rule.

### Efficient aggregation and approximation

The motion solve needs integrated force and torque, not every sub-grid stress value. Boundary
segments can aggregate node-scale information when the aggregation preserves the quantities
that govern the response:

- preserve total force and its moment arm for plate motion;
- preserve area/length-weighted material fractions and strength distributions when yield or
  threshold behavior depends on heterogeneity;
- retain blocked or yielding segment length and the relevant upper-tail load, rather than
  replacing a mixed contact with a mean that hides its strongest or weakest part;
- preserve physical material volume and provenance through deformation and remapping;
- report work/dissipation when a response law implies it, and reject approximations that add
  energy without an identified source.

Uniform linear drag may be aggregated by summing its response matrix. A thresholded or
nonlinear law may need several segment classes or an active-set solve; it must not be reduced
to a single average unless a comparison shows that the average preserves its outcome. Refine
only the segments whose force, yield, or deformation varies enough to change the answer.

Use real units at interfaces: metres, square metres, cubic metres, years, N, N/m, Pa, N·m,
and m/yr. A dimensionless approximation is allowed inside a process, but it must map to the
interface's physical response with stated bounds and a scale check.

### Timestep and resolution robustness

- Rates integrate over elapsed years, not per-step multipliers. A yielded deformation budget
  is rate × duration, subject to available material and capacity.
- A staggered/explicit response is acceptable if its error stays bounded as the timestep is
  refined. Processes with stiff drag, an intra-step yield transition, or another demonstrated
  instability use a local implicit/bounded solve or substeps; a global coupled solve is not
  assumed.
- Compare the same initial worlds at 2 and 5 Myr steps, and at node densities 1 and 4 where
  the interface aggregation is resolution-sensitive. Compare accumulated force response,
  relative motion, deformation, and ledgers, not just final maps. Record whether any mismatch
  comes from timestep lag, spatial aggregation, or the process law.
- Deterministic seed runs must remain reproducible. Any stochastic nucleation or weakening
  belongs to a physical event/time scale, not an incidental number of numerical substeps.

### Diagnosis gates

**Gate G1: process/interface inventory.** Before implementing a common interface, inventory
convergence, subduction, rifting, transform motion, basal drag, ridge push, merge/split, and
deformation. For each, record its current inputs, outputs, units, spatial support, state it
changes, feedback to plate motion, material owner, and cost. Identify which gaps actually
change an accepted geological metric. If a process already has the needed feedback, document
that and do not duplicate it.

**Gate G2: reduced coupled pilot.** Use one candidate contact process (RFC 0011's blocked
collision is the initial candidate) to compare its current implementation with a process
response interface on matched synthetic contacts and replayed worlds. Start with a staggered
update. Report force/torque response, local and mean relative motion, yield/deformation
allocation, material ledger closure, timestep sensitivity, and runtime. Compare with a local
inner solve on cases where the lagged update might be stiff. A pilot that only reproduces a
geometric event count does not pass.

**Gate G3: transferability.** Exercise the contract against at least one divergent/rifting and
one transform case before declaring it the general tectonic interface. If either requires
incompatible or meaningless shared fields, keep a narrower mechanical interface and document
process-specific extensions rather than forcing a universal API.

No full implementation plan is selected until G1–G3 are recorded. RFC 0011 can proceed as a
process-specific change if its own gates justify it; its implementation should not assume
that this RFC has already proven the shared interface.

### Cost and saves

The default cost is process responses over active boundary segments plus the current plate
motion solve. A local or global coupled solve adds work only for cases where the staggered
update fails its timestep or conservation checks. Work should scale with plates and active
segments, not with a fine global volume mesh. The pilot must measure its actual cost against
the current motion/deformation step at density 4. A general interface must add no more than
10% to median tectonic-step runtime on the acceptance worlds, unless a measured geological
improvement justifies a revised budget.

The interface itself should be derived/cacheable state and not serialized. Any new persistent
constitutive memory needs a default for old saves and tests covering its remap and topology
lifecycle. Existing surface material fields and ledgers remain authoritative.

## Detailed design (provisional touchpoints)

- `backend/app/torque.py`: current per-plate force gathering, torque/inertia calculation,
  implicit drag integration, and `shift_plate`; likely separate load assembly, motion update,
  and rigid rotation while preserving current force terms. A coupled solve is conditional on
  the gates, not a required replacement.
- `backend/app/world.py`: current step ordering and invalidation points; provisional location
  for the prepare → solve → deform → topology orchestration.
- `backend/app/quad_tectonics.py` and `backend/app/rheology.py`: boundary kinematics and
  material response for convergence, divergence, and transform motion.
- `backend/app/merge_split.py`: geometric split/merge triggers and topology changes; use the
  interface output rather than duplicating process force decisions.
- `backend/app/faults.py`: fault-local slip and deformation state; determine whether
  transform resistance belongs here or in a separate boundary-process module.
- `backend/app/tectonic_mechanics.py` (new, candidate): shared segment state and response
  protocols only if G1 shows they reduce duplication without erasing process-specific laws.
- Tests, if the gates select an implementation: equation-level force/torque cases; paired
  action/reaction; rigid-rotation limiting case; zero-load and zero-strength limits; timestep
  and resolution comparisons; material ledger closure; old-save loading for any new state.
- Docs: `docs/architecture.md`, `docs/simulation-model/foundations.md`, and
  `docs/simulation-model/tectonics.md` if the step flow or physical behavior changes.

## Assumptions and uncertainty

- **Rigid plate interiors remain useful at current scales.** Diffuse deformation exists and
  plate rigidity is an approximation [[Gordon, 1998]](#ref-gordon1998). This RFC localizes the
  first mechanical coupling to boundary zones; G3 can show when that boundary model is too
  narrow.
- **The correct process state is not yet known.** Crust type, thickness, keel strength,
  relief, and velocity are proxies already used by the model. G1 must show which are enough
  for the chosen metric and which missing variables (temperature, composition, damage,
  inherited fabric) materially alter the response.
- **The proposed three-layer state is a reduced representation.** Upper/lower crust and
  lithospheric mantle do not capture all stratigraphic or mantle compositional layering.
  Fractions, depletion, and coarse thermal states stand in for heterogeneous materials; their
  effect on force, yield, and melt response must be measured before adding more rock classes.
- **Sedimentary resources are not geomechanical state.** Surface coal and oil/gas inventories
  do not specify burial, lithification, pore pressure, or source-rock maturity. A future
  sediment weakening law needs an independently defined and conserved fluid/pressure state.
- **Force balance may be quasi-static or inertial at the chosen timestep.** Plate inertia is
  often negligible in geological torque balances [[Forsyth & Uyeda, 1975]](#ref-forsyth1975),
  while the current solver integrates angular acceleration. G1/G2 must distinguish a numerical
  stabilizer from a physical inertia term before choosing a general solve.
- **Averages can erase localization.** Segment aggregation is safe only when it preserves
  force/moment and the distribution needed by nonlinear yield laws. G2 and G3 explicitly test
  segment-level outcomes.
- **Coupling is not proof of physical validity.** A balanced solve can still use the wrong
  force law. Verification checks equations and invariants; validation requires geological
  comparisons, uncertainty ranges, and per-process evidence under design principle 9.
- **The project is a coupled Earth-like simulator, not a full geodynamic code.** This RFC
  prioritizes the mechanisms controlling plate-scale motion and material budgets. It does not
  claim to resolve grain-scale rheology, earthquake cycles, or mantle convection in 3-D.

## Acceptance criteria

- [ ] **G1** records the process inventory, physical quantities, ownership, feedback paths,
      and runtime for convergence, divergence/rifting, transform motion, and plate driving
      forces. Each identified gap is tied to a named design principle and a measurable
      geological or mechanical outcome.
- [ ] **G2** reports a staggered pilot on synthetic limiting cases and matched replay worlds;
      force and torque are dimensionally consistent, pair forces are equal and opposite,
      deformation budgets and material ledgers close to rounding, and the pilot stays within
      the runtime budget above.
- [ ] The pilot's accumulated response and key metrics agree within 10% between 2 and 5 Myr
      steps, except for explicitly quantized topology events whose event-time difference is
      reported and bounded. If the staggered update fails, a local solve/substep is measured
      before considering global coupling.
- [ ] **G3** reports whether the same interface can express a rift and a transform response;
      any process-specific extension is documented with a physical reason.
- [ ] Any pilot that consumes material properties identifies the upper-crust, lower-crust,
      and lithospheric-mantle state it reads; its composition fractions are bounded and its
      material transfers conserve the layer budgets. Temperature and any optional fluid state
      have an explicit update cadence and timestep-sensitivity result.
- [ ] Any implementation selected by the gates has equation-level tests, save compatibility
      for new persistent state, updated simulation docs, and measured runtime at density 4.
- [ ] No process is accepted solely because software tests pass; its force law and outcome
      are compared with a stated physical benchmark or a clearly labeled design hypothesis.

## Falsification criteria

- **Interfaces erase process physics.** If G3 requires collision, rift, and transform laws to
  share quantities that have no physical meaning for one of them, reject the universal
  contract and retain only shared geometry/unit conventions.
- **Aggregation changes threshold outcomes.** If the reduced segment model changes which
  contacts yield, where rifts localize, or which transform segments slip compared with a
  node-resolved reference beyond the acceptance tolerance, refine the sufficient statistics
  or abandon that aggregation.
- **Rigid interiors fail at the target scale.** If validation shows distributed intraplate
  deformation controls the target outcome, add a deformable plate-region model for that
  question rather than tuning boundary forces to compensate.
- **Coupled physics is too costly.** If a reduced solve cannot stay within the runtime budget
  without violating force balance, timestep robustness, or material conservation, narrow the
  coupled scope or choose a cheaper physically justified closure. Do not replace the missing
  response with an arbitrary geometric trigger.
- **The force law is not identifiable from available evidence.** If sensitivity runs show
  multiple plausible laws produce indistinguishable observables, leave the mechanism
  unresolved and report the uncertainty instead of choosing a single calibrated law.
- **The compact material state is insufficient.** If force, yield, or melt outcomes depend on
  vertical heterogeneity that the three layers and their coarse properties cannot preserve,
  add only the layer subdivision or material memory demonstrated to change those outcomes.

## References

- <a id="ref-forsyth1975"></a>Forsyth, D. W., & Uyeda, S. (1975). On the relative
  importance of the driving forces of plate motion. *Geophysical Journal of the Royal
  Astronomical Society*, 43(1), 163–200.
  [doi:10.1111/j.1365-246X.1975.tb00631.x](https://doi.org/10.1111/j.1365-246X.1975.tb00631.x)
- <a id="ref-gordon1998"></a>Gordon, R. G. (1998). The plate tectonic approximation:
  Plate nonrigidity, diffuse plate boundaries, and global plate reconstructions. *Annual
  Review of Earth and Planetary Sciences*, 26, 615–642.
  [doi:10.1146/annurev.earth.26.1.615](https://doi.org/10.1146/annurev.earth.26.1.615)
- <a id="ref-mulyukova2019"></a>Mulyukova, E., & Bercovici, D. (2019). The generation of
  plate tectonics from grains to global scales: A brief review. *Tectonics*, 38(12),
  4058–4076. [doi:10.1029/2018TC005447](https://doi.org/10.1029/2018TC005447)
- <a id="ref-lockner2002"></a>Lockner, D. A., & Stanchits, S. A. (2002). Undrained
  poroelastic response of sandstones to deviatoric stress change. *Journal of Geophysical
  Research: Solid Earth*, 107(B12), 2353.
  [doi:10.1029/2001JB001460](https://doi.org/10.1029/2001JB001460)
