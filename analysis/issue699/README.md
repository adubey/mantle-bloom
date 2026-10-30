# RFC 0010 Gate 3b: plate-scale rift closure

## Status

**Preregistration published before Gate 3b probes; method clarification appended 2026-10-30.** This record distinguishes the accepted
local-mechanism evidence in Gate 3a (#487) from the plate-scale closure now under test in
#699. No result in the Gate 3a coupons is counted as plate-scale validation.

The current whole-plate split violates design principle 1 (state the question and model
boundary) because flow-fit clustering selects an instantaneous plate-wide cut without
representing local lithospheric loading or failure. Its coin-flip breakup and whole-band
failed-rift thinning also do not express an accumulated local process (principle 3), and
craton-crossing splits conflict with principle 7. The observable invariant for this gate is:
**a plate remains one plate until locally failed surface adjacencies form a connected,
viable, kinematically divergent boundary; failure stays localized to the predeclared weak
or loaded region, and all material/topology changes remain attributable.**

## Candidate and scope

Use a cheap, grid-scale closure on production spherical plate meshes. Gate 3a supplies relative
local crack/weakness evidence, not a calibrated plate-scale law: seeded 3 MPa flaws initiated
under the tested 7 MPa host tension while the 14 MPa intact proxy stayed intact in 0/9 cases;
a moving tip advanced consistently in its controlled coupon and a 10x tougher band arrested it;
finite-source stress transfer increased tension ahead of the fixed seeded crack. Gate 3a did not
show spontaneous initiation or new damage ahead of a seed, and did not produce a spatial
plate-scale crack/weakness raster. Do not transfer these MPa or coupon lengths directly to plate
cells.

The candidate identifies local crack seeds and supported damaged adjacencies from the available
local weakness ranking and extension/loading signal. When segments need to link, a shortest-path
search over the production plate's exact cell-edge adjacency finds a low-cost route between
active tips. Physical edge length is the base cost; local resistance raises it, while inherited
weakness and locally supported extension lower it. Cratons/strong bands raise resistance. The
search only considers a bounded neighborhood of active tips and only advances an edge when the
local support rule permits damage. Thus the graph search proposes a cheap connection through
local cracks; it does not lay down a complete separator or create failure in unsupported cells.
If no supported route grows, the rift arrests. Split only after accumulated cuts make a
connected, viable, divergent boundary. Crack-tip behavior below cell scale stays an explicit
closure.

This is a local path-connection heuristic for plausible, affordable simulation, not new
crack-propagation science or a predictive geophysics claim. It differs from rejected #480
whole-plate static min-cut: endpoints arise from local segments, path search is local, edge
costs use local loading plus inherited weakness, and only locally supported edges are damaged.
The Gate 2 capacity anchors remain conditional scale references. No constitutive value is fitted
to historical split outcomes. The prior 0.6 weak-zone and 3.0 strong-zone multipliers remain
scenario values, not Earth calibration. Use production sparse-quad adjacency directly; the
territory-clipped Delaunay mesh is only a geometry diagnostic, not the rift graph.

The controlled screen uses generated cratons-first worlds for seeds 7, 11 and 2024, selecting
for each seed the largest connected continental plate by represented area before inspecting its
rift response. Mesh densities are 1, 4 (production), and 8; the cell spacing implied by the
model is 125, 62.5, and 44.19 km, respectively. Geometry, scenario masks, loading direction,
strength multipliers, and load schedule are matched across cases. An intact uniform-strength
case is the control for each seeded weak-corridor case. An actual craton mask and a separately
labelled 3.0x strong-band case test arrest/diversion; results from those two are not pooled.

If any selected plate is edge-disconnected, has under 95% coverage by a connected spherical
shell mesh, or cannot be mapped across resolutions without changing the prescribed physical
scenario, that case is recorded as ineligible with the reason. It is not silently replaced
based on its rift outcome. Whole-planet runs report all plates and the excluded share.

## Frozen predictions and thresholds

These are Gate 3b operational thresholds, frozen before probes. They do not revise RFC 0010's
existing historical acceptance targets.

| Question | Primary measure | Pass threshold |
|---|---|---|
| Weak-corridor localization | Applied-load multiplier at first failed adjacency, weak corridor / matched intact control | `<= 0.80` in each seed at density 4; at least 80% of failed edge length remains inside the prescribed corridor through `1.2 x` its initiation load, while the intact control has `<= 5%` failed edge length at that load |
| Strong/cratonic arrest | Whether a failed path crosses the predeclared 3.0x strong band; actual-craton response reported separately | No through-going path across the strong band by `2.0 x` the weak-path breakthrough load; at least 75% of tips touching the band arrest or turn along it rather than cross it |
| Segment linkage | Link load divided by the isolated single-segment propagation load at physical step-overs 0, 150, 300 and 600 km | 0 km links by `1.25x`; 150 km links by `1.50x`; 300 km stays unlinked through `1.50x`; 600 km stays unlinked through `2.00x`. Required ordering is monotonic with gap. The 150 km link-load change from density 4 to 8 must be `<= 20%`; density 1 is reported as a coarse-resolution ordering check, not required to resolve sub-cell geometry. |
| Transform/shear transfer | Sign and integrated magnitude of shear resultant in the predeclared transfer zone; topology; source-force and torque closure | Transfer-zone shear has the imposed sign, is nonzero and changes monotonically with imposed shear; no plate split occurs in the shear-only case. Force and torque residuals are each `<= 1%` of their applied source totals. This is a mechanics-transfer check, not a calibrated geological slip law. |
| Load/timestep partition | Failed-edge area intersection-over-union (IoU), initiation/link load | IoU `>= 0.80` and event-load drift `<= 10%` between one load increment and two half increments, and between 2 Myr and 1 Myr history steps where applicable. |
| Refinement | Per-case topology ordering and characteristic load thresholds | The qualitative weak/intact, arrest, and step-over ordering persists from density 1 to 4 to 8. Density 4 to 8 characteristic loads drift `<= 20%`; a threshold that does not converge is a failed gate result, not a tuning target. |
| Mechanical accounting | Elastic boundary work versus stored membrane energy plus explicitly booked failure dissipation; global force/torque residual | Work residual `<= 5%` of external boundary work; force and torque residuals `<= 1%`. If the candidate has no defensible failure-dissipation account, the work check fails and the closure must be revised before histories. |
| Material accounting | Continental-material ledger error; Hm source/sink accounts; daughter-plate material inventory | Continental ledger relative error `<= 1e-10` of initial continental inventory (the existing `continental_ledger.assert_closed` tolerance). Hm is not conserved: report each named source/sink and signed live change, with unattributed residual `<= 1%` of the absolute recorded gross change. No unbooked duplication or deletion at a split. |
| Runtime | Added rift-closure wall time per density-4 world step | `<= 2%` of matched baseline step time, measured over the same 20-step window after warm-up; report absolute seconds and per-plate scaling too. |

A miss on any predeclared local-mechanism or refinement threshold is reported as a miss. The
quantitative thresholds are operational screening cutoffs, not scientific constants; qualitative
ordering is the main claim. No
parameter sweep may be used to move the threshold after seeing the result. Any revised closure
or threshold requires a new dated preregistration and owner review before its acceptance run.

## Whole-planet comparison

After the controlled screen passes, compare the old split rule from `origin/main` with the
candidate closure on generated seeds 7, 11 and 2024. Generate matching initial worlds,
cratons-first layout, craton formation off, node density 4, 2 Myr steps, 200 Myr duration
(100 steps); repeat decisive comparisons at 1 Myr and density 8. Keep hash seed and software
versions fixed and record both source hashes. Run one world at a time. Report successful and
failed rift systems separately, plate count and daughter areas, split line-craton rate, band
craton exposure, craton loss, crust and Hm accounts, force/torque closure, boundary work,
runtime and event histories.

The 200 Myr history is an early scale screen; it cannot establish the RFC's 1 Gyr / multi-Gyr
craton-loss targets or a global breakup frequency from a small event count. Preserve the RFC
targets unchanged and label any underpowered historical comparison inconclusive.

## Reproduction inputs

- RFC branch base: `origin/main` commit `2baeac22d7c55a362063ce59b0770a398fe91d95`.
- Production schedule from the Gate 1 replay: 2 Myr; cratons-first layout; craton formation
  and post-PR-490 spherical relief disabled; seeds 7, 11, 2024.
- Controlled screen: generated worlds only, node density 1/4/8; no Downloads save is used.
  The matched seed-generated initial state is necessary for cross-arm equality and avoids
  selecting a saved trajectory based on its rift outcome.
- Local Gate 3a results are supporting mechanism evidence only. The 150 km flat-plate result
  failed its earlier refinement screen; this preregistration tests it anew on the sphere and
  treats failure to converge as falsification of this candidate.

## Method clarification (2026-10-30)

The issue's original preregistration comment described strongest-overload edge removal and
re-equilibration. The owner clarified that Gate 3b should use Gate 3a evidence to find a cheap
local path through cracks and weaknesses to form a rift, targeting a somewhat realistic simulation
rather than new science. Gate 3a's relevant evidence and limitations are summarized in
[PR #700](https://github.com/adubey/mantle-bloom/pull/700). This supersedes the earlier candidate
mechanics paragraph while preserving its scenarios and operational screening thresholds.

Before the first path diagnostic, the exact rule was frozen in [issue comment
6520271866](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520271866): on exact
production cell adjacency, use great-circle edge length times mean endpoint relative resistance
divided by positive opening projection as the Dijkstra cost; an edge is supported only when the
prescribed dimensionless load times its opening projection reaches its resistance. Grow one
adjacent edge per update and recompute. The first screen uses synthetic prescribed corridors and
a prescribed extension axis: Gate 3a supplies no plate-scale weakness raster. Dijkstra spans the
full exact plate graph between local endpoints (it is not clipped to a tip neighborhood), but
returns one route rather than a plate separator. This checks graph geometry, not a production
stress field or plate-scale validation. The geometry-only Delaunay mesh pre-screen is not rift
evidence. See [scope correction](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520392557).

## Results

### Geometry eligibility screen

The first geometry pre-screen selected the largest continental plate by represented area and
checked exact source-graph connectivity plus a territory-clipped spherical triangulation used
only to inspect shell coverage. Six of nine planned seed/resolution cases were eligible with
one exact source component and at least 99.99% area coverage. Seed 7 failed the deterministic
20-plate generation partition at densities 1, 4 and 8; these cases remain ineligible rather than
being replaced. Full per-case measurements are in `plate-mesh-screen.json`. The Delaunay shell
coverage is not used as the path graph.

### First local graph path screen

The corrected run used seed 11, density 4, a generated continental plate, endpoints selected
before reading material scenarios at 600 km great-circle separation, exact sparse-quad adjacency,
and the preregistered prescribed axis with `lambda=1.5`.

| Scenario | Result | Route length | Weak corridor share | Interpretation |
|---|---|---:|---:|---|
| Weak corridor (150 km wide) | Supported route, 13 nodes | 748.5 km | 100% | A low-cost connected route was found. |
| Intact control | No supported route | — | — | Negative at this one synthetic loading setting. |
| Weak corridor plus strong band (300 km wide) | No supported route | — | — | Arrest, not diversion, in this setup. |

This is a graph-algorithm screen only. Material masks and loading are prescribed, no production
local stress field is coupled, and no crack growth, topology change or plate split occurs. It
cannot establish the RFC's plate-scale or geological acceptance criteria. The initial attempt's
topology-diameter endpoints produced a roughly 20,255 km route; that endpoint-selection result is
retained separately in `local-path-screen-global-endpoints-ineligible.json` and excluded from
interpretation. The corrected result is in `local-path-screen.json`, and the result/limitations
were posted in issue comment
[6520330787](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520330787).

### Robustness subset

The same frozen graph rules were run at seed 11 densities 1/4/8 and seed 2024 density 4.
The actual-craton arm uses the existing `cratons.strength(craton_crust_m)` as a separate
resistance floor, `max(scenario resistance, 1 + 2 * strength)`; the factor of 2 is an uncalibrated
screening choice.

| Case | Weak only | Intact | Strong band | Actual-craton floor |
|---|---|---|---|---|
| Seed 11, density 1 | Connected; 805 km; 100% weak-mask path; 100% craton-node overlap | No path | No path | No path |
| Seed 11, density 4 | Connected; 749 km; 100% weak-mask path; 0% craton-node overlap | No path | No path | No path |
| Seed 11, density 8 | Connected; 795 km; 100% weak-mask path; 72% craton-node overlap | No path | No path | No path |
| Seed 2024, density 4 | No path | No path | No path | No path |

The weak/intact fork and arrest ordering repeat across seed 11's three mesh densities, but the
weak-only route's actual craton overlap changes sharply with resolution. Adding the existing
craton-strength floor arrests that route at each density. Seed 2024 has no supported route at the
fixed load in any arm. These are recorded as limitations/misses, not validation: the projection
and support rule are too restrictive or seed-sensitive for a general path mechanism, and a weak
mask can exploit cratons if the resistance map omits them. No tuning follows from these outcomes.
All robustness cases are synthetic graph diagnostics; no production stress, damage evolution,
segment linkage, split or geological history has been measured. Per-case results are in
`local-path-screen-seed11-density1-craton.json`,
`local-path-screen-seed11-density4-craton.json`,
`local-path-screen-seed11-density8-craton.json` and
`local-path-screen-seed2024-density4-craton.json`; the results were posted in issue comment
[6520377961](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520377961).

### Current design decision

Do not treat the prescribed-axis graph screen as a production rift mechanism. In production,
`torque.gather_per_node_force_field` supplies per-node source forces but its own API describes
them as inputs to an internal-load calculation, not as local stress. `thin_shell.solve_force_field`
can produce a reduced membrane diagnostic, but that solver is explicitly offline/diagnostic and
is not coupled to production damage or topology updates. RFC 0010 likewise calls local loading
unresolved: per-node forces currently reduce to whole-plate torque, and neither torque nor
`Plate.internal_stress` determines intraplate stress.

Before a production path rule, Gate 3b needs a local extension-demand field on the production
cell graph (or a documented shell solve), a force-transfer/accounting rule, and source-force,
torque and work checks. Then the Gate 3a weak/intact/tip-arrest ordering can inform actual edge
damage and linkage. Until those pieces exist, the synthetic graph results are feasibility
screens only; they neither pass nor falsify plate-scale rifting. The blocker and source references
were posted in issue comment
[6520406763](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520406763).


### Local-load transfer feasibility

The preregistered offline shell diagnostic did not produce a usable local demand field. At seed
11, density 4, the solve was interrupted after approximately 120 seconds while SciPy LSMR
was still iterating; no JSON result, demand quantiles, or force/torque closure measurements
were produced. This is a censored runtime outcome, not a convergence failure.

The predeclared reduced-resolution follow-up kept the same seed, plate-selection rule, force
assembly, material thickness, modulus and clipped-shell method at density 1. It failed while
assembling the mesh with `ValueError: thin-shell mesh contains a degenerate triangle`, before
the equilibrium solve; no output JSON was produced. This records a mesh-eligibility failure,
not a solver timing result. The preregistration and both outcomes are in issue comments
[6520456345](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520456345) and
[6520542634](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520542634).

Together with the solver's diagnostic-only status, these outcomes do not support putting the
current FEM solve in a production per-step route. They also do not justify relabeling per-node
source forces or plate torque as local stress. Gate 3b therefore remains at the phenomenological
path-screen stage: a candidate production load transfer would need its own explicitly approximate
force-to-edge rule and conservation/accounting checks, or a different validated low-cost solver.
No local damage or production rift update has been implemented from this diagnostic.

### Coarse shell-response follow-up

The preregistered coarse follow-up reduced the same seed-11, density-4 plate to deterministic
512- and 1,024-medoid meshes, aggregated existing per-node forces and thickness, and solved the
offline linear membrane diagnostic. The clipped meshes contained isolated triangle slivers, so
the method amendment retained each mesh's largest edge-connected triangle component and reported
the discarded components. The result is compactly recorded in
`coarse-shell-response-seed11-density4.json`.

Both meshes covered over 99% of represented plate area. Coarsening changed net force by 1.16%
and 0.30% of summed source-force magnitudes at 512 and 1,024 representatives, respectively;
torque residuals were 1.07% and 0.40% of summed source-torque magnitudes. Solver equilibrium
residuals were 1.02% and 2.88%, both above the preregistered 1% cutoff. The two spatial rankings
nonetheless agreed at Spearman 0.911, with 0.726 area IoU for the top demand decile. Each solve
took under one second. This supports feasibility for a coarse diagnostic field but misses the
predeclared force-transfer cutoff at 512 and equilibrium cutoff at both resolutions. The plotted
field is therefore labeled provisional; it is not a production extension, stress, damage or rift
field. No numerical parameters were tuned after seeing the miss. The [preregistration](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6527325898),
[mesh amendment](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6527356236) and
[result](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6527479690) are recorded
on issue #699.

The updated plate plot colors the coarse principal tensile-resultant response behind the craton
mask and crack-path candidates. Color values are clipped at the plotted mesh's 99th percentile to
show spatial contrast, then nearest-medoid mapped onto the plate cells for display. The shell
field is from generation time, while the route paths use the state after one 2 Myr tectonic step;
the overlay is therefore a spatial illustration, not a matched-state coupled diagnostic. The
figure is generated by `plot_gate3b_plate_paths.py --shell-json` and does not add validation
beyond the checks above.


### Data-backed divergent-boundary linkage screen

The preregistered initial-state case had no divergent nodes on the selected plate. The amended
rule used exactly one 2 Myr tectonic step after generation (climate/biomes disabled) before
selecting the largest connected continental plate. Seed 11 at density 4 then had 213 divergent
nodes in 38 exact-adjacency components. The closest qualifying pair of distinct components was
612 km apart. The fixed 500–700 km selection rule was applied without reading resistance fields.

The route used exact plate adjacency and great-circle edge length times mean endpoint resistance.
Inherited non-cratonic belts used relative resistance 0.6; intact cells used 1.0. The craton arm
floored resistance at `1 + 2 * cratons.rift_resistance(craton_crust_m)`. These values are routing
ranks, not calibrated strength or support thresholds.

| Case | Divergent nodes / components | Endpoint gap | Route length, belt arm | Belt fraction | Craton fraction |
|---|---:|---:|---:|---:|---:|
| Seed 11, density 1 | 404 / 41 | 614 km | 953 km | 0% | 0% |
| Seed 11, density 4 | 213 / 38 | 612 km | 612 km | 50% | 58% |
| Seed 11, density 8 | 130 / 29 | 604 km | 971 km | 0% | 0% |
| Seed 2024, density 4 | 0 / 0 | ineligible | — | — | — |

At density 4, the inherited-belt route cost was 500,151 relative metres versus 612,372 for the
intact control, an 18% reduction with the same geometric route. The craton floor increased that
route's cost to 1,110,101 but did not divert or arrest it; 58% of the route length crossed craton
cells. At densities 1 and 8, the selected route did not encounter belt or craton cells, so all
three resistance arms returned the same route. Seed 11 route length was not resolution-stable,
and seed 2024 supplied no eligible divergent pair after its fixed single step. These outcomes do
not establish robust weak-path preference or craton protection. No tuning or case substitution
followed.

The exact per-case results are `boundary-linking-seed11-density1-poststep.json`,
`boundary-linking-seed11-density4-poststep.json`,
`boundary-linking-seed11-density8-poststep.json`, and
`boundary-linking-seed2024-density4-poststep.json`. The runs and scope are summarized in issue
comments [6520975606](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520975606)
and [6520999485](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6520999485).

This advances the graph screen from prescribed crack endpoints to actual divergent boundary
segments and model belt/craton fields. It still only proposes a shortest connection: it does not
check local force support along interior edges, grow damage over time, split the plate, or test
conservation and work. The seed and resolution misses mean this candidate is not ready for a
production rift update.


### Craton barrier sensitivity

The first craton-protected route arm was still only a soft path cost: `rift_resistance` is at
most 0.75, so the chosen `1 + 2 * rift_resistance` formula capped a cell at 2.5 times
intact resistance. Dijkstra has no route-support cutoff and returns a path whenever the graph is
connected, so the route could still traverse the craton. This is why the original plotted route
looked nearly straight and crossed the shield; it was not a viable split boundary.

A preregistered sensitivity then removed all edges incident to craton cells. With the original
endpoint pair, the source itself was cratonic, so that arm had no path from its seed. A second
preregistered variant selected the same nearest-pair rule from divergent nodes outside craton
cells (endpoint gap 598 km). The soft actual-belt route was 801 km and crossed craton cells for
77% of its length. The finite craton penalty shifted the route slightly and reduced craton overlap
to 38%, but still crossed the shield. When craton-touching edges were prohibited, no craton-free
path connected these two off-craton divergent segments. Thus the current graph geometry offers no
ordinary around-craton connection for this endpoint pair; it does not justify pushing a rift
through the craton.

Per-case outputs are `boundary-linking-seed11-density4-craton-barrier.json` and
`boundary-linking-seed11-density4-noncraton-tips.json`. The preregistration and result are in
issue comments [6521203014](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6521203014)
and [6521224977](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6521224977). A craton-edge veto is a qualitative ordinary-rift sensitivity, not a claim of absolute physical impermeability.

### Geometry correction: boundary linkage is not an intraplate rift

The data-backed linkage screen above selected divergent components from the existing plate
boundary. Its endpoints and much of its route follow the plate edge, so it tests a connection
between boundary faults; it is not an intraplate rift, does not demonstrate plate breakup, and is
excluded from Gate 3b rift evidence. The craton-resistance and no-path observations remain valid
only for that boundary-linkage geometry.

Gate 3b needs to seed finite cracks inside continental plate material, using local weakness and
extension evidence, then propagate between interior crack tips on supported adjacencies. A route
that joins existing edge faults, or merely minimizes distance between edge segments, cannot meet
that criterion. The next diagnostic must explicitly keep its nucleation and growth corridor inside
the plate and check that a completed cut leaves two material-retaining sides. Until a local
interior loading signal exists, synthetic interior paths remain graph-feasibility screens rather
than a physical rift prediction.

### Interior weak-belt candidate: preregistered first trial (2026-10-30)

The plate-edge linkage screen is excluded from rift evidence. Before running the next candidate,
freeze an interior-only graph trial on seed 11, density 4: generate the same cratons-first world,
advance one 2 Myr tectonic step with climate/biomes disabled, then select the largest continental
plate by represented area. Define interior as nodes outside the production boundary-force band
(`dist_to_neighbor > boundary_force_reach_rad`, or no neighbor), and remove craton nodes from
candidate seed segments. Find exact-adjacency connected components in the remaining inherited
`cratons.belt_mask`; select the closest pair of distinct components whose nearest nodes are 500–700
km apart, nearest 600 km, with node-index tie breaks. If no such pair exists, record ineligible;
do not substitute endpoints after observing routes.

Use those component nodes as candidate interior segment tips. The tip-to-tip great-circle direction
is an illustrative prescribed opening axis, not a derived local stress field. For each edge, compute
the prior frozen tangent-plane opening projection `D_e=max(0, |n_e dot x_hat|)` and mean endpoint
relative resistance (belt 0.6, intact 1.0, synthetic strong band 3.0). Route cost is great-circle
edge length times mean resistance divided by positive opening projection; an edge is supported only
when `lambda * D_e >= R_e`, with the frozen illustrative `lambda=1.5`. Remove every edge incident
to an actual craton node in the ordinary-rift arm: no supported path means arrest, not craton
crossing. Advance at most one supported edge from an active tip per update, then recompute.

This one-trial candidate checks whether real inherited weak-belt geometry can connect interior
segments while respecting the craton barrier. It does not validate load transfer, damage physics,
segment birth, split viability, step-over refinement, or Gate 3 acceptance. The prescribed axis and
load multiplier remain graph falsification inputs only. The candidate must be replaced or stopped
if it cannot produce an interior supported route without crossing craton; no parameter tuning or
endpoint replacement follows this trial.

The route graph is further bounded to nodes within twice the selected tip separation of the
spherical midpoint between the tips. This local-search radius is fixed before reading the route;
it prevents the candidate from becoming a plate-wide separator search.

#### Implementation correction before interpreting the first result

The first implementation mistakenly required an all-supported path before beginning growth. The
frozen rule instead lets Dijkstra propose a low-cost route over positive-opening links, checks
support only on the next edge from either active tip, advances one supported edge, and recomputes.
The first output, saved as `interior-belt-path-seed11-density4-all-supported-precheck.json`, is
retained as a support-connectivity precheck, not as the candidate routing result. The code is corrected to implement the frozen
stepwise rule; the same seed, tip-selection rule, local radius, resistance values, load multiplier,
and craton veto are retained for the candidate rerun.

#### First interior weak-belt candidate result

The corrected stepwise trial selected two distinct interior inherited-belt components on seed 11,
density 4, with tips 598.8 km apart and a fixed 1,197.6 km local-search radius. The plate had
33,856 cells, 6,218 interior non-cratonic belt cells, and 37 eligible belt components. At the
illustrative `lambda=1.5`, the inherited-belt arm and matched intact arm both advanced 14 supported
links (12 from the source tip, 2 from the target tip), then arrested. Their paths were identical;
the remaining tip gap had grown to 1,526.5 km. The supported-edge share within the local search
neighborhood was 49.93% in both arms. The craton-edge veto was active throughout, so neither arm
crossed craton.

This is not evidence of a useful weak-path advantage: the real inherited-belt map did not change
the selected growth sequence relative to intact material, and the tips diverged instead of linking.
It does show the stepwise support rule can arrest a pair without a forced plate split. This single
roughly-600-km step-over trial is not a pass/fail of the full step-over/refinement gate and does not
provide local stress, crack nucleation, force transfer, viable breakup or geological validation.
The initial all-supported-connectivity precheck is preserved separately as
`interior-belt-path-seed11-density4-all-supported-precheck.json`; the corrected candidate output is
`interior-belt-path-seed11-density4.json`.

#### Next frozen candidate case: 150 km interior step-over (before run)

Use the same generated-world schedule, seed 11, density 4, plate-selection rule, actual inherited
belt/craton fields, local boundary exclusion, opening projection, resistance values, illustrative
`lambda=1.5`, stepwise tip growth, and craton-edge veto as the 600 km trial. This is a separate
predeclared case, not replacement of the 600 km pair. Select the closest pair of distinct eligible
interior belt components whose nearest nodes lie 100–200 km apart, nearest 150 km with node-index
tie breaks. No eligible pair means record ineligible without widening the interval. Keep the
local search radius at twice the selected gap. A linked path by `lambda=1.5` is the preregistered
short-step prediction; otherwise record arrest. Compare the same pair against the intact control.

#### 150 km step-over candidate result

The fixed selection found a 150.35 km gap between eligible interior inherited-belt components 10
and 11. At `lambda=1.5`, both the inherited-belt and intact controls advanced eight supported
links (five from the source and three from the target) before arrest. Their edge sequences were
identical, and both stopped with the tips 67.57 km apart. Each path bent along the local graph,
with cumulative growth lengths of 334.13 km and 200.13 km; neither crossed craton. The local
supported-edge shares were 54.17% for inherited belts and 52.08% for intact material.

The predeclared 150 km prediction (link by `lambda=1.5`) was missed. The small supported-edge
share difference did not alter propagation or connect the segments. Together with the 600 km
case, this first candidate does arrest without crossing cratons but has not demonstrated weak-path
localization or the required short/long step-over ordering. Do not tune the support multiplier or
replace endpoints from these results; the missing production local-loading signal and candidate
linkage rule need revision before any full gate run. Output: `interior-belt-path-seed11-density4-stepover150.json`.

Reproduction commands:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_gate3b_interior_belt_paths.py \
  --out analysis/issue699/interior-belt-path-seed11-density4.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_gate3b_interior_belt_paths.py \
  --target-gap-km 150 --minimum-gap-km 100 --maximum-gap-km 200 \
  --out analysis/issue699/interior-belt-path-seed11-density4-stepover150.json
```
