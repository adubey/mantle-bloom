# Gate 3b geological-scale evaluation plan (opened 2026-10-30)

## Decision this plan implements

Gate 3b closure evidence must be produced on Mantle Bloom's actual spherical plate meshes over
geological time. The existing crack coupons remain component verification; they cannot be scaled up
by making a local mesh larger or running it longer. This is a proposed revision to RFC 0010's
accepted Gate 3 criteria. Gate 3a is closed as a local evidence phase. Gate 3b changes the scale of
the pass evidence and does not claim that a plate-scale rift closure is already implemented or
validated.

## Scale and run schedule

The first full-history comparison spans **200 Myr**, matching the existing broad Earth breakup
context. At the app's default `100,000 years/step`, that is 2,000 simulation steps per history.
Use whole-planet saved starting states and the actual plate mesh, not a rectangular cutout.

| Stage | World and duration | Step and resolution comparison | Purpose |
|---|---|---|---|
| A. Implementation screen | One controlled spherical world with a weak corridor, an intact control, a strong/cratonic band, and two rift segments with named step-over distances | Default mesh and step; repeat the decisive subset at node density 4 and half step | Verify nucleation ordering, corridor localization, arrest, linkage, shear transfer, ledger closure, and cost before a long history |
| B. Geological history | Matched starting states for seeds 7, 11, and 2024; old split rule as baseline and candidate rift closure as treatment; at least 200 Myr each | Production mesh and 100 kyr step | Measure rift-system formation, successful splits, daughter plates, craton intersections, budgets, and runtime |
| C. Robustness history | A representative seed selected before seeing its Stage B result | Node density 4 and 50 kyr step, same starting state and physical weak/strong fields | Check whether decisive event order and rates persist as the mesh and step change |
| D. Mantle/magma attribution | Matched candidate-closure histories with only one driver changed at a time, where the required field and budget exist | Same histories and event definitions as the control | Measure whether mantle/plume or magma forcing changes initiation, growth, linkage, and split outcomes; do not require every rift to be plume-driven |

Stages A and B must use the same production failure and segment rules. Stage A cannot substitute for
Stage B. Stage D is conditional on a separately budgeted mantle/magma source; an unavailable field
must be reported as a model limitation, not replaced with an unbalanced force.

## Required measurements

Track a rift-system identity from its first failed surface edge through abandonment or plate split.
Define one system as a connected, evolving set of rift edges; when two systems link, retain both
parent identities and record the link time. For every system, report:

- initiation time and position, host crust type, weak-belt/craton overlap, and connected length;
- whether it arrested, linked to another segment, or formed a successful plate split;
- duration, connected length history, daughter-plate areas, and any transform/transfer shear;
- source-by-source force and torque residuals, crust/mantle volume ledgers, and boundary work;
- runtime per step and added runtime relative to the matched baseline.

Report counts by 200 Myr for initiated rifts, abandoned rifts, linked systems, successful split
systems, and split events; do not combine these into a single “split count.” The broad estimate of
6–8 major continental breakup systems per 200 Myr is context only until event grouping and the
model's represented geology are compared consistently. Keep the existing matched-baseline targets
for split frequency, plate count and craton protection unless the owner separately revises them.

## Proposed acceptance tests and unresolved tolerances

| Question | Pass evidence | Tolerance status |
|---|---|---|
| Does the rift closure prefer weak lithosphere? | Paired weak-corridor and intact-control whole-sphere cases initiate in the weak corridor first; report initiation and split hazards by material class | Numeric minimum contrast needs owner approval before the run |
| Does strong lithosphere matter? | A stronger band or craton arrests or diverts growth more often than the matched weak path; report crossings and connected length | Preserve the accepted 50% craton-line/exposure reduction target against the matched old rule, evaluated over the agreed histories |
| Do segments link physically? | Two independent spherical rift systems remain distinct at large step-over, link at small step-over, and transfer shear at the transform without a plate-wide cut | Use physical step-over distances held fixed across meshes; the 150 km case is mandatory because the scratch probe did not converge there. Numeric convergence tolerance needs approval |
| Does topology change follow a viable rift? | Every new plate split has a connected boundary history and viable daughter plates; failed rifts do not split the plate | Exact daughter-size and boundary-connectivity thresholds need approval before testing |
| Is behavior stable with resolution and step size? | The same qualitative ordering holds at node density 4 and half step; event counts and growth measures are compared per 200 Myr | Numeric allowable drift needs approval before testing |
| Is the model geologically plausible? | Report split-system counts against the 6–8/200 Myr context; report failed rifts separately | Context, not a hard threshold |
| Is the solve usable and conservative? | Crust/mantle and force/torque budgets close; added runtime is measured | Keep existing 2% per-step overhead target unless owner revises it with measured trade-off |

The scratch flat-plate probe is not calibration evidence: it used scenario values, lacked a sphere,
thinning, heat and material accounting, and failed refinement on the 150 km step-over. Its useful
contribution is identifying step-over convergence and shear transfer as must-test risks.

## What is required before Stage A

1. Implement or select the candidate grid-scale rift closure on actual plate-surface adjacency,
   with a documented stress averaging length tied to lithosphere thickness and a strongest-first
   failure update. Keep fracture-tip physics explicitly sub-grid where it cannot be resolved.
2. Define how failed edges alter the spherical membrane equilibrium, plate-boundary kinematics,
   thinning and crust/mantle material ownership. Verify conservation on the sphere.
3. Freeze the weak-belt, strong-band, seed, step-over, driver, resolution and tolerance table before
   running the candidate as an evaluation.
4. Run Stage A as a feasibility/falsification screen. If step-over linkage or shear transfer fails
   refinement, report that result and revise or reject the closure; do not tune toward the desired
   outcome.

No Stage A or geological-history run is claimed in this plan. The immediate deliverable is the
scale-correct design and preregistered criteria; implementation and the actual whole-sphere runs
follow only after the unresolved acceptance tolerances are settled.
