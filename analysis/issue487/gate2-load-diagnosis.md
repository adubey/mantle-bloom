# RFC 0010 diagnosis gate 2: local loading available to rift damage

**Date:** 2026-10-19
**Scope:** Read-only force/motion trace plus a diagnostic source-force adapter and reduced thin-shell equilibrium solver. Production motion, deformation, topology, and save behavior are unchanged.

## Result

The source-force interface preserves each existing torque contribution as a per-node vector in N. `thin_shell.solve_force_field` transfers the combined tangent load through a linear triangular membrane equilibrium solve and returns local in-plane membrane stress increments in N/m. It removes the area-weighted rigid-rotation load handled by plate motion; applied radial force is assigned to shell-normal support. The comparison is defined as `N_ext = max(0, lambda_max(Delta_N))`, the largest positive principal value of the incremental membrane-resultant tensor, divided by `N_cap = integral(Delta_sigma_yield(z) dz)`, the column's additional extensional differential strength to normal-fault yield from its assumed `K0` background. Both are resultants in N/m. The interpretation assumes vertical compression remains the major principal stress, an elastic stress increment about a zero-load reference, and the same representative column geometry at each node. This maps units and failure mode; it does not add prestress to the shell or feed local yielding back into equilibrium. Synthetic tests verify equilibrium residual, rigid-rotation removal, linear force scaling and the named capacity pairing. Diagonal stiffness scaling lets the iterative solve handle generated plate meshes with thousands of nodes.

This establishes an offline load-transfer and strength-comparison interface, not a production stress
field. The solver assumes linear plane-stress elasticity, free in-plane displacement, prescribed
nodal loads, no faults, bending or plasticity, and caller-supplied effective thickness and modulus.
For an isolated single-plate solve, free displacement is the baseline because the inputs contain no
displacement anchor to an external reference. Edge clamps remain artificial sensitivity probes. A shared-node two-plate manufactured case now
verifies interface transfer for uniform isotropic tension, while event-level paired remeshes measure
sensitivity to adding a force-contributing neighbor. These probes do not specify a geologic contact
law or calibrate actual interplate displacement and traction continuity. Gate 2 covers source-connected candidates; disconnected-source events
are counted separately and excluded from its coverage denominator.

The strength envelope has conditional anchors and an uncertainty bracket, but no spatially
calibrated field. The reference uses depth-varying freshwater hydrostatic pressure
`P_f = rho_water * g * z`; fixed `P_f/P_lith=0.42` is only an overpressure comparison. Friction is
bracketed by `mu=0.6-0.85` from Byerlee's rock-friction experiments. `K0`, cohesion, geotherm, creep
law and local strain rate remain assumptions. Six unit tests check the hydrostatic pressure
profile, Mohr-Coulomb normal-fault increment and dry/wet power-law inversion against analytic
values. The shell's curved-sphere manufactured solution checks load-to-resultant mapping for its
stated idealized case, not the real force field or interplate traction transfer.

The exact 22,896-node matched-load screen was re-scored without another world replay. The hydrostatic
reference classifies 1,563 nodes over capacity (6.83%), compared with 1,718 (7.50%) for fixed
`P_f/P_lith=0.42`. This is threshold sensitivity, not an observed rupture rate. Broader load and
boundary validation and locally constrained capacity inputs remain open.

## Synthetic solver validation

A shallow spherical cap spanning ±0.05 rad in each tangent-plane coordinate was loaded by
equivalent edge tractions for a uniform uniaxial membrane resultant of `1e10 N/m`. The thin-shell
solver used `E = 70 GPa`, `H = 100 km`, and a uniform supplied capacity of `2e10 N/m`.

| Node grid | Surface triangles | Mean recovered tensile resultant | Relative error | Equilibrium residual / load norm |
|---:|---:|---:|---:|---:|
| 3 × 3 | 8 | `9.984535e9 N/m` | −0.1547% | `9.42e-13` |
| 5 × 5 | 32 | `9.984318e9 N/m` | −0.1568% | `9.18e-9` |
| 9 × 9 | 128 | `9.984202e9 N/m` | −0.1580% | `2.01e-7` |

Across these refinements the recovered mean varies by 0.0034% of the analytic load. This checks
the manufactured linear membrane case and confirms that the partial-cap mesh excludes the convex
hull's artificial underside. This grid series by itself does not establish accuracy for general
plate outlines, other load distributions, alternative edge conditions or a natural lithosphere
strength profile.

An irregular convex cap with 67 nodes and 111 surface triangles was then loaded with the same
manufactured uniaxial edge traction. With free edges its mean resultant was `9.988472e9 N/m`
(0.1153% low). Clamping the six nodes on the western edge changed the mean to
`9.706457e9 N/m` (2.82% below the free-edge result); both solves had relative free-node residuals
below `2.3e-6`, and the clamped displacements were exactly zero. This confirms the solver responds
to the selected boundary condition on an irregular convex patch. The clamped case has no separate
closed-form benchmark, so this is a sensitivity check, not an accuracy claim. That check did not
cover nonconvex outlines, other edge placements or actual plate meshes.

The solver now accepts explicit triangle connectivity, allowing the mesh to preserve a concave
outline rather than filling it with the point-set convex hull. An L-shaped spherical patch was
loaded with the same manufactured uniform uniaxial membrane resultant:

| Grid | Nodes | Triangles | Free-edge mean resultant | Relative error | Residual / load norm |
|---:|---:|---:|---:|---:|---:|
| 9 × 9 | 65 | 96 | `9.984677e9 N/m` | −0.1532% | `1.96e-6` |
| 17 × 17 | 225 | 384 | `9.984629e9 N/m` | −0.1537% | `1.67e-6` |

The two free-edge means differ by 0.00048% of the analytic load. On the 17 × 17 mesh, clamping
the western edge gave `1.002079e10 N/m` mean resultant and clamping the southern edge gave
`8.374740e9 N/m`; the shift between those constrained cases was 16.5%. Relative free-node
residuals were `1.76e-6` and `1.11e-6`, and constrained displacements were exactly zero. These
clamped results have no separate closed-form benchmark, so they show sensitivity to constraint
placement, not accuracy or a preferred physical boundary condition. The concave case validates
explicit-mesh load transfer for one manufactured load; broader loads and physically justified
edge conditions remain open.

## Incremental extensional stress on a curved shell

The shell stores an incremental membrane-resultant tensor `Delta_N` in N/m, with tension positive
and zero load as its reference. Its scalar extension measure is
`N_ext = max(0, lambda_max(Delta_N))`. The column diagnostic returns
`N_cap = integral(min(Delta_sigma_b, Delta_sigma_d) dz)`, also in N/m, where each branch is
additional horizontal differential stress to failure from the stated background ratio `K0`.
Therefore `N_ext / N_cap` compares incremental extensional load against remaining normal-faulting
capacity. The interpretation requires `0 < K0 <= 1`, vertical compression as the major principal
stress, and a shared depth/area convention. Shell loads are not total prestress, and this ratio is
not valid for thrust-fault regimes, large rotations, or yielded material that redistributes load.
The new result properties are `max_incremental_extensional_resultant_n_per_m`,
`extensional_capacity_resultant_n_per_m` and `extensional_failure_ratio`; tensile-named attributes
remain compatibility aliases.

An independent curved-sphere solution checks a broad patch under uniform isotropic tension `N`
with its analytic conormal boundary traction. Radial pressure `2N/R` balances spherical curvature
and is supported by the shell-normal pressure term omitted from the in-plane solve. This gives
`lambda_max(Delta_N) = N` everywhere without using the solver's plane-stress element equations to
construct the reference. The reproducible mesh sequence is:

| Theta × phi nodes | Triangles | Mean error | p95 absolute error | Maximum absolute error | Residual / load |
|---:|---:|---:|---:|---:|---:|
| 5 × 9 | 64 | 2.94% | 6.17% | 12.97% | `4.20e-7` |
| 9 × 17 | 256 | 1.46% | 3.23% | 4.24% | `3.05e-6` |
| 17 × 33 | 1,024 | 1.18% | 2.75% | 2.91% | `8.65e-6` |

The mean and p95 errors fall monotonically with refinement. The fine mesh's p95 error is below
3%, while the maximum local error remains about 2.9%. The equilibrium residual rises but remains
below `1e-5` of load. This validates the stated extensional invariant on a curved shell for one
smooth analytic stress state; it does not validate arbitrary high-gradient production loads.
Reproduce it as part of `bin/debug/validate_thin_shell_plate_meshes.py`.

## Edge-condition decision

Use free in-plane displacement with the supplied loads as the single-plate diagnostic baseline.
The force field contains no prescribed displacement against a fixed mantle frame, so clamping a
selected side invents an external anchor. On the two original broad generated meshes, west/east/
south/north clamps, opposite-side clamp pairs and full-perimeter clamps were compared with free
edges:

| Generated plate | West+east clamp | South+north clamp | Full perimeter clamp |
|---|---:|---:|---:|
| #6, oceanic with hole | −38.53% | −39.55% | −100.00% |
| #8, continental, single loop | −46.57% | −14.49% | −100.00% |

Each value is the change in mean extensional resultant relative to the free-edge solve. Single-side
clamps also produce shifts from −0.58% to −23.52%, depending on location and outline. These checks
show that arbitrary constraints can dominate the inferred load. The free-edge condition is selected
for the current isolated-plate diagnostic because it adds no unsupported displacement anchor; it is
still a simplified isolated-body assumption, not a claim about natural plate boundaries. Resolving
tractions and displacement continuity across real plate boundaries require a coupled mesh. The
new shared-node benchmark and two event probes exercise that mesh numerically, but leave contact
force ownership, fault slip and the interface law unspecified; those physical aspects remain open.

## Representative generated plate meshes

`bin/debug/validate_thin_shell_plate_meshes.py` reproduces the probe below. It generates seed 7
with 12 plates and node density 0.5, selects the largest hole-bearing and largest one-loop plates,
then Delaunay-triangulates their cell centres in a tangent projection. Eight interior probes per
candidate triangle are checked with exact plate containment. The accepted connectivity is retained,
while node positions and cell areas are scaled to a shallow cap for a controlled manufactured
uniaxial traction. This tests representative generated node distribution and boundary topology
without claiming a uniform stress field on a broad spherical plate.

| Generated plate | Nodes | Territory loops | Hull triangles with outside centroids | Explicit triangles | Recovered resultant | Error | Residual / load |
|---|---:|---:|---:|---:|---:|---:|---:|
| #6, oceanic (largest with hole) | 2,390 | 2 | 118 / 4,757 (2.48%) | 4,548 | `9.991120e9 N/m` | −0.0888% | `8.29e-5` |
| #8, continental (largest single-loop) | 1,988 | 1 | 94 / 3,956 (2.38%) | 3,749 | `9.991270e9 N/m` | −0.0873% | `2.19e-5` |

The explicit meshes retained the generated two-loop and one-loop boundaries. They remove unused
outside points and rejected elements, so their active node counts are 2,389 and 1,987. This is a
positive topology and solver-conditioning check on two representative generated meshes. It does
not establish the fidelity of eight-probe clipping for arbitrarily narrow features, nor stress
accuracy on the original curved plates, real production forces or alternative edge conditions.
Compact measurements are saved in
[`representative-mesh-summary.json`](representative-mesh-summary.json).

### Original curved mesh and edge sensitivity

The same explicit meshes were solved on their original spherical coordinates and cell areas.
The load is a manufactured uniaxial edge-traction pattern defined in each plate's gnomonic
tangent projection, then projected into the local spherical tangent plane at each edge. It is a
controlled numerical probe, not a production torque field or an analytic stress benchmark. The
projected spans are 3.39 × 2.28 rad for plate #6 and 2.13 × 2.33 rad for plate #8, so these are
broad curved patches. Each clamp fixes only boundary nodes in a band 12% of the larger projected
span wide at the named west/east/south/north extreme.

| Plate | Free mean (`1e10 N/m`) | West change | East change | South change | North change | Residual / load across cases |
|---|---:|---:|---:|---:|---:|---:|
| #6, oceanic with hole | 1.3625 | −1.47% (17 nodes) | −0.58% (11) | −7.62% (35) | −23.52% (59) | `5.01e-5`–`8.50e-5` |
| #8, continental, single loop | 1.3213 | −1.20% (25 nodes) | −1.54% (37) | −8.86% (26) | −7.20% (32) | `2.13e-5`–`5.27e-5` |

The solver converges to a small free-node residual on both original curved meshes. Constraint
placement changes the mean and upper-tail stress distribution, with the direction of the largest
change depending on plate shape. These numbers show conditioning and boundary sensitivity; because
this broad-plate traction has no closed-form curved-shell solution, the resultant's difference
from `1e10 N/m` is not an accuracy error. The free edge is selected for the isolated diagnostic;
the clamp cases quantify unsupported displacement anchors, while physical interplate coupling
remains unmodeled.

### Explicit phenomenological capacity alternative

The available state cannot support a physical depth-resolved capacity. As a testable design
alternative, the diagnostic now evaluates a directly integrated threshold in N/m:

`T_cap(i) = T_ref * (1 - lambda * w_i)`

`T_ref` is a tunable reference capacity in N/m. `w_i` is a caller-specified weak-zone weight in
`[0, 1]`; `lambda` is the maximum fractional weakening in `[0, 1)`. This is deliberately
phenomenological: it assumes no depth profile and does not derive strength from the existing
Moho-temperature or Mohr-Coulomb proxies. The probe uses `lambda = 0.5` and a synthetic band with
`w_i = 1` when the projected x-coordinate is within 8% of the plate's x-span around its median,
and 0 elsewhere. It sweeps `T_ref = 1e10`, `2e10`, and `4e10 N/m` against the free-edge curved
mesh resultant.

At `T_ref = 2e10 N/m`, the weak-band and background failure fractions are 98.1% versus 10.1% for
plate #6, and 90.0% versus 2.5% for plate #8. Across the capacity sweep, failed-node counts
decrease monotonically (2,387 -> 892 -> 0 for #6; 1,935 -> 587 -> 17 for #8). This validates the
rule's units, threshold ordering and intended localization under the manufactured load. The band,
thresholds and weakening are test inputs, not calibrated parameters or evidence that a natural rift
will follow a weak belt. The alternative needs explicit acceptance and replay validation before
it can replace the physical capacity requirement.

## Source trace and units

| Path | Available quantity | Reduction / limitation for local rifting |
|---|---|---|
| `torque.gather_boundary_force_inputs` | Per-own-node position, crust/mantle thickness, nearest-neighbor direction and distance, nearest neighbor plate/node, and neighbor angular velocity. | Geometry/material sampling only. It has no stress tensor, adjacency force, or load-transfer operator. Nearest-neighbor direction is not an intraplate tensile axis. |
| `torque.gather_per_node_force_field` / `PerNodeForceField` | Per-node slab-pull, slab-drag, ridge-push, basal-drag and contact-drag force vectors (N); total external force and area-average applied traction (Pa). | Preserves the current spatial external forcing before torque aggregation. The vector sum reproduces the current torque terms, but `total_traction_pa` is not internal tensile stress. |
| `thin_shell.solve_force_field` / `solve_thin_shell` | Incremental membrane-resultant tensor (N/m), largest positive principal increment, supplied additional normal-fault capacity (N/m), and extensional failure ratio. | Solves reduced linear elastic equilibrium on outward-facing convex-hull facets by default, or caller-supplied triangles for explicit outlines. A broad curved sphere with uniform isotropic tension recovers the analytic principal resultant with 1.18% mean error and 2.75% p95 absolute error on a 17 × 33 grid. Free in-plane displacement is the single-plate diagnostic baseline; side, paired and perimeter clamps are sensitivity probes because no displacement anchor is present in the inputs. Radial load is excluded and rigid rotation is removed. Interplate coupling and production constitutive feedback remain outside this diagnostic. |
| `slab_pull_force_n` | A scalar force in N on selected subducting boundary nodes: positive column buoyancy (Pa) times sunk length and represented boundary length. | A boundary-local driving force; not interior traction. `slab_pull_torque` crosses it with the planet-radius lever arm and sums it. |
| `ridge_push_force_vectors_n` / `ridge_push_torque` | Boundary-local force magnitude from isostatic relief, density contrast, gravity, and represented boundary length (N). | Direction is the local away-from-neighbor tangent; it is now retained as a vector for load transfer before moment aggregation. |
| `basal_drag_force_vectors_n` / `basal_drag_coefficients` | Node-area-integrated basal drag. The underlying force per area has traction units (Pa), using mantle-flow and rigid-plate velocities. | Distributed basal force is retained node-wise and transferred by the diagnostic shell solve; it is not itself in-plane stress. |
| `collision_drag_coefficients` / `contact_drag_pa_s_per_m` | Boundary resistance based on relative plate velocity; coefficient in Pa·s/m times velocity gives Pa traction. | Applies to classified interplate collision/jammed contacts. It does not describe intraplate rift loading or damage. |
| `predict_plate_motion` / `integrate_omega` | Sum of explicit torques plus affine drag terms, integrated against whole-plate inertia to update one `omega` vector. | This is a plate-scale motion solve; it returns no spatial stress or strain-rate field. Rigid rotation alone has no internal strain. |
| `boundary_closing_rate_m_per_s` / `rheology.plastic_strain_rate_per_myr` | Boundary-normal relative plate velocity (m/s), mapped by an effective viscosity proxy to scalar Pa, then a calibrated plastic response. | Defined only for boundary-classified nodes. It is not a 2-D stress tensor or vertically integrated force, and cannot be reused as an intraplate rift field without a new physical rule. |
| `Plate.internal_stress` / `accumulate_plate_stress` | Dimensionless, time-integrated plate radius/overlap accumulator with decay. | Plate-level split eligibility relaxation, not a measured force/stress and not spatially locating damage. |

### Existing thermal and strength proxies

The current state provides ingredients for a future capacity model, but not an integrated tensile
capacity field:

| Source | Existing behavior | Why it does not yet define capacity |
|---|---|---|
| `orogeny.steady_moho_temperature_c(hc_m, hm_m)` | Computes steady Moho temperature from a one-dimensional conductive proxy with crustal heat production and fixed surface/base temperatures. | Gives one temperature at the Moho, not a calibrated temperature profile for every depth and material layer. |
| `orogeny.moho_temperature_c(..., lag_c)` and `moho_thermal_lag_c` | Applies a scalar per-column thermal lag, updated when tectonic thickness changes and relaxed over time. | The lag adjusts the Moho estimate; it does not evolve a depth-resolved geotherm. |
| `rheology.yield_excess(sigma_n_pa)` | Uses cohesion `2e7 Pa`, friction angle 30°, and a scalar effective-viscosity boundary-normal stress proxy. | Represents boundary-local Mohr-Coulomb yielding and deformation; it is not a tensile criterion or depth-integrated column strength. |
| Per-node `Hc`, `Hm`, `crust_type_code`, craton/restite and thermal-lag fields | Preserve thickness, broad crust identity and a few tectonic/material tracers. | Do not specify vertical composition, temperature, water/pore pressure or strain-rate-dependent brittle and ductile strength. |

A physical extensional capacity requires a stated failure criterion and a depth-integrated vertical
yield envelope, in N/m. It is differential strength against normal-faulting extension, not intact
rock tensile strength. The existing Moho-temperature and boundary-yield proxies cannot be
multiplied by column thickness and labeled physical capacity without new constitutive assumptions
and calibration. A phenomenological cohesion-times-thickness capacity would be a separate
alternative and should be labeled and tested as such.


## Candidate physical capacity: extensional strength envelope

The initial physical hypothesis should be a column strength envelope, not an intact-rock
tensile-strength threshold. Under extension, shallow fractured lithosphere is commonly limited
by frictional sliding and normal faulting; deeper hot lithosphere is limited by thermally
activated creep. Classic strength-envelope work combines those mechanisms and emphasizes their
dependence on effective pressure, temperature, strain rate and rock type
[[Brace & Kohlstedt, 1980]](https://doi.org/10.1029/JB085iB11p06248)
[[Kohlstedt et al., 1995]](https://doi.org/10.1029/95JB01460). This is a candidate first-order
column model, not a claim that two homogeneous layers capture natural rifts.

For each node, define depth `z` downward from the surface and evaluate the following inputs:

1. **Temperature:** construct a full crust-plus-mantle geotherm from surface temperature, crustal
   heat production and lithosphere-base temperature. The existing steady-state proxy documents a
   parabolic crustal profile and linear sub-Moho profile, but the cross-section helper currently
   holds temperature constant below the Moho and the stored thermal lag is only one Moho
   temperature offset. A first capacity run should use the documented steady-state profile and
   bracket it against a clearly labeled lagged-geotherm approximation. The latter needs an
   explicit depth interpolation rule; copying the Moho offset unchanged through the column is not
   a resolved transient heat solution.
2. **Lithostatic and pore pressure:** compute `P_lith(z) = integral_0^z rho(s) g ds` from the
   existing crust/lithosphere densities. Set `P_f(z) = lambda_f P_lith(z)` and use effective
   compression `sigma' = sigma - P_f`. Use dry rock (`P_f=0`) and a connected freshwater hydrostatic reference (`P_f=rho_water*g*z`) as distinct cases. The fixed-ratio tests (`P_f=lambda_f*P_lith`, including `lambda_f=0.42`) represent imposed fluid pressure fractions; 0.42 is not hydrostatic for the model's rock densities. None is inferred local fluid state. The world has no pore-pressure field; overpressure, drainage and fluid migration are
   omitted. Pore pressure materially changes brittle strength, so a single dry or hydrostatic
   choice could misstate rift resistance [[Brace & Kohlstedt, 1980]](https://doi.org/10.1029/JB085iB11p06248).
3. **Brittle extension:** use a Mohr-Coulomb normal-faulting limit with compressive stress positive.
   For vertical maximum compression `sigma_1 = sigma_v'` and horizontal minimum compression
   `sigma_3 = sigma_h'`, failure occurs at

   `sigma_h_fail'(z) = (sigma_v'(z) * (1 - sin(phi)) - 2 * C * cos(phi)) / (1 + sin(phi))`.

   Given initial horizontal compression `sigma_h0' = K0 * sigma_v'`, additional extension to
   failure is `Delta_sigma_b(z) = max(0, sigma_h0' - sigma_h_fail')`. With neutral baseline
   `K0 = 1`, this reduces to

   `Delta_sigma_b(z) = 2 * (sin(phi) * sigma_v'(z) + C * cos(phi)) / (1 + sin(phi))`.

   Here `sigma_v' = P_lith - P_f`, `C` is cohesion and `phi` is the friction angle. This assumes
   vertical stress is the maximum principal compression and a horizontal stress is the minimum;
   another stress orientation can change the failure mode. The existing `C = 20 MPa` is an inherited boundary-deformation proxy, not a measured intact-rock or fault cohesion. Bracket friction using coefficients `mu=tan(phi)=0.6–0.85` from Byerlee's laboratory synthesis; this is a frictional strength range, not a locally measured fault value. Cohesion and K0 remain design assumptions. A tensile cutoff for intact, unfractured rock and fault-orientation anisotropy are not
   represented in this first pass.
4. **Ductile extension:** for each layer and chosen rock scenario, invert a laboratory power-law
   creep relation at a stated strain rate, for example
   `epsilon_dot = A * sigma^n * d^(-m) * f_H2O^r * exp(-(E + P V)/(R T))`.
   The required stress is `sigma_d(z)`; define its remaining capacity as
   `Delta_sigma_d = max(0, sigma_d - abs(sigma_v' - sigma_h0'))`, assuming the same principal
   orientation for the background state. Laboratory flow laws distinguish crustal mineralogy
   from olivine mantle and show substantial wet/dry and parameter uncertainty
   [[Kohlstedt et al., 1995]](https://doi.org/10.1029/95JB01460)
   [[Hirth & Kohlstedt, 2003]](https://doi.org/10.1029/138GM06). The current state has no layerwise
   composition, grain size, water content or local strain-rate field, so those inputs must be
   scenario parameters. The first comparison should report at least weak/wet and strong/dry
   brackets; it must not select one silently.
5. **Column integration:** at each depth take the lower *additional* stress required for brittle
   failure and ductile flow, `Delta_sigma_y(z) = min(Delta_sigma_b(z), Delta_sigma_d(z))`, and integrate
   `N_y = integral_0^(Hc+Hm) Delta_sigma_y(z) dz`, yielding N/m. This is the classic integrated
   extensional differential-strength scale. It can be compared with the positive maximum
   principal increment of the in-plane shell resultant only as an *additional* extension measured
   from a zero-load reference, with the column's background horizontal stress specified by `K0`
   and vertical compression remaining the major principal stress. This restriction excludes
   thrust-fault regimes and does not include prestress or yielded load redistribution.

**Diagnostic result.** The standalone implementation is
[`bin/debug/diagnose_gate2_column_envelope.py`](../../bin/debug/diagnose_gate2_column_envelope.py);
its complete assumptions, reference profiles, parameter sources and 1,728 scenario results are in
[`gate2-column-envelope.json`](gate2-column-envelope.json). Using a 0.125 km depth interval, the
illustrative hydrostatic-reference capacities are 25.91 TN/m for the 35 km crust + 100 km mantle-lid
continental column, 0.516 TN/m for the 11 + 8 km young-ridge oceanic column, and 6.20 TN/m for
the 11 + 40 km mature-reference oceanic column. These are outputs of assumed material states, not
Earth-calibrated strength estimates. Across the explicit normal-faulting scenario sweeps, the
unweighted p10/median/p90 capacities are 6.17/13.98/29.62 TN/m (continental),
0.149/0.310/0.570 TN/m (young ridge), and 2.17/4.09/7.29 TN/m (mature oceanic). These percentiles
summarize a designed grid, not probability distributions. At 0.25 km spacing, reference capacities
shift by 0.12%, 0.92%, and 0.01%, respectively. In the continental sweep, median capacity changes
from 8.25 TN/m at `K0=0.6` to 23.14 TN/m at `K0=1.0`; wet versus dry mantle changes it from 10.73
to 17.98 TN/m, and the fixed pore-pressure-ratio cases change it from 11.84 TN/m at 0.42 to 16.71 TN/m when dry. Missing
background stress and material state therefore dominate the interpretation.

The diagnostic reconstructs a steady 1-D geotherm from current model constants; it does not use
the stored Moho thermal lag as a transient depth profile. It applies one quartz-rich upper crust and
one effective lower-crust law to continents, one effective mafic law to oceanic crust, dry/wet
olivine endmembers to the mantle, a depth-varying freshwater hydrostatic reference pressure, and a
normal-faulting Mohr-Coulomb branch with `K0` in [0.6, 1]. The 1,728-case factorial still samples
fixed `P_f/P_lith` fractions of 0 and 0.42 to test fluid-pressure sensitivity; 0.42 is not called
hydrostatic. The laws are laboratory/geodynamic parameterizations extrapolated to
geological strain rates. The mafic-granulite law is the central fit to reconstituted, nominally dry
granulite with 0.16–0.28 wt% dissolved water; its reported fit uncertainties are not sampled, so
this is not a water-free endmember. Semibrittle flow, localization, inherited faults, anisotropy,
evolving fluids, compositional mixing, thermal transients and 3-D stress transfer are omitted. In
particular, the dry-olivine creep inversion at the continental Moho exceeds its experimental regime
by orders of magnitude; the minimum brittle/ductile branch controls the reported envelope, but that
does not validate the extrapolation. See the script and JSON for full assumptions and per-case
outputs.

**Recommendation.** Keep the shell as a load-transfer diagnostic and compare
`max(0, lambda_max(Delta_N))` against the depth-integrated additional extensional capacity from the
column envelope. The shell exposes that pairing with explicit names and compatibility aliases.
The column diagnostic reports geotherm, pore-pressure, rock-law, strain-rate and `K0` scenarios,
branch crossover and integrated-capacity ranges. The independent curved-sphere benchmark passes the
stated refinement criteria, and free in-plane displacement is selected as the isolated-plate
diagnostic boundary condition; clamps are sensitivity probes only. A conditional scale-calibration
audit now brackets representative continental columns against independent strength-scale literature;
see [`gate2-capacity-calibration.md`](gate2-capacity-calibration.md). Do not wire this strength into
production or calibrate damage from legacy candidate cuts, which are not rupture labels. The
per-node `Hc`, `Hm`, crust type, craton share, belt history and Moho thermal lag are available,
and the matched replay now uses those fields for a thickness-aware diagnostic capacity map. The
steady and interpolated-lag screens differ by only eight over-capacity nodes across 22,896 samples,
but remain conditional on assumed temperature, rheology, fluid pressure and background stress.
Depth-resolved temperature, composition, fluid, strain-rate and fabric state remain missing. Keep
the map diagnostic and complete broader load/boundary validation before using a capacity threshold
in production; see [`gate2-capacity-state-audit.md`](gate2-capacity-state-audit.md).

**Known realism limits.** The current columns collapse a vertically layered, anisotropic, fluid-
and rate-dependent lithosphere into two thicknesses and broad crust identity. The temperature
proxy omits transient heat transport and the diagnostic thermal view does not resolve sub-Moho
temperature. The density profile omits topography-dependent lateral pressure, dynamic pressure and
fluid pressure. A single friction angle/cohesion omits inherited faults, fabric, damage and frictional
rate/state behavior. Power-law creep extrapolated from laboratory samples to geological conditions
has uncertain water, grain-size and mineral-mixture effects. A one-dimensional column ignores
along-strike segmentation, flexure, bending, faults and 3-D stress transfer. The reduced shell is
linear elastic with a uniform modulus/thickness and currently uses selected edge constraints; it
does not feed yielded regions back into load transfer. Each is a deliberate first-model
simplification that can make predicted rift thresholds and locations unrealistic; sensitivity
brackets and replay behavior can expose dependence on these choices but cannot turn missing state
into measured geology.

The complete trace is in `backend/app/torque.py` (`gather_boundary_force_inputs`, force terms, `gather_per_node_force_field`, `predict_plate_motion`, `integrate_omega`), `backend/app/thin_shell.py` (`solve_force_field`, `solve_thin_shell`), `backend/app/rheology.py` (`normal_closing_rate_m_per_s`, `plastic_strain_rate_per_myr`), and `backend/app/merge_split.py` (`accumulate_plate_stress`, `maybe_split_plate`). `test_spatial_force_field_preserves_each_existing_torque_term` verifies the vector adapter's moment; `test_thin_shell.py` checks synthetic equilibrium, rigid-mode removal and linearity. No production damage path consumes this diagnostic.

## Mesh-connectivity safeguard

The thin-shell solver requires triangle elements in one plate solve to form one connected mesh.
A disconnected mesh has additional rigid modes; the plate-wide rigid-rotation projection does not
remove those modes, so the least-squares residual can reflect mesh topology as well as load transfer.
For a source-connected plate with a mesh-only split, the replay first tries a spherical-hull mesh
clipped by centroid and edge-midpoint probes, then retains the dominant triangle component only if
at least 99% of triangles remain. It rechecks shared-edge connectivity and requires at least 95%
coverage of both the candidate cut band and its craton subset. A source-disconnected plate remains
excluded whole; no separate-island force boundary condition is supported.

A seed 2024 replay through 460 Myr compared exact edge-sharing `plate.adjacency()` with the
accepted triangle mesh for the two late plate 9 events. Both source plates have seven cell
components, so there is no single continuous source shell to repair into. The mappings are:

| Time | Exact source component sizes | Mesh components and source mapping | Result |
|---:|---|---|---|
| 432 Myr | `1848, 20, 19, 14, 1, 1, 1` | `1844→1848`, `19→20`, `19→19`, `12→14`; the three singleton cells have no triangles | Mesh keeps the four substantial source regions separate, but leaves 8 cells out |
| 460 Myr | `1760, 17, 17, 14, 1, 1, 1` | `1773→(1760+17)`, `17→17`, `12→14`; three singleton cells have no triangles | Delaunay clipping also bridges two source regions that exact edge adjacency says are separate |

These are the historical 55.52% and 51.47% residual cases. At 432 Myr, the residual is not a
least-squares fitting bug: the source shell consists of separate regions with no modeled material
connection. At 460 Myr, the mesh additionally introduces a false connection between source
components. Bridging those regions would invent a load path; solving each island independently
would require an unsupported force/traction boundary condition. The replay excludes any event whose exact source adjacency has multiple components. Those 56
source-disconnected candidates are not “mesh failures”; the plate cells themselves have no single
edge-connected material domain for this solver. The old partial replay numbers in the paragraph
above (25 solves and six residual failures) predate the final shared-edge/mesh repair treatment and
are historical only. The full regenerated replay matches all Gate 1 production cuts and ordered
candidate identities for seeds 7, 11 and 2024. It solves all 74 source-connected candidates and
repairs the ten former edge-disconnected meshes (six by spherical-hull outline probes; four by
retaining a ≥99% dominant triangle component). All 74 pass stop code 1 and residual/load ≤1%
(maximum 0.4514%) on 22,896 cut-band samples. Minimum mesh coverage is 99.67% for all cut-band
nodes and 99.44% for craton-band nodes.

The Gate 2 population is explicitly limited to candidates whose source cells form one exact
edge-connected component. Require valid solves for at least 90% of this eligible population overall
and 80% in each seed, event kind, nonempty time bin and nonempty event-kind × time-bin stratum. This
is a project evidence threshold, not a physical constant. Current connected-source coverage is
74/74 = 100% overall and in every nonempty stratum, so this completeness requirement passes. The
other 56/130 candidates (43.1%) have disconnected source-cell components and are excluded from the
Gate 2 denominator because the isolated shell solve cannot define their component load paths. Their
frequency is uneven: connected-source events cover 33/69 splits (47.8%) versus 41/61 failed-rift
events (67.2%); by time the rates over all events are 51/75 (68.0%), 14/34 (41.2%) and 9/21
(42.9%). Keep those exclusions visible and do not generalize the conditional capacity scores to
disconnected-source events or the full event population. Compact per-seed results and source/mesh
audit are in [`gate2-matched-replay.json`](gate2-matched-replay.json) and
[`gate2-connectivity-audit.json`](gate2-connectivity-audit.json). Gate 2 remains open for a defensible spatial source/capacity mapping and broader load/boundary
validation. Per-node thickness and proxy fields now feed the diagnostic capacity screen; its
scenario assumptions and results are recorded in [`gate2-capacity-state-audit.md`](gate2-capacity-state-audit.md)
and [`gate2-spatial-capacity-replay.json`](gate2-spatial-capacity-replay.json). These results do
not resolve the missing local thermal, material and stress state, and do not close Gate 2. A
one-factor sensitivity sweep on exact matched node loads finds the strongest capacity dependence
in `K0`, mantle water/rheology and pore pressure. Relative node rankings remain stable, but the
binary exceedance share ranges from 3.47% to 11.23% across the tested `K0` endpoints. The tested
envelope is not joint and does not identify these fields. Results and limitations are in
[`gate2-capacity-sensitivity.md`](gate2-capacity-sensitivity.md) and
[`gate2-capacity-sensitivity.json`](gate2-capacity-sensitivity.json).

Reproduce the partial source-filtered replay with:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
  bin/debug/replay_gate2_capacity.py --seed 2024 --years 4.6e8 \
  --out /tmp/issue487-g2-seed2024-source-filter.json
```

## Coupled-boundary load validation

The manufactured benchmark in [`benchmark_gate2_coupled_shell.py`](../../bin/debug/benchmark_gate2_coupled_shell.py)
uses a mildly curved spherical cap split into two equal-material regions. The isolated runs put
analytic uniform isotropic membrane traction on all edges, including equal-and-opposite traction
on the cut. The joined run shares interface nodes and loads only the outer boundary. The analytic
shell-normal Laplace pressure `p=2N/R` is included in the global force balance audit but remains
outside the in-plane shell solve.

| Check | Result | Acceptance limit |
|---|---:|---:|
| Global force balance including analytic normal-pressure reaction / boundary load | `9.44e-6` | `< 1e-3` |
| Boundary torque / (boundary load × radius) | `2.92e-17` | `< 1e-8` |
| Joined in-plane equilibrium residual / load | `2.48e-7` | `< 1e-5` |
| Maximum interior interface reaction mismatch / target resultant | `2.19e-7` | `< 1e-5` |
| Maximum joined-versus-isolated regional mean `N_ext` difference | `0.0094%` | `< 3%` |

The mesh has 119 nodes, 192 triangles and seven shared interface nodes; two interface corners are
excluded from the local reaction comparison because they also meet external edges. All three solves
stop with code 1. The joined and isolated regional mean resultants recover `1e10 N/m` within
`0.0290%` and `0.0196%`, respectively. The compact result is
[`gate2-coupled-shell-benchmark.json`](gate2-coupled-shell-benchmark.json). This verifies the
linear shared-node membrane implementation on its manufactured case; it does not validate a
natural plate contact, fault slip, damage, or a strength law.

The initial event comparison has been superseded after an audit found that its nodewise relative-change metric divided by each isolated-node N_ext, including zeros. The quoted 26.2% and 51.2% relative-change numbers were invalid; the event mean N_ext values were not affected. The corrected diagnostic now compares one connected-source split and one failed-rift per seed, checks exact event keys against the Gate 1 sidecars, measures mutual pair-specific contact-force and torque closure, and reports stress differences normalized by the mean isolated cut-band N_ext. Results and the precise force ownership assumptions are in gate2-interface-sensitivity.md and gate2-interface-sensitivity.json.

The six short event solves all achieve 100% source cut-band coverage and joined residual/load below 0.14%. The corrected mean absolute cut-band stress difference from the isolated result ranges from 5.3% to 97.0%. Four sampled boundaries have nonzero mutual contact drag; its pair force closure residual is 2.0–21.4% and torque closure is 4.6–21.6%. In the other two events the selected pair has no contact-drag force, so those cases do not test reciprocal traction. Removing mutual contact drag in the tied-interface variant changes cut-band stress by 0–71.8% of the current-contact joined mean, without a consistent direction.

These remain short-run sensitivity probes, not physical coupled predictions. The tied mesh is a perfectly bonded/no-slip idealization, while the torque engine's contact drag is a velocity-dependent force rather than a finite-slip shell interface law. A six-event pair-level rigid-rotation audit finds 3.0–95.1% mean absolute cut-band stress sensitivity between the joint three-mode projection and a separate-group projection. The separate-group variant retains bonded joint displacements, so it is not a valid free-slip alternative. Both variants converge; the sensitivity reflects unresolved interface assumptions, not numerical failure. See [gate2-rigid-projection-audit.md](gate2-rigid-projection-audit.md). Gate 2 remains open until a mechanically consistent interface and load-transfer formulation is validated, alongside source/capacity justification.

## Gate disposition

- **Current diagnostic interface:** `max(0, lambda_max(Delta_N))`, the largest positive principal shell stress increment, is paired with the column's additional normal-fault capacity from an assumed `K0` background; both are in N/m. The physical envelope has conditional reference anchors and a broad scenario bracket, but no
spatial calibration. A broad curved-sphere analytic benchmark converges to 1.18% mean error and 2.75% p95 error at 17 × 33 nodes. Free in-plane displacement is selected for the isolated-plate diagnostic because the input has no external displacement anchor. Single, paired and full-boundary clamps show substantial sensitivity and are not adopted as physical boundary conditions. The direct-threshold phenomenological rule has a monotonic capacity and weak-band localization check, but is a competing design only. Within the narrowed connected-source scope, the matched comparison passes the 90% overall / 80% stratum coverage gate (74/74); 56 disconnected-source candidates are excluded and reported separately, so results do not generalize to them or the full population. A local strength field remains underidentified. A shared-node manufactured benchmark passes. A corrected six-event split/failed-rift sample has 5.3–97.0% band-scale N_ext sensitivity to adding a neighbor; mutual contact forces are not consistently reciprocal, and the bonded interface is not a validated contact/slip law. Earlier 26.2%/51.2% nodewise-relative summaries were invalid due to zero denominators and are superseded. Gate 2 is not closed.
- **`internal_stress`:** retain only as a plate-scale eligibility/background signal if desired; it cannot trigger or weight local damage.
- **No calibrated strength field:** do not assign `Y_t`, `eta_rift`, or failure-strain parameters from these diagnostics. Do not infer stress from applied traction, rigid rotation, flow-cluster residuals or `Plate.internal_stress`.
- **Required next work:** justify or constrain local `K0`, pore pressure and mantle water/rheology independently of legacy event outcomes, then validate load and boundary sensitivity within the connected-source scope. Keep disconnected-source event counts and strata visible, and do not generalize findings to those events. The 90% overall / 80% stratum coverage rule applies to connected-source candidates. A defensible multi-component load-transfer model would be a separate scope expansion, not a prerequisite for the current Gate 2 population. Gate 3 accumulation and segment-growth calibration remains downstream. The six-event audit has identified mutual contact drag as internal to a selected pair, measured force/torque imbalance, and compared current pair contact loads with a tied-interface case that omits only that mutual term. Gate 2 is scoped to the single connected source plate and its free-edge shell diagnostic; coupled neighbor displacement and a finite-slip law are outside its closure criteria. The six-event pair audit is sensitivity evidence only: contact drag has measurable force/torque imbalance, and projection choices produce materially different stress under different kinematic assumptions. Do not present those pair results as physical coupled stresses. A future coupled-stress model is a separate scope expansion requiring consistent interface kinematics, rigid-mode treatment and traction ownership.

The RFC now records this outcome and keeps the local-cut damage design conditional. Gate 1 is closed by the explicit baseline/exclusion decision recorded in `analysis/issue487/README.md`; the archived-stream provenance limitation remains documented.
