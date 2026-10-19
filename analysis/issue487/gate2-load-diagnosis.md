# RFC 0010 diagnosis gate 2: local loading available to rift damage

**Date:** 2026-10-19
**Scope:** Read-only force/motion trace plus a diagnostic source-force adapter and reduced thin-shell equilibrium solver. Production motion, deformation, topology, and save behavior are unchanged.

## Result

The source-force interface preserves each existing torque contribution as a per-node vector in N. `thin_shell.solve_force_field` transfers the combined tangent load through a linear triangular membrane equilibrium solve and returns local tensile stress resultants in N/m. It removes the area-weighted rigid-rotation load handled by plate motion; applied radial force is assigned to shell-normal support. The result includes a ratio to caller-supplied depth-integrated tensile capacity in N/m. Synthetic tests verify equilibrium residual, rigid-rotation removal and linear force scaling. Diagonal stiffness scaling lets the iterative solve handle generated plate meshes with thousands of nodes.

This establishes a diagnostic load-transfer interface and validates generated plate topology in a controlled manufactured case, not a production stress field. The solver assumes outward convex-hull triangular elements by default, linear plane-stress elasticity, free edges and no faults, bending, plasticity or thermal strength envelope. A manufactured uniaxial load on one shallow spherical cap recovers its analytic resultant across three refinements; representative generated cell-centre meshes also recover it after topology-preserving mapping to a shallow cap. Original broad curved plates, physical boundary conditions, production load distributions and the strength provider remain unvalidated. Gate 2 therefore remains open; Gate 3 calibration must wait for those checks or a deliberate revision to a clearly labeled phenomenological rift model.

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
from `1e10 N/m` is not an accuracy error. The free edge remains the diagnostic default, not a
physically selected boundary condition.

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
| `thin_shell.solve_force_field` / `solve_thin_shell` | Nodal membrane stress resultants (N/m), maximum tensile principal resultant, supplied integrated strength (N/m), and failure ratio. | Solves reduced linear elastic equilibrium on outward-facing convex-hull facets by default, or caller-supplied triangles for explicit outlines. Shallow convex and L-shaped concave manufactured uniaxial loads are recovered within 0.16% over refinements. On the concave mesh, west versus south edge clamps shift the mean resultant by 16.5% with small free-node residuals; these are sensitivity checks, not accuracy benchmarks. Radial load is excluded and rigid rotation is removed. Representative plate meshes, physically justified boundary conditions, constitutive and omitted-process sensitivity, and a strength provider remain open. |
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

A physical tensile capacity would require a stated failure criterion and a calibrated vertical yield
envelope, then integration through the modeled column, `T_c = integral(sigma_yield(z) dz)`, in
N/m. The existing Moho-temperature and boundary-yield proxies could inform a proposed model, but
they cannot be multiplied by column thickness and labeled physical capacity without new
constitutive assumptions and calibration. A phenomenological cohesion-times-thickness capacity
would be a separate alternative and should be labeled and tested as such.

The complete trace is in `backend/app/torque.py` (`gather_boundary_force_inputs`, force terms, `gather_per_node_force_field`, `predict_plate_motion`, `integrate_omega`), `backend/app/thin_shell.py` (`solve_force_field`, `solve_thin_shell`), `backend/app/rheology.py` (`normal_closing_rate_m_per_s`, `plastic_strain_rate_per_myr`), and `backend/app/merge_split.py` (`accumulate_plate_stress`, `maybe_split_plate`). `test_spatial_force_field_preserves_each_existing_torque_term` verifies the vector adapter's moment; `test_thin_shell.py` checks synthetic equilibrium, rigid-mode removal and linearity. No production damage path consumes this diagnostic.

## Gate disposition

- **Current diagnostic interface:** spatial tensile resultants and comparison to supplied integrated strength are available. Manufactured shallow-cap, L-shaped concave, generated-topology shallow-cap and original curved-mesh load cases have been exercised. Original curved-mesh solves converge, and four clamp placements show up to 23.5% mean-resultant sensitivity. The direct-threshold phenomenological rule has a monotonic capacity and weak-band localization check, but it is a candidate design only. Physical edge conditions, curved-shell accuracy against an independent benchmark and acceptance/replay calibration of the alternative remain open, so Gate 2 is not closed.
- **`internal_stress`:** retain only as a plate-scale eligibility/background signal if desired; it cannot trigger or weight local damage.
- **No calibrated strength field:** do not assign `Y_t`, `eta_rift`, or failure-strain parameters from these diagnostics. Do not infer stress from applied traction, rigid rotation, flow-cluster residuals or `Plate.internal_stress`.
- **Required next work:** validate the curved-shell formulation against an independent manufactured solution and broader mesh refinement; choose physically justified boundary conditions or approve the sensitivity-tested free-edge diagnostic; decide whether to adopt the labeled phenomenological capacity or continue toward a physical profile; then test the chosen load/capacity rule on replay candidates. Gate 3 accumulation and segment-growth calibration remains downstream of those choices.

The RFC now records this outcome and keeps the local-cut damage design conditional. Gate 1 is closed by the explicit baseline/exclusion decision recorded in `analysis/issue487/README.md`; the archived-stream provenance limitation remains documented.
