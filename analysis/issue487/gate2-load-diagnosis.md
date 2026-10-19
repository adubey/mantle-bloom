# RFC 0010 diagnosis gate 2: local loading available to rift damage

**Date:** 2026-10-19
**Scope:** Read-only force/motion trace plus a diagnostic source-force adapter and reduced thin-shell equilibrium solver. Production motion, deformation, topology, and save behavior are unchanged.

## Result

The source-force interface preserves each existing torque contribution as a per-node vector in N. `thin_shell.solve_force_field` transfers the combined tangent load through a linear triangular membrane equilibrium solve and returns local tensile stress resultants in N/m. It removes the area-weighted rigid-rotation load handled by plate motion; applied radial force is assigned to shell-normal support. The result includes a ratio to caller-supplied depth-integrated tensile capacity in N/m. Synthetic tests verify equilibrium residual, rigid-rotation removal and linear force scaling.

This establishes a diagnostic load-transfer interface, not a validated production stress field. The solver assumes outward convex-hull triangular elements, linear plane-stress elasticity, free edges and no faults, bending, plasticity or thermal strength envelope. A manufactured uniaxial load on one shallow spherical cap recovers its analytic resultant across three refinements; other plate shapes, boundary conditions and the material strength provider remain unvalidated. Gate 2 therefore remains open; Gate 3 calibration must wait for those checks or a deliberate revision to a clearly labeled phenomenological rift model.

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
closed-form benchmark, so this is a sensitivity check, not an accuracy claim. Nonconvex outlines,
other edge placements and actual plate meshes remain untested.

## Source trace and units

| Path | Available quantity | Reduction / limitation for local rifting |
|---|---|---|
| `torque.gather_boundary_force_inputs` | Per-own-node position, crust/mantle thickness, nearest-neighbor direction and distance, nearest neighbor plate/node, and neighbor angular velocity. | Geometry/material sampling only. It has no stress tensor, adjacency force, or load-transfer operator. Nearest-neighbor direction is not an intraplate tensile axis. |
| `torque.gather_per_node_force_field` / `PerNodeForceField` | Per-node slab-pull, slab-drag, ridge-push, basal-drag and contact-drag force vectors (N); total external force and area-average applied traction (Pa). | Preserves the current spatial external forcing before torque aggregation. The vector sum reproduces the current torque terms, but `total_traction_pa` is not internal tensile stress. |
| `thin_shell.solve_force_field` / `solve_thin_shell` | Nodal membrane stress resultants (N/m), maximum tensile principal resultant, supplied integrated strength (N/m), and failure ratio. | Solves a reduced linear elastic equilibrium on outward-facing convex-hull facets. The shallow-cap manufactured load is recovered within 0.16% over 3 × 3 through 9 × 9 nodes. An irregular convex cap also passes; clamping its west edge changes the mean result by 2.82% with a small free-node residual. Radial load is excluded and rigid rotation is removed. Nonconvex outlines, other edge placements, constitutive and omitted-process sensitivity remain to be validated. |
| `slab_pull_force_n` | A scalar force in N on selected subducting boundary nodes: positive column buoyancy (Pa) times sunk length and represented boundary length. | A boundary-local driving force; not interior traction. `slab_pull_torque` crosses it with the planet-radius lever arm and sums it. |
| `ridge_push_force_vectors_n` / `ridge_push_torque` | Boundary-local force magnitude from isostatic relief, density contrast, gravity, and represented boundary length (N). | Direction is the local away-from-neighbor tangent; it is now retained as a vector for load transfer before moment aggregation. |
| `basal_drag_force_vectors_n` / `basal_drag_coefficients` | Node-area-integrated basal drag. The underlying force per area has traction units (Pa), using mantle-flow and rigid-plate velocities. | Distributed basal force is retained node-wise and transferred by the diagnostic shell solve; it is not itself in-plane stress. |
| `collision_drag_coefficients` / `contact_drag_pa_s_per_m` | Boundary resistance based on relative plate velocity; coefficient in Pa·s/m times velocity gives Pa traction. | Applies to classified interplate collision/jammed contacts. It does not describe intraplate rift loading or damage. |
| `predict_plate_motion` / `integrate_omega` | Sum of explicit torques plus affine drag terms, integrated against whole-plate inertia to update one `omega` vector. | This is a plate-scale motion solve; it returns no spatial stress or strain-rate field. Rigid rotation alone has no internal strain. |
| `boundary_closing_rate_m_per_s` / `rheology.plastic_strain_rate_per_myr` | Boundary-normal relative plate velocity (m/s), mapped by an effective viscosity proxy to scalar Pa, then a calibrated plastic response. | Defined only for boundary-classified nodes. It is not a 2-D stress tensor or vertically integrated force, and cannot be reused as an intraplate rift field without a new physical rule. |
| `Plate.internal_stress` / `accumulate_plate_stress` | Dimensionless, time-integrated plate radius/overlap accumulator with decay. | Plate-level split eligibility relaxation, not a measured force/stress and not spatially locating damage. |

The complete trace is in `backend/app/torque.py` (`gather_boundary_force_inputs`, force terms, `gather_per_node_force_field`, `predict_plate_motion`, `integrate_omega`), `backend/app/thin_shell.py` (`solve_force_field`, `solve_thin_shell`), `backend/app/rheology.py` (`normal_closing_rate_m_per_s`, `plastic_strain_rate_per_myr`), and `backend/app/merge_split.py` (`accumulate_plate_stress`, `maybe_split_plate`). `test_spatial_force_field_preserves_each_existing_torque_term` verifies the vector adapter's moment; `test_thin_shell.py` checks synthetic equilibrium, rigid-mode removal and linearity. No production damage path consumes this diagnostic.

## Gate disposition

- **Current diagnostic interface:** spatial tensile resultants and comparison to supplied integrated strength are available. Manufactured shallow-cap refinement and one irregular convex outline/free-versus-clamped edge comparison pass. Nonconvex geometry, broader boundary sensitivity and a material strength provider are not established, so this is not ready for damage calibration.
- **`internal_stress`:** retain only as a plate-scale eligibility/background signal if desired; it cannot trigger or weight local damage.
- **No calibrated strength field:** do not assign `Y_t`, `eta_rift`, or failure-strain parameters from these diagnostics. Do not infer stress from applied traction, rigid rotation, flow-cluster residuals or `Plate.internal_stress`.
- **Required next work:** derive a depth-integrated column strength field with explicit units and ownership; check nonconvex/representative mesh and broader edge sensitivity, then compare load/capacity on additional synthetic geometries and replay candidates. Otherwise revise/reject the local damage premise and specify a phenomenological rift rule with its model boundary and validation criteria.

The RFC now records this outcome and keeps the local-cut damage design conditional. Gate 1 is closed by the explicit baseline/exclusion decision recorded in `analysis/issue487/README.md`; the archived-stream provenance limitation remains documented.
