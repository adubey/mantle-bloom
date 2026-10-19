# RFC 0010 diagnosis gate 2: local loading available to rift damage

**Date:** 2026-10-19
**Scope:** Read-only trace of the force, motion, split-trigger, and boundary-deformation paths named by RFC 0010. No model equations or runtime behavior were changed.

## Result

The current code does **not** provide a defensible spatially varying intraplate tensile load for local rift damage. It has spatially sampled force terms and boundary-normal kinematics, but lacks the constitutive/load-transfer rule that maps those quantities to tensile stress across an interior cell adjacency. The plate-motion solve reduces the force field to a whole-plate angular acceleration. Consequently this gate cannot populate `sigma_t` or compare it with depth-integrated lithosphere strength. Gate 3 calibration must wait for either a force-interface extension and a repeated diagnosis, or a deliberate revision to a clearly labeled phenomenological rift model.

This is a diagnosis result, not evidence that local rifting is physically impossible. It says the current interface does not identify the needed field.

## Source trace and units

| Path | Available quantity | Reduction / limitation for local rifting |
|---|---|---|
| `torque.gather_boundary_force_inputs` | Per-own-node position, crust/mantle thickness, nearest-neighbor direction and distance, nearest neighbor plate/node, and neighbor angular velocity. | Geometry/material sampling only. It has no stress tensor, adjacency force, or load-transfer operator. Nearest-neighbor direction is not an intraplate tensile axis. |
| `slab_pull_force_n` | A scalar force in N on selected subducting boundary nodes: positive column buoyancy (Pa) times sunk length and represented boundary length. | A boundary-local driving force; not interior traction. `slab_pull_torque` crosses it with the planet-radius lever arm and sums it. |
| `ridge_push_torque` | Boundary-local force magnitude from isostatic relief, density contrast, gravity, and represented boundary length (N). | Direction is the local away-from-neighbor tangent, but only its torque is returned; there is no transfer into an intraplate stress field. |
| `basal_drag_coefficients` | Node-area-integrated basal drag. The underlying force per area has traction units (Pa), using mantle-flow and rigid-plate velocities. | Distributed basal traction is sampled node-wise, but is not an in-plane stress: the code integrates its moment over the plate. No shell equilibrium or elastic/viscous transfer maps it to extensional stress across links. |
| `collision_drag_coefficients` / `contact_drag_pa_s_per_m` | Boundary resistance based on relative plate velocity; coefficient in Pa·s/m times velocity gives Pa traction. | Applies to classified interplate collision/jammed contacts. It does not describe intraplate rift loading or damage. |
| `predict_plate_motion` / `integrate_omega` | Sum of explicit torques plus affine drag terms, integrated against whole-plate inertia to update one `omega` vector. | This is a plate-scale motion solve; it returns no spatial stress or strain-rate field. Rigid rotation alone has no internal strain. |
| `boundary_closing_rate_m_per_s` / `rheology.plastic_strain_rate_per_myr` | Boundary-normal relative plate velocity (m/s), mapped by an effective viscosity proxy to scalar Pa, then a calibrated plastic response. | Defined only for boundary-classified nodes. It is not a 2-D stress tensor or vertically integrated force, and cannot be reused as an intraplate rift field without a new physical rule. |
| `Plate.internal_stress` / `accumulate_plate_stress` | Dimensionless, time-integrated plate radius/overlap accumulator with decay. | Plate-level split eligibility relaxation, not a measured force/stress and not spatially locating damage. |

The complete trace is in `backend/app/torque.py` (`gather_boundary_force_inputs`, force terms, `predict_plate_motion`, `integrate_omega`), `backend/app/rheology.py` (`normal_closing_rate_m_per_s`, `plastic_strain_rate_per_myr`), and `backend/app/merge_split.py` (`accumulate_plate_stress`, `maybe_split_plate`). Existing unit tests verify force/torque formulas and boundary rheology, but do not define an intraplate equilibrium/load-transfer law.

## Gate disposition

- **Current force interface:** insufficient for local tensile stress or a load-versus-strength calculation.
- **`internal_stress`:** retain only as a plate-scale eligibility/background signal if desired; it cannot trigger or weight local damage.
- **No supported stress field:** do not assign `sigma_t`, `Y_t`, `eta_rift`, or failure-strain calibration values from the existing torque outputs. Do not distribute whole-plate torque by node area or infer stress from flow-cluster residuals without an explicit derivation.
- **Required next decision:** before gate 3, either (a) extend the interface with a documented spatial equilibrium/load-transfer rule, explicit force/traction/strength units and ownership, then compare local tensile load to integrated strength on synthetic plates; or (b) revise/reject the local damage premise and specify a phenomenological rift rule with its model boundary and validation criteria.

The RFC now records this outcome and keeps the local-cut damage design conditional. Gate 1's archived-stream provenance question remains separate and open.
