# Rifting and thin-shell diagnostic scripts

These offline diagnostics support RFC 0010's Gate 2 assessment of local loading and lithosphere
capacity, plus Gate 3's average-stress synthetic and history screens. Gate 2 is closed for the connected-source,
single-plate diagnostic. Its capacity results are conditional scenarios; they do not select a
production rupture threshold or validate coupled neighbor stresses. Gate 3 also includes a
cumulative sensitivity sweep over node-matched, 10 Myr post-event histories and a predeclared
one-hop growth screen on a third independent seed. These test damage sensitivity and spatial
coincidence under a frozen law; they do not calibrate physical parameters or simulate propagation
feedback. See the
[RFC](../../rfc/0010-rifts-follow-weak-lithosphere.md) and
[issue 487 findings](../../analysis/issue487/README.md).

Run scripts from the repository root with the worktree's `backend/.venv/bin/python`. Several tools
replay a full world and use substantial memory. Run one replay at a time and do not run it alongside
a test suite. Keep full saves and captures in scratch space; the compact outputs and findings belong
in `analysis/issue487/`. Each script's module docstring gives its complete usage and arguments.

### `audit_gate2_assumption_basis.py`

Calculates independent input-scale quantities for heat flow, strain rate, pore pressure and the
reference frictional K0 limit. It writes `gate2-assumption-basis.json`; it does not read event loads
or fit assumptions to old rift outcomes.

```bash
python3 bin/debug/audit_gate2_assumption_basis.py \
    --output analysis/issue487/gate2-assumption-basis.json
```

### `benchmark_gate2_coupled_shell.py`

Checks a manufactured two-region curved shell with a shared interface. The audit compares global
force and torque balance, interface reactions and recovered extensional resultants. It validates
the stated numerical benchmark, not a geological contact or finite-slip law.

### `calibrate_gate2_capacity.py`

Compares named column-capacity scenarios with published lithosphere strength scales and the matched
replay loads. The literature comparisons are order-of-magnitude checks with different spatial and
physical definitions; this is not a fit to legacy split locations or a calibrated global strength
field.

### `calibrate_gate3_average_stress.py`

Converts conditional Gate 2 resultant and capacity anchors to average Pa using the same column
thickness, then sweeps design choices for failure strain, damage viscosity and residual strength.
It checks sub-yield behavior, timestep partitioning, weak-column versus thick-column response, and
local segment nucleation, tip advance, arrest and linkage on a one-dimensional synthetic chain.
The parameter values are sensitivities, not physical calibrations or production defaults. The
two-threshold history command applies the proposed law retrospectively to the saved Gate 2
histories; those inputs are not a held-out dataset and retain their original 2 Myr cadence.

```bash
python3 bin/debug/calibrate_gate3_average_stress.py \
    --output analysis/issue487/gate3-average-stress.json
```

### `sweep_gate3_candidate_loads.py`

Extracts compact per-event median and p95 `N_ext/N_cap` ratios and average-strength scales from
the temporary Gate 2 replay outputs, then holds each snapshot constant over a duration and sweeps
failure strain, damage viscosity and residual strength. Its `history` command instead consumes the
node-matched JSONL from `replay_gate2_capacity.py --history-out`, integrates average stress and
strength over observed 2 Myr intervals, and sweeps the same damage parameters plus strength scales.
The `history-two-threshold` command applies the revised onset/completion law to the same histories.
Both history commands cover only 10 Myr after each connected-source candidate and are retrospective
sensitivities, not held-out validation or physical calibration.

```bash
python3 bin/debug/sweep_gate3_candidate_loads.py sweep \
    --snapshots analysis/issue487/gate3-candidate-load-snapshots.json \
    --output analysis/issue487/gate3-candidate-load-sweep.json

python3 bin/debug/sweep_gate3_candidate_loads.py history \
    --history /path/to/combined-gate3-histories.jsonl \
    --output analysis/issue487/gate3-cumulative-history-sweep.json

python3 bin/debug/sweep_gate3_candidate_loads.py history-two-threshold \
    --history /path/to/combined-gate3-histories.jsonl \
    --history-role heldout \
    --output analysis/issue487/gate3-two-threshold-history.json
```

### `evaluate_gate3_spatial_history.py`

Applies the frozen two-threshold reference law to high-cadence held-out node histories and
measures connected damage patches on each candidate band's source-cell adjacency graph. It uses
pre-topology coordinates and inherited-belt/craton markers, while taking accumulated damage from
post-topology histories. Patch size and great-circle diameter describe only the sampled candidate
band; the diagnostic has no tip propagation, damage feedback or full-plate linkage rule.

```bash
python3 bin/debug/evaluate_gate3_spatial_history.py \
    --history /path/to/issue487-g3-holdout-seed8675309-100myr.jsonl \
    --output analysis/issue487/gate3-spatial-seed8675309.json
```

### `evaluate_gate3_growth_history.py`

Applies the same fixed law to a held-out history that captured each candidate band plus a one-hop
source-cell adjacency collar. It classifies sampled threshold crossings as nucleation, adjacency to
one prior failed component, or linkage between multiple prior components. It reports collar outcomes,
connected patch span, simultaneous failures and a finite-window quiet label. The terms “tip advance”
and “window-arrested” are spatial/time-window diagnostics only: independent node damage has no stress
redistribution, triggering, production cuts, or propagation feedback.

```bash
python3 bin/debug/evaluate_gate3_growth_history.py \
    --history /path/to/issue487-g3-growth-seed27182818-100myr.jsonl \
    --output analysis/issue487/gate3-growth-seed27182818.json
```

### `test_gate3_magma_process_zone.py`

Checks a provisional Gaussian nonlocal tensile-stress driver at a pressurized dike tip. The
script recalibrates dike pressure to the same mean aperture on each mesh, then runs matched
magma-off/on near-tip seeded damage cases. The default caps nonlinear damage runs at 64 x 40 and
includes a 128 x 80 elastic screen; use `--mesh-nx-levels 128 --damage-max-nx 128` for the
fine-mesh damage pair. The 25 km process-zone length is inherited from Gate 3 damage regularization;
it is not a physical calibration. Use `--process-zone-length-km` to vary it. The sparse spatial
neighbor construction preserves the Gaussian weighting while avoiding large Python neighbor lists.
See the RFC 0010 calibration report and its fine-mesh/length follow-up for results and limitations.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/test_gate3_magma_process_zone.py \
    --output analysis/issue487/gate3-magma-process-zone.json
```

### `test_gate3_spatial_pressure_growth.py`

Screens a finite-compliance magma source feeding a mode-I dike with a linear center-to-tip
pressure profile. The same profile is projected into the crack-tip stress intensity and an
elliptical elastic opening mode; source depletion changes volume, mean aperture, slot flow and the
leading tensile stress ahead of the tip. This is a reduced analytical diagnostic: the pressure
profile is prescribed, slot resistance uses mean aperture, and it does not solve source-side host
stress relief or a lithosphere-scale fracture process zone. Its parameter grid is not calibrated.
See the RFC 0010 Gate 3 calibration report for the supercritical-seed result and limits.

```bash
backend/.venv/bin/python bin/debug/test_gate3_spatial_pressure_growth.py \
    --output analysis/issue487/gate3-spatial-pressure-growth.json
```

### `audit_gate3_pressure_toughness_basis.py`

Compares uniform and linear center-to-front pressure profiles at the same observed dike opening
scale. It reports equivalent elastic net pressure, tip stress intensity, energy release rate and
critical initial flaw lengths across separate laboratory and dike-scale effective-toughness
scenarios. Opening-derived pressure is not measured magma overpressure; the output is a conditional
parameter-basis audit, not a calibration.

```bash
backend/.venv/bin/python bin/debug/audit_gate3_pressure_toughness_basis.py \
    --output analysis/issue487/gate3-pressure-toughness-basis.json
```

### `test_gate3_distributed_dike_flow.py`

Maps spatial pressure values through a small-strain finite-element elastic crack and iterates them
with distributed Poiseuille flow. It uses a quasi-steady kinematic front-speed condition and
compares the resulting opening, volume, pressure, tip stress intensity and event-scale volume/flow
duration against the Afar 2005 geometry. The reported seismicity migration rate is only a speed
proxy; the dike geometry is fixed and the reservoir depletion is an equivalent compliance budget,
not a transient propagation simulation.

```bash
backend/.venv/bin/python bin/debug/test_gate3_distributed_dike_flow.py \
    --output analysis/issue487/gate3-distributed-dike-flow.json
```

### `simulate_gate3_transient_dike_front.py`

Runs a parameter sweep with evolving crack half-length, distributed lubrication pressure,
finite-reservoir depletion, flow-limited front speed, and a toughness arrest condition. It compares
event duration and volume with the Afar 2005 dike, and reports the leading LEFM tensile stress 1 km
ahead of the tip. The elastic opening remains a single-mode projection; the near-tip stress is not
a spatial host stress solution and is not fed back into propagation.

```bash
backend/.venv/bin/python bin/debug/simulate_gate3_transient_dike_front.py \
    --output analysis/issue487/gate3-transient-dike-front.json
```

### `prototype_gate3_transient_host_stress.py`

Transfers two snapshots of the representative transient dike's distributed pressure into the
resolved host-elastic cutout. It compares stress-tensor increments against intact controls under
zero and 5 MPa nominal tensile background loading. A near-tip stress fit replaces the analytical
`K` at the two zero-prestress snapshots and estimates its effect on one transient path. This is a
post-processed two-snapshot sensitivity, not iterative host-flow feedback; the prestress is a
sensitivity scenario and damage does not evolve.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_transient_host_stress.py \
    --output analysis/issue487/gate3-transient-host-stress.json
```

### `sweep_gate3_near_arrest_feedback.py`

Sweeps initial flaw length, reservoir compliance and effective toughness using the cutout-derived
tip-`K` correction. It compares analytical propagation, central feedback and a fit-window
sensitivity envelope. The correction is interpolated from two spatial snapshots, so this is a
near-threshold sensitivity test rather than a fully coupled or mesh-converged propagation model.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/sweep_gate3_near_arrest_feedback.py \
    --host-output analysis/issue487/gate3-transient-host-stress.json \
    --output analysis/issue487/gate3-near-arrest-feedback.json
```

### `converge_gate3_tip_k_cutout.py` and `rerun_gate3_threshold_feedback.py`

Compares tip-stress-intensity fits at 0.5 and 0.25 km cutout spacing, then independently computes
a finite-strain domain J with the applied dike-face pressure correction at four contour radii.
The companion rerun uses fine-grid domain-J equivalent K as its central correction and keeps the
full stress-fit window extrema as a sensitivity envelope. The J-to-K conversion assumes local
plane-strain mode-I behavior. The coarse parent and fine cutout receive the same transient
distributed pressure profile; output includes their aperture, pressure-resultant and pressure-work
mismatches. This remains post-processed feedback: the host cutout is not solved at every
propagation step.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/converge_gate3_tip_k_cutout.py \
    --output analysis/issue487/gate3-tip-k-cutout-convergence.json
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/rerun_gate3_threshold_feedback.py \
    --cutout-output analysis/issue487/gate3-tip-k-cutout-convergence.json \
    --output analysis/issue487/gate3-domain-j-threshold-feedback.json
```

### `sweep_gate3_parent_cutout_transfer.py`

Measures how parent-to-cutout aperture, pressure-work and integrated-load differences change as
the coarse parent mesh is refined from 64 x 130 to 512 x 520. It holds the local mesh at 0.25 km,
uses two short seed states and a longer arrest control, and applies the same transient distributed
pressure profile on both meshes. `--parent-nx`, `--parent-ny` and `--snapshot-labels` select a
single refined mesh and subset of cases. It is a transfer-convergence diagnostic, not propagation
validation.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/sweep_gate3_parent_cutout_transfer.py \
    --output analysis/issue487/gate3-parent-cutout-transfer.json
```

### `prototype_gate3_coupled_dike_propagation.py`

Runs one fixed case-A history with FE aperture feeding cubic-law flow resistance and FE volume
feeding finite-reservoir storage at each sampled length. It iterates aperture, pressure and source
pressure with under-relaxation, then uses matched parent/cutout mechanics and pressure-corrected J
for tip advance or arrest. This converges the sampled feedback states, but the case reaches only
0.848 km³, below the 1–2 km³ comparison range, so it does not close Gate 3. To continue an
interrupted run whose output has `coupling_failed_to_converge`, pass it with `--resume-from`.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_coupled_dike_propagation.py \
    --output analysis/issue487/gate3-coupled-aperture-feedback.json
```

### `refine_gate3_arrest_bracket.py`

Warm-starts an inserted total-length state from a saved coupled-history aperture profile and compares
the baseline fixed-point settings with a tighter tolerance and lower relaxation. It prints iteration
residuals and checkpoints each completed case. Use `--target-half-length-cells` to select a
parent-grid-aligned half-length. `--resolution-scale 2` uses a half-resolution
parent and cutout mesh for a provisional, faster check; results from that run must not be combined
with the high-resolution arrest trajectory as if they were a mesh-converged bracket. The default
resolution matches the 512 x 520 coupled history and can be used to refine its arrest bracket.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/refine_gate3_arrest_bracket.py --target-half-length-cells 121 \
    --input analysis/issue487/gate3-coupled-aperture-feedback.json \
    --output analysis/issue487/gate3-arrest-bracket-60-5km.json
```

### `prototype_gate3_resolved_dike_tip.py`

Runs a source-only elastic cutout with 0.5 km spacing around a complete 65 km pressurized dike.
The cutout inherits displacement from matched coarse dike-on/off solutions only on its outer
perimeter; the internal crack faces remain free. It reports stress-increment contours at
sensitivity thresholds and checks the local mean aperture against the coarse calibration.
Observed seismic/fault-activation dimensions are comparison scales, not assigned damage lengths
or calibrated stress thresholds. This does not evolve damage or establish dike propagation.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_resolved_dike_tip.py \
    --output analysis/issue487/gate3-resolved-dike-tip.json
```

### `prototype_gate3_resolved_dike_damage.py`

Applies the existing Gate 3 finite-strain damage law to a 1 km resolved dike cutout under a
90%-of-yield tectonic prestress and a finite 12-hour pressure episode, with an intact prestress-only
control. This is a compatibility screen, not propagation evidence: the present free-crack model
opens by tens of metres under that prestress, far beyond the pressure calibration's 4 m target.
See the issue 487 calibration report before interpreting peak stresses.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_resolved_dike_damage.py \
    --output analysis/issue487/gate3-resolved-dike-damage.json
```

### `prototype_gate3_cohesive_dike_interface.py`

Sweeps a linear bilateral normal-spring seam in both the coarse parent and 1 km local dike
cutout, recalibrates magma pressure at a fixed prestress to a 4 m aperture, then runs one
12-hour damage step for a selected seam stiffness. The seam values are sensitivity probes, not
independently calibrated interface properties. The output shows how pressure-aperture fit and
damage response change with that assumption; it is not evidence of propagation.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_cohesive_dike_interface.py \
    --output analysis/issue487/gate3-cohesive-dike-interface.json
```

### `prototype_gate3_frictional_fault.py`

Splits a resolved dike cutout on a 2-D in-plane fault interface and runs matched dike-only,
fault-only and fault-plus-dike equilibria over 50%, 80% and 95% Coulomb prestress scenarios.
The effective normal stress, friction, interface penalty and 5 MPa dike load are design inputs,
not Miyakejima measurements. The script supports `--spacing-km` for mesh sensitivity. This does
not reproduce the mapped 3-D strike-slip geometry or predict dike arrest; use it to inspect the
interface law and numerical sensitivity only.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_frictional_fault.py --spacing-km 1.0 \
    --output analysis/issue487/gate3-frictional-fault.json
```

### `prototype_gate3_coupled_dike_fault_history.py`

Solves a staged 2-D dike and preloaded Coulomb fault in one finite-strain equilibrium model. A
declining pressure schedule is supplied per dike-length stage; the default spans 16 to 2.5 MPa,
historical same-edifice feeder-dike analogues rather than measured 2000-event pressures. Cumulative
fault slip shifts the Coulomb spring's shear-stress offset between stages, while host displacement
is re-equilibrated for each changed crack geometry. The energy estimate uses mean pressure times the
change in integrated opening volume, less the host elastic-energy change and fault-interface work.
The reduced geometry, partial history carry, pressure assumptions, and FE energy estimate do not
reproduce the published displacement-discontinuity model or establish Miyakejima arrest. Supply
paired `--half-lengths-km` and `--pressure-history-mpa` lists to define a scenario.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_coupled_dike_fault_history.py \
    --spacing-km 1.0 \
    --half-lengths-km 10 15 20 25 30 \
    --pressure-history-mpa 16 12.6 9.3 5.9 2.5 \
    --output analysis/issue487/gate3-coupled-dike-fault-source-history.json
```

### `prototype_gate3_finite_reservoir_tip_advance.py`

Keeps one potential dike seam in a fixed finite-strain mesh. Cohesive springs tie the portions
ahead of the two active tips; each growth step frees a longer interval, carries the previous host
displacement and Coulomb interface jumps forward, and solves reservoir pressure from the opening
volume against a finite 1.04 km³ transfer budget. Its stepwise pressure work minus host energy
change is an upper-bound energy screen because explicit frictional dissipation is omitted. The
10 km strike width, 5 MPa initial reservoir pressure, fault prestress and 0.2 MJ/m² reference
fracture energy are stated scenarios. These runs are a reduced mechanism test and do not reproduce
the event pressure history or establish Miyakejima arrest.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_finite_reservoir_tip_advance.py \
    --spacing-km 1.0 --growth-step-km 2.0 \
    --output analysis/issue487/gate3-finite-reservoir-tip-advance-1km.json
PYTHONPATH=backend:bin/debug backend/.venv/bin/python \
    bin/debug/prototype_gate3_finite_reservoir_tip_advance.py \
    --spacing-km 0.5 --growth-step-km 2.0 \
    --output analysis/issue487/gate3-finite-reservoir-tip-advance-0p5km.json
```

### `prototype_gate3_mapped_3d_fault_transfer.py`

Uses the finite rectangular dike and left-lateral fault patches inferred for Miyakejima in
Ozawa et al. (2004), Table 1, in a homogeneous full-space elastic Kelvin model. It integrates
dike-to-fault normal and strike-shear stress increments over both finite patches while refining
midpoint quadrature. The inferred patches intersect, so stresses near their intersection are
singular in this representation. The script excludes 1, 2 or 4 km around the dike patch from the
fault average and reports those lengths as uncalibrated sensitivities. It estimates elastic
Coulomb tendency only; it does not solve slip, source evolution or arrest, and omits the free
surface and crustal structure.

```bash
backend/.venv/bin/python bin/debug/prototype_gate3_mapped_3d_fault_transfer.py \
    --output analysis/issue487/gate3-mapped-3d-fault-transfer.json
```

### `audit_gate3_miyakejima_reference_profiles.py`

Packages a 30 km terrain/bathymetry trace along the preferred full-event Hughes (2010)
inversion geometry. The sampled elevations come from the public NOAA ETOPO 2022 15-arc-second
grid; each 1 km profile point records the nearest native grid-cell coordinate and its offset
from the ideal trace. It also records an approximate digitization of Hughes Figure 46's observed
GPS displacement perpendicular to the dike. That curve is an auxiliary surface-deformation
check, not an aperture target. The trace is the straight inversion geometry, not a digitization
of Maccaferri's two-segment propagation path, and does not reproduce the paper's load conversion.

```bash
python3 bin/debug/audit_gate3_miyakejima_reference_profiles.py \
    --output analysis/issue487/gate3-miyakejima-reference-profiles.json
```

### `diagnose_gate2_column_envelope.py`

Computes brittle Mohr–Coulomb normal-fault strength and temperature/rate-dependent ductile flow
stress through representative lithosphere columns, then integrates the weaker additional strength
with depth. It reports a scenario envelope for missing local temperature, fluid, composition,
strain-rate, friction and background-stress state. It is a diagnostic rheology model, not production
rheology.

### `replay_gate2_assumption_screen.py`

Post-processes captured Gate 1 node loads without another world replay. It crosses the reference
stress, dry/wet mantle rheology and two pore-pressure profiles; its output records scenarios,
classification changes and input hashes. The near-critical K0 case is derived from the reference
friction angle and is a limiting sensitivity, not a measured stress estimate or probability bound.

### `replay_gate2_capacity.py`

Replays seeds 7, 11 and 2024 against the pinned Gate 1 source with diagnostic helpers overlaid. It
checks that candidate identities match Gate 1, builds a shell load for each craton-crossing event,
and compares physical column scenarios with a labeled phenomenological capacity alternative. The
capacity scores do not feed back into the replay. Use one seed at a time when memory is limited and
check the output's baseline count-match fields before interpreting results.
Its optional `--history-out` captures accepted candidate cut-band node IDs, pre-topology spatial
metadata and source-cell links, then follows them after topology across the requested number of
steps. Use `--step-years 2.5e5 --history-steps 40` for the 0.25 Myr, 10 Myr held-out Gate 3
screen. For RFC 0010 Gate 3, run against the pinned Gate 1 source with diagnostic helpers overlaid;
later production model revisions can change candidate populations and invalidate comparison to
the stored Gate 2 events.

### `replay_gate2_coupled_event.py`

Replays one event for a source plate and its dominant force-contributing neighbor, then compares
pair-level stress under current contact loads and a bonded solve without the pair's mutual contact
drag. It also reports force/torque closure and rigid-mode projection sensitivity. These cases are
mechanical sensitivity probes, not a validated interface law or physical coupled-stress prediction.

### `sensitivity_gate2_capacity.py`

Re-scores captured node loads over one-factor-at-a-time capacity alternatives. It varies the
assumptions recorded in the script and compares each result against the same load samples. The
result is a sensitivity study, not a joint uncertainty distribution or fitted rupture label.

### `validate_thin_shell_plate_meshes.py`

Builds representative meshes from generated plates and checks clipping/coverage, manufactured
loads on the original curved mesh and a shallow-cap mapping, localized and paired clamp sensitivity,
and an independent spherical-shell solution. Its capacity response sweep checks a phenomenological
rule's behavior only; the chosen values are not physical strength estimates.

### `prototype_gate3_2d_viscoplastic_section.py`

Solves one small-strain, quasi-static `x-z` section with local force equilibrium, a reduced
plane-strain stress mapping, and Gate 2's depth-dependent brittle and creep branches. Its default
command compares uniform intact rock with a prescribed initially thinned corridor on 10 km and
5 km grids. This is a controlled mechanics diagnostic; it has no evolving fault, crack, magma,
gravity, or transient temperature field.

```bash
backend/.venv/bin/python bin/debug/prototype_gate3_2d_viscoplastic_section.py \
    --output analysis/issue487/gate3-2d-viscoplastic-section.json
```

### `run_gate3_2d_viscoplastic_transient.py`

Runs a short paired sequence on the 10 km grid. After each quasi-static solve, it integrates
elementwise vertical strain into local lithosphere thickness for the next snapshot. The mesh is
fixed and material is not advected; strong element-scale thinning therefore needs refinement and
a physical localization length before geological interpretation. This is an exploratory
feedback diagnostic, not a rift-growth validation. Final 2-D maps of remaining thickness and
horizontal/vertical strain rate are saved for profile comparison. `--nx` and `--nz` select the
horizontal and vertical cell counts; for example, `--nx 60 --nz 32` repeats it at 5 km spacing.

```bash
backend/.venv/bin/python bin/debug/run_gate3_2d_viscoplastic_transient.py \
    --output analysis/issue487/gate3-2d-viscoplastic-transient.json
```

### `verify_gate3_material_advection.py`

Checks marker trajectories and deformation-gradient updates in manufactured affine extension
and simple shear against exact solutions. The extension position and area errors should decrease
with timestep refinement; simple shear should preserve area. This is the first kinematic
verification step only: it does not implement particle-to-grid interpolation/remapping, material
transfer conservation, or coupling to the Gate 3 stress solve.

```bash
backend/.venv/bin/python bin/debug/verify_gate3_material_advection.py \
    --output analysis/issue487/gate3-material-advection.json
```

### `verify_gate3_marker_grid_transfer.py`

Advects labeled markers through manufactured pure-shear extension and simple shear. At each step
it deposits mass, volume, per-material inventories and history to a fixed grid using bilinear
partition-of-unity weights, samples a prescribed affine grid velocity back at the markers, and
advects them. It checks per-step conservation, history reconstruction and trajectory convergence.
This verifies the transfer cycle only: velocity is prescribed, not stress-solved, and remapped
state is not coupled to Gate 3 mechanics.

```bash
backend/.venv/bin/python bin/debug/verify_gate3_marker_grid_transfer.py \
    --output analysis/issue487/gate3-marker-grid-transfer.json
```

### `run_gate3_section_marker_coupling.py`

Runs uniform and inherited-weak-corridor controls. Each step solves the reduced FE section,
advects Lagrangian material markers with its velocity, updates marker deformation gradients from
the FE velocity gradient, and projects marker-carried layer codes and vertical stretch back to the
elements for the next constitutive solve. The current default is 30 x 16 cells over 1 Myr with four
steps. Markers start away from the fixed boundaries; elements without markers use their prior
Eulerian layer and strain state. The diagnostic has no moving boundary, fracture, transient heat,
gravity, magma or mantle flow, and it is not a Gate 3 pass.

```bash
backend/.venv/bin/python bin/debug/run_gate3_section_marker_coupling.py \
    --output analysis/issue487/gate3-section-marker-coupling.json
```

### `run_gate3_section_marker_moving_mesh.py`

Runs uniform and inherited-weak-corridor controls on a Lagrangian moving FE mesh. Every triangle
starts with three material markers; mesh nodes and markers move with the stress-solved FE velocity,
and marker deformation gradients plus layer codes feed into the next constitutive solve. It checks
element coverage, marker/mesh position agreement, area closure, force residual and whether the
seeded weak corridor localizes extension. The default is 30 x 16 cells, 1 Myr and four steps. Use
`--nx`, `--nz`, `--steps` and `--duration-myr` for resolution and timestep comparisons. The
optional `--max-viscosity-pa-s` varies the effective-viscosity upper bound; its default remains
`1e28 Pa s`.

The output reports per-layer FE-realized equivalent stress separately from the flow-law stress
target, which may differ when the effective viscosity hits its bound. It identifies the FE stress
peak's location, temperature, pressure, strain rate and controlling branch. If a constitutive Picard
solve does not converge, the driver stops without advancing that step's mesh or markers and records
the completed duration and failed state. Convergence requires both a small log-viscosity update
and a final-response free-node force residual below `0.25%`; if the viscosity fixed point meets its
tolerance but the force residual does not, Picard iteration continues. `--max-picard-iterations`
can raise the solver iteration ceiling for a sensitivity run; an incomplete run is not a valid
full-duration comparison.

The diagnostic has no ALE remeshing or smoothing. The corridor is imposed, and it does not include
fracture, faults, transient heat, gravity, magma or mantle flow. Passing its internal numerical
checks is not a Gate 3 pass or geological validation.

```bash
backend/.venv/bin/python bin/debug/run_gate3_section_marker_moving_mesh.py \
    --output analysis/issue487/gate3-section-marker-moving-mesh.json
backend/.venv/bin/python bin/debug/run_gate3_section_marker_moving_mesh.py \
    --nx 60 --nz 32 --steps 4 \
    --output analysis/issue487/gate3-section-marker-moving-mesh-5km.json
backend/.venv/bin/python bin/debug/run_gate3_section_marker_moving_mesh.py \
    --nx 60 --nz 32 --steps 128 --duration-myr 1 --max-picard-iterations 800 \
    --corridor-only --progress --checkpoint-output /tmp/gate3-section-final.npz \
    --output /tmp/gate3-section-run.json
backend/.venv/bin/python bin/debug/run_gate3_section_marker_moving_mesh.py \
    --nx 60 --nz 32 --steps 128 --duration-myr 1 --max-picard-iterations 800 \
    --max-viscosity-pa-s 3e28 --corridor-only --progress \
    --output /tmp/gate3-section-high-viscosity-run.json
```

The optional checkpoint stores the final mesh and marker deformation/material state in compressed
NumPy format. Use it with `audit_gate3_viscosity_ceiling_same_state.py` to vary only the effective
viscosity upper bound while holding the final state and loading fixed. The check reports clipped
element counts, realized and target stress, and the free-node force residual. This is a static
constitutive sensitivity, not a rerun with the new ceiling and not a Gate 3 pass. The current
`1e28 Pa s` bound clips 250 elements in the saved 5 km, 128-step endpoint; `3e28 Pa s` removes those
clips and matches the stress target in that state. The compact result is
[`gate3-viscosity-ceiling-same-state.json`](../../analysis/issue487/gate3-viscosity-ceiling-same-state.json).

A complete 5 km, 128-step corridor trajectory at `3e28 Pa s` converges with no elements hitting
the upper bound. It raises the final realized stress from `200 MPa` to the `406 MPa` target, but
final thinning and localization are effectively unchanged from the `1e28 Pa s` run. The ceiling
therefore affected stress reporting in the static audit but did not materially limit geometry
evolution over this 1 Myr trajectory. See the [compact comparison](../../analysis/issue487/gate3-section-marker-moving-mesh-5km-128step-3e28-summary.json).

```bash
backend/.venv/bin/python bin/debug/audit_gate3_viscosity_ceiling_same_state.py \
    --checkpoint /tmp/gate3-section-final.npz \
    --output analysis/issue487/gate3-viscosity-ceiling-same-state.json
```

### `audit_gate3_rate_floor_same_state.py`

Repeats a static equilibrium solve on one saved marker/mesh endpoint while varying the
equivalent-strain-rate floor. It reports how many elements' raw rates fall below the floor, stress
target versus realized stress, viscosity-bound counts, and free-node force residual. Lower floors
can require larger viscosities, so use `--maximum-viscosity-pa-s` to check whether the ceiling is
confounding the comparison. In the current 5 km endpoint, `1e-21 s-1` with a `1e32 Pa s` ceiling
converges without viscosity clipping and lowers peak target stress by about 7%; at `1e-22 s-1` the
current Picard solver does not converge. These are numerical sensitivities, not physical
calibrations or a Gate 3 pass.

```bash
backend/.venv/bin/python bin/debug/audit_gate3_rate_floor_same_state.py \
    --checkpoint /tmp/gate3-section-final.npz \
    --maximum-viscosity-pa-s 3e28 \
    --output analysis/issue487/gate3-rate-floor-same-state.json
backend/.venv/bin/python bin/debug/audit_gate3_rate_floor_same_state.py \
    --checkpoint /tmp/gate3-section-final.npz \
    --maximum-viscosity-pa-s 1e32 \
    --output analysis/issue487/gate3-rate-floor-same-state-1e32.json
```

### `audit_gate3_section_stress_basis.py`

Reads the peak-stress record from the moving-mesh audit and evaluates 108 combinations of the Gate 2
scenario assumptions for local mantle brittle-plastic and dry-olivine creep strength. It varies base
temperature, `K0`, friction coefficient, cohesion and pore-pressure profile while holding geometry,
pressure and strain rate at the sample location fixed. It reports the branch target and FE-realized
stress after viscosity bounds, plus the designed-grid percentile range. This is a local constitutive
sensitivity screen: it does not re-solve FE equilibrium, assign probabilities to the scenarios, or
calibrate tectonic stress.

```bash
backend/.venv/bin/python bin/debug/audit_gate3_section_stress_basis.py \
    --input analysis/issue487/gate3-section-marker-moving-mesh-stress-audit.json \
    --output analysis/issue487/gate3-section-stress-basis-sensitivity.json
```
