# Gate 3 controlled pressure-source comparison

**Date:** 2026-10-25

**Disposition:** The corrected comparison isolates pressure on an existing crack. Pressure raises tip/collar tensile stress in both matched pairs, but its source-flank effect changes sign with the damage seed. The earlier claim of consistent pressure-driven flank relief is withdrawn. Gate 3 remains open.

## Question and control

Under the same extension loading, crack geometry, boundary conditions, and initial damage, does applying the calibrated dike pressure change source-flank and tip/collar stresses?

This follow-up corrects two issues in the merged #570 evidence. First, its seeded no-dike table entry used the 2.5 Myr endpoint while labeling it 0.2 Myr. Second, the pressure-on case used duplicated, disconnected crack faces while the comparison case used intact mesh. The pressure-specific pairs below now use the same cracked mesh and toggle only the pressure load. The zero-pressure crack is an artificial control that isolates pressure; it is not a physical no-dike host. Separate intact-host runs provide context for the added crack geometry; they are not the pressure control.

All cases fix horizontal displacement on both vertical edges and vertical displacement at the bottom-left and bottom-right corners. The pressure pairs also share extension rate, mesh topology, and seed state. This meets the source-attribution expectation in design principle 9: matched controls isolate the contribution being measured, and equal-time outcomes are reported.

The calculation uses a 400 × 130 km plane-strain strip, a 32 × 20 mesh, 4 mm/yr extension, 0.02 Myr steps, and a requested 2.5 Myr duration. Pressure is calibrated in the elastic mesh to the 2005 Afar dike's approximate 65 km length and 4 m mean opening; peak opening is checked against the observed ~8 m bound. The source holds pressure constant and is neither a finite-reservoir nor a plume-on/off experiment.

## Results at equal elapsed time

All stress values below are sampled at 0.2 Myr.

| Geometry | Pressure | Tip seed | Source-flank mean σxx | Tip/collar peak tensile stress | Stop / failures outside seed |
|---|---|---|---:|---:|---|
| Intact host | Off | No | 61.93 MPa | 61.93 MPa | Completes 2.5 Myr; 0 |
| Cracked mesh | Off | No | 27.79 MPa | 299.31 MPa | Stops at 0.24 Myr on stress guard; 0 |
| Cracked mesh | On | No | 34.47 MPa | 325.70 MPa | Stops at 0.20 Myr on stress guard; 0 |
| Intact host | Off | Yes | 60.83 MPa | 69.73 MPa | Completes 2.5 Myr; 0 outside 16 inherited failures |
| Cracked mesh | Off | Yes | 37.18 MPa | 328.12 MPa | Stops at 0.20 Myr on stress guard; 0 outside 16 inherited failures |
| Cracked mesh | On | Yes | 35.78 MPa | 339.13 MPa | Stops at 0.20 Myr on stress guard; 0 outside 16 inherited failures |

The matched pressure effect is seed-dependent. Without the seed, pressure increases source-flank mean σxx by 6.68 MPa and the tip/collar peak by 26.40 MPa. With the seed, pressure decreases source-flank mean σxx by 1.40 MPa and increases the tip/collar peak by 11.01 MPa. The original intact-host versus pressurized-crack comparison mixed pressure and crack-geometry effects; its source-flank relief cannot be attributed to pressure alone.

The elastic pressure-source check remains: pressure is 3.961 MPa, mean opening is 4.000 m, maximum opening is 5.762 m, equilibrium residual is `4.54e-13`, source/support force and moment residuals are `2.60e-13` and `5.19e-13`, and pressure-work/stored-energy mismatch is `8.05e-7`. The minimum deformation determinant in the seeded pressure-on case is positive (`0.9983`).

Independent constitutive checks remain green for affine stretch, the small-strain limit, elastic stress objectivity, the consistent tangent, implicit damage update, and plastic-volume preservation ([verification JSON](gate3-finite-strain-verification-review.json)). The rotation check holds plastic deformation at identity; plastic-flow update objectivity under finite rotations is not tested.

## Interpretation and next gate work

For this prototype, crack geometry causes most of the measured tip concentration and source-flank stress reduction relative to intact host. Adding pressure to that same crack consistently increases the tip/collar peak, but does not consistently reduce source-flank stress. No pressure-on or cracked seeded run produces failures outside the inherited seed, and the stress guard stops the cracked runs early. These results do not establish causal rift advance or support the earlier claim of a consistent pressure-driven stress-sign pattern.

Gate 3 remains open. The next mechanics step is to couple finite-reservoir depletion, dike-tip advance, and regularized host-rock stress in one stable calculation, then test failure beyond the inherited seed. A plume comparison still requires an explicit spatial-temporal plume source and matched source-on/off cases.

## Reproduction

```bash
backend/.venv/bin/python bin/debug/test_gate3_magma_intrusion.py \
  --output analysis/issue487/gate3-magma-intrusion-review.json
```

The compact result is [`gate3-magma-intrusion-review.json`](gate3-magma-intrusion-review.json). It identifies the two matched pressure pairs and the intact-host context cases. The diagnostic imports [`prototype_gate3_finite_strain.py`](../../bin/debug/prototype_gate3_finite_strain.py), which uses the shared local mechanics in [`prototype_gate3_local_feedback.py`](../../bin/debug/prototype_gate3_local_feedback.py).

Run the independent constitutive checks with:

```bash
backend/.venv/bin/python bin/debug/verify_gate3_finite_strain.py \
  --output analysis/issue487/gate3-finite-strain-verification-review.json
```
