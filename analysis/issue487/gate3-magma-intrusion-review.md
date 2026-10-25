# Gate 3 controlled local magma-source comparison

**Date:** 2026-10-24
**Disposition:** The prescribed dike produces the predicted local stress relief and tip stress concentration, but the diagnostic does not demonstrate propagation. Gate 3 remains open.

## Question and control

Under the same extension loading and inherited damage, does a localized pressurized dike change source-flank and tip/collar stresses in the predicted directions?

The diagnostic compares intact tectonic loading with and without a dike, then repeats the pair with an identical near-tip damage seed. It uses a 400 × 130 km plane-strain strip, a 32 × 20 mesh, 4 mm/yr extension, 0.02 Myr steps, and a requested 2.5 Myr duration. The pressure is calibrated in the elastic mesh to the 2005 Afar dike's approximate 65 km length and 4 m mean opening; the simulated maximum opening is an independent check against the observed ~8 m bound. This is a prescribed constant-pressure source, not a finite-reservoir or plume-on/off experiment.

## Result

| Matched case at 0.2 Myr | Source-flank mean σxx | Tip/collar peak tensile stress | Stop / failures outside seed |
|---|---:|---:|---|
| Intact, no dike | 61.93 MPa | 61.93 MPa | Continues; 0 |
| Intact, dike | 34.47 MPa | 325.70 MPa | 5× strength guard; 0 |
| Near-tip seed, no dike | 60.83 MPa | 69.73 MPa | Continues; 0 |
| Near-tip seed, dike | 35.78 MPa | 339.13 MPa | 5× strength guard; 0 |

The no-dike cases complete 2.5 Myr. In the seeded pair, both cases fail the same 16 inherited seed elements; neither fails an element outside the seed. The dike-on cases stop at 0.2 Myr when the pointwise peak tensile stress crosses the numerical guard, before outside-seed failure can test advance.

The source pressure is 3.961 MPa. The elastic dike solution gives 4.000 m mean and 5.762 m maximum opening. Its equilibrium residual is `4.54e-13`; source/support force and moment residuals are `2.60e-13` and `5.19e-13`; pressure-work and stored-energy mismatch is `8.05e-7`. The minimum deformation determinant in the seeded dike-on run is positive (`0.9983`).

## Interpretation and next gate work

The preregistered stress-sign prediction is supported: opening relieves extension in the source flanks and transfers a strong tensile concentration to the tips. The pointwise tip stress guard makes the propagation result inconclusive; it is not evidence of geological rift advance. Constant pressure also means this run has no finite-reservoir magma-volume closure. The source is prescribed independently of plume attribution, so it does not test RFC 0007 plume forcing.

Gate 3 remains open. The next useful mechanics step is to couple finite-reservoir depletion, dike-tip advance and regularized host-rock stress in a single stable calculation, then test whether it causes failure beyond the inherited seed. A separate plume comparison requires an explicit spatial-temporal plume source and should change only that source between matched on/off cases.

## Reproduction

```bash
backend/.venv/bin/python bin/debug/test_gate3_magma_intrusion.py \
  --output analysis/issue487/gate3-magma-intrusion-review.json
```

The result is [`gate3-magma-intrusion-review.json`](gate3-magma-intrusion-review.json). The diagnostic imports [`prototype_gate3_finite_strain.py`](../../bin/debug/prototype_gate3_finite_strain.py), which uses the shared local mechanics in [`prototype_gate3_local_feedback.py`](../../bin/debug/prototype_gate3_local_feedback.py). The mesh refactor reuses the base element gradients because the duplicated crack-face nodes retain identical reference coordinates; the controlled rerun reproduces the previous result values to floating-point roundoff.
