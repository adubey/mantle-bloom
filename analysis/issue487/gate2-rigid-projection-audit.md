# Gate 2 pair-level rigid-rotation projection audit

Date: 2026-10-22. This audit isolates the projection choice on six matched connected-source split and failed-rift events (seeds 7, 11 and 2024). Each comparison uses the same joint mesh, nodal force field, thicknesses, stiffness and source cut-band nodes.

## What the two projections mean

The current solver removes the area-weighted load component associated with the joint mesh’s three global rigid rotations. The joint mesh is a single continuous membrane, so that is the mathematically consistent nullspace projection for its bonded kinematics.

The audit variant independently removes the source and neighbor groups’ rigid-rotation loads, then solves the same continuous joint mesh. This tests sensitivity to letting each plate group shed its own torque before coupling. It is not a free-slip solve: the displacement mesh remains continuous, and there is no contact law or interface reaction equation. Neither result establishes the physical interface behavior.

## Results

| Seed / event | Source + neighbor | Pair-wide load removed | Separate-group load removed | Mean absolute cut-band ΔN_ext / pair-wide mean | p95 / pair-wide mean | Residual pair-wide / separate (%) |
|---|---:|---:|---:|---:|---:|---:|
| 7 / failed_rift @ 2 Myr | 1 + 19 | 2.7% | 19.3% | 95.1% | 715.4% | 0.010% / 0.007% |
| 7 / split @ 4 Myr | 3 + 14 | 15.4% | 27.8% | 36.6% | 71.6% | 0.019% / 0.015% |
| 11 / failed_rift @ 2 Myr | 1 + 3 | 6.9% | 15.5% | 22.7% | 38.6% | 0.131% / 0.039% |
| 11 / split @ 16 Myr | 1 + 3 | 1.3% | 2.8% | 3.0% | 6.4% | 0.034% / 0.027% |
| 2024 / failed_rift @ 2 Myr | 8 + 7 | 10.5% | 12.2% | 4.0% | 14.0% | 0.029% / 0.025% |
| 2024 / split @ 4 Myr | 8 + 7 | 1.2% | 2.5% | 3.1% | 8.9% | 0.018% / 0.018% |

All six events had 100% source cut-band coverage. Both projection variants stopped with code 1 for every case, and all residuals were below 0.132%; the stress differences are not explained by solver non-convergence. Mean absolute cut-band stress differences range from 3.0% to 95.1% of the pair-wide cut-band mean. Seed 7 is especially sensitive: 95.1% for the 2 Myr failed rift and 36.6% for the 4 Myr split. The separate-group projection removes 2.5–27.8% of raw pair-load norm across cases, compared with 1.2–15.4% under the joint projection.

The comparison uses mean absolute difference because per-node relative changes are unstable where tensile resultant is zero. The p95 column is likewise normalized by the pair-wide cut-band mean, not by each node’s baseline.

## Audit conclusion

The audit confirms the implementation matches its stated joint-mesh assumption, but that assumption materially affects the coupled stress result in some events. The separate-group variant is not a validated alternative because it removes plate-group torques while retaining bonded joint kinematics. The present pair solve also applies each side’s current pair-specific contact drag even though the six-event audit found those forces are not consistently action-reaction balanced. A physically interpretable coupled stress requires a stated interface model and a matching formulation for kinematics, traction ownership and torque balance.

Gate 2 remains open. Do not select a finite-slip law or use the coupled values as projection-independent physical predictions based on this sensitivity audit. See the machine-readable event metrics in [gate2-rigid-projection-audit.json](gate2-rigid-projection-audit.json) and the replay implementation in [replay_gate2_coupled_event.py](../../bin/debug/replay_gate2_coupled_event.py).
