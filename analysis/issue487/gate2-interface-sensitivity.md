# Gate 2 interface-force ownership and coupling sensitivity

## Question

A source-only shell solve applies the source plate's current nodal force field on a free-edge mesh. Adding a neighboring plate changes the mesh and applies that plate's loads too. The joined mesh shares displacement degrees of freedom across the pair, which is an idealized bonded/no-slip interface. The earlier pair solve also applied the torque engine's contact-drag loads on both plates. This audit asks whether those pair-specific forces are reciprocal and how much the interface assumption changes the source cut-band stress.

## Force ownership

The per-node fields retain slab pull, slab drag, ridge push, basal drag, and contact drag. In a pair solve, slab pull, slab drag, ridge push, basal drag, and contact drag from other neighbors are external to the selected pair. The contact-drag force between the selected two plates is an interaction internal to the pair: it should be represented once as an interface exchange.

The diagnostic now extracts contact drag attributed specifically to the selected source-neighbor pair using each node's nearest-neighbor plate id. It reports the sum of both sides' contact-force vectors and moments relative to the sum of nodal magnitudes. A zero ratio can mean either reciprocal forces or no pair contact drag; the reported L1 force total distinguishes those cases.

Two coupled cases are compared:

1. **Current pair contact loads:** the shared-node shell mesh with both plates' existing force fields, including their pair-specific contact drag.
2. **Tied interface without pair contact loads:** the same mesh and other loads, but remove only mutual contact drag. The tied displacement constraint supplies the interface reaction. Loads from other neighbors stay applied.

Neither is a calibrated tectonic contact law. The model's contact drag is a velocity-dependent force in the plate-motion torque calculation; it is not a shear stiffness or a frictional slip law for the static thin-shell displacement solve. The tied case is an idealized bonded endmember. It does not imply that real plate boundaries are bonded, nor that its stress is an upper bound.

## Matched event sample

One connected-source split and one connected-source failed-rift event were replayed for each Gate 1 seed. Each event key (seed, event index, kind, world time and source plate id) was checked against the saved Gate 1 sensitivity sidecar. All six runs had 100% matched source cut-band coverage. Full replay reports are kept in /private/tmp; this compact file contains the recorded metrics. Run the existing replay script with the seed/time/plate/kind rows below (the script docstring gives the full environment invocation); use /private/tmp/gate2-pair-s<seed>-<fail|split>.json as the output path.

The stress-change metric is the mean absolute nodewise difference divided by the mean isolated source cut-band N_ext. This shared band-scale denominator avoids dividing by a zero-tension value at an individual node. The earlier seed-7 reports instead divided by each isolated node value, including zero values; their 26.2% and 51.2% relative-change summaries are invalid and superseded. The stress means themselves remain valid.

| Seed / event | Source + neighbor | Mutual pair contact force L1 | Force closure | Torque closure | Current-contact mean absolute change | Tied/no-pair-contact change | Joined residual/load |
|---|---:|---:|---:|---:|---:|---:|---:|
| 7 failed rift, 2 Myr, plate 1 | 19 | 2.23e19 N | 21.4% | 21.6% | 97.0% | 118.7% | 0.010% |
| 7 split, 4 Myr, plate 3 | 14 | 4.70e18 N | 2.0% | 4.6% | 44.8% | 45.4% | 0.019% |
| 11 failed rift, 2 Myr, plate 1 | 3 | 0 N | n/a | n/a | 22.3% | 22.3% | 0.131% |
| 11 split, 16 Myr, plate 1 | 3 | 0 N | n/a | n/a | 12.8% | 12.8% | 0.034% |
| 2024 failed rift, 2 Myr, plate 8 | 7 | 7.43e18 N | 9.9% | 10.2% | 6.8% | 5.8% | 0.029% |
| 2024 split, 4 Myr, plate 8 | 7 | 8.09e18 N | 14.9% | 15.0% | 5.3% | 6.9% | 0.018% |

The paired contact-drag force is absent in both sampled seed-11 events, so those cases cannot validate action-reaction closure. For the four cases with nonzero mutual contact drag, the force closure residual ranges from 2.0% to 21.4%; torque closure ranges from 4.6% to 21.6%. Independent nearest-node/front-length discretizations therefore do not guarantee reciprocal pair forces and moments. This is a force ownership/discretization problem to solve before treating current contact loads as a physically closed pair interface.

Coupled N_ext is sensitive to adding the neighbor, with the corrected band-scale difference from 5.3% to 97.0%. Removing pair contact drag changes the result by 0% to 71.8% of the current-contact joined mean, and does not move stress monotonically toward or away from the isolated result. This is sensitivity evidence; it is not a physical prediction or a chosen law.

## Gate 2 disposition

The repeatable split/failed-rift comparison is now complete for the first connected-source event of each kind in each seed. It broadens the two earlier seed-7 split probes, but it remains a six-event short-run sample and does not validate the loads over the full 1 Gyr candidate population. The idealized tied interface and current velocity-based drag do not define a finite-slip traction law. The contact force imbalance and large event-to-event sensitivity keep Gate 2 open.

A later coupled-load model should discretize the shared contact once, assign an equal-and-opposite traction to both sides, and make opening, closing and tangential slip assumptions explicit. A finite frictional or viscous interface also needs a physically or empirically constrained slip-rate/state law and boundary geometry; this diagnostic does not choose its coefficients. Before selecting that law, separately audit the plate-wide rigid-rotation projection, since it removes a load mode already handled by the motion solver and its pair-level effect may contaminate local stress comparisons.
