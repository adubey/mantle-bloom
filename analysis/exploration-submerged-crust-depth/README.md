# Exploration: fixed submerged continental depth and seed 528019447

## Preregistration

**Question.** Does setting submerged continental crust to a fixed target of -666 m change the deepest disconnected below-sea basin in the cratons-first, synthetic-relief world generated with seed 528019447?

**Hypothesis.** The deep basin is a submerged hole in continental crust that `synthetic_relief.drain_land` preserves. If so, setting each submerged continental target to -666 m should raise the basin floor from the observed baseline near -1,422 m to approximately -666 m. The below-sea connected-component footprint should remain nearly unchanged because those nodes remain below sea level. If the floor is instead carried by oceanic crust, its deepest point should remain below -666 m.

**Controls.** Use the same source commit (`6abd9947e857a738e8af8635193fdc91b43a1aba`), seed, and Generate dialog settings in both arms: continental fraction 0.70, land fraction 0.29, cratons-first on, synthetic relief on, node density 4, climate density 4, fluid density 2. Run control and perturbed arms in separate fresh generations, sequentially. Do not step the world.

**Perturbation.** In an isolated scratch copy only, set `z_sea = -666.0` for `cloud.continental & ~land` in `synthetic_relief._target_elevation`; leave oceanic cells and land targets unchanged.

**Primary metric.** For each arm, report each disconnected below-sea component's area and minimum elevation; identify the component with greatest overlap against the baseline largest component. Also report the deepest component's continental/oceanic node counts and the total land area / land fraction.

**Prediction and rough size.** If the deepest baseline nodes are submerged continental crust, the deepest elevation in the matching component should rise by roughly 756 m (from -1,422 m to near -666 m), within 10 m of the imposed target, while its footprint changes by less than 1%. Land fraction should be unchanged. If its floor stays below -676 m or the component is mainly oceanic, this target does not control that floor.

**Validity checks.** Sea level stays at 0 m; land fraction and land mask stay unchanged; oceanic crust elevations stay unchanged; no elapsed time or tectonic step is advanced. An arm violating these checks is invalid for the hypothesis.

**Baseline evidence from the prior ten-world run.** Seed 528019447 at these settings had two disconnected below-sea components totaling 2,278,782 km²; its lowest below-sea node was -1,422.22 m. This earlier diagnostic did not retain per-node crust-type composition. The control arm will be regenerated before comparing.

## Results

Pending preregistration publication and paired run.

## Probe outcome

The paired runs support the hypothesis.

| Arm | Largest matching component area (km²) | Nodes | Continental / oceanic nodes | Floor (m) | Total disconnected below-sea area (km²) | Area land fraction |
|---|---:|---:|---:|---:|---:|---:|
| Control | 2,270,390.76 | 638 | 638 / 0 | -1,422.22 | 2,278,782.05 | 0.28541105 |
| Submerged continental target = -666 m | 2,270,390.76 | 638 | 638 / 0 | -666.00 | 2,278,782.05 | 0.28541105 |

The largest component's floor rose by 756.22 m and landed within 0.001 m of the imposed target. Its node membership and accounting area were unchanged. The second, 2-node component also retained its 8,391.30 km² area but its floor moved from -40.05 m to -666.00 m. Land area, land-node fraction, sea level and both below-sea components' topology were unchanged. No oceanic nodes were present in either closed component.

**Interpretation.** For this seed, the 1,422 m floor is carried entirely by submerged continental crust. The synthetic relief's variable submerged-crust target creates that depth, while `drain_land` preserves the enclosed submerged hole rather than filling it. The perturbation demonstrates that replacing the submerged continental target applies directly to the basin.

**Limits.** This is one seed at one resolution and one imposed target. It is a causal mechanism probe, not evidence that -666 m is a suitable production rule or that the same result generalizes across seeds. The source change existed only in `/private/tmp/mb-depth-666-probe`; no model implementation was changed in the shared checkout.

## Classification clarification

The basin's 638 nodes were all continental and below sea level, but they were not classified as land by the synthetic relief land mask, nor as open ocean by `connected_ocean_mask`. In `drain_land`, they are `hole = ~land & ~ocean`: an enclosed submerged hole. The drainage pass intentionally excludes sea nodes from its land-raising return, so the synthetic submerged-crust target survives there. The `-666 m` arm replaced the target for every submerged continental node, not only this basin; its exact effect on the basin confirms which field controls the floor, not that a basin-specific cap is a suitable fix.

The temporary `-666 m` assignment has been removed from the scratch source copy. The shared model source was never changed.

## Preregistered follow-up: Caspian-scale isolated-hole bounds

**Question.** Can the oversized isolated submerged continental hole in seed 528019447 be bounded to roughly a Caspian-sized footprint and floor without changing the rest of the generated coastline or basin connectivity?

**Working bounds.** Use 371,000 km² as the maximum below-sea footprint of any one isolated hole, matching the Caspian Sea's approximate surface area. Use -1,050 m as a rounded lower floor bound relative to world sea level: the Caspian surface is about 28 m below global sea level and its maximum water depth is about 1,025 m. These are exploratory Earth-scale ceilings, not calibrated Mantle-Bloom parameters.

**Control.** The already-published seed 528019447 control at the UI defaults: largest isolated hole 2,270,390.76 km², 638 continental nodes, floor -1,422.22 m, land area fraction 0.285411. Sea level 0 m. The control commit is `6abd9947e857a738e8af8635193fdc91b43a1aba`.

**Perturbation.** In a disposable source copy, act only on below-sea continental components disconnected from the main ocean. For a component above the area bound, keep a connected compact core around its deepest node up to the area limit and raise excess nodes to `MIN_LAND_M`; clamp the retained component's floor to -1,050 m. Leave land targets, open-ocean-connected nodes and oceanic crust unchanged. This is a test of the explicit size/depth policy, not a proposed production implementation.

**Primary metrics.** Largest and total isolated-hole area, component count, deepest floor, land area fraction, land-node fraction, and unchanged area/elevation outside the selected hole.

**Prediction.** For seed 528019447, the largest hole should shrink from 2.27 million km² to at most 371,000 km² (about 84% smaller), with its floor no lower than -1,050 m. Sea level and open-ocean connectivity should remain unchanged. The land area fraction should increase by approximately the excess raised area divided by planetary surface area; report the measured change rather than treating it as zero.

**Validity checks.** Sea level remains 0 m; no tectonic step advances; crustal backing reproduces the edited targets; oceanic elevations and all originally dry land targets remain unchanged. Report any failure, especially if selecting the compact core fragments it into multiple disconnected holes or causes other holes to exceed the cap.

**Results.** Pending the paired scratch run.
