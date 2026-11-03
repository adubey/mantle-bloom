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
