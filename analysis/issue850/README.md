# RFC 0018 spectral weighting review

## Scope and controls

This is a read-only comparison for issue #850. The frozen candidate and resistance
implementation at `a18603104773ddda9912c7e926f8bd868cb82554` was loaded without editing it
(source SHA-256 `1f635bd51e595f4c5aaedb95f359f05391dcb4334720a287574e7c4e617c0d13`). Only the
edge-cost-to-Laplacian-conductance mapping changed. Resistance factors, candidate score,
candidate proposal rule and breakup prior were not tuned. No world was stepped and no G2
seed was loaded.

The snapshot comparison used all 20 checksummed development snapshots in the frozen G0
manifest (source commit `63290b053de15b72280fe28fe0086ff37affdba2`): 68 eligible opportunities
and 812 candidate partitions per mapping. Each opportunity ran once under each mapping,
with identical deterministic starting-vector logic. The timing run started with no other
registered heavy workload. The compact per-candidate fingerprints, anchor costs and mapping
results are in [`spectral-weighting.json`](spectral-weighting.json).

## What the mapping means

For adjacency weights `w_ij`, the graph Laplacian is `L = D - A`, and its quadratic form is
`xᵀLx = Σ w_ij (x_i - x_j)²`. A larger adjacency weight penalizes separating those cells;
the Fiedler mode gives a low-cost way to divide the graph. RFC 0018 defines `edge_cost` from
structural resistance and edge length. With the current `w = 1 / edge_cost`, a low-resistance
edge becomes a high-weight connection, so the spectral partition tends to keep it together.
That direction conflicts with the stated goal of proposing cuts along weak, low-resistance
paths. `w = edge_cost` makes low-cost edges cheaper to cross in the spectral cut.

Three mappings were compared:

| Mapping | Conductance | Interpretation |
|---|---|---|
| Current | `1 / edge_cost` | Inverse resistance |
| Direct | `edge_cost` | Cost weighted cut |
| Uniform | `1` | Ignores resistance in the eigensolve |

## Synthetic controls

The control used 100 deterministic 40-node cycles. Each had two antipodal planted weak edges
with costs in `[0.05, 0.06)` and background costs in `[0.9, 1.1)`. The same graph was solved
with all three mappings using the balanced partition from the second normalized-Laplacian
eigenvector. A dense symmetric eigensolve avoids sparse-start-vector ambiguity on these small
controls. Each graph and mapping was solved twice; all 300 edge cuts repeated exactly.

| Mapping | Both planted weak edges crossed | Any planted weak edge crossed |
|---|---:|---:|
| Current inverse cost | 0/100 | 0/100 |
| Direct cost | 100/100 | 100/100 |
| Uniform | 8/100 | 8/100 |

This controlled result shows the direction of the mapping clearly: the current inverse
mapping avoided the weak links, while direct cost crossed them in every graph.

## Frozen G0 snapshot results

All mappings completed on the same 68 opportunities and each produced 812 candidate
partitions without diagnostic errors. Candidate partition fingerprints had zero overlap
between current and direct cost on every opportunity. Their selected partitions differed on
all 68 opportunities, and their selected anchor edge was the same on only 9/68. Direct cost
selected anchors lower in the initial edge-cost distribution (median percentile 11.5%, mean
19.3%) than current inverse cost (median 20.1%, mean 29.4%). Uniform selected anchors had
median 19.8% and mean 31.2%.

| Mapping | Base solves | Candidate refinements | Total sparse eigensolves | Runtime |
|---|---:|---:|---:|---:|
| Current inverse cost | 68 | 6,413 | 6,481 | 311.6 s |
| Direct cost | 68 | 6,004 | 6,072 | 259.6 s |
| Uniform | 68 | 812 | 880 | 50.6 s |

The current and direct runs each made one base partition solve per opportunity. The remaining
solves are candidate-specific resistance/refinement passes, up to eight per candidate. Direct
cost reduced the refinement count by 6.4% and observed diagnostic runtime by 16.7% compared
with the current mapping. Uniform weights were much faster because each candidate converged
after its first refinement; this is a different candidate method, not an equivalent
optimization. These snapshot timings do not predict full-simulation G3 runtime.

No candidate passed the frozen probe's `usable_candidate` screen under any mapping. That
historical screen includes the superseded `1.96 ×` residual-opening cutoff, so this result is
not a test of the accepted RFC 0018 `O/R/S` prior or its current G1 gate. It is reported only
to show that this geometry comparison did not create a candidate that passed the old screen.
The selected-candidate median net opening was `+0.087 cm/yr` under current inverse cost and
`−0.059 cm/yr` under direct cost. This is a development-snapshot observation, not a measured
breakup rate.

## Recommendation and uncertainty

The current `1 / edge_cost` mapping does not match the weak-path intent in the synthetic
control and produces entirely different candidate partitions on the G0 snapshots. Direct
cost has the expected weak-link direction and selects lower-ranked edge costs, but its
selected snapshot cuts have lower median opening and no candidates pass the obsolete screen.
This does not establish direct cost as the final method; it establishes that the mapping is a
material design choice that cannot be treated as a harmless implementation detail.

The production prototype and same-trace requirements in [#846](https://github.com/adubey/mantle-bloom/issues/846)
should be rebuilt only after the RFC owner settles this mapping through the open revision
tracker [#848](https://github.com/adubey/mantle-bloom/issues/848). If the mapping changes:

- Repeat G0 geometry, candidate availability, score and runtime evidence on the development
  controls because the candidate partitions and selected edges all changed here.
- Repeat G1 controlled probability and structural checks against the final candidate method.
- Re-establish G2 with a preregistered held-out sample for the revised method. The ten
  completed G2 seeds remain excluded from tuning and revised-method validation; preserve the
  unused seed registry for the owner’s explicit decision.
- Run G3 after the method and failed-rift semantics settle, against a paired production
  baseline with the required held-out seeds and runtime controls.

Limitations: this compares the frozen historical candidate probe on young development
snapshots, not a production implementation or full-length simulations. It changes only the
conductance mapping, but the changed partitions also change the candidate-specific resistance
refinements. The stored candidate fingerprints are orientation-invariant hashes, not full
geological traces. No physical resistance factors were tuned, and no held-out G2 evidence was
used.
