# RFC 0018 implementation review: runtime and failed-rift coverage

## Finding

The existing G2 graph method shows a substantial runtime increase that appears inconsistent
with G3's 2% median-runtime limit. The first ten paired G2 runs have median wall times of
505.2 seconds for the production baseline and 8,363.5 seconds for the graph arm, a 16.25×
ratio. Individual paired ratios range from 14.07× to 23.58×. Exact per-seed measurements are
in [`g2-runtime-review.json`](g2-runtime-review.json).

This is a strong warning, not a definitive G3 failure. G2's report does not record sufficient
runtime-isolation metadata to treat the wall-time ratio as a controlled performance
comparison. G3 specifically requires paired production runtime, and should not be marked
failed from these values alone.

## Production-semantics gap

The G2 callback uses the frozen mass-balanced candidate generator and `O/R/S` prior, but it
returns immediately when the event draw does not select breakup. It does not apply a failed
rift or record a graph-path craton cut on that branch. Production issue #846 requires a failed
breakup attempt to thin the selected graph path and use the same trace for craton diagnostics.
Therefore the G2 callback does not fully exercise production semantics. Applying failed-rift
thinning may alter subsequent world evolution and the actual runtime.

The G2 result still establishes its previously defined non-worsening outcome for the
experimental callback on the ten evaluated seed pairs. It does not establish the G3 runtime
limit or production behavior for failed rifts. The G2 seeds remain held out and must not be
used to tune an amended method.

## Reproduction and provenance

- G2 report and seed registry: `analysis/issue797/README.md`,
  `analysis/issue797/g2-summary.json`, and `analysis/issue797/g2-seed-preregistration.json`.
- Protocol: `RFC0018-G2-experimental-v1`; 200 Myr, 1 Myr per step; first ten seeds of the
  frozen 20-seed registry.
- Frozen candidate method source commit: `5da20be8cbf7760960a2557bcba915c8fb828292`.
- G2 runner SHA-256: `ebb0978027631ee8e54b26300497a85ef7da1b234d7cb33f930c0c06240da54c`.
- Seed registry SHA-256: `214623b962984be9e95c6ced1f87b470b105a6d628f96cdb7cfae5d5844eba2f`.
- Runtime summary calculations use the paired per-seed `wall_s` values in the compact G2
  summary; the ratios are `graph_wall_s / baseline_wall_s`.

## Revision question

Determine whether the runtime limit, candidate-search implementation, and failed-rift
semantics can be made mutually consistent without changing the frozen G0/G1/G2 model. If they
cannot, propose revised runtime criteria or a revised method and identify which diagnosis
gates and held-out evidence must be repeated. Do not run on or tune against the G2 seeds.
