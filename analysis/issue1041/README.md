# Issue #1041 — Hydrology Routing Benchmark

## Pre-registration

**Question:** Does the Rust flow-direction kernel improve the hydrology routing workload without changing flow targets or spill masks?

- **Control:** detached baseline worktree at `334e1eeb` (the fetched `origin/main` used when this task-specific baseline was created; current `origin/main` is `d332eadd`, whose intervening commits do not change simulation code).
- **Candidate:** `codex/issue-1041-biome-kernels` at the rebased implementation commit recorded with the results below.
- **Young case:** seed 0, density 4/4/2, diagnostic wind, one warmup and three measured steps.
- **Mature case:** `mantle-bloom-seed47235302-67900000y.bloom` (67.9 Myr), diagnostic wind, one warmup and three measured steps.
- **Primary metrics:** mean step wall time and the combined `Ordering the flow network` plus `Routing flow downstream` phase times; peak RSS is secondary. Compare each case on the same machine with identical inputs, running arms sequentially.
- **Expected direction:** lower flow-routing time and no material increase in total step time. Treat a total-step change below 2% as inconclusive given three-step sampling and simulation variance. A routing-phase decrease of at least 10% is the threshold for calling this slice a measurable speedup.
- **Parity invariant:** exact flow-target arrays and `should_spill` masks; a saved-world replay must preserve the existing output/ledger tolerances before claiming end-to-end parity.

The earlier exploratory young baseline run on `334e1eeb` is retained in the local scratch outputs; the registered baseline and candidate runs below are the comparison used for claims.
