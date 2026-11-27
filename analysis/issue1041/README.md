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

## Registered Benchmark Results

The comparison ran sequentially on the same 16 GB Apple Silicon host, with no other tracked
heavy run active. `profile_phases.py` used one warmup and three measured steps per run. Young
profiles use seed 0 at density 4/4/2; mature profiles load the named 67.9 Myr save. Both used
diagnostic wind. Peak RSS is `/usr/bin/time -l` maximum resident set size.

| Case | Baseline total step | Candidate total step | Baseline routing + order | Candidate routing + order | Peak RSS baseline → candidate |
|---|---:|---:|---:|---:|---:|
| Young, one registered run | 14.263 s | 14.276 s | 1.711 s | 1.717 s | 1,051 → 941 MB |
| Mature, run 1 | 16.757 s | 15.966 s | 2.031 s | 1.833 s | 979 → 1,007 MB |
| Mature, replicate | 13.925 s | 13.984 s | 1.463 s | 1.499 s | 1,076 → 1,045 MB |
| Mature, mean of two runs | 15.341 s | 14.975 s | 1.747 s | 1.666 s | 1,027 → 1,026 MB mean |

The mature routing/order phase mean fell **4.6%**, below the preregistered 10% threshold for
calling a subsystem speedup. The young routing/order phase was unchanged (+0.4%). Mature total
step time fell 2.4% across two runs, but individual runs varied by about 2.8 seconds, so that
difference is inconclusive. Peak RSS was effectively unchanged across the mature runs; the
young single-run reduction is not replicated and should not be treated as an established memory
saving.

A `cProfile` pass measured `_compute_flow_direction` and its Python bridge at about 20.4 ms →
4.6 ms per call in the young case and 22.6 ms → 5.5 ms per call in the mature case (about
4.0–4.2x faster). This kernel is called three times per measured step. Its local improvement
is real in the profile, but the full routing phase remains dominated by Python routing and DAG
ordering, outside this PR's current port.

**Decision:** the isolated flow-direction kernel is faster, but the registered whole-routing
metric did not reach its 10% threshold. Treat the current change as a small kernel-level
improvement; do not claim a material step-time or peak-RSS win. The next routing slice should
port `_routing_order` and `route_downstream` to address the remaining measured cost.

Reproduce with `bin/debug/profile_phases.py` (`--seed 0` for young; `--load <local-save>` for
mature), `--scenarios step --steps 3 --warmup 1 --wind-model diagnostic --cprofile <scratch.prof>
--json <scratch.json>`, wrapped in `/usr/bin/time -l`. Control commit: `334e1eeb`; candidate
commit: `934dbd94`. The control worktree predates current `origin/main` (`d332eadd`) only by
non-simulation RFC closeout commits. Scratch profiles and full run output are not checked in.
