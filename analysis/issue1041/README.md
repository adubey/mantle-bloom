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

## Glacier-march glue and `ice_step.rs` profiling (item 5)

**Question:** Where does the glacier march's remaining time go — Python glue around
`mantle_bloom_rust.step_ice`, the sea-level re-solve, or inside the Rust kernel — and what
exact-parity optimization does that support? Principle #2 (keep the planet moving at a
useful pace); parity invariant is bit-identical march output (the existing Rust-vs-Python
`test_glacier_rust.py` / `test_ice_sheet.py` parity tests are the check).

**Method:** `bin/debug/profile_phases.py --seed 0 --scenarios step --steps 3 --warmup 1
--wind-model diagnostic` (young Density-4 world, 124,046 nodes), plus per-section Rust
timing gated on `MANTLE_BLOOM_ICE_PROFILE=1` (new in this change; prints setup / components /
signed-distance / domain-build / solve+scatter per `step_ice` call to stderr). Baseline arm
is `origin/main` (`069d1d38`) measured in the same worktree by stashing the change; arms ran
sequentially on a shared 16 GB host, so total-step wall times are noisy and phase times are
the primary metric.

**Baseline findings (per step, mean of 3):**

| Phase | Baseline | Notes |
|---|---:|---|
| `Solving march ice` (238 substeps) | 7.56 s | of which Rust `step_ice` 7.40 s (28 ms/call) |
| `Re-solving march sea level` (237 re-solves) | 1.37 s | Rust `connected_ocean_mask` 0.89 s, Rust sea-level solve 0.37 s, Python glue the rest |
| Inside one `step_ice` call (~28 ms) | — | signed distance 8.8 ms, domain build 8.3 ms, plastic solve+scatter 7.9 ms, setup 1.5 ms, components/seed scan 1.2 ms |

Domains average ~93k nodes over ~2 plastic bodies per call, rebuilt from scratch every
substep. A cross-substep topology cache on `MarchGraph` was tried and reverted: domain
membership churns between substeps (27% hit rate over 927 body builds), and the snapshot
compares plus cache-entry cloning cost more than the gathers they saved (domain build
8.3 → 10.2 ms/call). The sea-level re-solve count and the Dijkstra reach are `MarchCost`
owner decisions (production `trend-1k` re-solves every substep); not changed here.

**Change:** chunked parallel gathers in `plastic::build_local_domain` (one rayon task per
2,048-node block, concatenated in node order) and a chunked parallel `touches_out` seed scan
in `step_ice`, both bit-identical to the serial builds, plus the gated
`MANTLE_BLOOM_ICE_PROFILE` section timing.

**Results (same profile, candidate vs baseline):**

| Phase | Baseline | Candidate | Delta |
|---|---:|---:|---:|
| `Solving march ice` | 7.56 s | 7.24 s | −4.2% |
| `Marching glaciers` (Rust `step_ice`) | 7.40 s | 7.09 s | −4.2% |
| `Re-solving march sea level` | 1.37 s | 1.33 s | — (noise) |
| Mean total step | 13.43 s | 12.86 s | noisy, shared host |

Per-call section means: components/seed scan 1.21 → 0.80 ms (parallel `touches_out`);
domain build unchanged within noise (the gathers are memory-bandwidth bound, so chunking
alone does not move them). Parity: `cargo test` 25 passed; `test_glacier_rust.py` + `test_ice_sheet.py` 56 passed,
including the Rust-vs-Python step-ice parity tests (the march calls `step_ice` repeatedly
on one graph, so the parallel gathers are exercised across substeps).

**Next levers, in measured order:** (1) the signed-distance Dijkstra (8.8 ms/call) — needs a
warm-start/bounded variant, which changes results and is an owner/`MarchCost` decision;
(2) domain-build gathers (8.3 ms/call) — needs partial topology reuse keyed on the stable
subset of a body's nodes, not whole-domain identity; (3) the per-substep
`connected_ocean_mask` (3.7 ms/call) inside the sea-level re-solve.
