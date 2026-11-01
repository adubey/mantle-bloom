# Issue #731: exported lake overflow ownership

## Replay preregistration

Compare the four 100 kyr steps in RFC 0016 windows A, B and C on the same recreated
checkpoints. The control is the current #729 river-ledger implementation before this fix; the
treatment is that implementation with #731's lake-owner and terminal-overflow routing changes.
The three input checkpoint SHA-256 values are A `9882568e681a0f9f207e3302453b3a92f18865d5640d7baa8d9d8eb09ab4633e`,
B `759d83107ae87d34d825b4836addc922d414089774f90739ea00e6a9c9dff5b5`, and C
`38cc6ff2054cde2db3a65ac0a3e08a09e854175735512b9c8c40ba7a76c6190c`. These are the event-
validated replacements documented in RFC 0016's G0 report. Each run uses the matching
`target-anchor-window-<a|b|c>.json`, disables the ice-age cycle, and advances four 100 kyr steps.

Primary metric: `stranded_overflow_m3 / (rain_m3 + meltwater_m3)` for each step, with a
decision threshold of at most 0.01% in every window. The rain-only land-ledger residual must
remain at or below 0.1% of input. The fix is expected to transfer some formerly stranded water
to a parent lake or the ocean, so ocean export and internal transfer may increase. Target-lake
volume and event ΔV may move substantially; record the change without tuning the fix to preserve
the control event. Also report lake volume, target event ΔV, ocean export, stranded water,
internal transfer, flush passes, step time, and peak RSS. Flow-cycle count and whether spill edges
still close cycles are separate diagnostic outcomes; a remaining cycle is acceptable only when
its water has an explicit owner or settlement.

Controls are the same checkpoint bytes, anchors, code/runtime dependencies, forcing
(`ice_age_period_years = 0`), step count, step length, central Raymond velocity and one-year
residence cap. Run baseline and treatment arms sequentially, after checking compute budget for
each. Treat runtime as indicative if another workload overlaps. The prior G0 control measured
5.30–6.10% stranded overflow per step in A, 1.11–5.25% in B, and 0.000003–0.34% in C; the
expected treatment direction is a reduction below threshold. No replicate is preregistered:
these are large deterministic routing changes, and event movement is the decision signal rather
than a small effect near replay noise.

## Replay results

The baseline is PR #729 at `b003a3f` before this issue's routing changes. The treatment is that
PR snapshot plus the bounded #731 code. The exact temporary source snapshots and per-step data
are identified in [`replay-comparison.json`](replay-comparison.json); full outputs stay in the
private capture directory. All six runs were sequential and compute-budget reported no competing
heavy run. The first baseline A attempt completed one step but could not write its output in the
restricted file context; it was discarded and rerun successfully in the authorized capture
context.

| Window | Max stranded overflow, control → treatment | Max land residual, treatment | Target event ΔV, control → treatment (km³) | End target volume change (km³) | Final flow cycles / spill-edge cycles, control → treatment | Median step time, control → treatment (s) | Peak RSS, treatment (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 6.096% → 0% | 0.00000468% | +495,623.987 → +495,623.987 | +112.244 | 63/63 → 63/63 | 7.96 → 7.97 | 5.45 |
| B | 5.251% → 0% | 0.004348% | +131,625.097 → +131,624.955 | −26,620.195 | 47/47 → 51/51 | 7.66 → 7.83 | 5.70 |
| C | 0.342% → 0% | 0.000000018% | −261,586.996 → −261,586.996 | −0.069 | 14/14 → 14/14 | 10.13 → 10.26 | 6.90 |

Every treatment step meets the 0.01% stranded-overflow threshold and the 0.1% land-ledger
residual threshold. The registered event changes are negligible: A and C match to the reported
precision, and B differs by 0.142 km³ from the recreated control event. The events do not need
re-baselining. End-state lake volumes still move: global lake volume changes by −4,512 km³ in A,
+52,477 km³ in B and +18,863 km³ in C. The endpoint target volume also shifts in B, so retain the
comparison in the RFC record even though its registered event ΔV is preserved.

Directed cycles remain in every window. At each final step, every counted cycle includes at
least one edge matching `spill_target`; all A and C cycles and all but four B treatment cycles
meet that check. This remains a legitimate graph state. Exported lake water in those cycles is
credited to the ocean when it has no receiving lake body; `route_downstream` reports cycle
settlement separately. Each step uses two river reservoir sweeps: runoff, then exported lake
overflow. Aggregate internal transfers change by +1.8% in A, +3.6% in B and −0.02% in C. Ocean
export rises as expected when overflow that formerly settled on land reaches the ocean.

The first treatment experiment extended the zero-duration flush to twice the leaf count and
repeated the final owner handoff. It reduced stranded overflow but increased gross internal
transfers by roughly an order of magnitude and made A/B about three times slower. Those runs
were superseded. The bounded implementation allows one final owner handoff plus one follow-through;
it preserves zero stranded overflow while bringing transfer totals and runtime close to control.
Both experiment sets remain in private captures and their checksums.

The stress test ran a generated world for 192 Myr. A first 24-step run used 8 Myr per step and
found one 0.103% land residual, just above the 0.1% tolerance; stranded overflow remained zero.
The retained 48-step run used 4 Myr per step over the same horizon and passed both limits. This
also exercises remapping and repeated routing on a mature world; a 4 Myr stress timestep is still
coarser than the 100 kyr replay timestep.

## Remaining acceptance work

PR #729 merged as `cd93d339` while this issue was being worked. The A/B/C replays used temporary
integration snapshots of that PR plus #731; the final issue branch must now merge refreshed
`origin/main`, wire the river-router cycle-settlement field and direct ledger invariant test into
the merged API, and add the mature-world stress test. Re-run the required pre-push checks and
update code provenance against that merged branch before opening the PR. The temporary replay
integration is evidence for the behavior, not final branch CI.
