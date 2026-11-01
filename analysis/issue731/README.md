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

Run outputs and full saves stay in the private capture directory; only compact summaries and
the comparison JSON will be checked in after both arms finish.
