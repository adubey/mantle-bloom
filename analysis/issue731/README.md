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

The final A/B/C comparison uses the four preregistered 100 kyr steps and checkpoint hashes
above. The control backend is detached `origin/main` at `cd93d339` (after PR #729); the treatment
uses the issue branch model code committed in `bbe59e16`. All six arms ran sequentially after a
compute-budget check, with no competing heavy run. Per-step measures and backend provenance are
in [`replay-comparison.json`](replay-comparison.json); full outputs remain in the private capture
folder.

An initial refresh accidentally resolved the replay runner's default backend to the task checkout
for both arms. Those control attempts were discarded. The controls were rerun from the detached
`origin/main` worktree; their output records confirm `cd93d339`. Treatment output records were from
the issue branch source tree, which matches the committed model code in `bbe59e16`.
The corrected command, with local paths replaced, was:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/replay_rfc0016_river_water.py \
  <checkpoint.mbworld> --backend <arm-worktree>/backend \
  --anchors analysis/rfc0016/target-anchor-window-<a|b|c>.json \
  --myr 0.4 --step-years 100000 --checkpoint-myr 0.1 --out <scratch-output>
```

| Window | Max stranded overflow, control → treatment | Max land residual, control → treatment | Target event ΔV, control → treatment (km³) | End target volume change (km³) | End global lake volume change (km³) | Final flow cycles / spill-edge cycles, control → treatment | Median step time, control → treatment (s) | Peak RSS, control → treatment (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 6.096% → 0% | 0.00000468% → 0.0000000000103% | +495,623.987 → +495,623.987 | +124.815 | −4,198.023 | 63/63 → 63/63 | 7.95 → 8.09 | 5.37 → 5.32 |
| B | 5.251% → 0% | 0.004373% → 0.00000000000510% | +131,625.097 → +131,624.955 | +20.384 | +79,994.726 | 47/47 → 49/49 | 7.70 → 7.93 | 6.16 → 5.15 |
| C | 0.342% → 0% | 0.0000000177% → 0.0000000000109% | −261,586.996 → −261,586.996 | ≈0 | +0.013 | 14/14 → 14/14 | 10.11 → 10.66 | 6.85 → 6.51 |

Every treatment step is below both thresholds: stranded overflow is zero, and the remaining land
residual is well below 0.1% of input. Routed external supply that ends at a land sink without an
active lake body is reported as `land_sink_settlement_m3`. It is a named terminal sink, not lake
storage or ocean export. This category explains the largest prior stress-run residual (0.1327% at
step 36); after naming that destination, the 48-step, 192 Myr mature-world stress test passes both
the residual and stranded-overflow limits. In the A/B/C runs, maximum treatment land-sink
settlement is 0.00000468% of input in A, 0.004348% in B and 0.0000000177% in C.

The registered events remain present. A and C match the current `origin/main` control to the
reported precision. B moves by 0.142 km³ (about 0.00011% of its event ΔV), while its event remains
intact; retain this small change in the RFC 0016 G1 decision record. End target volumes move by +124.815
km³ in A and +20.384 km³ in B. Global lake volume changes by −4,198 km³ in A, +79,995 km³ in B
and +0.013 km³ in C. Ocean export rises by 17.7%, 9.2% and 0.27%, respectively; aggregate
internal transfers rise by 1.8%, 3.6% and 0.055%. Each step uses two river reservoir sweeps: the
runoff pass and exported-overflow pass. Step times are indicative and were measured sequentially;
there was no competing replay workload.

Directed flow cycles remain in all windows. Every counted final cycle includes an edge matching
`spill_target`, and the cycle's exported overflow is credited to the ocean when no receiving lake
body owns it. These cycles remain a legitimate graph state because their water has a named
settlement.


## Review-fix validation preregistration

The review fixes keep global ocean credit separate from spatial river-mouth discharge and
identify terminal land-sink inputs from the initial active-body membership map. They route distinct capped-pair handoffs through their owners, escalating returning surplus
to an untried ancestor. Deferred transfers use normal flush passes; an unresolved owner transfer
fails explicitly instead of exporting through its parent. Repeat the treatment on the same A/B/C
checkpoint hashes, forcing, horizon and thresholds above; the unchanged control remains
`cd93d339`. Expected direction: zero stranded overflow, unchanged global ocean credit, and no
spatial mouth at an arbitrary ocean node. Lost input at an owned lake must remain in the land
residual. Record any event, lake-volume or runtime movement rather than tuning it away.
