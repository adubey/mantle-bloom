# Issue #520, part 2: river lines drawn across lakes

A spilling sink's `flow_target` is the node beyond its lake's rim (or, for a lake united within
the step, its lowest sink). The map and the River Inspector drew that edge as one straight river
segment, up to 1,669 km long across an inland sea. The branch draws each spill hop along its
spill path over neighbour edges (`hydrology.flow_segments`), stops at visible lake water, and
ends the upstream network there (`mouth_type` `"lake"`). Routing is unchanged.

Parts 1 (perched-lake false merges) and 3 (floor drains flipping) aren't addressed here.

## Measurement

`bin/debug/measure_river_hops.py`, one 100 kyr step from each save, `origin/main` `f3a1e46`
against the branch. "Long" is a drawn segment longer than 3× the median flow-edge length
(86 km and 93 km); "non-nbr" joins nodes that aren't k-NN neighbours.

| Save | State | Spill hops (river) | Map long, main → branch | Map non-nbr | Inspector long | Inspector non-nbr |
|---|---|---|---|---|---|---|
| seed308161711, 121.5 Myr (evidence) | loaded | 90 (18) | 4 → 0 | 9 → 0 | 6 → 1 | 18 → 1 |
| | +1 step | 100 (36) | 12 → 0 | 26 → 0 | 14 → 0 | 36 → 0 |
| seed374464662, 132.6 Myr | loaded | 24 (5) | 1 → 0 | 3 → 0 | 2 → 0 | 5 → 0 |
| | +1 step | 28 (8) | 1 → 0 | 5 → 0 | 2 → 0 | 8 → 0 |

- The one segment left on the loaded evidence save is a 457 km hop by a sink of a lake united
  within the step. A cache saved before this change has no `should_spill`, so hops are found
  from `flow_target == spill_target`, which misses those pooling hops. The cache is rebuilt
  on the next step.
- Spill-hop counts, river-network counts after the step, and the hops' lengths match between
  the two trees: routing, discharge and lake inflow don't move.
- The longest drawn segment after the step is 111 km (evidence save) and 97 km (other save),
  against 1,029 km and 343 km on `origin/main`.

Saves are local (see the shared handoff note for #520). Commands:

```bash
cd backend
.venv/bin/python ../bin/debug/measure_river_hops.py <save> --myr 0.1 --checkpoint-myr 0.1 \
    --center 43 153 --render-dir <scratch>/branch --out <scratch>/branch.json
.venv/bin/python ../bin/debug/measure_river_hops.py <save> --myr 0.1 --checkpoint-myr 0.1 \
    --backend <origin-main-worktree>/backend --center 43 153 --render-dir <scratch>/main \
    --out <scratch>/main.json
```

## Request cost

The first version recomputed every spill path, in plain Python, each time `flow_segments` was
called: twice per map render and once per river on `/world/rivers`. The branch now works the
paths out once per `HydrologyFields` (`_SpillPathCache`), bounds each search, and runs it as
numba kernels. `measure_river_hops.py --timing 7` on the evidence save, one 100 kyr step on, run
one configuration at a time. Median ms; "cold" clears the path cache first.

| | `03e5dc5` (before) | paths, no cache | cache | cache + numba |
|---|---|---|---|---|
| Elevation map render, cold | 995 | 2,469 | 1,694 | 536 |
| Elevation map render, warm | 981 | 2,356 | 1,005 | 574 |
| `/world/rivers` payload, cold | 76 | 2,469 | 792 | 82 |
| `/world/rivers` payload, warm | 50 | 2,442 | 78 | 45 |
| `flow_segments` over every river node, cold | — | 754 | 719 | 43 |

- Absolute times drift between runs: the same render on `03e5dc5` took 0.3 s in an earlier
  set. Compare within a column; cache + numba's render coming in below the baseline is that
  drift.
- All three branch versions draw the same segments.
- The kernels compile in 0.85 s on first use and load from numba's cache in 0.15 s per process.
- The search's climb ceiling and node cap never stopped a search on this save.

## Part 1 replay preregistration

**Question.** Does checking each child's main basin prevent the false parent merge on sea 26,
without materially changing the unrelated frozen lake 13?

**Comparison.** Replay the 100 kyr evidence save for 10 steps on baseline `origin/main`
`9ae6eaf40fee013f9989b8ad236837927be0496e` and on the issue branch
`d0bd84932009ae35c1f153a969d9b7b46ed838bf`. This baseline already contains the sibling #521
split-reset change and part 2's river rendering fix; both arms use the same save, step length,
tracked bodies (sea 26 and lake 13), and replay script. Save SHA-256 is
`c557251402f9d2df3d2f67b93c2ad528674cba32055791ed089ca56d6214d6e0`.

**Expected direction and threshold.** Primary metric is sea 26's wet-node extent. The issue's
historical baseline was 1,713 nodes at +3 and 1,020 at +8 on `59d2ef3`, before sibling #521
landed. Since #521 changes replay numbers, use the measured current `origin/main` counts as the
comparison; expect the branch to reduce sea 26's wet extent by at least 25% at both +3 and +8.
Secondary metrics are per-step visible lake area and volume, plus merge/split event counts;
expect fewer false parent transitions and no false rise in standing water. Lake 13's size is a
control expected to remain broadly unchanged because its inflation is attributed to the frozen
rule and split reset, not this perched-water check.

**Variability and decision rule.** Both arms start from the identical serialized world and take
the same fixed 100 kyr steps, so the replay is expected to be deterministic. Keep all 10
per-step observations. If sea 26 lands near the 25% reduction threshold, or event counts suggest
a marginal effect, rerun both arms from the same save before interpreting the comparison.

**Run controls.** Run baseline first, then branch, one replay at a time and without unit or
stress tests running. Record exact commits, command, endpoint values, runtime, and peak RSS in
the results section after both runs complete.

## Replay correction preregistration

The first branch replay at `d0bd84932009ae35c1f153a969d9b7b46ed838bf` produced a large step-2
departure from current `origin/main`: 39 vs 14 merges, 5 vs 53 splits, and 1.57M vs 0.52M km³
of visible lake water. The first implementation also used the persisted-water helper to inspect
one floor cell during the same-step merge check, even though the just-resolved child levels are
available there. The correction keeps the helper for persisted state and restores the direct
child-level check for a merge happening in the current step.

**Expected direction and threshold.** Rerun the same baseline save and branch after this change.
At step 2, expect at least half of the first trial's branch-to-baseline gaps to close: no more
than 27 merges, at least 29 splits, global visible lake volume no more than 1.05M km³, and lake
13 volume no more than 529,000 km³. These are diagnostic thresholds for the correction, not a
replacement for the primary sea 26 acceptance measurement above.

**Variability and controls.** Reuse the completed baseline run from `9ae6eaf40fee013f9989b8ad236837927be0496e`; it starts from the same saved state and is unchanged by the correction. Run the corrected branch once at its recorded commit, after all tests have stopped. Retain the first branch run as an exploratory result and report why the correction was made.

## Part 1 replay results

**What was replayed.** The 128,340-node save `seed308161711-121.5Myr`, SHA-256
`c557251402f9d2df3d2f67b93c2ad528674cba32055791ed089ca56d6214d6e0`, for 10 steps of
100,000 years. Current baseline was `origin/main` `9ae6eaf40fee013f9989b8ad236837927be0496e`;
the corrected branch was `14ae639a5f31e30f44b8ee1c71a5d5afcb1e1789`. #521's split-reset removal
was already on the baseline. Each replay ran alone, without unit or stress tests. The first
branch attempt at `88f4abf945070a80d2c9e4dc8750b29440c00b85` is retained in the JSON; it used
the new helper in the current-step merge check as well. The correction was motivated by that
run's event and water-volume shift. The corrected run's step-2 measurements were effectively the
same, so that hypothesis was not supported: it still had 39 merges, 5 splits, 1.57M km³ of
visible lake water, and 1.04M km³ in lake 13. All four preregistered correction thresholds
failed.

**Measures.** The preregistered current-main comparison does not meet the 25% sea-26 reduction
criterion at +8: baseline has 6 wet nodes and the corrected branch has 68. Both have 0 at +3.
The original issue's 1,713 (+3) and 1,020 (+8) nodes were measured on `59d2ef3`, before #521;
the corrected branch is far below those historical extents, but current main already removed
that expansion by landing #521 first. At +10 the branch has 340 sea-26 wet nodes versus 1,685
on current main. Across the 10 steps, the branch has 207 merge and 291 split events, versus 46
and 765 on current main. Lake 13 also changes materially (for example, 671 vs 1,260 wet nodes
at +8), despite the expectation that this control would remain broadly unchanged. The issue
therefore shows a strong topology and water-balance interaction with the #521 baseline; the
global event counts alone do not identify which transitions are false.

| Step | Sea 26 wet nodes (main → branch) | Sea 26 volume km³ | Visible lake area km² | Visible lake volume km³ | Merges (main → branch) | Splits (main → branch) |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 → 0 | 0 → 0 | 2,376,497 → 2,376,497 | 735,164 → 735,708 | 2 → 2 | 83 → 73 |
| 2 | 42 → 42 | 455 → 455 | 2,838,633 → 3,630,305 | 524,447 → 1,574,278 | 14 → 39 | 53 → 5 |
| 3 | 0 → 0 | 0 → 0 | 4,115,398 → 3,311,601 | 1,908,866 → 1,176,986 | 3 → 29 | 98 → 20 |
| 4 | 0 → 0 | 0 → 0 | 2,054,957 → 3,493,854 | 583,388 → 1,850,286 | 2 → 20 | 72 → 5 |
| 5 | 30 → 27 | 395 → 373 | 1,825,698 → 4,268,335 | 463,101 → 1,760,200 | 4 → 26 | 80 → 18 |
| 6 | 15 → 13 | 144 → 130 | 1,744,844 → 5,609,712 | 455,449 → 1,964,652 | 2 → 19 | 87 → 25 |
| 7 | 0 → 50 | 0 → 488 | 1,493,123 → 4,907,186 | 428,551 → 1,864,956 | 3 → 18 | 84 → 40 |
| 8 | 6 → 68 | 64 → 594 | 6,257,890 → 4,776,587 | 2,753,476 → 2,026,376 | 7 → 22 | 63 → 23 |
| 9 | 2 → 0 | 19 → 0 | 5,896,513 → 6,414,787 | 2,745,636 → 1,995,511 | 4 → 16 | 75 → 42 |
| 10 | 1,685 → 340 | 150,153 → 3,512 | 8,647,382 → 4,776,826 | 777,997 → 1,246,655 | 5 → 16 | 70 → 40 |

**Commands.** On each arm, run from `backend/` with its own environment:

```bash
source .venv/mantle-bloom-env.sh
.venv/bin/python -I ../bin/debug/measure_lake_split_reset.py <local-save> \
    --steps 10 --track 26 13 --out <scratch-output>
```

The baseline took 72.13 s of summed step runtime (peak RSS 1,118.62 MB); the corrected branch
took 89.17 s (1,097.19 MB). The exploratory branch attempt took 88.69 s (1,135.05 MB). The
complete compact per-step measurements and run labels are in
`analysis/issue520/replay-part1-10steps.json`; full replay logs and saves remain local.
