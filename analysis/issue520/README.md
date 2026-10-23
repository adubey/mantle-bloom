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
