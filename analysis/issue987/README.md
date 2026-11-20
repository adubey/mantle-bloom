# Issue 987: the flank-recipient search no longer grows with recipient distance cubed

RFC 0018's damage pass plans, for every plate with a rift line, where each straining band cell's thinned volume goes
(`rift_strength.flank_plan`). A donor takes the unbroken, same-crust-type cells outside the band within `flank_rows`
hops, widening with no limit until it has some. The old search expanded a donor x node matrix one hop at a time, the
whole disk of every donor, until the farthest donor settled: about `donors x hops^3`. As a widening broken corridor
pushed recipients away, animation steps with braking off grew to hours.

Design principle 2 (bounded, vectorized work over repeated passes over the same graph) was the one violated. The
fix keeps the plan exactly and changes only how it is found:

1. Every donor searches its `flank_rows` hops. Most settle here, at the old cost.
2. For the rest, one multi-source search from the recipients gives each cell's hop distance to the nearest. That is
   the donor's reach. Each donor then searches breadth-first only through cells `k` hops out whose own distance is
   at most `reach - k`, which are the cells on shortest paths to the nearest recipients.

Cost is now about `donors x hops^2` at worst (the corridor), and never more than `donors x plate cells`.

## Evidence

Script: [`bin/debug/time_flank_plan.py`](../../bin/debug/time_flank_plan.py), run in this branch and in a detached
`origin/main` worktree. Seed 106533960, node density 4 (124k nodes), plates only, braking off, `k_cap = 0.14` (the
#888 D1 arm strength, so cells break). Compact results: [`flank-ab-d4.json`](flank-ab-d4.json).

**Same model.** After 10 steps of 1 Myr (935 broken cells), the rift and column state digests match:
`5447bc37033da80e` on both. All 18 benchmarked plans below are bit-identical (reach, held, hops).

**Stepping.** The flank search took 0.022-0.035 s per step on the branch. On `origin/main` it took 0.020-0.046 s,
except 0.23 s and 0.12 s on the two steps where `hops_max` reached 40. The whole damage pass took 0.69-0.81 s per
step on the branch and 0.70-1.14 s on `origin/main`.

**Receding recipients.** These are the three largest plates with a line, after the 10 steps. Every unbroken cell
within K hops of the band is marked broken, and the timing is for one `flank_plan` call:

| Plate (cells, donors) | K hops | Broken cells | `hops_max` | Held | origin/main | Branch | Speed-up |
|---|---|---|---|---|---|---|---|
| 7 (16,859, 382) | 0 | 0 | 4 | 0 | 0.002 s | 0.004 s | 0.6x |
| | 10 | 2,981 | 19 | 0 | 0.022 s | 0.007 s | 3x |
| | 20 | 5,254 | 27 | 0 | 0.065 s | 0.010 s | 7x |
| | 40 | 8,866 | 78 | 0 | 1.54 s | 0.025 s | 61x |
| | 80 | 13,398 | 101 | 0 | 1.93 s | 0.048 s | 40x |
| 4 (15,797, 355) | 0 | 0 | 4 | 0 | 0.002 s | 0.004 s | 0.5x |
| | 10 | 2,900 | 37 | 0 | 0.112 s | 0.007 s | 17x |
| | 40 | 10,256 | 73 | 0 | 0.578 s | 0.028 s | 21x |
| | 80 | 15,291 | 175 | 0 | 4.51 s | 0.078 s | 58x |
| 8 (11,318, 335) | 0 | 0 | 12 | 0 | 0.011 s | 0.004 s | 3x |
| | 40 | 10,392 | 74 | 0 | 0.641 s | 0.023 s | 28x |
| | 80 | 10,982 | 0 | 335 | 7.57 s | 0.001 s | 7,000x |

The last row has no recipients left: the old search walked the whole plate component hop by hop before holding
every donor, and the new one sees the infinite distance at once. When some donors need phase 2, the near case
(K = 0) costs about 2 ms more per plate, for the plate-wide distance search.

**240 steps (branch only).** The same world and settings, stepped 240 times. A frame is 10 steps. Compact
results: [`steps240-d4.json`](steps240-d4.json). Peak RSS was 3.3 GB.

| Steps | Median step | Median damage pass | Median flank search | Flank per plate | `hops_max` | Held, max | Broken cells | Plates |
|---|---|---|---|---|---|---|---|---|
| 1-40 | 4.7 s | 0.87 s | 0.038 s | 1.5 ms | 40 | 4 | 7,558 | 35 |
| 41-80 | 6.9 s | 1.24 s | 0.062 s | 1.3 ms | 44 | 128 | 19,260 | 55 |
| 81-120 | 7.6 s | 1.47 s | 0.072 s | 1.1 ms | 54 | 100 | 27,775 | 69 |
| 121-160 | 8.5 s | 1.65 s | 0.086 s | 1.1 ms | 85 | 92 | 35,757 | 89 |
| 161-200 | 9.7 s | 1.90 s | 0.099 s | 1.1 ms | 51 | 163 | 41,538 | 98 |
| 201-240 | 10.5 s | 1.99 s | 0.099 s | 1.0 ms | 65 | 133 | 43,012 | 102 |

The flank search's cost per plate stays flat (it even falls) as broken cells grow from 0 to 43,000 and `hops_max`
reaches 85. On `origin/main`, one plan on a plate of this size took 0.6-7.6 s with recipients 40-80 hops away (the
benchmark above).

The issue's whole-step target is **not met** in this run: p95 step time is 16.5 s, 4.4x the first frame's median of
3.8 s. That growth comes from other terms:

- The plate count grows from 19 to 102 under `k_cap = 0.14`. The damage pass costs a steady 20-33 ms per plate,
  most of it `prepare_plate`, so its total grows with the plate count.
- The rest of the step grows from about 3.9 s to 8.5 s. This run doesn't isolate that growth. The fault-count
  growth the issue noted is one candidate.

Both are outside this change. The default `k_cap` strains this seed far less (#987's run 4), so this arm is a stress
case, not the UI's default.

## Limits

- Solution 2 of the issue (capping the widening, which changes RFC 0018 section 4) is not part of this change. The
  plan's semantics are as the RFC specifies them.
- With `flank_rows >= 3`, the old loop held every donor on a plate whose components it exhausted before hop
  `flank_rows`. The new search gives them their recipients, as the per-donor reference test expects. The default,
  and every setting in the repo, is 2, where the results are identical.
- Absolute step times come from a machine shared with other agents' runs.
