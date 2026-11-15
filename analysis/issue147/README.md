# #147: step profile after the glacier march, and three result-identical speedups

Measured at `f25174b2` (`origin/main`) with
[`bin/debug/profile_phases.py`](../../bin/debug/profile_phases.py) on a 6-core Apple Silicon
machine with 8 GB, load 2-4. Inputs are the UI defaults: node density 4, climate density 4,
fluid density 2, diagnostic wind, 100 kyr steps.

## Where a step goes now

The profiled runs carry cProfile overhead, so read the shares, not the seconds.

| Phase | Seed 0, new world (8 steps) | Seed 443034896 at 67.4 Myr (8 steps) |
|---|---:|---:|
| Glacier march: shaping ice sheets (`ice_sheet._shape`) | 38.5% | 34.0% |
| Glacier march: re-solving sea level (`eustasy._solve_sea_level_connected`) | 15.9% | 28.9% |
| Glacier march: own work (signed distance, graph restriction, sea components) | 6.9% | 9.5% |
| Resolving lakes | 3.2% | 13.1% |
| Everything else (faults, plate motion, gaps, deformation, climate, routing) | about 35% | about 14% |

The march (RFC 0015) is 72% of a step on the new world and over 90% on the mature one. A
mature step runs 24 to 205 s depending on how many substeps the march takes: the 8 steps
averaged 172 substeps, each solving about 7 ice bodies and re-solving sea level once. For
comparison, #801 measured 6.3 s and 14.8 s per step on the same two worlds before the march.

## What changed here

All three give the same results as before; only the work done to reach them changed.

| Change | Before | After |
|---|---|---|
| Sea-level bisection: one flood from the ocean seed gives every candidate level's connected-ocean mask (`eustasy._seeded_flood_levels`), instead of a union-find over all nodes per candidate (about 18 candidates per solve) | 69.1 s | 11.1 s |
| Lake spill order: a heap instead of re-sorting the ready list after every release (`lakes._ordered_spill_bodies`) | 9.1 s | 1.3 s |
| Lake merge check: one cell's level on scalars, stopping at the first child short of the saddle (`lakes._already_merged`, inside `governing_bodies`) | 11.3 s | 1.8 s |
| **Three steps from the 67.4 Myr save, total** | **237.2 s** | **164.0 s** |

Seconds are summed over the same three steps (35, 25 and 177 s before; 28, 20 and 117 s
after), unprofiled, old and new run back to back.

Equivalence: after those three steps every array in the exact save is equal between old and
new (3,671 members, only the manifest's bytes differ), and the solved sea level matches to
the last bit at each step. The same holds for one step of a seed-0 world at node density 1,
and for sea level over four steps at density 2. The flood is also checked against the
union-find mask in `test_flood_levels_reproduce_the_seeded_mask_at_every_level`.

## Second pass: march bookkeeping (after #904)

Five more changes, again with the same results. Function times are cProfile cumulative
seconds over the same three steps of the 67.4 Myr save (222 s profiled in total before),
so they carry profiler overhead; the last row is unprofiled.

| Change | Before | After |
|---|---:|---:|
| Inside distance of the signed-distance field: searched on the positive set and its boundary instead of from every node outside it (`ice_sheet._inside_distance_local`, already proven exact in #822, now the only path) | 8.7 s | 3.7 s |
| Below-sea components: a union-find over the neighbor table instead of a sparse graph built per call (`hydrology._below_sea_components`) | 9.1 s | 3.8 s |
| Graph restriction: CSR rows filtered directly instead of a round trip through COO (`ice_sheet._restrict`; one more call per substep than before, for the inside distance) | 5.5 s | 0.6 s |
| Median edge length: found once per march instead of once per substep (`ice_sheet.median_edge_m`) | 2.5 s | 0.0 s |
| Sea-level flood: the undirected edge lists are built once per march (`eustasy.flood_adjacency`) | 12.0 s | 9.8 s |
| **Three steps from the 67.4 Myr save, unprofiled, old and new back to back** | **172.1 s** | **161.8 s** |

The per-step times were 28.6, 20.3 and 123.2 s before and 29.2, 21.4 and 111.2 s after, so
the gain is in the long march and within machine noise on the short ones.

Equivalence, checked with [`bin/debug/exact_step_ab.py`](../../bin/debug/exact_step_ab.py):
after those three steps all 2,961 arrays and 320,735 scalars of the exact save are equal
between `origin/main` (`849e30c1`) and this branch, as are the 1,450 arrays of a seed-0
world at density 1 after two steps. Only the stored wall-clock timer
`collision_polarity_stats.prepass_seconds` differs. Compare fields by name, as that script
does: the numbering of a save's `.npy` members changes from one process to the next, so
two runs of the same code differ in hundreds of members when paired by file name.

Tried and dropped: `np.unique` in `lakes._body_leaf_mapping` in place of the Python set.
Most bodies have a handful of members, where `np.unique` costs more than the loop (3.8 s
became 8.9 s).

## Third pass: the opening extent, found near last substep's margin

Measured at `0d4505ee`, where the trend-adaptive march is the default (#842), so the
three steps of the 67.4 Myr save take 84 s where the second pass measured 162 s.

Each body's substep solve has two searches: the extent whose profile holds the opening
volume (`_volume_extent`, a bisection over the body's whole domain), then the margin root
find that walks out from it. Counting profiles (`_shape` calls) per body solve over those
three steps (3,109 body solves) corrects the first pass's "13 of 16":

| Profiles per body solve | Before | After |
|---|---:|---:|
| Opening extent (`_volume_extent`) | 5.3 | 1.4 |
| Margin root find (walk, then Illinois to 100 m) | 11.1 | 11.1 |
| **Three steps, unprofiled, old and new back to back** | **83.6 s** | **63.8 s** |

The step times were 27.9, 19.7 and 36.0 s before and 23.8, 16.9 and 23.1 s after.

The bisection itself is unchanged. The march now passes each body the margin its nodes
ended the previous substep with, and `_volume_extent` first builds the profiles at the two
ends of the bisection cell that margin falls in. V(d) is nondecreasing, so every other
midpoint the bisection tests is then on a known side and needs no profile. A hint that
turns out wrong costs a few more probes, doubling outward, and the bisection still
returns the same extent. 93% of body solves had a hint.

Equivalence, with `bin/debug/exact_step_ab.py`: after those three steps all 2,961 arrays
and 320,735 scalars are equal between `origin/main` (`0d4505ee`) and this branch, and
likewise the 1,450 arrays of a seed-0 world at density 1 after two steps. This rests on
V(d) being nondecreasing at the extents compared. The profile sweep keeps that order by
construction, but summing a volume in floating point can break it by about 1e-9 of the
volume, so an opening volume that falls inside that gap could in principle send the
hinted search one cell away from the full bisection. None did in these 3,109 solves.

## Fourth pass (not merged): starting the root find at the last margin

`MarchCost.margin_start` (preset `trend-1k-margin`, off by default) starts each body's
margin root find at the margin its previous substep ended on, stepping out from 100 m and
doubling, in place of the opening extent rounded to an edge length and a whole-edge
stride. Three steps of the 67.4 Myr save, one run per arm:

| Arm | Profiles per body solve | Body solves | Substeps per step | Closing ice, km^3 |
|---|---:|---:|---|---:|
| `trend-1k` (default) | 12.5 | 3,109 | 24, 18, 90 | 163,201 |
| `trend-1k-margin` | 5.9 | 6,799 | 24, 70, 90 | 178,483 |
| `trend-1k`, `ICE_ROOT_TOLERANCE_M` 90 in place of 100 | 15.0 | not recorded | 24, 90, 90 | 176,623 |

The option does build fewer profiles per solve, but it is not a change "within the
existing tolerance": the ice closes 9.4% larger, and the second step takes 70 substeps to
settle where the default took 18, so the three steps were no faster (93 s against 84 s).

The third row is the control that explains it. Tightening the root tolerance from 100 m to
90 m, with nothing else changed, moves the closing ice by 8.2% and one node by 1.9 km of
ice. On this world the march's result depends on which of several sign changes of the
budget function a body's root find lands on, so any change to the search, however small,
gives a different trajectory of about this size. A single three-step comparison therefore
cannot say whether `margin_start` is better, worse or the same. It needs the kind of
evidence #842's round 3 used: many steps on several saves, both arms, judged by ledger
closure, step robustness and swings.

## What is left, largest first

Shares are of the profiled three steps after the second pass (192 s), before the third
pass and the trend-adaptive march.

1. **Ice shape solve (about half).** 83% of it is the
   ordered sweep, `_plastic_sweep`, which is already close to what the hardware allows: an
   indexed heap in place of `heapq` tuples bought 15%, and solving bodies on threads topped
   out at 1.5x. The cost is how many sweeps run. After the third pass a body takes about
   12.5 profiles per substep, 11 of them in the margin root find, and each profile sweeps
   about 3 times while it drops floating cells. The root find walks out from the opening
   extent one edge length (about 60 km) at a time, then narrows that bracket to 100 m. Its
   function is not monotone, so nothing can be inferred as in the third pass: starting it
   from a narrow bracket around last substep's margin would save roughly half of those
   profiles, but ends on a different 100 m bracket, so results move within the existing
   tolerance. That needs a replay comparison and belongs with the march-cost work (#822,
   #842). The mass
   balance inside each evaluation (`SurfaceClimate.rates`, 19 s) is numpy on small arrays;
   its `arcsin` rules out a bit-identical numba port.
2. **Lake march (about 17%).** `_settle_prognostic_slice` is Python loops over small numpy
   arrays: `budget` (0.3 M calls, 7 s), `mean_depth` (0.4 M), `_body_leaf_mapping` (4 s),
   and a full walk of the lake forest in `governing_bodies` several times a slice. A body's
   leaf ids never change within a solve, so `_body_leaf_mapping` could cache them per body.
3. **Sea-level flood (5%).** What remains is the flood itself, 17 ms per solve at 128 k
   nodes: about 70% of the nodes are seeds and all go through a heap of tuples.
4. **Outside distance (under 2%).** `_outside_distance` is one scipy Dijkstra over the
   eligible graph per substep, and it needs the predecessors.

## Reproduce

```bash
backend/.venv/bin/python bin/debug/profile_phases.py --seed 0 --scenarios step --warmup 2 --steps 8
backend/.venv/bin/python bin/debug/profile_phases.py \
    --load ~/Downloads/mantle-bloom-seed443034896-67400000y.mbworld \
    --scenarios step --warmup 2 --steps 8
```
