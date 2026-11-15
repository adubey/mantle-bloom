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

## What is left, largest first

Shares are estimates of a mature step after these changes, scaled from the profile above.

1. **Ice shape solve (roughly half).** 83% of it is the
   ordered sweep, `_plastic_sweep`, which is already close to what the hardware allows: an
   indexed heap in place of `heapq` tuples bought 15%, and solving bodies on threads topped
   out at 1.5x. The cost is how many sweeps run. Each body takes about 16 shape evaluations
   per substep, 13 of them in `_volume_extent`'s full-range bisection that recovers last
   substep's extent from the opening volume, and each evaluation sweeps about 3 times while
   it drops floating cells. Cutting the evaluations (a bracket carried from the previous
   substep, or a secant search) would move results within the existing tolerances, so it
   needs a replay comparison and belongs with the march-cost work (#822, #842).
2. **Lake march (roughly 10%).** `_settle_prognostic_slice` is Python loops over small numpy
   arrays: `budget` (0.5 M calls), `mean_depth` (0.7 M), `_body_leaf_mapping`, and a full
   walk of the lake forest in `governing_bodies` several times a slice.
3. **March bookkeeping (roughly 14%).** `_signed_distance`, two `_restrict` calls per substep
   that round-trip the graph through COO, and `hydrology._below_sea_components`.
4. **Sea-level flood (7%, measured).** The symmetric adjacency is rebuilt for every solve
   although the node graph is fixed for the whole march.

## Reproduce

```bash
backend/.venv/bin/python bin/debug/profile_phases.py --seed 0 --scenarios step --warmup 2 --steps 8
backend/.venv/bin/python bin/debug/profile_phases.py \
    --load ~/Downloads/mantle-bloom-seed443034896-67400000y.mbworld \
    --scenarios step --warmup 2 --steps 8
```
