# Issue 808: plate outline rebuilds

Save: seed 443034896, 67.4 Myr, 127,936 nodes, 8 plates (B2: `mantle-bloom-seed443034896-67400000y.bloom`).
`bin/debug/profile_phases.py --scenarios step --steps 3 --warmup 1`, run alone, origin/main
(b7885e2e) vs this branch.

## Cause

`_stitch_xyz_loops` rebuilt a k-d tree over the whole growing outline once per extra boundary
loop. Plate 2 has 679 loops (12.9 K vertices) and plate 8 has 296, so each call was quadratic
(0.53 s and 0.14 s per call with the loops cached). Everything else in a cold rebuild is
~8-15 ms per plate. Plates in the 18-28 plate saves have few loops, which is why the cost
only shows on this save.

## Change

* Stitch with k-d trees over geometrically growing vertex groups plus linked insertion
  (O(V log V)).
* Cache the stitched local outline with the other topology caches, so a rotation re-poses it.

## Results

| | main | branch |
|---|---|---|
| Step wall time (mean of 3) | 12.73 s | 10.61 s |
| "Outlining plate territory" | 2.63 s (20.6%) | 0.54 s (5.1%) |
| Calls per step | 29.3 | 29.3 |

`bin/debug/bench_plate_outline.py` (per plate, ms):

| Plate | Cells | Loops | Cold main | Cold branch | Warm main | Warm branch |
|---|---|---|---|---|---|---|
| 2 | 21,905 | 679 | 586 | 125 | 531 | 0.04 |
| 8 | 24,848 | 296 | 179 | 70 | 143 | 0.03 |

`get_neighbours` is identical for every plate on the save (exact equality). The outline's seam
positions can differ (the keyhole may join a different equidistant vertex), so the outline
digest changed; the enclosed territory did not.

## Not done

The recompute count is unchanged (29.3 per step); the remaining cost is the cold rebuild
(`_probe_neighbour_indices` is ~80% of it on ordinary plates).
