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
