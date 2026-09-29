# Issue #228 Phase 4: gap filling, volcanism, overlap tracking and relattice on quad worlds

Seed 7, `node_density` 1, 1 Myr steps. `main` is `99d73e2`.

Reproduce:

    cd backend
    .venv/bin/python ../bin/debug/measure_quad_passes.py --seed 7 --steps 60                                   # main-seed7-60.json (quad + lines)
    .venv/bin/python ../bin/debug/measure_quad_passes.py --seed 7 --steps 150 --every 25 --surfaces quad       # *-seed7-150-quad.json
    .venv/bin/python ../bin/debug/ablate_quad_passes.py --seed 7 --variant no-gapfill                          # one change reverted

`uncovered` and `overlapping` come from a 200k-point Fibonacci sample tested against every
plate's `contains_batch`. That is exact for quad plates, whose cells are their territory.

## What was wrong on `main`

| after 60 Myr, quad | main | lines (reference) |
|---|---:|---:|
| sphere uncovered | 2.6% | (polygon test, not exact) |
| gap sweep's uncovered fraction | 0.0% | 0.0% |
| overlap tracker recall vs containment | 0.61 | — |
| volcano nodes (cumulative) | 4,427 | 8,798 |
| eruptions, last 10 Myr | 1,084 | 3,076 |

- **Gap filling never fired.** The sweep calls a point covered when a node lies within 1.5
  spacings. Boundary advance keeps a new cell's centre 1.3 spacings from a neighbour's
  nodes, which leaves a seam about one cell wide between plates. 93% of uncovered sample
  points had two plates within 1.5 spacings, and none was farther than 1.28 spacings from a
  node. All of it was under `MIN_GAP_NODES`: 166 clusters, the largest 276 cells.
- **Overlap tracking missed about a third of real overlap.** Two plates' lattices sit in
  different frames, so their nodes are rarely within the 0.5-spacing proximity tolerance
  even when one plate's centres lie inside the other's cells.
- **Quad ridges were volcanically silent.** In the line engine, every node a row claim or gap
  fill creates is seeded thin and erupts (`seed_and_erupt_new_nodes`). Quad's `_open_rift`
  seeds magmatic cells at the full column, so they only erupted if stretching melted them.
- **Relattice was skipped by accident.** `relattice_continental_plates` duck-types, and quad
  plates happen not to have the method.

## Changes

- **Overlap tracking.** `PlateSurface.territory_is_exact` (true on quad). When every plate
  is exact, `compute_node_overlap` flags a node whose centre another plate contains.
  Candidate pairs come from bounding caps, so a buried plate is still seen.
- **Gap filling.** On exact worlds the sweep uses containment. Every cluster below
  `MIN_GAP_NODES` is pooled and grown into by its adjacent plates, once per plate, with no
  standoff from neighbour nodes and without logging events. `MIN_GAP_NODES` still gates
  spawns and the logged per-cluster growth.
- **Volcanism.** A new cell that is mostly magmatic (stretch share < 0.5) ignites as a
  volcano without changing its column. The quad eruption roll is keyed by each node's stable
  ID (`_node_uniforms`), so inserting or removing cells no longer reshuffles other draws.
- **Relattice.** None on quad, by design; documented and tested.

Line worlds are unchanged. A hash of two 9-step line worlds (seed 3, density 0.5, one with
diagnostics) matches `main` byte for byte.

## Results, quad

| seed 7 | main 60 Myr | branch 60 Myr | main 150 Myr | branch 150 Myr |
|---|---:|---:|---:|---:|
| sphere uncovered | 2.6% | 0.7% | 3.0% | 2.2% |
| sphere doubly covered | 0.6% | 1.1% | 0.6% | 0.7% |
| overlap recall | 0.61 | 1.00 | 0.66 | 1.00 |
| volcano nodes | 4,427 | 6,583 | 10,713 | 14,973 |
| eruptions, last reporting window | 1,084 / 10 Myr | 1,629 / 10 Myr | 3,422 / 25 Myr | 4,426 / 25 Myr |
| one-cell-thin cells | 0.36% | 0.62% | 0.59% | 0.64% |
| hole loops | 0 | 1 | 0 | 0 |
| crust volume (10⁹ km³) | 12.67 | 12.69 | 9.53 | 9.10 |
| land fraction | 0.420 | 0.408 | 0.289 | 0.273 |

- **Seams.** Uncovered area follows the 4-step gap-fill cadence. Right after a pass it is
  about 0.7%, the same as at generation (0.75%): whole cells on two unaligned lattices can't
  tile a boundary exactly. Advance reopens seams to about 2% before the next pass. Removing
  advance's standoff would close them at the source, but that changes deformation and is left
  as a follow-up.
- **Crust volume is within run-to-run spread.** At 150 Myr, branch vs main is −4.5% on seed 7,
  +3.9% on seed 11 (9.62 vs 9.31) and +1.6% on seed 23 (5.04 vs 4.96). Reverting any single
  change on seed 7 lands between 9.29 and 9.45.
- **Relattice isn't needed.** Over 150 Myr the quad lattice stays unrefined, under 1% of
  cells are one cell thin, and at most one transient hole appears. The line lattice reaches
  15% one-node lines by 60 Myr (Phase 0a).

## Cost

Measured sequentially on an otherwise quiet machine, 40 steps each, seed 7:

| s/step | main | branch |
|---|---:|---:|
| mean | 1.21 | 1.30 |
| gap-fill steps (every 4th) | 1.54 | 1.92 |
| other steps | 1.10 | 1.09 |

Pooling the seam clusters took one gap-fill pass from 4.2 s to 0.45 s. Before pooling, each
of about 200 clusters ran its own frontier walk and rebuilt every plate's KD-tree. Overlap
tracking costs about 0.02 s per step.
