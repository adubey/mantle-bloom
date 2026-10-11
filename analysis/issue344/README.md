# Issue #344: quad merges kept out of third plates' cells

Part of #344 (the scheduled stress run's three failures). This covers
`test_quad_merge_of_neighbouring_plates_in_a_stepped_world`, which bisects on main's
first-parent history to `4c06455` (#259, gap filling every step on quad worlds).

## What was wrong

`quad_merge._transfer` activated a new survivor-lattice cell unless its centre was nearer a
third plate's node than either fusing plate's. Where the absorbed plate already overlapped a
third plate, its own nodes won that test, so the rotated lattice put new cells' centres inside
the third plate.

The first fix (`66cde3f`) also dropped a candidate whose *centre* another live plate contains.
The PR #346 review pointed out that a cell centred just outside can still lie across the third
plate's edge. The current fix (`8356842`) samples each candidate's whole footprint (9 × 9 points,
corners and edges included) and drops it if any sample is in another live plate. The dropped
cells' volume is restored across the kept new cells, so the merge still conserves every
extensive field. The seam they leave is closed by the next step's gap filling.

On the stress test's world (seed 5, 6 plates, 3 × 1 Myr), merging the largest plate with its
largest neighbour, overlap with the other plates (m²):

| Sampled at | before merge | after, `main` | after, centre guard | after, footprint |
|---|---:|---:|---:|---:|
| Cell centres (the test's measure) | 8.15e11 | 9.64e11 (+18%) | 7.88e11 | 7.88e11 |
| 4 × 4 sub-cells per cell | 1.617e12 | 1.678e12 | 1.592e12 | 1.155e12 |
| 8 × 8 sub-cells per cell | 1.654e12 | 1.719e12 (+4%) | 1.628e12 | 1.173e12 |

The footprint check keeps 3,028 new cells of the absorbed plate's 3,142, against 3,125 for
`main`. The merged plate covers 98.4% of the pair's area, against 99.8%.

An adversarial case is in `unit_tests/test_quad_merge.py`
(`test_merge_plates_grows_no_new_cell_across_a_third_plates_edge`): the absorbed plate stops
short of a third plate on its own lattice. The centre guard grows 1.78e9 m² of overlap there,
and the footprint check none. A scan of third-plate placements (3 lattice rotations, edge gaps
from −0.6 to +1.0 spacing) finds no layout where the footprint check grows sub-cell overlap.
The check can still miss a third-plate cell reaching less than 1/8 of a cell into a candidate
between two samples.

## Replay

Seed `997271774` from the 87.4 Myr save, 431 × 100 ky steps to 130.5 Myr (save SHA-256
`08ad9461…ed762b`, as in [#320](../issue320/README.md)). Three runs, one after the other:

- **main**: `origin/main` at `6b4e102`.
- **centre**: the first fix, `66cde3f` (recorded as `ace2033`, which adds only this write-up).
- **footprint**: the current fix. It ran on a working tree identical to `8356842`, so its JSON
  records `ace2033`.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  /path/to/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 --checkpoint-myr 2.5 --out footprint.json
python3 bin/debug/compare_crust_transfer.py main=base.json centre=centre.json footprint=footprint.json
```

The compact comparison is [`compare-seed997271774.json`](compare-seed997271774.json); the
full replay JSONs are not checked in. `main` reproduces #320's final column exactly.

### At 130.5 Myr

| Measure | main | centre | footprint |
|---|---:|---:|---:|
| Land area | 141.7 M km² | 143.2 M km² | 141.1 M km² |
| Continental area | 237.0 M km² | 236.2 M km² | 238.0 M km² |
| Emergent continental area | 137.7 M km² | 139.2 M km² | 137.1 M km² |
| Land elevation p10 / p50 / p90 / p99 | 167 / 1150 / 6692 / 8842 m | 147 / 1145 / 6700 / 8828 m | 185 / 1209 / 6731 / 8867 m |
| Continental Hc p50 / p90 / p99 | 33.6 / 70.9 / 84.0 km | 33.5 / 72.4 / 84.0 km | 33.5 / 70.7 / 84.0 km |
| Continental Hm-cap fraction (change) | 6.13% (+0.18 pp) | 5.58% (−0.36 pp) | 5.84% (−0.11 pp) |
| Continental material change | −0.31% | −0.38% | −0.31% |
| Collision loss of continental material | 0.0076%/Myr | 0.0074%/Myr | 0.0073%/Myr |
| Runtime / peak RSS | 2,071 s / 2,874 MiB | 2,251 s / 2,584 MiB | 4,781 s\* / 2,343 MiB |

\*Not a clean timing: another agent's full unit suite (10 workers) ran alongside it for about
half the replay, and swap filled. `main` and `centre` each ran alone.

### Reading it

- **The footprint run is identical to `main` to step 38 (91.2 Myr),** earlier than the centre
  run (step 111), because the stricter check changes more merges. It tracks `main` within
  ~1.3 M km² of land and ~1.1 M km² of continental area at every 2.5 Myr checkpoint:

  | Myr | 91.2 → 112.4 | 114.9 | 119.9 | 124.9 | 127.4 | 130.5 |
  |---|---|---:|---:|---:|---:|---:|
  | Land, main / footprint (M km²) | within 0.3 | 143.5 / 143.1 | 143.0 / 144.1 | 143.1 / 141.8 | 141.6 / 140.4 | 141.7 / 141.1 |
  | Continental, main / footprint (M km²) | within 0.1 | 246.1 / 245.9 | 243.3 / 243.5 | 240.4 / 240.9 | 238.9 / 239.7 | 237.0 / 238.0 |

- **Every area gap is inside run-to-run noise** (about ±2% on continental area at 43 Myr,
  #315). Footprint against `main`: land −0.4%, continental +0.4%, emergent −0.4%.
- **The Hm-cap fraction ends 0.3 pp below `main`,** between `main` and the centre run. With
  one save and three diverged trajectories, that's not resolved as an effect of the fix.
- **The ledgers close the same way in all three runs:** the continental and mobile-cover
  balance errors agree to 1e-8 relative. (The −2.04e8 km³ continental balance error is the same
  on `main` and predates this change.)
- **Runtime:** the footprint check adds one `contains_batch` over ~81 points per candidate per
  other plate, per merge. Merges are rare, so the replay's 2× wall time is most likely from
  the concurrent suite, not the check. A clean timing would need a rerun alone.
