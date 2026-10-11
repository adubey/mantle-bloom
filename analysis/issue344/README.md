# Issue #344: quad merges kept out of third plates' cells

Part of #344 (the scheduled stress run's three failures). This covers
`test_quad_merge_of_neighbouring_plates_in_a_stepped_world`, which bisects on main's
first-parent history to `4c06455` (#259, gap filling every step on quad worlds).

## What was wrong

`quad_merge._transfer` activated a new survivor-lattice cell unless its centre was nearer a
third plate's node than either fusing plate's. Where the absorbed plate already overlapped a
third plate, its own nodes won that test, so the rotated lattice put new cells' centres inside
the third plate. On the stress test's world (seed 5, 6 plates, 3 × 1 Myr), merging the
largest plate with its largest neighbour:

| Overlap with the other plates (m²) | before merge | after, `main` | after, fix |
|---|---:|---:|---:|
| Sampled at cell centres (the test's measure) | 8.15e11 | 9.64e11 (+18%) | 7.88e11 |
| Sampled at 8 × 8 sub-cells per cell | 1.654e12 | 1.719e12 (+4%) | 1.628e12 |

The 6 new cells inside plate 0 had 0.68–0.75 coverage by the absorbed plate and centres
0.2–0.5 spacings from its nodes, against 0.5–0.6 from plate 0's. The fix rejects any
candidate cell whose centre another live plate contains (quad territory is exact, as
`gaps.py` already treats it). The rejected cells' volume is restored across the new cells as
before, so the merge still conserves every extensive field.

## Replay

Seed `997271774` from the 87.4 Myr save, 431 × 100 ky steps to 130.5 Myr (save SHA-256
`08ad9461…ed762b`, as in [#320](../issue320/README.md)). `main` is `origin/main` at
`6b4e102`; `fix` is this branch. The two ran one after the other, each alone.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  /path/to/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 --checkpoint-myr 2.5 --out fix.json
python3 bin/debug/compare_crust_transfer.py main=base.json fix=fix.json
```

The compact comparison is [`compare-seed997271774.json`](compare-seed997271774.json); the
full replay JSONs are not checked in. `main` reproduces #320's final column exactly.

### At 130.5 Myr

| Measure | main | fix |
|---|---:|---:|
| Land area | 141.7 M km² | 143.2 M km² |
| Continental area | 237.0 M km² | 236.2 M km² |
| Emergent continental area | 137.7 M km² | 139.2 M km² |
| Land elevation p10 / p50 / p90 / p99 | 167 / 1150 / 6692 / 8842 m | 147 / 1145 / 6700 / 8828 m |
| Continental Hc p50 / p90 / p99 | 33.6 / 70.9 / 84.0 km | 33.5 / 72.4 / 84.0 km |
| Continental Hm-cap fraction (change) | 6.13% (+0.18 pp) | 5.58% (−0.36 pp) |
| Continental material change | −0.31% | −0.38% |
| Collision loss of continental material | 0.0076%/Myr | 0.0074%/Myr |
| Runtime / peak RSS | 2,071 s / 2,874 MiB | 2,251 s / 2,584 MiB |

### Reading it

- **The runs are identical to step 111 (98.5 Myr),** the first merge the fix changes. They
  track within ~1 M km² of each other after that:

  | Myr | 98.5 → 112.4 | 114.9 | 119.9 | 124.9 | 130.5 |
  |---|---|---:|---:|---:|---:|
  | Land, main / fix (M km²) | within 0.3 | 143.5 / 144.1 | 143.0 / 144.4 | 143.1 / 142.1 | 141.7 / 143.2 |

  The land gap changes sign along the way, so it's trajectory divergence, not a trend.
- **Every area gap is inside run-to-run noise** (about ±2% on continental area at 43 Myr,
  #315): land +1.0%, continental −0.35%, emergent +1.1%.
- **The Hm-cap fraction ends 0.5 pp lower.** That's one save, and the run diverges at a merge,
  so treat it as unresolved, not as an effect of the fix.
- **The ledgers close the same way in both runs:** the continental and mobile-cover balance
  errors are equal to 1e-12 relative. (The −2.04e8 km³ continental balance error is the same
  on `main` and predates this change.)
- **Runtime is 9% higher on the fix run.** Both ran alone, but other apps were open on the
  machine. The new containment check is one `contains_batch` per other plate per merge, which
  is too small to account for that. The diverged trajectory and machine load are the likelier
  causes. Median step time is 4.1 s against 4.6 s.
