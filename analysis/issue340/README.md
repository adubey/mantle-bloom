# Issue #340: cratons lose crust only to margin attrition

Replays of the issue's two saves, split by the site that booked each craton loss
(`bin/debug/split_craton_rework.py`), plus a 300 Myr, three-seed rerun of the calibration
behind the old `CRATON_RETREAT_DELAY_YEARS` (`bin/debug/compare_cratons.py`).

| Save | SHA-256 | Span | Steps |
|---|---|---|---|
| `mantle-bloom-seed997271774-87400000y.mbworld` | `08ad9461…ed762b` | 87.4 → 130.5 Myr | 431 × 100 ky |
| `mantle-bloom-seed910211954-220200000y.mbworld` | `bd487409…b060d2` | 220.4 → 240.6 Myr | 202 × 100 ky |

Runs:

- **main**: `origin/main` at `6b4e102`.
- **steps 1+2** (seed 997271774 only): `4e1901d`, craton immobile in orogenic flow, anatexis
  and delamination.
- **branch**: `f721e73`, all of #340 plus the two merge fixes the replays found (below).
- **no braking** (seed 997271774 only): `f721e73` with `CRATON_COLLISION_FRICTION_GAIN = 0`.

Every replay ran alone (one at a time), from a snapshot of its commit's `backend/app` and
`bin/debug`.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \
  ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 --out branch.json \
  --baseline main.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \
  --rescore branch.json --baseline main.json
backend/.venv/bin/python bin/debug/compare_cratons.py --seeds 7,11,2024 --variants on \
  --years 300e6 --workers 1 --out cratons.jsonl
```

The compact comparison is [`replays.json`](replays.json). The full replay JSONs are not
checked in.

## Acceptance

"Non-attrition loss" is the craton booked to every sink account over the replay, less what
the `cratons.attrite` site booked, per Myr. It counts rifting, subduction other than
attrition, delamination, collisional reworking, erosion, topology removal, audit clips and
unattributed loss. The target is at most 20% of main's rate on the same save and span.

| Save | Non-attrition loss, branch / main | Ratio | Craton-edge relief, branch / main |
|---|---:|---:|---:|
| seed 997271774 | 0.29 / 4.41 M km³/Myr | **6.5%** | 1,732 / 1,610 m |
| seed 910211954 | 1.04 / 5.95 M km³/Myr | **17.4%** | 3,575 / 3,482 m |

Both ledgers close: the balance error is 0.0 km³ on every run.

Craton-edge relief is the mean relief of non-cratonic continental crust within 300 km of a
craton, at every craton edge, not only edges facing a convergent front.

## What was lost, by site (M km³)

| Sink and site | 997 main | 997 branch | 910 main | 910 branch |
|---|---:|---:|---:|---:|
| Orogenic flow (`relax_orogens`) | 172.4 | 0 | 86.5 | 0 |
| Anatexis | 3.9 | 0 | 0.13 | 0 |
| Suture accretion and crust transfer | 0.35 | 0 | 4.6 | 0 |
| Shed roots (`delaminated_m3`) | 0.81 | 0 | 6.5 | 0 |
| Merge and topology audit clips | 4.2 | 0 | 0.65 | 0 |
| Faults-phase clip | 0.07 | 0.22 | 0.06 | 0.14 |
| Rifting | 5.1 | 5.3 | 20.8 | 19.5 |
| Erosion | 3.1 | 6.8 | 1.0 | 1.3 |
| Trench attrition (allowed) | — | 0 | — | 1.0 |
| **Craton at end** (start 2,544.2 / 2,047.5) | 2,354.3 | 2,531.9 | 1,928.8 | 2,027.4 |

What remains is mostly real rifting: rift opening, divergent thinning and, on seed 997271774,
one failed rift (`merge_split.maybe_split_plate`, 4.9). Erosion doubles on seed 997271774,
where craton edges stand higher. Seed 910211954 forms 1.9 M km³ of new craton on the branch
against 1.5 on main.

## Two merge paths that destroyed craton

Steps 1+2 alone took seed 997271774 to 29% of main's rate. A plate-merge clip grew from 4.2 to
44.2 M km³, because cratons the front no longer drained were left for a merge to stack. Scratch
probes (taps on the step's functions, not checked in) found two paths, both now fixed
in `quad_merge.py`:

- **Stacking past the Hc cap.** At 114.6 Myr on main, a merge stacked two cratons on 412
  cells. Their summed craton (median 86 km) exceeded the 84 km cap, the overflow Hc was placed
  elsewhere without its craton, and the merge audit clipped 4.2 M km³ to `delaminated_m3`. The
  craton above the capped column now goes with the overflow, keeping its oldest date.
- **The crust-type vote.** At 240.4 Myr on seed 910211954, defragmentation accreted a stranded
  "terrane" onto plate 45. The fragment was 29 oceanic-coded cells at the 84 km cap carrying
  8–40 km of continental material. The merge votes crust type by crust volume, so it retyped
  26 cratonic columns oceanic, and the topology audit stripped 4.4 M km³ to `unattributed_m3`.
  A column holding craton now stays continental. Main does the same on a smaller scale in
  that window (0.6 M km³). Oceanic-coded columns at the Hc cap are a separate anomaly, not
  addressed here.

## Side effects at the end of each replay

| Measure | 997 main | 997 branch | 997 no braking | 910 main | 910 branch |
|---|---:|---:|---:|---:|---:|
| Cratonic area (M km²) | 66.61 | 69.07 | 68.47 | 64.62 | 66.01 |
| Cratonic area under Hc > 50 km (M km²) | 21.8 | 26.5 | 25.8 | 31.6 | 33.7 |
| Land area (M km²) | 141.7 | 146.0 | 143.8 | 156.9 | 158.9 |
| Emergent continental area (M km²) | 137.7 | 141.3 | 139.8 | 122.0 | 121.4 |
| Land elevation p10 / p50 / p90 / p99 (m) | −218 / 765 / 6,307 / 8,457 | −267 / 728 / 6,382 / 8,476 | −224 / 768 / 6,301 / 8,434 | −369 / 1,834 / 5,522 / 8,316 | −277 / 1,971 / 5,478 / 8,135 |
| `ceiling_overflow_m3` added (km³) | 0 | 0 | 0 | 779,676 | 223,860 |
| Mean plate surface speed (cm/yr) | 2.59 | 2.93 | 2.67 | 4.09 | 3.54 |
| Plate merges | 1 | 2 | 2 | 0 | 0 |
| Runtime (s) | 2,307 | 2,596 | 2,052 | 5,383 | 3,856 |

- **Cratons don't spread into belts.** Cratonic area never grows past its start, on either
  save. What does grow is the cratonic area under orogenic crust: orogens now pond on top of
  cratons rather than draining them.
- **Land.** The land and emergent-continental gaps are within about ±3% of main. Run-to-run
  noise on continental area is about ±2% at 43 Myr (#315). Seed 997271774's +3.0% land on the
  branch falls to +1.5% without braking.
- **Ceiling overflow.** Ponding against cratons sends no more crust to the Hc cap. Seed
  910211954 books 71% less; seed 997271774 books none on any run.
- **Plate speed.** On seed 910211954, with no merges, the branch's plates end 14% slower. On
  seed 997271774 all runs agree to 0.02 cm/yr until the first merge at about 115 Myr. After
  it the topologies differ, so the speeds there can't isolate the braking.
- **What braking does.** Turning it off changes craton loss by under 1 M km³. It accounts
  for about 40% of the craton-edge relief gain (1,653 m without it, 1,732 with it, 1,610 on
  main).

## The old calibration, rerun

`CRATON_RETREAT_DELAY_YEARS`' comment recorded that never consuming cratons "kept no more
continent than 60 Myr did". Rerun at the same settings (300 Myr, 5 Myr steps, seeds 7, 11 and
2024), main against the branch:

| Seed | Craton kept at 300 Myr, main / branch | Genuine continental area at 300 Myr, main / branch | Mean plate rotation, main / branch |
|---|---:|---:|---:|
| 7 | 69% / 96% | 144.2 / 145.0 M km² | 10.6 / 9.0 × 10⁻⁹ rad/yr |
| 11 | 59% / 88% | 80.7 / 89.5 M km² | 8.7 / 6.0 × 10⁻⁹ rad/yr |
| 2024 | 74% / 106% | 95.4 / 106.3 M km² | 6.7 / 5.9 × 10⁻⁹ rad/yr |

The conclusion changes. The old comparison changed only boundary consumption, while orogenic
flow, melt and delamination went on draining craton, so cratons that escaped retreat were
lost anyway and kept no extra continent. Now that the craton is immobile everywhere, the
branch keeps 0.6%, 11% and 11% more genuine continental area at 300 Myr. Plates rotate
12–31% slower, with the same number of merges per seed. Seed 2024 ends with more craton than
it started with, from formation (279 M km³).

Two sinks grow over 300 Myr on the branch and are worth watching:

- **Erosion.** Seed 11 books 292 M km³ of craton eroded on the branch, against 99 on main.
- **`subducted_m3`.** It stays large on the branch: 107, 300 and 93 M km³. `compare_cratons.py`
  doesn't split it by site, so this run can't say how much is attrition and how much is
  no-outlet subduction. At 5 Myr steps, attrition removes about 1.6 km of a 125 km cell's
  craton per step.

## Not addressed here

- **Attrition placement.** Attrition acts on the craton plate's own retreating edge at an
  oceanic-override contact (#177), not on an overriding forearc. Real subduction erosion
  removes the upper plate's forearc. It also uses the nonaccreting-margin rate everywhere.
- **Partition rule.** Anatexis treats craton as spread through its column, not as a basal
  layer (`docs/simulation-model.md`).
- **Double-counted jams.** A jam is resolved only by a merge. Until then both overlapping
  cratons exist, counted twice in area.
