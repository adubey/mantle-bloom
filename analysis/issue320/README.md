# Issue #320: crust transfer at polarized collision fronts

Replay of seed `997271774` from the 87.4 Myr save, 431 × 100 ky steps to 130.5 Myr (save
SHA-256 `08ad9461…ed762b`, as in [issue #317](../issue317/README.md)). Three runs:

- **#314**: the checked-in [`issue314/cascade-seed997271774.json`](../issue314/cascade-seed997271774.json)
  (metrics only).
- **main**: `origin/main` at `4d9a28a` (#319 merged: one-sided retreat, Hm sink, crust kept on
  the lower plate).
- **transfer**: this branch, default partition (scrape 0.70, underthrust 0.27, loss 0.03 of
  the crust below the mobile cover; 200 km underthrust reach).

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  /path/to/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 --checkpoint-myr 2.5 --out transfer.json
python3 bin/debug/compare_crust_transfer.py main=base.json transfer=transfer.json \
  --baseline314 analysis/issue314/cascade-seed997271774.json
```

The compact comparison is [`compare-seed997271774.json`](compare-seed997271774.json); the
full replay JSONs are not checked in.

## At 130.5 Myr

The transfer column is the replay after the PR #328 review fixes (shared root-shedding
allowance per upper plate, no oceanic hand-off, approach aimed at the upper plate, cover
booked where its crust lands).

| Measure | #314 | main | transfer |
|---|---:|---:|---:|
| Land area | 141.5 M km² | 146.1 M km² | 141.7 M km² |
| Continental area | 235.1 M km² | 232.4 M km² | 237.0 M km² |
| Emergent continental area | 137.8 M km² | 141.5 M km² | 137.7 M km² |
| Land elevation p10 / p50 / p90 / p99 | 183 / 1086 / 6796 / 8043 m | 197 / 1207 / 6779 / 8728 m | 167 / 1150 / 6692 / 8842 m |
| Continental Hc p50 / p90 / p99 | 33.3 / 77.6 / 84.0 km | 33.6 / 74.7 / 84.0 km | 33.6 / 70.9 / 84.0 km |
| Continental Hm-cap fraction (change) | 13.62% (+7.67 pp) | 6.15% (+0.20 pp) | 6.13% (+0.18 pp) |
| Continental material change | — | −0.36% | −0.31% |
| Collision loss of continental material | — | 0.0021%/Myr | 0.0076%/Myr |
| Runtime / peak RSS | — | 2,467 s / 2,692 MiB | 2,659 s* / 2,562 MiB |

\*Run alongside the full unit-test suite, so not a clean timing. The first transfer run,
before the review fixes, ran alone and took 2,253 s against main's 2,467 s.

Collision loss is delamination + no-outlet + the lower-crust loss share, per Myr, as a share
of the opening inventory. 0.0076%/Myr is about half of Earth's ~0.015%/Myr orogenic
delamination and at the low end of #276's 0.01–0.05%/Myr seed range. By account: 25.1 M km³
lower-crust loss share, 2.4 M km³ delaminated roots, 0 no-outlet.

Land and emergent area end ~4.4 M km² (3%) below `main`, back on the #314 replay's values;
continental area ends 4.6 M km² above `main`. Run-to-run noise on continental area is about
±2% at 43 Myr (#315). `main` keeps consumed crust on the lower plate's own survivors, which
thickens it into land behind the suture; the transfer moves it onto the upper plate's
frontal belt and under its front, where part of it goes into the root rather than the
surface. The Hm-cap change matches `main` within noise (the pre-fix run ended at −0.19 pp).

### Reading the land gap

Against `main` the transfer costs ~4.4 M km² of land at 130.5 Myr, but the trajectory shows
`main` is the run that drifts:

| Myr | #314 | main | transfer |
|---:|---:|---:|---:|
| 87.4 | 148.8 | 148.8 | 148.8 |
| 107.4 | 147.2 | 148.7 | 148.2 |
| 114.9 | 142.2 | 144.0 | 143.5 |
| 119.9 | 142.0 | 145.2 | 143.0 |
| 124.9 | 142.9 | 147.0 | 143.1 |
| 127.4 | 141.6 | 147.8 | 141.6 |
| 130.5 | 141.5 | 146.1 | 141.7 |

- **The runs agree until ~117 Myr.** All three stay within ~1 M km² of each other through
  117 Myr. They share a ~5 M km² drop at 114.9 Myr, so that drop isn't from crust handling.
- **Then `main` climbs ~4 M km² while #314 and the transfer stay together.** A likely
  cause, not yet tested directly: #319's interim treatment subducts the lower plate's Hm
  but piles its crust on the lower plate's own margin. That gives near-sea-level margin
  cells crust with no dense mantle under it, and they emerge. #319 called that state
  temporary, and this change removes it.
- **The transfer piles crust onto highlands.** Its crust goes onto the upper plate's
  orogen and under its front, onto already-high columns: median land elevation is lower
  (1150 vs 1207 m) and p99 higher (8842 vs 8728 m). It also keeps more continental area
  (+4.6 M km²).
- **The gap is about the size of the noise.** This is one seed, and the gap opens only in
  the last 12 Myr. It is about the size of #315's run-to-run noise (~±2%, ~4 M km²), and
  the pre- and post-review transfer runs differ by up to 1 M km² at some checkpoints.
- **Open checks:**
  - loss = 0, and underthrust = 0, to see whether the partition is the lever;
  - where `main`'s late land appears;
  - two more seeds.

## Transfer at polarized fronts

1,848 fronts over 420 steps moved 1,074 M km³ of lower-plate crust (the rest of the 1,829
M km³ donated was at fronts without polarity, still on the old path):

| Share | Volume |
|---|---:|
| Scraped (incl. cover) | 806.2 M km³ |
| Underthrust attempted / placed | 240.9 / 163.2 M km³ |
| Underthrust overflow, thrust up | 77.7 M km³ |
| Placed on the upper plate by staged placement | 883.9 M km³ |
| Lost (lower crust subducted) | 26.8 M km³ |
| Handed to other overriders / no outlet | 0 / 0 |

Cost: 1.6 ms mean per front and 7 ms mean per step with fronts, against ~4.7 s per step
overall. The slowest front took 91 ms and the slowest step 109 ms, once each, under the
concurrent test load. Seeds come from each upper plate's cached KD-tree, and one cell graph
per upper plate serves all its fronts.

## Closure

Over the replay the continental-material, craton and mobile-cover ledgers' errors changed by
+1.0 × 10⁻⁶, 0 and +6.5 × 10⁻⁶ km³ (the save opens with a fixed offset in two of them, the
same in every run). The Hm source/sink ledger closed to −6.4 × 10⁻⁵ km³. On `main` it is off by
−2.6 × 10⁵ km³ over the same replay; the transfer path books polarized donors' Hm through one
call per front instead of #319's survivor-spreading path.

## Not addressed here

`boundary_retreat` still put 18.4 M km² of continental area newly at the Hm cap (main: 22.2
M km²). That is Hm placement at fronts without a polarity, which still go through the old
same-plate accretion.
