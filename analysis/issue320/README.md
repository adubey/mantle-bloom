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

| Measure | #314 | main | transfer |
|---|---:|---:|---:|
| Land area | 141.5 M km² | 146.1 M km² | 141.4 M km² |
| Continental area | 235.1 M km² | 232.4 M km² | 237.2 M km² |
| Emergent continental area | 137.8 M km² | 141.5 M km² | 137.6 M km² |
| Land elevation p10 / p50 / p90 / p99 | 183 / 1086 / 6796 / 8043 m | 197 / 1207 / 6779 / 8728 m | 181 / 1172 / 6739 / 8875 m |
| Continental Hc p50 / p90 / p99 | 33.3 / 77.6 / 84.0 km | 33.6 / 74.7 / 84.0 km | 33.5 / 70.5 / 84.0 km |
| Continental Hm-cap fraction (change) | 13.62% (+7.67 pp) | 6.15% (+0.20 pp) | 5.76% (−0.19 pp) |
| Continental material change | — | −0.36% | −0.46% |
| Collision loss of continental material | — | 0.0021%/Myr | 0.0077%/Myr |
| Runtime / peak RSS | — | 2,467 s / 2,692 MiB | 2,253 s / 2,639 MiB |

Collision loss is delamination + no-outlet + the lower-crust loss share, per Myr, as a share
of the opening inventory. 0.0077%/Myr is about half of Earth's ~0.015%/Myr orogenic
delamination and inside #276's 0.01–0.05%/Myr seed range at its low end. By account
(transfer): 25.6 M km³ lower-crust loss share, 2.4 M km³ delaminated roots, 0 no-outlet.

Land and emergent area end ~4.7 M km² (3%) below `main`, back on the #314 replay's values;
continental area ends 4.8 M km² above `main`. Run-to-run noise on continental area is about
±2% at 43 Myr (#315). `main` keeps consumed crust on the lower plate's own survivors, which
thickens it into land behind the suture; the transfer moves it onto the upper plate's
frontal belt and under its front, where part of it goes into the root rather than the
surface.

## Transfer at polarized fronts

1,948 fronts over 424 steps moved 1,082 M km³ of lower-plate crust (the rest of the 1,836
M km³ donated was at fronts without polarity, still on the old path):

| Share | Volume |
|---|---:|
| Scraped (incl. cover) | 809.6 M km³ |
| Underthrust attempted / placed | 245.3 / 170.3 M km³ |
| Underthrust overflow, thrust up | 75.0 M km³ |
| Placed on the upper plate by staged placement | 884.6 M km³ |
| Lost (lower crust subducted) | 27.3 M km³ |
| Handed to other overriders / no outlet | 0 / 0 |

Cost: 2.0 ms mean and 11 ms max per front; 8.4 ms median and 29 ms max per step, against
~4.7 s per step overall. Seeds come from each upper plate's cached KD-tree, and one cell graph
serves the underthrust and the staged placement.

## Closure

Over the replay the continental-material, craton and mobile-cover ledgers' errors changed by
−1.0 × 10⁻⁶, 0 and +6.4 × 10⁻⁷ km³ (the save opens with a fixed offset in two of them, the
same in both runs). The Hm source/sink ledger closed to −7.7 × 10⁻⁵ km³. On `main` it is off by
−2.6 × 10⁵ km³ over the same replay; the transfer path books polarized donors' Hm through one
call per front instead of #319's survivor-spreading path.

## Not addressed here

`boundary_retreat` still put 19.0 M km² of continental area newly at the Hm cap (main: 22.2
M km²). That is Hm placement at fronts without a polarity, which still go through the old
same-plate accretion.
