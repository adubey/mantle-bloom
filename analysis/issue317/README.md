# Issue #317: Hm ratchet baseline

## Reproduction

Run the saved cascade state from 87.4 Myr for 43.1 Myr (431 × 100 ky) and compare it
with the supplied 130.5 Myr reference save:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  /Users/adubey/Downloads/mantle-bloom-seed997271774-87400000y.mbworld \
  --myr 43.1 --checkpoint-myr 2.5 \
  --compare /Users/adubey/Downloads/mantle-bloom-seed997271774-130500000y.mbworld \
  --out analysis/issue317/cascade-seed997271774.json
```

The saved input SHA-256 is `08ad9461d9585e92be93a58d81e326c7e4b685ae7d462135d01b0e9931ed762b`.
The run used commit `4d0a8d2fcfac2c19d206112a2270cf4fcabd175a`. The complete 431-step,
2.5-Myr-checkpoint record is in [`cascade-seed997271774.json`](cascade-seed997271774.json).
All volume and area values below use each quad cell's exact accounting area.

## Result

| Measure | 87.4 Myr | 130.5 Myr | Change |
|---|---:|---:|---:|
| Continental area | 268.552 M km² | 235.133 M km² | −33.419 M km² |
| Continental Hm-cap area | 15.963 M km² | 32.016 M km² | +16.053 M km² |
| Continental Hm-cap fraction | 5.944% | 13.616% | +7.672 pp |
| Mean continental Hm | 107.810 km | 119.389 km | +11.579 km |
| Cratonic Hm-cap fraction | 6.984% | 11.627% | +4.643 pp |
| Non-cratonic Hm-cap fraction | 5.577% | 14.400% | +8.633 pp |

The supplied 130.5 Myr comparison save has 12.051% continental cap coverage (27.979 M km²),
so the replay endpoint is 1.565 percentage points higher. It is not an identical trajectory
endpoint: the comparison has 232.166 M km² continental area and 145.888 M km² land, versus
235.133 M km² continental area and 141.535 M km² land in the replay.

## Attribution and closure

The Hm inventory falls by 1.392 billion km³ overall. The debug-only ledger records 4.743
billion km³ of gross sources and 6.135 billion km³ of gross sinks. By effective node type:

| Scope | Live Hm change | Sources | Sinks | Ledger residual |
|---|---:|---:|---:|---:|
| All nodes | −1,391.724 M km³ | 4,743.063 M km³ | 6,134.787 M km³ | −0.000003 km³ |
| Continental nodes | −880.380 M km³ | 218.222 M km³ | 1,098.602 M km³ | +0.000005 km³ |
| Oceanic nodes | −511.343 M km³ | 4,556.193 M km³ | 5,067.537 M km³ | +0.000008 km³ |

These residuals are floating-point noise; the largest absolute per-step ledger residual was
0.000017 km³ for all nodes and 0.000012 km³ for continental nodes. The phase-budget residual
is a separate audit and is not added to the source/sink ledger. The aligned-node
reclassification column is zero for this replay; plate snapshots now preserve the owning
plate type on both sides of topology operations so inherited nodes are classified against
the correct before/after default.

There were 4,282 connected donor fronts: 7,431.277 M km³ donor Hm, 6,596.365 M km³ placed
on survivors, and 834.912 M km³ unplaced/delaminated. Convergent Hm/Hc clipping stayed at
zero throughout the replay.

## Runtime profile

The replay took 2,524.66 seconds (42.1 minutes), with 2,306.91 MiB peak RSS. Instrumented
totals were 155.54 s in boundary classification (7,564 calls), 14.62 s in spatial graph
builds (25,011), 8.19 s in spatial-index builds (19,916), 1.98 s in graph traversal (4,288),
and 1.85 s in Hm spreading (4,282).
