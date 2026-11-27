# RFC 0015 scale check

This report compares the four RFC 0015 D8 replay saves with the Last Glacial Maximum
(LGM) reference of about 52 million km³ of extra grounded ice and 134 m of sea-level
lowering [[Lambeck et al., 2014]](../../rfc/0015-frozen-water-in-physical-units.md#ref-lambeck2014).
The RFC criterion is report-only: the comparison measures scale and does not set a new
pass/fail threshold.

## Method

Each save was replayed for 10 Myr in 1 Myr steps on `origin/main` commit
`c4584b20b42b353272bd958ad90c37565a680c49`, with ice ages disabled and the glacier march
enabled. The Rust glacier solver was available for every run. The reproducing diagnostic
is [`bin/debug/measure_rfc0015_scale.py`](../../bin/debug/measure_rfc0015_scale.py); its
compact checkpoint outputs are alongside this report. The saves are local and not
available for download; their names and SHA-256 checksums are:

| Seed | Opening age (Myr) | Save | SHA-256 |
|---:|---:|---|---|
| 649658380 | 6.4 | `mantle-bloom-seed649658380-6400000y.bloom` | `06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f` |
| 997271774 | 87.4 | `mantle-bloom-seed997271774-87400000y.bloom` | `07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8` |
| 628977350 | 226.7 | `mantle-bloom-seed628977350-226700000y.bloom` | `e5acdb30a642b1c6d3cccc1480874acf807a16ff0b5258ef9958ee6e89af2303` |
| 295569319 | 180.1 | `mantle-bloom-seed295569319-180000000y.bloom` | `1af7fcb73f6b7c727fb36f2037af712088bae692509b7842b08c6a6b9e4d9c27` |

Grounded volume is the sum of grounded ice thickness times each cell's accounting area.
The volume percentage is relative to the LGM's *extra* 52 million km³, not an anomaly
from a modeled present-day baseline. Flat-ocean equivalent converts grounded ice to
water-equivalent depth with `RHO_ICE / RHO_WATER`, then divides by connected ocean area.
It is a scale comparison, not the simulated sea level. The latter is shown separately
because topography, lakes, and water inventory also affect it.

Visible ice is depth ≥ `GLACIER_VISIBLE_DEPTH_M` (10 m). Margin ice is visible ice with at
least one neighboring node below that threshold; interior ice is the rest. Speed quantiles
are accounting-area-weighted quantiles of the model's calibrated mass-balance-flux sliding
speed estimate (`hydro.ice.sliding_speed_m_per_yr`). This is not velocity integrated by a
physical ice-flow solver. The maximum is a node maximum, not an area-weighted statistic.

## Endpoint comparison

| Seed | Endpoint age (Myr) | Grounded ice (million km³) | LGM volume | Flat-ocean equivalent (m) | LGM sea-level scale | Solved sea level (m) | Ocean area (million km²) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 649658380 | 16.4 | 2.12 | 4.1% | 5.6 | 4.2% | −104.0 | 348.7 |
| 997271774 | 97.4 | 7.60 | 14.6% | 19.5 | 14.6% | −523.8 | 357.0 |
| 628977350 | 236.7 | 64.39 | 123.8% | 138.1 | 103.1% | +5,101.0 | 427.4 |
| 295569319 | 190.1 | 28.71 | 55.2% | 75.8 | 56.5% | −262.9 | 347.5 |

Two endpoint states are below the LGM volume scale (2.12 and 7.60 million km³); the
two held-out states are above it (28.71 and 64.39 million km³). The high-volume held-out
seed's flat-ocean equivalent is close to the LGM's 134 m, while the other is about 57% of
that level. Their solved sea levels differ sharply from those equivalents, including
+5.1 km on seed 628977350. That solved value is dominated by this model world's state and
is not an Earth-like sea-level estimate.

The 1 Myr trajectories are variable, especially on the held-out saves:

| Seed | Grounded ice range over run (million km³) | Range as % of LGM volume |
|---:|---:|---:|
| 649658380 | 2.12–3.56 | 4.1–6.8% |
| 997271774 | 5.65–13.19 | 10.9–25.4% |
| 628977350 | 18.82–65.82 | 36.2–126.6% |
| 295569319 | 7.75–29.57 | 14.9–56.9% |

Both held-out runs jump well above their saved opening volume during the first simulated
Myr. The two high-ice trajectories repeatedly report glacier-depth-cap and drainage-route
warnings; these diagnostics limit interpretation of those large values. This report
therefore records the scale outcome and its limits, rather than claiming the model has
matched LGM climate or ice-sheet dynamics.

## Ice extent and speed distribution at endpoints

| Seed | Visible ice (million km²) | Margin share | Interior speed p90 / p95 / max (m/yr) | Margin speed p90 / p95 / max (m/yr) |
|---:|---:|---:|---:|---:|
| 649658380 | 3.00 | 36.2% | 57.9 / 120.2 / 447.9 | 0.0 / 0.0 / 318.0 |
| 997271774 | 12.90 | 35.7% | 80.9 / 142.0 / 1,122.0 | 0.0 / 0.0 / 3,487.3 |
| 628977350 | 71.51 | 30.0% | 77.7 / 147.1 / 3,839.2 | 365.7 / 633.4 / 6,298.1 |
| 295569319 | 37.39 | 20.1% | 54.7 / 115.4 / 1,749.4 | 0.0 / 0.0 / 4,245.1 |

The high maximum estimates are concentrated at a few nodes and should not be read as
representative flow rates. The p90 and p95 values give a more stable description of the
interior and margin fields, but remain estimates from the model's current mass-balance
proxy.

## Reproduction

Run once per save from the repository root:

```sh
backend/.venv/bin/python bin/debug/measure_rfc0015_scale.py <save> \
  --myr 10 --step-years 1000000 --checkpoint-myr 1 --out <output.json>
```

The outputs report every 1 Myr checkpoint, the opening and closing state, ocean area,
grounded and total ice, flat-ocean equivalent, margin share, speed summaries and peak RSS.
