# RFC 0007 tuning and sensitivity runs (issue #584)

Harness: [tune_plumes_584.py](../../bin/debug/tune_plumes_584.py).
One world at a time. Baselines seeds 7/11/2024 to 1 Gyr, regenerated seed
456862899 to 3 Gyr, single-knob sensitivity cases on seed 7 (500 Myr).
Controls: default `generate_world`, cratons on, 5 Myr steps.
Starting knobs are the code defaults
(`PLUME_HEAD_RATE_PER_MYR` 0.1, log-uniform volumes, `PLUME_HEAT_MYR` 20,
`PLUME_CRATON_WEAKENING` 0.75, `PLUME_HEAD_STRESS` 1.0, docking thresholds
5/15 km, `DOCK_WIDTH_M` 3 km, `DOCK_LOWER_CRUST_LOSS` 0.2,
`PLATEAU_PRESERVED_FRACTION` 0.9).
Raw per-case JSON: [tuning.json](tuning.json) (21 entries).

Determinism: re-running baseline_7 reproduces the landed result bit-for-bit
(119 LIPs, 45 splits, identical ledgers), so the deterministic head schedule
(#579) holds and knob effects are shape-only.

## Baselines (starting knobs)

| case | LIPs/Gyr | vol/Gyr | splits |
| ---- | -------- | ------- | ------ |
| 7 | 119 | 1.6e9 km³ | 45 |
| 11 | 110 | 1.3e9 km³ | 47 |
| 2024 | 93 | 0.9e9 km³ | 67 |
| regen 3 Gyr | 95 | 1.2e9 km³ | 90/Gyr |

Head counts sit inside the Poisson(100) 95% interval on all seeds, and the
3 Gyr regen matches the 1 Gyr rate, so the rate is stable over deep time.

Realized head-volume deciles (min, p50, max) in km³: default log-uniform
~1.1e6 / ~7e6 / ~4.8e7; lognormal 1.0e6 / 6.4e6 / 5.0e7; power-law (α=1.5)
1.1e6 / 3.6e6 / 4.2e7. Variants share support, so spread is shape-only (<2x
in total volume).

Active hotspots (heads + tails) at every 100 Myr checkpoint stay 24–35,
inside the RFC's 20–50 band (t=0 reads 0 by construction: no plume has been
emplaced yet). Buoyant-ledger `unattributed_m3` is 0.0 on every case: the
ledger closes.

## Sensitivity spreads (seed 7, 500 Myr)

| knob | LIPs/Gyr | splits | heated frac |
| ---- | -------- | ------ | ----------- |
| rate 0.05 / 0.2 | 62 / 238 | 25 / 29 | 0.04 / 0.17 |
| Cox on | 60 | 22 | 0.07 |
| lognormal / power-law | 124 / 124 | 28 / 24 | 0.09 / 0.08 |
| dock easy / hard | 124 / 124 | 26 / 24 | 0.08 / 0.12 |
| loss 0.4 / plateau 0.5 | 124 / 124 | 26 / 26 | 0.12 / 0.13 |
| heat 5 / 50 Myr | 124 / 124 | 27 / 26 | 0.05 / 0.16 |
| weaken 0.5 / 1.0 | 124 / 124 | 25 / 28 | 0.10 / 0.09 |
| stress 0→3 | 124 | 26–27 | 0.06–0.13 |

LIP counts are identical across same-seed cases (deterministic schedule);
rate and Cox move counts linearly. Split counts barely move under any knob:
25–29 except the 2024 baseline outlier (67). A majority of split parents is
fully unheated on every case.

## Juvenile continental by source (km³/Gyr)

| case | flood basalt | hotspot | plateau |
| ---- | ------------ | ------- | ------- |
| baselines | 2.1–6.0e8 | 4.9–8.8e8 | 1.3e4–2.1e5 |
| sens max | 1.3e9 (rate 0.2) | 1.4e9 (rate 0.2) | 1.7e6 (rate 0.2) |

Flood basalt tracks continental head volume; hotspot tails dominate the
juvenile budget. Converted plateau is 3–4 orders of magnitude below both on
every case (see falsifiers).

## Falsification checks

| check | result | verdict |
| ----- | ------ | ------- |
| dock share ≥10% | 0.53–0.75 everywhere | pass |
| converted spread <5x | 16x baselines, 668x all | fail |
| converted 0.3–1.5 km³/yr | best case 1.7e-3 | fail |
| rift concentration, 5–50 Myr | heated frac ≤0.16, unheated majority | weak |
| heads → splits | 26–27 splits for stress 0–3 | no effect |

Docking works (about 70% of trench-arriving plateau volume docks rather than
subducts), but almost none of it converts: `_convert_buoyant_on_continental`
fires so rarely that the plateau→continent pathway contributes ~1e-4 km³/yr
against the 0.3–1.5 km³/yr sanity range, and no knob comes within 2 orders of
magnitude. The stress sweep is flat and agrees with the #582 comparator
(head-following 0.19 vs matched plume-off 0.52): nothing supports the
half-of-heads breakup target either.

## Knob decisions

Keep every starting knob: rate 0.1, log-uniform volumes, docking thresholds,
loss 0.2, preserved 0.9, heat 20 Myr, weakening 0.75, stress 1.0. No sweep
gives a reason to move any of them, and the flat responses mean the values
are load-bearing only through the pathways that already work.

Two RFC revisions are proposed (not applied here):
1. Replace the 0.3–1.5 km³/yr converted-plateau sanity range or rework the
   conversion pathway — the mechanism as built under-produces by ~4000x.
2. Drop the half-of-heads breakup target and fix the knob table's
   `PLUME_HEAD_STRESS` 1.5 to the code's 1.0.

Step-partition robustness (1 vs 10 Myr steps) was not run: outside #584's
scope, still open against the RFC acceptance list.
