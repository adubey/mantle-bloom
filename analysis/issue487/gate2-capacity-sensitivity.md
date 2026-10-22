# Gate 2 capacity assumption sensitivity

## Result

The 74 connected-source candidates were re-run on the pinned Gate 1 source to capture exact
cut-band node loads and state because the prior replay artifact retained summaries only. The
capture matched all 380 production cuts and the ordered candidate time/plate identities for seeds
7, 11 and 2024. All 74 eligible candidates passed the existing solver and mesh criteria, yielding
22,896 accepted node samples.

The post-processing study changes one capacity assumption at a time using the endpoint values from
the already documented factorial grid. The reference is `K0=0.8`, pore-pressure ratio 0.42,
1330 C base temperature, strain rate `1e-15 s-1`, 20 MPa cohesion, 30 degree friction, plagioclase
lower crust and dry olivine mantle. The stored Moho-lag comparison is evaluated separately at that
reference. These equally designed scenarios are sensitivity probes, not likelihood weights or a
joint uncertainty distribution.

| Capacity scenario | Nodes above capacity | Share of 22,896 | Nodes changing classification vs reference |
|---|---:|---:|---:|
| Reference, steady profile | 1,746 | 7.63% | — |
| Reference, stored Moho lag interpolated | 1,738 | 7.59% | 0.03% |
| Base temperature 1,250 C | 1,485 | 6.49% | 1.14% |
| Base temperature 1,400 C | 1,954 | 8.53% | 0.91% |
| Strain rate `1e-16 s-1` | 2,094 | 9.15% | 1.52% |
| Strain rate `1e-14 s-1` | 1,495 | 6.53% | 1.10% |
| Pore-pressure ratio 0.0 | 1,123 | 4.90% | 2.74% |
| `K0=0.6` | 2,854 | 12.47% | 4.84% |
| `K0=1.0` | 1,017 | 4.44% | 3.18% |
| Cohesion 1 MPa | 2,059 | 8.99% | 1.37% |
| Friction `atan(0.6)` | 1,718 | 7.50% | 0.12% |
| Mafic-granulite lower crust | 1,768 | 7.72% | 0.10% |
| Wet-olivine mantle | 2,448 | 10.69% | 3.07% |

The rank ordering of `N_ext/N_cap` is highly stable: Spearman correlation with the reference is
at least 0.9938 in every scenario. The binary threshold result is less stable. Across this finite
one-factor envelope, 1,009 nodes (4.41%) exceed capacity under every tested scenario, and 3,016
(13.17%) exceed capacity under at least one. These are OAT bounds only; interacting weak settings
could produce a different joint envelope.

`K0` has the largest single-factor effect: the over-capacity fraction ranges from 4.44% at `K0=1`
to 12.47% at `K0=0.6`. The wet-olivine mantle case gives 10.69%; removing the pore-pressure
reduction gives 4.90%. The Moho-lag profile changes only eight nodes, consistent with the stored
lag being near zero for most samples. The result is conditional on the current stress-state and
material assumptions, especially the missing spatial `K0`, mantle water and pore-pressure fields.

In the reference case, continental nodes have a 1.68% exceedance fraction, while mature-oceanic
and active-divergence young-ridge-proxy nodes have 22.71% and 22.35%, respectively. The latter is
not an ocean-age result: active divergence remains a coarse oceanic-age proxy. Seed 11 has the
highest reference exceedance fraction (11.87%), compared with 3.53% for seed 7 and 5.29% for seed
2024. These are screening summaries, not observed rupture rates or labels for fitting a threshold.

## Reproduction

The capture-only replay writes per-event JSONL node samples outside the repository. Run it on the
pinned Gate 1 archive with the diagnostic force helpers overlaid:

```bash
PYTHONHASHSEED=0 PYTHONPATH=/path/to/gate1/backend backend/.venv/bin/python \
  /path/to/gate1/bin/debug/replay_gate2_capacity.py --seed 7 --years 1e9 \
  --out /tmp/g2-seed7-capture.json \
  --node-samples-out /tmp/g2-seed7-samples.jsonl --capture-only
```

Repeat for seeds 11 and 2024, then run `bin/debug/sensitivity_gate2_capacity.py` with all three
sample files, capture reports and the accepted Gate 1 replay reports. The script verifies cut
counts, candidate identities, connected-source coverage, per-node array lengths and the expected
22,896-sample population before scoring.

Compact results and source hashes are in
[`gate2-capacity-sensitivity.json`](gate2-capacity-sensitivity.json). Full per-node capture files
remain under `/private/tmp/` and are not checked into the repository.
