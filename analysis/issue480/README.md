# Issue #480: production cut-routing measurements

These are diagnostic results for the implementation in progress. They do **not** pass
issue #480's 50% per-seed acceptance gate, and should not be presented as a completed
acceptance replay.

## Method

The production-path probe generated cratons-first worlds for seeds 7, 11 and 2024, with
craton formation disabled, at 2 Myr steps to 1 Gyr. It exercised the actual min-cut in
`maybe_split_plate` and recorded its production counters. Separately, the probe's existing
`prototype_mincut_split._evaluate` shadow calculator measured routed-vs-flow results on each
flow cut in the same evolving world. The line-crossing reductions in the table below are
from that shadow calculator, not a second production counter. The full probe outputs are in
`/tmp/issue480-production-1gyr.json`; they are not checked in.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
  --layout cratons-first --seeds 7,11,2024 --years 1e9 --step 2e6 \
  --every 100e6 --out /tmp/issue480-production-1gyr.json
```

The reference calculator used the RFC defaults: 1,000 km corridor and 0.05 cost floor. The
same-state flow comparison isolates the cut routing but does not model how a full flow-only
replay would diverge later. The earlier RFC 0009 G3 comparison is in
[`analysis/rfc0009/g3.json`](../rfc0009/g3.json), from model commit
`5339d0ca0ee29986852da1f1c33f750a2cb8a9a4`. This implementation branch starts at
`68f8bf3`; #477's intervening initial-relief changes mean the older G3 run is not a strict
same-commit baseline.

## Results

| Seed | Same-state flow line crossings → routed | Line-crossing reduction | Band craton share | Split line crossings / splits | Ledger error (m³) |
|---:|---:|---:|---:|---:|---:|
| 7 | 51 → 45 | 11.8% | 19.8% | 29 / 76 | -1,536 |
| 11 | 60 → 48 | 20.0% | 16.1% | 17 / 59 | 1,024 |
| 2024 | 42 → 30 | 28.6% | 11.9% | 18 / 72 | 1,536 |

The line-crossing reduction is against the flow line on the same events, among cuts whose
flow band crossed craton. Band craton share is the area-weighted share across actual split
and failed-rift bands. The actual production `split_line_cuts_through_craton / split_cuts`
ratios were 38.2%, 28.8% and 25.0%. For comparison, RFC 0009 G3's historical band shares
were 32.4%, 35.3% and 30.6%, respectively. The current band shares are about 39%, 54% and
61% lower; seed 7 misses 50%. The shadow line-crossing reduction misses 50% for all three
seeds. A no-routing baseline was not rerun at branch-base commit `68f8bf3`, so these numbers
are diagnostic and do not close issue #480's matched-baseline gate.

## Min-cut rejection and admissibility margins

The earlier table reports the production-probe aggregate and the same-state line-crossing
comparison. A separate margin-instrumented run evaluated admissibility only for candidate
flow cuts whose measured band contained craton. It used seeds 7, 11 and 2024, cratons-first,
1 Gyr, 2 Myr steps, a 1,000 km corridor and cost floor 0.05. Its attempt counts differ from
the earlier production-probe aggregate despite the nominally matching controls; keep the
series separate until that discrepancy is diagnosed. Per-cut outputs are in `/tmp` and are
not checked in; the compact margins are in
[`mincut-gate-margins-1gyr.json`](mincut-gate-margins-1gyr.json).

| Seed | Flow-band craton candidates | Admissible | Rejected | Size failures | Residual failures | Pole failures | Connectivity failures | B admitted | Line-clear / admitted | Zero-band / admitted |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 64 | 29 | 35 | 27 | 9 | 0 | 1 | 0.690 | 13 / 29 | 3 / 29 |
| 11 | 41 | 29 | 12 | 11 | 5 | 1 | 0 | 0.747 | 9 / 29 | 0 / 29 |
| 2024 | 32 | 20 | 12 | 6 | 3 | 4 | 0 | 0.643 | 8 / 20 | 0 / 20 |
| **Pooled** | **137** | **78** | **59** | **44** | **17** | **5** | **1** | **0.703** | **30 / 78** | **3 / 78** |

Gate-failure columns overlap: one rejected candidate may fail several checks. Size was the
leading reason. Among size failures, the short side was short by 12, 101 and 76.5 nodes at
the median for seeds 7, 11 and 2024, respectively. Residual failures exceeded the 1.10 limit
by modest amounts (seed maxima 1.175, 1.190 and 1.121); the one connectivity failure had a
smallest-side largest-component share of 0.726. B admitted is routed craton area in the
2.5-spacing band divided by flow-cut craton area for admitted candidates only. It leaves
70.3% of that exposure pooled (29.7% reduction). Forty-eight of 78 admitted lines still
touched craton; 27 of the 30 line-clear cuts still had craton in the wider band.

## Short tuning checks

Seed 7 was run to 100 Myr at the RFC settings shown below. These checks are too short to
select new production constants.

| Corridor | Cost floor | L | B |
|---:|---:|---:|---:|
| 1,000 km | 0.05 | 63.6% | 0.972 |
| 1,500 km | 0.05 | 54.5% | 0.907 |
| 2,500 km | 0.05 | 27.3% | 0.904 |
| 1,000 km | 0.25 | 61.5% | 0.866 |

The exploratory 1,500 and 2,500 km corridors and 0.25 floor did not reach the required
50% band reduction. The 1,000 km, 0.05 defaults remain in place pending a better-supported
change.

A matched 100 Myr seed 7 cratons-first probe at 2 Myr steps and a 1,000 km corridor also
compared floors 0.01 and 0.025 against 0.05. All three had 12 flow-band
craton-crossing cuts. B was 0.9191, 0.9182 and 0.9188 respectively, or about 8.1% less
band craton area than the flow cuts. Lowering the floor therefore did not approach the
50% band-reduction target; 0.025 was only marginally better than 0.05, and 0.01 was
marginally worse. The line-clear admissible share L was 25.0%, 16.7% and 16.7%. These
short diagnostic runs are not the matched 1 Gyr acceptance replay. Compact results are in
[`floor-sweep-100myr.json`](floor-sweep-100myr.json); full outputs remain in `/tmp`.

### Cut line by crust type

The earlier coastline-overlap check counted only links whose endpoints had different
effective continental classifications; that did not answer whether a rift line runs over
oceanic or continental crust. The corrected probe classifies every cut-crossing link as
continental-continental, oceanic-oceanic, or mixed, then measures each class's share of cut
length. In the matched seed 7, 100 Myr run at 2 Myr steps, all 16 cut events were on plates
with both crust types. The selected cut was oceanic-oceanic over **32.2%** of pooled cut
length (10 of 16 events had at least one such link); the flow cut was oceanic-oceanic over
**34.4%**. Two selected cuts were entirely oceanic-oceanic links. So yes, rifts do occur over
oceanic crust in this run; routing slightly reduced, rather than increased, the oceanic share.
This is expected to be possible because the cut has no continent-only restriction. The result
is a short diagnostic, not the matched 1 Gyr acceptance replay. Compact results are in
[`cut-crust-types-100myr.json`](cut-crust-types-100myr.json).

## Rift-loss resistance

At the original 0.75/0.75 resistances, seed 456862899 was replayed for 3 Gyr at 5 Myr
steps with craton formation set to 200 Myr. Loss across all three rifting paths
(`rifted_divergent_m3`, `rifted_stretch_m3`, and `rifted_failed_m3`) averaged 11.85%, 11.32%
and 9.83% of live craton per Gyr in the first, second and third Gyr. Stretch donation was the
largest path. At 0.97/0.97, the same replay averaged 2.45%, 2.08% and 1.79%, missing the
2% target in the first two Gyr. At 0.99/0.99, it averaged 1.61%, 1.12% and 0.97%, meeting
the per-Gyr ≤2% target for this seed. Its ledger error was 5,120 m³, wall time 37.56 minutes,
and peak RSS 1,254 MiB. This single-seed result does not establish the loss gate across
other seeds or conditions.

For a matched old-cut comparison, the same seed and schedule were rerun with the production
min-cut forced to fall back to the original flow-defined great-circle cut and with the
original 0.75/0.75 resistances. The old-flow run averaged 10.52%, 9.49% and 6.56% per Gyr.
The routed implementation at those same 0.75/0.75 resistances had averaged 11.85%, 11.32%
and 9.83%: 12.6%, 19.3% and 49.9% higher than flow-only in the three corresponding windows.
This indicates routing alone did not reduce rift loss in this seed; the 0.99 resistance
change is what brings the candidate below 2%.

Comparing old-flow/0.75 with routed/0.99 combines the algorithm and resistance changes. The
observed reductions are 84.7%, 88.2% and 85.3% per corresponding Gyr window. These are
matched-seed full-run comparisons, not paired event-by-event counterfactuals: after the first
different cuts, the world trajectories diverge. Within the 0.99 trajectory itself, loss fell
30.5% from Gyr 1 to 2 and a further 13.7% from Gyr 2 to 3 (40.0% from Gyr 1 to 3). Those
period changes describe evolution over time, not an algorithm effect.

The interval rate is calculated from rifted-volume change divided by the average live craton
volume at the interval endpoints, annualized to one Gyr. Full aggregates are in
[`comparison.json`](comparison.json); the old-flow output is in
`/tmp/issue480-seed456-3gyr-oldflow.json` and is not checked in.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py \
  --layout cratons-first --seeds 456862899 --years 3e9 --step 5e6 \
  --every 100e6 --formation-myr 200 --out /tmp/issue480-seed456-3gyr-r099.json
```

The full replay output remains in `/tmp/issue480-seed456-3gyr-r099.json` and is not checked
in; the compact window aggregates and resource measurements are in
[`comparison.json`](comparison.json).

## Run notes

The 1 Gyr runs took 26.8–30.5 minutes each, in series. Existing glacier-depth warnings
appeared during the replays; all simulations completed. `git diff --check` and Ruff passed.
Unit and stress tests have not been run. The line-crossing and band-area cut gates remain
unmet, so this work is not ready to claim issue #480 acceptance.
