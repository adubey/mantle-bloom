# RFC 0015 D8 and G6

## Replay plan

Compare D8 with the pre-D8 main commit `57651bf2425d8d8b71a2a3230a8814ae74f7d7bd`.
This baseline includes the lake tier-cap fix in PR #1030, which landed after D8 began.
For each arm record the exact implementation commit and save SHA-256. Use
`bin/debug/replay_d8_meltwater.py` with the same save, controls and step schedule in both
arms. The diagnosis saves are seeds 997271774 and 649658380 with ice ages off. Replay each
for 10 Myr at 0.25 Myr and 1 Myr steps. Report seed 374464662 with ice ages at an 8 Myr
period for 18 Myr at both step sizes, and with a 1 Myr period at 0.25 Myr steps. The ice-age
arms do not gate the RFC. The save inventory and original hashes are in
[`analysis/issue548`](../issue548/README.md). Converted v5 saves have different checksums;
record their checksums as well as their original hashes. The script enables the glacier
march in both arms and records the saved march controls and ice-age period.

The converted v5 inputs have these SHA-256 checksums. The private B2 catalog links each
to the corresponding original checksum recorded by #548.

| Set | Seed | Converted save SHA-256 |
|---|---:|---|
| Diagnosis | 997271774 | `07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8` |
| Diagnosis | 649658380 | `06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f` |
| Diagnosis, ice ages | 374464662 | `a9eed5dae3b685c3ed0350a1035bf3adbc66afbe958b64e639139583da52801a` |
| Held out, ice ages | 447276134 | `22b60fee062d10eea11838b3006beed5fd228e4e22b28e0d556602ff7b99a30d` |
| Held out | 628977350 | `e5acdb30a642b1c6d3cccc1480874acf807a16ff0b5258ef9958ee6e89af2303` |
| Held out | 295569319 | `1af7fcb73f6b7c727fb36f2037af712088bae692509b7842b08c6a6b9e4d9c27` |

Each step records the two source budgets, their residuals, filled hydraulic depressions,
their throughput, unresolved drainage, surface melt routed from thick warm-based ice, and
the count and ten largest discharges of ice-margin river networks. Their headwaters must
be on an ice cell adjacent to ice-free ground, or on ice-free ground adjacent to ice. A filled pit is a
warm-based local-minimum plateau whose hydraulic head was raised. Throughput is summed on
entry to each connected region of raised head. A region can contain multiple pits. The
baseline has no D8 source budgets or filled-pit diagnostics, so these values
are unavailable there. The networks are computed on both arms.

The expected direction is more glacier water reaching a margin, but no size of the change
is assumed. Routing, precipitation and river storage can change together over a replay.
Report the full trajectories and runtime, rather than infer a causal effect from one final
value. No variability estimate is available for these measures before the replay. G6 is a
reporting gate. After G6, evaluate the fixed candidate on held-out seeds 628977350 and
295569319, with the same 10 Myr and two-step schedule. The RFC requires each source budget
to close to round-off and at least half of the ice-derived source to reach a lake or ocean.
Below one fifth on a held-out save falsifies the routing premise. Seed 447276134 is the
report-only ice-age held-out arm.

Example for each arm, substituting the pinned backend and local save:

```bash
backend/.venv/bin/python bin/debug/replay_d8_meltwater.py <save> --backend <backend> \
  --myr 10 --step-years 250000 --checkpoint-myr 1 --out <output.json>
```

## Results

The reviewed D8 head `0cf10cd596a63d3db5dc1703271d90c941620064` and the pinned baseline
were replayed sequentially from the converted seed 649658380 save for 10 Myr at 0.25 Myr
steps. The corrected per-step data and ten largest network discharges are in
[`baseline`](g6-baseline-649658380-025-corrected.json) and
[`D8`](g6-d8-649658380-025-reviewed.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 75–129 | 70–114 |
| Ice-margin networks, final step | 97 | 88 |
| Sum of step wall times | 1,684 s | 1,605 s |
| Peak RSS | 1,041 MB | 1,074 MB |

The D8 glacier source (ice release and rain on ice) totalled 3.112 × 10¹⁸ m³; 99.55%
of that total was delivered to lakes or the ocean during the steps, and the minimum
per-step delivery share was 98.69%. The bare-snow source totalled 1.767 × 10²⁰ m³,
with 96.17% delivered to lakes or the ocean in aggregate and at least 95.17% in each
step. The largest absolute source-ledger residual relative to its input was 2.09 × 10⁻¹⁶
for glacier water and 3.23 × 10⁻¹⁶ for bare snowmelt. Unresolved glacier drainage totalled
1.39% of glacier-source input; the highest per-step fraction was 7.68%. This water remains
in ground routing, so the unresolved amount is a diagnostic and not an additional ledger
destination. The maximum number of filled pits in a step was 15; surface melt on routed
warm-based ice thicker than 1 km totalled 1.124 × 10¹⁶ m³.

The network counts differ, but this single trajectory does not isolate routing from changes
in lake, river and ice evolution or establish variability. The other diagnosis seed,
ice-age reporting arms and the held-out gate remain open; these results do not
complete G6. Both runs were sequential and the machine's budget check reported no other
heavy runs at the start of each arm.

The same seed was then replayed for 10 Myr at 1 Myr steps, using the same baseline and
D8 code (the D8 head `ffb54ee2e6f1f09714a1513bdfef22b0a9266f0f` adds only the preceding
analysis). Per-step results are in [`baseline`](g6-baseline-649658380-100.json) and
[`D8`](g6-d8-649658380-100.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 78–109 | 75–101 |
| Ice-margin networks, final step | 94 | 87 |
| Sum of step wall times | 1,090 s | 1,126 s |
| Peak RSS | 1,033 MB | 1,025 MB |

D8's 3.849 × 10¹⁸ m³ glacier source delivered 99.33% to lakes or ocean in aggregate
(at least 98.49% per step). The 1.767 × 10²⁰ m³ bare-snow source delivered 95.52%
(at least 93.85% per step). Their maximum absolute relative residuals were 1.68 × 10⁻¹⁶
and 4.35 × 10⁻¹⁶ respectively. Unresolved drainage was 3.79% of total glacier source,
with a maximum of 12.22% in one step; it remained in ground routing. The maximum filled
pit count was 21, and warm-based surface melt on ice thicker than 1 km totalled
5.854 × 10¹⁶ m³. On the second D8 step two nodes reached the 5 km depth guard; the solver
reported zero shed meltwater. The two step sizes produce different warm-thick melt totals
and unresolved shares, and further seeds are needed before attributing those differences.
The arms ran sequentially with no other heavy run reported at either start.

Seed 997271774 was replayed next from the converted 87.4 Myr save with the same 0.25 Myr
steps. Its SHA-256 is `07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8`.
Per-step results are in [`baseline`](g6-baseline-997271774-025.json) and
[`D8`](g6-d8-997271774-025.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 203–279 | 176–284 |
| Ice-margin networks, final step | 223 | 262 |
| Sum of step wall times | 2,069 s | 2,138 s |
| Peak RSS | 1,432 MB | 1,293 MB |

D8's 6.730 × 10¹⁹ m³ glacier source delivered 89.33% to lakes or ocean in aggregate,
with at least 83.31% in any step. The 2.118 × 10²⁰ m³ bare-snow source delivered 92.40%
in aggregate and at least 90.02% per step. Maximum absolute relative source-ledger residuals
were 3.32 × 10⁻¹⁶ and 3.25 × 10⁻¹⁶. Unresolved drainage was 2.31% of total glacier input,
at most 3.73% in a step; it remained in ground routing. This is a second diagnosis save,
but the timestep and ice-age reporting arms and held-out gate remain open.

The same save was also replayed at 1 Myr steps. Per-step outputs are in
[`baseline`](g6-baseline-997271774-100.json) and [`D8`](g6-d8-997271774-100.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 223–282 | 214–285 |
| Ice-margin networks, final step | 273 | 252 |
| Sum of step wall times | 1,722 s | 1,694 s |
| Peak RSS | 1,323 MB | 1,387 MB |

D8's 6.363 × 10¹⁹ m³ glacier source delivered 91.48% to lakes or ocean in aggregate
(at least 86.67% per step). The 2.082 × 10²⁰ m³ bare-snow source delivered 93.25%
(at least 92.08% per step). Maximum absolute relative source-ledger residuals were
1.87 × 10⁻¹⁶ and 2.44 × 10⁻¹⁶. Unresolved glacier drainage was 5.29% of total input,
with at most 14.12% in one step; unresolved water remains in ground routing. This step
size also reported repeated nodes at the 5 km ice-depth guard in both baseline and D8;
on D8 one step shed about 1,690 km³ as meltwater while other reported steps shed zero.
The guard behaviour warrants review alongside the remaining G6 arms.

The first baseline arm on seed 649658380 at 0.25 Myr steps is superseded: review found its
network diagnostic included headwaters in the interior of the ice sheet. Its 40-step output
is preserved in [`g6-baseline-649658380-025.json`](g6-baseline-649658380-025.json), but is
not G6 evidence. The corrected baseline and comparison are reported above.
The corresponding D8 arm was started and interrupted at the owner's request for reviews;
it wrote no completed checkpoints. A one-step development smoke replay on seed
649658380 at 0.25 Myr used branch commit `809d123351e727f5789648c368c0a05817f4894e` and
converted save SHA-256 `06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f`.
The glacier and bare-snow source residuals were 16 m³ of 1.88 × 10¹⁷ m³ and −650 m³ of
4.13 × 10¹⁸ m³, respectively. No drainage was unresolved. This development run used the
earlier connected-region count for pits and is not a G6 measurement.

An initial baseline replay from `eb3cbc4febbfa7a2df51c30d01916ab091a9f647` on seed
649658380 at 0.25 Myr steps was stopped after 1.5 Myr when #1030 landed. Its six completed
steps remain in the local run output. They are not used in the comparison because the lake
tier-cap fix changes lake routing. The replay plan above pins the updated baseline before
the comparison is restarted.
