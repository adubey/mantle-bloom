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
value. No variability estimate is available for these measures before the replay. G6 reports the
fixed candidate on held-out seeds 628977350 and 295569319 with the same 10 Myr and two-step
schedule. The RFC requires each source budget to close to round-off and at least half of the
ice-derived source to reach a lake or ocean. Below one fifth on a held-out save falsifies the
routing premise. The prior candidate passed both held-out seeds at both step sizes; those
results are historical after the outlet fix described below. Seed 447276134 is a separate
report-only ice-age arm; it is not part of the gate.

Example for each arm, substituting the pinned backend and local save:

```bash
backend/.venv/bin/python bin/debug/replay_d8_meltwater.py <save> --backend <backend> \
  --myr 10 --step-years 250000 --checkpoint-myr 1 --out <output.json>
```

## Review-fix replay preregistration

The Codex review found that a warm cell could keep a higher outlet seed instead of the
steepest route on the filled hydraulic surface. This changes D8 routing, so the prior D8
diagnosis and held-out candidate measurements are not evidence for the revised code. The
revised D8 arm will be replayed on diagnosis seeds 649658380 and 997271774 and held-out
seeds 628977350 and 295569319, each for 10 Myr at 0.25 and 1 Myr steps. The pinned baseline
remains `57651bf2425d8d8b71a2a3230a8814ae74f7d7bd`; its existing paired records are retained.
Use the same converted saves and controls listed above, and record each revised candidate
commit, save checksum, and glacier solver backend (`rust` or `python`). Run one arm at a time.

The primary criterion remains at least 50% of glacier-source input delivered to a lake or
ocean for each held-out seed at both step sizes; below 20% falsifies the routing premise.
Both glacier and bare-snow source residuals must close to round-off at every step. The
glacier ledger now starts from ice release plus rain on ice before margin routing, so its
residual also measures loss or duplication in the ice-to-margin handoff. No numeric change
in delivery is assumed after correcting outlet selection. Report unresolved drainage
separately because water in this category may later be delivered by ground routing.

## Results

**Prior candidate only; superseded by the review-fix reruns below.** Seeds 628977350 and
295569319 exceeded the 50% glacier-source delivery criterion at both step sizes on the
earlier candidate. That evidence remains useful historically but does not close G6 for the
revised routing. The separate report-only ice-age held-out seed 447276134 was also recorded
at both step sizes on the earlier candidate.

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

Seed 374464662 is the report-only 8 Myr ice-age-period arm. Baseline and D8 were run for
18 Myr at 0.25 Myr steps from the converted 132.6 Myr save (SHA-256
`a9eed5dae3b685c3ed0350a1035bf3adbc66afbe958b64e639139583da52801a`). Per-step data are
in [`baseline`](g6-baseline-374464662-025-age8m.json) and
[`D8`](g6-d8-374464662-025-age8m.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 139–481 | 150–460 |
| Ice-margin networks, final step | 200 | 194 |
| Sum of step wall times | 4,213 s | 4,349 s |
| Peak RSS | 1,055 MB | 1,054 MB |

D8's 7.020 × 10¹⁹ m³ glacier source delivered 97.22% to lakes or ocean in aggregate
(at least 91.04% per step). The 2.621 × 10²⁰ m³ bare-snow source delivered 95.54%
(at least 94.69% per step). Both source ledgers closed within 2.80 × 10⁻¹⁶ relative
residual. Unresolved glacier drainage was 1.24% of total source and at most 7.33% in a
step. Both arms repeatedly reported a node at the 5 km depth guard with zero shed meltwater.
This is a report-only ice-age arm, not a gate result.

The same report-only save was replayed for 18 Myr at 1 Myr steps. Both arms reached
150.6 Myr. The 1 Myr D8 replay had previously stopped after 14 steps; because the replay
script does not persist a resumable world state, the D8 arm was rerun from the original
save and completed all 18 steps. Per-step data are in [`baseline`](g6-baseline-374464662-100-age8m.json)
and [`D8`](g6-d8-374464662-100-age8m.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 220–438 | 220–418 |
| Ice-margin networks, final step | 263 | 295 |
| Sum of step wall times | 2,891 s | 3,111 s |
| Peak RSS | 1,055 MB | 1,031 MB |

D8's 8.369 × 10¹⁹ m³ glacier source delivered 97.64% to lakes or ocean in aggregate
(at least 93.32% per step). The 2.473 × 10²⁰ m³ bare-snow source delivered 95.56%
(at least 94.24% per step). Maximum absolute relative source-ledger residuals were
2.87 × 10⁻¹⁶ and 2.25 × 10⁻¹⁶. Unresolved glacier drainage was 2.19% of total source,
with at most 7.93% in a step; it remains in ground routing. The maximum filled-pit
count was 73 and warm-based surface melt on ice thicker than 1 km totalled
1.304 × 10¹⁸ m³. The 5 km depth guard recurred; the 146.6 Myr step shed 120 km³ as
meltwater. These are report-only ice-age results, not a G6 gate result. The arms ran
sequentially with no other heavy run reported at either start.

The same 18 Myr ice-age save was also replayed at 0.25 Myr steps with a 1 Myr ice-age
period. Both arms reached 150.6 Myr after 72 steps. The per-step data are in
[`baseline`](g6-baseline-374464662-025-age1m.json) and
[`D8`](g6-d8-374464662-025-age1m.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 121–480 | 103–474 |
| Ice-margin networks, final step | 472 | 412 |
| Sum of step wall times | 4,102 s | 4,096 s |
| Peak RSS | 1,047 MB | 1,054 MB |

D8's 1.020 × 10²⁰ m³ glacier source delivered 97.91% to lakes or ocean in aggregate
(at least 92.06% per step). The 2.320 × 10²⁰ m³ bare-snow source delivered 95.00%
(at least 91.30% per step). Maximum absolute relative source-ledger residuals were
3.59 × 10⁻¹⁶ and 3.27 × 10⁻¹⁶. Unresolved glacier drainage was 0.56% of total source,
with at most 4.01% in a step; it remains in ground routing. The maximum filled-pit count
was 82 and warm-based surface melt on ice thicker than 1 km totalled 1.615 × 10¹⁸ m³.
One late step hit the 5 km depth guard and shed about 66,600 km³ as meltwater. These are
report-only ice-age results, not a G6 gate result. The baseline and D8 arms ran
sequentially.

The first held-out save, seed 628977350, was replayed for 10 Myr at 0.25 Myr steps. Both
arms reached 236.7 Myr after 40 steps. Per-step data are in
[`baseline`](g6-baseline-628977350-025-heldout.json) and
[`D8`](g6-d8-628977350-025-heldout.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 1,409–2,267 | 1,223–1,994 |
| Ice-margin networks, final step | 1,580 | 1,468 |
| Sum of step wall times | 4,814 s | 4,913 s |
| Peak RSS | 1,334 MB | 1,459 MB |

D8's 2.861 × 10²⁰ m³ glacier source delivered 98.64% to lakes or ocean in aggregate
(at least 97.53% per step). The 1.069 × 10¹⁹ m³ bare-snow source delivered 99.21%
(at least 98.51% per step). Maximum absolute relative source-ledger residuals were
2.79 × 10⁻¹⁶ and 2.23 × 10⁻¹⁶. Unresolved glacier drainage was 0.023% of total source,
with at most 0.38% in a step. This clears the RFC's 50% held-out delivery criterion and
is well above its 20% falsification threshold. One late step reported 65,200 km³ shed at
the 5 km depth guard. Seed 295569319 and the other step-size arm are reported below.
The baseline and D8 arms ran sequentially.

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


The same held-out seed 628977350 was also replayed for 10 Myr at 1 Myr steps. Both arms
reached 236.7 Myr after 10 steps. Per-step data are in
[`baseline`](g6-baseline-628977350-100-heldout.json) and
[`D8`](g6-d8-628977350-100-heldout.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 1,535–2,138 | 1,325–1,886 |
| Ice-margin networks, final step | 1,567 | 1,325 |
| Sum of step wall times | 3,584 s | 3,615 s |
| Peak RSS | 1,511 MB | 1,394 MB |

D8's 2.258e+20 m³ glacier source delivered 98.92% to lakes or ocean in aggregate
(at least 98.29% per step). The 1.052e+19 m³ bare-snow source delivered 99.38%
(at least 98.57% per step). Maximum absolute relative source-ledger residuals were
1.62e-16 and 2.13e-16. Unresolved glacier drainage was 0.022% of total
source, with at most 0.165% in one step. The maximum filled-pit count was
264 and warm-based surface melt on ice thicker than 1 km totalled 5.265e+15 m³.
This clears the RFC's 50% held-out delivery criterion and exceeds the 20% falsification
threshold. The final 1 Myr step reported no-route drainage and lake depth-cap overflow,
and depth-guard logs during the run included shed meltwater; ledger closure remains
round-off.

Seed 628977350 clears the held-out gate at both step sizes. Seed 295569319 results are reported below. The baseline and D8 arms ran sequentially.


Held-out seed 295569319 was replayed for 10 Myr at 0.25 Myr steps. Both arms reached
190.1 Myr after 40 steps. Per-step data are in
[`baseline`](g6-baseline-295569319-025-heldout.json) and
[`D8`](g6-d8-295569319-025-heldout.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 582–693 | 574–675 |
| Ice-margin networks, final step | 616 | 594 |
| Sum of step wall times | 2,542 s | 2,590 s |
| Peak RSS | 2,562 MB | 2,264 MB |

D8's 9.032e+19 m³ glacier source delivered 98.47% to lakes or ocean in aggregate
(at least 97.82% per step). The 2.169e+20 m³ bare-snow source delivered 95.33%
(at least 92.58% per step). Maximum absolute relative source-ledger residuals were
1.60e-16 and 3.09e-16. Unresolved glacier drainage was 1.18% of total source,
with at most 2.90% in a step. The maximum filled-pit count was 144 and
warm-based surface melt on ice thicker than 1 km totalled 4.596e+17 m³. This clears the
RFC's 50% held-out delivery criterion. The final step had no-route glacier drainage; the
measured unresolved share remains small relative to source. The 5 km depth guard reported
zero shed meltwater.


The same held-out seed 295569319 was also replayed for 10 Myr at 1 Myr steps. Both arms
reached 190.1 Myr after 10 steps. Per-step data are in
[`baseline`](g6-baseline-295569319-100-heldout.json) and
[`D8`](g6-d8-295569319-100-heldout.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 584–716 | 575–692 |
| Ice-margin networks, final step | 611 | 582 |
| Sum of step wall times | 1,842 s | 1,914 s |
| Peak RSS | 2,267 MB | 2,466 MB |

D8's 8.907e+19 m³ glacier source delivered 98.37% to lakes or ocean in aggregate
(at least 98.01% per step). The 2.113e+20 m³ bare-snow source delivered 94.74%
(at least 92.39% per step). Maximum absolute relative source-ledger residuals were
1.49e-16 and 3.06e-16. Unresolved glacier drainage was 2.21% of total source,
with at most 3.69% in a step. The maximum filled-pit count was 183 and
warm-based surface melt on ice thicker than 1 km totalled 6.431e+17 m³. This clears the
RFC's 50% held-out delivery criterion and remains well above the 20% falsification
threshold. Repeated no-route glacier drainage is included in the unresolved diagnostics;
the source ledgers close to round-off. The 5 km depth guard repeatedly logged nodes and
zero shed meltwater.

Seed 295569319 clears the held-out criterion at both step sizes. Together with seed
628977350, both preregistered held-out saves pass at 0.25 and 1 Myr steps. The separate
ice-age seed 447276134 remains report-only; it does not gate RFC 0015.


The separate held-out ice-age save, seed 447276134, was replayed for 18 Myr with an 8 Myr
ice-age period at 0.25 Myr steps. Both arms reached 134.4 Myr after 72 steps. Per-step data
are in [`baseline`](g6-baseline-447276134-025-age8m.json) and
[`D8`](g6-d8-447276134-025-age8m.json). This report-only arm does not gate RFC 0015.

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 207–943 | 203–892 |
| Ice-margin networks, final step | 527 | 487 |
| Sum of step wall times | 5,063 s | 5,091 s |
| Peak RSS | 1,096 MB | 1,081 MB |

D8's 1.006e+20 m³ glacier source delivered 98.73% to lakes or ocean in aggregate
(at least 96.42% per step). The 4.437e+20 m³ bare-snow source delivered 95.00%
(at least 92.80% per step). Maximum absolute relative source-ledger residuals were
2.59e-16 and 2.93e-16. Unresolved glacier drainage was 0.48% of total source,
with at most 2.84% in a step. The maximum filled-pit count was 115 and
warm-based surface melt on ice thicker than 1 km totalled 4.216e+17 m³. No-route fluxes
were logged throughout. Both arms repeatedly hit the 5 km depth guard; baseline shedding
rose from about 0.6 to 3.4 million km³ per step, while D8 shedding rose to about 4.7 million
km³ per step. This forcing is outside the RFC's gate scope. The machine's unit-test worker
overlapped the baseline arm; the replay arms themselves ran sequentially.

The 1 Myr step-size report is recorded below.


The same report-only ice-age save was replayed for 18 Myr at 1 Myr steps. Both arms
reached 134.4 Myr after 18 steps. Per-step data are in
[`baseline`](g6-baseline-447276134-100-age8m.json) and
[`D8`](g6-d8-447276134-100-age8m.json).

| Measure | Baseline | D8 |
| --- | ---: | ---: |
| Ice-margin networks, range over steps | 299–912 | 296–897 |
| Ice-margin networks, final step | 572 | 639 |
| Sum of step wall times | 4,293 s | 3,980 s |
| Peak RSS | 1,025 MB | 1,041 MB |

D8's 1.194e+20 m³ glacier source delivered 98.53% to lakes or ocean in aggregate
(at least 97.05% per step). The 4.227e+20 m³ bare-snow source delivered 94.90%
(at least 92.57% per step). Maximum absolute relative source-ledger residuals were
2.24e-16 and 3.14e-16. Unresolved glacier drainage was 0.71% of total source,
with at most 1.95% in one step. The maximum filled-pit count was 146 and
warm-based surface melt on ice thicker than 1 km totalled 1.108e+18 m³. The baseline's
5 km depth guard shed about 2.1–6.6 million km³ on the final three steps; D8's final step
shed about 2.2 million km³, after about 85,500 km³ on the previous step. Both arms ran
without other heavy runs. These ice-age results are report-only and do not gate RFC 0015.
