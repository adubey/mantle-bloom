# RFC 0015 D8 and G6

## Replay plan

Compare D8 with the pre-D8 main commit `eb3cbc4febbfa7a2df51c30d01916ab091a9f647`.
For each arm record the exact implementation commit and save SHA-256. Use
`bin/debug/replay_d8_meltwater.py` with the same save, controls and step schedule in both
arms. The diagnosis saves are seeds 997271774 and 649658380 with ice ages off. Replay each
for 10 Myr at 0.25 Myr and 1 Myr steps. Report seed 374464662 with ice ages at an 8 Myr
period for 18 Myr at both step sizes, and with a 1 Myr period at 0.25 Myr steps. The ice-age
arms do not gate the RFC. The save inventory and original hashes are in
[`analysis/issue548`](../issue548/README.md). Converted v5 saves have different checksums;
record their checksums as well as their original hashes. The script enables the glacier
march in both arms and records the saved march controls and ice-age period.

Each step records the two source budgets, their residuals, filled hydraulic depressions,
their throughput, unresolved drainage, surface melt routed from thick warm-based ice, and
the count and ten largest discharges of ice-margin river networks. A filled pit is a
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

G6 and held-out runs have not yet been made. A one-step development smoke replay on seed
649658380 at 0.25 Myr used branch commit `809d123351e727f5789648c368c0a05817f4894e` and
converted save SHA-256 `06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f`.
The glacier and bare-snow source residuals were 16 m³ of 1.88 × 10¹⁷ m³ and −650 m³ of
4.13 × 10¹⁸ m³, respectively. No drainage was unresolved. This development run used the
earlier connected-region count for pits and is not a G6 measurement.
