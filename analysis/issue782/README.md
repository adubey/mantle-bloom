# Issue 782: RFC 0015 gate G2 part 2 — melt-factor and snow-phase sensitivity

Gate G2 (part 2) of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) for
[#782](https://github.com/adubey/mantle-bloom/issues/782). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #8, each process follows its own physics — the melt factors and snow
phase are the least-constrained numbers in D1's mass balance, so their influence is
measured before their values are fixed.

## Preregistration

Written before any run in this report. The only data seen are G0/G1/#765 runs at the
RFC starts (`DDF_snow` 4, `DDF_ice` 8, `T_rs` 1.0 C).

### Tree, saves and arms

D1+march prototype: gate branch `claude/issue-766-g1-march` at `884fe85` (D1 physical
rates, march arm m3 with longest substep 250 yr). Its `G0Config` carries the knobs as
`ddf_snow_mm` (4.0), `ddf_ice_mm` (8.0) and `rain_snow_c` (1.0), overridable per run
with `--set KEY=VALUE`.

Run: march arm (m3) with `--set f_snow=True` (see amendment), 0.25 Myr steps, draw 0
(each world's own noise), 10 Myr:

| Save | Arms |
|---|---|
| Seed 997271774 at 87.4 Myr | base + `ddf_snow_mm` {3, 5} + `ddf_ice_mm` {5.5, 8.2} + `rain_snow_c` {0, 2} |
| Seed 649658380 at 6.4 Myr | same 7 arms |
| Seed 374464662 at 132.6 Myr, ice ages 8 Myr | none: G1's 250 yr cost run (`cost-s374`, end 629 M km³) is the report-only base; knob effects are read on the two stable saves |

The full factorial is out of budget and out of scope; one-at-a-time around the starts
is the preregistered design. Each run alone (peak ~1.2 GB).

```bash
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 \
  --checkpoint-myr 0.25 --myr 10 --set march_max_substep_years=250 --set f_snow=True \
  [--set <knob>=<value>] [--ice-age-period-years 8e6] --out <scratch-output>/<arm>.json
```

### Amendment (before any reported result): every arm sets `f_snow=True`

The launched sweep ran arm (m3) as defined, which uses the old `is_frozen` snow
switch (`f_snow=False`). Under that switch `rain_snow_c` (`T_rs`) is never read
(`ice_sheet.py:205-208`), so the trs0/trs2 arms would repeat base bit-for-bit, and
the `DDF_snow` arms would measure the old switch's snow rather than D1's `f_snow`
partition. G2 measures D1's knobs, so every arm adds `--set f_snow=True` (arm (e)'s
fix composed with the march). The `f_snow`+march combination is previously unrun;
its settling behavior is observed, not assumed.

Discarded (wrong switch; peeked but not reported as results): s997 base, snow3,
snow5, ice55 with end volumes 13.36, 13.70, 11.23, 13.37 M km³. Their logs stay in
the captures dir under the same run labels.

### Reading

Primary metric: ice volume at 10 Myr per arm. Ice area reported. Per knob, the
max/min volume across its tested values at fixed others; a knob whose range moves
volume by more than 2x is flagged weakly constrained, per the RFC. Outcome: the fixed
starting values (kept starts get reasons).

**Expected:** higher `DDF` and lower `T_rs` mean less ice; rough sizes unquantified —
this gate measures them.

### Not tested

Held-out saves (447276134, 628977350, 295569319) are never used here. The ice-age
runaway is G4's problem; this gate does not re-litigate it.

## Results

All 14 arms ran to 10 Myr alone (peak RSS ~1.2 GB) on gate branch `884fe85`.
`g2.json` is the compact result; raw logs are `g2f-<save>-<arm>.log` in
`~/mantle-bloom-captures/issue-782/` (private path).

### End volumes and areas at 10 Myr

| Arm | s997 vol (M km³) | s997 area (M km²) | s649 vol | s649 area |
|---|---|---|---|---|
| base (4 / 8 / 1.0) | 12.56 | 13.43 | 2.12 | 2.85 |
| snow3 | 13.24 | 14.11 | 2.15 | 2.98 |
| snow5 | 11.84 | 12.79 | 2.09 | 2.88 |
| ice55 | 12.68 | 14.31 | 2.12 | 3.02 |
| ice82 | 12.10 | 13.10 | 2.05 | 2.91 |
| trs0 | 12.34 | 13.48 | 2.13 | 2.96 |
| trs2 | 11.76 | 13.31 | 2.10 | 2.98 |

### Per-knob max/min volume (flag rule: >2x)

| Knob | s997 | s649 | Flag |
|---|---|---|---|
| `DDF_snow` {3, 4, 5} | 13.24 / 11.84 = 1.12 | 2.15 / 2.09 = 1.03 | none |
| `DDF_ice` {5.5, 8, 8.2} | 12.68 / 12.10 = 1.05 | 2.12 / 2.05 = 1.03 | none |
| `T_rs` {0, 1.0, 2} | 12.56 / 11.76 = 1.07 | 2.13 / 2.10 = 1.01 | none |

Max step-to-step swing after 2 Myr is 1.7–5.1% on every arm, far under the 25%
limit. Mean step times 16–102 s (the s649 trs arms march more substeps; cost is
not a gate here).

### Outcome: keep all three starts

- **`DDF_snow` 4** stays: ±1 moves s997 volume +5.4% / −5.7%, s649 ±2%.
- **`DDF_ice` 8** stays: the full sourced range {5.5, 8.2} plus the start moves
  ≤5% on both saves. Bare-ice melt rarely binds: most melt is of snow.
- **`T_rs` 1.0 C** stays: the full {0, 2} range moves ≤7%, and the knob is live
  (trs0 and trs2 differ on both saves, so the dead-switch failure mode is
  excluded). Note the ordering on s997 (trs2's 11.76 below trs0's 12.34) runs
  against the naive more-snow→more-ice reading; the gap is within the save's
  bistable scatter (E2 split 12.2–13.1 under identical forcing), so no
  mechanism is claimed.

No knob is flagged weakly constrained. The ice-age arm is G1's 250 yr cost run
(end 629 M km³, runaway) as preregistered; it does not enter the flag rule.
