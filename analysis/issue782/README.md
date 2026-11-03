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

Run: march arm (m3), 0.25 Myr steps, draw 0 (each world's own noise), 10 Myr:

| Save | Arms |
|---|---|
| Seed 997271774 at 87.4 Myr | base + `ddf_snow_mm` {3, 5} + `ddf_ice_mm` {5.5, 8.2} + `rain_snow_c` {0, 2} |
| Seed 649658380 at 6.4 Myr | same 7 arms |
| Seed 374464662 at 132.6 Myr, ice ages 8 Myr | none: G1's 250 yr cost run (`cost-s374`, end 629 M km³) is the report-only base; knob effects are read on the two stable saves |

The full factorial is out of budget and out of scope; one-at-a-time around the starts
is the preregistered design. Each run alone (peak ~1.2 GB).

```bash
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 \
  --checkpoint-myr 0.25 --myr 10 --set march_max_substep_years=250 [--set <knob>=<value>] \
  [--ice-age-period-years 8e6] --out <scratch-output>/<arm>.json
```

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

Pending.
