# Issue #588: C-default in-step lake fill-spill-merge

## Preregistration

Registered before the diagnostic replay. Compare the #588 branch with a task-specific
detached `origin/main` worktree over ten 100 kyr steps on each of the three diagnosis saves
registered on #538. Use the same save bytes and the ice-age control `never`. The held-out
acceptance saves remain reserved.

| Seed | Age | Density | SHA-256 |
|---:|---:|---:|---|
| 910211954 | 220.2 Myr | 4 | `bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2` |
| 635056485 | 452.8 Myr | 4 | `ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045` |
| 317739385 | 200.8 Myr | 6 | `6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd` |

**Expected direction:** the C path should remove split-reset storage jumps and close each
step's lake budget. It will retain finite-step storage and may change lake volume and area
relative to the balanced-level baseline; no directional prediction is registered for those
trajectory measures. The default 64 integration slices may increase runtime; the issue does
not pre-register an allowable runtime change.

**Primary invariant:** on the branch, absolute liquid-ledger residual divided by external
lake supply is at floating-point tolerance per step. Report the maximum normalized residual
and any independent unledgered loss estimate. The baseline has no lake ledger, so its
unledgered loss estimate is measured from opening storage + external supply − evaporation −
exported overflow − closing storage; positive unexplained deficit is reported as loss, and a
negative value as creation.

**Secondary measures:** global visible lake volume and area per step, large-body flicker
reversals using the RFC's definition (consecutive volume-change signs reverse, with both
changes exceeding 50% of the prior volume, for bodies above 10,000 km²), breached-pit counts,
floor-drain toggles, mean step runtime and peak RSS. These diagnose changes and are not
standalone pass/fail thresholds for this issue.

Runs will be serial, with no tests or other replay/probe running. Record the backend commits,
save checksums, commands, runtime and RSS. Full outputs stay in scratch; check in only compact
summary metrics.
