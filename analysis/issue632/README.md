# RFC 0014 D2(c): production G5 replay

## Question and governing principles

Does the owner-selected 2 m physical ice ceiling behave as a bounded, conservative lake-ice
reservoir when production climate, terrain, routing, and liquid-plus-ice ledgers advance
through long steps?

The prior behavior violated **design principle 3** (represent seasonal processes in aggregate
across long steps) by applying sublimation in proportion to open water, inverting the intended
ice-covered loss. It also violated **principle 6** (account for material through every process)
when the phase-specific lake ledger did not expose that loss correctly. The fix weights
sublimation by ice-covered fraction and books it as loss from the ice reservoir.

The observable invariants are: physical ice thickness ≤2 m; ice water equivalent does not
exceed local liquid-plus-ice inventory; liquid water remains beneath ice when present; deep
bedfast fraction (water depth >2 m) ≤1%; and each production step's combined lake-water ledger
residual ≤0.1% of positive external supply. The empirical cap represents the owner's
non-permafrost scope. This model does not calculate ground heat or permafrost.

## Production diagnosis

Three preregistered diagnosis saves completed 10 Myr replays at 100 kyr steps (101 snapshots
per save), with ice_age_period_years = 0. See [production-g5-summary.json](production-g5-summary.json)
for checkpoint metrics, inputs, and provenance. Replays used backend commit
`855d5fe574f18c9efa0779ed94677eb6d65dcbc5`.

| Seed | Runtime | Ice water-equivalent range | Liquid-under-ice area max | Shallow bedfast max | Deep bedfast max | Max ledger residual | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 910211954 | 1,796 s | 0–2,721 km³ | 62.09% | 4.64% | 0% | 2×10⁻⁹% | 5,693 MB |
| 635056485 | 906 s | 0–26.5 km³ | 17.49% | 0% | 0% | 2×10⁻⁹% | 3,644 MB |
| 317739385 | 1,856 s | 0–1,771 km³ | 49.00% | 3.63% | 0% | 1×10⁻⁹% | 5,919 MB |

All snapshots met the physical cap, nonnegative-reservoir, local-inventory, deep-bedfast, and
ledger checks. Ice volume changed over the run (0.1-to-10 Myr changes: −14.36%, +36.17%, and
−58.99%, respectively); this is descriptive behavior under evolving forcing, not a plateau
gate. The plateau is only expected in the fixed-forcing prototype. Shallow bedfast area is
reported separately. Cumulative sublimation is recorded in the JSON. A few computed physical
thickness values exceed 2 m by floating-point roundoff at about 10⁻¹⁵ m; they are 2 m within
numerical precision.

## Held-out acceptance

Three reserved #538 saves completed ten 100 kyr production steps each. Results and checksums
are in [heldout-g5-summary.json](heldout-g5-summary.json); full replay outputs remain in
scratch storage.

| Seed | Start age | Wet lake nodes | Max physical ice | Deep bedfast max | Liquid-under-ice area max | Max ledger residual | Zero-supply steps | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 373936961 | 244.4 Myr | 275–568 | 2.00 m | 0% | 49.88% | 1.02×10⁻¹⁰% | 0 | 4,714 MB |
| 997271774 | 130.5 Myr | 218–549 | 2.00 m | 0% | 25.63% | 1.10×10⁻¹⁰% | 0 | 2,035 MB |
| 456862899 | 19,191.8 Myr | 0–1 | 0 m | 0% | 0% | 1.82×10⁻⁹% | 7 | 2,874 MB |

All five coded acceptance checks pass on all three saves: physical ice ≤2 m, nonnegative
reservoirs, ice within local water inventory, deep bedfast ≤1%, and per-step ledger residual
≤0.1% of positive external supply. The first two saves had lake ice during the replay and
retained liquid under ice. The extreme-age third save had no lakes for seven steps and one
liquid-only node for three; its lake-ice and under-ice checks are non-informative. Therefore
held-out evidence is positive but comes from two informative ice-bearing saves.

The first held-out load exposed an incompatible collision-polarity cache in older saves;
persistence now drops cached frames without stable cell keys so production rebuilds them. A
focused old-pickle regression check passes. Glacier depth-limit warnings during replay are
outside the lake-ice metrics.

## Verification and reproduction

Final checks passed: ruff on the changed hydrology files, `./bin/quality_check.sh`, and
`./bin/affected_test.sh --base origin/main` (1,876 passed). The held-out and production replays
were each run sequentially, with no replay overlapping tests or another replay.

Run one save at a time with the matching local save substituted for `<local-save>`:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_production.py \
  <local-save> --myr 10 --step-years 100000 --checkpoint-myr 0.1 \
  --out analysis/issue632/production-g5-seed<seed>.json
```

For held-out acceptance, use a preregistered #538 save and
`bin/debug/gate_lake_ice_heldout.py` with `--myr 1 --step-years 100000 --checkpoint-myr 0.1`.
