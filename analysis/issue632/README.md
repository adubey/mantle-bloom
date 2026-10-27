# RFC 0014 D2(c): production G5 replay

## Pre-registration

Before running the production replay, the expected direction and metrics are registered here.
The replay checks whether the owner-selected 2 m non-permafrost ice cap behaves correctly when
the full world, liquid routing, and lake ledger advance together.

### Inputs and controls

- Diagnosis saves only: seeds 910211954, 635056485, and 317739385, with the SHA-256 values
  registered in #538 and repeated in `analysis/issue541/README.md`.
- Held-out acceptance saves registered in #538 are used only for the final acceptance replay;
  no thresholds or implementation choices are tuned against them.
- Replay 10 Myr per save in 100 kyr steps, with the ice-age control set to `never`
  (`ice_age_period_years = 0`). Run one save at a time.
- Use the production `world.step_world` path, including evolving climate, terrain, liquid
  routing, and the lake liquid-plus-ice ledger. The 2 m physical cap is not retuned.

### Metrics and expected direction

- Maximum physical lake-ice thickness must remain at or below 2 m.
- Water-equivalent lake ice must not exceed the combined local lake-water inventory.
- Deep-water bedfast area (total water depth >2 m) must remain at or below 1% of deep-water
  lake area. Report shallow bedfast area separately.
- Report cap volume at 0.1 Myr intervals, the change from 0.1 to 10 Myr, and integrated
  sublimation. The static-prototype plateau is a reference; evolving production forcing means
  volume trajectory is descriptive rather than a separate pass/fail threshold.
- Every step's combined liquid-plus-ice ledger residual must be at most 0.1% of that step's
  external lake supply. Report steps with zero supply separately.
- Report liquid volume under ice, so the run shows whether storage and routing continue while
  the cap is present.

The owner-selected cap and ledger criteria are the primary checks. A failed bound or ledger
criterion requires a code correction before held-out acceptance. No diagnosis result will be
used to alter the 2 m cap.

## Results

All three preregistered diagnosis saves completed 10 Myr production replays at 100 kyr
steps (101 snapshots per save), with ice_age_period_years = 0. Results are in
[production-g5-summary.json](production-g5-summary.json); it records each checkpoint's
lake area, liquid and ice volumes, liquid-under-ice area, physical ice maximum, inventory
excess, shallow/deep bedfast fractions, cumulative sublimation, ledger residual, and peak
RSS. Inputs and checksums matched the pre-registration.

| Seed | Runtime | Ice water-equivalent range | Liquid-under-ice range | Shallow bedfast range | Deep bedfast max | Max ledger residual | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 910211954 | 1,857 s | 0–5,992 km³ | 0–71.0% | 0–70.6% | 0% | 3×10⁻⁹% | 4,330 MB |
| 635056485 | 1,016 s | 0–258 km³ | 0–19.0% | 0–92.6% | 0% | 2×10⁻⁹% | 3,651 MB |
| 317739385 | 1,940 s | 0–1,755 km³ | 0–55.5% | 0–65.0% | 0% | <0.5×10⁻⁹% (rounds to 0) | 5,753 MB |

Every checkpoint passed the preregistered bounds: maximum physical ice thickness was 2 m,
water-equivalent ice never exceeded local combined lake-water inventory, deep bedfast area
was zero, and stepwise lake-ledger residuals were far below the 0.1% tolerance. Liquid
storage remained under ice during every frozen checkpoint. Shallow bedfast fractions varied
widely and are reported separately; they do not count toward the deep-water bound.

Ice volume was strongly time-dependent, not a plateau. The 0.1-to-10 Myr changes were
+49.5%, +163.5%, and −30.0% for seeds 910211954, 635056485, and 317739385, respectively.
This is descriptive production behavior under evolving climate and terrain, not a failed
gate. Cumulative sublimation was recorded in the JSON for each run. Glacier depth-limit
warnings appeared during the broader world simulation; they concern glacier handling and
did not correspond to any lake-ice gate violation.

An earlier attempt at seed 910211954 was discarded after the ledger exposed an implementation
bug. The production limiter was corrected and committed as ae745ec; the reported result is
the complete rerun on that corrected commit.

## Held-out acceptance

The three reserved #538 saves completed ten 100 kyr production steps each with the ice-age
control set to `never`. Results are summarized in
[heldout-g5-summary.json](heldout-g5-summary.json); only compact metrics are checked in, while
full replay outputs remain in scratch storage.

| Seed | Start age | Wet lake nodes across 10 steps | Max ice | Deep bedfast max | Liquid under ice max | Max ledger residual | Zero-supply steps | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 373936961 | 244.4 Myr | 369–982 | 2.00 m | 0% | 50.94% | 3.39×10⁻¹⁰% | 0 | 4,052 MB |
| 997271774 | 130.5 Myr | 500–823 | 2.00 m | 0% | 41.80% | 6.98×10⁻¹¹% | 0 | 2,074 MB |
| 456862899 | 19,191.8 Myr | 0–1 | 0 m | 0% | 0% | 1.86×10⁻⁹% | 7 | 2,518 MB |

All five coded acceptance checks passed on all three saves: the 2 m physical ice bound,
nonnegative reservoirs, ice within local water inventory, deep bedfast area at most 1%, and
per-step ledger residual at most 0.1% of positive external supply. The first two saves had
lake ice throughout. The extreme-age third save had no lake water in seven steps, then one
liquid-only lake node for three steps; its lake-ice and under-ice checks are therefore
non-informative, though its water ledgers close when supply occurs. Keep this limitation in
view when interpreting aggregate held-out acceptance.

The first held-out load exposed a stale pre-#633 saved collision-polarity cache without stable
per-node keys. Persistence now drops such incompatible cached frames so the next production
step rebuilds them. A focused old-pickle regression check passes. The first save also emitted
the existing glacier 5 km depth-limit warning during replay; this is separate from the lake
ice metrics.

## Reproduction

Run one save at a time with the matching local save substituted for `<local-save>`:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_production.py \
  <local-save> --myr 10 --step-years 100000 --checkpoint-myr 0.1 \
  --out analysis/issue632/production-g5-seed<seed>.json
```

For held-out acceptance, substitute a preregistered #538 save and use
`bin/debug/gate_lake_ice_heldout.py` with `--myr 1 --step-years 100000 --checkpoint-myr 0.1`.
