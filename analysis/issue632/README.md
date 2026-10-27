# RFC 0014 D2(c): production G5 replay

## Pre-registration

Before running the production replay, the expected direction and metrics are registered here.
The replay checks whether the owner-selected 2 m non-permafrost ice cap behaves correctly when
the full world, liquid routing, and lake ledger advance together.

### Inputs and controls

- Diagnosis saves only: seeds 910211954, 635056485, and 317739385, with the SHA-256 values
  registered in #538 and repeated in `analysis/issue541/README.md`.
- Held-out acceptance saves remain untouched and will not be used to tune this implementation.
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

Pending the production replay.

## Reproduction

Run one save at a time with the matching local save substituted for `<local-save>`:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_production.py \
  <local-save> --myr 10 --step-years 100000 --checkpoint-myr 0.1 \
  --out analysis/issue632/production-g5-seed<seed>.json
```
