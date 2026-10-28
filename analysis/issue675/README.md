# RFC 0014 follow-up: large-lake flicker mechanisms

## Scope and provenance

This investigation extends the final RFC 0014 acceptance analysis in [#673](https://github.com/adubey/mantle-bloom/issues/673) and tracks the large-lake follow-up in [#675](https://github.com/adubey/mantle-bloom/issues/675). It does not change model behavior or the accepted flicker threshold.

The registered seed-997271774 save (SHA-256 `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7`) was replayed for 0.8 Myr with ice-age control `never`. The source files for lake solving, hydrology, erosion and world stepping match the #673 model-source checksums. The replay worktree was at `f22bda1785a8792453e4513e95d2655d41267eb5`; its measured model base is `6ebf8cd19ee0098822ff9eae6685a22e20f9db32`. The local save is not available for download. The compact results and eight geographic anchor vectors are in this directory.

Three dynamic schedules were compared: 100 kyr (8 steps), 25 kyr (32 steps), and 10 kyr (80 steps). At each step, `trace_lake_anchor.py` finds the nearest grid cell to each anchor and records all connected wet components containing those cells. Anchors are locations, not persistent lake IDs. The maximum nearest-cell distance was 2.9 km in the longest replay.

## Results

### Timestep changes the connected basin at the target location

The eight anchors identify the shallow 59.62 km³ component at 131.2 Myr in the registered 100 kyr replay. The same locations give very different component sizes under the three dynamic step schedules:

| Elapsed age | 100 kyr | 25 kyr | 10 kyr |
|---:|---:|---:|---:|
| 130.55 Myr | next checkpoint at 130.6: 20,090 km³ / 931,713 km² | 19,870 km³ / 931,713 km² | 67,807 km³ / 1,310,953 km² |
| 131.2 Myr | 59.62 km³ / 31,119 km² | 77,620 km³ / 1,446,251 km² | 84,985 km³ / 1,489,253 km² |
| 131.3 Myr | 157.59 km³ / 65,796 km² | 79,663 km³ / 1,462,044 km² | 86,642 km³ / 1,512,189 km² |

All size entries are volume / flooded area. In the 10 kyr run, the anchor component grows from 19,903 km³ / 931,713 km² at 130.54 Myr to 67,807 km³ / 1,310,953 km² at 130.55 Myr. In the 25 kyr run, connectivity oscillates during the transition: 36,262 km³ / 1,122,339 km² at 130.6 Myr, 67,465 km³ / 1,314,945 km² at 130.65, 45,921 km³ / 1,233,657 km² at 130.7, and 67,088 km³ / 1,314,945 km² at 130.725.

The 100 kyr trajectory also alternates wet and dry at the anchor location: it maps to no wet component at 130.7 Myr, to 49.1 km³ / 23,381 km² at 130.8, 39.7 km³ / 23,381 km² at 130.9, no wet component at 131.0, 147.8 km³ / 61,995 km² at 131.1, 59.62 km³ / 31,119 km² at 131.2, and 157.59 km³ / 65,796 km² at 131.3.

The standard #588 large-body reversal detector counted zero qualifying reversals in the 10 kyr run. That detector asks for a sign reversal in consecutive volume changes above a threshold; it does not flag the sustained connectivity change seen by geographic matching. This is not a pass of the RFC 0014 flicker criterion.

At 131.2 and 131.3 Myr, the 10 and 25 kyr schedules agree within about 10% in component volume and 3% in area, while the 100 kyr trajectory remains a much smaller, intermittently connected basin. The transition timing is not monotonic in step size: the 10 kyr schedule merges earlier than 25 kyr, and the 100 kyr trajectory follows another path. Because every world process updates at each chosen step, this is evidence of dynamic trajectory and topology sensitivity, not a lake-solver-only or single-mechanism test.

All global lake ledgers remain closed. Maximum normalized residuals were `1.93e-13` for 10 kyr and `5.58e-13` for 25 kyr, below the 0.1% requirement.

### Evaporation scale

The current temperature-only evaporation proxy averaged 1.445 m/year on the shallow seed-997271774 footprint immediately before the 100 kyr rebound (10th–90th percentile 1.440–1.451 m/year). Its 59.62 km³ storage over 31,119 km² corresponds to about 1.33 years of no-inflow storage at that rate. For the deep seed-373936961 component immediately before its step-5 drop, the proxy averaged 0.817 m/year; its 136,570 km³ over 185,010 km² corresponds to about 904 years of no-inflow storage.

These are diagnostic turnover ratios, not drought forecasts: catchment supply continues in the replay. The saves do not carry wind, humidity, radiation, or deep-water temperature, so this does not test wind-enhanced evaporation or vertical thermal moderation. The rates do show why shallow storage is sensitive to the interaction of supply, evaporation and event timing.

### Sediment scale and placement

At the current modeled deposition rate, the prior step's silt on the shallow lake's eight-cell footprint is about 56 km³, roughly 57% of its 97.97 km³ rebound. Deposition is applied after that step's water solve; the next solve still expands the connected area and raises storage. On the deep event's 44-cell footprint, prior-step deposition is about 416 km³, roughly 0.3% of the 132,450 km³ component drop. These are equal-area footprint-volume screens, not sediment-only counterfactuals.

For the shallow case, the solver's one-hop spill target before expansion received about 0.04 m of silt in the prior step. The reported target after expansion was a different dry cell. The one-hop target is not a surveyed geomorphic sill, so these figures do not establish connector aggradation or sediment causality.

### Subsurface storage screens

Earlier capacity-only estimates remain upper-bound scale checks, not observed or dynamically exchanged water:

- Shallow-basin pore capacity is about 18.6 km³ for a 1 m sediment layer at 0.3 porosity, or 248 km³ for 10 m at 0.4 porosity. The smaller case is about 19% of the rebound; the larger exceeds it but is not all empty, connected or exchangeable.
- A 60 km nearshore band, 3 m head rise and specific yield 0.2 gives about 32 km³ of shallow-basin bank storage, about one-third of the rebound. Actual capacity and return times depend on lakebed and aquifer geology.

No porewater or shoreline-aquifer storage was added to the model in this investigation. The scratch per-body transfer tap was discarded because hierarchy flush passes recount gross transfers; the global ledger cannot identify event-basin ancestry.

## Interpretation and next test

The evidence supports timestep-sensitive topology at the shallow-lake location. It does not show that evaporation, groundwater exchange, thermal stratification, or sediment is irrelevant. It does show that mechanism tests must hold the integration schedule and non-lake forcing fixed, and must track basin lineage geographically with basin-level inflow and overflow budgets. The RFC's existing flicker threshold remains open and unchanged. The next causal sequence is:

1. Freeze the timestep and non-lake terrain/climate inputs around a geographically tracked event window; add basin-lineage accounting that closes with the whole-world ledger.
2. Test a bounded shallow-basin porewater and shoreline-storage exchange range one at a time, with equal-and-opposite reservoir transfers.
3. Test connector-local sediment using the actual sill/floor clearance and deposition field, not a lake-wide median.
4. Add wind and vertical heat-storage forcing only with explicit meteorological and mixed-layer inputs; keep Lake Winnipeg as a cold, smaller geometry/wind analog rather than coefficient calibration.

These are exploratory gates, not RFC acceptance-criterion changes. Only the owner may revise the existing flicker threshold.

## Reproduction

Use the registered seed-997271774 `.mbworld` save, which is local and not available for download. The geographic trace script accepts the save and the checked-in anchor coordinates:

```bash
backend/.venv/bin/python bin/debug/trace_lake_anchor.py <save> \
  --anchors analysis/issue675/shallow-lake-anchor.json \
  --backend <backend-at-6ebf8cd-or-source-identical> \
  --myr 0.8 --step-years 10000 --checkpoint-myr 0.01 \
  --out /tmp/seed997271774-anchor-10k.json
```

Repeat with `--step-years 25000 --checkpoint-myr 0.025` and `--step-years 100000 --checkpoint-myr 0.1`. The script sets the ice-age control to `never`, writes each step to the output JSON, and reports geographic anchor-to-component matches and normalized ledger residuals. The model-source SHA-256s used for the registered runs are recorded in `timestep-and-mechanism-screen.json`.
