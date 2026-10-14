# Issue #378: craton rifting by path (RFC 0003 phase 1)

The gate in [RFC 0003](../../rfc/0003-rifting-around-cratons.md): before steering split cuts
or giving failed rifts craton resistance, find out which rifting path destroys cratons.

## Method

`bin/debug/compare_cratons.py --seeds <S> --variants on --years 1e9 --every 100e6`, seeds 7,
11 and 2024, default density, run one at a time on this branch. The branch only adds
diagnostics: nothing it books feeds back into the model, so these runs are origin/main's
model. Numbers are in [rift_paths.json](rift_paths.json).

- `rifted_divergent_m3`: thinning in `lithosphere_plate.deform_columns`.
- `rifted_stretch_m3`: donors to new rift cells, `quad_tectonics._open_rift`.
- `rifted_failed_m3`: failed rifts, `merge_split.maybe_split_plate`.
- A cut "goes through craton" when any node within the failed-rift band
  (`FAILED_RIFT_BAND_MULT` = 2.5 spacings) of the cut holds craton.

## Results (0–1 Gyr)

| Seed | Rifted, % of start craton | Divergent | Stretch | Failed | Splits through craton | Failed rifts through craton | Craton share of band area (split / failed) |
|---|---|---|---|---|---|---|---|
| 7 | 4.0% | 32% | 51% | 17% | 5 of 33 (15%) | 4 of 29 (14%) | 11% / 18% |
| 11 | 4.3% | 28% | 58% | 14% | 6 of 37 (16%) | 5 of 31 (16%) | 7% / 6% |
| 2024 | 10.7% | 41% | 50% | 10% | 4 of 32 (12%) | 5 of 28 (18%) | 4% / 9% |

Path columns are shares of `rifted_m3`. The craton ledgers balance to within 3e-6 km³.

## Findings

- **Stretch donation is the largest path** (50–58%), then divergent thinning (28–41%).
  Failed rifts take 10–17%.
- **Solution 2 alone can't meet the target.** Failed rifts cost 0.4–1.0% of the live craton
  per Gyr. Removing that path entirely would leave 3.4–9.7%/Gyr, against the RFC's 2%/Gyr
  target.
- **Few cuts cross craton.** Only 12–18% of split and failed-rift cuts have craton in their
  band. Yet divergent and stretch losses, which run along every divergent boundary, are 83–90%
  of the sink. That points toward the RFC's first falsification case: much of the loss
  happens at divergent boundaries that later run through cratons, not at the moment a cut is
  made. These counters can't prove it, because they don't link a divergent boundary to the
  cut that created it.
- Seed 2024 rifts 2.5× more than the others, with a similar share of cuts through craton.
  Its divergent path is 4.4%/Gyr, against 1.2–1.3%/Gyr on seeds 7 and 11.

## What this means for the next phase

Solution 2 (failed-rift resistance) is a cheap consistency fix and worth doing, but it acts
on 10–17% of the sink. Before building solution 1 (steering the cut), measure how much
divergent and stretch loss happens at boundaries born from a craton-crossing cut. One way:
tag nodes inside each cut's band with the cut's time, and book later divergent and stretch
loss by whether the node was tagged. If that share is small, steering the cut won't reach the
target, and the work should go toward boundary migration or solution 3, as the RFC's
falsification section says.

Not run: seed 456862899 to 3 Gyr (the RFC's acceptance run for solutions 1 and 2, not the
gate).
