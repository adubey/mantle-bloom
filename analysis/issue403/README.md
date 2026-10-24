# RFC 0001–0007 baseline after #358

Issue [#403](https://github.com/adubey/mantle-bloom/issues/403). Model baseline:
`2562f6381b5777b021f68a958a196bf99d0d25aa` (origin/main when this investigation
started). `backend/app/` and the pinned dependency constraints are unchanged.
The [provenance manifest](provenance.json) separates known model commits from
publication commits and the old save's unknown generating commit.

## Reproduce

Run from the repository root after `./bin/setup.sh`. Run one world at a time:

```sh
backend/.venv/bin/python bin/debug/rebaseline_rfcs.py \
  --out analysis/issue403/tuning.json --save-dir /tmp/mantle-bloom-issue403
backend/.venv/bin/python bin/debug/rebaseline_rfcs.py --seeds 456862899 --years 3e9 \
  --out analysis/issue403/long.json --save-dir /tmp/mantle-bloom-issue403
backend/.venv/bin/python bin/debug/ledger_history.py \
  /tmp/mantle-bloom-issue403/seed456862899.mbworld \
  --json analysis/issue403/long-history.json
```

`rebaseline_rfcs.py` calls `compare_cratons.run_job` with the same inputs as
`compare_cratons.py --seeds 7,11,2024 --variants on --years 1e9 --every 100e6
--step 5e6 --workers 1`. It adds read-only taps for individual collision merges,
`cratons.attrite`, and continental inventories, saving the final world without a
second simulation. The regenerated seed 456862899 uses the same defaults to
3 Gyr. No simulation controls are overridden. JSON records controls, Python
version, constraints checksum, baseline commit, and save checksums. Saves and
full logs stay outside Git; the JSON is sufficient to review every quoted rate.

All runs use 5 Myr steps, default density and the 2 Myr year-based split
cooldown. This investigation does not test sensitivity to other step sizes or
to the cadences still counted in steps.

## Definitions and limitations

- Craton sink shares use the **live craton volume at the window start** and
  physical window duration, not a fixed initial volume across multiple Gyr.
  Rifting path accounts are subdivisions of `rifted_m3`, not additional sinks.
  `unattributed_m3` is a signed audit residual, not an identified physical loss.
  Gross destruction and net live-volume change differ when new craton forms.
- `collision_merges` measures each call to
  `_merge_plates_with_craton_audit`; its ledger delta attributes loss directly
  to that merge. Tiny floating-point remap residuals are not geological loss.
  It does not capture the placement inputs needed for RFC 0001's H1–H4 gate.
- `attrited_craton_km3` is the direct return from `cratons.attrite`.
  The remaining `subducted_m3` is **not** automatically trench retreat: collision
  partition and suture-placement remainder also write that account. Continental
  material attrition can differ from craton attrition; do not subtract this
  craton measure from the continental ledger.
- Continental source/sink rates are terminal ledger fluxes. Juvenile additions
  are not a measure of gross arc production or net preservation. Numerical and
  declined marine-deposition sinks are included in total terminal flux.
  Balance errors are reported separately and are not added to physical sinks.
- `land_fraction_node` is the raw node-count share in `compare_cratons`;
  saved-history `land_fraction` is the area-weighted, connectivity-aware stats
  measure. They are separate metrics. Continental mean inventories in
  `ledger_history` are sample means over the uniformly spaced 5 Myr snapshots.
- `ledger_history.py` now includes each Gyr endpoint (the previous `-1` dropped
  the last sample's loss), and normalizes partial windows by their measured
  duration. Old saved-summary percentages remain historical as published;
  current numbers use the corrected summarizer.
- The old private `mantle-bloom-seed456862899-19191800000y.mbworld` is **not
  available for download**. Its generating commit was not recorded; `d795a09`
  is the first RFC publication commit, not its simulation baseline. Per #403,
  it used 0.1 Myr steps through 21.8 Myr then 10 Myr steps: a pre-#358 20-step
  cooldown of 2 Myr initially and **200 Myr for most of the run**. Its loss
  figures document discovery of the problem; they do not calibrate current
  main. The fresh 3 Gyr run is a new experiment, not a replay of its topology.

## Findings

The tuning runs reproduce the earlier post-#358 craton baseline. Delamination
is the largest pooled craton sink (erosion is larger on seed 11 individually),
and direct merge taps account for the delamination ledger. Of 55, 27 and 36
merges, 27, 7 and 9 lose more than 1 km³ of craton; this reporting floor excludes
floating-point remap noise. The largest individual losses are 387.96, 280.87
and 157.48 M km³, or 16.66%, 12.54% and 14.08% of global live craton immediately
before the merge. These identify substitute capture candidates for RFC 0001,
not the H1–H4 placement diagnosis.

Direct attrition accounts for the tuning runs' `subducted_m3`; the contact-owner
question in RFC 0002 remains open. Rifting is 8.97–15.62% of window-start live
craton per Gyr, and the existing phase 1b selection of steering together with
stronger resistance stands. No new F measurement is inferred from cut counts.

Continental inventory falls 73.0–86.0% in the first Gyr, with terminal fluxes
9.88–18.32 km³/yr. Seed 7's collision-placement remainder rivals deep subduction
in the early run; seed 11 remains deep-subduction dominated. A sediment-first
diagnosis cannot be selected from the pooled ledger. Continental balance errors
are −0.1013%, +0.0866% and +0.0509% of initial inventory; they are reported,
not silently treated as closure or as geological losses.

Formation's last positive samples are 355, 470 and 335 Myr. Craton share starts
at 39.9%, 50.7% and 37.6%, then reaches 85.2–91.7% at 500 Myr while genuine
continental area shrinks. RFC 0004's 30–65% band fails without welded growth;
the ratio alone cannot diagnose a formation ratchet.

| RFC | Refreshed conclusion / gate | Follow-up |
|---|---|---|
| 0001 | Merge loss remains large; placement diagnosis still required; recommendation unchanged | Existing H1–H4 gate |
| 0002 | Attrition is measurable; proposed margin ownership must be resolved | [#405](https://github.com/adubey/mantle-bloom/issues/405) |
| 0003 | Preserve published post-#358 F result and selected 1 + 3 branch | Existing [#389](https://github.com/adubey/mantle-bloom/issues/389), #388 / PR #392 |
| 0004 | Revisit source priority and share-based acceptance / falsification | [#409](https://github.com/adubey/mantle-bloom/issues/409) |
| 0005 | Rank collision and trench terminal sinks before selecting a recycling fix | [#408](https://github.com/adubey/mantle-bloom/issues/408) |
| 0006 | Existing elevation-driven reset can create continent; recovery gate must isolate arcs | [#406](https://github.com/adubey/mantle-bloom/issues/406) |
| 0007 | Isolate plume / docking / reset sources and compare breakup against background splits | [#407](https://github.com/adubey/mantle-bloom/issues/407) |

The tables below include every craton sink account and terminal continental
account; zeros are retained. Full-precision numerators, window-start stocks,
cut counters, checkpoint plate counts, and per-Gyr topology are in the JSON
files.

### Regenerated seed 456862899: 3 Gyr stress run

The regenerated run is a fresh current-main experiment. Its first-Gyr craton
sink mix is comparable in scale to the tuning runs, but the dominant sink
changes over time: there are no merge-attributed losses after the first Gyr,
while the direct attrition tap rises to 32.82% of each Gyr's starting live craton
in the final interval. The 3 Gyr result therefore supports ordering the merge,
attrition and erosion diagnoses by time window rather than treating the first
Gyr as a steady-state rate.

| Gyr | Start live craton (M km³) | Delaminated | Subducted | Rifted | Eroded | Collision reworked | Signed unattributed |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0–1 | 1,448.98 | 27.65% | 15.45% | 9.06% | 21.34% | 1.81% | +0.20% |
| 1–2 | 355.42 | 0.00% | 21.69% | 7.14% | 3.81% | 0.93% | −0.70% |
| 2–3 | 238.60 | 0.00% | 32.82% | 11.81% | 9.94% | 2.07% | −0.03% |

The first Gyr has 12 merges, including 6 with more than 1 km³ of craton loss;
the largest removes 181.00 M km³ (13.77% of live craton immediately before
that merge, at 125 Myr). Across the full run there are 15 merges, 261 splits
and 231 failed rifts; plate count changes from 9 to 58. The merge-loss tap
matches the first-Gyr `delaminated_m3` ledger delta. Placement inputs were not
captured, so this still does not pass RFC 0001's H1–H4 gate.

| Gyr | Deep subduction | Collision subduction | Collision lower crust | Rift thinning | Other terminal sinks | Total terminal flux (km³/yr) |
|---|---:|---:|---:|---:|---:|---:|
| 0–1 | 5.500 | 0.810 | 0.018 | 0.444 | 0.012 | 6.784 |
| 1–2 | 1.733 | 0.000 | 0.001 | 0.136 | 0.000 | 1.870 |
| 2–3 | 1.427 | 0.001 | 0.002 | 0.123 | 0.002 | 1.554 |

Continental inventory falls from 3,864.20 to 190.38 M km³ (−95.07%), while
genuine continental area falls from 82.33 to 2.72 M km². Node-count land
fraction falls from 14.8% to 0.5%. The residual is −9.119 M km³ (−0.2360% of
the initial inventory); it remains a reported balance error, not a physical
sink. Craton formation totals 0.722 M km³ and last has a positive sample at
255 Myr. Craton area reaches 105.0% of genuine continental area at 3 Gyr as
the continental-area denominator collapses, so RFC 0004's area-share gate
cannot be read as a formation-rate signal on this run.

The long run does not measure the RFC 0006 source-attribution gate, RFC 0007's
matched plume-off breakup baseline, or RFC 0001's placement inputs. Those
remain explicit diagnosis work tracked by [#406](https://github.com/adubey/mantle-bloom/issues/406),
[#407](https://github.com/adubey/mantle-bloom/issues/407), and the H1–H4 gate.
