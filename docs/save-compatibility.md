# Save compatibility and legacy line-backed worlds

Issue #248 (Phase 5b of #228). This document is the compatibility policy for `.mbworld` saves
as the plate surface moves from `PlateWithLines` (latitude rows of nodes) to
`PlateWithSparseQuadPatch` (cube-sphere cells). It covers which saves load, what happens to a
line-backed save, and what a conversion is allowed to change.

Code: [`persistence.py`](../backend/app/persistence.py) (envelope, versions, errors),
[`legacy_conversion.py`](../backend/app/legacy_conversion.py) (the converter; its module
docstring describes the method in detail). Tests: `backend/unit_tests/test_persistence.py`,
`backend/unit_tests/test_legacy_conversion.py`.

## 1. Decision

**Line-backed saves move to sparse quads by a one-way conversion.** A converted world is an
ordinary quad world: it saves as one, steps as one, and keeps no line state. On the nine
representative saves in §5, conversion loses no field or provenance and passes every
invariant audit. It keeps crust volume within 1.2% and land within about 1 pp; seven saves
meet every tolerance and two miss one narrowly. There is no quad-to-line path, and no
long-lived line loader for legacy saves after #251.

Until the cutover (#250), line worlds remain a supported surface, so conversion is opt-in:

| stage | a line-backed save loads as | how to convert |
|---|---|---|
| now (before #250) | a line world, unchanged | `POST /world/load?convert_lines=true`, `persistence.load_world_bytes(data, convert_lines=True)`, or `bin/debug/convert_legacy_saves.py --write-converted` |
| after #250 (quads are the default) | flip `convert_lines` to default on: a quad world, converted on load | (automatic) |
| after #251 (`PlateWithLines` deleted) | a quad world, converted on load through `legacy_conversion.legacy_unpickler` | (automatic) |

The converter reads line state structurally (§3), never through `PlateWithLines` methods, so
it outlives the line classes. #251 removes the line *engine*, not `legacy_conversion.py`,
`LegacyRecord` and `legacy_unpickler`; it switches `load_world_bytes` to unpickle a
line-backed save with `legacy_unpickler` and always convert it.
`test_conversion_reads_saves_whose_line_classes_are_retired` proves conversion through the
stand-ins gives exactly the same world as through the real classes.

## 2. Save format

A save is a pickled envelope. `persistence.SAVE_FORMAT_VERSION` is 3.

| version | shape | loads? |
|---|---|---|
| 1 | a bare pickled `World` (every save before the envelope existed) | yes; may be a line world |
| 2 | `{"format": "mantle-bloom-world", "version": 2, "world": World}` | yes |
| 3 | v2 plus `"surface": "lines" \| "quad" \| "empty"`, checked against the unpickled plates | yes |
| > 3 | written by a newer build | no: `UnsupportedSaveVersionError` |

Each `PlateWithSparseQuadPatch` also versions its own pickled state
(`QUAD_SURFACE_FORMAT_VERSION`, now 1). A quad plate from a newer build makes the whole save
`UnsupportedSaveVersionError`.

**Errors.** `load_world_bytes` raises one of two `SaveFormatError`s (a `ValueError`) and never
loads part of a save:

- `CorruptSaveError`: not a readable mantle-bloom save. That covers bytes pickle can't read
  (garbage, truncation, a class this build doesn't have), another program's envelope,
  anything other than a `World` inside, a declared surface that doesn't match the plates, a
  world mixing line and quad plates, or a quad plate whose cells break the leaf-topology
  invariants.
- `UnsupportedSaveVersionError`: an envelope version outside 1–3 (including a non-integer), or
  a quad plate in another surface format version.

`POST /world/load` returns `400` with `invalid or incompatible world file: <reason>` for
either.

**Fields a save predates.** A `World` field added after a save was written is backfilled on
load (`persistence._backfill_added_fields`: a `default_factory` field gets its empty default;
a plain-default field falls through to the class attribute). An `ElevationLine` field added
since reads as its registry default (`ElevationLine.__getattr__`, matching
`surface_fields.SURFACE_FIELDS`). The eustatic water budget, if missing or in its pre-#257
units, is re-snapshotted from the save's own hypsometry and sea level, so the shoreline
doesn't move. Attributes an older build set on `World` that this build no longer declares
(found in real saves: `gap_fill_algorithm`, `volcanic_field_plate_ids`, the old
`ocean_water_column_m`) are inert; the last is dropped.

**Fields from a newer build.** A save from a newer build has a newer envelope version and is
refused. Conversion separately refuses line state it doesn't recognise (§3) instead of
dropping it.

## 3. What conversion keeps, drops and recomputes

| state | conversion |
|---|---|
| plate id, frame, omega, nominal crust type, age, internal stress | kept exactly |
| every `SURFACE_FIELDS` field | transferred by its `RemapClass` (§4) |
| line positions (`phi`, `theta` per row) | replaced by cells; topology is now the cell set |
| `_leading_row_retreat_years` (a row-end retreat tracker, in 33 of 51 real saves) | dropped: it has no quad meaning |
| topology/geometry revisions, outline/row-lookup/KD-tree caches | dropped (derived) |
| `collision_progress`, `overlap_progress`, `pinned_omegas` (plate-id keyed) | kept |
| faults and fault systems (plate id + plate-local trace) | kept: frames are unchanged |
| collision evidence and front records (plate id + plate-local points) | kept: frames are unchanged |
| `collision_polarity_frame` (this step's per-node masks) | dropped; the next step's prepass rebuilds it |
| magma parcels, earthquakes, gap/stranded-basin tracks, removed-points and corner-notch logs (world-space) | kept |
| `stats_history`, `phase_budget`, `events` | kept; the conversion is logged as an event |
| climate, hydrology, erosion and node-index caches | dropped; recomputed on the next use |
| `ocean_water_volume_m3` | re-snapshotted against the converted hypsometry at the same sea level |
| `World.surface_conversion` | new: the conversion summary, kept through later saves |

**Refused rather than dropped.** A plate attribute or line field not in the tables above,
non-finite field values or line coordinates, field arrays whose length disagrees with their
line, and malformed plate motion (a frame that isn't a finite proper rotation, a non-finite
or misshapen omega, a non-finite internal stress, a negative or fractional age, an unknown
crust type) all raise `ValueError`. That's the "no silent field loss" guarantee for state added by a build this
converter predates.

**Not carried, by design.** The line engine's random streams are keyed by line index
(docs/plate-surface-baseline.md §4.4), so a converted world's future differs from the line
world's even where its state matches. Comparisons are statistical, as in #247/#249.

## 4. Method

`legacy_conversion.py`'s docstring has the detail. In short:

- **Node areas.** Each node stands for the nominal footprint at the world's spacing, shared
  with the same-plate nodes stacked within half a spacing of it (issue #230's stacked rows).
  These deduplicated areas weight every average and are the reference for volume deltas.
- **Territory.** Each plate gets a level-0 lattice in its own frame at the world's spacing. A
  cell is active when its centre lies in some node's footprint (0.6 spacings), or in a
  pinhole enclosed by nodes on every side within 1.5 spacings. The plate must also win the
  area-weighted vote of the 8 nearest nodes of any plate. So no cell centre is claimed twice,
  interleaved overlap rows don't fragment into salt-and-pepper cells, edges facing open sphere
  move by about half a cell either way, and real gaps stay gaps.
- **Fields.** A node goes to its own cell when that is active, otherwise to the nearest
  active cell. A cell combines its own plate's nodes by remap class, with the same rules
  `coarsen_cells` uses: area-weighted mean for extensive thicknesses and intensive values,
  votes for categories (structural `elev_change_reason` codes win ties), OR for volcano
  provenance, max for countdowns and channel depth, and earliest for history. Elevation is
  isostasy from the new column plus the carried residual. A cell no node reaches copies its
  nearest node, but not its volcano.
- **Overlaps.** A node whose target cell is on another plate is crust the line world held
  twice. It stacks onto that cell the way `quad_merge` stacks a suture, every field by its
  remap class. Volume adds and provenance and history combine as above. Composition is voted
  by crust volume, and other categories by area with ties keeping the cell's own. Intensive
  and clock fields blend by area (soil contents by soil depth, channel width by channel
  depth), and so does the elevation residual.
- **Caps.** Thickness that stacking pushes past the Hc/Hm caps diffuses outward across the
  receiving plate's cells, through full ones, for up to 12 rings (~750 km at the default
  density), so a converted collision thickens a belt instead of one column. Hc/Hm are then
  clamped into their caps (the suture-accretion limit), shifting elevation isostatically.
  This clamp is the only place volume leaves. It also clamps any column a pre-#256 save
  holds past the caps, which the line engine itself does at the end of its next step (#262);
  the report separates the two (`source_over_cap` versus `clamped`).

## 5. Tolerances and results

`ConversionReport` (in `World.surface_conversion` as a summary, in full from
`convert_world_to_quad`) records each quantity below. `bin/debug/convert_legacy_saves.py`
converts saves, audits the result with the #247 hard-invariant audit, round-trips it, steps
the line and converted worlds side by side, and checks these tolerances:

| check | tolerance | why |
|---|---|---|
| Hc transfer delta: quad volume against the line volume clipped to the caps, not counting suture delamination | ≤ 1% | crust is the conserved quantity |
| Hm transfer delta, same basis | ≤ 1.5% | Hm sits at its caps more often |
| suture delamination: Hc (Hm) the clamp removes beyond the save's own over-cap volume | ≤ 1% (1.5%) | overlap crust stacked past the cap even after spreading |
| worst plate Hc volume, plates ≥ 500 nodes not touched by an overlap (< 1% of nodes stacked to or from) | ≤ 5% | each edge moves up to half a cell; an overlap legitimately moves crust between plates |
| land share of the covered sphere at the save's sea level, sampled | ≤ 1 pp | |
| continental share of the covered sphere, sampled | ≤ 1 pp | |
| uncovered sphere, quad minus line | ≤ +0.5 pp | conversion must not open gaps |
| multiply-covered sphere (quad) | ≤ 2.0% (warn above 1.5%) | rotated lattices overlap at edges by construction; the #247 harness's C3 levels |
| sea level at conversion | unchanged | the water budget is re-snapshotted |
| #247 hard-invariant audit, at conversion and after each step | no violations | |
| save/load round trip | identical authoritative state hash | |
| provenance: every volcano node's cell is a volcano, every active countdown is kept, no creation time gets later | 0 violations | write-once and provenance fields |

Volumes use the deduplicated node areas of §4. Area shares are sampled the same way on both
sides (each of 200k sphere points reads its nearest node or cell centre within 1.5
spacings): per-node area sums would count a line overlap twice.

Exact on the synthetic shapes in `test_legacy_conversion.py` (isolated plates with a hole, a
concave notch, the antimeridian, a pole, and two disconnected fragments): every point well
inside a shape converts to its own plate, everything more than 2.5 spacings from a node stays
open, the hole, both fragments and both polar and antimeridian plates survive, nothing is
stacked or clamped, and crust volume is within 2%.

### Representative saves

Nine of the 51 real saves, chosen to span age (0–1,064 Myr), density (one at 0.5, the rest
at 4), fragmentation (the 352.4 Myr #228 reproducer and the most fragmented save, 363 Myr with
10,947 lines) and legacy damage (saves holding up to 7.6% of their crust past today's caps).
Each converted world was audited, round-tripped and stepped twice (100 kyr) beside its line
original. Results are in [`analysis/issue248/`](../analysis/issue248/): the inventory of all
51 saves, one JSON per conversion, `summary.md`, and before/after renders.

| save | Myr | line nodes | quad cells | Hc Δ | Hm Δ | Hc over cap in save | Hc delaminated | worst untouched plate Hc Δ | stacked nodes | land line→quad | continental line→quad | uncovered line→quad | multiply covered line→quad | sea level after steps, line / quad | result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|---|---|
| seed745352920-0y | 0.0 | 130578 | 129210 | +0.00% | +0.00% | 0.00% | 0.00% | 0.08% | 0.11% | 0.290 → 0.289 | 0.766 → 0.766 | 0.85% → 0.30% | 1.52% → 0.31% | -0 / -3 | pass |
| seed125087475-11600000y | 11.6 | 130649 | 124570 | -0.01% | +0.11% | 0.00% | 0.01% | 0.56% | 0.66% | 0.331 → 0.332 | 0.655 → 0.655 | 1.98% → 0.99% | 0.26% → 0.52% | -222 / -258 | pass |
| seed19335323-33200000y | 33.2 | 16557 | 15386 | -0.63% | -0.36% | 0.29% | 0.63% | 0.32% | 2.04% | 0.251 → 0.254 | 0.292 → 0.293 | 5.09% → 1.60% | 1.17% → 1.57% | -871 / -874 | pass (warn: multiply_covered) |
| seed343559903-104700000y | 104.7 | 110465 | 104430 | -1.15% | -0.16% | 0.14% | 1.15% | 0.22% | 3.00% | 0.204 → 0.211 | 0.698 → 0.698 | 20.80% → 20.49% | 2.25% → 0.75% | -327 / -1082 | FAIL: crust_delamination |
| seed331006609-239700000y | 239.7 | 163177 | 130845 | -0.18% | -0.12% | 3.85% | 0.17% | 0.42% | 0.97% | 0.269 → 0.279 | 0.368 → 0.370 | 1.96% → 0.77% | 0.34% → 0.55% | -1943 / -2036 | FAIL: land_fraction_delta |
| seed804913535-352400000y | 352.4 | 158589 | 129602 | -0.04% | -0.06% | 0.00% | 0.04% | 0.48% | 1.14% | 0.329 → 0.331 | 0.420 → 0.421 | 2.32% → 0.83% | 0.47% → 0.71% | -820 / -816 | pass |
| seed875551829-363100000y | 363.1 | 229650 | 129889 | +0.00% | +0.02% | 0.03% | 0.00% | 0.13% | 1.20% | 0.132 → 0.134 | 0.042 → 0.042 | 4.21% → 1.64% | 0.12% → 1.03% | -3256 / -3251 | pass |
| seed579428537-626500000y | 626.5 | 222068 | 130346 | +0.04% | +0.08% | 7.62% | 0.00% | 0.27% | 1.45% | 0.195 → 0.196 | 0.043 → 0.043 | 4.28% → 2.01% | 0.24% → 1.03% | -3161 / -3158 | pass |
| seed896200538-1063700000y | 1063.7 | 151091 | 126176 | -0.39% | -0.01% | 1.61% | 0.39% | 1.33% | 1.90% | 0.136 → 0.138 | 0.177 → 0.177 | 11.55% → 4.13% | 0.57% → 0.89% | -2486 / -2475 | pass |

Every save passes the exact checks: no audit violations at conversion or after either step,
identical round trips, and no provenance violations. Sea level is unchanged at conversion. It
then tracks the line world within a few tens of metres, except at 104.7 Myr: that save is 20%
open sphere, which quad gap filling covers with new ocean floor on the first step (§6).

Two saves miss a tolerance, narrowly, and they are recorded as failures rather than
re-tuned away:

- **104.7 Myr, Hc delamination 1.15% (limit 1%).** 3% of this save's nodes sit in overlaps,
  the most of any save, and their crust lands on columns already near the cap. The overflow
  spreads 12 rings before the rest is removed.
- **239.7 Myr, land share +1.03 pp (limit 1 pp).** The per-plate conservation ratio thins
  or thickens a few small plates by up to 7.5%, which moves coasts. Crust volume itself is
  within 0.2%.

The policy decision does not depend on these: conversion keeps every save's crust within
1.2% and its land within about 1 pp, with no field or provenance lost.

Rerun:

```sh
backend/.venv/bin/python bin/debug/inventory_legacy_saves.py ~/Downloads/*.mbworld --out analysis/issue248/save_inventory.json
backend/.venv/bin/python bin/debug/convert_legacy_saves.py SAVE.mbworld [...] --out analysis/issue248/conversions --steps 2 --render --jobs 3
backend/.venv/bin/python bin/debug/convert_legacy_saves.py --out analysis/issue248/conversions --recheck   # after a tolerance change
```

### Inventory of real saves

All 51 `.mbworld` files on hand (`analysis/issue248/save_inventory.json`) are version 1 (bare
pickles), hold only `LithospherePlate`s, and load in this build. Older ones lack up to 44
`World` fields (all backfilled) and up to four line fields (`crust_type_code`,
`node_created_years`, `elev_change_reason`, `overlap_onset_years`; all defaulted). No line
carries a field this build doesn't know. The only line-only plate state is
`_leading_row_retreat_years` (33 saves), which conversion drops (§3).

## 6. Known limits

- Edges move by up to half a cell, so small plates (long edges for their area) see larger
  relative area changes than big ones.
- A notch or channel narrower than about three spacings is enclosed by nodes on both sides,
  so it fills. In line saves these are almost always row stubs, the artefact #228 removes.
- Where line plates overlapped, the overlap is resolved at conversion: the vote picks one
  owner and the other plate's crust stacks onto it, capped. Plates deep in an overlap can
  gain or lose a large share of their own crust that way, though the world total is kept.
- Real voids in a line save stay voids, and the quad engine's every-step gap filling (#259)
  fills them on the first step. A save with a large void will see sea level move on that
  step as new ocean floor appears. The report's `coverage.after.uncovered` shows how much to
  expect.
- `corner_notch_log` entries describe line-row events; they are kept as history but have no
  quad counterpart going forward.
