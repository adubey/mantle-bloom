# Issue 894: why the damage-band prototype peels edge slivers

Exploration for [#894](https://github.com/adubey/mantle-bloom/issues/894). It measures the
RFC 0018 D1 prototype (`bin/debug/rift_damage_step.py`, PR #890) without changing it, and
proposes a direction. Nothing here is a model change or acceptance evidence.

## Conclusion

The prototype does not fail because loads sit at plate edges. It fails for three reasons:

- **A. The failure test is per cell, but break-up is a property of a whole section.** The load
  is nearly uniform across a plate at about a tenth of rock strength. Only the weakest few
  percent of cells exceed their strength, and those are thin oceanic cells and plate margins.
  Damage therefore follows weak rock, not the line the plate is being pulled apart along.
- **B. Softening spreads damage instead of concentrating it.** A cell's load is set by the
  plate's forces and does not fall when a neighbour yields. Softening lowers the strength of
  every cell within 125 km, so neighbours cross the threshold and the damaged patch grows
  outward. Nothing limits the stretch across a section or makes a band grow along its strike.
- **C. The path and split rules reject the cuts that remain.** A third of path searches are
  dropped for overlapping legs, stray pieces of 1–27 cells fail the two-component rule, and the
  700-cell daughter minimum rules out a 70/30 split on any plate under 2,333 cells.

The total pull across straight sections of a plate does peak across the plate's middle, and
stays there for tens of Myr. That is the signal a rift algorithm can use to choose a large
split.

## Runs

All runs are plates-only (`simulate_climate_biomes = False`), 1 Myr steps, default knobs,
development seeds only, prototype at `4ce4f5cf`.

| Probe | Seeds | Length | Command |
|---|---|---|---|
| 1 | 106533960 | 60 Myr | `probe_rift_edge_peeling.py --seed 106533960 --years 60 --out probe1.json --shell-at 1 30` |
| 1b | 106533960 | 60 Myr | the same with `--maps DIR --map-at 10 30 60` (adds path outcomes and band shape) |
| 2 | 1097624328, 1097710202, 1397817354, 1439649077 | 1 Myr | `probe_rift_edge_peeling.py --seed S --years 1 --out probe2-seedS.json` |

Probe 1 and 1b gave identical counters (567 ligament passes, 0 splits, 4,481 withheld
thinnings), matching the run reported in the issue. Predictions were posted on the issue before
each run. Summaries: `probe1-seed106533960-summary.json` (from
`summarize_rift_edge_peeling.py`) and `probe2-sections-summary.json`. Raw outputs and maps are
outside the repository under the `issue894` captures directory.

"Hops" is the number of cell-to-cell steps from the nearest plate-boundary cell. One hop is
about 62.5 km.

## Findings

### 1. Load is flat; strength is low at edges

Seed 106533960, seven snapshots pooled and cell-weighted:

| Hops from edge | 0 | 1 | 2 | 3–4 | 5–9 | 10+ |
|---|---:|---:|---:|---:|---:|---:|
| Median load, proxy C (TN/m) | 0.50 | 0.50 | 0.50 | 0.50 | 0.53 | 0.58 |
| Median load, thin shell (TN/m) | 0.25 | 0.32 | 0.41 | 0.50 | 0.67 | 0.68 |
| 99th percentile load, proxy C (TN/m) | 1.26 | 1.09 | 1.05 | 1.02 | 0.96 | 0.89 |
| Median strength (TN/m) | 4.3 | 4.9 | 5.6 | 6.7 | 8.2 | 9.7 |
| Median load / strength | 0.14 | 0.12 | 0.11 | 0.10 | 0.10 | 0.16 |
| Cells loaded past strength | 5.8% | 5.2% | 4.7% | 3.8% | 2.4% | 3.0% |
| Cells with an inherited weak mark | 31% | 27% | 19% | 10% | 8% | 8% |
| Share of all cells | 12% | 10% | 10% | 17% | 31% | 20% |
| Share of straining cell-steps | 20% | 10% | 7% | 11% | 25% | 27% |

The prediction that edge load would be at least twice interior load was wrong. Strength is
2.2 times lower at the edge than at 10+ hops. On the four probe 2 seeds the load is also flat
or rises inward (edge to interior medians 0.43 to 0.69, 0.22 to 0.34, 0.42 to 0.44 and 0.33 to
0.54 TN/m) and strength is 1.3 to 2.7 times lower at the edge.

Slab pull supplies 0.28–0.32 TN/m along the tension direction at every depth. Basal drag
supplies 0.22 TN/m at the edge and 0.47 TN/m at 10+ hops. Ridge push is negligible.

### 2. Strain is mostly oceanic, and it reaches interiors

- 66% of straining cell-steps are oceanic (78% at the snapshots). Broken bands that reach the
  path search are 89% oceanic on average.
- 22% of straining cells are within 3 hops of a divergent boundary, against 12% of all cells.
- 27% of straining cell-steps are at 10+ hops, above that zone's 20% share of cells. The
  prediction of under 1% was wrong. More time does not turn this strain into lines.

### 3. Damage forms patches, not lines

Maps at 10, 30 and 60 Myr show broken cells as wide patches: about 15 by 8 cells (roughly
900 by 500 km) in the continental interior of plate 13, and broad oceanic patches on plates 9
and 16. Bands of 50 or more cells have a median elongation of 3 (long axis over short axis),
and the long axis lies a median 31° from the tension direction, where a rift would strike at
90° to it. This is the #699 failure mode in a milder form.

### 4. Path and split rules

Outcomes of all 1,680 path searches:

| Band size (cells) | Searches | Legs overlap, skipped | Ligament rule fails | Ligament rule passes |
|---|---:|---:|---:|---:|
| 2–9 | 1,157 | 362 | 504 | 291 |
| 10–49 | 294 | 137 | 26 | 131 |
| 50+ | 229 | 70 | 14 | 145 |

- Most searches start from tiny bands on the boundary: the median band has 3 cells at 1.5 hops,
  and its path ends are 7% of the plate's span apart.
- 569 searches (34%) are dropped because the two Dijkstra legs share cells, including all 27 on
  plate 13's continental interior band. The decile endpoint rule assumes a line-shaped band.
- 42 of the 567 ligament-passing paths would leave a second piece of at least 15% of the plate:
  plate 16 at 7–22 Myr (for example 1,265 / 640 cells, a 67/33 split) and plate 9 at 18–47 Myr
  (about 1,000 / 270 cells, along a 70-cell path whose ends are 80% of the span apart). They
  were rejected for stray pieces of 1–27 cells and for the daughter minimum.
- `SPLIT_MIN_NODES = 700` is applied to each daughter, as production does. A 70/30 split then
  needs a parent of at least 2,333 cells. At 30 Myr, 4 of 18 plates are that large; the median
  plate has about 1,700 cells.

### 5. Whole sections are pulled hardest across the middle

For each plate and each of proxy C's six directions, the pull across every straight section
was summed and compared with the section's summed strength.

| Measure | Seed 106533960 (79 plate-snapshots, 60 Myr) | Probe 2 seeds (1 Myr, plates of 1,400+ cells) |
|---|---|---|
| Section of largest total pull has at least 20% of the plate on its smaller side | 72% | 67%, 67%, 80%, 67% |
| Median smaller-side share at that section | 33% | 39%, 35%, 36%, 27% |
| Best interior section's load over its intact strength, median | 0.13 | 0.40, 0.22, 0.10, 0.40 |
| Section of largest pull-to-strength ratio is an edge sliver (under 10% share) | 84% | not summarized |

- The largest ratio sits at plate tips because a short chord has almost no strength to sum. A
  ratio test alone therefore also peels slivers.
- Probe 2's prediction for the interior ratio (0.05 to 0.3) was missed on two seeds, which came
  out higher at 0.40.
- The pull axis is persistent: between consecutive 10 Myr snapshots it stayed within 30° in 47
  of 56 cases, and seven plates kept their pull peak at 29–46% share for 50 Myr.

## Verdicts on the issue's questions

| Question | Verdict |
|---|---|
| 1. Load geometry | No. Load is flat. Strength is low at edges. |
| 2. Double-counting at ridges | Partly. Strain is mostly oceanic; ridge proximity is a smaller effect. |
| 3. Path construction | Yes, a prototype defect: overlapping legs drop 34% of searches. Tiny bands on the boundary explain the 2-cell paths. |
| 4. Sparse inheritance | Supported. 8% of interior cells carry a weak mark, against 31% at the edge. |
| 5. Timescale | No. Interior strain is common; it does not localize. |
| 6. Process mismatch | Partly. Margin and oceanic weakness, not slab rollback, drives the peeling. |
| 7. Strain spreading within a step | Not tested. Finding 3 suggests spreading across steps matters more. |

## Earth comparison

- Rifting takes a long slow phase, then accelerates abruptly before rupture as rift strength
  decays (Brune et al. 2016, *Nature* 536, 201–204, doi:10.1038/nature18319). In the proposal
  below this comes from the intact part of a rift line shrinking.
- In the western branch of the East African Rift, rifting localizes in Proterozoic mobile belts
  at the margin of the Tanzanian craton (Koptev et al. 2018, *Geophysical Research Letters*;
  read as an abstract). Rift lines should follow weak belts and pass around cratons.
- Lithospheric buoyancy forces are reported as the primary driver of East African block motion,
  with mantle tractions needed along the rift (Rajaonarison, Stamps and Naliboff 2021,
  *Geophysical Research Letters*; read as an abstract). The model has basal drag and plumes but
  no buoyancy force from uplift.
- Back-arc opening, which does peel margins, depends on the mode of subduction and is not a
  consequence of subduction alone (Uyeda and Kanamori 1979, *JGR* 84(B3), 1049–1061,
  doi:10.1029/JB084iB03p01049). Margin peeling is a real but separate process.
- Buck 1991 (*JGR*, "Modes of continental lithospheric extension") is cited by the RFC and the
  issue; its details were not re-verified here.

## Proposed direction: rift lines

Choose the line first from whole-plate pull, then let it grow from the bottom up. This reuses
the RFC's load proxy, strength, damage fields, thinning and ligament rule. It replaces the
per-cell strain test (step 3) and the band-first path search (step 6).

1. **Pull axis.** From proxy C's six sorted running sums, take the total pull across every
   straight section and find the direction and position where it peaks. This costs nothing
   beyond proxy C (about 2 ms per plate).
2. **Corridor and line.** The corridor is the strip of sections whose total pull is within a
   set fraction of the peak. Inside it, one Dijkstra search finds the cheapest path from one
   side of the plate to the other, with cost equal to strength times length. Belts, old rift
   scars, plume-heated and thinned cells are cheap; cratons are expensive. The path wiggles
   around strong blocks, so it is not a straight line.
3. **Section load.** `Λ` is the pull across the line divided by the summed strength of its
   intact cells, which is the RFC's ligament rule read as a ratio.
4. **Nucleation.** When `Λ` passes a calibrated threshold, a rift starts at the weakest stretch
   of the line, or at a plume on it.
5. **Tip growth.** Each step the two tips advance along the line at a rate that rises with
   `Λ`. As the intact part shortens, `Λ` rises and growth speeds up: slow, then fast. Only line
   cells strain, to the RFC's width `w_0`, and they thin with the existing flank transfer.
6. **Failure.** If the pull falls or the tips reach rock too strong to cross, growth stops. The
   warm part heals and the scar stays cheap for a later attempt.
7. **Break-up.** When both tips reach the plate boundary or another rift, the plate splits
   along the line. A simple boundary-to-boundary path gives two daughters; small stray pieces
   join the adjacent daughter.

Cost per plate is one path search on a corridor when a rift is active, plus proxy C. The
capacity table, disk averages and per-cell strain pass are needed only on line cells.

### What this does not settle

- Evidence is five development worlds, and the growth rule is untested. This is a hypothesis.
- 20–33% of plates have their pull peak at an edge. They would not rift under this rule.
- Lines may still run through oceanic crust on mostly oceanic plates. The RFC's oceanic
  strength factor is the existing answer.
- Six straight directions are a coarse description of pull on large curved plates.
- A daughter minimum relative to parent area, and attaching strays, are changes to kept
  production rules and need an owner decision.
- The model has no buoyancy load from plume uplift, which is a main driver in East Africa.
