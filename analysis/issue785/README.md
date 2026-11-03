# Bound isolated submerged basins in initial synthetic relief (#785)

## Implementation verification

The cap was checked on ten fresh cratons-first synthetic-relief worlds using Generate-dialog
defaults: continental fraction 0.70, land fraction 0.29, node, climate and fluid densities
4, 4 and 2. The worlds were not stepped. Seeds were preregistered in the issue comment:
528019447, 7, 11, 2024, 123456789, 987654321, 42, 314159265, 876543210 and 135792468.
Seed 528019447 was generated a second time for a determinism check.

The diagnostic counts each below-sea continental component disconnected from the open-ocean
mask. Every measured component stayed at or below 371,000 km² and at or above -1,050 m
(allowing 1e-6 m for isostatic floating roundoff). Repeated seed 528019447 produced the same
elevation hash. Sea level was 0 m for every world.

| Seed | Holes | Total hole area (km²) | Largest hole (km²) | Deepest floor (m) | Generated land area (%) |
|---:|---:|---:|---:|---:|---:|
| 528019447 | 2 | 374,080 | 365,688 | -1,050 | 28.91 |
| 7 | 1 | 3,634 | 3,634 | -33 | 29.16 |
| 11 | 10 | 233,407 | 145,795 | -154 | 28.21 |
| 2024 | 0 | 0 | 0 | — | 28.07 |
| 123456789 | 0 | 0 | 0 | — | 29.04 |
| 987654321 | 1 | 25,834 | 25,834 | -51 | 28.91 |
| 42 | 2 | 71,932 | 67,533 | -230 | 29.44 |
| 314159265 | 5 | 50,765 | 30,453 | -116 | 27.26 |
| 876543210 | 5 | 531,692 | 366,715 | -293 | 29.62 |
| 135792468 | 0 | 0 | 0 | — | 29.48 |

Across the ten worlds, total isolated-hole area averaged 129,134 km² (median 50,765 km²),
with a maximum of 531,692 km². The sum can exceed the per-hole cap when a world has multiple
holes; seed 876543210 has five, with its largest at 366,715 km². The maximum total is 1.43
times the Caspian Sea surface area and 2.66 times the ~200,000 km² below-sea footprint
reported for Earth's Caspian Depression. No standardized global aggregate for all Earth's
landlocked below-sea floors was found. The often cited ~20% of land for Earth's endorheic
basins describes entire drainage catchments, including high ground, so it is not comparable
to this below-sea-floor measure.

For seed 528019447, the prior unbounded control had 2,270,391 km² in its largest hole and
land area fraction 0.285411. The implementation produced a 365,688 km² largest hole and land
area fraction 0.289146, an increase of about 0.374 percentage points (about 1.91 million
km²) because the excess area is raised to 5 m. Other seeds' land-area fractions are reported
as generated; this run did not generate paired unbounded controls for all ten seeds.

The tested Earth comparison is only a scale reference, not validation of the chosen cap.
The below-sea area of the Caspian Depression is reported at about 200,000 km²
([Geology.com](https://geology.com/below-sea-level/)); the global endorheic catchment area
covers about a fifth of land ([global endorheic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6267997/)).
The cap itself uses the Caspian Sea's approximate 371,000 km² water-surface area, as chosen
for this implementation.

## Reproduction

From the repository root, run:

```bash
backend/.venv/bin/python bin/debug/measure_initial_basin_bounds.py
```

The script uses the fixed seeds above and generates one extra copy of seed 528019447. It
reports hole areas and floors, land fractions and an elevation hash for determinism.
