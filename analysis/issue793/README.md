# Even continent sizes at the start: concentration 2 against 32

Findings for [#793](https://github.com/adubey/mantle-bloom/issues/793).

**Principle:** 10, because a starting layout whose smallest continent is a sliver and whose
largest holds half the continental area gives the planet's history little to work with.

**The change.** `CONTINENT_SIZE_CONCENTRATION` goes from 2 to 32. The Dirichlet draw that splits
continental area between the 4-7 continents stays unequal and every continent is substantial.

## Measures

`bin/debug/continent_layout_measures.py --seeds 40 --concentration 2 32`, seeds 1-40, land
fraction 0.29, node density 1 (`layout-measures.json` has the summaries and the core share on
seeds 7, 11 and 2024).

| Measure over 40 seeds | Concentration 2 | Concentration 32 |
|---|---|---|
| Smallest continent's share of continental area (minimum) | 0.5% | 9.3% |
| Largest continent's share (maximum) | 59.5% | 32.2% |
| Seeds outside 5-40% | 25 | 0 |
| Cratons per world | 27-58 | 30-62 |
| Cratons within 300 km of the continental edge, median (range) | 38% (15-55%) | 41% (11-65%) |
| Belt share of continental area | 11.0-20.4% | 12.3-21.5% |
| Core share of continental area, mean (range) | 39.7% (37.7-40.0%) | 39.8% (38.5-40.0%) |
| Core share on seeds 7, 11, 2024 | 40.0, 37.7, 39.6% | 40.0, 38.5, 39.8% |
| Nuclei that joined a craton, mean | 82% | 69% |
| Cratons of 10 M km² or more per world, median | 1 | 0 |
| Largest craton | 20.9 M km² | 18.2 M km² |
| Craton area 5th percentile; share under 1 M km² | 0.51 M km²; 31% | 0.51 M km²; 30% |

## Reading

- The split is the only thing the change moves directly. The core share is fixed by
  `CRATON_SHARE_OF_CONTINENT` and does not change, and the 5th percentile craton size and the
  share of small cratons are the same.
- Continents are smaller at the top, so the largest craton is smaller and a median world has
  no craton of 10 M km² or more (1 before). The median share of cratons near the edge moves
  3 points and the share of nuclei that join a craton 13 (82% to 69%). Belt share stays
  inside the RFC's 10-30%.
- The comments on the constants in `continent_layout.py` and the RFC 0009 knob table carry
  these figures. Several of the figures they held before (cratons 20-55, belts 11-20%, edge
  share median 35%) no longer match the layout at concentration 2, so they were out of date
  already.
- This is the layout at t=0. No replay of a saved world was run; whether the continents stay
  distinct as plates merge and rift is not measured here.
