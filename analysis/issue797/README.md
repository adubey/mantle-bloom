# RFC 0018 G2 evidence

G2 compares the production great-circle split method with the experimental graph-cut
method on paired held-out worlds. A tiny split has a smaller daughter below 7% of the
parent plate's area. For seed `s`, `D_s = T_baseline,s - T_graph,s`; negative values mean
more tiny splits with the graph method, and zero is a tie.

## Current gate and amendment

The 2026-11-08 owner-directed sample-size amendment evaluates the first ten seeds, in order,
from the unchanged frozen 20-seed registry. The gate requires all ten pairs to be complete,
no more than one pair with `D_s < 0` (preserving the prior 10% tolerance), and at least 50
combined actual splits. The complete amendment record is in
[`g2-sample-size-amendment.json`](g2-sample-size-amendment.json); the frozen registry is
[`g2-seed-preregistration.json`](g2-seed-preregistration.json).

This amendment was made after all first-ten outcomes were observed. It is post-result evidence,
not a preregistered confirmatory analysis. It supports non-worsening on the evaluated pairs,
not positive efficacy, and the smaller sample has less precision.

| Seed | Baseline splits | Graph splits | Tiny baseline | Tiny graph | `D_s` |
|---:|---:|---:|---:|---:|---:|
| 228191914 | 22 | 5 | 0 | 0 | 0 |
| 2058294723 | 20 | 14 | 0 | 0 | 0 |
| 2069607991 | 15 | 10 | 0 | 0 | 0 |
| 1305368628 | 9 | 9 | 0 | 0 | 0 |
| 857560705 | 15 | 15 | 0 | 0 | 0 |
| 640625437 | 20 | 10 | 0 | 0 | 0 |
| 2137727007 | 17 | 10 | 0 | 0 | 0 |
| 1471988496 | 16 | 9 | 0 | 0 | 0 |
| 1249837662 | 11 | 8 | 0 | 0 | 0 |
| 1416208252 | 12 | 8 | 0 | 0 | 0 |
| **Total** | **157** | **98** | **0** | **0** | **0** |

Across the ten pairs, 255 actual split events meet the 50-event requirement. No pair has
`D_s < 0`; the mean is 0 and the 95% paired seed-bootstrap interval is [0, 0]. Both tiny
shares are zero, so the relative share reduction is undefined and the percentage-point
difference is 0. Under the amended criteria, all three gate conditions pass, conditional on
adoption of the amendment.

Seed `205240614` has a completed baseline arm but no graph arm and is outside the first ten
registry seeds. It is unpaired and excluded from this gate. The unmodified full interim
summaries, including any runs beyond the amended sample, remain in the private run directory.
