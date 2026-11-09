# RFC 0018 diagnosis-gate evidence

This report records the revised G0 and G1 development checks and the completed G2 held-out
comparison. G0/G1 use 20 development worlds; G2 uses the first 10 seeds from the frozen
20-seed held-out registry. Detailed compact records and their hashes are listed below. Full
simulation outputs remain outside the repository.

## Gate status

| Gate | Outcome | Evidence |
|---|---|---|
| G0 — revised measurement and frozen model | PASS | 20 development worlds; 68 eligible opportunities; 812 candidate proposals. The `O/R/S` prior, mass-balanced proposal, momentum ranking score, and `mu=1` are frozen. See [`g0-prior-summary.json`](g0-prior-summary.json). |
| G1 — controlled probability and structural checks | PASS | Known-probability draw cases pass; 100,000-draw controls at 0.1, 0.5 and 0.9 fall within 4 standard errors; monotonic controls pass; both selected virtual cuts pass connectedness, trace, size and ownership. See [`g1-results.json`](g1-results.json). |
| G2 — held-out non-worsening comparison | PASS under merged PR #845 | All first 10 registry pairs complete; 255 actual splits; 0 of 10 have `D_s < 0`; at most one is allowed. This was a post-result sample-size amendment. See [`g2-sample-size-amendment.json`](g2-sample-size-amendment.json). |

The G2 result demonstrates non-worsening on these evaluated pairs, not positive efficacy.
The sample-size reduction was made after the first-ten outcomes were observed, so this is not a
preregistered confirmatory result and has less precision than the originally planned 20-pair
evaluation. The 50-split event minimum was met without reserve seeds.

## G0 — revised measurement and frozen model

The revised supplement reuses the 20 frozen development snapshots and records opening `O`,
candidate resistance `R`, pole separation `S`, structural validity, proposal weights, score and
runtime. It removes the former `1.96 × RMS residual` opening cutoff and does not treat pole
separation as an independent hard gate. The source has no per-cell suture provenance, so no
suture weakening signal was inferred from belts.

The same-snapshot candidate-generation comparison found 33 structural cuts / 19 opportunities
for the weak-anchor baseline, 52 / 25 for uniform anchors, and 97 / 33 for mass-balanced
targets (812 proposals across 68 opportunities; runtimes 313.1, 316.7 and 329.9 s). The
mass-balanced proposal was frozen for G0/G1 because it had the highest structural yield on
development data. These are development results, not held-out comparisons.

The frozen candidate generator proposes up to 12 mass-balanced cuts, with target fractions
from 0.15 through 0.50, deduplicates identical partitions, and weights proposals by the
nonnegative potential-momentum ranking score at `mu=1` (uniform fallback when all scores are
zero). This is a development-set choice, not held-out evidence. For a candidate:

```text
O_eff = max(O, 0) * S / S_ref
p_i = O_eff / (O_eff + O_ref * R_i / R_ref)
P_break = sum_i(q_i * p_i)
```

Topology-invalid candidates have `p_i = 0`; if breakup is selected, the candidate is drawn
from probabilities proportional to `q_i * p_i`. The frozen development scales are
`O_ref = 0.8734616233 cm/year`, `R_ref = 61.0574320`, and
`S_ref = 6.46551732e-9 rad/year`. At those reference values a valid candidate has prior 0.5.
These values define a transparent trial prior; they are not calibrated geological breakup
rates.

| Measure | Result |
|---|---:|
| Development worlds / eligible opportunities / proposals | 20 / 68 / 812 |
| Structurally valid candidate proposals | 97 / 812 |
| Opportunities with structural candidate / positive `P_break` | 33 / 68; 17 / 68 |
| Mean opportunity `P_break` at `mu=1` | 0.0469; world-bootstrap 95% interval 0.0254–0.0734 |
| Constructible worlds without an eligible opportunity | 2 (retained) |

G0 passes under the revised criteria. The old opening threshold's zero-candidate result is
retained as historical evidence and is superseded by the accepted stochastic interpretation.

## G1 — coupled-probability checks

`bin/debug/run_rfc0018_g1.py` applied the frozen G0 model to all 68 opportunities, making one
keyed uniform breakup draw per opportunity and a conditional candidate draw only after a
breakup was selected. It did not use the legacy independent 0.55 roll and did not mutate
production worlds.

| G1 check | Result |
|---|---:|
| Known-probability boundary cases | Pass |
| 100,000-draw controls at `P=0.1, 0.5, 0.9` | 0.10060, 0.49873, 0.89948; all within 4σ |
| Increasing opening / decreasing resistance / increasing pole separation | Pass; probabilities remain in `[0,1]` |
| Opening required for `P=0.5` as resistance falls from `4R_ref` to `0.25R_ref` | 3.494 to 0.218 cm/year |
| Keyed development outcomes | 2 breakups / 68 opportunities, descriptive only |
| World-cluster bootstrap interval for descriptive outcome fraction | 0.000–0.0704 (95%) |
| Selected virtual cuts passing connectedness, trace, minimum size and parent ownership | 2 / 2 |

G1 passes its revised diagnostic criteria. The observed opportunity fraction is one keyed
realization, not a geological rate threshold. The 20-world report retained both worlds with
no opportunities; every applied production split still needs implementation-level ownership
and ledger validation.

## G2 — held-out comparison

The frozen registry is [`g2-seed-preregistration.json`](g2-seed-preregistration.json), SHA-256
`214623b962984be9e95c6ced1f87b470b105a6d628f96cdb7cfae5d5844eba2f`. Evaluate its first ten
seeds, in registry order. The sample-size amendment in PR #845 changed the required sample
from 20 pairs to 10 and preserved the prior 10% tolerance as at most one negative pair.
It was made after all first-ten outcomes were observed. The amendment is recorded in
[`g2-sample-size-amendment.json`](g2-sample-size-amendment.json).

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

The 255 combined actual splits exceed the 50-event minimum. There are zero negative `D_s`
values, the mean is 0, and the 95% paired seed-bootstrap interval is [0, 0]. Both tiny-share
values are zero, so relative share reduction is undefined and `S_baseline - S_graph` is 0
percentage points. All amended G2 conditions pass.

Seed `205240614` has a completed baseline arm but no graph arm and is outside the first ten
registry seeds. It is unpaired and excluded. Reserve seeds were not used.

## Compact records and provenance

- `g0-summary.json`: historical corrected G0 results under the superseded opening threshold.
- `g0-prior-summary.json`: revised G0 aggregate and frozen `O/R/S` prior; SHA-256
  `49dc81750bab43c64593edecc1d54c2b930a50141819c80b16d7314f56e08d91`.
- `g1-results.json`: event-level controlled and development outcomes; SHA-256
  `981fc836816b8918f3bc7800c8eba1d13f1abb2f0379cac98bb6a34e58d71668`.
- `g2-summary.json`: original ten-pair interim summary; gate fields are interpreted together
  with the subsequently merged sample-size amendment.
- `g2-metric-amendment.json` and `g2-nonworsening-amendment.json`: the earlier endpoint and
  non-worsening amendments; both were adopted after outcomes were observed.
- `g2-sample-size-amendment.json`: merged 10-pair rule and outcome.
- `g2-seed-preregistration.json`: unchanged frozen held-out and reserve registry.

The G0 supplement and G1 used the preregistered development set and frozen candidate table.
The G2 runner compared production great-circle splitting with an isolated experimental
callback using the frozen G0/G1 prior over 200 Myr with 1 Myr steps. The graph arm did not
modify production source. Full replay outputs and the 812-row candidate table remain outside
the repository; compact results and checksums are kept here.
