# Issue #589: G1-qualified endpoint fast path

## Current decision

Keep `lakes.solve_lakes_prognostic` as the only production solver. The refreshed 10 kyr screen
found a subset of isolated leaf bodies that meet both G1 bounds, but they contain only 0.0047%
of pooled measured starting lake volume. On those same bodies, endpoint evaluation and budget-
bound checks took 2.72× one production C integration slice, before ledger bookkeeping. There is
no measured runtime savings to justify activating this path.

This follows issue #589's fail-closed condition: keep the optimization disabled when it fails
the runtime gate. The current implementation remains C-default, as RFC 0014's owner chose after
the earlier G1 diagnosis.

## Refreshed G1 evidence

`gate_lake_endpoint.py` captures current production C liquid forcing and the following D2(c)
lake-ice phase for one frozen-terrain 10 kyr step on each checksum-registered diagnosis save.
Every active body above 10,000 km² is measured, plus a seeded 5% sample of smaller bodies.
Merged bodies, leaves that can receive hierarchy transfers, spilling bodies, invalid cases and
unsampled small bodies are assigned to C fallback. For the remaining isolated leaf bodies, the
screen compares endpoint equilibrium with a direct C trajectory, then applies the same lake-ice
phase to both liquid outcomes. This is an offline measurement; it does not create a production
eligibility cache.

The 1,000-substep direct C trajectory was checked against the prior 4,000-substep result for 47
common bodies on seed 910211954: maximum relative final-volume difference was 5.7e-9 and no
eligibility classifications changed. All screened eligible bodies met the 1% liquid transient
bound; their largest secondary combined liquid-plus-ice difference was 0.0053%. The per-save
volume shares, runtime measurements and checksums are in
[`g1-d2c-summary.json`](g1-d2c-summary.json). Full per-body outputs and saves remain in scratch
space.

| Seed | Eligible bodies | Measured starting-volume coverage | Eligible starting-volume share | Endpoint / C slice time |
|---:|---:|---:|---:|---:|
| 910211954 | 228 | 100.000% | 0.00032% | 2.68× |
| 635056485 | 39 | 100.000% | 0.02390% | 3.13× |
| 317739385 | 96 | 99.999% | 0.00720% | 2.65× |
| **Pooled** | **363** | **>99.999%** | **0.00465%** | **2.72×** |

The body-time comparison is endpoint root solving and budget-bound checks versus one production
`_integrate_body_volume` call on the same eligible bodies. It excludes ledger bookkeeping, so
the measured ratio understates endpoint-path cost. Captured full 10 kyr C step times are
reported in the JSON for context; a density-4 origin/main versus endpoint implementation replay
was not run because the candidate path already costs more per eligible body and does not meet
the runtime-savings gate.

## Prior evidence and D2(c) correction

- The earlier G1 replay reported 46.71% pooled union-flagged lake-volume share and supported C
  as default, with endpoint equilibrium only considered as an optimization for fast bodies in
  the [#538 owner decision](https://github.com/adubey/mantle-bloom/issues/538#issuecomment-6394597033).
- That replay used zero sublimation and said to rerun if D2 changed the loss law
  ([#538 G1 result](https://github.com/adubey/mantle-bloom/issues/538#issuecomment-6394305581)).
- Current D2(c), merged through PR #660, applies sublimation to the separate lake-ice reservoir
  after liquid routing. The refreshed diagnostic captures this phase after the production C
  solve and compares both liquid-only and combined liquid-plus-ice outcomes.

## Scope left inactive

No production solver change is made. Endpoint budget booking, nested and asymmetric routing
support, randomized hierarchy comparisons, runtime gains and held-out acceptance were not
implemented or run because the diagnosis shows no runtime case for enabling the optimization.
The production lake water ledger and C fallback remain unchanged.
