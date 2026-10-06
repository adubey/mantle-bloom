# Issue #249 parity campaign: final report

The #249 campaign rerun on `main` at `2bba51d`, after the fixes for #253, #255, #256, #257
and #259. The first campaign, in `analysis/issue249-campaign/`, is the reference.

**Recommendation: no-go for the #250 cutover**, until #254 and #268 are resolved. Every hard
invariant holds in every run, the #228 artifacts are gone, and quad is faster than lines at
the #147 profile settings. Two problems remain, both on the quad side:

- **#254:** the quad deep ocean floor sits about 300–560 m shallower at 400 Myr. As a result, quad sea
  level ends up about 900 m higher by 400 Myr, and quad has half the line land fraction.
- **#268:** at node density 0.5, quad double-claims 1.9–2.2% of the sphere. That reaches the
  2.0% C3 ceiling and crosses it at several checkpoints.

## Sets and commands

All commands run from `backend/` with `.venv/bin/python ../bin/debug/surface_parity.py`. Each
set's `provenance.json` records the exact command, commit and platform.

| set | what | command | verdict |
|---|---|---|---|
| `long` | seeds 1–5, density 1, to 400 Myr | `paired --preset long --seeds 1,2,3,4,5 --jobs 5` | fail |
| `repro-equivalent` | the 352.4 Myr reproducer's seed on both surfaces, to 352 Myr, with renders | `paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render` | fail |
| `stress` | seeds 1–2, density 0.5, to 800 Myr | `paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2` | fail |
| `repro-lines` | the real 352.4 Myr line save, continued 1 Myr under every audit | `run --preset issue147 --seed 804913535 --surface lines --from-world <save> --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render` | pass |
| `issue147` | the #147 profile world (seed 0, density 4, climate 4, fluid 2, 100 kyr steps, to 60 Myr), sequential, quiet machine | `paired --preset issue147 --seeds 0 --jobs 1` | fail (K4 only) |

The save for `repro-lines` (`mantle-bloom-seed804913535-352400000y.mbworld`) isn't in the
repo. There is no line-to-quad converter yet (#248), so `repro-equivalent` stands in for the
quad side of the reproducer.

The `run_*.sh` scripts are the exact invocations used. `summary.md` is
`bin/debug/summarize_parity_campaign.py ../analysis/issue249-final`, and has the full gate
and metric tables. As in the reference campaign, the per-step timings and logs aren't
committed; `phase-means.json` keeps the per-phase means.

The `long`, `repro-equivalent` and `stress` sets ran with several worlds in parallel, and
`long` shared the machine with unrelated test runs. Their timings only mean something as
quad/line ratios within a set. Only `issue147` is a clean timing run.

## #228 quality gates

| #228 acceptance criterion | result | evidence |
|---|---|---|
| No folded quads or invalid neighbours in long stress runs | **pass** | H1–H3 pass at every audit (every 10 steps, plus checkpoints) in every set, including 800 Myr at density 0.5. |
| Boundary loops and intervals are reproducible derived views that can't silently drift | **pass** | H6 (derived caches equal fresh rebuilds, exact KD-tree parity) and H7 (revisions) pass at every audit in every set. |
| Holes, concavity, fragments, antimeridian and polar plates represented correctly | **pass** | H1/H2 pass with up to 13 hole loops per world (`issue147`). X1 antimeridian, pole and plate-boundary pass everywhere. M6 warns for 1–3 extra fragments, which is not a correctness failure. |
| Crustal volume conserved within a documented tolerance | **pass on quad**; K1/K2/K4 fail as paired gates (see below) | K3 (Σ cell area = sphere − uncovered + overlap) passes everywhere. Quad keeps more Hc than lines in every set. |
| Persistent fields survive remeshing with appropriate semantics | **pass** | H4 (fields finite and sized), H11 (caps) pass in every set. H8/H9 save/load and continuation are identical. |
| The 352.4 Myr reproducer no longer shows stretched-line/row-stub artifacts | **pass** | `repro-equivalent/renders/`: the 352 Myr line elevation map has swirling row streaks across the oceans, and the line plate map has one-row stubs. The quad maps have neither. Stacked nodes: 0 on quad against 32% on lines. M2 (row alignment) passes. |
| Coverage and overlap invariants at least as strong as lines | **pass at density 1 and 4; fail at density 0.5 (#268)** | Uncovered sphere 1.5% on quad against 8.2% on lines (long, 400 Myr); C1, C2, C4 pass everywhere. C3 is under 1.6% at densities 1 and 4, but reaches 2.23% at density 0.5. |
| Performance measured against #147, no unaccepted dominant-phase regression | **pass, with one accepted regression** | See the next section. |
| Rendering/projection remain downstream and representation-neutral | **pass** | Renders are produced by the same path for both surfaces. X1 HEALPix lookups pass, except where samples are too few (holes). |
| Climate and hydrology stable | **pass** | S1 and every S2 jitter gate pass in every set. Quad sea-level jitter is lower than lines after #257/#259. |
| Ensemble parity of world statistics (P1) | **fail: elevation p05 and sea level (#254)** | Long set, 5 seeds: land fraction, Hc, continental Hc, plates, elevation p50/p95 pass. Elevation p05 is 372 m shallower on quad at 60 Myr (limit 255 m); sea level is 703 m higher at 240 Myr (limit 573 m). |

### Conservation gates K1, K2, K4

These gates fail when quad's and lines' drift differs by more than 10 pp, whichever surface
drifts more. Broken down by direction:

- **K1 (total Hc):** in every set and every failing checkpoint, lines lose more. For example,
  the long set at 400 Myr averages −20% on quad against −42% on lines. This fits #254's
  finding that line endpoint growth thins the line ocean floor. It is a line baseline loss,
  not a quad regression.
- **K4 (Hm):** lines drift more almost everywhere, in either direction. In `issue147` both
  surfaces grow Hm during the cooling transient, by 23% on quad and 42% on lines. This is
  also a line baseline finding.
- **K2 (continental Hc):** mixed. Lines lose more early on, but late in a run quad loses more
  on three seed-runs:

  | run | quad | lines |
  |---|---:|---:|
  | long seed 1, 400 Myr | −38% | −26% |
  | long seed 5, 400 Myr | −38% | −12% |
  | stress seed 1, 800 Myr | −78% | −58% |

  The long-set means at 400 Myr are −53% on quad against −44% on lines. Before #253 they were
  −69% against −46%. P1 continental Hc passes. #253 closed most of the gap; what's left is a
  late-run excess on some seeds, which this report records as a tuning difference rather
  than a blocker.

## Performance against #147

`issue147` ran one world at a time on a quiet machine. The #147 profile ran under cProfile,
so compare its column by ratio only.

| s/step | lines | quad | quad / lines | #147 profile (cProfile) |
|---|---:|---:|---:|---:|
| step total | 8.451 | 3.767 | **0.45** | 10.26 |
| deform + topology | 1.562 | 1.364 | **0.87** | — |
| — deform | 1.339 | 0.545 | 0.41 | 2.58 |
| — topology | 0.142 | 0.085 | 0.60 | 0.20 |
| — gap fill | 0.039 | 0.665 | 16.9 | — |
| — overlap tracking | 0.041 | 0.069 | 1.66 | — |
| magma transport | 4.721 | 0.446 | 0.09 | — |
| faults | 0.628 | 0.685 | 1.09 | 1.16 |
| shift | 0.402 | 0.315 | 0.78 | 0.81 |
| climate, erosion, hydrology | 0.832 | 0.742 | 0.89 | 1.07 |
| sea level | 0.213 | 0.179 | 0.84 | 0.24 |

- The #147 dominant phases are no slower on quad. Deform is 0.41× and the deform + topology
  bucket is 0.87×.
- **Accepted regression: gap fill.** It costs 0.66 s/step on quad (17× lines) and takes most
  of the deform saving. It runs every step since #259, which fixed sea-level jitter, and the
  bucket as a whole is still faster than lines. It's the obvious next target if quad
  deformation needs to get faster.
- Magma transport, the largest #147 hotspot, costs 4.7 s/step on lines but 0.45 s/step on
  quad at this density. That accounts for most of the 0.45× step total. At density 1 the two
  surfaces are about equal. Why quad is cheaper at density 4 hasn't been investigated, so
  don't count on that saving until it's explained.
- Faults are 9–21% slower on quad across the sets. That is within R1/R2 tolerance.

Indicative ratios from the parallel sets: long step total 0.68× and deform + topology 0.45×;
stress 0.72× and 0.46×.

## Changes since the reference campaign

| | reference (`f345cc5`) | final (`2bba51d`) |
|---|---|---|
| H11 field caps | warn (line side) | pass (#256 clamp) |
| C3 multiply covered, long | fail | warn (1.58% worst) |
| C3 multiply covered, stress | fail (2.01% worst) | fail (2.23% worst), #268 |
| S2 sea level / land stability | warn | pass |
| land fraction at 400 Myr, long mean (quad / lines) | 0.060 / 0.137 | 0.077 / 0.136 |
| continental Hc drift at 400 Myr, long mean (quad / lines) | −69% / −46% | −53% / −44% |
| sea level at 400 Myr, long mean (quad / lines) | −1395 / −2058 m | −1120 / −2006 m |
| step total, long (quad / lines) | 0.61× | 0.68× |

## Issues

| issue | status | effect on cutover |
|---|---|---|
| #253 quad retreat destroys continental crust | fixed | residual late-run excess on some seeds (K2 above), tuning |
| #254 quad ocean floor ~430 m shallower | open | **blocker**: P1 elevation p05 and sea level fail, and land fraction is about half of lines by 400 Myr. The latest finding on #254 points at line-side under-thinning. If that holds, the quad floor may be the more correct one. Either way, the cutover needs either a fix or an explicit decision to accept the higher sea level. |
| #255 quad double-claims 1.2–1.9% | fixed (2.0% ceiling) | holds at density 1 and 4 |
| #256 Hm leaves its caps | fixed | H11 passes everywhere |
| #257 nominal-area sums | fixed | S2 passes |
| #259 gap-fill cadence | fixed | sea-level jitter gone; gap-fill cost accepted above |
| #268 C3 at density 0.5 | open | **blocker** if density 0.5 is a supported setting, otherwise a tolerance decision |
