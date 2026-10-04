# Issue #249 parity campaign: final report (campaign 3)

The #249 parity campaign, rerun on `main` at `fa7ecdf` after the #268 seam-overlap fix (PR #270).
It supersedes `analysis/issue249-final/` (PR #269, run at `2bba51d`). The first campaign, in
`analysis/issue249-campaign/`, remains the reference.

Since PR #269, two other changes have landed:

- #248 added a line-to-quad converter, so this campaign also continues the real 352.4 Myr save on
  quad (`repro-saves`), not only the seed-equivalent run.
- #267 made sparse quads the default for new worlds. The cutover that #249 gates has already
  happened, so this report checks that cutover after the fact.

**Recommendation: go, conditional on resolving #254.** The cutover can stand; no quad-side problem
found here is serious enough to revert it. Every hard invariant holds in every run, the #228
artifacts are gone, #268 is resolved, the real reproducer save converts and runs cleanly, and quad
is faster than lines. What's left:

- **#254 (open):** the quad deep ocean floor sits 355 m shallower at 60 Myr. As a result, quad sea
  level ends up about 600–900 m higher, and by 400 Myr quad has about half the line land fraction.
  The latest finding on #254 points at line-side under-thinning. If that holds, closing #254 means
  accepting a wetter world rather than fixing quad. Either way, the issue needs a decision.
- **Late-run continental crust loss (K2):** quad loses more continental Hc than lines after about
  240 Myr. The gap is slightly wider than in campaign 2. This is a tuning difference, not a
  blocker; see below.

## Sets and commands

All commands run from `backend/` with `.venv/bin/python ../bin/debug/surface_parity.py`.
`run_campaign.sh` and `run_perf_clean.sh` are the exact scripts used, and each set's `provenance*.json` records the command,
commit and platform.

| set | what | command | verdict |
|---|---|---|---|
| `long` | seeds 1–5, density 1, to 400 Myr | `paired --preset long --seeds 1,2,3,4,5 --jobs 5` | fail |
| `repro-equivalent` | the 352.4 Myr reproducer's seed on both surfaces, to 352 Myr, with renders | `paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render` | fail |
| `stress` | seeds 1–2, density 0.5, to 800 Myr | `paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2` | fail |
| `repro-saves` | the real 352.4 Myr line save, and its #248 quad conversion, each continued 1 Myr under every audit | `run --preset issue147 --seed 804913535 --surface lines\|quad --from-world <save> --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render` | warn (M6 only) |
| `issue147-clean` | the #147 profile world (seed 0, density 4, climate 4, fluid 2, 100 kyr steps, to 60 Myr), one world at a time on a quiet machine | `paired --preset issue147 --seeds 0 --jobs 1` | fail (K4, line side) |

The line save (`mantle-bloom-seed804913535-352400000y.mbworld`) isn't in the repo. Its quad
conversion comes from `bin/debug/convert_legacy_saves.py <save> --steps 0 --write-converted`; the
log is `convert.log`.

`summary.md` is `bin/debug/summarize_parity_campaign.py ../analysis/issue249-campaign3`.
Per-step timings and logs aren't committed; `phase-means.json` keeps the per-phase means.

Only `issue147-clean` is a clean timing run. The correctness sets ran several worlds in parallel, so
their timings only mean something as quad/line ratios within a set. `run_campaign.sh` also ran an
`issue147` set, but the desktop app was competing with it for CPU, so it was discarded and rerun as
`issue147-clean`.

## #228 quality gates

| #228 acceptance criterion | result | evidence |
|---|---|---|
| No folded quads or invalid neighbours in long stress runs | **pass** | H1–H3 pass at every audit in every set, including 800 Myr at density 0.5 and the converted save audited every step. |
| Boundary loops and intervals are reproducible derived views that can't silently drift | **pass** | H6 (derived caches equal fresh rebuilds, exact KD-tree parity) and H7 (revisions) pass at every audit in every set. |
| Holes, concavity, fragments, antimeridian and polar plates represented correctly | **pass** | H1/H2 pass everywhere. The converted save starts with 38 hole loops and 39 extra fragments; within 1 Myr it is down to 2 hole loops. M6 warns for 2–6 extra fragments in the generated sets. |
| Crustal volume conserved within a documented tolerance | **pass on quad**; K1/K2/K4 fail as paired gates (see below) | K3 (Σ cell area = sphere − uncovered + overlap) passes everywhere. The #248 conversion changes Hc by −0.04% and Hm by −0.06%. On the 1 Myr continuation of the save, quad Hc drifts −0.14% and line Hc −0.85%. |
| Persistent fields survive remeshing with appropriate semantics | **pass** | H4 and H11 pass in every set. H8/H9 save/load and continuation are identical, including on the 130k-cell converted save. |
| The 352.4 Myr reproducer no longer shows stretched-line/row-stub artifacts | **pass** | `repro-equivalent/renders/`: the 352 Myr line maps show row streaks and stubs; the quad maps don't. Quad has zero stacked nodes. On the converted save (`repro-saves/renders/`), conversion keeps the streaky relief the line world had already built up, but it adds no new artifacts, and the row-stub island cluster near the top centre is gone after 1 Myr. |
| Coverage and overlap invariants at least as strong as lines | **pass** | C1, C2, C4 pass everywhere. C3 passes at densities 1 and 4 (worst about 1.4%) and warns at density 0.5 (worst 1.79%, under the 2.0% ceiling). In campaign 2 it failed at 2.23%; #268 is resolved. |
| Performance measured against #147, no unaccepted dominant-phase regression | **pass, with one accepted regression** | See the performance section. |
| Rendering/projection remain downstream and representation-neutral | **pass** | Both surfaces render through the same path. X1 HEALPix lookups pass, except for warnings on the few hundred samples that fall in holes. |
| Climate and hydrology stable | **pass** | S1 and the S2 sea-level and land jitter gates pass everywhere. S2 air temperature warns on `repro-equivalent` (0.34 against a 0.30 limit); every other set passes it. |
| Ensemble parity of world statistics (P1) | **fail: elevation p05 and sea level (#254)** | Long set, 5 seeds: land fraction, Hc, continental Hc, plates and elevation p50/p95 pass. Elevation p05 is 355 m shallower on quad at 60 Myr (limit 255 m), and sea level is 629 m higher at 240 Myr (limit 573 m). Both margins are slightly smaller than in campaign 2 (372 m and 703 m). |

### Conservation gates K1, K2, K4

These gates fail when quad's and lines' drift differs by more than 10 pp, whichever surface drifts
more. Broken down by direction:

- **K1 (total Hc):** in every failing case, lines lose more. The long set at 400 Myr shows the
  same pattern as campaign 2. This fits #254's finding that line endpoint growth thins the line
  ocean floor, so it is a line baseline loss.
- **K4 (Hm):** lines drift more in almost every case, whether the drift is up during the cooling
  transient or down late in the run. This is also a line baseline finding.
- **K2 (continental Hc):** lines lose more before about 240 Myr; after that, quad loses more on
  several seeds.

  | seed means | reference | campaign 2 | campaign 3 | lines |
  |---|---:|---:|---:|---:|
  | long, continental Hc drift at 400 Myr | −69% | −53% | −57% | −44% |
  | long, land fraction at 400 Myr | 0.060 | 0.077 | 0.067 | 0.136 |
  | stress, continental Hc drift at 800 Myr | −94% | −84% | −88% | −73% |
  | stress, land fraction at 800 Myr | 0.029 | 0.056 | 0.046 | 0.168 |

  The worst seeds are long seed 1 (−54% against −26%), long seed 4 (−68% against −41%) and stress
  seed 1 at 800 Myr (−87% against −58%). Campaign 3 is about 4 pp worse than campaign 2 on the seed
  means, but the seeds move in both directions: any change to gap filling sends quad down a
  completely different trajectory. Five seeds can't separate an effect of #270 from noise. P1
  continental Hc passes. This report records K2 as a late-run tuning difference that needs its own
  investigation, not as a cutover blocker.

## Performance

`issue147-clean` ran one world at a time on a quiet machine. The #147 profile ran under cProfile,
so compare its column by ratio only.

| s/step | lines | quad | quad / lines | campaign 2 quad / lines | #147 profile (cProfile) |
|---|---:|---:|---:|---:|---:|
| step total | 8.518 | 4.068 | **0.48** | 0.45 | 10.26 |
| deform + topology | 1.583 | 1.360 | **0.86** | 0.87 | — |
| — deform | 1.358 | 0.535 | 0.39 | 0.41 | 2.58 |
| — topology | 0.144 | 0.086 | 0.60 | 0.60 | 0.20 |
| — gap fill | 0.040 | 0.671 | 16.9 | 16.9 | — |
| — overlap tracking | 0.042 | 0.068 | 1.60 | 1.66 | — |
| magma transport | 4.733 | 0.791 | 0.17 | 0.09 | — |
| faults | 0.636 | 0.621 | 0.98 | 1.09 | 1.16 |
| shift | 0.411 | 0.313 | 0.76 | 0.78 | 0.81 |
| climate, erosion, hydrology | 0.844 | 0.766 | 0.91 | 0.89 | 1.07 |
| sea level | 0.215 | 0.180 | 0.84 | 0.84 | 0.24 |

- The line timings are within 1% of campaign 2's in every phase. That confirms the machine was
  quiet, and that nothing merged since `2bba51d` changed the line path's cost.
- The #147 dominant phases are no slower on quad. Deform is 0.39× lines and the deform +
  topology bucket is 0.86×.
- **Accepted regression: gap fill.** It costs 0.67 s/step on quad (17× lines) and takes most of
  the deform saving. It has run every step since #259, which fixed sea-level jitter. #270 didn't
  change its cost (0.665 s/step before, 0.671 s/step now), and the bucket as a whole is still
  faster than lines. It's the obvious next target if quad deformation needs to get faster.
- Magma transport, the largest #147 hotspot, costs 4.7 s/step on lines but 0.79 s/step on quad.
  That accounts for most of the 0.48× step total. In campaign 2 it was 0.45 s/step; the quad
  trajectory differs after #270, and magma transport cost depends on the world's state. At
  density 1 the two surfaces cost about the same. Why quad is cheaper at density 4 hasn't been
  investigated, so don't count on that saving until it's explained.

The ratios from the parallel correctness sets agree: step total 0.73–0.77× and deform + topology
0.49–0.55× on long, stress and repro-equivalent.

## Changes since campaign 2 (PR #269)

| | campaign 2 (`2bba51d`) | campaign 3 (`fa7ecdf`) |
|---|---|---|
| C3 multiply covered, long | warn (1.58% worst) | pass (about 1.4% worst) |
| C3 multiply covered, stress | fail (2.23% worst) | warn (1.79% worst) |
| uncovered at 400 Myr, long | 1.2–1.7% | 1.6–1.8% |
| P1 elevation p05 / sea level | 372 m / 703 m over the limits of 255 m / 573 m | 355 m / 629 m |
| continental Hc drift at 400 Myr, long mean (quad / lines) | −53% / −44% | −57% / −44% |
| land fraction at 400 Myr, long mean (quad / lines) | 0.077 / 0.136 | 0.067 / 0.136 |
| M6 extra fragments, worst | 1–3 | 6 |
| real 352.4 Myr save on quad | not possible | converts and runs; warn on M6 only |

Lines results are identical in both campaigns: no merged change touches the line path, so the line
runs act as a fixed control.

## Issues

| issue | status | effect on cutover |
|---|---|---|
| #253 quad retreat destroys continental crust | fixed | some late-run excess remains (K2 above); tuning |
| #254 quad ocean floor shallower | open | **needs a decision**: P1 elevation p05 and sea level fail. Either fix the floor, or confirm the line side is the one that's wrong and accept a higher sea level. |
| #255 / #268 quad double-claimed area | fixed | C3 within its ceiling at every density tested |
| #256 Hm leaves its caps | fixed | H11 passes everywhere |
| #257 nominal-area sums | fixed | S2 passes |
| #259 gap-fill cadence | fixed | sea-level jitter is gone; gap-fill cost is accepted (performance) |
| #248 legacy conversion | fixed | the 352.4 Myr save converts within tolerance and continues cleanly |
