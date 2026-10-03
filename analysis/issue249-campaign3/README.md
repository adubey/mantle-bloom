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
`run_campaign.sh` is the exact script used, and each set's `provenance*.json` records the command,
commit and platform.

| set | what | command | verdict |
|---|---|---|---|
| `long` | seeds 1–5, density 1, to 400 Myr | `paired --preset long --seeds 1,2,3,4,5 --jobs 5` | fail |
| `repro-equivalent` | the 352.4 Myr reproducer's seed on both surfaces, to 352 Myr, with renders | `paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render` | fail |
| `stress` | seeds 1–2, density 0.5, to 800 Myr | `paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2` | fail |
| `repro-saves` | the real 352.4 Myr line save, and its #248 quad conversion, each continued 1 Myr under every audit | `run --preset issue147 --seed 804913535 --surface lines\|quad --from-world <save> --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render` | warn (M6 only) |
| `issue147` | the #147 profile world (seed 0, density 4, climate 4, fluid 2, 100 kyr steps, to 60 Myr), one world at a time | `paired --preset issue147 --seeds 0 --jobs 1` | fail (K4, line side); R1 warn from contention |

The line save (`mantle-bloom-seed804913535-352400000y.mbworld`) isn't in the repo. Its quad
conversion comes from `bin/debug/convert_legacy_saves.py <save> --steps 0 --write-converted`; the
log is `convert.log`.

`summary.md` is `bin/debug/summarize_parity_campaign.py ../analysis/issue249-campaign3`.
Per-step timings and logs aren't committed; `phase-means.json` keeps the per-phase means.

**All timings in this campaign are indicative.** The correctness sets ran several worlds in
parallel. `issue147` ran one world at a time, but the mantle-bloom desktop app was running alongside
it the whole time, using about 75% of a core. Treat its numbers as quad/line ratios. The clean
absolute comparison with #147 is campaign 2's (`analysis/issue249-final/README.md`).

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
| Performance measured against #147, no unaccepted dominant-phase regression | **pass** (indicative) | See the performance section. |
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

The `issue147` run in this campaign doesn't give clean timings. The desktop app slowed its quad
half much more than its lines half. The table compares each phase with the same run in campaign 2,
which was one world at a time on a quiet machine:

| s/step | lines, c2 | lines, c3 | slowdown | quad, c2 | quad, c3 | slowdown |
|---|---:|---:|---:|---:|---:|---:|
| step total | 8.451 | 10.187 | 1.21× | 3.767 | 8.380 | 2.22× |
| deform | 1.339 | 1.755 | 1.31× | 0.545 | 1.487 | 2.73× |
| topology | 0.142 | 0.167 | 1.18× | 0.085 | 0.157 | 1.84× |
| gap fill | 0.039 | 0.044 | 1.11× | 0.665 | 1.165 | 1.75× |
| faults | 0.628 | 0.764 | 1.22× | 0.685 | 1.203 | 1.76× |
| shift | 0.402 | 0.564 | 1.40× | 0.315 | 0.975 | 3.10× |
| climate, erosion, hydrology | 0.832 | 0.997 | 1.20× | 0.742 | 1.454 | 1.96× |
| magma transport | 4.721 | 5.488 | 1.16× | 0.446 | 1.471 | 3.30× |

#270 changed only gap filling, yet every quad phase slowed down by 1.5–3.3×, including shift and
climate, which #270 doesn't touch. So the R1 warning here (deform + topology at 1.45× lines) comes
from contention, not code. Gap filling itself slowed down less than the untouched phases, so this
run gives no sign that #270 made it slower.

**The #147 comparison therefore rests on campaign 2** (`analysis/issue249-final/README.md`, PR #269,
`2bba51d`), which was run cleanly:

- whole step 0.45× lines, deform + topology 0.87×, deform alone 0.41×;
- gap fill is the accepted regression, at 0.66 s/step (17× lines);
- most of the whole-step saving comes from magma transport, at 0.09× lines, and that saving hasn't
  been explained yet.

The ratios from the parallel correctness sets agree. Step total is 0.73–0.77× lines and
deform + topology 0.49–0.55× on long, stress and repro-equivalent, against 0.68–0.72× and
0.45–0.46× in campaign 2.

For clean absolute numbers on `fa7ecdf`, rerun the quad half alone on a quiet machine:
`run --preset issue147 --seed 0 --surface quad --out ../analysis/issue249-campaign3/issue147`,
then `compare` that directory.

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
