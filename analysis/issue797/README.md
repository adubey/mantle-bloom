# RFC 0018 G0 — development measurement

## Scope and status

RFC 0018 is Accepted by the project owner's instruction dated 2026-11-03. Tracking issue
[#797](https://github.com/adubey/mantle-bloom/issues/797) owns G0–G3. This worktree is based
on `origin/main` at `63290b053de15b72280fe28fe0086ff37affdba2` and is dedicated to issue #797.

G0's corrected diagnostic run completed on 2026-11-04. It did **not pass G1**: none of 68
eligible opportunities had a candidate whose net opening exceeded the frozen residual
envelope. The current source also lacks the RFC's per-node suture signal, so this run cannot
establish a pass for the full stated resistance composition. The development seeds remain
preregistered below; no held-out seeds have been selected for G2. G0 must
record every eligible opportunity, candidate, rejection, daughter
mass share, cut length and continuity, opening and uncertainty, score, and runtime. Keep all
20 primary seeds, including constructible worlds with no eligible split opportunity, in the
report. If the RFC's inconclusive-result rule applies, use only the 30 separately predeclared
expansion seeds to reach 50 total. Never use G2 held-out seeds for tuning.

## Principle and measure

The current great-circle split choice conflicts with design principle 8 (let each process
follow its own physics): it selects a geometry from mantle-flow disagreement without using
local weakness or checking outward opening. G0 measures whether graph-cut candidates can
produce connected, dominant edge-to-edge traces with outward opening above a fit-derived
uncertainty, while recording mass-weighted daughter shares and candidate scores. This is a
development gate; it does not establish lower tiny-breakup rates or held-out performance.

## Recovered development evidence

The RFC's three-world pilot is recorded in
[`../exploration-graph-rift-algorithm/README.md`](../exploration-graph-rift-algorithm/README.md)
and its compact JSON. It made one production-checked split in three worlds (40.2% smaller
daughter, +0.013 cm/yr net opening) and two no-splits. These seeds and results are development
evidence and are excluded from G0 and G2.

The probe script was found in the existing #770 worktree. Its SHA-256,
`c6be15d742663a9dd508e9d032394bd2970661ce7658f3efef55280ee1aa089b`, matches the hash in
all three compact pilot records. It was written against model commit
`96b5b8c833ad98696b7a2759ff93fcfb4ba93d64`; it is not in `origin/main` and its inputs are
local saves that the RFC says are unavailable for download. The script is therefore a
source for reconstructing the diagnostic, not current-source G0 evidence. Reconcile it with
the accepted RFC's finalized resistance composition and momentum score before any gate run.

## G0 seed preregistration

Selected before any G0 run, on 2026-11-03. Seed IDs are 31-bit positive integers derived by
SHA-256 from the UTF-8 label
`mantle-bloom/rfc0018/g0/v1/{group}/{index:02d}`: interpret the first four digest bytes as a
big-endian unsigned integer and mask to 31 bits. Reject zero, duplicates, and IDs found in
the issue #770 JSON evidence or the three RFC 0018 pilot records. The #770 artifact scan used
the evidence tree at commit `96b5b8c833ad98696b7a2759ff93fcfb4ba93d64`; it found 28 unique
historical IDs including the pilot seeds. The selected 55 IDs are unique and do not overlap
that set. Candidate outcome data from these selected seeds has not been inspected; the failed
initialization batch's baseline split counts were observed and are excluded from tuning.

| Group | Index labels | Seeds | Use |
|---|---|---|---|
| Primary | `primary/00`–`primary/19` | 1629738962, 1097710202, 591453549, 969845080, 763144949, 870996356, 494568917, 1854897345, 954941394, 600370317, 1643146570, 146046339, 1397817354, 1630166320, 604997962, 704718994, 1097624328, 487506346, 106533960, 1439649077 | The 20 initial G0 development worlds |
| World-generation replacements | `replacement/00`–`replacement/04` | 1227913777, 303885793, 2045331463, 157959777, 1208093067 | In order, only when a primary seed cannot construct a valid baseline world before stepping. A valid world with no eligible opportunity remains in the sample and is not replaced. |
| Expansion | `expansion/00`–`expansion/29` | 58813960, 1533142689, 139601706, 986078511, 1862603747, 1195410342, 651078283, 295468763, 1348581301, 591259338, 1343302795, 1083351774, 483800127, 1419040997, 1725816656, 986489725, 165063606, 1248030031, 1871340641, 1404866888, 637746806, 847386583, 543848301, 913181899, 988144353, 846786010, 1061237420, 1048620058, 882673937, 618341104 | Thirty additional worlds, used only if the RFC's G1/G0 uncertainty rule requires expanding to 50 total. |

The G0 and G2 seed sets remain disjoint: no G2 list has been selected yet, and none of these
development IDs may be reused for held-out acceptance runs.

## Corrected G0 results

The hash-frozen diagnostic completed all 20 constructible primary worlds at the first
age-qualified 2 Myr snapshot. The 19-world 200 Myr initialization attempt is excluded: its
historical suture-mask calls failed at every candidate before any candidate was measured.
Raw corrected output is preserved outside the repository at
`~/mantle-bloom-captures/rfc0018-g0/rfc0018-g0-results-corrected.json`; its SHA-256 is
`da65c6bceff7bd4881a69b79034e4a19a47b6bb23196f3d13ae4c51f20634086`. Compact aggregate data
is in [`g0-summary.json`](g0-summary.json).

| Measure | Result |
|---|---:|
| Constructible worlds | 20 / 20 |
| Eligible opportunities | 68 |
| Candidate cuts measured | 812 |
| Candidates meeting connectedness, trace, daughter-size and pole-separation checks | 31 |
| Candidates with positive net opening | 398 / 812 |
| Candidates with net opening above `1.96 ×` fit-residual RMS | 0 / 812 |
| Opportunities with at least one fully usable candidate | 0 / 68 |
| World-bootstrap 95% interval for usable-opportunity fraction | 0.000–0.000 |
| Total runtime / median per world | 430 s / 22.4 s |
| Peak process RSS | 627 MB |

The best candidate's net-opening-to-threshold ratio was 0.494; among the 31 candidates that
passed the other daughter/trace checks, the maximum was 0.481. Thus no candidate was close to
the frozen opening gate. Seventeen opportunities had at least one structurally usable
candidate. For those, the momentum score with `mu=1` and the prototype score selected the same
candidate in 12/17 opportunities; `mu=0.5`, `1`, and `2` selected the same top candidate as
`mu=1` in 17/17, 17/17, and 16/17, respectively. These rankings are descriptive only; no
score or `mu` is frozen because G1 failed.

The 95% world-bootstrap interval is wholly below the two-thirds threshold, with 68
opportunities (above the ten-opportunity minimum). Under the RFC's rule this is a clear G1
failure: do not expand to 50 and do not tune resistance weights. Revise the opening-
uncertainty definition before rerunning gates. Separately resolve the missing suture signal in
the source/RFC design. The run used the fixed suture multiplier with no active cells and must
not be represented as validation of that term.

## Frozen G0 protocol

- Generate each world with `world.generate_world(seed)` and no overrides: the current
  production defaults for plate count, land, mantle centers, node/climate/fluid density,
  cratons-first layout, relief, tilt and climate remain in force. A world is constructible
  when generation returns a non-null world with at least one plate. Constructible worlds
  with no split opportunities are retained. Use listed replacements only for a generation
  failure before stepping.
- Advance each world to 2 Myr in 1 Myr `step_world` calls, with default production movement
  and climate controls. At the first age-qualified split-decision step, measure every plate
  that passes the current residual trigger from the same pre-split state. The 2 Myr horizon
  reaches the existing minimum plate age; G0 measures split geometry and does not need the
  long-history horizon reserved for G2. The production path still runs for that step after
  measurement, but no later trajectory is included in G0. Record wall runtime and maximum RSS.
- An opportunity is a plate that passes the existing minimum-node and minimum-age checks and
  whose current mantle-flow residual is at least the existing size/stress-adjusted residual
  trigger in `maybe_split_plate`, before k-means, cut refinement, daughter pole, and minimum
  daughter-size checks. Record later production gate outcomes separately. This operational
  denominator captures the point where production first attempts a split decision.
- Build each shared-edge length from the exact spherical endpoints of its edge segment. For
  a coarse cell adjacent to two fine cells, count the two half-edge segments separately;
  aggregate segments that share the same node pair. Thus cut length/opening weights use the
  mesh interface, not center-to-center distance.
- For a candidate partition, fit each daughter rotation to its own mantle-flow vectors. On
  every crossed edge oriented A-to-B, project the difference between the two fitted surface
  velocities onto the edge tangent for signed opening. Fit residuals are projected the same
  way: `d_e = R * (residual_B - residual_A) dot tangent`. Define the fit envelope as the
  shared-edge-length-weighted RMS of `d_e`; the candidate passes the opening gate only when
  length-weighted net opening exceeds `1.96 * RMS`. This is a conservative residual envelope,
  not a calibrated confidence interval. Record both values in cm/yr. This rule is frozen
  before looking at any selected-seed outcomes.
- Evaluate the RFC's candidate-specific motion-resistance term against each proposed side
  assignment; keep all tabulated structural resistance values and the 12-candidate count
  fixed. Compare the RFC potential-momentum ranking with the prototype's reduced-mass times
  net-opening ranking, and evaluate `mu = 0, 0.5, 1, 2`. Do not let the diagnostic mutate the
  production split path.
- Sampled anchors are weighted by spectral separation divided by edge resistance. A stable
  SHA-256 exponential-race key over seed, plate ID, elapsed simulation time and the two stable
  cell keys selects up to 12 anchors without replacement. Candidate-specific fit, resistance
  and spectral cut are iterated up to eight times, stopping early when the side assignment is
  unchanged. The spectral solver uses a deterministic starting vector; record the iteration
  count. The diagnostic uses the prototype's 250 km fault
  trace reach and defines the thinning signal as areal column mass below the plate median.
- The current source has no per-node suture or pseudo-suture mask (only dated belt history).
  Do not infer a suture from every belt. Keep the RFC's suture multiplier fixed but inactive
  for this diagnostic and record `suture_signal_available=false`. This means candidate scores
  diagnose the available current-source signals; the full stated resistance composition needs
  a source representation or explicit RFC revision before G0 can be called a pass.
- The dominant trace is measured by joining crossed-edge midpoints within 0.75 nominal cell
  spacing, weighting each component by exact shared-edge length. Its edge contacts are
  counted as connected components among plate-boundary cells. Record this operational trace
  rule with the measured share and contact count.

The source commit is `63290b053de15b72280fe28fe0086ff37affdba2`; the corrected diagnostic
script SHA-256 is `3e5a489f3725285115b159251d0c6b52cec8e2b14b8a2752d3605f8f932fff15`. The first
instrumentation batch completed 19 worlds, but
all 27,355 candidate probes failed at the first resistance calculation because the historical
prototype's `cratons.suture_mask`/`pseudo_suture_mask` APIs do not exist at this source commit.
Those trajectories are baseline data only and are excluded from G0 candidate measurements.
The seed schedule, candidate rules and uncertainty formula remain unchanged. The corrected
run uses a 2 Myr first-age-qualified snapshot instead of 200 Myr: RFC G0 asks for geometry
and score measurability, and its long paired-history horizon belongs to G2. At the diagnostic
step, all eligible plates are measured from the same pre-split state, avoiding irrelevant
long-history cost and omission after the production path's first successful split. A corrected,
hash-frozen rerun is required.

## Next steps

1. Revise the RFC's fit-uncertainty rule and specify how the missing suture signal is
   represented before changing any weights.
2. Rerun the preregistered G0/G1 gates only after that revision; keep the held-out G2 sample
   unselected and separate.
