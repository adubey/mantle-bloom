# Mass-weighted graph-cut rift prototype

This diagnostic implements a first, offline version of the proposed graph-cut rift selector.
It reads saved rift-onset worlds, leaves them unchanged, and evaluates candidate cell cuts.
It does **not** replace `merge_split.maybe_split_plate`, create daughter plates, or replay the
post-split simulation. The default great-circle behavior is untouched.

The probe follows the simulation design principles of area-weighting finite cells (#4),
giving rifting a process-specific opening and weakness score (#8), and retaining per-seed
evidence and uncertainty (#9).

## Prototype algorithm

1. Use the plate's exact edge-sharing cell adjacency and cluster mantle-flow velocities into
   two regimes.
2. Assign dimensionless edge resistance from edge length, column areal mass (`Hc` and `Hm`),
   cratonic crust, residual from either fitted flow regime, and available weakening signals:
   thinning, generic belts, sutures, volcanoes, plume heating, and active or historic fault
   traces. Coefficients are exploratory heuristics, not calibrated strengths. There is no
   per-cell stress field in this first pass.
3. Approximate low-cost balanced partitions with a weighted Fiedler vector. Sample weak edges
   that cross the two flow clusters; each sampled edge sets a threshold that the resulting
   graph cut crosses.
4. Keep candidates whose two sides are each at least 95% connected in the plate adjacency
   graph. Approximate the cut trace by clustering crossed-edge midpoints; require a dominant
   trace to reach the plate edge at two locations. These are topology-quality checks, not a
   daughter-mass ratio rule.
5. Fit a rigid Euler rotation to each proposed daughter, measure signed opening normal to
   every crossed edge, and compute `reduced_mass * length_weighted_mean_signed_opening`.
   Reject net-closing cuts and rank the remaining cuts by that mass-weighted opening score.

## Saved-world pilot

The three snapshots are the saved production rift-onset states used in the earlier #770
diagnostics. Seeds were not selected as a validation sample for this graph algorithm. The
maps show land, water, cratons, belts, pseudo-sutures, the parent plate edge in gold, and the
selected graph cut in red. A rejected cut is shown as a diagnostic candidate, not as a
resulting split.

| Seed / time / plate | Candidate result | Smaller daughter mass | Net across-cut opening | Opening edge length | Cut trace |
|---|---|---:|---:|---:|---|
| 1173383782 / 3 Myr / plate 1 | Rejected: net closing | 22.1% | -3.85 cm/yr | 3.8% | 1 component; 17,289 km |
| 1183015207 / 15 Myr / plate 8 | Accepted by prototype gates | 28.4% | +0.23 cm/yr | 74.9% | 2 components; dominant trace 97.0%, 15,921 km |
| 1486436815 / 34 Myr / plate 4 | Accepted by prototype gates | 18.4% | +0.11 cm/yr | 59.1% | 2 components; dominant trace 95.9%, 3,120 km |

“Accepted” here means the offline geometry and opening gates passed. It is not evidence that
the full tectonic simulation can realize or sustain the split. The 17,289 km closing cut and
the secondary trace fragments show that candidate topology still needs review. The reduced
mass score helped rank cuts but did not itself prevent the 18.4% daughter, so it should remain
a ranking term rather than be described as a guarantee against small fragments.

## Evidence limits

The three input saves and rendered images are local files and are not included in this
proposal. The scripts ran in the issue #770 worktree and are not part of this `origin/main`
proposal branch. The saved-world names, measured values and image checksum are retained here
so the result can be checked without treating it as acceptance evidence.

## Applied graph-cut split replay

The follow-up replay used the same three saved onset worlds, advanced each for up to 20 Myr,
and applied the graph cut as an actual plate-mask split on the first successful event. It
monkeypatched the split hook only inside the diagnostic process; it did not change the
production great-circle implementation or the baseline generator. Production minimum-node,
plate-age, and daughter-pole-separation checks were retained. The stochastic 0.55 split
success gate was also retained. In this probe, a failed rift does not thin an aulacogen band.

| Seed / target plate | Evaluations | Result | Split time | Smaller daughter mass | Net opening | Cut length |
|---|---:|---|---:|---:|---:|---:|
| 1173383782 / 1 | 17 | Split | 21 Myr | 40.2% | +0.013 cm/yr | 17,002 km |
| 1183015207 / 8 | 20 | No split | — | — | — | — |
| 1486436815 / 4 | 20 | No split | — | — | — | — |

The first split produced daughter plates 1 and 9, with about 50.7% of crossed-edge length
opening. The 40.2% smaller daughter is a fairly balanced result without an explicit mass
fraction floor. The opening is only slightly positive, however, so this event should be
treated as a diagnostic demonstration rather than a calibrated physical success. In the
second world, proposed daughters repeatedly failed the production minimum-size or pole
separation checks; the third failed the minimum-size check throughout the replay. The
previous offline probe's “accepted” candidates therefore did not necessarily become valid
splits when production topology checks were applied.

The rendered split image shows land, water, cratons, orogenic belts, pseudo-sutures, the
daughter plate edge in orange, and the graph-cut trace in red/white. The image and full
per-step replay log remain outside the repository; the compact results and image checksum
are in `applied-split-replay.json`.
