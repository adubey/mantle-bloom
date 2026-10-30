# Gate 3 acceptance matrix

**Disposition: Gate 3a (local-scale evidence) is closed by the owner on 2026-10-30. Gate 3b
(plate-scale geological breakup) remains open and is tracked by [#699](https://github.com/adubey/mantle-bloom/issues/699).**
Closing 3a completes the bounded local investigation and records its evidence; it does not claim the
coupons are a validated predictive rift model. The millimeter-to-kilometer results are supporting
physics evidence for 3b, whose pass requires controlled histories on real spherical plate meshes at
production resolution over geological time. Its run matrix and tolerances must be preregistered.

## What worked locally, and what did not

| Piece tested | What worked | What did not, or remains unproved |
|---|---|---|
| **Weak-crack initiation** | At about `7 MPa` host tension, the `3 MPa` weak-flaw scenario initiated in the screening cases; the `14 MPa` intact-basalt proxy stayed intact in `0/9` flaw cases. | Flaws were supplied to the test and the 3 MPa strength is a sensitivity choice, not a measured natural value. This supports a weak-versus-intact response, not spontaneous rift initiation. |
| **Local crack-tip growth** | A refined cohesive-tip coupon gave the same `80 mm` advance on `2.5 mm` and `5 mm` meshes; cohesive work differed by `0.014%`. A stronger boundary also arrested a local tip in its controlled coupon. | The tip was driven by an imposed crack-tip load. It does not show that a plate-scale rift starts, finds its path, or grows under a naturally evolving tectonic load. |
| **Finite magma-source stress transfer** | On both `32×20` and `48×40` meshes, the dike relieved opening-direction stress behind it and raised tension ahead by similar amounts (`0.41–0.45 MPa`). The modeled source-volume relation closed closely. | The crack was fixed and pre-seeded: there was no new damage ahead. The pressure-work comparison still differed by `7.5–7.9%` and omitted tectonic boundary work and dissipation. So this verified the direction of stress transfer, not propagation. |
| **Combined host, flaw, magma and mantle-replacement screen** | The pieces could be connected in a single reduced scenario, with volume replacement and source accounting reported. | The 65–66 km endpoint was prescribed. The weak case used an uncalibrated 3 MPa strength, and the simplified flow model did not feed tip changes back into host stress. It was an interaction check, not a predicted rift. |
| **Longer extension and lateral mantle-load screen** | Extending the run to `30–100 Myr` produced more failed cells; adding an idealized outward mantle traction increased their count. | Failure stayed fragmented, with many disconnected components rather than one continuous, advancing rift. Counts changed with mesh, and the imposed traction was not a solved convection field. |
| **Kilometer-scale cohesive seam** | A local fine-mesh test could resolve crack-tip response. | At `1 km` spacing, the estimated `14 MPa` process zone was only `0.0294 m`; all `136` interface points separated together. The full-section mesh could not resolve gradual tip growth. |
| **Host-to-tip transfer over long times** | Short local plumbing and conservation checks were useful for debugging. | The affine transfer error was vastly larger than the local displacement; non-affine transfer still had millimeter errors on a sub-meter patch, and longer host solves failed numerically. This route did not provide a reliable geological-timescale coupling. |

In short, Gate 3a established several local pieces of mechanics and showed their limits. It did **not** produce a validated, plate-scale rift: spontaneous initiation, sustained growth and arrest in an evolving host, and successful plate splitting remain Gate 3b questions.

| RFC requirement | Best existing evidence and measured result | Status against criterion |
|---|---|---|
| **Gate 3a: local-scale mechanism evidence** | Constitutive and conservation coupons; weak-seed versus intact-control initiation (3 MPa sensitivity initiates; 14 MPa control did not in 9 tested flaws); local finite-strain, layered brittle–ductile and moving-tip growth/arrest diagnostics; matched finite-source stress transfer. See the committed [magma-intrusion review](gate3-magma-intrusion-review.md) and the RFC 0010 evidence history. | **Closed by owner as a bounded evidence phase.** Limitations remain explicit: outcomes are not calibrated natural-rift parameters; the recent finite-source test uses fixed seeded crack geometry and has no outside-seed damage; its pressure-work comparison omits boundary and dissipative terms. This is not a predictive plate-rifting validation. |
| **Gate 3b: plate-scale geological breakup** — actual spherical mesh; weak/strong lithosphere; connected boundary before split; successful and failed rifts; daughter plates; refinement and budgets | Existing local evidence includes a 1.5× thinned weak corridor that takes `1.87×` the domain-mean extension, plus a moving-tip coupon where a 10× tougher band arrests advance. These use 5–10 km and 2.5–5 mm meshes, respectively, and do not test a complete production plate or geological history. | **Open; tracked by [#699](https://github.com/adubey/mantle-bloom/issues/699).** First preregister a matched whole-sphere run matrix and tolerances. |
| **Weak versus intact initiation fork** | At the `7 MPa` coupon stress, the `3 MPa` weak-seed sensitivity initiates while the `14 MPa` intact-basalt proxy initiates in `0/9` tested flaw cases. The strengths are scenarios, not calibrated values; this is a controlled sensitivity result, not a natural-rift calibration. | **Evidence supports the fork.** Treat the intact result as a valid negative and the weak case as a sensitivity; do not tune strength further to force initiation. |
| **Gate 3b conditional mantle/magma attribution** — matched geological-duration whole-sphere histories, compare rift and breakup outcomes with/without represented plume or magma forcing; balance source and force/torque budgets | Local finite-source tests show dike-normal stress relief behind and increased tension ahead at 32×20 and 48×40, but the 100-year test has fixed crack geometry, no ahead-of-seed damage, and a 7.5–7.9% pressure-work difference before omitted energy terms. | **Local component evidence only.** No plate-scale causal attribution has been measured. Plumes/magma need not cause every modeled rift. |
| **Host-to-tip coupling method and tectonic timescale** | The affine transfer has a `0.143 m` fit residual against a `4–8 μm` local signal; the non-affine synthetic transfer has about `1 mm` boundary error on a `0.4 m` patch. Year-scale tests are plumbing checks, and longer host solves fail numerically. | **Retire affine transfer and stop duration/radius sweeps.** Neither failed transfer nor an unreachable Myr integration is a Gate 3 failure criterion. Nested cutout/submodelling is a named follow-up only if host-derived coupling is needed for a later validation. |
| **Criteria 4–7** — split frequency and count context, craton protection, ledgers, save compatibility and runtime | No scale-corrected whole-sphere replay has yet been run against these requirements. | **Open.** These must be measured as part of the plate-scale Gate 3 evaluation, not deferred until after it. |

## Gate 3b scope (opened 2026-10-30)

The millimeter-to-kilometer crack and host coupons remain useful for constitutive verification,
local energy accounting and mechanism comparison. They cannot close Gate 3 because their domain and
time span are far below the spherical production model. The 1–30 year coupling runs remain
software-plumbing checks, not tectonic outcomes. Do not continue local mesh/radius/duration sweeps
as Gate 3 closure work, and do not tune material strength to force a coupon outcome.

The next Gate 3b action is to preregister a whole-sphere geological replay matrix on actual plate
meshes, including a matched baseline, the proposed rift closure, representative weak corridors and
strong/cratonic regions, and separate mantle/plume/magma controls where those fields exist. It must
state the world seeds and initial states, represented duration, production and refined time steps,
node densities, event grouping, outcome metrics, tolerances, material/force/torque budgets and
runtime limit. The 6–8 breakup-systems-per-200-Myr Earth reconstruction remains a broad plausibility
reference, not a hard threshold, until the model's event grouping is fixed. No new scale-corrected
pass is claimed by the existing coupon results.
