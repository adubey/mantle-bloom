# Issue 1048: rainfall-pattern calibration packet

The packet covers the previously exposed development seed `910211954`, 221.6 Myr,
density 4, with the `ridge_025` arm, two transport solvers and both diagnostic and CFD
winds. It is calibration evidence, not a fresh evaluation or an amendment to G3's failed
result. Cold-land rainfall remains descriptive. The [independent preliminary
reviews](../rfc0017/g3b-independent-calibration-review.md) identify the missing evidence
that this packet is intended to address; a repeat independent review has not yet been done.

`bin/debug/rfc0019_g3b_calibration_packet.py CAPTURE OUT` recomputes the candidate and
control from the captured fields without loading or stepping a world. The generator is
restricted to the already-exposed pinned input save `seed910211954-221.6myr.bloom`,
SHA-256 `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67`.
The packet was generated at code commit `420705ad` using verified capture inputs; its
`packet.json` records individual input and source checksums. The capture and full packet
are stored outside Git; see the private `rfc-0019-packet-1048.md` note for recovery.
The packet JSON SHA-256 is `fe61fb5747d969d01b891198bd6be2c4c1f16f2bc7c4337e8cd51e76e5f9adb8`;
the file manifest SHA-256 is `fd3c7bfb7fc28496a4eb1001c8c76cd70f3d96790e0e7ab6cd471b4fa8f571c5`.

All reference and candidate maps share a linear 0–3469.8875 mm/yr scale, computed once
as the unweighted pooled 99th percentile of positive current-climate native cells across
the two winds. Magenta marks over-range values; unsaturated diagnostic views use separate
scales and must not be compared by color. The packet includes terrain and wind context,
area-weighted 5-degree latitude profiles, 15-degree latitude and upwind-ocean exposure
stratified coast-distance tables with area support, and geolocated multi-cell ridge
profiles. Upwind-ocean exposure is sampled every 50 km to 500 km along the local wind;
it is a geographic proxy, not a complete description of orographic shelter. Coastal
strata pool landmasses and do not prove that an interior and coast are on the same
continent. Ridge points overlap and are not independent ranges. Reviewers must consider
these limitations before assigning labels; the packet is not a frozen evaluation protocol.

The recomputed multi-cell ridge summaries match the previously stored calibration JSON
exactly in all four cases. Candidate hard gates pass in each case, with maximum global
relative water-ledger residual `2.64e-16`. The earlier review record has been corrected:
the crossing counts were mistakenly called eligible profile counts.

| Wind | Solver | Eligible multi-cell profiles | Median windward minus lee (mm/yr) | Native cells over fixed scale |
| --- | --- | ---: | ---: | ---: |
| Diagnostic | Finite volume | 812 | -0.387 | 13,384 |
| Diagnostic | Wind arrows | 729 | -0.082 | 16,671 |
| CFD | Finite volume | 832 | -5.546 | 13,849 |
| CFD | Wind arrows | 817 | -8.539 | 16,579 |

Next, independent reviewers should apply the draft rubric to this packet and resolve
whether its maps and stratification suffice. Before any fresh evaluation, settle the
landmass/range comparability rules and freeze the candidate, supported modes, masks,
seed selection and reviewer procedure for owner approval.
