# G3b independent calibration review, 2026-11-27

This is calibration on an already-exposed development world, not fresh evaluation or an
owner-approved protocol. Two independent Codex reviewers examined seed 910211954 at
221.6 Myr, density 4, `ridge_025`, in both finite-volume and wind-arrow solvers and
both diagnostic and CFD winds. Each received the draft and the same capture, without
seeing the other's review. Their observations below are provisional because the available
packet does not meet the draft's required display and supporting evidence.

| Pattern | Reviewer 1, all four combinations | Reviewer 2, all four combinations | Resolution for calibration |
| --- | --- | --- | --- |
| Rain belts | Provisionally recognizable: broad equatorial and southern wet regions span several rows; flag a narrow northern stripe. | Inconclusive: maps suggest broad bands, but latitude profiles by surface and hemisphere and circulation context are missing. | Inconclusive. Visible bands are suggestive, but the stipulated profile and circulation checks are needed to distinguish a broad belt from artifacts. |
| Rain shadows | Provisionally weak or absent: three-cell signed contrasts are negative and fewer than half of profiles are windward-wetter. | Provisionally weak or absent on the same evidence; maps lack geographic overlays and per-range evidence. | Provisional weak/absent, pending a compliant packet. A crest peak alone does not establish a multi-cell shadow. |
| Dry interiors | Inconclusive: no matched latitude/exposure evidence or labeled geography. | Inconclusive for the same reason. | Inconclusive; aggregate coast-distance means alone cannot establish a comparable interior. |

The 752–849 overlapping eligible profiles have median windward-minus-lee contrasts of
−0.08 to −8.54 mm/yr and 36.6–45.8% windward-wetter profiles across the four combinations.
They are not independent mountain ranges. These measurements describe the calibration
world and are not new automatic thresholds. See the draft and
[g3b-rain-shadow-transect-calibration-seed910211954-d4.json](g3b-rain-shadow-transect-calibration-seed910211954-d4.json).

Both reviewers identified a material packet mismatch: the four existing calibration
PNGs use logarithmic scales and different upper limits (approximately 6,863–7,723 mm/yr),
whereas the draft calls for one shared linear scale. They lack co-registered terrain and
wind context; the provided packet also lacks the prescribed latitude profiles and
latitude/exposure-stratified interior comparisons. Therefore these ratings do not validate
the proposed review procedure. Rebuild a compliant calibration packet and repeat the
independent review before freezing it. Do not represent the numerical summaries alone as
independent map-based adjudication.

## Proposed clarifications, pending calibration and owner agreement

- A rating is inconclusive if essential geography, wind direction, support, or corroborating
  profile evidence needed to distinguish a pattern from an artifact is missing or unreliable.
  Weak/absent requires enough applicable evidence to assess the pattern. A coherent visible
  feature is not automatically disqualified by a missed historical G3 numerical cutoff.
- Reviewers record independent initial ratings and evidence before seeing each other's
  judgments. They discuss disagreements against frozen rules; record both original ratings
  and the rationale. If disagreement cannot be resolved on the existing packet, the final
  rating is inconclusive, not a forced majority or a new threshold.
- Show maps on one fixed projection and one numerical linear precipitation scale, with
  over-range values marked, plus co-registered geography and winds. Include area-weighted
  latitude profiles, geographic eligibility and support, multi-cell transects with range
  context where feasible, and coast-distance comparisons stratified by latitude and
  exposure. Define reproducible masks, grouping, and stratification before evaluation.
- Freeze the candidate and exact code, solver/wind modes, density, generation and snapshot
  procedure, display scale and units, masks, support rules, reviewer identities, seed
  selection and replacement rule, and artifact adjudication before opening fresh candidate
  rainfall. At least five evaluable independent seeds per pattern and mode are proposed;
  fresh seed IDs and inputs have not yet been selected or approved.

G3 remains formally failed. Cold-land precipitation remains descriptive, and water-ledger
closure and nonnegative water remain separate invariants. No fresh evaluation has begun.
