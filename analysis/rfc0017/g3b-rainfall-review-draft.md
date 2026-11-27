# RFC 0019 G3b: draft rainfall-pattern review

This working draft translates the owner's 2026-11-25 direction into a repeatable
review. It does not replace the registered G3 protocol or retroactively change its
failed outcome. Numerical map scales, geography eligibility and reviewer assignment
must be fixed before inspecting the new evaluation worlds.

## Purpose and scope

Assess whether the chosen atmosphere produces broadly recognizable rainfall patterns
across independent synthetic worlds. G3a addresses downstream landscape diversity;
water accounting, runtime and other engineering requirements remain separate checks.
Cold-land precipitation is a reported diagnostic, not a G3b pass requirement.

Use at least five independent seeds. Freeze the candidate, code, parameters, input
hashes, wind modes, density, geographic masks, map scales and review rules before
opening evaluation rainfall. Previously exposed G3 held-outs and G4 worlds may inform
development and interpretation but are not fresh evaluation evidence. Do not choose
or discard a seed based on its candidate rainfall. Record missing or incompatible
inputs and the reason before selecting a replacement.

## Evidence and per-world assessment

For each seed and supported wind mode, produce a candidate precipitation map and
terrain/ocean map on a common projection, with a common fixed rainfall color scale
across all worlds and arms. Supply latitude profiles and the geographic measurements
below; also show the legacy map as context, not a target. Label units, seed, save age,
wind mode, density and code version. Predefine land and ocean masks, coast-distance
bands, eligible ridge crossings and area support using geography and winds alone.
Report support area and crossing count so that missing geography is visible.

For development maps, set the common linear scale once from the pooled 99th percentile
of positive precipitation in the matched current-climate maps. Use that same scale for
every candidate and world; mark values above it as over-range and include an unsaturated
version in the diagnostic bundle. Do not normalize each map separately, since that can
hide large differences in amount. Use the established 0.1% planetary-area support rule
for region pairs and the 100-crossing / 25-distinct-crest support rule for rain shadows.
These are geography screens, not rainfall thresholds. Choose and record the map projection
and orientation once for the whole evaluation set.

For each applicable pattern, record **recognizable**, **weak or absent**, or
**inconclusive**, with a short explanation referencing map features and measurements:

| Pattern | Evidence to inspect | Recognizable signal |
| --- | --- | --- |
| Rain belts | Area-weighted latitude profiles by hemisphere and surface; belt width and neighboring rainfall; map. | One or more broad, coherent wet regions with a visible contrast to adjacent areas, consistent with the circulation. Profiles distinguish broad belts from one-row spikes. |
| Rain shadows | Map and distribution of matched windward/lee rainfall at eligible mountain crossings, with sample count and wind direction. | Exposed and sheltered sides show a repeated geographic contrast across multiple crossings. Report the signed median contrast and fraction of crossings in the expected direction. |
| Dry interiors | Rain by coast distance, latitude and exposure; map and land area in each band. | A coherent inland region on an eligible large landmass is drier than nearby, comparable coastal land. Stratification shows whether latitude or exposure explains the difference. |

Do not require the original fixed 1.5 belt ratio, 60% ridge-crossing threshold, or
0.6 inland/coast ratio for this draft review. Report those numbers as descriptive
comparisons when available. In particular, distinguish the absolute windward/lee
pattern from the earlier metric that measured change relative to a matched control.
If a world lacks adequate geography or the measurement cannot be trusted, mark that
pattern inconclusive; never count it as a pass or failure. Record cold-coast and
cold-inland rainfall, frozen precipitation and their geographic coverage separately.

Two reviewers independently classify each map set as recognizable, weak/absent, or
inconclusive and cite the map and measurements supporting the rating. Resolve different
ratings against the definitions above before pooling worlds; preserve both initial ratings
and the reason for the resolution. Pool only independent seeds. A pattern is recognizable
for this sample when a strict majority of its evaluable worlds receives the final
recognizable rating; require at least five evaluable worlds for a sample decision. Report
counts and individual worlds, not only the majority label.

## Across-world decision

For each pattern and wind mode, count independent seeds assessed as recognizable and
weak/absent. A recognizable majority means more than half of the evaluable worlds;
at least five evaluable independent worlds are needed for that pattern's decision.
Fewer than five is inconclusive and calls for additional suitable worlds, not a pass.
Multiple snapshots, grid densities and wind modes from one seed are not independent
worlds. Report results by mode, with maps and per-world judgments available for review.

Separately flag repeated artifacts, including narrow wet stripes, rainfall concentrated
in a small number of cells, or widespread reversal of terrain-related rain patterns.
Record their location, extent and number of affected worlds. A flagged artifact calls
for an explicit review decision; this draft does not assert a numerical failure
threshold. Report whether an artifact recurs in a majority of evaluable worlds, alongside
its per-world evidence. An artifact flag prompts review of its landscape effects; it does
not automatically fail G3b. Freeze reviewers and resolution procedure before the fresh
evaluation; do not change a definition after seeing evaluation maps.

## Calibration examples for borderline evidence

Use these examples to train reviewers on the labels, not as new candidate rankings or
acceptance evidence. Both density-4 held-out worlds have already been exposed. The older
G3 screen reported five of sixteen inland/coast combinations just above its `0.6` limit
(range `0.613–0.835`). Under this draft, a ratio of `0.613` is not an automatic failure:
reviewers should ask whether a coherent inland region is visibly drier than comparable
coastal land and whether latitude/exposure stratification supports that reading. If the
map suggests a dry interior but the available aggregate ratio cannot separate it from
latitude or exposure, label it **inconclusive** until the stratified diagnostic is
available. If the map and stratified values agree, label the pattern **recognizable**
even though the historical cutoff was narrowly missed. This does not change the recorded
G3 failure under its original preregistration.

Use the `ridge_025` candidate record for calibration, not the separate saturation-bounded
column probe in `results.json`; they are different arms. In the candidate record, the
diagnostic finite-volume distance means for `seed507012855` are about `696`, `252`, and
`225 mm/yr` across the 0–500, 500–1,500, and over-1,500 km bands. For `seed762517225` they
are about `769`, `470`, and `251 mm/yr`. Under CFD winds they are about `417`, `325`, and
`307 mm/yr` for the first seed, and `529`, `538`, and `274 mm/yr` for the second. These
are descriptive bands, not matched-latitude tests: the first two diagnostic examples
support inland drying, while the CFD examples need maps and latitude/exposure-stratified
values because the decline is weak or nonmonotone. None of these existing exposed worlds
counts toward a new evaluation.

The two apparently conflicting rain-shadow summaries came from different arms and
different masks, so they should not be read as one result. In the saturation-bounded
column probe, `windward_leeward_area_weighted_pair` divides area-weighted rain means over
separate extreme-ascent and extreme-descent land masks (the top and bottom deciles of
terrain-ascent projection). `pair_fraction_windward_wetter` instead counts only one-cell
downwind neighbors that also fall in the opposite extreme mask, then reports the
unweighted fraction whose windward rain is greater. The latter is a narrow adjacency
sample, not a census of ranges; the former can be high because of broad regional
differences unrelated to a particular ridge. Neither is a stand-alone rain-shadow verdict.

The frozen `ridge_025` output has a more mechanism-specific `ridge_transects` diagnostic:
it follows qualifying transport faces over a ≥250 m rise to the next downwind lower land
cell. On `seed507012855`, the finite-volume diagnostic-wind arm reports 1,479 crossings
over 324 distinct crest **cells**; median absolute crest-minus-lee rain is `12.1 mm/yr`,
median contrast gain over the paired control is `30.1 mm/yr`, and 90.1% of crossings have
crest rain increase with lee rain decrease. This is a strong local signal. On
`seed762517225` with the same mode and solver, 1,859 crossings over 435 crest cells have
only `0.003 mm/yr` median absolute contrast, despite a `2.2 mm/yr` median contrast gain
and an 88.9% share with the expected paired change. That is a borderline example: the
candidate adds the expected contrast at many crossings, but its absolute median amounts
are nearly equal and low. Also, “distinct crests” here means distinct grid cells, not
independent mountain ranges, and the gain is relative to control rather than absolute
rainfall climatology. For G3b, use these transects as diagnostics; rate recognizability
from the candidate's absolute map and fixed, geographically supported multi-cell
windward-to-lee transects, with per-range evidence where feasible.

The existing ridge-rainout diagnostic now emits `ridge_native_cell_profiles` for each arm.
For each eligible coarse transport crossing, it centers a transect on the highest native
land cell in that crossing's coarse crest cell and samples offsets from three native cells
upwind through three downwind. It deduplicates profiles by coarse crest cell and flow
direction. The compact summary reports per-offset sample counts and candidate/control
rain quantiles, plus per-profile averages for the three-cell windward and lee sides when
at least two land cells are available on each side. `--include-profile-samples` retains
every individual transect; the default output omits those rows. The profiles are fixed
across paired control and ridge arms because the crossing geometry is shared; they are
not independent mountain ranges. The exposed capture files contain saved terrain, winds
and source inputs but not the selected candidate's native rain array, so this diagnostic
recomputes the offline candidate from those captured fields. It does not load or step a
world. These already-opened worlds remain calibration only and do not count as fresh
evaluation.

The recomputations on both exposed density-4 worlds show a clear crest-local peak but
little median windward-minus-lee contrast when averaged across three native cells on each
side. Candidate median rainfall at one cell upwind / crest / one cell downwind was `79.7 /
98.1 / 84.0 mm/yr` for seed 507012855 finite-volume diagnostic wind, `125.9 / 164.4 / 139.8`
for finite-volume CFD, `38.3 / 51.6 / 47.5` for wind arrows diagnostic, and `70.6 / 92.2 /
82.2` for wind arrows CFD. For seed 762517225, the corresponding triplets are `14.7 / 26.0 /
26.1`, `23.6 / 41.5 / 30.1`, `25.1 / 30.8 / 23.6`, and `13.5 / 23.3 / 19.2 mm/yr`.
The wider paired profiles are:

| Seed | Solver and wind | Eligible profiles | Median 3-cell windward-minus-lee (mm/yr) | Profiles with windward wetter | Median contrast gain vs control (mm/yr) |
|---|---|---:|---:|---:|---:|
| 507012855 | Finite volume, diagnostic | 383 | -0.08 | 45.2% | 0.01 |
| 507012855 | Finite volume, CFD | 377 | -6.41 | 41.6% | 0.18 |
| 507012855 | Wind arrows, diagnostic | 330 | -0.07 | 48.8% | 0.01 |
| 507012855 | Wind arrows, CFD | 368 | -3.67 | 42.4% | 0.09 |
| 762517225 | Finite volume, diagnostic | 473 | -0.0003 | 46.1% | 0.00008 |
| 762517225 | Finite volume, CFD | 481 | -0.0036 | 39.1% | 0.00028 |
| 762517225 | Wind arrows, diagnostic | 425 | -0.011 | 44.7% | 0.00001 |
| 762517225 | Wind arrows, CFD | 473 | -0.080 | 35.5% | 0.00029 |

The common-scale maps show localized ridge changes, not a broad displacement. Across these
two worlds the paired 3-cell contrast is near zero or negative in every solver/wind mode,
while one-cell crest rain is higher than its immediate upwind and lee neighbors. This
supports a crest-local enhancement, but does not establish a robust multi-cell dry-lee
pattern. The high treatment-relative `crest_wetter_and_lee_drier_share` in the older
diagnostic therefore must not be treated as an absolute recognizable-shadow result.
These are two already-exposed calibration worlds, not fresh evaluation seeds; profiles
share crest cells and are not independent mountain ranges. The per-world compact outputs
are [seed 507012855](g3b-rain-shadow-transect-calibration-seed507012855-d4.json) and
[seed 762517225](g3b-rain-shadow-transect-calibration-seed762517225-d4.json).

A third calibration used the registered development save `seed910211954-221.6myr.bloom`
(SHA-256 `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67`) at density 4.
It is explicitly marked `heldout: false`; it is method calibration only and does not count
toward evaluation. All four solver/wind combinations have adequate crossing support. For
the frozen `ridge_025` arm, the median of the per-profile three-cell windward-minus-lee
contrast remains negative in every combination, even though about half the profiles gain
contrast relative to control:

| Solver and wind | Eligible profiles | Median candidate windward-minus-lee (mm/yr) | Profiles where windward-side mean exceeds lee-side mean | Profiles with positive contrast gain vs control |
|---|---:|---:|---:|---:|
| Finite volume, diagnostic | 812 | -0.39 | 37.7% | 50.9% |
| Finite volume, CFD | 832 | -5.55 | 37.5% | 52.4% |
| Wind arrows, diagnostic | 729 | -0.08 | 45.8% | 50.3% |
| Wind arrows, CFD | 817 | -8.54 | 36.6% | 51.8% |

The pooled median rainfall by native-cell offset for `ridge_025` peaks at the selected crest
(offset 0) in each combination: `35.3 mm/yr` for finite-volume diagnostic winds,
`170.5 mm/yr` for finite-volume CFD, `20.9 mm/yr` for wind-arrows diagnostic winds, and
`140.1 mm/yr` for wind-arrows CFD. In the CFD profiles, the downwind side stays wetter than
the upwind side in the aggregate offset medians despite the crest-local peak. This is a
useful calibration example: a crest peak and roughly 50% treatment-relative contrast gain
do not imply a repeated absolute windward-wet / lee-dry pattern. The common-scale maps also
show that changes are localized around ridges while broad rain-belt structure is retained.
For each profile, the windward and lee values are arithmetic means over the valid native
cells at offsets `-3..-1` and `+1..+3`, respectively, with at least two valid land cells on
each side. Thus the table's percentage is the share of eligible profiles whose windward-side
mean exceeds its lee-side mean; it is neither a percent rainfall increase nor a fraction of
mountain land area. Profiles overlap at coarse crest cells and are not independent ranges,
so the profile fractions describe this world rather than provide seed-level uncertainty. See the compact
[development calibration output](g3b-rain-shadow-transect-calibration-seed910211954-d4.json);
the full captured inputs, logs, maps and `SHA256SUMS` remain in the private capture store.

**Provisional calibration assessment:** the `ridge_025` rain-shadow signal is weak or absent
in the three examined calibration worlds. The multi-cell absolute contrast is near zero or
negative across the reported solver/wind combinations, and windward-side means exceed lee
means in fewer than half of eligible profiles. The repeated crest-local rainfall peak does
not provide the missing broad windward-to-lee contrast. Treat this as evidence that the
current ridge-rainout treatment is not yet a successful rain-shadow candidate; it is not a
formal G3b evaluation because these worlds were used for calibration and have not received
the frozen independent-review procedure.

A post-hoc flux trace on the development capture separates the mechanism from the rainfall
pattern. It re-solves the same stationary control and `ridge_025` budgets without stepping a
world or changing parameters, groups eligible uphill faces by unique crest-to-lee reservoir
link, and measures incoming advection, intercepted crest rain, crest-to-lee net transport,
lee column water and lee rainfall. The paired changes below compare the same links. Rates
are medians per coarse link; percentages are the fraction of links where the candidate value
is lower than control. The rainfall column uses actual placed precipitation after overflow
redistribution, averaged over the receiving lee reservoir.

| Solver and wind | Links | Incoming advected vapor, control → candidate (kg/s) | Intercepted crest rain (kg/s) | Crest-to-lee transport Δ (kg/s; lower links) | Lee column water Δ (kg/m²; lower links) | Lee rain Δ (mm/yr; lower links) |
|---|---:|---:|---:|---:|---:|---:|
| Finite volume, diagnostic | 828 | 78,647 → 71,302 | 9,674 | -4,151; 90.0% | -0.0102; 91.1% | -1.35; 93.8% |
| Finite volume, CFD | 849 | 608,632 → 508,866 | 66,782 | -105,211; 99.5% | -0.1275; 98.9% | -6.25; 93.4% |
| Wind arrows, diagnostic | 752 | 35,762 → 33,970 | 4,530 | -2,816; 93.9% | -0.0083; 94.7% | -0.52; 94.0% |
| Wind arrows, CFD | 842 | 521,993 → 434,121 | 54,961 | -69,876; 99.2% | -0.0653; 98.2% | -2.85; 93.0% |

The treatment is therefore active: it captures water at the ridge, and most paired links
deliver less moisture and rain to the lee reservoir. The median lee-rain reduction is still
only `0.52–6.25 mm/yr` after saturation overflow is spatially placed, and it does not
overcome the broader absolute rainfall gradients in the native-cell transects. The
component decomposition shows that the paired median bulk-sink rain falls by `0.33–5.12
mm/yr`; most lee reservoirs have no change in intercepted rain (median change 0), while
overflow placement has a zero median change and a negative mean change because the effect is
concentrated in fewer links. The combined local land-source proxy (lake/sea evaporation,
river evaporation, and vegetation) is fixed between arms; its median equivalent rate is
`46–47 mm/yr` for CFD lee reservoirs and `2.8–3.7 mm/yr` for diagnostic-wind lee reservoirs.
That local source could replenish downwind columns in CFD mode, but its equivalent rate is
not rainfall and does not prove the cause of the weak absolute shadow. Global relative
closure is at most `2.7e-16` and summed local closure at most `2.7e-15`, so the observed
limitation is pattern strength rather than a water-ledger leak. This is descriptive evidence
from one calibration world; links and profiles overlap and are not independent ranges. See
the compact [flux trace](g3b-ridge-flux-transect-seed910211954-d4.json) and [placed-rain
decomposition](g3b-ridge-rain-budget-decomposition-seed910211954-d4.json).

For rain belts, a high latitude-band mean by itself is not enough. A broad band should
occupy a visible span of adjacent latitude rows/bins and have a neighboring contrast on
the same fixed scale. A single bright row, especially if it repeats at a grid boundary,
is an artifact candidate and does not establish a belt. When map and area-weighted profile
disagree, first check the mask, area weights, clipping and whether the apparent band is
only one or two cells wide; resolve to **inconclusive** if those checks do not settle it.
These examples motivate the qualitative labels and diagnostics above; they do not create
new numeric thresholds.

## Existing evidence and next work

The earlier frozen `p=0.25` ridge candidate passed the original belt, ridge and
exposed-cold-coast checks on two held-out worlds, while five of sixteen combinations
missed the original inland/coast limit. See the [held-out result](g3-coarse-transport-plan.md#held-out-pattern-screen-for-frozen-ridge_025-candidate).
The subsequent G4 candidate also used a coast-focused sink; its results do not establish
the pattern performance of the ridge-only version or vice versa. Existing worlds can
help set map scales, eligibility and the review rubric, but cannot satisfy the fresh
evaluation count.

Before evaluation, apply this rubric to development worlds and resolve any ambiguous
examples, then freeze the candidate, supported wind modes, reviewers, masks, projection,
map scale, and at least five fresh independent seed IDs and save hashes. If fewer than
five selected worlds support a pattern, add suitable worlds chosen by geography alone.
Compare the finalized protocol with the original G3 criteria in an RFC revision before
using it as an acceptance gate; this draft is not a new G3b result.
