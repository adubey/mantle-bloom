# Issue #547: post-topology plate rates

## Current scope

On 2026-10-25, the project owner narrowed #547 to the **merge rate assignment**. This issue
removes the angular-momentum blend and rate clamp from merges while retaining split daughter
pole assignments until RFC 0010 changes split selection. The earlier replay and traces below
tested removing both merge and split assignments; they are historical evidence for coupling,
not acceptance evidence for the merge-only implementation. The merge-only matched metric check
is now complete; see the results below. On 2026-10-29, RFC 0012 changed the timestep
criterion: report four-seed mean metrics for each arm and step, plus per-seed spread, topology
differences, and per-run ledger closure. The former 10% timestep tolerance is no longer a
pass/fail gate. The older combined-change replay remains historical. The implementation PR rebases the merge
change on current main `2cefff3b`; the four-seed matrix below was run earlier against baseline
`f04cbc07` and candidate `d553c3cf`. Current main has since changed its material-overflow tracer
handling. No new 24-run matrix was run for this PR, so those metrics characterize the recorded
matched pair rather than a fresh comparison against `2cefff3b`.

## Merge-only matched replay (recorded comparison)

**Question.** Does removing the merge angular-momentum blend and rate clamp change motion,
deformation, ledgers, or geological outcomes when split daughter pole assignments remain?

**Controls.** The matrix has 24 runs: baseline commit
`f155b7eb36727f452139955b13350593db57cbcd` and candidate commit
`45b035c34ed0c370e52d1386e17c14318e893d5a`, paired across seeds `326494039`, `748935664`,
`580014012`, `111065103` and steps 0.1, 1, and 2 Myr. Each run starts from the same generated
world, uses node density 4, `PYTHONHASHSEED=0`, and advances 50 Myr with 10 Myr checkpoints.
Replay outputs remain in scratch storage. The compact matched summary is
[merge_only_replay_summary.json](merge_only_replay_summary.json); the checked-in summarizer is
[`bin/debug/summarize_rfc0012_merge_replay.py`](../../bin/debug/summarize_rfc0012_merge_replay.py).

**Matched outcomes.** Mean solved rate and accumulated driving torque changed by at most 2.28%
and 2.43%, respectively; peak solved rate was unchanged in all 12 pairs. Contact response was
less stable: mean absolute closing rate changed by as much as +35.45%, mean contact force by
+32.97%, and peak absolute closing rate by −31.76%. The largest mean-closing change is
`580014012` at 2 Myr (0.0448 to 0.0607 cm/yr). The largest force change is `748935664` at 2 Myr
(6.106 to 8.119 × 10¹² N/m). These are real differences in contact metrics, not just percentage
changes around a near-zero value.

| Seed | Step (Myr) | Mean solved rate Δ | Mean closing Δ | Mean contact force Δ | Topology events (baseline/candidate) |
|---:|---:|---:|---:|---:|---:|
| 326494039 | 0.1 | 0.00% | 0.00% | 0.00% | 17 / 17, same sequence |
| 326494039 | 1 | +0.07% | +9.41% | +1.62% | 10 / 9, different |
| 326494039 | 2 | 0.00% | 0.00% | 0.00% | 7 / 7, same sequence |
| 748935664 | 0.1 | +0.48% | +0.60% | −0.11% | 36 / 32, different |
| 748935664 | 1 | +0.00% | −0.00% | +0.02% | 15 / 15, same sequence |
| 748935664 | 2 | +0.31% | −8.15% | +32.97% | 14 / 15, different |
| 580014012 | 0.1 | −2.28% | −0.55% | +2.40% | 11 / 11, same sequence |
| 580014012 | 1 | +0.01% | +0.11% | −0.03% | 11 / 11, same sequence |
| 580014012 | 2 | +1.74% | +35.45% | −10.95% | 10 / 10, same sequence |
| 111065103 | 0.1 | 0.00% | 0.00% | 0.00% | 9 / 9, same sequence |
| 111065103 | 1 | 0.00% | 0.00% | 0.00% | 7 / 7, same sequence |
| 111065103 | 2 | 0.00% | 0.00% | 0.00% | 8 / 8, same sequence |

Three pairs have different topology sequences. Their first mismatch is at 35.6 Myr for
`748935664` / 0.1 Myr (baseline accretes a stranded fragment; candidate sheds two nodes), at
40 Myr versus 47 Myr for `326494039` / 1 Myr (fragment accretion versus a later split), and at
22 Myr versus 24 Myr for `748935664` / 2 Myr (different plates split). Endpoint plate count
differs in one pair: 19 baseline versus 21 candidate for `748935664` / 2 Myr. The other nine
event sequences and eleven endpoint plate counts match.

Endpoint land fraction moves by at most 0.227 percentage points (`580014012` / 0.1 Myr).
Area-weighted Hc/Hm quantiles move by at most 0.521% (Hc p99, `748935664` / 0.1 Myr); Hm p99
moves at most 0.345%. The largest Hm boundary-retreat budget difference is −9.77 × 10¹⁶ m³
(−2.82%) for `748935664` / 2 Myr. Thus global thickness and land metrics remain close even in a
few pairs whose contact response and event history diverge.

The largest absolute Hm ledger closure residual is 1.022 × 10¹⁵ m³, or 0.00375% of inventory;
continental-node Hm closure is 1.017 × 10¹⁵ m³ (0.00810%) and oceanic-node closure is 1.005 ×
10¹⁴ m³ (0.00107%). The largest continental material-balance residual is 1.241 × 10¹⁴ m³;
craton and mobile-cover residuals are 512 and 9,472 m³. These maxima span both arms and all
checkpoints. The Hm residuals are small relative to inventory but are not at floating-point
roundoff scale in every coarse-step run.

Runtime changes range from −7.7% to +11.2% (mean +0.4%); peak RSS ranges from 1.18 to 1.90
GiB. These do not suggest a consistent performance change.

**Timestep characterization under the current acceptance criterion.** For each arm, the
following are the relative differences between the four-seed arithmetic mean metric at 1 or 2
Myr and the corresponding mean at 0.1 Myr. Per-seed values and spread remain in the JSON summary.
These values characterize sensitivity; no 10% pass/fail limit applies.

| Continuous metric | Arm | 1 vs 0.1 Myr | 2 vs 0.1 Myr |
|---|---|---:|---:|
| Mean solved rate | Baseline | +9.2% | +8.2% |
| Mean solved rate | Candidate | +9.7% | +9.7% |
| Mean absolute closing | Baseline | +35.8% | +118.8% |
| Mean absolute closing | Candidate | +39.3% | +103.8% |
| Mean contact force | Baseline | −14.4% | −22.5% |
| Mean contact force | Candidate | −14.1% | −20.7% |
| Peak absolute closing | Baseline | −7.6% | −9.3% |
| Peak absolute closing | Candidate | −4.8% | −13.6% |
| Accumulated driving torque | Baseline | +2.9% | +5.6% |
| Accumulated driving torque | Candidate | +3.0% | +5.3% |

The comparison does not show a consistent timestep-invariance improvement: the candidate is
closer on some contact metrics at 2 Myr, while worse on mean closing at 1 Myr and peak closing
at 2 Myr. Contact response remains strongly step-sensitive in both arms. The full comparison
and topology differences are recorded above and in the compact JSON.

**Assessment.** The merge-only change preserves many global metrics but changes contact response
and topology in some matched pairs. The focused merge trace confirms that a merge-rate change can
alter the next solve's classification. Under the current RFC contract, timestep differences are
reported as a limitation rather than used to fail acceptance. These results do not claim
unchanged trajectories or show that the candidate is consistently more timestep-invariant than
baseline. The old 10% threshold below records the preregistration used for the runs; RFC 0012's
updated acceptance wording supersedes it. No tuning was applied to these outcomes.

> **Historical preregistration note (2026-10-29):** the 10% timestep thresholds in the
> original preregistration and dated sections below were the criteria in force when those
> runs were designed. RFC 0012 now treats timestep sensitivity as characterization, not a
> pass/fail invariance gate. See the current four-seed averages above and [#547](https://github.com/adubey/mantle-bloom/issues/547).

### Updated matched replay preregistration on current main (2026-10-25)

The preceding 24-run matrix used baseline `f155b7eb36727f452139955b13350593db57cbcd`. Current
`origin/main` is `f04cbc0772eca6d2484f7477b20004deb368206c` and adds plume evolution to the world
step. To keep the #547 comparison current, repeat the same matrix against this new baseline and
the merge-only candidate built directly on it. This is a new matched comparison; results from
the earlier matrix remain historical and will not be combined with it.

**Question and prediction.** Does retaining the survivor's last solved rate after a merge,
without angular-momentum blending or a rate floor, change motion, deformation, ledgers, or
geological outcomes? The direction and size are not assumed. The focused trace predicts that a
merge can change the next solve's classification through the lagged rate, while the physical
solve then derives the updated rate from post-merge loads.

**Controls.** Compare baseline `f04cbc0772eca6d2484f7477b20004deb368206c` with the candidate
commit recorded below, on generated worlds for seeds `326494039`, `748935664`, `580014012`, and
`111065103`, node density 4, and 50 Myr horizons at 0.1, 1, and 2 Myr. Set
`PYTHONHASHSEED=0`, checkpoint every 10 Myr, and run one simulation at a time with no unit or
stress test running concurrently. Record commit, runtime, peak RSS, checkpoints, all metrics,
topology sequence and event times, and ledger closure. Full outputs stay in scratch; check in a
compact JSON summary and the command template.

**Measures and decision rule.** Compare accumulated driving and residual torque, mean and peak
solved rate, relative contact closing/slip and force, per-step Hc/Hm deformation, Hc/Hm
distributions, land fraction, topology event time/order, and continental/craton/mobile-cover/Hm
ledger closure. Apply RFC 0012's 10% continuous-metric timestep criterion separately to each
arm; report topology differences separately. Any >10% baseline/candidate continuous difference
must be traced to its first divergent step before acceptance or tuning. The result is diagnostic
if either arm misses the timestep tolerance. Ledger closure must be reported, including any
residual larger than rounding scale.

**Reproduction.** Reuse `bin/debug/replay_rfc0012_merge_rate.py`, selecting the baseline with
`MANTLE_BLOOM_REPLAY_BACKEND`; run each seed/step/arm serially and write compact output to a
scratch path. Compare the matched summaries with
`bin/debug/summarize_rfc0012_merge_replay.py`. No metric tuning will be made during this rerun.

### Current-main matched replay results (2026-10-25)

The preregistered current-main matrix completed all 24 runs. Baseline was
`f04cbc0772eca6d2484f7477b20004deb368206c`; candidate was
`d553c3cfc05db4372d0bb33ec7e58aa2390b793b`. All other controls matched the preregistration.
The compact results are in
[merge_only_replay_current_main_summary.json](merge_only_replay_current_main_summary.json).

**Matched outcomes.** Six of 12 event sequences differ; endpoint plate count differs in one
pair. The largest absolute candidate/baseline changes in solved-rate metrics were 4.07% for
mean rate and 3.08% for peak rate. Driving torque changed by at most 0.86%; residual torque by
9.58%. Contact response has a larger outlier: for seed `748935664` at 2 Myr, mean absolute
closing changed by −25.91% and peak closing by −32.95%. Mean contact force changed by at most
7.61%. Endpoint land fraction changed by at most 0.208 percentage points, and Hc/Hm quantiles
by at most 1.625%. The contact-closing differences exceed the preregistered 10% trace threshold. A follow-up trace
of that pair is recorded below; the overall matrix remains diagnostic because the separate
10% timestep tolerance is not met in either arm.

| Seed | Step (Myr) | Mean solved rate Δ | Mean closing Δ | Mean contact force Δ | Topology same? | Events (baseline/candidate) |
|---:|---:|---:|---:|---:|:---:|---:|
| 326494039 | 0.1 | −0.24% | −6.10% | −0.03% | No | 13 / 11 |
| 326494039 | 1 | +1.00% | +0.73% | −0.35% | No | 10 / 11 |
| 326494039 | 2 | −0.89% | −9.03% | +7.61% | No | 7 / 7 |
| 748935664 | 0.1 | +0.10% | +0.51% | +0.18% | No | 37 / 33 |
| 748935664 | 1 | +0.22% | −0.10% | +2.17% | No | 16 / 15 |
| 748935664 | 2 | +4.07% | −25.91% | −0.81% | No | 15 / 15 |
| 580014012 | 0.1 | 0.00% | 0.00% | 0.00% | Yes | 9 / 9 |
| 580014012 | 1 | 0.00% | 0.00% | 0.00% | Yes | 9 / 9 |
| 580014012 | 2 | 0.00% | 0.00% | 0.00% | Yes | 9 / 9 |
| 111065103 | 0.1 | 0.00% | 0.00% | 0.00% | Yes | 8 / 8 |
| 111065103 | 1 | 0.00% | 0.00% | 0.00% | Yes | 6 / 6 |
| 111065103 | 2 | 0.00% | 0.00% | 0.00% | Yes | 9 / 9 |

**Trace of the contact-closing outlier.** The initial state hashes match. Both arms merge plates
4 and 12 at 14 Myr; baseline applies the old blended rate, while candidate retains plate 4's last
solved rate. At the first subsequent solve, at 16 Myr, boundary-class counts differ on all 14
plates, but both solves converge. The first topology-path difference follows at 22 Myr: baseline
splits plate 4, while candidate does not split then and next splits plate 16 at 24 Myr. This is
consistent with the lagged-rate → classification → topology pathway. It does not establish that
this event alone explains the full 50 Myr contact-closing delta. The compact trace is in
[merge_closing_trace_summary.json](merge_closing_trace_summary.json); raw traces remain in scratch.

The separate timestep tolerance also fails in both arms. Across the 1 and 2 Myr runs relative
to 0.1 Myr, mean closing changes by as much as 192%, contact force by 41.5%, mean solved rate
by 39.7%, peak closing by 100.3%, and accumulated driving torque by 16.2%. Those step-size
effects occur in both arms and overlap topology changes, so they cannot be attributed solely to
the merge-rate change. Maximum Hm ledger closure is 0.0000121 of the corresponding inventory;
coarse-step continental balance residuals remain material and are reported in the compact JSON.

**Assessment.** The current-main results again show that most global motion and endpoint metrics
stay close, while topology and contact response diverge in some pairs. The selected merge-rate
rule is supported by the focused handoff trace, but this broader matrix is not an unchanged-
trajectory result. The preregistered >10% contact-closing difference still needs a per-step trace
before #547 can be accepted. No tuning was applied.

### Design decision: kinematics presented to the first solve after a merge

**Principle 8** requires merge topology and physical motion to meet at a clear interface. The
merge changes plate identity, geometry, and material; it does not supply a new physical rate.
RFC 0012's staggered step evaluates classifications and rate-dependent loads from the
pre-motion state, solves plate motion, then applies topology. Its specified lagged response is
an approximation; it does not require a nonlinear classification/solve loop each time topology
changes.

Options considered:

1. **Keep the surviving plate's last solved rate as the next-step kinematic estimate.** This
   retains a rate that came from the prior physical solve and avoids inventing an impulse or
   minimum speed. The next solve samples this rate for classifications and load linearization,
   then replaces it with the rate derived from the merged plate's current loads.
2. **Set the survivor's rate to zero.** This avoids importing the absorbed plate's motion, but
   injects an arbitrary stop at the topology event and can change convergence/divergence labels.
3. **Blend the two rates by angular momentum and clamp.** This preserves an inertial quantity
   that RFC 0012 explicitly removes from geological-time motion and can impose a rate floor.
4. **Iterate classification and motion until they agree.** This could address the dependency
   between class-dependent loads and solved rates, but it introduces a nonlinear active-set
   problem with convergence, tie-breaking, and runtime requirements. It is beyond #547's
   merge-only deliverable and is not established by the current replay.

**Selected rule for #547:** retain the survivor's last solved rate through the topology
operation; on the first following motion pass, classify from that lagged rate and solve from
post-merge loads. Do not blend with the absorbed rate, apply a rate floor, or iterate the
classification/solve in this issue. This is a numerically explicit lagged interface, not a claim
that the retained rate is the post-merge equilibrium. It follows RFC 0012's provisional
staggered coupling and keeps the issue's scope limited to removing the unphysical merge
assignment.

**Check against the traced merge.** The matched seed-`580014012` trace exercises this exact
rule. The prior survivor rate is 0.377 cm/yr in both arms; the baseline merge changes it to
0.500 cm/yr through blending and clamping, while the candidate retains 0.377 cm/yr. At the next
solve, classification changes on four plates (including the merged plate), while the merged
plate's solved rate differs by only 0.009 cm/yr. Thus the rule is implementable and the rate
handoff is causal for classification in this event. It does **not** make classification
invariant across the old and new merge behavior; the resulting contact-response and topology
differences are reported above and remain subject to #547's acceptance decision.

The 24-run matched replay already tests this selected rule because its candidate is precisely
the merge-only implementation retaining the survivor rate. Repeating those same runs would add
no evidence. The matrix found 3/12 different topology sequences, contact metrics beyond 10% in
some pairs, and did not meet RFC 0012's 10% timestep criterion in either arm. This rules out
claiming acceptance from the current evidence; it does not justify restoring the blend or
introducing a classification iteration. Any future iteration should be a separately scoped RFC
decision with a convergence invariant and a new preregistered replay.

### Matched merge-rate handoff trace

The preregistered diagnostic was posted to [#547](https://github.com/adubey/mantle-bloom/issues/547#issuecomment-6400942895).
It tests whether the rate written at a merge changes the next step's boundary classification
or solve. Principle 8 applies: each process should follow its own physics and connect through
clear interfaces. A merge's historical momentum blend and floor feed a later physical
classification before the quasi-static solve replaces plate rates.

**Controls and validity.** Compare baseline `f155b7eb36727f452139955b13350593db57cbcd` and
merge-only candidate `45b035c34ed0c370e52d1386e17c14318e893d5a`; generate seed `580014012` at
node density 4; use `PYTHONHASHSEED=0`, 2 Myr steps, and stop at 18 Myr. The first collision
merge is the same in both arms at 16 Myr: plates 8 and 13 merge into surviving plate 8. Initial
state SHA-256 fingerprints match, and all nine topology events through 18 Myr match, including
the subsequent split at 18 Myr. Continental, craton, mobile-cover, and Hm ledger residuals at
18 Myr remain at rounding scale. The matched trace is
[merge_rate_handoff_trace.json](merge_rate_handoff_trace.json); the instrumentation-only runner
is [`bin/debug/trace_rfc0012_merge_rate_handoff.py`](../../bin/debug/trace_rfc0012_merge_rate_handoff.py).

**Outcome.** Before the merge, survivor 8 has a rate of 0.377 cm/yr in both arms. At the merge,
baseline's momentum calculation and floor set it to 0.500 cm/yr; the candidate retains
0.377 cm/yr. The vector difference is 0.613 cm/yr. On the next solve (step 9, advancing 16–18
Myr), boundary-class counts differ on plates 0, 8, 10, and 12, before the solver overwrites
their rates. On the merged plate 8, convergent nodes change 103→112, divergent 16→0, transform
1,699→1,706, and subducting 9→21. Across all plates, convergent counts change 1,661→1,670,
divergent 1,003→975, transform 11,285→11,304, and contested stays 290. Driving-torque norm
changes −0.50%; maximum blocked-contact force changes −0.37%. The largest per-plate solved-rate
change is 0.088 cm/yr (plate 12); survivor 8's solved rate changes 0.009 cm/yr.

This supports the predicted causal pathway for this event: the merge assignment changes the
survivor's stored rate, and that rate changes boundary classes at the first following solve.
The event sequence remains matched through that solve. This one seed and one merge do not
establish how often the effect explains the 50 Myr replay differences, nor do they satisfy the
issue's timestep acceptance check.

**Reproduction.** Run the baseline and candidate commands one at a time with the corresponding
`<model-root>` and scratch output:

```sh
PYTHONHASHSEED=0 PYTHONPATH=<model-root>/backend \
  <model-root>/backend/.venv/bin/python \
  <trace-script> --model-root <model-root> --seed 580014012 --step-myr 2 \
  --stop-myr 18 --keep-plate 8 --absorb-plate 13 --out <scratch-output>.json
```

The baseline raw trace is labeled `merge-handoff-baseline-v2.json`, SHA-256
`dd1b101482157aad02c1bb5be8b87e8d28cbed94f78463dcec66ff8659fe544b`; the candidate is
`merge-handoff-candidate-v2.json`, SHA-256
`31f6b4b0eb789e9270c4cdebf94d16eb74a73e3be12761a0fb266eaa21e19b8e`. An initial instrumentation
attempt captured the solve on the merge step instead of the following step; it is excluded.
The corrected step selector was rerun serially on both arms. Both corrected invocations emitted
an `astropy_healpix.neighbours_nested` invalid-value runtime warning; the state
fingerprints, event histories, and ledger validity checks above matched.

## Historical matched replay preregistration (merge and split changes together)

**Question.** Does removing merge angular-momentum blending and split mantle-pole rate
assignments change physical motion, deformation, material closure, or geological outcomes,
once each resulting plate is advanced by the quasi-static load solve?

**Expected direction and size.** Expect no change in the next motion solve or in the subsequent
world trajectory. The solver is quasi-static and replaces every live plate's rate from the
current loads and resistance before rotating it; topology runs after the motion and deformation
phases. The removed state assignments should therefore affect only the stored rate between a
topology event and the next solve. The primary comparison is exact equality of per-step solved
rates, forces, deformation, event ordering, and final inventories for each matched baseline and
implementation run. Continuous step-size comparisons use RFC 0012's 10% tolerance. Report any
nonzero trajectory difference and its first checkpoint without tuning the model.

**Controls.** Generate the same world independently in each arm using `generate_world(seed,
node_density=4.0)`. Use seeds `326494039`, `748935664`, `580014012`, and `111065103`; advance to
50 Myr with step sizes 0.1, 1, and 2 Myr. This matches the RFC 0012 timestep-check horizon and
provides 500, 50, and 25 steps per run, respectively. Compare baseline
`402f00944251f383e7baf5ef7b418e8d606cde63` with the issue #547 implementation commit recorded
below. Run baseline and implementation serially, paired by seed and timestep, with only one
simulation active at a time. Set `PYTHONHASHSEED=0` in both arms. Record checkpoints every
10 Myr. Exclude 5 Myr.

**Metrics.** Record accumulated driving and residual torque response; mean and peak solved
plate rate; relative contact closing/slip; per-step Hc/Hm deformation and area-weighted Hc/Hm
distributions; land fraction; plate count; topology event time/order; and continental-material,
craton, mobile-cover, and Hm-ledger closure. Report endpoint values and checkpoint trajectories.
Check every continuous metric across 0.1/1/2 Myr against the RFC 0012 10% tolerance, and report
topology event-time/order differences separately. Record runtime and peak RSS for every run.

**Decision.** Because the proposed change removes state that G1 found unused by the next physical
solve, any change in a solved-load metric or a material ledger is unexpected and requires finding
the first divergent step before this issue is accepted. A state-rate difference immediately
after a topology event is expected; it must not feed subsequent motion or deformation.

## Results

Implementation replay commit: `9804addb2b1f29c53b4e22c473db18220c74a248`
(`Fix Python hash seed in replay controls`). The final code commit is
`b3f4e7db3d329dfd053c6c035ddd1df428bce785`; it adds a load-response regression test without
changing replayed model code.

### Preregistration amendment

The initial preregistration used a 500 Myr horizon, matching G1. A baseline 0.1 Myr pilot was
interrupted after more than an hour, before its first 50 Myr checkpoint, while it was routing
hydrology; it produced no comparison metrics. The 500 Myr, 24-arm matrix would take multiple
days at that measured pace. Before any replay metrics were collected, the horizon was changed
to 50 Myr, the existing RFC 0012 timestep-check horizon. Seeds, step sizes, metrics, comparison
rules, and one-simulation-at-a-time execution remain unchanged.

The first 50 Myr pair (seed `326494039`, 0.1 Myr) is retained as a discarded instrumentation
pilot, not as acceptance evidence. The baseline and implementation reached the same four split
events at 2.0, 4.1, 8.0, and 11.4 Myr, but their later event logs and continuous metrics diverged.
The runner took a shallow copy of cumulative `phase_budget` entries, so its deformation deltas
were invalid; it also recorded all world events under the topology-event field. Both defects are
fixed before repeating the pair. The other metrics from that pilot are provisional: mean solved
plate rate changed by -3.6%, mean blocked-contact closing by -6.9%, peak blocked-contact closing
by -12.2%, and mean blocked-contact force by -3.8%; Hc/Hm distribution shifts were at most 1.1%.
These are reported only to preserve the exploratory run, not to accept or reject the change.

The corrected 12-pair matrix is complete. All arms used `PYTHONHASHSEED=0`, node density 4,
50 Myr, 10 Myr checkpoints, and the same generated-world seeds. Full outputs and logs remain in
scratch storage; [corrected_replay_summary.json](corrected_replay_summary.json) contains the
matched metric deltas, topology differences, run time, peak RSS, and ledger closure errors.

| Seed | Step (Myr) | Topology events (base/candidate) | Largest relative motion-metric change |
|---:|---:|---:|---|
| 326494039 | 0.1 | same (10/10) | peak blocked-contact closing −12.2% |
| 326494039 | 1 | different (13/12) | mean blocked-contact closing −30.2% |
| 326494039 | 2 | different (10/12) | mean blocked-contact force −31.7% |
| 748935664 | 0.1 | different (15/14) | peak blocked-contact closing +16.6% |
| 748935664 | 1 | different (22/21) | mean blocked-contact closing −24.3% |
| 748935664 | 2 | different (24/21) | mean blocked-contact closing −39.2% |
| 580014012 | 0.1 | same (12/12) | peak blocked-contact closing +127.2% |
| 580014012 | 1 | different (13/11) | mean solved plate rate −25.9% |
| 580014012 | 2 | same (16/16) | mean blocked-contact closing −53.2% |
| 111065103 | 0.1 | same (9/9) | mean blocked-contact closing +184.0% |
| 111065103 | 1 | different (21/12) | peak blocked-contact closing +206.6% |
| 111065103 | 2 | different (13/11) | peak blocked-contact closing +31.8% |

Eight of the 12 pairs have different topology event sequences. Every pair has a different saved
10 Myr checkpoint; the runner did not save sub-checkpoint state or per-step motion metrics, so
these runs do not identify the exact first divergent solve. The earliest event differences in
the unequal sequences occur between 2.3 and 8 Myr. The largest relative percentages can be
amplified by small baseline contact values; absolute values and all metric deltas are in the
JSON summary. Final land-fraction deltas range from −0.01134 to +0.00160. Maximum endpoint
area-weighted Hc/Hm quantile shift is 4.83%. Material-ledger closure errors are recorded per
pair in the JSON summary.

These results reject the preregistered expectation of an unchanged trajectory. The existing
code confirms several couplings: `torque.classify_boundary_nodes` classifies boundary loads
from `Plate.omega`; `merge_split.maybe_split_plate` compares mantle flow with that stored rate;
and merge eligibility and probability use stored relative closing rates. `step_world` calls
motion/deformation before topology changes, so rates changed by topology remain inputs to the
next solve and topology pass. Removing those assignments therefore changes physical loads and
later topology, even though rates are subsequently solved. The implementation is not accepted
pending a step-level trace and a revised design that decouples topology bookkeeping from
physical classification without restoring an unsupported angular-momentum blend.

The 50 Myr 0.1 Myr baseline pilot that ran before the replay instrumentation fix is exploratory
only and excluded from this table. It was discarded because its deformation snapshots were
shallow copies of cumulative phase budgets and its event filter included non-topology events.

## Rate-classification coupling trace preregistration

**Question.** Does the post-topology stored rate change the next solve's boundary classification
and load assembly, explaining the matched replay divergence?

**Hypothesis and prediction.** `torque.classify_boundary_nodes` classifies contacts from stored
`Plate.omega`. For seed `326494039` at 0.1 Myr steps, the existing matched traces show the same
first split at 2.0 Myr in both arms. Baseline assigns split poles to the daughters while the
candidate retains the parent's rate. If this rate-dependent classification is causal, the first
solve after the split (2.1 Myr) will show at least one changed per-plate boundary-class count;
the associated torque or solved rate may also differ. The direction and size of the continuous
change are not predictable in advance.

**Primary metric and decision threshold.** Compare each plate's convergent, divergent, transform,
and contested node counts on the 2.1 Myr solve. Any count difference supports the classification
pathway. Secondary measures are each plate's pre-solve angular velocity, closing-rate summary,
driving-torque components, solved rate, blocked-contact response, and topology events. If all
class counts match and every per-plate load and solved rate agrees within 1e-8 relative (with
1e-12 absolute tolerance), this pathway is not supported for the first split; later split-gate
and merge-decision rate dependencies remain alternatives.

**Controls.** Run baseline `402f00944251f383e7baf5ef7b418e8d606cde63` and candidate
`9804addb2b1f29c53b4e22c473db18220c74a248` independently from `generate_world`, with seed
`326494039`, node density 4, `PYTHONHASHSEED=0`, 0.1 Myr steps, and a 2.1 Myr horizon. Run one
simulation at a time. An instrumentation-only scratch script wraps boundary classification and
records per-step state without changing model calculations. Use the selected checkout's Python
and code on `PYTHONPATH`:

```sh
PYTHONHASHSEED=0 PYTHONPATH=<model-root>/backend <model-root>/backend/.venv/bin/python \
  <scratch>/trace_rfc0012_rate_coupling.py --model-root <model-root> --out <scratch-output>
```

**Validity checks.** Initial generated-world state and events through the 2.0 Myr topology step
must match exactly. Record final plate count, topology events, and material-ledger closure; if
initial state or conservation validity differs, do not use the trace to support the hypothesis.
The result is a causal diagnostic only, not an acceptance replay or candidate-fix test.

### Trace outcome

The preregistered prediction was supported. Baseline and candidate classifications, motion
diagnostics, and topology events match through the 2.0 Myr step. That step splits plate 11 into
plates 11 and 19. After topology, their stored rates differ: baseline uses the fitted daughter
poles; candidate leaves both daughters at the former parent rate. On the first following solve
at 2.1 Myr, boundary-class counts differ on 16 of 20 plates. For plate 11, convergent/divergent/
transform counts change from `22/311/236` to `29/354/186`; for plate 19 they change from
`333/136/181` to `334/196/120`. Aggregate counts change from `11564/7760/27675/1047` to
`11490/7727/27786/1046` (convergent/divergent/transform/contested).

The first post-split solve's total driving-torque norm changes by only +0.20%, while the maximum
solved rate changes from 4.336 to 4.852 cm/yr (+11.9%). The residual torque is still tiny in
absolute terms (`1.52e12` vs `2.29e12 N m`, compared with `~6.2e26 N m` driving torque). There
are no blocked contacts at this step. Ledger closure remains at rounding scale: continental
balance errors are +2,048 and −1,024 m³, craton errors are zero, and Hm closure errors are at
most 7,680 m³. The complete per-plate trace and closure values are in
[rate_coupling_trace.json](rate_coupling_trace.json).

This directly demonstrates that the changed post-topology rate alters the next physical
classification and load solve. The late-run divergence is therefore not only a replay or
instrumentation artifact. It does not establish that classification is the only later pathway:
split gates and merge decisions also read stored rates. The issue's current design cannot meet
its unchanged-trajectory expectation by removing assignments alone; any further model design
change needs the project owner's direction under RFC 0012.

The trace used the instrumentation-only scratch script `trace_rfc0012_rate_coupling.py` with
the same baseline and candidate controls above, but a 2.1 Myr horizon. The script wrapped
`torque.classify_boundary_nodes` to record counts and rates, then recorded the ordinary motion
diagnostics after each step. Its command was:

```sh
PYTHONHASHSEED=0 PYTHONPATH=<model-root>/backend <model-root>/backend/.venv/bin/python \
  <scratch>/trace_rfc0012_rate_coupling.py --model-root <model-root> --out <scratch-output>
```

## Later topology gate trace preregistration

**Question.** Does the stored-rate-dependent split gate explain the first different topology
event in a matched replay?

**Hypothesis and prediction.** In the 50 Myr, 0.1 Myr pair for seed `748935664`, the first
topology mismatch is a split of Plate 6 at 2.3 Myr in the candidate versus 2.4 Myr in baseline.
The split gate compares the mantle-flow velocities with `Plate.omega` using
`_fit_residual_rms`, then applies clustering, pole-separation, minimum-size, and deterministic
rift-outcome gates. If the stored rate is the cause of the one-step difference, the Plate 6
RMS residual relative to its threshold (or a downstream gate input) will differ at 2.3 Myr,
with the candidate clearing all gates at that step while baseline does not; baseline should
clear the corresponding gates at 2.4 Myr. Direction and magnitude of residual differences are
not assumed.

**Primary metric and decision threshold.** Record every Plate 6 split attempt from 2.2 through
2.4 Myr in both arms: stored omega, RMS residual, size/stress relaxation and resulting residual
and pole-separation thresholds, clustering/refinement outcome, fitted pole separation, minimum
node gate, deterministic rift outcome, and whether a split or failed rift occurs. The primary
comparison is the first timestep and gate at which the split outcomes differ. Any difference
at or before the residual gate supports a stored-rate contribution; matching gate inputs and
outcomes through 2.4 Myr falsifies this specific explanation. This trace does not test merge
decisions or establish that the split gate is the only later coupling.

**Controls.** Run baseline `402f00944251f383e7baf5f7b418e8d606cde63` and candidate
`9804addb2b1f29c53b4e22c473db18220c74a248` independently from `generate_world`, with seed
`748935664`, node density 4, `PYTHONHASHSEED=0`, 0.1 Myr steps, and a 2.4 Myr horizon. Run
one simulation at a time. An instrumentation-only scratch script wraps the split gate to
record values and decisions without changing model calculations.

**Validity checks.** The generated initial states and topology event history must match through
2.2 Myr. Record final plate count, topology events, and material-ledger closure; if initial state
or conservation validity differs, do not use the trace to support the hypothesis. This is a
causal diagnostic only, not an acceptance replay or candidate-fix test. Preregistration was
committed and posted to issue #547 before either arm was run.

### Trace outcome

The validity checks passed. Independently generated plate fields and mantle centers were equal;
Plate 6's gate snapshots through 2.2 Myr matched exactly; and both arms split Plates 0, 2, and 4
at 2.0, 2.1, and 2.2 Myr. Candidate and baseline were run serially through 2.4 Myr.

At the Plate 6 decision for the 2.3 Myr topology step, candidate `Plate.omega` was
`[6.534, -1.086, -1.891]e-9 rad/yr`. Its residual RMS was `6.853e-9 m/s`, above the
`6.456e-9 m/s` threshold, and the split cleared clustering, refinement, pole separation, and
minimum-size gates. Baseline omega was `[6.790, -0.743, -2.461]e-9 rad/yr`; its residual RMS
was `6.272e-9 m/s`, below the same `6.456e-9 m/s` threshold, so it did not cluster or split.
At the following step, baseline residual RMS rose to `7.061e-9 m/s`, above its `6.456e-9 m/s`
threshold, and Plate 6 split at 2.4 Myr. The candidate's split therefore occurred one step
earlier, as predicted. Size/stress relaxation was effectively identical (`0.31449`), as was the
residual threshold; the changed stored rate and resulting flow-residual comparison are the first
recorded gate difference.

Final plate count was 12 in both arms. Continental, craton, mobile-cover, and Hm ledger closure
errors remained at rounding scale (largest absolute Hm closure error was `8,192 m³`). This
supports the split-gate pathway as a contributor to the topology divergence. It does not isolate
the counterfactual effect of changing only omega from other small trajectory differences, nor
does it rule out later merge-decision coupling. The compact gate values and closure results are
in [split_gate_trace.json](split_gate_trace.json); full outputs and the instrumentation-only
runner remain in scratch storage.

Baseline replay commit: `402f00944251f383e7baf5f7b418e8d606cde63`. Candidate model code is
unchanged from `9804addb2b1f29c53b4e22c473db18220c74a248`; the trace ran with preregistration
commit `854f3733eaa18ad1b0b17bb7199dc8fdcd12a45f` at HEAD. Both arms used seed `748935664`, node
density 4, `PYTHONHASHSEED=0`, 0.1 Myr steps, and were run serially to 2.4 Myr.

Baseline replay commit: `402f00944251f383e7baf5ef7b418e8d606cde63`. Candidate checkout at run:
`ea595b210e0f37a95f6a3f2ecb7328a1d1ce0b30`; its model code is unchanged from
`9804addb2b1f29c53b4e22c473db18220c74a248` (the added commits only preregister and record the
probe). Both runs used seed `326494039`, node density 4, 0.1 Myr steps, `PYTHONHASHSEED=0`, and
were run serially. Runtime was 132.5 seconds baseline and 226.2 seconds candidate; the
instrumentation-only runtime difference is not a model performance comparison.

## Merge-only matched replay preregistration

**Question.** Does removing the merge angular-momentum blend change motion, deformation,
material closure, or geological outcomes when split-rate assignments remain unchanged?

**Prediction.** A merge-only candidate leaves the survivor's stored rate in place until the
next motion pass, while baseline blends the survivor and absorbed plate's angular momentum and
clamps the resulting rate. Because boundary classification can read stored rates before the
simultaneous motion solve, post-merge classifications or solved loads may differ. No direction
or size is assumed. The run will record the first changed classification/load after each merge
where captured, along with the complete preregistered replay metrics. The prior combined-change
matrix cannot answer this question.

**Metrics and decision threshold.** Compare accumulated driving and residual torque response,
mean and peak solved plate rate, relative contact closing/slip, per-step Hc/Hm deformation and
area-weighted Hc/Hm distributions, land fraction, topology event time/order, and
continental-material, craton, mobile-cover, and Hm-ledger closure. Compare continuous metrics
across 0.1/1/2 Myr steps with RFC 0012's 10% tolerance; report topology differences separately.
If any continuous metric exceeds tolerance, identify its first divergent step and cause before
acceptance or tuning. Ledger closure must remain at numerical-rounding scale.

**Controls.** Use `generate_world(seed, node_density=4.0)` for seeds `326494039`, `748935664`,
`580014012`, and `111065103`; compare baseline `f155b7eb36727f452139955b13350593db57cbcd`
with merge-only model implementation `f7f4c0ae2bedfbadfc57bcba53661ba24adaae3c` (the current
candidate checkout may include later test or diagnostic-only commits). Run to 50 Myr at 0.1,
1, and 2 Myr steps, serially with only one simulation active. Set `PYTHONHASHSEED=0`; record checkpoints
every 10 Myr. Exclude 5 Myr. Save compact results and exact commands under this directory; keep
full outputs in scratch storage.

**Runner and decision rule.** The candidate's
`bin/debug/replay_rfc0012_merge_rate.py` runner can target either arm through
`MANTLE_BLOOM_REPLAY_BACKEND`; use the matching worktree's Python environment for each run.
Example baseline command:

```sh
MANTLE_BLOOM_REPLAY_BACKEND=<baseline-root>/backend PYTHONHASHSEED=0 <baseline-root>/backend/.venv/bin/python <candidate-root>/bin/debug/replay_rfc0012_merge_rate.py --seed <seed> --step-myr <0.1|1|2> --myr 50 --checkpoint-myr 10 --out <scratch-output>.json
```

Compare checkpoint trajectories and endpoints for the preregistered measures; continuous
metrics must stay within 10% across step sizes. Any larger difference must be traced to its
first divergent step before acceptance. Report topology sequence changes separately and
require all four ledger balance errors to remain at numerical-rounding scale. The runner stores
one compact JSON summary per arm; full output and saves remain outside Git. This preregistration
and runner are committed before the first matrix arm.

**Status.** Preregistered for the narrowed candidate. The matrix has not started.
Before the matrix, run one 0.1 Myr instrumentation smoke on each arm to validate backend
selection, the runner's schema and ledger fields; exclude both from all baseline/candidate
comparisons and decision metrics. The initial smoke exposed that ordinary world events include
fault and earthquake messages, so they cannot serve as topology outcomes; discard those smoke
outputs. The runner now taps `merge_split.apply_topology_changes` directly. Repeat the smoke on
both arms with this corrected topology-event capture before starting the matrix.

### Focused timestep refinement: 0.05 Myr (2026-10-26)

The preregistered check was posted to [issue #547](https://github.com/adubey/mantle-bloom/issues/547#issuecomment-6408770512) before either 0.05 Myr run. It compares seed `748935664` at node density 4, `PYTHONHASHSEED=0`, 50 Myr, and 10 Myr checkpoints. The 0.1 Myr runs are reused from the current-main matrix: baseline `f04cbc0772eca6d2484f7477b20004deb368206c`, candidate model code `d553c3cfc05db4372d0bb33ec7e58aa2390b793b`. The 0.05 candidate ran from `fb516de9ee61cb54bb29a9a7050a9fa1f6c07481`, which adds only the report to the same model code. The two new runs ran serially and completed successfully.

| Continuous metric, 0.05 vs 0.1 Myr | Baseline | Candidate | 10% limit |
|---|---:|---:|:---:|
| Accumulated driving torque | +5.17% | +5.13% | pass |
| Accumulated residual torque | −7.06% | −2.99% | pass |
| Mean solved rate | −0.66% | −2.29% | pass |
| Peak solved rate | −3.08% | 0.00% | pass |
| Mean absolute contact closing | +46.08% | +38.67% | **fail** |
| Peak absolute contact closing | −3.39% | −3.39% | pass |
| Mean contact force | +28.22% | +13.49% | **fail** |

At every 10 Myr checkpoint, endpoint land fraction and area-weighted Hc/Hm quantiles stayed within 4.57% of the 0.1 Myr values. Cumulative deformation did not stay within tolerance: for example, divergent Hc deformation changed −27.7% baseline / −24.2% candidate, failed-rift Hc changed −349% / −260%, and failed-rift Hm changed −798% / −627%. These component budgets also reflect topology changes, so they should not be read as fixed-topology convergence tests.

The candidate/baseline mean-closing difference changes from +0.51% at 0.1 Myr to −4.59% at 0.05 Myr. Its absolute size grows from `1.24e-13` to `1.62e-12 m/s` (about 13×) and changes sign. Thus refinement does not bound the contact response at the preregistered tolerance. Peak closing remains equal between arms at both step sizes.

Topology differs with step size and between arms. In both arms, the Plate 2 split occurs at 2.05 Myr with 0.05 Myr steps versus 2.1 Myr with 0.1 Myr steps. At 50 Myr, baseline has 38 events and 21 plates at 0.05 Myr, versus 37 events and 20 plates at 0.1 Myr; candidate has 40 events and 20 plates at 0.05 Myr, versus 33 events and 20 plates at 0.1 Myr. At 0.05 Myr the first baseline/candidate event-log mismatch is 16.8 versus 18.0 Myr. At 0.1 Myr it is 34.4 versus 34.0 Myr.

Both new runs reached 50 Myr with exit code 0. Runtime / peak RSS was 5,137 s / 1.72 GiB baseline and 5,363 s / 1.53 GiB candidate. Material and Hm ledger closure remained at rounding scale; the largest Hm residual was `1.58e5 m³` against inventory `2.71e19 m³`. Repeated glacier-depth warnings were emitted in both arms; neither run failed. Full outputs remain in scratch. The compact checkpoint and metric comparison is [timestep_refinement_005_comparison.json](timestep_refinement_005_comparison.json).

**Assessment.** The prediction of bounded refinement is not supported. Torque, solved-rate, land, and thickness summaries stay within 10%, but contact closing, contact force, and several deformation budgets do not. Topology also shifts with step size. This is one-seed diagnostic evidence; it neither establishes full RFC acceptance nor supports metric tuning.

**Reproduction.** With the runner's 0.05 Myr option, run the baseline and candidate serially using the same command shape as the matched replay, changing `--step-myr` to `0.05` and writing each arm to a separate scratch output:

```sh
MANTLE_BLOOM_REPLAY_BACKEND=<arm-root>/backend PYTHONHASHSEED=0 <arm-root>/backend/.venv/bin/python bin/debug/replay_rfc0012_merge_rate.py --seed 748935664 --step-myr 0.05 --myr 50 --checkpoint-myr 10 --out <scratch-output>.json
```

### Follow-up refinement probe: 0.025 Myr (2026-10-26)

The follow-up question and readout method were recorded in this analysis before either run.
We asked whether continuous measures move toward stable values across 0.1→0.05→0.025 Myr.
The same seed, density, hash seed, horizon, and checkpoints were used.
The two 0.025 runs completed serially. This is a one-seed trend check, not a new pass/fail
threshold or proof of convergence.

| Measure | Baseline: 0.1→0.05; 0.05→0.025 | Candidate: 0.1→0.05; 0.05→0.025 | Adjacent absolute change shrinks? |
|---|---:|---:|:---:|
| Accumulated driving torque | +5.17%; −3.27% | +5.13%; −3.74% | Yes, both |
| Accumulated residual torque | −7.06%; +4.63% | −2.99%; +6.34% | Baseline only |
| Mean solved rate | −0.66%; −0.46% | −2.29%; +0.33% | Yes, both |
| Peak solved rate | −3.08%; 0.00% | 0.00%; 0.00% | Baseline contracts; candidate stays flat |
| Mean absolute contact closing | +46.08%; +3.27% | +38.67%; −1.85% | Yes, both |
| Peak absolute contact closing | −3.39%; +32.63% | −3.39%; +32.63% | No, both |
| Mean contact force | +28.22%; −14.03% | +13.49%; −6.15% | Yes, both |

The candidate/baseline gap in mean absolute closing is +0.51% at 0.1 Myr, −4.59% at 0.05,
and −9.32% at 0.025. The within-arm candidate average changes little at 0.025, while its
difference from baseline grows toward 10%; these are separate comparisons.

Mean absolute closing changes sharply from 0.1 to 0.05 Myr, then changes only 3.27% in the
baseline and −1.85% in the candidate at 0.025 Myr. That measure appears to be settling in
this seed. Mean contact force also has a smaller absolute step-to-step change at 0.025, though
the baseline still shifts 14.03% from 0.05 to 0.025. Peak closing does not settle: it rises
32.63% in both arms on the final refinement. Residual torque also fails to contract in the
candidate. So the refinement trend is measure-dependent; it does not support a general claim
that the simulation has settled.

Topology remains timestep-sensitive: baseline event counts are 37, 38, and 32 with final plate
counts 20, 21, and 20 at 0.1, 0.05, and 0.025 Myr; candidate event counts are 33, 40, and 27,
with 20 plates at all three resolutions. Cumulative deformation components do not all
contract either, including boundary-retreat Hm and convergent-deformation Hc/Hm. Endpoint
land fraction shifts by roughly 1.7–2.5% between adjacent resolutions; median Hm thickness
is unchanged to the reported precision. Ledger residuals stayed at rounding scale. The runs
took 9,705 s / 1.74 GiB peak RSS (baseline) and 9,557 s / 1.84 GiB (candidate); both reached
50 Myr successfully. Full outputs remain in scratch. The compact comparison is
[timestep_refinement_0025_comparison.json](timestep_refinement_0025_comparison.json).

**Assessment.** Lowering the timestep substantially reduces the drift in mean contact
closing, which is evidence that this particular average is approaching a stable value. Other
metrics still oscillate or move farther at 0.025 Myr. Because topology events also change,
three resolutions for one seed are not enough to establish convergence of the full simulation.

## Peak contact closing trace

To locate the 32.63% rise in peak absolute contact closing from 0.05 to 0.025 Myr, I repeated
the seed-748935664, 50 Myr baseline/candidate pairs at both resolutions with a step-level peak
record. All four runs completed serially. The candidate replay used diagnostic revision
`933b54a`; its model code is unchanged from candidate revision `fb516de`. The baseline model
revision is `f04cbc0`.

| Arm | Step | Peak time | Peak contact | Absolute closing | Yielded at capacity | Nearby topology |
|---|---:|---:|---|---:|:---:|---|
| Baseline | 0.05 Myr | 4.200 Myr | 4–11 | 4.78624×10⁻¹⁰ m/s | Yes | Plate 4 split into 4 and 11 at 4.100 Myr; 2 steps before peak |
| Candidate | 0.05 Myr | 4.200 Myr | 4–11 | 4.78624×10⁻¹⁰ m/s | Yes | Plate 4 split into 4 and 11 at 4.100 Myr; 2 steps before peak |
| Baseline | 0.025 Myr | 2.125 Myr | 4–9 | 6.34803×10⁻¹⁰ m/s | Yes | Plate 4 split into 4 and 9 at 2.025 Myr; 4 steps before peak. Plate 6 split at 2.100 Myr, 1 step before peak |
| Candidate | 0.025 Myr | 2.125 Myr | 4–9 | 6.34803×10⁻¹⁰ m/s | Yes | Plate 4 split into 4 and 9 at 2.025 Myr; 4 steps before peak. Plate 6 split at 2.100 Myr, 1 step before peak |

At each resolution, baseline and candidate have the same peak value, contact pair, and event
timing. The higher 0.025 Myr peak is a different contact after a different split sequence; it
is not a candidate-only increase. Both peaks occur at yielded contacts whose force equals the
reported capacity. The 0.05 peak is two steps after the split that created its contact. The
0.025 peak is four steps after the split that created its contact and one step after a separate
split. This is consistent with a topology-sensitive peak, but timing alone does not establish
causation. The trace localizes the rise for this one seed; it does not show that peak closing
has converged or that all topology-sensitive peaks behave the same way.

The repeated runs took 4,563 s / 1,514 MiB (baseline 0.05), 4,914 s / 1,647 MiB (candidate
0.05), 9,234 s / 1,808 MiB (baseline 0.025), and 10,454 s / 1,330 MiB (candidate 0.025).
Ledger residuals remained near rounding scale. Full outputs remain in scratch. The compact
trace is [peak_closing_trace_comparison.json](peak_closing_trace_comparison.json); the local
run command used the existing replay runner with `--seed 748935664 --step-myr <0.05|0.025>
--myr 50 --checkpoint-myr 10`, the baseline or candidate backend, and one output file per run.
