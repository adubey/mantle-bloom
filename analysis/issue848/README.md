# Issue 848: where the RFC 0018 candidate search spends its time

Issue [#848](https://github.com/adubey/mantle-bloom/issues/848) asks whether RFC 0018's
candidate method can meet G3's 2% runtime limit, and what to do if it cannot. This report is
the first step: a function-level profile of the search on mature plates, plus a check of how
selective a rift-candidate marker based on recent extension (#836) could be. It is a profile,
not a G3 measurement. No production code, gate or RFC text changed.

## Profiling world and seed registration

The G0 snapshots are 2 Myr old and the G2 summaries kept no saves, so this report uses one new
world. Its seed is **not part of any registry, development set or held-out sample, and must not
be counted as one.** Further uses of it are development evidence.

| Item | Value |
|---|---|
| Seed | `1091931222`, the first four SHA-256 bytes of `mantle-bloom/rfc0018/profile/v1/00`, big-endian, masked to 31 bits |
| Checked absent from | every number in `analysis/**/*.json`, `analysis/**/*.md` and `rfc/*.md` at `origin/main` b72e4461 |
| Layout | cratons-first, relief from synthetic history, both passed explicitly (`synthetic_relief` defaults off, so the G0 snapshots had none) |
| Stepping | 1 Myr steps to 200 Myr, plates-only (`World.simulate_climate_biomes = False`), production great-circle splitter |
| Result | 11 to 34 plates, 336 s total, median 1.65 s per step |
| Code | `origin/main` b72e4461 plus the scripts in this PR |
| Saves (private) | `world_50myr`, `world_100myr`, `world_200myr`; SHA-256 `be337d44`, `0b5127a4`, `fa6cb174` (full sums in the private capture) |

A first run with the defaults took about 385 s per Myr (4 steps, stopped), so the world is
plates-only. That is also the likeliest reason G2's baseline averaged 2.5 s per step, but the G2
settings are not recorded, so that is an inference.

## Rebuild of the frozen mass-balanced generator

The G1/G2 generator's code is lost, but its output survived: `mass-balanced-probe.json` (SHA-256
`c57a813c`, the `source_probe_sha256` in `g0-prior-summary.json`) and the 20 G0 snapshots
(`LEGACY_HASHES.json` maps each recompressed `.bloom` back to the SHA-256 recorded in the G0
evidence). `bin/debug/rebuild_rfc0018_mass_balanced.py` rebuilds the generator and the `O/R/S`
prior from the written spec and compares every candidate with that probe. Compact result:
[mass-balanced-rebuild.json](mass-balanced-rebuild.json).

| | Rebuilt | Recorded |
|---|---:|---:|
| Opportunities / candidates | 68 / 812 | 68 / 812 |
| Structurally valid candidates | 97 | 97 |
| Opportunities with a valid candidate / positive `P_break` | 33 / 17 | 33 / 17 |
| Mean `P_break` (`mu = 1`) | 0.046925976976340836 | 0.046925976976340836 |

All 812 candidates match the probe on daughter node counts, cut cost, structural validity,
momentum score and net opening. What the evidence fixed that the written spec did not:

- The number of targets is `min(12, edges that cross the two flow clusters and have a nonzero
  Fiedler gap)`, spaced by `linspace(0.15, 0.50, count)`; one plate had only 8 such edges.
- A target is a share of plate mass for the low-Fiedler side, taken at the **first** node where
  the cumulative mass reaches it (not the closest).
- The mass-quantile cut is re-applied at every refinement iteration; the seed-edge midpoint of the
  prototype is not used. Identical partitions are dropped, though none repeated here.
- The conductance is still `1 / edge_cost`, so this reproduces the inverse-weight method. The
  direct-cost revision of #857 changes this one line.

This validates the rebuild on the G0 development snapshots only. It says nothing about the
applying callback (how a chosen cut becomes a split, and failed-rift thinning), which is still
not rebuilt, so G2 behaviour is not reproduced.

## Which candidate generator the profile measured

**This is the original G0 prototype, not the generator G1 and G2 froze.** `measure_rfc0018_g0.py`
(commit 236b159c) proposes 12 variants from sampled weak seed edges. The G1/G2 method instead
proposes up to 12 mass-balanced cuts (target fractions 0.15 to 0.50, deduplicated, weighted by
the potential-momentum score; see `analysis/issue797/README.md`). That code is in the G2 runner
and the #846 worktree, neither of which I could reach, so it was not profiled. The two share the
per-variant refinement loop (one initial solve, then up to 8 per variant: up to 97 solves per
decision, the same figure #848 gives for the frozen implementation), which is why the solve
share is probably representative. The proposal rule, the `usable` test and the variant
statistics (including the 96% iteration-cap share) may differ, so treat those as unconfirmed for
the frozen method. Re-run this profile against the frozen generator once its code is reachable.

## Method

`bin/debug/profile_rfc0018_candidate_search.py` loads `world_100myr`, steps it 5 unmeasured
warm-up steps (101 to 105 Myr) then 15 measured steps (106 to 120 Myr), and at each
measured step measures every plate the G0 prototype treats as a split opportunity
(`_is_opportunity`), from the same pre-split state as the G0 runner. It calls the prototype's
`choose_cut` from `measure_rfc0018_g0.py` (commit 236b159c on `origin/codex/issue-797-evidence-rfc0018-g0`,
file SHA-256 `3e5a489f`) unchanged and wraps its helpers with timers. Both mappings ran on every
opportunity: `inverse` is the prototype's `conductance = 1 / edge_cost`, `direct` is `conductance =
edge_cost` from the #850 revision. Production splitting continued unchanged.
Summary: [profile-summary.json](profile-summary.json), built by `bin/debug/summarize_rfc0018_profile.py`.

## Results

150 opportunities over 15 steps (13 distinct plates, median 1,177 nodes, max 2,143), each
measured under both mappings. Times are inclusive, so rows overlap.

| | Inverse (as in G0 prototype) | Direct (#850) |
|---|---:|---:|
| Search time, all opportunities | 179 s | 145 s |
| Median per opportunity | 1.11 s | 0.80 s |
| Sparse solves, total / median per opportunity | 14,248 / 97 | 12,941 / 91 |
| `eigsh` share of search time | 84.7% | 82.3% |
| `_fiedler_vector` (matrix build plus `eigsh`) | 88.0% | 85.9% |
| Resistance refits (`_node_resistance`) | 5.5% | 6.1% |
| Edge extraction (`_shared_edges`) | 5.3% | 6.5% |
| Opening, trace, component checks | under 1% each | under 1% each |
| Variants using all 8 fit iterations | 1,722 of 1,800 (96%) | 1,371 of 1,800 (76%) |
| Usable variants / accepted selected cuts | 0 / 0 | 0 / 0 |

1. **The sparse eigensolve is the cost.** About 85% of search time is `eigsh`, and the search
   makes about 97 solves per opportunity, which is the cap of 1 + 12 variants x 8 iterations.
   Resistance refits, edge extraction and trace validation together are about 12%. Speeding up
   anything other than the solve count or cost cannot reach the target.
2. **The per-variant fit loop almost never settles early.** 96% of variants (inverse) ran the
   full 8 iterations. This count cannot tell a loop that converged on the eighth pass from one
   that hit the cap, so whether the later iterations change the cut is unmeasured. That is the
   most valuable next measurement: if they rarely change it, a smaller cap or a warm start
   would remove most solves. It would still change the method's outputs at some tolerance, so
   it needs its own check against the G0/G1 results.
3. **The direct mapping is somewhat cheaper (about 19% less time) but not a different order.**
   It does not by itself bring the cost near 2%.
4. **Scale check against G2.** The 150 opportunities took 179 s of search time over 15 steps,
   about 12 s per step against a plates-only step of about 1.7 s: total step time about 8x
   higher from the search alone, with profiler overhead included. G2 measured 16.25x on different worlds and an
   uncontrolled setup, so this is the same order of magnitude, not a reproduction.
5. **Usable variants were 0 of 1,800 under both mappings.** The `usable_candidate` test in
   that script still includes the superseded `1.96 x RMS` opening screen, so this is not a
   verdict on the revised method; the accepted `O/R/S` prior was not evaluated here.

### Rift-candidate marker: recent extension is not selective

For each opportunity the script records the steps since the plate last called
`quad_tectonics._open_rift`. A value of 1 means the plate was stretched in the previous step. The data suggest the
split pass runs before a step's own deformation (a value of 0 is rare, 5 of 150), but I did
not check the call order in the code.

| Plate extended within | Opportunities | Share of search time (inverse) |
|---|---:|---:|
| 1 step | 100 of 150 (67%) | 73.2% |
| 3 steps | 124 of 150 (83%) | 85.7% |
| 5 steps | 138 of 150 (92%) | 93.0% |
| Never since warm-up | 3 of 150 (2%) | n/a |

A marker that only asks whether a plate was stretched recently would keep two thirds or more
of the search. On a mature world almost every large plate has some extending edge nearly every
step, so a plate-level on/off flag cannot cut the call count. A useful marker would have to be
graded and local: how much extension, how concentrated, and whether it sits near an already
weak line. #836 has not established what magnitude or localization is physically meaningful, and
it cannot be added to the accepted method without G0 and G1 evidence, so this report does not
choose one.

## What this changes in the plan

- **Priority:** reduce the number or cost of eigensolves per decision. Candidates, from least
  to most invasive: (a) measure how often iterations 2 to 8 change a variant's cut; (b) warm-start
  each solve from the previous eigenvector; (c) reuse a plate's candidate set across steps while
  its inputs are unchanged; (d) lower the cap or variant count. (b) and (c) can keep outputs
  within solver tolerance. (a) is needed first, because (d) changes the distribution.
- **Marker:** keep it as a design option, but as a graded signal feeding resistance, not as a
  gate on the search.
- **Failed-rift semantics:** not addressed here; the event model in PR #857 is still unaccepted.
- **#866:** unrelated to this question; the outlet-incision result is inconclusive and does not
  concern stretching.

## Limits

One world, one 15-step window at 106 to 120 Myr, 13 distinct plates; no uncertainty is
reported because opportunities within a world are not independent. The world is a baseline-arm
history, so plate sizes differ from a graph-arm run. The timers add overhead and the helper
shares overlap. Nothing here is a G3 result: G3 needs paired, serialized, un-instrumented full
runs after the method settles.

## Reproduce

```bash
backend/.venv/bin/python bin/debug/generate_issue848_world.py --out ~/mantle-bloom-captures/issue848
backend/.venv/bin/python bin/debug/profile_rfc0018_candidate_search.py \
    --save ~/mantle-bloom-captures/issue848/world_100myr.mbworld \
    --g0-script <measure_rfc0018_g0.py from origin/codex/issue-797-evidence-rfc0018-g0> \
    --out ~/mantle-bloom-captures/issue848/profile.jsonl --warmup-steps 5 --steps 15
backend/.venv/bin/python bin/debug/summarize_rfc0018_profile.py \
    ~/mantle-bloom-captures/issue848/profile.jsonl analysis/issue848/profile-summary.json
```

The raw `profile.jsonl` (300 rows, SHA-256 `e8e5cbef`) stays in the private capture.
