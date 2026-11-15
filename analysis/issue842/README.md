# #842 round 3: the trend-adaptive march (`trend-1k`) against fixed-250 yr D2

Question: can the trend-adaptive march with a 1 kyr longest substep replace fixed-250 yr D2
as the production default (RFC 0015, D2)? The preregistration is `preregistration.md`, with
nine amendments. Scores are in `round3.json` (round 3b) and `round3_nov8.json` (the earlier
round 3). Harness: `bin/debug/run_842_round3.sh`, `score_842_round3.py`,
`score_842_round3_nov8.py`, `probe_822_march_cost.py`. Captures: `~/mantle-bloom-captures/issue-842-r3b/`
(round 3b, with `SHA256SUMS`) and `issue-842/` (the earlier round 3).

## Verdicts

Both round 3s, under every rule set. Amendments 6–9 were each set by the owner after
seeing results (6 and 7 after the earlier run, 8 and 9 after round 3b).

| Run | As written | Per-step thresholds (owner's first relaxations) | Amendments 6–8 | Amendments 6–9 (decides) |
|---|---|---|---|---|
| Earlier round 3 (2026-11-08/09, tree `a81c4ccc`) | Fail (rules 1, 2, 5) | n/a | Fail (5b: flips) | **Pass** (6 + 7) |
| Round 3b (this plan, 108 runs, one tree) | Fail | Fail | Fail (gap, 3a) | **Pass** |

Where the ice-age saves matter: round 3b's two failures under 6–8 were both on the ice-age
saves, and the trend march has no passing evidence there (below).

## Round 3b by threshold: old and new

| # | Criterion | Old threshold | New threshold (gating) | Result |
|---|---|---|---|---|
| 1 | Glacier ledger (E1, E5; first 10 steps) | ≤ 0.1% on every step | Mean ≤ max(0.1%, `d2` mean + 0.01%) | Old: ❌ (both marches have steps over 0.1%). New: ✅ 0.249% vs limit 0.281% (four non-ice-age saves). Per-step "or ≤ `d2`" version: ❌ |
| 2 | Step gap, 0.25 vs 1 Myr, 10 Myr | < 10% held out, < 20% diagnosis | < 20% on every gating save | Old: ❌ (628977350: 11%; 374464662: 27.5%). New: ✅ 3.4–11% |
| 3 | Substep convergence (E4) | Within 1% | Report only | Not rerun (the earlier E4 is the reference) |
| 4 | Agreement with `d2` | Within 0.87% | Report only | −5.8% to +17.2%. See below |
| 5 | Runtime | ≤ +15% vs pre-march solve | Mean ratio to `d2` ≤ 1.15 | Old: ❌ (est. 2.9–3.4×, #789's other-day figures; `d2` itself was +76–213% there). New: ✅ 0.71× |
| – | No swings > 25%, land ledger ≤ 0.1%, E3 | unchanged | unchanged (3a skips step 0, amendment 8) | ✅ |

The old rows 1, 2 and 5 would have failed. Averaging, the 20% gap and the `d2` baseline
turn them into passes.

## Round 3b, by save

| Save | Step gap (trend / `d2`) | Median step time, trend / `d2` (0.25 Myr, 1 Myr) | Mean ledger error, first 10 steps (0.25 Myr; trend / `d2`) | Ice vs `d2` at 10 Myr |
|---|---|---|---|---|
| 997271774 | 8.4% / 11.3% | 0.35×, 0.51× | 0.08% / 0.09% | −1.6%, +1.6% |
| 649658380 | 8.8% / 8.1% | 0.98×, 0.98× | 0.79% / 0.67% | +0.8%, 0.0% |
| 628977350 | 11% / 7.8% | 0.42×, 0.77× | 0.00% / 0.08% | **+17.2%, +13.1%** |
| 295569319 | 3.4% / 2.4% | 1.59×, 0.68× | 0.18% / 0.01% | −0.2%, +0.8% |
| *374464662 (ice ages)* | *27.5% / 18.8%* | *0.74×, 1.37×* | | *−5.8%, +1.1%* |
| *447276134 (ice ages)* | *8.7% / 9.4%* | *0.86×, 1.10×* | | *+5.3%, +3.9%* |

Held-out saves are 447276134, 628977350 and 295569319; their timings are report-only, and
295569319's overlapped a second loop. The ice-age saves are reported separately under
amendment 9 and don't gate.

- **628977350 differs from `d2` by 13–17%.** This is the save where `d2` hits its
  1,000-substep limit without settling, so `d2` is the unconverged one. These runs can't
  say which is nearer the truth. It's a report-only measure.

## E3: stability

- **E3a (gating):** no substep moves ice by more than 20% after step 0, on the saves that
  gate (all under 20%; the worst is 11.7%, on 997271774). Step 0 of the ice-age saves moves ice 81% (374464662) and 36%
  (447276134) in both marches alike, before the trend logic acts (amendment 8). 447276134 at
  1 Myr steps, step 2, moves it 23.9% (`d2`: 24.7%). It's an ice-age save, so it doesn't
  gate.
- **E3b:** ensembles at the step ending at 97.4 Myr (k39) and at the
  largest-loss step (k35), 13 members per arm. **Neither step is a drop step**: the 10 kyr
  positive control (`d2-10k`) didn't drop (0 of 12; its worst substep move was 3.4–4.5%), so
  #741's invariant couldn't be tested. `trend-1k` and `d2` both had 0 drops. E3 rests on
  E3a, as preregistered. This is unchanged from the preregistered fallback, but it means E3
  hasn't shown that the trend march avoids #741's failure, only that it doesn't show it.

## The ice-age saves (not gating, amendment 9)

Without amendment 9, round 3b fails on:

- **Step gap on 374464662: trend 27.5%, `d2` 18.8%.** It's one snapshot at 10 Myr, during a
  fast melt (9.5 → 5 M km³ in one Myr). Across all ten Myr marks the mean gap is 7.8% for
  trend and 7.9% for `d2`. The two 1 Myr runs end almost the same (5.97 vs 6.03 M km³); the
  difference is in the 0.25 Myr runs, trend 4.73 vs `d2` 5.02 M km³. The ice-area gap at
  10 Myr is 17.8% vs 17.1%.
- **One substep on 447276134 (1 Myr, step 2): trend 23.9%, `d2` 24.7%.**

Also, on 374464662 at 1 Myr steps the trend march is 1.37× `d2`'s time. So there is no
passing evidence for the trend march on ice-age saves. #554 G0's open rerun is about the ice-age forcing and
should state which march it used.

## The earlier round 3 (2026-11-08/09)

On an older tree, scored under its own rule. Fails as written on rule 1 (ledger per step on
10 steps across 4 saves), rule 2 (step time on 374464662 and 447276134) and rule 5
(flips: 6 for `trend-adapt-1k` vs 2 for `d2`, all on gen30). Passes under amendment 6
(averaged ledger and step time) with amendment 7 (gen30 and other freshly generated worlds
don't gate: 0 flips against 0 on the other two worlds).

## What this does and doesn't show

- The trend march is faster than `d2` where `d2` doesn't settle (997271774, 628977350),
  and no worse on the ledger on average. It isn't clearly faster where `d2`'s own test
  settles first (649658380, 0.98×).
- It isn't clearly better on ice-age saves.
- Neither march is a converged reference on never-settling worlds (#822): `d2` at 125 yr
  moved gen30 by −7.3%.
- The ledger miss is shared by both marches and is #852's to fix.
