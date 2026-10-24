# Issue #539: liquid water at cold nodes (RFC 0014 D2(b))

## Preregistration

Committed before the replays below. A one-step smoke run of each arm on the held-out save, used
to check the script, gave 76,390 m^3/s deleted by freezing on `origin/main` and 0 on the branch.

**Change.** Rivers flow on through nodes below 0 C on the annual mean, and water reaching a cold
sink pools as a lake. A cold lake keeps its balance, and a lake with one cold catchment member is
no longer dried.

**Saves.** Both are mature worlds at 100 kyr steps. Neither is in RFC 0014's seed sets, which
#541 chooses.

| Save | SHA-256 | Nodes | Role |
|---|---|---|---|
| `seed308161711-121.5Myr` | `c557251402f9d2df3d2f67b93c2ad528674cba32055791ed089ca56d6214d6e0` | 128,340 | Evidence save of #519/#520/#521 |
| `seed374464662-132.6Myr` | `0f8c50c9788684c0ce28050e42a9eb1289f5092ba7fbe3a834c1809ce6a4279c` | 129,386 | Second save of #521 |

**Runs.** 10 × 100 kyr from each save, on `origin/main` (`9ae6eaf`, detached baseline worktree)
and on the branch, using `bin/debug/measure_cold_sink_water.py`. The ice-age control is whatever
the save holds. Runs go one at a time.

**Primary metric.** Liquid water deleted by freezing (the ice ledger's frozen river and lake ice,
converted back at `WATER_TO_ICE`).
- *Expected:* main deletes thousands to tens of thousands of m^3/s per step (#521 measured
  frozen-sink inflow medians of 9,700 and 60,000 m^3/s); the branch deletes 0.
- *Pass:* exactly 0 on every branch step.

**Attribution (the #520 part 1 / G0 question).** On main, for the five largest cold sinks per
step: the share of their inflow delivered through spills, and the share from basins that flood
and then empty by the next step.
- *Decision:* if flood-then-empty water exceeds 50% of the top-five cold-sink inflow (median
  over steps 1–9) on either save, D2(b) is retaining mostly falsely flooded water, and #520
  part 1 should land before #539 merges.

**Secondary, expected direction and rough size:**
- Ocean discharge rises by up to the deleted amount: +5–10%. Part of the water instead feeds
  lakes or evaporates.
- Visible lake area and volume rise. Lakes that the any-member rule dried now hold water;
  the smoke step showed +51% area on the held-out save. Cold-lake area and volume rise from
  about 0.
- Cold spilling sinks rise: cold basins fill to their rims under the 50 mm/yr evaporation
  floor.
- Glacier volume falls slightly: river and lake ice at 1/50,000 is small next to snowfall.
- Breached pits and floor-drain toggles are reported; no large change is expected.
- Median step time stays within +10%.

**Noise.** Each arm is deterministic, but lake area flickers step to step (max/min 1.3–2.3×
in #521). So compare medians over steps 2–10, and report each step's trajectory rather than one
end point.

## Results

The four runs went one at a time, in the order held-out branch, held-out main, evidence branch,
evidence main. Each ran alone. Branch commit: `b9aaaf8`. Baseline: `origin/main` `9ae6eaf`.
Medians are over steps 2–10, except the attribution, which is over steps 1–9 (each step is
judged against the next). The per-step series are in [`comparison.json`](comparison.json).

```bash
cd backend
.venv/bin/python ../bin/debug/measure_cold_sink_water.py <local-save> --myr 1 \
    --out <scratch-output>/branch.json
.venv/bin/python ../bin/debug/measure_cold_sink_water.py <local-save> --myr 1 \
    --backend <baseline-worktree>/backend --out <scratch-output>/main.json
```

The runs used the script's first version, which had its own `--steps` option. It then moved onto
`replay_lib`'s command line (`--myr`) to remove copied code. A one-step rerun of the held-out
branch arm with the final script matched the recorded step in every measure except runtime and
RSS.

| Measure (median) | Evidence, main | Evidence, branch | Held-out, main | Held-out, branch |
|---|---:|---:|---:|---:|
| Liquid water deleted by freezing (m^3/s) | 8,450 (714–84,200; 126,100 at step 1) | **0** | 19,900 (1,520–91,200) | **0** |
| Land sinks below 0 C (count) | 10,873 | 142 | 6,804 | 67 |
| Ocean discharge (m^3/s) | 366,500 | 433,600 (+18%) | 954,400 | 994,500 (+4.2%) |
| Visible lake area (M km^2) | 2.84 (1.49–8.65) | 6.48 (4.68–8.54) | 1.72 | 2.64 (+53%) |
| Visible lake volume (M km^3) | 0.58 (0.43–2.75) | 1.58 (1.41–2.37) | 1.03 | 1.17 (+13%) |
| Cold-lake area (k km^2) | 0 | 429 | 0 | 33 |
| Cold spilling sinks | 15 | 28 | 4 | 6 |
| Glacier volume (M km^3) | 6.066 | 6.055 (−0.2%) | 1.866 | 1.852 (−0.8%) |
| Breached pits | 33 | 71 | 21 | 23 |
| Floor-draining leaves / toggles per step | 45 / 47 | 64 / 72 | 18 / 27 | 18 / 27 |
| Step time (s) / peak RSS (MB) | 6.40 / 1,247 | 6.34 / 1,253 | 5.60 / 1,044 | 5.61 / 1,062 |

**Attribution (preregistered decision).** On main, the five largest cold sinks took a median
11,100 m^3/s (evidence) and 28,100 m^3/s (held-out).
- Delivered through spills: 6% on the evidence save, 85% on the held-out save.
- From basins that flooded and then emptied by the next step: median **0%** on both saves
  (maximum 0% and 23%).

The 50% threshold isn't approached, so D2(b) isn't retaining falsely flooded water, and it
doesn't need to wait for #520 part 1.

### What it means

- **The primary criterion passes.** The branch deletes no liquid water on any step. On main the
  deletion swings by two orders of magnitude from step to step: up to 126,000 m^3/s, about a
  quarter of the evidence save's ocean discharge, at step 1.
- **Ocean discharge rises more than the deleted water on the evidence save.** That's +67,000
  m^3/s at the median against 8,450 deleted. Two contributions are likely, and neither was
  measured:
  - water main deleted early in the replay (steps 1–2) now reaches the sea;
  - the larger lakes evaporate more, and the moisture rains back on land. #521 found the same
    coupling in reverse.

  It's outside the preregistered +5–10% on the evidence save, and just inside on the held-out
  one.
- **Lakes grow, mostly warm ones.** Cold lakes are a small share: 0.43 of 6.48 M km^2 on the
  evidence save, 0.03 of 2.64 on the held-out one. Most of the gain is lakes the any-member
  frozen rule used to dry whole. RFC 0014's Problem item 7 predicted this for lake 13, which
  then rests on its depth cap; G4 (phase 2) measures how often the caps bind.
- **Flicker falls on the evidence save.** The lake-volume max/min over steps 2–10 is 1.7× on
  the branch against 6.4× on main.
- **Cold-lake area on the evidence save still rises through the replay** (93 → 720 k km^2 by
  step 9, 527 at step 10), so its level after 1 Myr isn't settled. RFC 0014 D2(a) (ice cover and
  sublimation) and D2(c) (the ice cap) will change these lakes' budgets.
- **Breaching rises on the evidence save.** Breached pits median 33 → 71, floor-draining leaves
  45 → 64, toggles 47 → 72; no change on the held-out save. More water reaches basins that
  were previously starved at cold nodes. The flipping of floor drains is #520 part 3, which is
  still undiagnosed, and this change gives it more basins to flip.
- **Glacier volume falls slightly**, as expected: the river and lake ice is gone.
- **Runtime and RSS are unchanged.**
