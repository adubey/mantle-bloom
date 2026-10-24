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
