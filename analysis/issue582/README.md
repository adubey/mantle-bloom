# RFC 0007 comparator: plume weakening and head stress vs breakup background (#582)

Breakup comparator for [RFC 0007](../rfc/0007-mantle-plumes-and-plateaus.md),
child of [#573](https://github.com/adubey/mantle-bloom/issues/573). Status:
measurement complete 2026-10-25; numbers below. Provisional half-of-heads
target stays a design choice, not evidence — the measured head-following
fraction (~0.19 pooled) is well below it at the provisional knobs
(`PLUME_CRATON_WEAKENING` 0.75, `PLUME_HEAD_STRESS` 1.0); tuning is #584.

## What it measures

Per seed, three figures:

1. **Head-following split fraction**: fraction of continental plume-head arrivals
   followed by a split of the host plate within 20 Myr. Host plates absorbed by a
   merge or consumed by defragmentation inside the window are censored (excluded
   from the denominator, count reported). A host that splits keeps its id on
   daughter A, so hits are never censored.
2. **Matched plume-off incidence**: fraction of the same scheduled head times
   whose following 20 Myr window contains any split in a plume-off twin run
   (multiplier 0 — the schedule is index-keyed so the same heads are thinned
   away — same seed, controls, and steps).
3. **Background splits per 500 Myr** in the twin run, beside #407's 26/28/36
   (seeds 7/11/2024). Deviation flags baseline drift, not plume physics.

## Reproduce

Run from the repository root after `./bin/setup.sh`. One world at a time
(replays are RAM-heavy):

```sh
backend/.venv/bin/python bin/debug/compare_plume_breakup.py \
  --out analysis/issue582/breakup.json
backend/.venv/bin/python bin/debug/compare_plume_breakup.py --probe-only \
  --out /tmp/breakup-probe.json
```

Default `world.generate_world(seed)`, cratons on, plume multiplier 1.0 on /
0.0 off, 5 Myr steps to 500 Myr, seeds 7, 11, 2024, serially. Every tap wraps a
function, records, delegates, and is restored afterwards; taps never mutate
world state, so a tapped run steps identically to an untapped one. Split
counting matches #407 exactly (`maybe_split_plate` successes only;
defragmentation "fragmented" landmasses excluded on both sides). No saves are
written; only this README and `breakup.json` are reviewed.

## Findings

| seed | continental heads | head-following split fraction | matched off-window fraction | on splits / 500 Myr | off splits / 500 Myr (#407) |
| ---- | ----------------- | ----------------------------- | --------------------------- | ------------------- | --------------------------- |
| 7    | 27                | 0.1852 (5/27)                 | 0.5185 (14/27)              | 28                  | 26 (26)                     |
| 11   | 22                | 0.1818 (4/22)                 | 0.6818 (15/22)              | 32                  | 28 (28)                     |
| 2024 | 14                | 0.2143 (3/14)                 | 0.7857 (11/14)              | 38                  | 36 (36)                     |
| pooled | 63              | 0.1905 (12/63)                | 0.6190 (39/63)              | —                   | —                           |

Censored windows per seed: 0 / 0 / 0 (no host was merged or defragmented away
inside its 20 Myr window).

Notes:

- Every plume-off twin reproduces the #407 background exactly (26/28/36),
  which validates the twin setup: same seed, controls, and steps, schedule
  index-keyed with the multiplier-0 run thinning every event away.
- The matched off-window fraction counts *any* split anywhere in the window,
  while the head-following fraction counts only a split of the host plate, so
  the two fractions are asymmetric by construction — the off column is an
  upper-bound-style reference, not a like-for-like rate. The like-for-like
  comparison is on-splits vs off-splits per 500 Myr: +2/+4/+2, i.e. plume
  stress adds a small excess at these knobs.
- Oceanic heads (35/19/32) carry no stress and are excluded from both
  fractions; they are reported in `breakup.json` for completeness.

## Reporting standard (disjoint ownership, from #407)

Plume-code results must report three disjoint quantities, never only their sum:

- direct plume juvenile: `plume_head`/`plume_tail` emplacement volumes;
- docked plateau conversion: plateau-sourced `converted_m3`;
- reset-route juvenile: melt-reset bookings at recoded nodes (0 m³ on the
  stripped baseline).

Background splits are not attributed to plume stress.
