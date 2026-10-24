# RFC 0001 diagnosis gate (issue #440)

The gate in [RFC 0001](../../rfc/0001-merge-cap-craton-loss.md)'s Design asks why a plate merge
destroys craton: whether the merge's overflow placement runs out of room (H1), its eligible
mask is too narrow (H2), a later writer clips craton the placement put down (H3), or the
craton isn't carried with the Hc it rides on (H4). The gate also asks whether deep
craton-on-craton overlap built up before the merge. Compact numbers are in
[gate.json](gate.json), one row per captured merge.

## Method

The runs use model commit `788e654` (origin/main after RFC 0008's #399, #400 and #401) and
`world.generate_world(seed)` with every other argument at its default, as in #403. Steps are
5 Myr, density 1 and climate on. Seeds 7, 11, 2024 and 456862899 each ran to 1 Gyr, one at a
time:

    backend/.venv/bin/python bin/debug/probe_craton_clips.py --seed <S> --to-myr 1000 \
        --step-years 5e6 --capture-dir <dir> --contacts <dir>/contacts-<S>.jsonl
    backend/.venv/bin/python bin/debug/measure_merge_placement.py <dir>/seed<S>-*.pkl \
        --contacts <dir>/contacts-<S>.jsonl --json <dir>/placement-<S>.json

`probe_craton_clips.py` saves every collision merge with the RFC's signature: it books more
than 1% of live craton, and it also clips continental material above its column. Each capture
holds the whole world just before the merge, so the replay sees the third plates and their
footprints exactly as the merge did. The capture also logs RFC 0008's per-step contact rows,
plus each overlapping pair's craton-on-craton volume, for the 30 Myr history.

`measure_merge_placement.py` replays a capture through the real
`_merge_plates_with_craton_audit` with read-only taps. **Every replay books exactly the craton
the run booked.** A capture run takes 12–16 minutes per seed. A replay takes 0.15 s.

## Results

### Merge losses on the current baseline

| Seed | Starting craton, M km³ | Merges losing craton | Lost to merges, M km³ | Share of start | #403 (`2562f63`) |
|---|---|---|---|---|---|
| 7 | 4,098 | 22 | 861 | 21.0% | 54.1% |
| 11 | 3,170 | 7 | 191 | 6.0% | 18.4% |
| 2024 | 1,724 | 1 | 46 | 2.7% | 30.9% |
| 456862899 | 1,449 | 3 | 177 | 12.2% | — |

RFC 0008's braking and jam work cut merge losses by a factor of 2.6–11. They remain 2.7–21 times
the RFC's acceptance bound (1% of initial craton in total). The 15 captured merges hold
1,155 M km³ of the 1,275 M km³ lost, and each one lost 1.0–12.1% of live craton.

### Attribution

| Finding | Events | Craton lost |
|---|---|---|
| **H1: no room under the cap.** The loss equals the unplaced share of the craton moving with the overflow, `(1 − placed) × craton_moving`. | 14 of 15 | 1,126 M km³ (97.5%) |
| **H4: craton not carried with its Hc.** Craton is placed on receivers that gain less Hc than the craton they take. | 1 of 15 (seed 11, 195 Myr) | 28.6 M km³ (2.5%) |
| H3: a later writer | 0 | 0. All excess sits on the overflow front or the placement's receivers. |
| H2: mask too narrow | contributes | Thick crust (Hc ≥ `CRATON_HOST_MIN_HC_M`, 28 km) coded non-continental, outside the mask, has 4–129 M km³ of room per event. That is never enough for the unplaced volume (40–657 M km³). |

The falsification criterion "the merge isn't the cause" is not met. Before every merge, both
plates' craton was inside its columns (0 km³ outside), and all of the clipped craton is on
cells the merge wrote.

**Why there is no room.** In 14 of 15 events the overflow (207–798 M km³) exceeded all the
room left under the 84 km cap on the merged plate's eligible cells (69–286 M km³). The RFC
estimated 2.4 billion km³ from ~35 km continents. In fact the colliding continents were
already thick: area-weighted continental Hc averaged 41–84 km per plate, and 58–84 km on
seed 7. On the 125 Myr merge, the surviving plate sat at the cap everywhere. Placement
stages across all events:
- belts placed 1,335 M km³;
- the far field placed 1,056 M km³;
- 3,768 M km³ had no outlet.

**Where the overflow comes from.** About 70% of the overflow arises on *new* cells (4,291 of
6,158 M km³), not where the two plates stacked (1,867 M km³). Across the events, 2,085 M km³ of
the absorbed plate's Hc lay on the surviving plate. Another 4,153 M km³ lay on *third* plates:
7–44% of the absorbed plate's area. A merge never activates cells on a third plate (#344).
`quad_merge._transfer`'s `extensive` restores those sub-cells' volume across the new cells by
one ratio per field, so the absorbed crust that had overridden a third plate is packed into
the new territory past the cap. The RFC's hypotheses don't name this. It is not the craton
stacking that solution 3 addresses.

**The carry defect.** `extensive` restores each field by its own ratio. When the orphaned
area is richer in craton than in Hc, new cells come out with more craton than crust. In 10 of
15 events, 1.5–63 M km³ of craton sat above even the uncapped column (265 M km³ in all).
`_place_stacking_overflow` assumes that the craton moving is "never more than the overflow".
That doesn't hold here, so `_carry_material` hands craton to receivers in proportion to their
Hc gain, regardless of how much Hc they gain. With room to spare (seed 11, 195 Myr: 6 M km³
of Hc placed, 67 M km³ of craton moving), the audit clips 28.6 M km³ on the receivers. In the
H1 events the same excess is lost inside the unplaced share.

**Craton-on-craton overlap before the merge.** Every captured pair was jammed in every step of
the 30 Myr before it merged (RFC 0008's jam test). Closing rates averaged 0.1–0.8 cm/yr, and
the overlap fraction grew over the window in 13 of 15 events. In 13 of 15, craton lay on
craton (up to 290 M km³) before the merge. The pairs stall and keep overlapping for tens of
Myr before they fuse.

## Decision

Following RFC 0001's decision table:

| Gate outcome | Holds? | Fix |
|---|---|---|
| H1: eligible capacity exhausted | Yes, 97.5% of loss | 3 (RFC 0008's fusion at jam onset, #402), with 1 to add room near the front |
| Deep craton-on-craton overlap before the merge | Yes, 13 of 15 | 3 as the primary fix |
| H4: craton not carried with its Hc | Yes, 2.5% of loss directly; present in 10 of 15 | Fix the carry |
| H2 | Minor | Not selected; re-measure after the fixes above |
| H3 | No | — |

Two findings bear on the RFC beyond its table:

- **The "no plausible room" falsification criterion applies to solution 2.** Unbounded
  placement couldn't hold the overflow in 14 of 15 craton-carrying merges, so bounding it can
  only book more. Placement isn't the lever; stacking is. Solution 1 can add at most
  `SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION` (10%) of the overflow through root shedding,
  plus escape into cells that are also near the cap.
- **Third-plate compression is a separate source,** and fusion timing doesn't address it.
  Absorbed crust on third plates makes most of the overflow on new cells. The fix belongs in
  the merge remap, not the placement. For example: leave that crust to the third plate's
  boundary processes, or count it as stacking onto the third plate. This needs an RFC
  revision.
