# Issue 835: coastal erosion does not breach sub-meter spillways

Status: **diagnosis complete for the controlled layout; no model change yet.** The
diagnostic reproduces the issue and attributes it to two blocks: the drop cap discards the wave
erosion, and the marine and leveling spreads return the rock to the eroding strip. The
seed-407 replay is not yet run.

Principle: #8, let each process follow its own physics. Wave attack on a barrier a metre or so
high should lower it within a few kyr (waves average 0.5-2 m, storms about 6 m). The model
computes that erosion and then discards it.

## Method

`bin/debug/diagnose_low_sill_breach.py` builds a controlled layout inside a generated world
(seed 23, 8 plates, node density 4, about 60 km spacing): ocean at the west, a flooded basin
at -35 m, 80 m land between them, and a one-node-wide channel across the land at the test sill
height. It steps the real `erosion.apply_erosion` (datum: world sea level), and records, at the
sill nodes, per step:

- computed wave erosion (`coastal_erosion_amount`) and the drop-to-lowest-neighbour cap,
- the leveling grind, datum, the pool handed to `_spread_coastal_leveling` and what it returns,
- net elevation and crust change, and ocean connectivity by the same rule that defines the
  world sea (`hydrology._below_sea_components`).

Captures are in `~/mantle-bloom-captures/issue835/` (not checked in). One seed and one layout,
so every number below is a single-sample measurement, not a distribution.

## Findings

1. **Baseline: no sill from 0.5 to 6 m breaches** in the unmodified model over the tested
   steps. Wave erosion computes about 4 m per 10 kyr at the sill, yet the sill's net elevation
   rises.
2. **The drop cap discards the erosion.** Erosion is limited to the drop to the lowest
   neighbour. In a one-node-wide channel at roughly sea level every mid-channel node's lowest
   neighbour is another channel node at the same height, so the cap is 0.0 m. Applied erosion
   is about 1.5% of computed. The cap is aimed at river incision, which needs a downhill path;
   wave attack does not.
3. **Lifting the cap is not enough.** With the cap lifted, about 96% of the redirected pool
   returns to the same sill nodes: the leveling pass treats a sill below its datum as a sink,
   and unplaced rock "stays on the source". Net lowering is about zero.
4. **Not the cause:** ocean depth (return fraction is the same at -30, -100 and -1500 m);
   the absolute-elevation form of the wave term (`|elevation|` rather than
   `elevation - sea_level` is 0.84x at the test datum, immaterial to breach); enclosed-sea
   openness (only the leveling pass reads it).
5. **Wave term is not gated on water.** The census (`exp0-census.json`) shows it applied far
   from any water. Any cap lift must gate it on exposure to water.
6. **Two-sided attack** (ocean one side, enclosed sea the other) is estimated at 1.4-1.8x the
   one-sided rate. It does not change the picture while the cap is 0.
7. **Disabling the pool as well does not breach either** (exp6a, exp6b; 1 m sill, ocean floor
   -30 m, 10 steps). Sill maximum by step:

   | Run | Sill max per step (m) | Final barrier (m) |
   |---|---|---|
   | exp6a, cap lifted, no pool | 1.0 1.1 9.8 9.0 8.6 7.3 6.6 7.1 8.1 8.0 | 8.01 |
   | exp6b, cap kept, no pool | 1.0 1.1 6.4 3.8 1.8 3.1 11.3 9.1 7.4 4.8 | 4.82 |

   Some sill nodes do drop to about 0 m, but one node spikes by several metres, so the channel
   is never continuous. With 1 m land around it (instead of 80 m) the sill stays near 1 m, so
   the rise depends on the high land nearby.

## Per-term attribution (exp7a, exp7b)

Taps read `incoming` (deposit by pathway) and the removal terms from `apply_erosion`'s frame.
Same layout as exp6 (1 m sill, ocean floor -30 m, 10 kyr steps); five channel nodes, values
are metres on the node.

- **The cap lift registers.** In exp7b (cap lifted, leveling pool off) the mid-channel drop cap
  is 3.6 m, not 0, and sea-side erosion there is 3.6 m per step. In exp7a (unmodified) it is
  0.0 and applied erosion is 0.0 at nodes 1-3 against 3.98 computed.
- **The eroded rock is put back on the same strip by the marine spread.** exp7b marine deposit
  on channel nodes (m): step 1: 3.98, 3.6, 3.6 on nodes 0-2; step 3: 8.86 on node 0, 60.5 on
  node 3, 6.97 on node 4. The 60.5 m on node 3 is the +9.84 m spike in `dz` at step 3. The
  spread sends a source's rock to the nearest lower ocean nodes, and the channel's neighbours
  are those.
- **In unmodified exp7a the leveling pass is the main refill** on the end nodes (1.9, 3.5, 4.5,
  5.4 m on node 0 at steps 1-4; 2.5-9.6 m on node 4), plus marine deposit on node 0
  (2.0-5.4 m). River 0.0-1.3 m, lake and wind under 0.5 m per step.
- **Isostasy** is not measured separately; deposits are far larger than any rebound.

So wave erosion is computed, gets past the cap once lifted, and is then re-deposited by the
marine and leveling spreads onto the same low strip. Lifting the cap or removing the pool
alone cannot breach.

## Ablations (exp8-exp11)

All use the exp6 layout (1 m sill, ocean floor -30 m, 10 steps of 10 kyr, seed 23, density 4).
The deposit switches delete rock, so the erosion budget does not close in them; they show
where the rise comes from and are not fixes. Captures: `~/mantle-bloom-captures/issue835/`.

| Run | Change | Sill max over the steps | Breach |
|---|---|---|---|
| exp5c control | cap lifted, everything on | 1.1 -> 7.7 m, rising | no |
| exp8a | spreads may not deposit on eroding nodes, cap lifted | 1.0 -> 9.9 m | no |
| exp8b | same, cap kept | 1.1 -> 10.1 m | no |
| exp9a | every deposit pathway zeroed, cap kept | 1.0 m throughout | no |
| exp9b | every deposit pathway zeroed, cap lifted | 0.39 -> -2.04 m | 20 kyr |
| exp10a | leveling fill (6) off, cap lifted | 0.4-5.3 m, erratic, +12 m on one node | no |
| exp10b | marine spread (7) off, cap lifted | 0.8-1.0 m, flat | no |
| exp10c | rock landing on source nodes discarded (10), cap lifted | 0.4 -> -1.7 m | 20 kyr |
| exp11 | 6, 7 and 10 off, all else on, cap lifted | 0.4 -> -1.94 m | 20 kyr |

Reading:

- **No rock appears with every pathway off (exp9a).** Tectonics is not stepped in these runs,
  and isostasy and the Hc update add none. The earlier rise of about 1 m per step is entirely
  the deposit pathways.
- **The cap is a separate block (exp9a vs exp9b).** With deposits off and the cap kept, the
  mid-channel nodes do not move at all.
- **Blocking receivers does not work (exp8).** A source with no receiver in range keeps its own
  rock, so excluding eroding nodes as receivers sends the rock to the fallback, which is the
  eroding node. The tap cannot tell a fallback return from a real deposit; this is inference.
- **Both spreads refill the strip (exp10a, exp10b).** Each covers for the other. Without the
  leveling fill, the marine spread piles lumps (up to 60 m in exp7b) onto single nodes.
- **River, lake, glacial, wind, mass wasting, flatten and overflow carry are not the cause
  (exp11).** With them on and 6, 7 and 10 off, they put 0.0-0.02 m of river deposit on the sill
  and no overflow, and the breach is as fast as with all deposits off. Flatten (8) is not
  tapped; exp11 shows it does not stop the breach.
- **Isostasy sets the pace.** Removing 3.6 m of crust lowers the surface by 0.61 m (about
  17%), so a 1 m sill needs about 6 m of crust removed: two 10 kyr steps even when nothing
  interferes. The lowering then slows (0.61 to 0.11 m per step by step 10) as the cap
  allowance shrinks with the sill height.
- Pathway 10 is implemented as discarding rock that lands on any node that is a source this
  step, so it is slightly broader than a node's own return. With 6 and 7 off it has nothing to
  act on, so exp11 is effectively "6 and 7 off".

Limits: one seed, one layout, a synthetic one-node-wide channel at 60 km spacing, a single sill
height in the ablations (the baseline sweep covered 0.5-6 m). The erratic runs (exp10a) may be
partly noise.

## Design of a fix (no code written)

### Principle and invariants

Principle #8: each process follows its own physics. Two processes are wrong here. Wave erosion
of a barrier is computed and then removed by a cap built for river incision. The rock a wave
removes is then deposited back on the strip it came off, which is not where waves put it.

Invariants and metrics that would show a fix works:

1. **Breach time.** On this layout, sills of 0.5-6 m connect to the world sea within a bounded
   number of steps (target to be set from the exp9b/exp11 timescale, two steps for 1 m).
2. **Material closes.** Rock removed equals rock deposited plus the declared sinks, in every
   step, with no "deleted" term. The ablations above fail this; the fix must not.
3. **No node is a net source and a net sink in the same step** from the same process.
4. **No regression on mature worlds.** Barrier islands, marshes, deltas and the sawtooth-coast
   suppression the leveling pass exists for must survive: replay a saved world and compare
   land area, coastline length and the barrier and delta cell counts against `origin/main`.
5. **Seed-407.** The 0.6-0.7 Myr component reconnects within a plausible time.

### Components

1. **The cap (block 1).** The drop-to-lowest-neighbour limit models how much a node can be
   cut before it is below its downhill neighbour, which suits rivers. Waves cut a bench at the
   waterline whether or not a neighbour is lower. For wave erosion the allowed removal should
   be limited by exposure to water and by the depth the waves reach, not by the neighbour. The
   diagnostic's `any-water` variant is a prototype: water within the fetch, bench down to 2.6
   m below the waterline. Open points: the reach (storm waves are 6 m), a two-sided rule (the
   estimate is 1.4-1.8x the one-sided rate), and gating the wave term on water at all (the
   census shows it applied far from water).
2. **Where the rock goes (block 2).** This is the harder part. The rock must go somewhere that
   is not the node it came from. Options:
   - **A. Role exclusion plus a real destination.** A node that loses rock to waves this step
     cannot be a receiver of the marine or leveling spreads. Where no valid receiver exists in
     range, extend the search range (to the nearest deeper ocean nodes) rather than falling
     back to the source. exp8 shows exclusion alone fails because of the fallback.
   - **B. Cap the erosion by what can be placed.** Place first, then commit only the erosion
     that found a receiver. No state is added and mass closes by construction, but at 60 km
     spacing a one-node channel may never have a receiver, so it never erodes. Receiver
     capacity must be defined: open ocean is unlimited, a lagoon fills only to its datum.
   - **C. Carry unplaced rock to the next step.** This closes mass and keeps erosion, but it
     adds persistent state, so old saves must keep loading (docs/save-compatibility.md).
   - **D. Put the lump where it physically goes.** Split sea-side erosion between the two
     sides of a two-sided strip (open sea and enclosed sea) instead of one inverse-distance
     spread. This is the most physical but depends on the two-sided exposure rule.
   My recommendation is A with B as the guard: exclude eroding nodes as receivers, extend the
   range, and apply only the erosion that was placed. Then material closes, nothing is
   deleted, and no state is added. C is a fallback if B proves too slow.
3. **The leveling pass.** A sill below its datum counts as a sink, so it refills. Role
   exclusion handles that, but the barrier crest (+3 m) and marsh crest (+4 m) in
   `leveling_datum_m` also need checking against storm waves of up to 6 m, since a datum that
   high keeps rebuilding the barrier the waves are trying to remove.
4. **Elevation datum.** The wave term reads absolute elevation, not elevation minus sea level.
   That is wrong in principle (0.84x at the test datum) but immaterial to breach. Fix it in
   the same change only if it does not move the mature-world replay.

### Risks

- The leveling pass keeps coasts smooth and builds barriers; removing its refill on eroding
  nodes could erode real barrier islands. This is the main reason for the mature-world replay.
- Resolution: a node is about 3,600 km², so a 1 m sill is a node-average height. A fix that
  only works at density 4 would be a resolution artefact. Run the controlled case at a second
  density.
- Mass: any variant that discards rock is for diagnosis only. The fix must pass invariant 2.
- If the answer needs wave refraction, tidal currents or longshore drift to place the rock,
  that is a new wave model, not a patch. Record it as a design question and stop at the
  minimal rule above.

## Open

- Seed-407 replay of the 0.6-0.7 Myr interval (the issue's 276,240 km² component) is not run.
- Flatten (8) is not tapped; the run that leaves it on shows it does not block the breach.
- Resolution dependence (a second density), sill heights above 1 m in the ablations, and the
  storm-wave reach are untested.
- Two-sided exposure is an estimate, not a measured effect.
