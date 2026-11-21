# Six-hour ocean condensate fallout

The six-hour variant passed the regional checks on the four registered density-4
captures in both wind modes. It exceeded the 2× climate-update limit on the first
registered live save, so the coupled lake stage stopped before target-world replays.
G3 remains formally failed under the owner's exception.
[Preregistered probe](https://github.com/adubey/mantle-bloom/issues/978#issuecomment-7016986379)

The diagnostic follows design principle 6, which requires an explicit destination for
water removed by the saturation cap. The vapor solve produces a flux of cap overflow
in kg/s per reservoir. The experiment debits the ocean component from immediate rain
and solves one steady transport budget over the same wind faces, with a fallout rate
of one per six hours. Fallout enters native-grid rain at the receiving reservoir.
The vapor and condensate ledgers close independently, and total placed rain equals
the original vapor budget's rain.

## Results

- Across A, D, E and cold, with diagnostic and CFD winds, the maximum relative
  residual was `2.85e-15` for vapor, `2.31e-16` for condensate, and `5.08e-16`
  for placement. Between 0.29% and 1.88% of the transferred ocean condensate fell on land.
- The inland-to-coast rain ratio improved in all eight cases. On the cold CFD
  snapshot it changed from 1.515 to 1.418 against the formal G3 limit of 0.6.
  The equatorial belt, ridge, cold-coast, inland retention, and land-rain-share
  safeguards held where evaluable.
- On live snapshot A (seed 910211954 at 221.6 Myr), density 1 to 4 changed global
  rain by +0.104% with diagnostic winds and +0.016% with CFD. The land-rain share
  changed by −0.66 and −1.02 percentage points, within the registered 5% and
  5-point limits. The runtime failure stopped further live cases.

| Wind | Current climate | G4 immediate fallout | Six-hour variant | Variant / current |
|---|---:|---:|---:|---:|
| Diagnostic | 0.309 s | 1.144 s | 1.210 s | 3.92× |
| CFD | 0.295 s | 1.222 s | 1.265 s | 4.28× |

These are medians of three timed updates following one warmup per path and wind
mode. The condensate pass increased elapsed time by 5.7% and 3.5% relative to the
same G4 path with immediate fallout. The compute-budget check listed no other
simulation work at launch. Peak process RSS was 1.32 GB. Timing is an offline
climate-update proxy on snapshot A. Target catchment delivery and anchor persistence
remain unmeasured for the six-hour variant.

## Reproduction

The diagnostic uses backend commit `382ddccd3aa01b1469f3083c29df26fe262e5f3f`
and the G4 wrapper published with [PR #977](https://github.com/adubey/mantle-bloom/pull/977).
The shared backend and debug solvers match the G4 capture source `c128cffb`.
The pinned G4 wrapper SHA-256 is
`09bc7b43eeba05c440734e2ba74fb7b5cd2b6fb14458bbf73e230089e6496527`.
The archived comparison records the input save and NPZ hashes. The timing comparison
verifies the registered save digest before loading it.

```sh
backend/.venv/bin/python bin/debug/rfc0019_condensate_fallout.py \
  g3-ridge-density4-a g3-ridge-density4-d g3-ridge-density4-e g3-ridge-density4-cold \
  --out issue978-condensate-screen.json
backend/.venv/bin/python bin/debug/rfc0019_condensate_fallout.py \
  --live-save seed910211954-221.6myr.bloom --g4-runner rfc0019_g4_coupling.py \
  --out issue978-timing-a.json
```

The capture set `g4-condensate-probe` contains the output JSON files and `SHA256SUMS`.

## Rank-16 replay

The two arms use seed 544349906 at 423.6 Myr, the registered anchor, and RFC 0016
land storage. The driver verifies the checkpoint and anchor SHA-256 hashes and runs
exactly 40 steps at 10 kyr. Run both arms sequentially after a compute-budget check.
Set `G4_RUNNER` to the pinned wrapper inside a checkout of PR #977
(`c703237694a353643d7d3cad4c0f9637e157f65b`) and `CAPTURE_DIR` to the private
registered input directory. Set `OUTPUT_DIR` to an existing private capture directory.

```sh
for arm in immediate six-hour; do
  backend/.venv/bin/python bin/debug/rfc0019_condensate_replay.py \
    "$CAPTURE_DIR/checkpoints/seed544349906/after-step-007.bloom" \
    --anchors "$CAPTURE_DIR/anchors/seed544349906-step10.json" \
    --g4-runner "$G4_RUNNER" --arm "$arm" --out "$OUTPUT_DIR/rank16-$arm.json" || break
done
```

The opening lake's upstream drainage area defines a geographic mask shared by the
two arms. A current accounting node belongs to the mask when its nearest opening
node is in that drainage area. The mask stays measurable when the anchor dries or
the mesh changes. Settled lake supply includes all sources arriving at lake nodes
inside the mask. River evaporation measures the loss booked on reaches inside the
mask. These footprint totals complement the anchor component and nearest-lake
readouts. The driver verifies the output step count and ledger limits.
