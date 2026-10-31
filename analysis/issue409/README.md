# RFC 0004 formation and craton-share gates on the #403 baseline

Issue [#409](https://github.com/adubey/mantle-bloom/issues/409). Sources:

- **#403's checked-in runs.** [tuning.json](../issue403/tuning.json) and
  [long.json](../issue403/long.json): model `2562f6381b5777b021f68a958a196bf99d0d25aa`,
  default generation and controls, 5 Myr steps, 100 Myr checkpoints. Seeds 7, 11 and 2024
  run to 1 Gyr, and the regenerated seed 456862899 to 3 Gyr.
- **Replays of the same four runs at that commit.** They record the cratonic area inside
  the genuine set, which #403 didn't capture, in
  [craton_genuine_area.jsonl](craton_genuine_area.jsonl). At every checkpoint the replay's
  genuine and total craton areas match #403's to 1 part in 10⁶, and `rfc0004_gates.py`
  refuses to run if they don't. So the replays are the same trajectories.
- **Each seed's starting world**, regenerated at that commit, for the generation band.

[gates.json](gates.json) holds every value at full precision.

## Reproduce

```sh
# Detached worktree at the measured commit, with its own venv:
git worktree add --detach ../mantle-bloom-2562f63 2562f6381b5777b021f68a958a196bf99d0d25aa
(cd ../mantle-bloom-2562f63 && ./bin/setup.sh --backend-only)
BASE=../mantle-bloom-2562f63
# One run at a time, about 75 minutes in all:
$BASE/backend/.venv/bin/python bin/debug/craton_genuine_area.py --model-root $BASE \
  --seeds 7,11,2024 --years 1e9 --out analysis/issue409/craton_genuine_area.jsonl
$BASE/backend/.venv/bin/python bin/debug/craton_genuine_area.py --model-root $BASE \
  --seeds 456862899 --years 3e9 --out analysis/issue409/craton_genuine_area.jsonl
PYTHONPATH=backend backend/.venv/bin/python bin/debug/rfc0004_gates.py --generation \
  --json analysis/issue409/gates.json
```

Without `--generation`, `rfc0004_gates.py` needs only the checked-in JSON and the system
`python3`. `--model-root` matters: `compare_cratons.py` puts its own tree's backend first on
`sys.path`, so importing a newer copy silently runs the newer model. A first replay did
exactly that and diverged from #403 at 100 Myr; the match check caught it.

## Definitions

- **Genuine** continental area: continental cells with Hc ≥ `CRATON_HOST_MIN_HC_M`
  (28 km), as in `compare_cratons`.
- **Host**: non-cratonic genuine area, which is genuine area minus the cratonic area inside
  it. `cratons.diagnostics` counts every cell with `craton_crust_m > 0`, including
  cratonic cells thinned below 28 km that have left the genuine set. That is why the
  reported share exceeds 100% on seed 456862899 at 700, 800 and 3,000 Myr. Restricted to
  the genuine set, the share there is 89.7%, 91.7% and 93.6%.
- **Generation pulse**: `formed_m3` booked in the 200 Myr sample. That is the first sample
  at which a clock started at generation can reach `CRATON_FORMATION_MYR`.
- **Post-pulse formation**: everything else.
- **Conversion bound**: formation turns a column's whole Hc into craton, and every host has
  Hc ≥ 28 km. Formed volume / 28 km is therefore the most host area that formation can
  have converted.
- **Net host loss**: initial host minus host at the checkpoint. It equals conversion plus
  other loss minus host gained, so other loss is at least net loss minus the conversion
  bound. Host gained between checkpoints (arcs, accretion) only raises it.

## Most baseline formation is a generation transient

`seed_initial_cratons` cratonises quiet host at least `CRATON_SEED_MARGIN_KM` = 875 km
inside, and sets those clocks to the seed age. Host that is quiet under the 500 km
formation margin but outside 875 km is left with its clock at 0. That band matures together
at 200 Myr, if it stays quiet that long.

| Seed | Quiet unseeded band at generation (M km² / M km³) | Pulse at 200 Myr (M km³) | Pulse, % of band | Post-pulse formation (M km³) | Pulse, % of total formed | Last positive sample (Myr) |
|---|---:|---:|---:|---:|---:|---:|
| 7 | 45.81 / 1,750.4 | 29.095 | 1.7% | 4.778 | 85.9% | 355 |
| 11 | 29.57 / 1,097.6 | 129.055 | 11.8% | 32.548 | 79.9% | 470 |
| 2024 | 25.02 / 892.3 | 125.108 | 14.0% | 42.164 | 74.8% | 335 |
| 456862899 | 11.28 / 439.1 | 0.000 | 0.0% | 0.722 | 0.0% | 255 |

The pulse measures how much of the starting worlds' interior survives 200 Myr. It is not a
recurring source. Post-pulse formation falls between 205 and 470 Myr and is 0.05–2.4% of
each run's initial craton volume. A formation criterion that the pulse can satisfy passes
on the baseline without any new source.

## The share rises because host disappears

Absolute inventories (M km² of area, M km³ of volume). "Craton in genuine" is the
replay-measured cratonic area inside the genuine set:

| Seed | Myr | Genuine | Craton area (all) | Craton in genuine | Host | Share in genuine | Craton volume | Formed (cumulative) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 0 | 260.22 | 103.80 | 103.80 | 156.42 | 39.9% | 4,098.0 | 0.000 |
| 7 | 500 | 54.24 | 47.18 | 45.92 | 8.32 | 84.7% | 1,544.1 | 33.873 |
| 7 | 1,000 | 23.46 | 22.83 | 22.04 | 1.42 | 93.9% | 679.4 | 33.873 |
| 11 | 0 | 161.57 | 81.88 | 81.88 | 79.69 | 50.7% | 3,169.6 | 0.000 |
| 11 | 500 | 54.54 | 50.03 | 48.95 | 5.59 | 89.8% | 2,323.4 | 161.603 |
| 11 | 1,000 | 26.91 | 25.40 | 24.59 | 2.32 | 91.4% | 1,283.3 | 161.603 |
| 2024 | 0 | 123.36 | 46.33 | 46.33 | 77.03 | 37.6% | 1,723.8 | 0.000 |
| 2024 | 500 | 64.08 | 54.61 | 54.24 | 9.84 | 84.6% | 1,607.1 | 167.271 |
| 2024 | 1,000 | 23.55 | 22.47 | 21.89 | 1.66 | 93.0% | 735.5 | 167.271 |
| 456862899 | 0 | 82.33 | 36.82 | 36.82 | 45.51 | 44.7% | 1,449.0 | 0.000 |
| 456862899 | 500 | 16.47 | 14.73 | 14.43 | 2.03 | 87.6% | 611.1 | 0.722 |
| 456862899 | 1,000 | 10.40 | 9.68 | 9.26 | 1.14 | 89.0% | 355.4 | 0.722 |
| 456862899 | 3,000 | 2.72 | 2.85 | 2.54 | 0.17 | 93.6% | 103.5 | 0.722 |

Every 100 Myr checkpoint is in `gates.json`. By 500 Myr:

| Seed | Craton-in-genuine area kept | Host kept | Net host loss (M km²) | Conversion bound (M km², % of net loss) | Other loss, at least |
|---|---:|---:|---:|---:|---:|
| 7 | 44.2% | 5.3% | 148.10 | 1.21 (0.8%) | 99.2% |
| 11 | 59.8% | 7.0% | 74.10 | 5.77 (7.8%) | 92.2% |
| 2024 | 117.1% | 12.8% | 67.20 | 5.97 (8.9%) | 91.1% |
| 456862899 | 39.2% | 4.5% | 43.47 | 0.03 (0.1%) | 99.9% |

Cratons survive better than the crust around them. They keep 8.3–9.2 times as large a
fraction of their starting area as the host does, and that is what raises the share. At
least 91.1% of each run's net host loss came from something other than formation:
destruction, or thinning below 28 km. On seed 2024, the one run whose craton area grew,
formation can explain at most 6.0 of the 67.2 M km² lost. The 30–65% band therefore fails
here from survival bias, not from formation.

## Which windows can test a source

A source needs host to act on. Three qualifiers for each 500 Myr window, each evaluated on
this baseline:

| Seed | Window (Myr) | Genuine at start | Genuine minimum | Host at start (M km², % of initial) | Original rule at start / throughout | Revised rule | Post-pulse formed (M km³) |
|---|---|---:|---:|---:|---|---|---:|
| 7 | 0–500 | 100.0% | 20.8% | 156.42 (100.0%) | yes / no | yes | 4.778 |
| 7 | 500–1,000 | 20.8% | 9.0% | 8.32 (5.3%) | no / no | no | 0 |
| 11 | 0–500 | 100.0% | 33.8% | 79.69 (100.0%) | yes / yes | yes | 32.548 |
| 11 | 500–1,000 | 33.8% | 16.7% | 5.59 (7.0%) | yes / no | no | 0 |
| 2024 | 0–500 | 100.0% | 51.9% | 77.03 (100.0%) | yes / yes | yes | 42.164 |
| 2024 | 500–1,000 | 51.9% | 19.1% | 9.84 (12.8%) | yes / no | **yes** | 0 |
| 456862899 | 0–500 | 100.0% | 20.0% | 45.51 (100.0%) | yes / no | yes | 0.722 |
| 456862899 | 500–1,000 … 2,500–3,000 | ≤ 20.0% | — | ≤ 2.03 (≤ 4.5%) | no / no | no | 0 |

- The **original rule** (genuine area ≥ 25% of initial) doesn't say whether it is tested
  at the window start or throughout. That choice changes the result for four windows.
  Tested at the start, it admits seed 11's 500–1,000 Myr window, which has only
  5.59 M km² of host left: surviving craton holds up the genuine area.
- The **revised rule** (host at the window start ≥ 10% of its initial value) qualifies the
  four 0–500 Myr windows and seed 2024's 500–1,000 Myr window. Any floor in
  (7.0%, 12.8%] selects the same windows.
- The 0–500 Myr windows pass on the baseline from post-pulse formation alone. The **only**
  qualifying window that starts after the formation shutoff, and so can test a restart, is
  seed 2024's 500–1,000 Myr, and its baseline formation there is zero. The 1–2 Gyr windows
  that the RFC's 2 Gyr runs would add can't qualify on this baseline: seed 456862899 keeps
  at most 3.2% of its host after 1 Gyr.

## Conclusions for RFC 0004

1. On `2562f63` the missing-source premise can be tested on one window of one seed.
   Collapse of the host pool, not the absence of a source, decides whether any formation
   path has material to act on after 500 Myr.
2. Sink repairs (RFC 0001's merge placement, then RFC 0005's #408 sink ranking and
   RFC 0002) should land before RFC 0004 is implemented or tuned, and the gate should be
   re-measured on that baseline.
3. Acceptance must exclude the generation pulse, qualify windows on host from the matched
   baseline arm, and compare ledger flows (formed volume and area) rather than end-state
   stocks or raw shares.
4. The ratchet test must measure conversion, the area formation takes from host, not the
   share. The share's numerator must be restricted to the genuine set.
5. The pre-#358 save's 500 Myr shutoff stays historical evidence only.

The revised gate is written into [RFC 0004](../../rfc/abandoned/0004-proto-cratons.md).

## Limitations

- Checkpoints are 100 Myr apart. Window boundaries fall on checkpoints, but the
  "throughout" qualifier can miss a dip between them.
- The pulse is identified from the 200 Myr sample in #403's largest-jump list. On every
  run that sample is either the largest jump or absent (seed 456862899).
- Formed area isn't recorded, so conversion is bounded from volume. The bound assumes
  every formed column sat at the 28 km host floor, so real conversion is smaller. The
  other-loss figure is a floor on net terms: host gained between checkpoints isn't
  measured, so gross loss is larger still.
