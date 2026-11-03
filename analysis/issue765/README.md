# Issue 765: D2 revision evidence — invariant under the 250 yr march cap

Evidence for [#765](https://github.com/adubey/mantle-bloom/issues/765), the RFC 0015 D2
revision. Tracking: [#544](https://github.com/adubey/mantle-bloom/issues/544).
Builds on [#741](https://github.com/adubey/mantle-bloom/issues/741) (drops are long-substep
artefacts) and [#766](https://github.com/adubey/mantle-bloom/issues/766) (G1 on the march:
τ/4 under 250 yr on every save; 250 yr cost runs).

**Principle:** #3, fast processes in aggregate over long steps. A fast body's state should be
set by its forcing, not by the integration step.

**Invariant (from #765):** on #741's two drop steps, an ensemble of 12 opening-ice
perturbations (±0.01–0.5%) drops in at most 2 of 12 members, with no substep moving the
volume by more than 20%. The march's 0.25 vs 1 Myr step gap is re-measured under the limit.

## Preregistration

Written before any run below. The only data seen are #741's (9 of 12 drop at a 10 kyr
longest substep, 0 of 12 at a 1 kyr cap, on both steps) and #766's (τ/4 under 250 yr on
every save; the 250 yr march settles in ~13–16 substeps and costs 0.32–0.51× the 10 kyr
march per step).

Tree: gate branch `claude/issue-766-g1-march` at `884fe85` (G0 prototype `ea9a187`,
drifting noise `a94047d`, #741's ensemble probes, `g1_march.py`). March arm (m3):
damped, adaptive, flip fraction 1%, hold window 6, settling tolerance 0.1% over
3 pair-means in a row. Save: seed 997271774 at 87.4 Myr, SHA-256
`08ad9461…ed762b` (matches #741/#766). Draws are #741's `--noise-salt-offset` 1
and 2. Each run alone; peak RSS 1.1–1.2 GB.

| Run | What | Reads |
|---|---|---|
| E1 | Draw 1, 97.4 Myr: 12-member opening-ice ensemble, this step's climate, longest substep 250 yr | Does the 1 kyr result (0 of 12) hold at the revised 250 yr cap? |
| E2 | Draw 2, 96.65 Myr: same 12-member ensemble at 250 yr | Same, on the gradual drop |
| G1 | Draw 1, full 10 Myr replay at 0.25 and 1 Myr steps, longest substep 250 yr | Step gap under the limit, against G0's as-run 17.0% on this draw |

**Method (E1, E2):** replay to the step as run (10 kyr longest substep, the opening #741's
drops start from), then re-march the step in-process from its real opening ice with the
longest substep capped at 250 yr:

```bash
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 \
  --checkpoint-myr 0.25 --noise-salt-offset <1|2> --myr <10|9.25> \
  --factorial-myr <97.4|96.65> --ensemble-only --ensemble-max-substep-years 250 \
  --capture-dir <scratch-output> --out <scratch-output>.json
```

**Reading:** a drop is a member closing below 11.5 M km³ (#741's threshold). The cap
removes the drop if at most 2 of 12 members drop (against E1's 9 of 12 at 10 kyr) and no
member has a substep that moves the volume over 20%. If 8 or more still drop with no such
substep, the drop is a real transition the long substeps only hasten. Anything between is
reported as it stands.

**Expected:** 0 of 12 on both steps (the 1 kyr cap already gives 0 of 12; 250 yr only
shortens substeps further), with the E8/E9 split states (about 12.3 vs 13.2 M km³ on draw
2) still visible. The G1 run is report-only: no threshold is set on the step gap, since no
under-limit gap has ever been measured.

**Not tested:** other draws' drops; the ice-age save (its runaway is G4's problem, and G0
isn't closed by this issue).

## Results

### E1: draw 1, 97.4 Myr, longest substep 250 yr

Replay to the step as run (10 kyr longest substep), then the 13-cell opening-ice
ensemble re-marched at a 250 yr cap. Command as preregistered; tree `884fe85`.

| Cell | Closing (M km³) | Worst substep move |
|---|---|---|
| control | 12.96 | 0.68% |
| ±0.01–0.5% (12) | 12.96–13.18 | 0.57–0.74% |

**0 of 12 drop** (drop = closing below 11.5 M km³), and no substep moves the volume by
more than 20% — the largest is 0.74%. The invariant holds on draw 1 at 250 yr, as at
1 kyr (E8: 0 of 12). Raw log and capture: `e1-250.log`, `e1-250/` in
`~/mantle-bloom-captures/issue-765/` (private path; run label `e1-250`).

### E2: draw 2, 96.65 Myr, longest substep 250 yr

Same method as E1 with `--noise-salt-offset 2 --myr 9.25 --factorial-myr 96.65`.
Tree `884fe85`.

| Cell | Closing (M km³) | Worst substep move |
|---|---|---|
| control | 12.92 | 0.45% |
| ±0.01–0.5% (12) | 12.21–13.13 | 0.33–0.60% |

**0 of 12 drop**, and no substep moves the volume by more than 20% — the largest is
0.60%. The invariant holds on draw 2 at 250 yr, as at 1 kyr (E9: 0 of 12). The members
still split between about 12.2 and 13.1 M km³, the same two states E9 saw (12.3 vs
13.2): the cap removes the drop, not the bistability. Raw log and capture:
`e2-250.log`, `e2-250/`, run label `e2-250`.

### G1 step gap: draw 1, 0.25 vs 1 Myr steps, longest substep 250 yr

Full 10 Myr replays from 87.4 Myr at both step sizes under the revised march
(`--set march_max_substep_years=250`), draw 1. Tree `884fe85`.

| Arm | Steps | End volume (M km³) | End area (M km²) | Max swing after 2 Myr | Mean step time |
|---|---|---|---|---|---|
| 0.25 Myr | 40 | 13.57 | 14.37 | 2.08% | 29.7 s |
| 1 Myr | 10 | 12.26 | 14.64 | 0.39% | 57.7 s |

**Step gap: 9.71%**, against G0's as-run 17.0% on this draw. The cap more than halves
the gap and brings this draw under the 10% "step robustness" bar (acceptance is on
held-out saves; as-run, the march exceeds 10% on every draw). Both arms stay far under
the 25% swing limit. Raw logs: `gap-025.log`, `gap-1.log`, run labels `gap-025`,
`gap-1`.

### Cost note

The 1 Myr-step march costs more per step (57.7 s) than the 0.25 Myr one (29.7 s):
settling takes the same ~4 kyr either way, and each step re-marches it. Step time is
set by substeps-to-settle, not by step length — the hold makes long steps cheap per
simulated year, not per step.
