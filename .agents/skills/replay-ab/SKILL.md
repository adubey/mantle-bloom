---
name: replay-ab
description: Measure a mantle-bloom change by replaying a saved world (.mbworld) on origin/main and on the branch, then writing it up under analysis/issue<N>/. Use when an issue or review asks for replay numbers, a before/after comparison, or whether a change moves land area, crust/mantle thickness, ledgers or runtime.
---

# Replay A/B against origin/main

1. **Pick the save and horizon.** Use the ones the issue names. Saves live in `~/Downloads`
   as `mantle-bloom-seed<seed>-<years>y.mbworld`. Record each save's SHA-256
   (`shasum -a 256`), since earlier analyses cite it. If the save isn't there, ask the user;
   don't swap in another one.

2. **Pick the replay script.** Reuse one from `bin/debug/` and read its docstring first.
   `attribute_hm_ratchet.py` is the general one: it reports checkpoints, Hc/Hm and crust
   ledgers, land and hypsometry, and runtime and peak RSS. Extend an existing script, or add
   a new one under `bin/debug/` with a usage docstring, only when nothing measures what you
   need.

3. **Baseline = a detached `origin/main` worktree** (see the `start-issue` skill). Record
   its commit. Never `git stash` to get a "before" run: the stash is shared by every
   worktree, including the user's own. Never use the primary checkout either, because its
   `main` lags `origin`.

4. **Run replays in series.** The machine is short on RAM: each replay takes 0.6–2.7 GB.
   Run at most two at once (the baseline and branch pair on a small save), and never
   alongside a unit or stress test run. A 43 Myr replay of the seed 997271774 save takes about
   40 minutes. Run it in the background, and write outputs to a scratch directory, not the
   repo:

   ```bash
   PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
     ~/Downloads/<save>.mbworld --myr 43.1 --checkpoint-myr 2.5 --out <scratch>/branch.json
   ```

5. **Allow for noise.** Run-to-run noise on continental area is about ±2% at 43 Myr (#315).
   A gap inside that range is not a result. When a gap matters, compare the trajectory
   across checkpoints, not just the end point, or replay a second save.

6. **Write it up** in `analysis/issue<N>/README.md`, using
   [`analysis/issue320/README.md`](../../../analysis/issue320/README.md) as the model:
   - what was replayed: the save, its SHA, the steps, and the commit of each run;
   - the exact commands;
   - a table of measures per run at the end point;
   - what the numbers mean, including the noise caveat;
   - runtime and RSS, flagged when a run didn't run alone.

   Check in a compact comparison JSON. Don't check in full replay outputs or saves. Put the
   headline table in the PR description too.
