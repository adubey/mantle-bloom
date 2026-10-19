---
name: replay-ab
description: Measure a Mantle Bloom change by replaying a saved world (.mbworld) on origin/main and on the branch, then writing results in the relevant analysis issue directory. Use when an issue or review asks for replay numbers, a before/after comparison, or whether a change moves land area, crust/mantle thickness, ledgers or runtime.
---

# Replay A/B against origin/main

1. **Pick the save and horizon.** Use the ones the issue names. Resolve local saves only
   from a provided path or `MANTLE_BLOOM_SAVES_DIR`; do not broadly scan personal
   directories. If the variable is unset or the named save is missing, ask the user; don't
   substitute another save. Record each save's SHA-256 (`shasum -a 256`), since earlier
   analyses cite it.

2. **Pick the replay script.** Reuse one from `bin/debug/` and read its docstring first.
   `attribute_hm_ratchet.py` is the general one: it reports checkpoints, Hc/Hm and crust
   ledgers, land and hypsometry, and runtime and peak RSS. Extend an existing script, or add
   a new one under `bin/debug/` with a usage docstring, only when nothing measures what you
   need.

3. **Baseline = a detached `origin/main` worktree** (see the `start-issue` skill). Record
   its commit. Never `git stash` to get a "before" run: the stash is shared by every
   worktree, including the user's own. Never use the primary checkout either, because its
   `main` lags `origin`.

4. **Run the replay.** Follow the concurrency limits in `AGENTS.md`; never run replay work
   alongside a unit or stress test. A 43 Myr replay of the seed 997271774 save takes about
   40 minutes. Run it in the background, and write outputs to a scratch directory, not the
   repo:

   ```bash
   PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
     <local-save> --myr 43.1 --checkpoint-myr 2.5 --out <scratch-output>
   ```

5. **Pre-register and interpret the comparison.** Before running, record the expected
   direction and rough size, primary metric, decision threshold, controls, and seed(s) in
   the analysis README. Estimate variability for the metric and run context with baseline
   replicates when needed. The approximately ±2% continental-area variation at 43 Myr in
   #315 applies only to a comparable measurement; don't generalize it to other metrics or
   runs. A gap within measured variability is inconclusive. When a gap matters, compare the
   trajectory across checkpoints or use a held-out save. Record every run, including
   discarded runs and why, to avoid selective reporting.

6. **Write it up** in `analysis/issue<N>/README.md`, using
   [`analysis/issue320/README.md`](../../../analysis/issue320/README.md) as the model:
   - what was replayed: the save, its SHA, the steps, and the commit of each run;
   - the command with local input and scratch paths replaced by neutral placeholders such
     as `<local-save>` and `<scratch-output>`;
   - a table of measures per run at the end point;
   - what the numbers mean, including the noise caveat;
   - runtime and RSS, flagged when a run didn't run alone.

   Check in a compact comparison JSON. Don't check in full replay outputs or saves. For
   private machine paths, follow the local-only evidence handoff procedure in
   [`AGENTS.md`](../../../AGENTS.md#workflow). Put the headline table in the PR description.
