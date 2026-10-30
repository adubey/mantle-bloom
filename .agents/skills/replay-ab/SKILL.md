---
name: replay-ab
description: Measure a Mantle Bloom change by replaying a saved world (.mbworld) on origin/main and on the branch, then writing results in the relevant analysis issue directory. Use when an issue or review asks for replay numbers, a before/after comparison, or whether a change moves land area, crust/mantle thickness, ledgers or runtime.
---

# Replay A/B against origin/main

1. **Pick the save and horizon.** Prefer an explicit local save path named by the issue or
   user. Otherwise use `MANTLE_BLOOM_SAVES_DIR`; when unset, source the private
   `backend/.venv/mantle-bloom-env.sh` created by `./bin/setup.sh`, which defaults to the
   frontend Save location `/Users/$USER/Downloads`. Precedence is explicit path, existing
   environment variable, then the setup default. Keep `$USER` literal in shared instructions
   and never write a resolved personal path into a report, issue, PR or commit. Do not scan
   other personal directories. If the named save is missing, ask the user; don't substitute
   another save. Record each save's SHA-256 (`shasum -a 256`), since earlier analyses cite it.

2. **Pick the replay script.** Reuse one from `bin/debug/` and read its docstring first.
   `attribute_hm_ratchet.py` is the general one: it reports checkpoints, Hc/Hm and crust
   ledgers, land and hypsometry, and runtime and peak RSS. Extend an existing script, or add
   a new one under `bin/debug/` with a usage docstring, only when nothing measures what you
   need.

3. **Baseline = a task-specific detached `origin/main` worktree** (see `start-issue`).
   Never share or reuse another task's baseline worktree. Record its commit. Never
   `git stash` to get a "before" run: the stash is shared by every worktree, including the
   user's own. Never use the primary checkout either, because its `main` lags `origin`.

4. **Pre-register the comparison.** Before running any replay, record the expected direction
   and rough size, primary metric, decision threshold, controls, and seed(s) in the analysis
   README, and commit the preregistration before the run. Estimate variability for each
   primary metric; the approximately ±2% variation in
   continental area at 43 Myr in #315 applies only to a comparable measurement. If the result
   is near measured variability, plan to replicate both baseline and branch runs and report
   each arm's spread.

5. **Run the replay.** Run the arms one after the other, and follow
   [`compute-budget`](../compute-budget/SKILL.md) before starting each. A 43 Myr replay of the seed 997271774 save takes about 40 minutes. Run it
   in the background, and write outputs to a scratch directory, not the repo:

   ```bash
   (
     source backend/.venv/mantle-bloom-env.sh
     PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
       "$MANTLE_BLOOM_SAVES_DIR/<save-name>" --myr 43.1 --checkpoint-myr 2.5 \
       --out "<scratch-output>"
   ) > "<scratch-output>.log" 2>&1 &
   ```

   If the issue or user supplied a specific save path, substitute it for the
   `$MANTLE_BLOOM_SAVES_DIR/<save-name>` argument. Keep environment setup and replay in
   the same shell invocation; shell variables do not carry between separate agent tool
   calls.

6. **Interpret the comparison.** If the gap is near measured variability, run the planned
   replicates on both arms and report the spread of each. A gap within measured variability
   is inconclusive. When a gap matters, compare the trajectory across checkpoints or use a
   held-out save. Record every run, including discarded runs and why, to avoid selective
   reporting.

7. **Write it up** in `analysis/issue<N>/README.md`, using
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
