---
name: pre-push-check
description: Run mantle-bloom's full verification (ruff, dead/duplicate-code checks, full backend unit suite, frontend lint and build when touched) and report the results. Use before pushing a branch, opening or updating a PR, or saying a task is done.
---

# Pre-push check

Run these from the worktree root, in order, and stop to fix anything that fails.

1. **Is the branch current?** Run `git fetch origin` and
   `git merge-base --is-ancestor origin/main HEAD`. If it's behind and touches the same
   files as recent merges, merge `origin/main` first. Ask before rewriting history that's
   already pushed.

2. **Lint:** `backend/.venv/bin/ruff check`. `--fix` is safe for unused imports. Delete a
   dead local only after checking that it isn't a sign of a missing step. A C901 finding
   means split the function, not add a `noqa`.

3. **Dead and duplicate code:** `./bin/quality_check.sh` (vulture, knip, jscpd; seconds).
   Fix findings rather than adding them to `vulture_whitelist.py` or the jscpd baseline;
   see the ratchet note in AGENTS.md.

4. **Full unit suite:** `./bin/unit_test.sh`, about 12 minutes on 10 cores.
   `affected_test.sh` doesn't count as this check. Don't run the suite while replays are
   going (RAM). You can skip it if this exact HEAD already passed in this session.

5. **Frontend**, only if `git diff --name-only origin/main...HEAD -- frontend` is non-empty:
   `cd frontend && npm run lint && npm run build`. The build runs `tsc`, which is the
   typecheck.

6. **Stress tests:** `./bin/stress_test.sh` is slow. Run it when the change touches the step
   pipeline, world generation or persistence (the areas `backend/stress_tests/` covers).
   Otherwise just say it wasn't run.

7. **Report** each check: pass or fail, the test counts, and failing output verbatim. Only
   call a failure pre-existing after reproducing it on the `origin/main` baseline worktree.
