## Summary

<!-- What changed and why, in model terms. When the model moved, give before/after replay
numbers in a table (see .agents/skills/replay-ab). -->

Closes #

## Verification

<!-- What you ran and what it showed. Failing output verbatim; call a failure pre-existing
only after reproducing it on origin/main. -->

- [ ] `backend/.venv/bin/ruff check`
- [ ] `./bin/quality_check.sh` (dead and duplicate code)
- [ ] `./bin/affected_test.sh --base origin/main`: N passed
- [ ] Frontend touched: `cd frontend && npm run lint && npm run build`
- [ ] Stress tests (`./bin/stress_test.sh`), if the step pipeline, world generation or persistence changed
- [ ] Old `.mbworld` saves still load, if pickled state changed (docs/save-compatibility.md)
- [ ] Docs updated, if a module, route or model behaviour was added (docs/architecture.md, api-reference.md, simulation-model.md)
