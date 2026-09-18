# Working in this repo

## Running backend tests

Prefer `./bin/affected_test.sh` over `./bin/unit_test.sh` while iterating on
`backend/app/` or `backend/unit_tests/` changes -- it runs only the unit tests
affected by what's changed on the branch (committed since `main` plus anything
uncommitted), via `bin/list_affected_tests.py`'s static import-graph walk. The full
suite is ~900 tests and takes several minutes; the affected subset is usually a small
fraction of that and finishes in seconds.

It falls back to the full suite on its own whenever the affected set can't be trusted
(e.g. a change to `app/__init__.py`, or something under an `app/` subdirectory like
`app/data/`), so it's safe to reach for by default.

Still run the full `./bin/unit_test.sh` before finishing a task or pushing --
`affected_test.sh` is a dev-loop shortcut, not a substitute for the full suite as a
final check.
