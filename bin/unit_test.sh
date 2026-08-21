#!/usr/bin/env bash
# Run the backend's pytest suite -- the only unit tests in this repo (the frontend has no
# test framework set up yet, see frontend/package.json's scripts). Extra arguments are
# forwarded to pytest, e.g. `./bin/unit_test.sh -k biome` or `./bin/unit_test.sh tests/test_hydrology.py`.
#
# Runs across all available cores via pytest-xdist (-n auto). --dist loadscope keeps every
# file's own tests on a single worker rather than freely interleaving them: test_main.py's
# tests share app/main.py's module-level `_state` global through a `TestClient(app)` fixture,
# and several of them (e.g. test_render_before_generate_returns_404) depend on running before
# any test in that same file has generated a world -- true today only because plain pytest
# runs a file's tests in source order. loadscope preserves that per-file order while still
# parallelizing freely *across* files (no other test file touches `_state`). A passed -n/
# --dist further down this command line overrides these, since pytest takes the last value
# given for either.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$REPO_ROOT/backend"
source .venv/bin/activate
python -m pytest tests/ -q -n auto --dist loadscope "$@"
