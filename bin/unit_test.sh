#!/usr/bin/env bash
# Run the backend's fast unit tests (backend/unit_tests/) -- every test under ~1s, meant to
# be run often during development. The suite's slow, full-simulation tests (many-step
# integration/determinism checks) live in backend/stress_tests/ instead -- see
# bin/stress_test.sh -- so this stays quick. Extra arguments are forwarded to pytest, e.g.
# `./bin/unit_test.sh -k biome` or `./bin/unit_test.sh unit_tests/test_hydrology.py`.
#
# Runs across all available cores via pytest-xdist (-n auto), with xdist's default `load`
# distribution: tests go to whichever worker is free, so one file's tests can be split across
# workers and run in any order. That relies on every test setting up its own state --
# test_main.py's `client` fixture resets app/main.py's module-level world/stop-signal globals
# around each test for exactly this (#334); it used to run under --dist loadscope, which kept
# each file on one worker in source order, and its largest files finished last (full suite
# 7:02 under loadscope vs 6:27 under load on the same machine). A passed -n/--dist further
# down this command line overrides these, since pytest takes the last value given.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$REPO_ROOT/backend"
source .venv/bin/activate
python -m pytest unit_tests/ -q -n auto "$@"
