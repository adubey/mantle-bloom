#!/usr/bin/env bash
# Run the backend's slow, full-simulation tests (backend/stress_tests/) -- many-step
# integration and determinism checks that individually take anywhere from a few seconds to
# several minutes (see backend/unit_tests/ via bin/unit_test.sh for the fast suite meant for
# everyday use). Extra arguments are forwarded to pytest, e.g. `./bin/stress_test.sh -k gaps`.
#
# Runs across all available cores via pytest-xdist (-n auto), with xdist's default `load`
# distribution -- see bin/unit_test.sh's own comment: tests don't depend on file order, since
# test_main.py's `client` fixture resets app/main.py's globals around each test (#334).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$REPO_ROOT/backend"
source .venv/bin/activate
python -m pytest stress_tests/ -q -n auto "$@"
