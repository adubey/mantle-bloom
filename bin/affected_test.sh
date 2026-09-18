#!/usr/bin/env bash
# Run only the backend unit tests affected by this branch's changes -- a faster loop than
# bin/unit_test.sh's full suite when you're iterating on a handful of files. "Affected" is
# computed by bin/list_affected_tests.py: the change set is the diff against main plus
# anything currently uncommitted, mapped to test files through app/'s own import graph (a
# module's tests, plus every module that transitively imports it). Pass --base <ref> to
# compare against something other than main; anything else is forwarded to pytest, e.g.
# `./bin/affected_test.sh -v`.
#
# Falls back to the full suite whenever the affected set can't be trusted (a file the import
# graph can't map, e.g. a new __init__.py) -- see list_affected_tests.py's docstring. Always
# review with bin/unit_test.sh before pushing; this is a dev-loop shortcut, not a substitute
# for the full suite.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# --base can appear anywhere, not just first -- everything else is forwarded to pytest
# untouched, in order (e.g. `./bin/affected_test.sh -v --base main` must not hand pytest a
# literal "--base main" it doesn't understand).
BASE_ARGS=()
PYTEST_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base)
      BASE_ARGS=(--base "$2")
      shift 2
      ;;
    *)
      PYTEST_ARGS+=("$1")
      shift
      ;;
  esac
done

# Word-split intentionally (no mapfile -- macOS ships bash 3.2, which lacks it): every path
# list_affected_tests.py prints is a plain unit_tests/test_*.py filename, never containing a
# space.
AFFECTED="$(python3 "$SCRIPT_DIR/list_affected_tests.py" "${BASE_ARGS[@]+"${BASE_ARGS[@]}"}")"
if [[ -z "$AFFECTED" ]]; then
  echo "No affected tests." >&2
  exit 0
fi

cd "$REPO_ROOT/backend"
source .venv/bin/activate
python -m pytest $AFFECTED -q -n auto --dist loadscope "${PYTEST_ARGS[@]+"${PYTEST_ARGS[@]}"}"
