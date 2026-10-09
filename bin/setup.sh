#!/usr/bin/env bash
# One-step setup for a fresh clone or a new worktree: creates backend/.venv (if missing) and
# installs the backend's runtime, test and lint deps at the versions pinned in
# backend/constraints.txt, then installs the frontend's deps from frontend/package-lock.json.
# Safe to re-run any time -- pip skips what's already satisfied, and npm ci only runs when
# package-lock.json is newer than the last install.
#
# Each worktree gets its own venv: don't symlink another tree's, since a branch that changes
# backend/constraints.txt would then upgrade the venv every other tree runs from.
#
# PYTHON overrides the interpreter used to create the venv (default python3, needs 3.10+).
# --backend-only skips the frontend.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

frontend=true
while [ $# -gt 0 ]; do
  case "$1" in
    --backend-only) frontend=false; shift ;;
    *)
      echo "Unknown argument: $1 (supported: --backend-only)" >&2
      exit 1
      ;;
  esac
done

cd "$REPO_ROOT/backend"
if [ ! -x .venv/bin/python ]; then
  "${PYTHON:-python3}" -c 'import sys; sys.exit(sys.version_info < (3, 10))' \
    || { echo "Need Python 3.10+ (set PYTHON=/path/to/python3.x)" >&2; exit 1; }
  echo "Creating backend/.venv with $("${PYTHON:-python3}" --version)..."
  "${PYTHON:-python3}" -m venv .venv
  .venv/bin/python -m pip install -q --upgrade pip
fi
echo "Installing backend deps..."
.venv/bin/python -m pip install -q -r requirements.txt -c constraints.txt

if [ "$frontend" = true ]; then
  cd "$REPO_ROOT/frontend"
  if [ ! -f node_modules/.package-lock.json ] || [ package-lock.json -nt node_modules/.package-lock.json ]; then
    echo "Installing frontend deps..."
    npm ci --no-audit --no-fund
  else
    echo "Frontend deps up to date."
  fi
fi

echo "Done. Tests: ./bin/affected_test.sh (or ./bin/unit_test.sh); app: ./bin/restart.sh"
