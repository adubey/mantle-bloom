#!/usr/bin/env bash
# Dead-code and duplicate-code checks, the same ones CI runs (complexity is part of the
# ordinary lint: ruff's C901 and oxlint's eslint/complexity, both capped at 20).
#   - vulture: unused Python code. Findings that predate the check live in
#     vulture_whitelist.py (see vulture.toml).
#   - knip: unused frontend files, exports and dependencies (frontend/knip.json).
#   - jscpd: copy-pasted blocks across backend/app, frontend/src and bin (.jscpd.json). Only
#     clones missing from .jscpd-baseline.json fail. Editing inside an existing clone changes
#     its fingerprint, so after deduplicating or deliberately keeping one, refresh the
#     baseline with `cd frontend && npm run dupes:baseline` and commit it.
# Run from anywhere; needs ./bin/setup.sh to have run.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

status=0
echo "== vulture (Python dead code)"
backend/.venv/bin/vulture --config vulture.toml || status=1
echo "== knip (frontend dead code)"
(cd frontend && npm run --silent deadcode) || status=1
echo "== jscpd (duplicate code)"
(cd frontend && npm run --silent dupes) || status=1
exit $status
