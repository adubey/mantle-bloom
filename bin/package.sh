#!/usr/bin/env bash
# Build the single-process desktop app into dist/mantle-bloom/ (and, on macOS,
# dist/mantle-bloom.app). Run on the OS you want to ship for -- PyInstaller does not
# cross-compile. See docs/packaging.md.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

echo "==> Building frontend (production bundle)"
(cd frontend && npm ci && VITE_API_BASE="" npm run build)

echo "==> Preparing backend build venv"
cd backend
if [ ! -d .venv ]; then
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
pip install -q -r requirements-build.txt
cd "$REPO_ROOT"

echo "==> Running PyInstaller"
rm -rf build dist
pyinstaller bin/mantle-bloom.spec --noconfirm --clean

echo
echo "==> Done. Artifacts in dist/:"
ls -1 dist/
echo
echo "Smoke-test it:  ./dist/mantle-bloom/mantle-bloom"
