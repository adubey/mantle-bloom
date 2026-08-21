#!/usr/bin/env bash
# Start (or restart) the backend and frontend servers in the background.
#
# The frontend runs its production build (`vite build` + `vite preview`) by default -- the
# same optimized, minified bundle a real deploy would serve, so this is what to reach for to
# check how the app actually behaves in prod. Pass --dev to run the Vite dev server instead
# (hot module reload, unminified, slower first paint) for active frontend development.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

dev=false
for arg in "$@"; do
  case "$arg" in
    --dev) dev=true ;;
    *)
      echo "Unknown argument: $arg (only --dev is supported)" >&2
      exit 1
      ;;
  esac
done

"$SCRIPT_DIR/stop.sh"

echo "Starting backend on :8000..."
(
  cd backend
  source .venv/bin/activate
  nohup uvicorn app.main:app --host 127.0.0.1 --port 8000 > /tmp/mantle-bloom-backend.log 2>&1 &
)

if [ "$dev" = true ]; then
  echo "Starting frontend dev server on :5173..."
  (
    cd frontend
    nohup npm run dev > /tmp/mantle-bloom-frontend.log 2>&1 &
  )
else
  echo "Building frontend for production..."
  (cd frontend && npm run build) > /tmp/mantle-bloom-frontend.log 2>&1
  echo "Starting frontend (production build) on :5173..."
  (
    cd frontend
    nohup npm run preview -- --port 5173 --strictPort >> /tmp/mantle-bloom-frontend.log 2>&1 &
  )
fi

for i in $(seq 1 30); do
  backend_up=$(curl -s -o /dev/null -w "%{http_code}" http://127.0.0.1:8000/docs || true)
  frontend_up=$(curl -s -o /dev/null -w "%{http_code}" http://localhost:5173/ || true)
  if [ "$backend_up" = "200" ] && [ "$frontend_up" = "200" ]; then
    echo "Backend:  http://127.0.0.1:8000 (logs: /tmp/mantle-bloom-backend.log)"
    echo "Frontend: http://localhost:5173 (logs: /tmp/mantle-bloom-frontend.log)"
    exit 0
  fi
  sleep 1
done

echo "Timed out waiting for servers to come up -- check the logs above." >&2
exit 1
