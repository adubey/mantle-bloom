#!/usr/bin/env bash
# Start (or restart) mantle-bloom in the background.
#
# Default: the single-process app (backend + production frontend bundle served from one
# origin by FastAPI, see backend/app/desktop.py). Builds the frontend, then runs it on
# :8000 -- one process, one port, the same thing bin/package.sh freezes into a binary.
#
# --dev: the two-process development setup instead -- uvicorn + the Vite dev server (hot
# module reload, unminified) on a separate port, for active frontend work. Only --dev needs
# a second port.
#
# --port overrides the default 8000 (e.g. to run a second instance alongside one already
# up); --frontend-port overrides the Vite dev port (5173, --dev only). Logs go to
# /tmp/mantle-bloom-<port>.log (and /tmp/mantle-bloom-frontend-<frontend-port>.log in --dev)
# so parallel instances don't clobber each other.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

dev=false
port=8000
frontend_port=5173
stay_awake=false
while [ $# -gt 0 ]; do
  case "$1" in
    --dev) dev=true; shift ;;
    --stay-awake) stay_awake=true; shift ;;
    --port|--backend-port)
      port="$2"; shift 2 ;;
    --frontend-port)
      frontend_port="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1 (supported: --dev, --stay-awake, --port PORT, --frontend-port PORT)" >&2
      exit 1
      ;;
  esac
done

"$SCRIPT_DIR/stop.sh" --port "$port" --frontend-port "$frontend_port"
if [ "$stay_awake" = true ]; then
  CAFFEINATE=caffeinate
else
  CAFFEINATE=
fi

if [ "$dev" = true ]; then
  BACKEND_LOGS="/tmp/mantle-bloom-${port}.log"
  FRONTEND_LOGS="/tmp/mantle-bloom-frontend-${frontend_port}.log"
  api_base="http://localhost:$port"

  echo "Starting backend (uvicorn) on :$port with logs ${BACKEND_LOGS}..."
  (
    cd backend
    source .venv/bin/activate
    FRONTEND_PORT="$frontend_port" nohup ${CAFFEINATE:-} uvicorn app.main:app --host 127.0.0.1 --port "$port" > "$BACKEND_LOGS" 2>&1 &
  )

  echo "Starting frontend dev server on :$frontend_port with logs ${FRONTEND_LOGS}..."
  (
    cd frontend
    VITE_API_BASE="$api_base" nohup ${CAFFEINATE:-} npm run dev -- --port "$frontend_port" --strictPort > "$FRONTEND_LOGS" 2>&1 &
  )

  for i in $(seq 1 30); do
    backend_up=$(curl -s -o /dev/null -w "%{http_code}" "http://127.0.0.1:$port/docs" || true)
    frontend_up=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:$frontend_port/" || true)
    if [ "$backend_up" = "200" ] && [ "$frontend_up" = "200" ]; then
      echo "Backend:  http://127.0.0.1:$port"
      echo "Frontend: http://localhost:$frontend_port  <-- open this"
      exit 0
    fi
    sleep 1
  done
  echo "Timed out waiting for the dev servers to come up -- check the logs above." >&2
  exit 1
fi

# Default: single-process app.
LOGS="/tmp/mantle-bloom-${port}.log"
echo "Building frontend (production bundle)..."
(cd frontend && VITE_API_BASE="" npm run build) > "$LOGS" 2>&1

echo "Starting mantle-bloom on :$port with logs ${LOGS}..."
(
  cd backend
  source .venv/bin/activate
  MANTLE_BLOOM_PORT="$port" MANTLE_BLOOM_NO_BROWSER=1 \
    nohup ${CAFFEINATE:-} python -m app.desktop >> "$LOGS" 2>&1 &
)

for i in $(seq 1 30); do
  if [ "$(curl -s -o /dev/null -w "%{http_code}" "http://127.0.0.1:$port/" || true)" = "200" ]; then
    echo "mantle-bloom: http://127.0.0.1:$port  <-- open this"
    exit 0
  fi
  sleep 1
done

echo "Timed out waiting for mantle-bloom to come up -- check ${LOGS}." >&2
exit 1
