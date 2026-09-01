#!/usr/bin/env bash
# Start (or restart) the backend and frontend servers in the background.
#
# The frontend runs its production build (`vite build` + `vite preview`) by default -- the
# same optimized, minified bundle a real deploy would serve, so this is what to reach for to
# check how the app actually behaves in prod. Pass --dev to run the Vite dev server instead
# (hot module reload, unminified, slower first paint) for active frontend development.
#
# --backend-port/--frontend-port override the default 8000/5173 ports (e.g. to run a second,
# independent instance alongside one already up on the defaults). The frontend is always
# wired to reach the backend at its actual port (VITE_API_BASE, baked in at build/dev-server-
# start time -- see frontend/src/api.ts), and the backend's CORS allowlist is pointed at the
# actual frontend port (FRONTEND_PORT -- see backend/app/main.py), so a non-default pair
# still talks to itself correctly.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$REPO_ROOT"

dev=false
backend_port=8000
frontend_port=5173
stay_awake=false
while [ $# -gt 0 ]; do
  case "$1" in
    --dev) dev=true; shift ;;
    --backend-port)
      backend_port="$2"; shift 2 ;;
    --frontend-port)
      frontend_port="$2"; shift 2 ;;
    --stay-awake)
      stay_awake=true; shift 2 ;;
    *)
      echo "Unknown argument: $1 (supported: --dev, --backend-port PORT, --frontend-port PORT)" >&2
      exit 1
      ;;
  esac
done

"$SCRIPT_DIR/stop.sh" --backend-port "$backend_port" --frontend-port "$frontend_port"
if [ "$stay_awake" = true ]; then
  CAFFEINATE=caffeinate 
else
  CAFFEINATE=
fi

echo "Starting backend on :$backend_port..."
(
  cd backend
  source .venv/bin/activate
  FRONTEND_PORT="$frontend_port" nohup "$CAFFEINATE" uvicorn app.main:app --host 127.0.0.1 --port "$backend_port" > /tmp/mantle-bloom-backend.log 2>&1 &
)

api_base="http://localhost:$backend_port"

if [ "$dev" = true ]; then
  echo "Starting frontend dev server on :$frontend_port..."
  (
    cd frontend
    VITE_API_BASE="$api_base" nohup "$CAFFEINATE" npm run dev -- --port "$frontend_port" --strictPort > /tmp/mantle-bloom-frontend.log 2>&1 &
  )
else
  echo "Building frontend for production..."
  (cd frontend && VITE_API_BASE="$api_base" npm run build) > /tmp/mantle-bloom-frontend.log 2>&1
  echo "Starting frontend (production build) on :$frontend_port..."
  (
    cd frontend
    nohup "$CAFFEINATE" npm run preview -- --port "$frontend_port" --strictPort >> /tmp/mantle-bloom-frontend.log 2>&1 &
  )
fi

for i in $(seq 1 30); do
  backend_up=$(curl -s -o /dev/null -w "%{http_code}" "http://127.0.0.1:$backend_port/docs" || true)
  frontend_up=$(curl -s -o /dev/null -w "%{http_code}" "http://localhost:$frontend_port/" || true)
  if [ "$backend_up" = "200" ] && [ "$frontend_up" = "200" ]; then
    echo "Backend:  http://127.0.0.1:$backend_port (logs: /tmp/mantle-bloom-backend.log)"
    echo "Frontend: http://localhost:$frontend_port (logs: /tmp/mantle-bloom-frontend.log)"
    exit 0
  fi
  sleep 1
done

echo "Timed out waiting for servers to come up -- check the logs above." >&2
exit 1
