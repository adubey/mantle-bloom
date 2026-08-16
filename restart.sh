#!/usr/bin/env bash
# Start (or restart) the backend and frontend dev servers in the background.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

./stop.sh

echo "Starting backend on :8000..."
(
  cd backend
  source .venv/bin/activate
  nohup uvicorn app.main:app --host 127.0.0.1 --port 8000 > /tmp/mantle-bloom-backend.log 2>&1 &
)

echo "Starting frontend dev server on :5173..."
(
  cd frontend
  nohup npm run dev > /tmp/mantle-bloom-frontend.log 2>&1 &
)

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
