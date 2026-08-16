#!/usr/bin/env bash
# Stop the backend and frontend dev servers started by restart.sh.
set -uo pipefail

for port in 8000 5173 5174; do
  pids=$(lsof -ti:"$port" -sTCP:LISTEN 2>/dev/null || true)
  if [ -n "$pids" ]; then
    echo "Freeing port $port (pid(s): $pids)"
    echo "$pids" | xargs kill 2>/dev/null || true
  fi
done
