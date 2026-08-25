#!/usr/bin/env bash
# Stop the backend and frontend dev servers started by restart.sh.
#
# --backend-port/--frontend-port target a non-default pair (matching the ports restart.sh was
# given) instead of the default 8000/5173, so restarting a custom-port instance doesn't kill a
# default-port instance running alongside it, and vice versa.
set -uo pipefail

backend_port=8000
frontend_port=5173
while [ $# -gt 0 ]; do
  case "$1" in
    --backend-port)
      backend_port="$2"; shift 2 ;;
    --frontend-port)
      frontend_port="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1 (supported: --backend-port PORT, --frontend-port PORT)" >&2
      exit 1
      ;;
  esac
done

# 5174 is Vite dev's own fallback port when its primary port is taken, so it's always worth
# checking regardless of which frontend port was requested.
for port in "$backend_port" "$frontend_port" 5174; do
  pids=$(lsof -ti:"$port" -sTCP:LISTEN 2>/dev/null || true)
  if [ -n "$pids" ]; then
    echo "Freeing port $port (pid(s): $pids)"
    echo "$pids" | xargs kill 2>/dev/null || true
  fi
done
