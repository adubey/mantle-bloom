#!/usr/bin/env bash
# Stop whatever bin/restart.sh started -- the single-process app, or the --dev two-process
# setup.
#
# --port/--frontend-port target a non-default pair (matching what restart.sh was given) so
# stopping a custom-port instance doesn't kill a default one running alongside it.
set -uo pipefail

port=8000
frontend_port=5173
while [ $# -gt 0 ]; do
  case "$1" in
    --port|--backend-port)
      port="$2"; shift 2 ;;
    --frontend-port)
      frontend_port="$2"; shift 2 ;;
    *)
      echo "Unknown argument: $1 (supported: --port PORT, --frontend-port PORT)" >&2
      exit 1
      ;;
  esac
done

# 5174 is Vite dev's own fallback port when its primary is taken, so it's always worth
# sweeping regardless of which frontend port was requested.
for p in "$port" "$frontend_port" 5174; do
  pids=$(lsof -ti:"$p" -sTCP:LISTEN 2>/dev/null || true)
  if [ -n "$pids" ]; then
    echo "Freeing port $p (pid(s): $pids)"
    echo "$pids" | xargs kill 2>/dev/null || true
  fi
done
