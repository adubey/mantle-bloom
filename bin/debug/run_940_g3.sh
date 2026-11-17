#!/usr/bin/env bash
# #940, RFC 0015 gate G3: run the glacial erosion comparison and knob sensitivity
# replays one at a time, resumably.
#
# Runs gate_g3_glacial_erosion.py on the two non-ice-age #548 diagnosis saves:
#   - 12 main runs: 2 saves x 3 arms (main, cook2020, humphrey1994) x 2 step sizes
#     (1.0 Myr x 10 steps, 0.25 Myr x 40 steps = 10 Myr each)
#   - 8 knob sensitivity runs (cook2020, 1.0 Myr x 10 steps): multiplier 0.25x and 4x,
#     both full (abrasion + flattening) and abrasion-only (flatten_multiplier = 1.0).
#
#   OUT=~/mantle-bloom-captures/issue-940 bin/debug/run_940_g3.sh
#
# Saves are read from $SAVES (default $MANTLE_BLOOM_SAVES_DIR, else ~/Downloads). Set
# NO_WAIT=1 to skip the compute-budget wait (an owner override).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="$ROOT/backend/.venv/bin/python"
OUT="${OUT:-$HOME/mantle-bloom-captures/issue-940}"
SAVES="${SAVES:-${MANTLE_BLOOM_SAVES_DIR:-$HOME/Downloads}}"
mkdir -p "$OUT"

SAVES_LIST="649658380|mantle-bloom-seed649658380-6400000y.bloom|2
997271774|mantle-bloom-seed997271774-87400000y.bloom|2"

run_one() {
  local name="$1"
  local save="$2"
  local gb="$3"
  shift 3
  if [[ -f "$OUT/$name.json" ]]; then
    echo "skip $name"
    return
  fi
  if [[ -z "${NO_WAIT:-}" ]]; then
    until python3 "$ROOT/bin/compute_budget.py" check --run sim --gb "$gb" >/dev/null 2>&1; do sleep 120; done
  fi
  echo "$(date '+%F %T') start $name"
  PYTHONPATH="$ROOT/backend" "$PY" "$ROOT/bin/debug/gate_g3_glacial_erosion.py" "$save" \
    "$@" --out "$OUT/$name.json" >"$OUT/$name.log" 2>&1
  echo "$(date '+%F %T') done $name"
}

while IFS='|' read -r seed file gb; do
  save="$SAVES/$file"
  [[ -f "$save" ]] || { echo "missing save: $save" >&2; exit 1; }
  for arm in main cook2020 humphrey1994; do
    run_one "$seed-$arm-1.0" "$save" "$gb" --arm "$arm" --steps 10 --myr 1.0
    run_one "$seed-$arm-0.25" "$save" "$gb" --arm "$arm" --steps 40 --myr 0.25
  done
  for mult in 0.25 4.0; do
    run_one "$seed-cook2020-1.0-mult$mult" "$save" "$gb" --arm cook2020 --steps 10 --myr 1.0 --multiplier "$mult"
    run_one "$seed-cook2020-1.0-mult$mult-flat1" "$save" "$gb" --arm cook2020 --steps 10 --myr 1.0 --multiplier "$mult" --flatten-multiplier 1.0
  done
done <<<"$SAVES_LIST"

echo "$(date '+%F %T') all runs done"
