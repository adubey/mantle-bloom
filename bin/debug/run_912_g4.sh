#!/usr/bin/env bash
# #912, RFC 0015 gate G4: run the six replays, one at a time, resumably.
#
# Runs measure_flotation_debit.py on the three #548 diagnosis saves at 0.25 and 1 Myr steps
# (10 steps each), writing <seed>-<myr>.json and .log under $OUT. A run whose JSON exists is
# skipped. Each save's SHA-256 is checked first, and each run waits for compute_budget.py.
#
#   OUT=~/mantle-bloom-captures/issue-912 bin/debug/run_912_g4.sh
#
# Saves are read from $SAVES (default $MANTLE_BLOOM_SAVES_DIR, else ~/Downloads). Set
# NO_WAIT=1 to skip the compute-budget wait (an owner override).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="$ROOT/backend/.venv/bin/python"
OUT="${OUT:?set OUT to the capture directory}"
SAVES="${SAVES:-${MANTLE_BLOOM_SAVES_DIR:-$HOME/Downloads}}"
mkdir -p "$OUT"

# seed, save file, SHA-256, expected peak GB, extra arguments. Young save first.
RUNS="649658380|mantle-bloom-seed649658380-6400000y.mbworld|ebbfb5da073ae27b3baa20e53ef958cedd9153f4ce38ca549b742c561273f7b3|3|
997271774|mantle-bloom-seed997271774-87400000y.mbworld|08ad9461d9585e92be93a58d81e326c7e4b685ae7d462135d01b0e9931ed762b|6|
374464662|mantle-bloom-seed374464662-132600000y.mbworld|0f8c50c9788684c0ce28050e42a9eb1289f5092ba7fbe3a834c1809ce6a4279c|6|--ice-age-period-years 8e6"

while IFS='|' read -r seed file sha gb extra; do
  save="$SAVES/$file"
  got="$(shasum -a 256 "$save" | cut -d' ' -f1)"
  [[ "$got" == "$sha" ]] || { echo "checksum mismatch for $file: $got" >&2; exit 1; }
  for myr in 0.25 1; do
    name="$seed-$myr"
    if [[ -f "$OUT/$name.json" ]]; then
      echo "skip $name"
      continue
    fi
    if [[ -z "${NO_WAIT:-}" ]]; then
      until python3 "$ROOT/bin/compute_budget.py" check --run sim --gb "$gb" >/dev/null 2>&1; do sleep 120; done
    fi
    echo "$(date '+%F %T') start $name"
    # shellcheck disable=SC2086
    PYTHONPATH="$ROOT/backend" "$PY" "$ROOT/bin/debug/measure_flotation_debit.py" "$save" \
      --steps 10 --myr "$myr" $extra --out "$OUT/$name.json" >"$OUT/$name.log" 2>&1
    echo "$(date '+%F %T') done $name"
  done
done <<<"$RUNS"
echo "$(date '+%F %T') all runs done"
