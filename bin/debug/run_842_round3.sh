#!/usr/bin/env bash
# #842 round 3: run every preregistered probe, one at a time, resumably.
#
# Implements analysis/issue842/preregistration.md's run table with probe_822_march_cost.py,
# writing <name>.json (+ .npz, .log) under $OUT. A run whose JSON exists is skipped, so the
# script can be restarted. Before each run it waits for compute_budget.py to say RUN; timed
# runs (E1, E2, E5, E6) also wait until no other run is listed, since overlapping runs'
# times aren't comparable (compute-budget skill, section 5).
#
#   OUT=~/mantle-bloom-captures/issue-842 bin/debug/run_842_round3.sh [phase ...]
#
# Phases, in order by default: e12 e3 e4 e6 e5. Saves are read from $SAVES (default
# $MANTLE_BLOOM_SAVES_DIR, else ~/Downloads) by #548's file names; gen5 and gen30 come
# from #822's capture directory ($GEN_DIR).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
PY="$ROOT/backend/.venv/bin/python"
PROBE="$ROOT/bin/debug/probe_822_march_cost.py"
OUT="${OUT:?set OUT to the capture directory}"
SAVES="${SAVES:-${MANTLE_BLOOM_SAVES_DIR:-$HOME/Downloads}}"
GEN_DIR="${GEN_DIR:-$HOME/mantle-bloom-captures/issue-822}"
GB="${GB:-3}"
mkdir -p "$OUT"

# macOS ships bash 3.2, which has no associative arrays.
save_for() {
  case "$1" in
    997271774) echo "$SAVES/mantle-bloom-seed997271774-87400000y.mbworld" ;;
    649658380) echo "$SAVES/mantle-bloom-seed649658380-6400000y.mbworld" ;;
    374464662) echo "$SAVES/mantle-bloom-seed374464662-132600000y.mbworld" ;;
    447276134) echo "$SAVES/mantle-bloom-seed447276134-116300000y.mbworld" ;;
    628977350) echo "$SAVES/mantle-bloom-seed628977350-226700000y.mbworld" ;;
    295569319) echo "$SAVES/mantle-bloom-seed295569319-180000000y.mbworld" ;;
    *) echo "no save for seed $1" >&2; return 1 ;;
  esac
}
ICE_AGES=" 374464662 447276134 "
DIAGNOSIS="997271774 649658380 374464662"
HELD_OUT="447276134 628977350 295569319"

wait_budget() {  # $1 = "timed" to also require an otherwise idle machine
  local out
  while true; do
    if out="$(python3 "$ROOT/bin/compute_budget.py" check --run sim --gb "$GB" 2>&1)"; then
      if [[ "$1" != timed ]] || grep -q "nothing else is running" <<<"$out"; then
        return
      fi
    fi
    sleep 120
  done
}

probe() {  # probe <name> <timed|untimed> <save path> [probe args...]
  local name="$1" timed="$2" save="$3"
  shift 3
  if [[ -f "$OUT/$name.json" ]]; then
    echo "skip $name"
    return
  fi
  wait_budget "$timed"
  echo "$(date '+%F %T') start $name"
  PYTHONPATH="$ROOT/backend" "$PY" "$PROBE" "$save" --out "$OUT/$name.json" "$@" >"$OUT/$name.log" 2>&1
  echo "$(date '+%F %T') done $name"
}

ice_args() { [[ "$ICE_AGES" == *" $1 "* ]] && echo "--ice-age-period-years 8e6" || true; }

replays() {  # 10 Myr at both step sizes, both arms
  local seed arm
  for seed in "$@"; do
    for arm in d2 trend-1k; do
      # shellcheck disable=SC2046
      probe "$seed-$arm-0.25" timed "$(save_for "$seed")" --preset "$arm" --steps 40 --myr 0.25 $(ice_args "$seed")
      # shellcheck disable=SC2046
      probe "$seed-$arm-1" timed "$(save_for "$seed")" --preset "$arm" --steps 10 --myr 1 $(ice_args "$seed")
    done
  done
}

e3_steps() {  # the step ending at 97.4 Myr (opens after 39 steps) and the largest loss after 2 Myr
  "$PY" - "$OUT/997271774-d2-0.25.json" <<'EOF'
import json, sys
steps = json.load(open(sys.argv[1]))["per_step"]
vols = [s["ice_closing_km3"] for s in steps]
losses = sorted(((vols[i - 1] - vols[i]) / vols[i - 1], i) for i in range(8, len(vols)))
pick = next(i for loss, i in reversed(losses) if i != 39)
print(39, pick)
EOF
}

e3() {
  local k base arm scale extra
  for k in $(e3_steps); do
    base="$OUT/e3-open-k$k.mbworld"
    if [[ ! -f "$base" ]]; then
      probe "e3-rerun-k$k" untimed "$(save_for 997271774)" --steps "$k" --myr 0.25 \
        --write-save-after "$k" --write-save "$base"
    fi
    for arm in trend-1k d2 d2-10k; do
      extra=(--preset "${arm%-10k}")
      [[ "$arm" == d2-10k ]] && extra+=(--substep-years 10000)
      for scale in 1.0 0.995 0.998 0.999 0.9995 0.9998 0.9999 1.0001 1.0002 1.0005 1.001 1.002 1.005; do
        probe "e3-k$k-$arm-$scale" untimed "$base" "${extra[@]}" --steps 1 --myr 0.25 --scale-opening-ice "$scale"
      done
    done
  done
}

e4() {
  local seed sub myr
  for seed in $DIAGNOSIS; do
    for myr in 0.25 1; do
      for sub in 1000 500 250 125; do
        # shellcheck disable=SC2046
        probe "e4-$seed-sub$sub-$myr" untimed "$(save_for "$seed")" --preset d2 --substep-years "$sub" \
          --steps 1 --myr "$myr" $(ice_args "$seed")
      done
    done
  done
}

e6() {
  local arm
  for arm in d2 trend-1k; do
    probe "e6-gen30-$arm" timed "$GEN_DIR/gen30-p8.mbworld" --preset "$arm" --steps 1 --myr 1
    probe "e6-gen5-$arm" timed "$GEN_DIR/gen5-d025.mbworld" --preset "$arm" --steps 1 --myr 0.25
  done
}

phases=("$@")
[[ ${#phases[@]} -eq 0 ]] && phases=(e12 e3 e4 e6 e5)
for phase in "${phases[@]}"; do
  case "$phase" in
    e12) replays $DIAGNOSIS ;;
    e3) e3 ;;
    e4) e4 ;;
    e6) e6 ;;
    e5) replays $HELD_OUT ;;
    *) echo "unknown phase $phase" >&2; exit 2 ;;
  esac
done
echo "$(date '+%F %T') all phases done"
