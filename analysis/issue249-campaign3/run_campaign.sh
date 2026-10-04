#!/bin/sh
# Third #249 campaign, on main after the #268 seam-overlap fix (PR #270).
# Correctness sets run in parallel; the issue147 timing run goes last, alone.
# SAVE: the real 352.4 Myr line save (not in the repo). CONVERTED: scratch dir for its quad conversion.
set -u
SAVE=${SAVE:-$HOME/Downloads/mantle-bloom-seed804913535-352400000y.mbworld}
CONVERTED=${CONVERTED:?set CONVERTED to a scratch directory}
cd "$(dirname "$0")/../../backend"
OUT=../analysis/issue249-campaign3
P=".venv/bin/python ../bin/debug/surface_parity.py"

$P paired --preset long --seeds 1,2,3,4,5 --jobs 5 --out $OUT/long > $OUT/long.log 2>&1
$P paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render --out $OUT/repro-equivalent > $OUT/repro-equivalent.log 2>&1 &
$P paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2 --out $OUT/stress > $OUT/stress.log 2>&1 &
wait

# The 352.4 Myr reproducer on both surfaces: the line save as is, and its #248 quad conversion.
.venv/bin/python ../bin/debug/convert_legacy_saves.py "$SAVE" --out "$CONVERTED" --steps 0 --write-converted > $OUT/convert.log 2>&1
QSAVE="$CONVERTED/$(basename "$SAVE" .mbworld)-quad.mbworld"
$P run --preset issue147 --seed 804913535 --surface lines --from-world "$SAVE" --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render --out $OUT/repro-saves > $OUT/repro-lines.log 2>&1
$P run --preset issue147 --seed 804913535 --surface quad --from-world "$QSAVE" --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render --out $OUT/repro-saves > $OUT/repro-quad.log 2>&1
$P compare $OUT/repro-saves > /dev/null 2>&1

$P paired --preset issue147 --seeds 0 --jobs 1 --out $OUT/issue147 > $OUT/issue147.log 2>&1
touch $OUT/campaign.done
