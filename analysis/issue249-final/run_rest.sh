#!/bin/sh
# Second half of run_correctness.sh, restarted after the first attempt was cut off.
set -u
cd "$(dirname "$0")/../../backend"
OUT=../analysis/issue249-final
P=".venv/bin/python ../bin/debug/surface_parity.py"
$P paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render --out $OUT/repro-equivalent > $OUT/repro-equivalent.log 2>&1 &
$P paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2 --out $OUT/stress > $OUT/stress.log 2>&1 &
wait
touch $OUT/rest.done
