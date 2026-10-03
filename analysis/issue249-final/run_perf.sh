#!/bin/sh
# #147 profile comparison: sequential (--jobs 1), run on a quiet machine.
set -u
cd "$(dirname "$0")/../../backend"
OUT=../analysis/issue249-final
.venv/bin/python ../bin/debug/surface_parity.py paired --preset issue147 --seeds 0 --jobs 1 --out $OUT/issue147 > $OUT/issue147.log 2>&1
touch $OUT/perf.done
