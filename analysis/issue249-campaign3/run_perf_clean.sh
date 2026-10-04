#!/bin/sh
# Clean #147 comparison for fa7ecdf: both surfaces, one world at a time, quiet machine.
set -u
cd "$(dirname "$0")/../../backend"
OUT=../analysis/issue249-campaign3
.venv/bin/python ../bin/debug/surface_parity.py paired --preset issue147 --seeds 0 --jobs 1 --out $OUT/issue147-clean > $OUT/issue147-clean.log 2>&1
touch $OUT/perf-clean.done
