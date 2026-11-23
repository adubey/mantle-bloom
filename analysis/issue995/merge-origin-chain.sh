#!/bin/bash
# Baseline rerun with merge and daughter-id recording (#995): which big plates split soon after becoming big, and how they formed.
D=$HOME/mantle-bloom-captures/issue995/merge-origin
mkdir -p $D
for seed in 591453549 106533960 1097624328 390382636 1760721944; do
  out=analysis/issue995/rift-merge-$seed.json
  [ -f $out ] || python3 bin/debug/big_plate_rift_probe.py --settings analysis/issue995/settings.json --setting brake-o0.25-c1 --seed $seed --years 200 --every 5 --scratch $D --out $out > $D/$seed.log 2>&1
done
echo finished
