#!/bin/bash
# Time-dependent strength ramp arms (#995): label early late tau_myr, on the five baseline seeds.
D=$HOME/mantle-bloom-captures/issue995/ramp
mkdir -p $D
for arm in "e3-l0.5-t50 3 0.5 50" "e3-l0.5-t100 3 0.5 100" "e1-l0.5-t50 1 0.5 50"; do
  set -- $arm
  for seed in 591453549 106533960 1097624328 390382636 1760721944; do
    out=analysis/issue995/rift-state-ramp-$1-$seed.json
    [ -f $out ] || python3 bin/debug/big_plate_rift_probe.py --settings analysis/issue995/settings.json --setting brake-o0.25-c1 --seed $seed --years 200 --every 5 --ramp $2 $3 $4 --scratch $D --out $out > $D/$1-$seed.log 2>&1
  done
done
echo finished
