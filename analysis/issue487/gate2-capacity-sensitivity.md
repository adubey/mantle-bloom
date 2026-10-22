# Gate 2 capacity assumption sensitivity

## Result

The exact 22,896 node loads from the 74 connected-source candidates were held fixed and re-scored.
The Gate 1 event identities and cut counts were validated against the pinned baseline; this is a
capacity-model change only, so no 1 Gyr world replay was needed.

The reference uses the steady conductive geotherm, `K0=0.8`, hydrostatic freshwater pressure
`P_f=rho_water*g*z` with `rho_water=1000 kg/m3`, `mu=0.6` (`phi=atan(mu)`), strain rate
`1e-15 s-1`, 20 MPa cohesion, plagioclase lower crust and dry olivine mantle. Only friction's
0.6–0.85 range is directly anchored in Byerlee's rock-friction experiments. The other values are
scenario choices or inherited model constants, not local observations. In particular, a fixed
`P_f/P_lith=0.42` is retained as an overpressure comparison and is no longer labeled hydrostatic.

| Capacity scenario | Nodes above capacity | Share | Changed vs reference |
|---|---:|---:|---:|
| Reference: hydrostatic water profile | 1,563 | 6.83% | — |
| Stored Moho-lag interpolation | 1,558 | 6.80% | 0.03% |
| Base temperature 1,250 C | 1,283 | 5.60% | 1.22% |
| Base temperature 1,400 C | 1,797 | 7.85% | 1.02% |
| Strain rate `1e-16 s-1` | 1,929 | 8.43% | 1.60% |
| Strain rate `1e-14 s-1` | 1,297 | 5.66% | 1.16% |
| Dry, zero pore pressure | 1,080 | 4.72% | 2.13% |
| Fixed lithostatic pore-pressure ratio 0.42 | 1,718 | 7.50% | 0.69% |
| `K0=0.6` | 2,572 | 11.23% | 4.41% |
| `K0=1.0` | 795 | 3.47% | 3.35% |
| Cohesion 1 MPa | 1,880 | 8.21% | 1.38% |
| Friction coefficient 0.85 | 1,344 | 5.87% | 0.96% |
| Mafic-granulite lower crust | 1,592 | 6.95% | 0.13% |
| Wet-olivine mantle | 2,230 | 9.74% | 2.91% |

The node ranking is highly stable: Spearman correlation of `N_ext/N_cap` against the reference is
at least 0.9936 across the one-factor variants. The threshold classification is more sensitive.
Across this finite OAT envelope, 795 nodes (3.47%) exceed capacity in every variant, and 2,701
(11.80%) exceed it in at least one. These bounds are not probabilities and omit interactions.
`K0` remains the strongest tested single-factor effect; wet mantle rheology also shifts the
exceedance fraction. The fixed 0.42 lithostatic ratio increases exceedance relative to the
hydrostatic profile, as expected for the stronger assumed pore pressure.

At the hydrostatic reference, continental, mature-oceanic and young-ridge-proxy exceedance fractions
are 1.57%, 20.32% and 18.06%, respectively. The young-ridge category still means an active
divergence proxy, not oceanic age. Reference fractions by seed are 3.18% (seed 7), 10.69% (seed 11)
and 4.64% (seed 2024). None of these are observed rupture rates: the replay candidates were
selected by the old algorithm, not independently labeled as failures.

The capacity assumptions are therefore **partly justified and numerically validated, but not
identified**. The hydrostatic pressure and Mohr-Coulomb equations pass independent analytic checks;
the shell load transfer has separate manufactured-solution checks documented in
[`gate2-load-diagnosis.md`](gate2-load-diagnosis.md). But the simulation has no spatial K0, fluid,
water, grain-size, local strain-rate or fabric fields, and its elastic shell does not solve
interplate traction continuity. The threshold comparison remains an offline scenario screen; Gate 2
remains open.

## Reproduction

The three raw node-sample files remain outside the repository under `/private/tmp/`. Re-score them
with the accepted Gate 1 identity references:

```bash
backend/.venv/bin/python bin/debug/sensitivity_gate2_capacity.py \
  --samples /private/tmp/g2-gate2-sensitivity-seed7.jsonl /private/tmp/g2-gate2-sensitivity-seed11.jsonl /private/tmp/g2-gate2-sensitivity-seed2024.jsonl \
  --replay-report /private/tmp/g2-gate2-sensitivity-seed7.json /private/tmp/g2-gate2-sensitivity-seed11.json /private/tmp/g2-gate2-sensitivity-seed2024.json \
  --identity-reference /private/tmp/issue487-g2-spatial-g1baseline-seed7.json /private/tmp/issue487-g2-spatial-g1baseline-seed11.json /private/tmp/issue487-g2-spatial-g1baseline-seed2024.json \
  --output analysis/issue487/gate2-capacity-sensitivity.json
```

The compact output is [`gate2-capacity-sensitivity.json`](gate2-capacity-sensitivity.json). The
source now distinguishes hydrostatic water pressure from a fixed lithostatic-pressure fraction and
records the equations in `gate2-column-envelope.json`.
