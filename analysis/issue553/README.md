# RFC 0014 G4: depth-cap bindings

Final diagnosis-set replay for issue #553, run after D1c (#668) and D3 were merged.
The code under measurement is the C-default production solver on commit
`70b3e904899f22f4cda2d3134ca2e8fb8e898070`. The earlier endpoint/candidate comparison
in the separate #553 worktree predates D1c and is superseded for G4 by this report.

## Method

`bin/debug/gate_lake_depth_caps.py` replays each #538 diagnosis save for ten 100 kyr steps,
with ice-age control set to `never`. The save SHA-256 values below match the preregistered
checksums. The probe wraps the production `_advance_prognostic_body` and
`_integrate_body_volume` path. A binding event is one active body in one step that produces
positive overflow while limited by its artificial lake-tier or sea-tier depth cap. Rim
spills and floor drains are excluded. Body identity is the sorted set of leaf sink node IDs;
the JSON preserves the complete identity and event details.

`surplus_m3` is the positive overflow volume routed from that capped body. It can include
water passed from upstream bodies, and the same water can bind a second cap downstream; sums
across body events are therefore gross per-body surplus, not a unique global water total.
`source_supply_m3` records direct external supply to that body's members during the active
integration call; it does not include stored water or internal transfers. All event-level
fields and per-step counts are in [`g4-final.json`](g4-final.json).

## Results

| Seed | Lake cap body-steps (steps with a bind) | Lake overflow surplus (m³) | Sea cap body-steps (steps with a bind) | Sea overflow surplus (m³) |
|---:|---:|---:|---:|---:|
| 910211954 | 34 (10/10) | 4.5591e20 | 12 (7/10) | 2.2244e20 |
| 635056485 | 2 (2/10) | 5.6545e17 | 0 (0/10) | 0 |
| 317739385 | 13 (8/10) | 1.0229e20 | 1 (1/10) | 2.8775e17 |

The largest lake-tier body-step surpluses were 2.7078e20 m³ on seed 910211954 (step 1,
1.191 million km²), 5.6536e17 m³ on seed 635056485 (step 5, 36,369 km²), and 9.9171e19 m³
on seed 317739385 (step 7, 7.077 million km²). The largest sea-tier event was 1.8372e20 m³
on seed 910211954 (step 1, 561,076 km²). The complete sink-node identities and per-body
surplus values are recorded in the JSON.

Both caps bind repeatedly on two of the three saves, and seed 910211954 has substantial
overflow under both tiers. Seed 635056485 has only two lake-cap events and no sea-cap event.
G4 is measured. The owner chose the D4 option to retain both ceilings only as numerical
guards, with an aggregate warning when a cap binds. These thresholds are not validated
physical lake-depth limits; all overflow remains booked through the hierarchy and ledger.

## Reproduction

Run from the repository root after `./bin/setup.sh`:

```bash
backend/.venv/bin/python bin/debug/gate_lake_depth_caps.py <save> --myr 1 \
    --step-years 100000 --out <scratch>/seed.json
```

| Seed | Save basename | Save SHA-256 |
|---:|---|---|
| 910211954 | `mantle-bloom-seed910211954-220200000y.mbworld` | `bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2` |
| 635056485 | `mantle-bloom-seed635056485-452800000y.mbworld` | `ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045` |
| 317739385 | `mantle-bloom-seed317739385-200800000y.mbworld` | `6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd` |

The checked-in compact result SHA-256 is `02e8b2c3b0356626be2613c1a79ec8e772b95b485c30c1e90c384d48d42320a5`.
The per-save probe output hashes, before consolidation into `g4-final.json`, are:

| Seed | Probe output SHA-256 |
|---:|---|
| 910211954 | `4737ccc57c42f9eaa6788c0cd2d16b6ec12df03a2212dbe002fceadb83ff353d` |
| 635056485 | `9b74102ea2c76ef14ee58f05f782457f8cc397d44f6977cf545bd34ce4d5be45` |
| 317739385 | `04afdb9405916ae5451fdcbd7b583c971f06d2e46a913b65b7e33edd6111f839` |
