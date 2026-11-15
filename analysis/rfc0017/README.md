# RFC 0019: atmospheric moisture diagnosis

## G0: archived G1b evidence (partial)

The owner accepted [RFC 0019](../../rfc/0019-conservative-atmospheric-moisture.md) on
2026-11-14; [#898](https://github.com/adubey/mantle-bloom/issues/898) tracks the work.
This is an audit of existing outputs, not a new simulation or a completed G0.

The `rfc0016-g1b-outputs.tar.gz` archive from #771 has SHA-256
`90277c33c631239b480fe55060cb1a27deb531ea3f851d8467002131077828de`.
All 84 JSON files pass the archive's `SHA256SUMS` checks. Their recorded source commit is
`0bda513df0c9b0baeb3cd708be142f82eb332f1d`; the original run commands, checkpoint
hashes, controls and limitations are documented in [RFC 0016 G1b](../rfc0016/README.md#g1b-result).
The raw outputs remain private.

For each window, calculate each step's central-treatment / control ratio of
`land_rain_only.rain_m3`, then average the ratios after the first step, when the
humidity feedback has not yet acted. This reproduces #771's rounded values:

| Window | Step size | Mean ratio after first step | All steps, for context |
|---|---:|---:|---:|
| A | 100 kyr | 1.3807 | 1.2855 |
| D | 10 kyr | 1.4802 | 1.4682 |
| E | 10 kyr | 1.5078 | 1.4951 |

This confirms the arithmetic behind the reported 38–51% rise conditional on excluding the
initial step. It does not independently reproduce the simulation, establish a global water
budget, or establish that the rise persists on current source. The averaging convention
matters: including the initial, uncoupled step gives the lower ratios shown above.

Still required for G0: check exact input bytes and current-source compatibility, reproduce
`measure_climate_fields.py` for the named worlds, independently reproduce G1b where its
pinned code and inputs permit, trace the complete source/sink ownership, and establish
area-weighted totals. Do not mark the gate passed on this archive audit alone.
