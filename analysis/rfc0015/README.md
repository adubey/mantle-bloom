# RFC 0015: how saved worlds represent frozen water

Evidence for [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md), Problem items 5,
8, 9 and 10, and D6's sea-ice scale check. [`bin/debug/measure_frozen_water.py`](../../bin/debug/measure_frozen_water.py)
reads a save without stepping it. It measures land ice against the Biome view's
classification, the precipitation hydrology counts as snow, and the ice that sea level debits.
The numbers are in [`frozen-water.json`](frozen-water.json).

Run on `8e2c1ba` (no model changes from `origin/main`), one save at a time (peak RSS
2.2 and 2.5 GB):

```bash
cd backend
.venv/bin/python ../bin/debug/measure_frozen_water.py <local-save> --out <scratch-output>
```

| Measure | Seed 997271774, 87.4 Myr (`08ad9461…ed762b`) | Seed 649658380, 6.4 Myr (`ebbfb5da…73f7b3`) |
|---|---:|---:|
| Land under visible ice (≥ 10 m) | 6.0% | 2.8% |
| …classified EF (Ice Cap) | 93% | 33% |
| …classified ET (Tundra) | 5% | 60% |
| …warmest month ≥ 0 °C at the bed but < 0 °C at the ice surface | 3% | 37% |
| Mean ice thickness on iced land | 391 m | 853 m |
| EF land without ice | 39% | 3% |
| Land with an annual mean below 0 °C | 27% | 13% |
| …its precipitation that falls in months above 0 °C | 28% | 44% |
| Ice debited from the ocean (km³ water equivalent) | 5,184,194 | 3,267,068 |
| …floating | 0.4 | 71 |
| …grounded below sea level | 3,551 | 2,803 |
| Ocean poleward of 60°, share of ocean | 14% | 14% |
| …lowest sea-surface temperature | +1.2 °C | +1.2 °C |
| …at the −2.0 °C floor, or below −1.8 °C in the annual mean or warmest month | 0% | 0% |
| …mean air temperature over it | 8.4 °C | 8.7 °C |
| …coldest month below −1.8 °C | 52% | 49% |
| …warmest month below −1.8 °C | < 0.01% | < 0.01% |

**Sea-ice balance by latitude** (RFC 0015 D6, `sea_ice_balance`; both saves are at Earth
tilt, so both give the same rows). `T_si` is the land-temperature construct at sea level,
FDD and TDD are degree-days below −1.8 °C and above 0 °C, growth is Stefan's law, and the
state uses an ocean heat flux of 2 W m⁻²:

| Latitude | Ocean baseline | `T_si` | FDD | TDD | Growth (m/yr) | State, `DDF_ice` 4 or 7 mm |
|---:|---:|---:|---:|---:|---:|---|
| 50° | 16.7 °C | 6.3 °C | 0 | 2,333 | 0 | ice-free |
| 55° | 14.7 °C | 2.6 °C | 289 | 1,472 | 0.58 | seasonal |
| 60° | 12.5 °C | −1.3 °C | 979 | 832 | 1.08 | seasonal |
| 65° | 10.3 °C | −5.5 °C | 1,926 | 359 | 1.51 | seasonal |
| 70° | 8.1 °C | −9.6 °C | 3,045 | 68 | 1.90 | perennial |
| 80° | 4.2 °C | −16.7 °C | 5,437 | 0 | 2.54 | perennial |
| 89° | 1.4 °C | −21.7 °C | 7,283 | 0 | 2.94 | perennial |

Both saves are local and not available for download. The map measures use the Biome view's
own grid, area-weighted by cos(latitude). Months above 0 °C come from the model's sinusoidal
annual cycle (`biomes.seasonal_temp_amplitude`), with precipitation assumed even through the
year. Neither save has ice-age forcing, so the ice-age sea-ice cap is absent here. Polar-ocean
shares use the same grid; the air temperature over the ocean is the ocean's own in both
non-CFD climate modes.
