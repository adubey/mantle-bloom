# RFC 0015: how saved worlds represent frozen water

Evidence for [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md), Problem items 5,
7 and 8. [`bin/debug/measure_frozen_water.py`](../../bin/debug/measure_frozen_water.py)
reads a save without stepping it. It measures land ice against the Biome view's
classification, the precipitation hydrology counts as snow, and the ice that sea level debits.
The numbers are in [`frozen-water.json`](frozen-water.json).

Run on `origin/main` `9d6fa07`, one save at a time (peak RSS 2.2 and 2.5 GB):

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

Both saves are local and not available for download. The map measures use the Biome view's
own grid, area-weighted by cos(latitude). Months above 0 °C come from the model's sinusoidal
annual cycle (`biomes.seasonal_temp_amplitude`), with precipitation assumed even through the
year. Neither save has ice-age forcing, so the ice-age sea-ice cap is absent here.
