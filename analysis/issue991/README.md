# Issue 991 · Frozen-river cleanup

RFC 0015 D5 retires the old liquid-to-glacier conversion and its zero-valued
ledger terms. Principles 3, 5, 6, and 8 govern this boundary. Streams deliver
liquid water to cold lake sinks, the lake budget closes on that delivery,
and the glacier receives ice only through its own mass balance or D4's cap
handoff. The affected tests cover older v4 and v5 saves with cached ice ledgers.

## Precipitation accounting

The fixed-phase all-land test supplies 800 mm/yr for 100 kyr to a twelve-cell
continental row. The four cold cells stay at −40 °C and the others start at
50 °C, with no seasonal swing. All rainfall enters the liquid ledger. The
sum of snow accumulated as ice and bare-ground snowmelt, converted to water
equivalent, equals the snowfall input within a relative tolerance of 10⁻⁹.

A second check isolates the limit of that result when the glacier march can
change the phase of precipitation. It uses the same twelve-cell row and
precipitation with seed 0, sea level −1,000 m, the glacier march enabled,
lake substeps disabled, and three 100 kyr intervals. Each interval carries
forward glacier depth, river water, and the glacier-water account. The first
four nodes have mean −20 °C and zero seasonal amplitude. The remaining nodes
have mean −3 °C and amplitude 10 °C. The fixture follows
`_glacier_ledger_run` in `backend/unit_tests/test_hydrology.py`, with
`sea_level_m=-1000.0` and the mixed-phase temperatures above.

On both detached `origin/main` at `9cb6d373` and the D5 edits, opening rain plus
diagnostic frozen input matches total precipitation within 52 m³ on an input
of 1.75 × 10¹⁶ m³ per interval. The glacier solver's snow input, measured as
`(accumulation_m3 + snow_melted_m3) * RHO_ICE / RHO_WATER`, differs from the
opening snowfall used by the precipitation account. The differences were
+29.8206%, −1.0004%, and −1.0004% across the three intervals on both backends.
The runoff partition reads the opening surface, while the glacier rates read
the changing trial surface. Exact physical closure across those two partitions
remains unverified for mixed-phase cells. The D5 cleanup removes the historical
cold-sink loss path without changing either partition.
