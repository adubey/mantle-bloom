# RFC 0008 evidence: craton sinks on current main

Background numbers for [RFC 0008](../../rfc/0008-craton-collision-jamming.md) and the 2026-10-14
revision of [RFC 0002](../../rfc/0002-craton-trench-attrition.md). Not the RFC's gate, which needs
its own probe.

`bin/debug/compare_cratons.py --seeds 7,11,2024 --variants on --years 1e9 --every 100e6 --workers 1`
on origin/main at 1b4f686, default density, 5 Myr steps. Numbers are in
[craton_sinks_1gyr.json](craton_sinks_1gyr.json).

This baseline includes #358, which counts the split cooldown in years (2 Myr) instead of steps.
At 5 Myr steps the old cooldown was 20 steps, 100 Myr, so these runs split about 3.5 times as
often as the runs behind #378, #385 and #388 (104–115 splits and 27–55 merges per Gyr, against
32–37 splits and 9–15 merges). Numbers from before #358 at 5 Myr steps aren't comparable.
