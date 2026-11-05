# Issue 552: ice in an enclosed basin below sea level

Fix: PR #812. Water floats ice only where `is_ocean` (`lithosphere.floating_water_depth_m`).

## Replay

`bin/debug/measure_enclosed_calving.py`, run 2 Myr in 100 kyr steps from one ice-age save
(seed 997271774 at 87.4 Myr, RFC 0015's diagnosis set), on `origin/main` (`d02484e6`) and on the
branch. `replay.json` has the compact results, including the save's SHA-256. Windows are 0.5 Myr.

| Window end (Myr) | Calved main (m³) | Calved branch (m³) | Off-ocean calved, main (m³) | Off-ocean share, main | Ice volume main / branch (m³) | Sea level main / branch (m) |
|---|---:|---:|---:|---:|---|---|
| 87.9 | 2.07e10 | 2.03e10 | 3.91e8 | 1.9% | 5.662e15 / 5.662e15 | -411.1 / -411.1 |
| 88.4 | 2.23e11 | 0 | 0 | 0 | 5.698e15 / 5.698e15 | -418.7 / -418.7 |
| 88.9 | 6.33e11 | 2.00e10 | 0 | 0 | 5.712e15 / 5.712e15 | -427.2 / -427.4 |
| 89.4 | 1.13e11 | 2.10e11 | 0 | 0 | 5.740e15 / 5.739e15 | -435.0 / -434.9 |

Branch off-ocean calving is 0 in every window. The ice-ledger residual stays at 1-2 m³ on both.

## Findings

- On this save the bug is small. `origin/main` calved 3.9e8 m³ from non-ocean nodes, once, in the
  first window (at most 2 nodes), which is under 1e-7 of the ice volume. After that, none.
- Ice volume and sea level agree to better than 0.2 m. The fix moves neither on this save.
- Total calved volume differs between the runs after the first window, in both directions. The
  off-ocean calving on `origin/main` is zero in those windows, so the fix does not explain
  it. The states diverge from step 1, and a calved volume of 1e11 m³ is 2e-5 of the ice volume. I haven't
  separated that divergence from noise, for example with a second baseline run.
- Not measured: enclosed pits holding ice, and their depth. That is issue #813.

## Limits

One save, one ice-age window of 2 Myr. A cold, pit-rich world could show more. The off-ocean
share uses the end-of-step `is_ocean`.
