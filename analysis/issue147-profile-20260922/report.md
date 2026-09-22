# Issue #147: 60-frame, 60-million-year animation profile

## Run and method

- Clean `main` at commit `3dc562bfdb9718437707f64e7d469a3d5d0d5a8c`; fault-cache PR #215 is excluded. Python 3.14.6, Darwin 25.6.0 arm64.
- Seed 0; node and climate density 4.0; fluid density 2.0; diagnostic wind; Elevation & Biome view (`combined`), Eckert IV, 2200 × 1222 pixels.
- Exactly 60 frames, each with 10 `step_world` calls at 100,000 years, followed by the normal PNG render and H.264 encoding calls. Frame 1 advances before rendering; frame 60 ends at **60 Myr**. The app's normal animation function renders frame 1 before stepping, so its 60-frame run would end at 59 Myr.
- `cProfile` covered all 600 steps and 60 renders. Timers around major functions supplied per-frame wall seconds. These timings include profiler overhead and any host contention; they are not an unprofiled runtime estimate.
- World generation took 1.60 seconds, outside the frame totals. Plates grew from 19 to 36 and nodes from 130,587 to 139,705. The completed MP4 is 4,457,105 bytes.

## Headline

| Phase | Total seconds | Share of run | Mean seconds/frame |
| --- | ---: | ---: | ---: |
| Scheduled transport and other step work | 2435.7 | 38.8% | 40.60 |
| Plate deformation | 1547.8 | 24.6% | 25.80 |
| Fault processing | 694.1 | 11.0% | 11.57 |
| Climate, erosion, hydrology | 644.0 | 10.2% | 10.73 |
| Plate shifts | 485.8 | 7.7% | 8.10 |
| Sea level | 144.2 | 2.3% | 2.40 |
| Rendering | 123.5 | 2.0% | 2.06 |
| Topology changes | 120.8 | 1.9% | 2.01 |
| Volcanism and resource formation | 85.0 | 1.4% | 1.42 |
| PNG decode and H.264 encode | 3.0 | 0.0% | 0.05 |

**Total:** 6,285.1 seconds (104.8 minutes): 6,157.4 seconds stepping, 123.5 seconds rendering, and roughly 4.2 seconds encoding and other work. The frame timers sum to 6,284.1 seconds; encoder flush and bookkeeping account for roughly 1.0 second outside frame intervals.

### Low-level hotspots

`cProfile` cumulative times include child calls and cannot be added together.

| Function | Calls | Cumulative seconds | Own seconds | Meaning |
| --- | ---: | ---: | ---: | --- |
| `magma_transport.run_magma_transport` | 150 | 2,339.8 | 146.6 | Scheduled every four steps; 37.2% of the entire run and about 96% of the high-level scheduled/other bucket. |
| `magma_transport._weighted_destination_pairs` | 150 | 2,149.5 | 1,286.0 | Dominant subroutine; builds and weights parcel-to-continental-node candidate pairs. |
| `LithospherePlate.deform` | 18,320 | 1,547.1 | 132.2 | Second largest independent high-level phase. |
| `faults.update_faults` | 600 | 694.1 | 1.4 | Includes boundary-fault generation (~210.2 s) and fault relief (~210.9 s). |
| `erosion.apply_erosion` | 600 | 644.0 | 10.6 | Includes climate (~249.9 s) and hydrology (~142.6 s). |
| `LithospherePlate.shift` | 18,320 | 485.7 | 0.2 | Plate motion and torque work. |
| `faults.fault_tangent_components` | 364,252 | 347.5 | 244.5 | Repeated tangent geometry during rift stretching; PR #215 targets only the candidate-list filter. |
| `eustasy.update_sea_level` | 600 | 144.2 | 0.0 | Includes 17,587 seeded connectivity checks (~134.3 s). |
| `render_png` | 60 | 123.5 | 0.2 | About 2% of the run; biome soft blending takes ~69.9 s. |

`threading.join` shows 1,086.5 cumulative seconds over 2.58 million calls. These waits occur inside parallel SciPy KD-tree queries and overlap worker activity, so this is not an additional serial cost.

### Time trend

- Frames 1–10 averaged 52.8 seconds/frame; frames 21–30 averaged 104.2; frames 31–40 averaged 112.4. Growth in scheduled transport and maintenance accounts for most of the change.
- Magma transport and gap work run every four steps: three times in each even frame and twice in each odd frame. The scheduled/other bucket averaged 47.9 seconds in even frames and 33.3 seconds in odd frames.
- Frames 55–58 were unusually slow across plate movement, faults, climate, and rendering at once; frame 55 reached 380.1 seconds. This broad inflation is consistent with host contention, but host load was not recorded, so the cause is unproven. The raw values are retained below.

## Every frame, seconds

Each frame advances one million simulated years. `Motion` is plate shift plus deformation. `Climate/erosion` includes hydrology. `Geo/sea` combines volcanism, resource formation, and sea-level solving. `Scheduled/other` includes magma transport, gap work, reconciliation, statistics, and step overhead. The low-level profile attributes 2,339.8 of its 2,435.7 total seconds to magma transport; its exact per-frame split was not timed. `Render/encode` includes PNG render, decode, and H.264 encode.

| Frame / end Myr | Plates | Motion | Faults | Topology | Climate/erosion | Geo/sea | Scheduled/other | Render/encode | Total |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 19 | 15.6 | 6.8 | 0.3 | 7.7 | 3.2 | 3.4 | 2.5 | **39.4** |
| 2 | 20 | 16.9 | 6.8 | 0.6 | 8.0 | 2.8 | 7.5 | 1.9 | **44.4** |
| 3 | 24 | 18.4 | 7.8 | 0.5 | 8.1 | 3.0 | 6.7 | 1.9 | **46.4** |
| 4 | 24 | 18.6 | 8.1 | 0.8 | 8.3 | 2.9 | 11.5 | 1.8 | **52.0** |
| 5 | 25 | 19.4 | 8.1 | 0.7 | 8.2 | 3.0 | 8.2 | 1.9 | **49.6** |
| 6 | 25 | 19.4 | 8.1 | 1.0 | 8.3 | 3.1 | 14.0 | 1.8 | **55.8** |
| 7 | 27 | 20.7 | 8.7 | 0.9 | 8.6 | 3.1 | 11.1 | 1.9 | **54.9** |
| 8 | 27 | 20.5 | 8.6 | 1.3 | 8.3 | 3.0 | 16.9 | 1.7 | **60.4** |
| 9 | 28 | 21.5 | 8.7 | 1.2 | 8.4 | 3.2 | 13.4 | 1.7 | **58.1** |
| 10 | 28 | 21.7 | 8.7 | 1.7 | 8.5 | 3.2 | 21.1 | 1.7 | **66.5** |
| 11 | 28 | 22.3 | 8.7 | 1.4 | 8.4 | 3.2 | 16.1 | 1.7 | **61.7** |
| 12 | 28 | 22.3 | 9.0 | 1.8 | 8.6 | 3.2 | 25.7 | 1.7 | **72.3** |
| 13 | 28 | 23.1 | 8.9 | 1.4 | 8.6 | 3.3 | 18.7 | 1.7 | **65.8** |
| 14 | 28 | 22.5 | 8.9 | 1.8 | 8.6 | 3.2 | 29.5 | 1.7 | **76.3** |
| 15 | 28 | 23.5 | 8.9 | 1.4 | 8.7 | 3.3 | 21.9 | 1.8 | **69.4** |
| 16 | 29 | 23.3 | 9.0 | 1.7 | 8.6 | 3.2 | 36.5 | 1.7 | **84.0** |
| 17 | 29 | 24.7 | 9.2 | 1.3 | 8.7 | 3.3 | 27.1 | 1.8 | **76.2** |
| 18 | 30 | 24.7 | 9.4 | 1.8 | 8.7 | 3.3 | 47.6 | 1.7 | **97.3** |
| 19 | 31 | 25.8 | 9.6 | 1.5 | 8.9 | 3.4 | 34.2 | 1.8 | **85.3** |
| 20 | 31 | 25.7 | 9.8 | 1.9 | 8.9 | 3.3 | 53.3 | 1.7 | **104.8** |
| 21 | 31 | 26.9 | 9.9 | 1.6 | 8.9 | 3.6 | 36.9 | 1.8 | **89.6** |
| 22 | 31 | 26.2 | 9.9 | 2.1 | 9.1 | 3.4 | 54.8 | 1.8 | **107.2** |
| 23 | 31 | 27.0 | 9.9 | 1.7 | 8.9 | 3.5 | 35.7 | 1.7 | **88.3** |
| 24 | 31 | 26.2 | 9.9 | 2.0 | 9.2 | 3.4 | 57.8 | 1.7 | **110.2** |
| 25 | 31 | 26.7 | 9.9 | 1.6 | 8.9 | 3.5 | 41.0 | 1.7 | **93.2** |
| 26 | 31 | 25.8 | 10.1 | 2.1 | 9.0 | 3.3 | 63.2 | 1.7 | **115.2** |
| 27 | 31 | 26.8 | 9.9 | 1.7 | 9.0 | 3.4 | 43.0 | 1.7 | **95.6** |
| 28 | 31 | 25.8 | 10.0 | 2.2 | 9.0 | 3.4 | 66.0 | 1.7 | **118.0** |
| 29 | 31 | 26.7 | 10.2 | 1.8 | 9.0 | 3.4 | 46.2 | 1.7 | **99.0** |
| 30 | 31 | 26.0 | 10.0 | 2.2 | 9.1 | 3.3 | 73.6 | 1.7 | **126.0** |
| 31 | 31 | 27.2 | 10.0 | 1.8 | 9.2 | 3.4 | 51.1 | 1.7 | **104.4** |
| 32 | 31 | 26.2 | 10.0 | 2.3 | 9.2 | 3.3 | 69.5 | 1.7 | **122.1** |
| 33 | 31 | 27.6 | 9.9 | 1.8 | 9.1 | 3.4 | 48.6 | 1.7 | **102.1** |
| 34 | 31 | 26.5 | 9.9 | 2.3 | 9.2 | 3.4 | 62.4 | 1.7 | **115.2** |
| 35 | 31 | 27.7 | 10.0 | 2.0 | 9.4 | 3.5 | 41.6 | 1.7 | **95.8** |
| 36 | 31 | 27.0 | 10.2 | 2.3 | 9.3 | 3.5 | 66.2 | 1.7 | **120.2** |
| 37 | 31 | 28.5 | 10.1 | 1.9 | 9.3 | 3.5 | 44.8 | 1.7 | **99.8** |
| 38 | 31 | 27.7 | 10.2 | 2.4 | 9.6 | 3.5 | 73.9 | 1.7 | **129.1** |
| 39 | 31 | 28.6 | 10.4 | 1.9 | 9.4 | 3.4 | 49.9 | 1.7 | **105.4** |
| 40 | 31 | 27.3 | 10.2 | 2.4 | 9.5 | 3.5 | 74.9 | 1.7 | **129.6** |
| 41 | 31 | 28.5 | 10.1 | 1.9 | 9.5 | 3.5 | 45.7 | 1.7 | **100.9** |
| 42 | 32 | 28.5 | 10.5 | 2.4 | 9.8 | 3.4 | 59.7 | 1.7 | **116.2** |
| 43 | 32 | 29.3 | 10.3 | 1.8 | 9.8 | 3.6 | 37.5 | 1.7 | **94.0** |
| 44 | 32 | 28.4 | 10.5 | 2.4 | 9.7 | 3.5 | 53.0 | 1.7 | **109.1** |
| 45 | 32 | 29.6 | 10.3 | 1.9 | 9.6 | 3.6 | 35.1 | 1.7 | **91.8** |
| 46 | 32 | 28.7 | 10.9 | 2.3 | 9.8 | 3.5 | 49.4 | 1.7 | **106.3** |
| 47 | 32 | 33.5 | 11.7 | 1.9 | 12.2 | 4.1 | 39.7 | 2.8 | **105.9** |
| 48 | 34 | 49.1 | 16.5 | 3.1 | 15.3 | 5.1 | 72.0 | 3.0 | **164.0** |
| 49 | 34 | 43.3 | 14.1 | 2.0 | 12.7 | 4.5 | 38.5 | 2.2 | **117.3** |
| 50 | 34 | 32.3 | 12.3 | 2.2 | 10.8 | 3.7 | 39.2 | 2.0 | **102.6** |
| 51 | 34 | 32.0 | 11.5 | 1.8 | 10.0 | 3.8 | 22.1 | 1.7 | **82.8** |
| 52 | 36 | 30.1 | 11.1 | 2.0 | 9.8 | 3.6 | 32.9 | 1.8 | **91.4** |
| 53 | 36 | 32.1 | 11.6 | 1.7 | 10.3 | 3.7 | 22.0 | 2.1 | **83.5** |
| 54 | 36 | 31.3 | 11.6 | 2.2 | 10.7 | 4.2 | 33.7 | 6.4 | **100.1** |
| 55 | 36 | 181.2 | 47.9 | 7.7 | 39.6 | 12.3 | 83.7 | 7.7 | **380.1** |
| 56 | 36 | 152.7 | 34.0 | 5.3 | 30.4 | 9.1 | 81.4 | 4.6 | **317.5** |
| 57 | 36 | 109.2 | 28.1 | 4.0 | 26.0 | 7.4 | 50.1 | 4.7 | **229.4** |
| 58 | 36 | 85.2 | 22.0 | 4.2 | 20.7 | 6.4 | 56.4 | 2.7 | **197.6** |
| 59 | 36 | 42.9 | 14.5 | 2.3 | 12.8 | 4.8 | 25.2 | 2.4 | **104.9** |
| 60 | 36 | 36.3 | 12.4 | 2.6 | 12.0 | 3.9 | 32.9 | 2.0 | **102.1** |

## Artifacts

- [Full per-frame CSV](frames.csv): every measured major function, node and plate counts, full precision seconds.
- [Cumulative profile](profile-cumulative.txt) and [internal-time profile](profile-internal.txt): top 120 functions.
- `profile.pstats`: complete machine-readable profile (local worktree only).
- `animation.mp4`: completed rendered video (local worktree only).
- `metadata.json`: settings and run totals.
- [Profiling driver](../../bin/debug/profile_animation_60m.py): code used for this run.
