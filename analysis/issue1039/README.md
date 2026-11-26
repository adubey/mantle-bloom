# Issue #1039: Rust Rendering & Pipelined Animation Benchmarks

Measurements for [#1039](https://github.com/adubey/mantle-bloom/issues/1039) comparing the Python/Numba rendering path (`MANTLE_BLOOM_DISABLE_RUST=1`) against the Rust `mantle_bloom_rust.render` path (`combined`, `biome`, `elevation`, plus direct in-memory `globe` resampling) and pipelined animation rendering (`MANTLE_BLOOM_RENDER_THREADS=1|2`).

## Methodology & Host Conditions

- **Command:** `backend/.venv/bin/python bin/debug/bench_render_rust.py --densities 1,4 --steps-per-frame 1,10 --num-frames 4` (commit/branch `gemini/issue-1039-rust-parallel-render`).
- **Isolation:** Each arm runs in a fresh subprocess so `resource.getrusage(RUSAGE_SELF).ru_maxrss` records uncontaminated per-arm peak RSS (MB) and cold Numba/Rust initialization cost.
- **World setup:** `seed=0`, `node_density=climate_density=D` (`D in {1, 4}`), `fluid_density=1` (`D=1`) or `2` (`D=4`), `2` warmup steps of `100,000` yr, rendered at `2200x1222`. Animation runs use `behrmann/combined`, `step_years=100,000`, `num_frames=4`.
- **Concurrent host load (`compute-budget` Section 5):** Run under user override on a 16 GB, 10-core Apple Silicon Mac (`7.0 GB` available, memory pressure `2`, load `4.0` at start). One background replay (`pid 31294`, `issue1019`, `0.2–0.4 GB` RSS) overlapped earlier arms and exited before the final `D=4, steps_per_frame=10` parallel arm finished. Peak RSS (`ru_maxrss`) is unaffected by concurrent processes; wall-clock ratios should be read with that background load noted.

## Single-Frame Render Performance (`2200x1222`, `seed=0`)

| Density | Case | Python Cold (s) | Rust (2t) Cold (s) | Cold Speedup | Python Warm Mean (s) | Rust (2t) Warm Mean (s) | Warm Speedup | Python Peak RSS (MB) | Rust (2t) Peak RSS (MB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 (`1/1/1`) | `behrmann/combined` | 0.7517 | 0.1859 | 4.04x | 0.3106 | 0.1728 | 1.80x | 823.2 | 769.3 |
| 1 (`1/1/1`) | `behrmann/biome` | 0.2625 | 0.1420 | 1.85x | 0.2958 | 0.1385 | 2.14x | 823.2 | 769.3 |
| 1 (`1/1/1`) | `behrmann/elevation` | 0.2872 | 0.1427 | 2.01x | 0.2855 | 0.1411 | 2.02x | 823.2 | 769.3 |
| 1 (`1/1/1`) | `globe/combined` | 0.8993 | 0.3259 | 2.76x | 0.7809 | 0.3316 | 2.35x | 823.2 | 769.3 |
| 4 (`4/4/2`) | `behrmann/combined` | 2.6063 | 1.1568 | 2.25x | 1.8348 | 0.9754 | 1.88x | 2782.6 | 1138.8 |
| 4 (`4/4/2`) | `behrmann/biome` | 2.0141 | 0.8644 | 2.33x | 1.6392 | 0.8567 | 1.91x | 2782.6 | 1138.8 |
| 4 (`4/4/2`) | `behrmann/elevation` | 0.4327 | 0.2901 | 1.49x | 0.4248 | 0.2889 | 1.47x | 2782.6 | 1138.8 |
| 4 (`4/4/2`) | `globe/combined` | 2.6416 | 1.1677 | 2.26x | 2.6226 | 1.1647 | 2.25x | 2782.6 | 1138.8 |

Key observations:
- At Density 4, Rust cuts single-frame peak RSS from **2,782.6 MB to 1,138.8 MB (-1,643.8 MB / -59%)**, driven by eliminating the intermediate `(N, 41)` `float64` biome-softening tensor (~420 MB+ plus NumPy/Numba allocator overhead) in favor of an inline `[f64; 41]` stack accumulator, plus avoiding Numba JIT compilation for render kernels.
- Cold-frame latency on `behrmann/combined` drops from `0.75s` to `0.19s` at Density 1 (**4.04x**) and `2.61s` to `1.16s` at Density 4 (**2.25x**); warm renders are **1.47x–2.35x** faster across all tested views.

## Animation Loop Performance (`behrmann/combined`, `2200x1222`, `num_frames=4`)

| Density | `steps_per_frame` | Mode (`label`) | Total (s) | Per-Frame (s) | Speedup vs Python | Peak RSS (MB) | RSS Delta vs `py_seq` (MB) |
|---|---:|---|---:|---:|---:|---:|---:|
| 1 (`1/1/1`) | 1 | Python sequential (`py_seq`) | 38.4929 | 9.6232 | 1.00x | 607.0 | — |
| 1 (`1/1/1`) | 1 | Rust sequential, 2 threads (`rust_seq_2t`) | 18.7697 | 4.6924 | 2.05x | 698.1 | +91.1 |
| 1 (`1/1/1`) | 1 | Rust pipelined, 2 threads (`rust_par_2t`) | 17.3525 | 4.3381 | 2.22x | 612.0 | +5.0 |
| 1 (`1/1/1`) | 1 | Rust pipelined, 1 thread (`rust_par_1t`) | 17.8339 | 4.4585 | 2.16x | 675.1 | +68.1 |
| 1 (`1/1/1`) | 10 | Python sequential (`py_seq`) | 360.2819 | 90.0705 | 1.00x | 635.4 | — |
| 1 (`1/1/1`) | 10 | Rust sequential, 2 threads (`rust_seq_2t`) | 165.3427 | 41.3357 | 2.18x | 675.5 | +40.1 |
| 1 (`1/1/1`) | 10 | Rust pipelined, 2 threads (`rust_par_2t`) | 164.0550 | 41.0137 | 2.20x | 644.7 | +9.3 |
| 1 (`1/1/1`) | 10 | Rust pipelined, 1 thread (`rust_par_1t`) | 166.5115 | 41.6279 | 2.16x | 623.3 | -12.1 |
| 4 (`4/4/2`) | 1 | Python sequential (`py_seq`) | 105.2836 | 26.3209 | 1.00x | 2668.9 | — |
| 4 (`4/4/2`) | 1 | Rust sequential, 2 threads (`rust_seq_2t`) | 46.4333 | 11.6083 | 2.27x | 1213.8 | -1455.1 |
| 4 (`4/4/2`) | 1 | Rust pipelined, 2 threads (`rust_par_2t`) | 44.3797 | 11.0949 | 2.37x | 1248.5 | -1420.4 |
| 4 (`4/4/2`) | 1 | Rust pipelined, 1 thread (`rust_par_1t`) | 45.7839 | 11.4460 | 2.30x | 1303.9 | -1365.0 |
| 4 (`4/4/2`) | 10 | Python sequential (`py_seq`) | 958.6926 | 239.6732 | 1.00x | 2603.0 | — |
| 4 (`4/4/2`) | 10 | Rust sequential, 2 threads (`rust_seq_2t`) | 407.6791 | 101.9198 | 2.35x | 1256.5 | -1346.5 |
| 4 (`4/4/2`) | 10 | Rust pipelined, 2 threads (`rust_par_2t`) | 395.4648 | 98.8662 | 2.42x | 1256.0 | -1347.0 |
| 4 (`4/4/2`) | 10 | Rust pipelined, 1 thread (`rust_par_1t`) | 395.9561 | 98.9890 | 2.42x | 1274.0 | -1329.0 |

## Recommendation on Defaults (Coordinating with #1035)

1. **Keep pipelined parallel Rust rendering enabled by default** (`parallel_render=None` in `stream_animation_mp4`) for eligible views (`combined`, `biome`, `elevation` across `behrmann`, `eckert4`, `globe`), with `MANTLE_BLOOM_PARALLEL_RENDER=0` and `MANTLE_BLOOM_DISABLE_RUST=1` retained as clean opt-outs.
   - Pipelining adds negligible memory overhead over sequential Rust (`+34.7 MB` at `D=4, spf=1`, `-0.5 MB` at `D=4, spf=10`, `-86.1 MB` / `-30.8 MB` at `D=1`) while reducing wall time in every configuration and cutting Density-4 peak RSS by **~1.35–1.45 GB** relative to Python.
2. **Keep the bounded `rayon::ThreadPool` default** (`1` thread on `<= 4` logical CPUs, `2` threads on `> 4` CPUs via `MANTLE_BLOOM_RENDER_THREADS`):
   - `rust_par_1t` achieves virtually the same speedup (`2.16x–2.42x`) with even lower host CPU contention, making constrained (`<= 4`-core) containers/laptops safe while giving `2`-thread desktops extra overlap on single-frame and `spf=1` workloads.
3. **End-to-end step win note:** Sharing Rust `biomes.rs` (`grid_continentality`, `grid_dist_to_land_rad`, `biome_relative_shade_factor`, `smooth_biome_field*`) and `spill_path_rs` with `hydroclimate` yields a **2.05x–2.42x** full-animation speedup in addition to overlapping frame rendering with `step_world`.
