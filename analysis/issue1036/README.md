# Issue #1036: Accelerate Glacier March in Rust

## 1. Hotspot Confirmation

Profiling young-world steps on `origin/main` (`bf382b23`, `profile_phases.py --seed 0 --scenarios step --steps 3 --warmup 1 --wind-model diagnostic`) confirmed that the glacier march dominates stepping time at both densities (consistent with mature-world figures of 72–90%+ in `analysis/issue147/README.md`):

- **Density 4 (`4/4/2`, 124,046 nodes, 17 plates):** Glacier march takes **20.96 s / 26.47 s** per step (**79.2%**):
  - `Solving march ice`: `14.26 s`
  - `Solving connected-ocean sea level`: `5.09 s`
  - `Re-solving march sea level` (own time): `1.61 s`
  - Non-march step time: `5.51 s`
- **Density 1 (`1/1/1`, 30,712 nodes, 19 plates):** Glacier march takes **8.77 s / 10.91 s** per step (**80.4%**):
  - `Solving march ice`: `7.18 s`
  - `Solving connected-ocean sea level`: `1.17 s`
  - `Re-solving march sea level` (own time): `0.42 s`

## 2. Implementation Summary

- **Rust PyO3 extension (`backend/rust/`, module `mantle_bloom_rust`):**
  - `MarchGraph`: prebuilt compact `u32`/`i32` CSR and flood adjacency reused across all substeps of a march.
  - `step_ice` (`ice_step.rs`, `distance.rs`, `plastic.rs`, `climate.rs`): parallel setup with `rayon`, parallel signed-distance sweeps (`4`-ary `IndexedHeap`), local plastic domain solver with in-place scratch reuse and scalar climate evaluation, frozen-precip / ablation / grounding-line mass balance, releasing the GIL via `py.detach`.
  - `seeded_flood_levels`, `solve_sea_level_connected`, `below_sea_components`, `connected_ocean_mask` (`eustasy.rs`, `graph.rs`): priority-flood with Barnes plain-queue optimization, warm-bracket ceiling pruning, pre-integrated deep-ocean + coastal-slice connected-ocean sea-level solve, and `u32` union-find component labeling.
  - `march_substep_prep`, `apply_march_abrasion` (`climate.rs`): substep climate prep and abrasion kernels.
- **Python wiring & fallback (`backend/app/hydroclimate/_glacier_rust.py`, `ice_sheet.py`, `eustasy.py`, `hydrology.py`):**
  - Uses Rust automatically when compiled and no test monkeypatch overrides `_shape` or `_seeded_flood_levels` (with `MANTLE_BLOOM_DISABLE_RUST=1` override and transparent pure-Python fallback).
  - Kept the outer `_march_glacier` Python orchestration (`0.13 s` own time; PyO3 marshaling ~`0.002 s`) and `backend/constraints.txt` untouched.
- **Local, CI, and Railway builds:**
  - `bin/build_rust.py` builds via `cargo build --release --locked` and installs `mantle_bloom_rust.so` into the active venv without altering pinned Python packages (`bin/setup.sh` invokes `--skip-if-unchanged` unless `MANTLE_BLOOM_SKIP_RUST=1`).
  - `Dockerfile` compiles `mantle_bloom_rust.so` in a `rust-builder` stage (`cargo build --release --locked`) and copies it into the runtime image `platlib`, verified with `import mantle_bloom_rust, app.hosted` in CI (`docker` job).
  - Desktop packaging (`bin/package.sh`, `bin/mantle-bloom.spec`, `.github/workflows/package.yml`) includes `mantle_bloom_rust`.

## 3. Correctness & Equivalence Verification

- **Unit tests:** `backend/unit_tests/test_glacier_rust.py` (5 tests) plus `test_ice_sheet.py`, `test_eustasy.py`, `test_hydrology.py` (245 tests passed).
- **Same-input parity (`seed=0, d=1`, substep 97):** `exact_step_ab.py` measured a max ice depth difference of `3.99e-10 m` against the Python reference on identical substep inputs.
- **Multi-step sensitivity (`seed=0, d=1`, 3 steps):** `generate_world` matched to `4.5e-13 m` on `glacier_depth`, and substeps `0..96` matched to `<1e-4 m`. At substep 97, `1e-13` floating-point rounding flipped a discrete `trend-1k` `SettleWatch` convergence threshold (`4.5e-4 m`), branching substep counts (`245` vs `198` on step 0, yielding `-180.98 m` vs `-188.00 m` sea level on step 0 and `-122.86 m` vs `-116.71 m` after 3 steps). At Density 4 (`4/4/2`), substep count (`235.0` mean), lake tier bodies (`20, 23, 23, 25`), and gross overflow (`2.04e+08, 9.83e+08, 3.13e+08, 6.41e+08 km^3`) remained bit-for-bit identical across all 4 steps.

## 4. Before/After Benchmarks (Idle Machine, Back-to-Back)

Measured with `profile_phases.py --seed 0 --scenarios step --steps 3 --warmup 1 --wind-model diagnostic` and `/usr/bin/time -l`:

### Density 4 (`4/4/2`, 124,046 nodes, 17 plates)

| Metric | `origin/main` (`bf382b23`) | Rust (`issue-1036`) | Speedup |
|---|---:|---:|---:|
| Step wall time (mean of 3) | **26.47 s** (`26.78, 26.15, 26.47`) | **12.68 s** (`13.06, 12.38, 12.59`) | **2.09x** |
| Glacier march total | **20.96 s** | **7.18 s** | **2.92x** |
| — `Solving march ice` | 14.26 s | 5.97 s | **2.39x** |
| — `Solving connected-ocean sea level` | 5.09 s | 0.34 s | **14.97x** |
| — `Re-solving march sea level` (own) | 1.61 s | 0.87 s | **1.85x** |
| Non-march step time | 5.51 s | 5.50 s | 1.00x |
| Mean substeps / step | 235.0 | 235.0 | — |
| Peak RSS | 974.8 MB (`974,798,848` B) | 1,012.9 MB (`1,012,858,880` B) | +3.9% |

### Density 1 (`1/1/1`, 30,712 nodes, 19 plates)

| Metric | `origin/main` (`bf382b23`) | Rust (`issue-1036`) | Speedup |
|---|---:|---:|---:|
| Step wall time (mean of 3) | **10.91 s** (`11.27, 11.22, 10.24`) | **4.52 s** (`4.44, 4.46, 4.66`) | **2.41x** |
| Glacier march total | **8.77 s** | **2.40 s** | **3.65x** |
| — `Solving march ice` | 7.18 s | 2.12 s | **3.39x** |
| — `Solving connected-ocean sea level` | 1.17 s | 0.08 s | **14.6x** |
| — `Re-solving march sea level` (own) | 0.42 s | 0.20 s | **2.10x** |
| Mean substeps / step | 236.7 | 200.3 | — |
| Peak RSS | 441.3 MB (`441,286,656` B) | 504.4 MB (`504,446,976` B) | +14.3% |
