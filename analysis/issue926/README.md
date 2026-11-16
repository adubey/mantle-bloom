# Issue 926: Dependency and toolchain refresh — replay A/B

## Pre-registration

- **Why:** `bin/pin_deps.sh` bumped `pillow` (`12.3.0` -> `13.0.0`) and `certifi` (`2026.7.22` -> `2026.10.23`), and also refreshed `numpy` (`2.5.2` -> `2.5.5`), `scipy` (`1.18.0` -> `1.18.1`), `numba` (`0.67.0` -> `0.68.1`), and `llvmlite` (`0.49.0` -> `0.50.0`). Per `AGENTS.md`, backend pin updates require a before/after replay to verify whether simulation results moved.
- **Save:** `mantle-bloom-seed47235302-0y.bloom` (`sha256: f946ca1b2e309cec2a29409321b073233a1832d466276d9ba1b1f7c3de0d4ed5`), replayed for 10.0 Myr (`--myr 10.0 --checkpoint-myr 2.5`) using `bin/debug/attribute_hm_ratchet.py`.
- **Arms:**
  - **Baseline (`origin/main`, `ed121ebc`):** detached worktree with original `backend/constraints.txt` (`numpy==2.5.2`, `scipy==1.18.0`, `numba==0.67.0`).
  - **Branch (`gemini/issue-926-dep-refresh`):** updated `backend/constraints.txt` (`numpy==2.5.5`, `scipy==1.18.1`, `numba==0.68.1`).
- **Expected direction and rough size:** Zero difference (or at most floating-point ULP noise) across all simulation metrics, since only patch versions of `numpy`/`scipy`/`numba` changed and no simulation code was altered.
- **Primary metrics:** Land area (fraction / Mkm²), mean continental crust thickness (`Hc`), mean mantle lithosphere thickness (`Hm`), `Hm` cap saturation fraction, and `unattributed` `Hm` change at 10.0 Myr.
- **Decision threshold:** Any non-zero difference > 0.1% in land area, `Hc`, or `Hm` will be investigated and documented; exact or near-exact match confirms the pin refresh does not move simulation results.
- **Controls:** Same save file, same `--myr 10.0 --checkpoint-myr 2.5` arguments, sequential execution.

## Command

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  <local-save> --myr 10.0 --checkpoint-myr 2.5 --out <scratch-output>
```

## Results (10.0 Myr, step 100)

Both arms ran sequentially while another simulation (`bug/compare_era5_transport_proxy.py`, pid 81239) was active on the machine, so runtimes are flagged as non-isolated.

| Metric | Baseline (`ed121ebc`) | Branch (`fe8a008e`) | Delta |
|---|---|---|---|
| Plate count | 21 | 21 | 0 |
| Sea level (m) | -83.723009 | -83.723009 | 0 |
| Land area (km²) | 160,009,493.96 | 160,009,493.96 | 0 |
| Continental area (km²) | 357,796,131.82 | 357,796,131.82 | 0 |
| Emergent continental area (km²) | 159,576,243.68 | 159,576,243.68 | 0 |
| Continental Hc p50 / p90 (km) | 30.913162 / 39.855937 | 30.913162 / 39.855937 | 0 / 0 |
| Continental Hm mean (km) | 98.085376 | 98.085376 | 0 |
| Hm cap saturation (M km²) | 0.00 | 0.00 | 0 |
| Hm volume all (km³) | 39,824,050,265.09 | 39,824,050,265.09 | 0 |
| Hm closure residual (km³) | -1.49e-06 | -1.49e-06 | 0 |
| Runtime (s, overlapped) | 862.5 | 842.8 | -19.7 |
| Peak RSS (MB) | 962.2 | 1018.1 | +55.9 |

## Interpretation

Every simulation metric, continental/Hm ledger, and closure residual at 10.0 Myr is **bit-for-bit identical** (`exact_match_excluding_timing_rss: true` in [comparison.json](comparison.json)) between `origin/main` and the refreshed pins (`numpy==2.5.5`, `scipy==1.18.1`, `numba==0.68.1`, `llvmlite==0.50.0`, `pillow==13.0.0`, `certifi==2026.10.23`). The dependency refresh does not move simulation results.

