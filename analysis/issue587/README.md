# Issue #587 replay comparison

## Preregistration

**Registered before replay.** Replay seed 997271774 from the 87.4 Myr save for 43.1 Myr (431 steps at 100 kyr), comparing `origin/main` at `f04cbc0772eca6d2484f7477b20004deb368206c` with the RFC 0014 D1a branch at `244ea48f4a9ece9f358d2bde691d1623fed082a8`.

Save SHA-256: `08ad9461d9585e92be93a58d81e326c7e4b685ae7d462135d01b0e9931ed762b`.

**Expected direction and scale:** D1a changes lake-margin rain from dry Budyko runoff to direct lake supply and removes transit loss from that wet share, while glacier melt keeps its previous transit loss. The net change in land and continental area has no reliable directional prediction; expect it to remain within 2% after 43.1 Myr.

**Primary measure:** continental area at the endpoint. **Decision threshold:** a difference greater than 2% exceeds the comparable run-to-run continental-area variation reported for a 43 Myr replay in #315. Land area and crust/mantle ledgers are secondary context; report runtime and peak RSS, but do not interpret a single-run timing difference as a performance result.

**Controls:** same save bytes, seed, 431-step horizon, timestep, and replay script in both arms. Baseline is a task-specific detached worktree at the recorded `origin/main` commit. Runs execute one at a time with no tests or other replay/probe running. One run per arm; if the primary difference is near 2%, treat the result as inconclusive rather than expanding into extra runs.

## Reproduction

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py \
  <local-save> --myr 43.1 --checkpoint-myr 2.5 --out <scratch-output>
```

Both successful runs used a scratch `NUMBA_CACHE_DIR`. The first baseline launch failed because Numba could not locate a writable cache directory in its managed worktree; it was relaunched with that cache setting. Only the successful run is included below.

## Results

Both runs completed all 431 steps to 130.5 Myr. Continental area changed by **+0.0164%**, within the preregistered 2% threshold.

| Metric at 130.5 Myr | origin/main | RFC 0014 D1a | Change |
|---|---:|---:|---:|
| Continental area (km²) | 264,942,578 | 264,986,124 | +0.0164% |
| Land area (km²) | 152,129,138 | 151,873,845 | −0.1678% |
| Emergent continental area (km²) | 151,208,071 | 150,955,595 | −0.1670% |
| Hm-cap area (km²) | 15,331,121 | 15,155,673 | −1.144% |
| Hm ledger residual (km³) | −142,315.54 | 0.0000048 | — |
| Runtime (s) | 3,272 | 3,529 | +7.85% |
| Peak RSS (MiB) | 1,741 | 1,611 | — |

The Hm ledger residual is included as secondary diagnostic context; the primary preregistered area metric is within threshold. Both runs intermittently warned that one node reached `GLACIER_MAX_DEPTH_M` with zero meltwater shed. The single-run runtime difference is not a performance result.

The compact metric record is in [`replay-preregistration.json`](replay-preregistration.json).
