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

## Results

Replay results will be added here after both runs complete.
