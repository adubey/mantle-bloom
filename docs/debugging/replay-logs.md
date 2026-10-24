# Replay logs: catching a rare step as it happens

A rare event (a merge that clips craton, a ledger that stops closing) is hard to reproduce by
replaying a seed for hundreds of Myr. Started with `--replay-log`, the server saves the world
before any step whose health metrics go askew, so the step can be rerun offline (#442; how it
works: [architecture.md#replay-log](../architecture.md#replay-log)).

```bash
./bin/restart.sh --port 8442 --replay-log   # records under ~/.mantle-bloom/replay-logs/8442/
```

Set `MANTLE_BLOOM_REPLAY_LOG_DIR` to put them elsewhere. A record is a directory named
`<UTC time>-seed<seed>-<Myr before the step>myr-<triggers>` holding `before.mbworld` and
`record.json`. Find them with `GET /worlds/{world_id}/replay-logs`, the `?deb` Event Console
(`Replay log: <triggers> -> <name>`), or the `replay log:` warnings in
`/tmp/mantle-bloom-<port>.log`.

The starting triggers, each a constant in `replay_log.py`:

| Trigger | Fires when one step... |
|---|---|
| `craton_step_loss` | books more than `CRATON_STEP_LOSS_FRACTION` (1%) of the starting live craton to one `CRATON_SINK_ACCOUNTS` account |
| `craton_ledger_balance` | moves `cratons.balance_error_m3` by more than `CRATON_BALANCE_TOLERANCE` (1e-9) of live craton |
| `continental_ledger_open` | leaves `continental_ledger.assert_closed` failing when it passed before |
| `continental_area_jump` | changes genuine continental area by more than `CONTINENTAL_AREA_JUMP_FRACTION` (5%) |
| `land_fraction_jump` | changes the area-weighted land fraction by more than `LAND_FRACTION_JUMP` (0.05) |
| `plate_count_jump` | changes the plate count by more than `PLATE_COUNT_JUMP_FRACTION` (25%) |
| `hc_at_max` | takes the share of nodes with Hc at `MAX_CRUSTAL_THICKNESS_M` above `AT_CAP_FRACTION` (0.05%) |
| `hm_at_max` | takes the share of nodes with Hm at `MAX_MANTLE_LITHOSPHERE_THICKNESS_M` above 0.05% |
| `glacier_at_max` | takes the share of nodes with glacier ice at `GLACIER_MAX_DEPTH_M` (5000 m) above 0.05% |
| `step_exception` | raises; the record has no post-step state |

The three cap triggers count nodes the way the Stats panel's `hc_at_max_fraction` and
`hm_at_max_fraction` do, and fire on the step that crosses the threshold, not on every step
while the share stays above it. The thresholds are loose to start with: in 20 one-Myr steps
of a default seed-42 world none fired.
Each trigger is rate-limited to one record per 20 logged steps, and records past 2 GB in total
are removed oldest first.

`record.json` has the triggers and their values, `step_years`, `seed`,
`elapsed_years_before`, `steps_taken_before`, the world's scalar `controls` (flags, densities,
`/worlds/{world_id}/controls` settings), `git` (commit and dirty flag), `step_seconds`, `health_before`
and `health_after`, `fingerprint_after`, and any `annotations`.

Replay a record, checking it reproduces:

```bash
backend/.venv/bin/python bin/debug/replay_logged_event.py \
    ~/.mantle-bloom/replay-logs/8442/<record> [--save-after /tmp/after.mbworld]
```

It loads `before.mbworld` (which carries the controls), runs the recorded step, and compares
the fingerprint, plate count and craton ledger with the record. A `step_exception` record
reproduces when the replay raises the same exception type. It prints `REPRODUCED` or lists
the world attributes that differ, and exits 1 on a mismatch. A mismatch usually means the code
changed since the record was written: check out the record's `git` commit. It ends with a
`probe_craton_clips.py` command for the same step, since `before.mbworld` is an ordinary save
any `bin/debug` script loads.

To catch a new signature, register it in `replay_log.py` (or from an RFC's module):

```python
def _my_signature(outcome: replay_log.StepOutcome) -> dict | None:
    ...  # outcome.before / outcome.after are replay_log.Health; return values or None

replay_log.TRIGGERS["my_signature"] = _my_signature
```

Tests: [`unit_tests/test_replay_log.py`](../../backend/unit_tests/test_replay_log.py).
