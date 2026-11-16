# RFC 0020: UX Telemetry, B2 Log Export, and Online Experimentation

- **Status:** Draft
- **Date:** 2026-11-16
- **Related issues / RFCs:** None yet

## Problem

`mantle-bloom` currently logs only backend operational errors (checkpoint failures, worker startup timeouts, replay-log warnings, and animation stream errors in `hosted.py`, `hosted_worker.py`, `run_control.py`, and `main.py`). There is no structured record of how users actually experience the product:

1. **No product/UX visibility:** We cannot measure onboarding funnels (`demo` vs. `new_world` vs. `import`), time-to-first-step, which map views/projections or inspectors (`PlateInspector`, `RiverInspector`, `LakeInspector`, `PlatesAndFaults`, `MantleCrossSection`) are used, or how often hosted cold starts (`_IDLE_SECONDS = 90.0`) and step latencies cause abandonment.
2. **No durable log pipeline on Railway:** The Railway volume (`RAILWAY_VOLUME_MOUNT_PATH` / `MANTLE_BLOOM_CHECKPOINT_DIR`) is sized for world checkpoints, not unbounded log growth. Any telemetry written on Railway must be rotated, exported to Backblaze B2, and deleted locally without risking disk exhaustion.
3. **No off-box dashboarding:** Aggregation and visualization should run on a separate machine via cron rather than consuming Railway CPU/memory.
4. **No online experimentation or flag system:** We have no mechanism to split traffic across experimental arms, run overlapping experiments or Latin Hypercube (LHS) parameter sweeps online, or pass behavioral flags to the frontend and Python simulation/render binaries while preserving `.mbworld` determinism.

## Proposed solutions

### Alternatives considered

1. **Third-party SaaS analytics + feature-flag service (e.g., PostHog, LaunchDarkly):**
   - *Pros:* Off-the-shelf dashboards and flag UIs.
   - *Cons:* Adds external runtime dependencies, recurring cost, and data egress; awkward to wire into standalone Python simulation/CLI binaries and offline `.mbworld` replay determinism.
2. **Direct database writes on Railway (SQLite/Postgres):**
   - *Pros:* Queryable in place.
   - *Cons:* Adds write contention during simulation steps, consumes Railway volume/DB storage, and couples serving performance to analytics queries.

### Recommended solution

Build a self-contained, four-part pipeline aligned with `mantle-bloom`'s existing B2 tooling (`bin/b2_saves.py`):

1. **Structured raw event log** (`telemetry.py`): Append-only, buffered JSONL events emitted at action/request boundaries (kept strictly separate from downstream computed metrics).
2. **Two-phase B2 export & delete**: Rotate sealed `.jsonl.gz` chunks on Railway, upload to B2 with SHA-256 verification, `unlink()` only after confirmation, and enforce a local disk cap if B2 is unreachable.
3. **Off-box cron & static dashboard**: A scheduled job on an analysis machine incrementally syncs `telemetry/v1/` from B2 and runs DuckDB/Python aggregations to produce a static HTML dashboard and experiment readouts.
4. **Layered online experimentation & binary flags** (`experiments.py`): Orthogonal traffic layers (`sha256(f"{unit_id}:{layer_id}") % 10_000`), native support for Latin Hypercube design points as layer arms, deduplicated exposure logging, and a typed `Flags` object with CLI/env overrides and `World` pinning for simulation determinism.

## Design

### 1. Raw event stream vs. derived dashboard metrics

Raw logs are immutable, atomic records; funnels, percentiles, and experiment readouts are computed downstream on the dashboard machine.

Every event shares a compact envelope:

```json
{
  "ts": "2026-11-16T07:15:00.123Z",
  "event": "step_completed",
  "mode": "hosted",
  "session_id": "s_...",
  "world_id": "w_...",
  "exp_arms": {"render_cache": "treatment_v2"}
}
```

Core event catalog:
- **Lifecycle & creation:** `session_started` (`entry_path`, `user_type`), `world_created` (`method`: `generate|sketch|import|debug_scenario`, slider values—not raw names/seeds), `world_loaded`, `world_deleted`.
- **Simulation & compute:** `step_completed` / `step_failed` (`years_stepped`, `elapsed_years`, `duration_ms`, `simulate_climate_biomes`, `status`: `ok|error|client_disconnect`), `run_started` / `run_stopped`, `animation_started` / `animation_stopped`, `animation_recovered`, `worker_lifecycle` (`cold_start|wake_from_idle|crash_recovery`, `wait_ms`).
- **Rendering & exploration:** `map_rendered` (`view`, `projection`, `width`, `height`, `duration_ms`, `cache_hit`), `cross_section_rendered`, `inspector_opened`, `point_probed` (`endpoint`, `hit` boolean), `controls_changed` (`changed_keys`).
- **Exports & errors:** `world_exported` (`format`: `mbworld|hexgrid|animation`), `public_world_viewed`, `client_error` (`source`, `status_code`, `endpoint`).
- **Experiments:** `experiment_exposed` (`experiment_id`, `layer_id`, `arm_id`, `unit_id`, `params`), deduplicated in memory per `(unit_id, experiment_id)`.

High-frequency gestures (e.g., globe drag-rotation frames) are **not** logged per frame; if needed, they are summarized on idle/session end.

### 2. Log rotation, B2 export, and local deletion

- **Active vs. sealed files:** Processes write to `telemetry/active-<pid>-<ts>.jsonl`. A file is rotated and gzip-compressed to `telemetry/sealed-<instance>-<ts>.jsonl.gz` when it exceeds `MAX_CHUNK_BYTES` (default 10 MB), `MAX_CHUNK_AGE_S` (default 900 s), or on graceful shutdown.
- **Two-phase export & purge:** A background task in `hosted.py` scans sealed `.jsonl.gz` files, uploads each to `telemetry/v1/dt=YYYY-MM-DD/<filename>` via B2 (`b2_upload_file` with `X-Bz-Info-sha256`), and `unlink()`s the local file only after B2 returns HTTP 200 with matching checksum.
- **Disk safety valve:** If B2 is unreachable, total local `telemetry/` size is capped at `MAX_LOCAL_BYTES` (default 200 MB). When exceeded, the oldest sealed files are deleted first so telemetry never starves `CheckpointStore`.

### 3. Off-box cron & dashboard generation

A standalone stdlib-first script (`bin/telemetry_dashboard.py`, runnable via cron on another machine using `~/.config/mantle-bloom/b2.env`):
1. **Incremental pull:** Lists `telemetry/v1/` in B2 and downloads only keys not yet present in `~/.cache/mantle-bloom-telemetry/`.
2. **Aggregate:** Queries the local `.jsonl.gz` files (via DuckDB or Python) to compute:
   - Activation funnel (`session_started` → `world_created` → first `step_completed`) and time-to-first-step (`p50`/`p95`).
   - Feature adoption (`view × projection` matrix, inspector attach rates, probe miss rate, export conversion).
   - SLIs (`step_completed.duration_ms` by `elapsed_years`, `map_rendered` latency & cache hit rate, `worker_lifecycle` cold-start tax, `client_error` hotspots).
   - **Experiment readouts:** Per `(experiment_id, arm_id)` sample counts, primary/guardrail metric means and confidence intervals, and parameter-response tables for LHS layers.
3. **Output:** Writes a self-contained static HTML dashboard (plus JSON summary) to a configured output directory.

### 4. Online experimentation: overlapping layers & Latin Hypercube

- **Orthogonal layers:** Experiments are grouped into named **layers** (e.g., `ui_onboarding`, `render_cache`, `lhs_hydrology`). Within a layer, traffic is partitioned across `[0, 10_000)` buckets so arms do not overlap. Across layers, experiments overlap independently via:
  ```python
  bucket = int.from_bytes(hashlib.sha256(f"{unit_id}:{layer_id}".encode()).digest()[:4], "big") % 10_000
  ```
- **Latin Hypercube (LHS) support:** An LHS experiment is represented as a single layer containing `K` pre-generated design-point arms (`arm_0` .. `arm_{K-1}`), each carrying a `params` dict of parameter values (plus an optional `control` arm). Uniform bucket assignment distributes traffic evenly across all `K` strata.
- **Randomization unit:** Each experiment declares `unit`: `"session"` (UI/rendering/latency) or `"world"` / `"user"` (simulation behavior or multi-session retention).
- **Exposures:** `experiment_exposed` is emitted only when code actually reaches the flagged decision point, preventing unexposed traffic from diluting statistical power.

### 5. Binary flags & simulation determinism

- **Manifest (`experiments.json`):** A versioned JSON config (loaded from env path, local file, or B2) defining active layers, arms, traffic ranges, and default flag values.
- **Typed `Flags` object:**
  - `hosted.py` / `main.py` evaluate `Flags` for the current `(session_id, world_id)`, return active UI flags in session/world bootstrap responses to the frontend, and forward backend flags to `HostedWorker`.
  - Simulation and CLI binaries accept `Flags` with safe production defaults (so unit tests and local runs need zero config) plus `MB_FLAGS` env / `--flag key=value` overrides to reproduce any experimental arm locally.
- **Save compatibility & determinism:** Any flag that alters world generation or step physics must be randomized on `world_id` and recorded on `World` (e.g., `world.experiment_flags: dict[str, Any]`, defaulting to `{}` on old `.mbworld` saves via `getattr` / `__setstate__`) so reloading a checkpoint preserves the world's physics.

### Runtime cost & knobs

- **Runtime cost:** Buffered JSONL append is < 0.1 ms per event; gzip + B2 upload runs asynchronously on sealed chunks outside request paths.
- **Knobs (design choices, tunable via env):**
  - `MB_TELEMETRY_CHUNK_BYTES = 10_485_760` (10 MB)
  - `MB_TELEMETRY_CHUNK_AGE_S = 900` (15 min)
  - `MB_TELEMETRY_MAX_LOCAL_BYTES = 209_715_200` (200 MB)
  - `MB_EXPERIMENTS_BUCKETS = 10_000`

## Detailed design (provisional touchpoints)

- **New:** `backend/app/telemetry.py` — event schema, buffered JSONL writer, rotation, disk-cap enforcement, and B2 upload/delete worker.
- **New:** `backend/app/experiments.py` — manifest loader, orthogonal layer bucketing, LHS arm parameter resolution, `Flags` dataclass, and deduplicated `experiment_exposed` helper.
- **New:** `bin/telemetry_dashboard.py` — off-box B2 incremental sync, DuckDB/Python aggregations, experiment/LHS readout builder, and static HTML generator.
- **Existing:** `backend/app/hosted.py` — start telemetry exporter task, log `worker_lifecycle` cold/idle wakes, evaluate/forward `Flags`.
- **Existing:** `backend/app/hosted_worker.py` & `backend/app/main.py` — emit lifecycle/step/render/probe/error events, expose UI flags, accept `Flags` context.
- **Existing:** `backend/app/world.py` — optional `experiment_flags` field with backward-compatible loading for old `.mbworld` saves.
- **Existing:** `frontend/src/api.ts` & `frontend/src/App.tsx` — consume active UI flags and report `inspector_opened` / `client_error` events.
- **Docs:** `docs/architecture.md`, `docs/api-reference.md`, `docs/save-compatibility.md`.

Implementation is intended to be split into four dependency-ordered sub-issues:
1. Core event logger + Railway rotation + B2 export/delete (`telemetry.py`, `hosted.py`).
2. Experiment manifest, orthogonal/LHS bucketing, and `Flags` wiring (`experiments.py`, `main.py`, `world.py`).
3. Frontend UX event & flag integration (`api.ts`, `App.tsx`).
4. Off-box cron sync & static dashboard/experiment readout generator (`bin/telemetry_dashboard.py`).

## Assumptions and uncertainty

- Railway provides `RAILWAY_VOLUME_MOUNT_PATH` (or `MANTLE_BLOOM_CHECKPOINT_DIR`) with sufficient headroom for a 200 MB telemetry buffer alongside checkpoints.
- Traffic volume on hosted mode is low enough that 15-minute / 10 MB chunks keep B2 Class B/C API calls negligible.
- For LHS experiments with many strata (`K`), statistical power per stratum is limited by active user count; the dashboard should report per-parameter regression/response surfaces in addition to per-arm means.

## Acceptance criteria

- [ ] Unit tests verify JSONL rotation, SHA-256-verified B2 upload + local `unlink()`, and disk-cap eviction when B2 upload fails.
- [ ] Unit tests verify orthogonal layer hashing (independence across layers, determinism for a given `(unit_id, layer_id)`), LHS arm parameter resolution, and `experiment_exposed` deduplication.
- [ ] Unit tests verify `Flags` defaults, `MB_FLAGS` / `--flag` overrides, and old `.mbworld` save compatibility when `experiment_flags` is absent.
- [ ] Off-box `bin/telemetry_dashboard.py` generates a valid static HTML report from sample `.jsonl.gz` fixtures including funnel, latency, and experiment arm tables.
- [ ] `./bin/affected_test.sh --base origin/main`, `ruff check`, `./bin/quality_check.sh`, and frontend lint/build pass for each implementation PR.

## Falsification criteria

- If B2 upload/delete lags behind event generation under normal hosted load and repeatedly hits `MAX_LOCAL_BYTES`, the chunking/upload cadence or event verbosity (e.g., high-frequency probe/render events) must be reduced or sampled.
- If `telemetry.py` adds measurable latency (> 2 ms p95) to `/step` or `/render` requests, event emission must move to an async queue off the request path.

## References

- [`bin/b2_saves.py`](../bin/b2_saves.py) — existing B2 sync-and-purge pattern and credential conventions.
- [`docs/architecture.md`](../docs/architecture.md) — backend/frontend request flow and `HostedWorker` lifecycle.
- [`docs/save-compatibility.md`](../docs/save-compatibility.md) — `.mbworld` backward-compatibility requirements.
