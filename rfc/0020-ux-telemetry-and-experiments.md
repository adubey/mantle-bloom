# RFC 0020: Logging & Experimentation

- **Status:** Draft
- **Date:** 2026-11-16
- **Related issues / RFCs:** None yet

## Problem

We currently log only backend operational errors (checkpoint failures, worker startup timeouts, replay-log warnings, and animation stream failures). There is no structured record of how people use the product or how it performs for them:

1. **No product or UX visibility:** We can't see which entry path users take, how long until their first simulation step, which views and inspectors they open, what they click in lists, or whether hosted cold starts make them leave.
2. **No durable log pipeline:** Hosted storage is sized for world checkpoints, not unbounded logs. Any events written on the server must be rotated, exported off-box, and deleted locally without risking checkpoint storage.
3. **No off-box analysis:** Aggregation and dashboard generation should run on a separate machine via cron rather than consuming serving CPU or memory.
4. **No online experimentation:** We have no mechanism to gate behavior behind typed flags, run controlled experiments, or keep simulation replays and saved worlds reproducible when experimental behavior is active.

### Hard constraints

- **No paid logging or analytics services:** The owner will not pay for a third-party logging/analytics SaaS.
- **No meaningful increase in hosting storage or RAM cost:** The solution must not require a continuously running database container or large persistent volume.

## Proposed solutions

*(Note: Because this is a product/infrastructure RFC rather than a physical-simulation model change, the "what happens in nature" and scientific-literature sections of `rfc-authoring` do not apply.)*

This RFC is structured around **five architectural decisions**. Accepting the RFC records the chosen option for each decision.

### Decision 1: Where events are stored

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **1a. Append-only JSONL on volume → B2** *(Recommended)* | **$0** (within B2 free 10 GB storage & free uploads; <\$0.01/mo Railway volume at 50 MB cap) | Torn last line on crash; orphaned active file if worker dies; many small objects if uncompacted | Medium |
| **1b. SQLite table on volume → B2** *(Fallback)* | **$0** (same B2/volume footprint) | Write contention during concurrent requests; larger file overhead | Medium |
| **1c. Structured stdout → Railway logs** | **$0** (included in Railway plan) | Short retention window; no reliable bulk export for historical analysis | Low |
| **1d. Hosted free tier (e.g., PostHog)** | **$0** up to free-tier event cap, then paid | Hard monthly event cap; third-party client scripts; data leaves our control | Low |

**Why not a NoSQL store with cheap appends?** An append to a local file is already the cheapest append available, and a database adds cost or failure modes without solving any problem we have. A self-hosted store (MongoDB, Redis, ClickHouse) on Railway bills RAM and disk continuously, breaking the cost constraint. A hosted free tier (DynamoDB, Firestore, Atlas) puts a network call on the write path and imposes hard caps.

**Volume & scale estimate:** At ~300 bytes/event and ~200 events/session (~60 KB raw, ~6 KB gzipped per session), **1,000 sessions/day** produces ~60 MB/day raw (~6 MB/day gzipped, or **~2.2 GB/year**), well inside B2's 10 GB free storage tier. Above **~4,500 sessions/day** sustained, raw annual volume exceeds 10 GB and requires B2 lifecycle expiration on raw chunks after compaction.

**Recommendation:** **1a (JSONL files on volume → B2)**, with **1b (SQLite on volume → B2)** as the fallback if multi-process file handling proves fiddly.

### Decision 2: Where analysis runs

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **2a. Off-box cron + DuckDB → static HTML** *(Recommended)* | **$0** (runs on existing analysis machine; B2 downloads within 2,500/day free allowance) | Dashboard is periodic (e.g., daily/hourly), not real-time | Medium |
| **2b. In-place queries on Railway** | **$0** incremental, but risks CPU/RAM spikes | Competes with simulation workers for memory and CPU | Medium |

**Recommendation:** **2a (Off-box cron + DuckDB)** with daily compaction to Parquet so queries stay fast as object count grows.

### Decision 3: UI event model

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **3a. Element ontology + generic verbs** *(Recommended)* | **$0** | Requires upfront registry and lint enforcement | Medium |
| **3b. Ad-hoc per-feature event names** | **$0** | Time series breaks whenever UI labels or components change; easy to forget new controls | Low |

**Recommendation:** **3a (Element ontology + generic verbs)**. Stable IDs keep metrics intact across UI redesigns, and lint/CI enforcement prevents unlogged or silently broken controls.

### Decision 4: Experimentation scope

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **4a. Typed flags + plain A/B split** *(Recommended)* | **$0** | Cannot run overlapping multi-layer or multi-arm sweeps without extension | Low–Medium |
| **4b. Orthogonal layers + Latin Hypercube (LHS)** | **$0** | Severely underpowered at low traffic; requires multiple-comparison adjustment and per-experiment re-salting | High |
| **4c. None yet (logging only)** | **$0** | Cannot test UI or caching changes online | Low |

**Power calculation:** For a baseline 20% conversion metric, detecting a 5 percentage-point lift (20% → 25%) at 80% power and $\alpha = 0.05$ requires **~1,100 sessions per arm** (~2,200 total for a 2-arm A/B test). With $K = 10$ LHS arms, each arm receives $1/K$ of traffic (~11,000 sessions total). At early traffic (e.g., 100–300 sessions/week), a 2-arm test already takes ~2–5 months, and a 10-arm LHS would take years to reach significance.

**Recommendation:** **4a (Typed flags + plain A/B split)**, keeping **simulation physics out of online experiments** so saved worlds are never contaminated. Orthogonal layers and LHS remain a clean future extension once traffic justifies them (requiring per-experiment salt to avoid carryover and Benjamini–Hochberg / Bonferroni correction for multiple comparisons).

### Decision 5: Identity, privacy, and retention

- **Guest & session identity:** Mint an anonymous random `session_id` (stored in a dedicated cookie/session storage, separate from the existing view-state cookie) with a bounded session lifetime. Never log GitHub logins, IP addresses, free text, custom world names, or raw seeds.
- **Scope:** Telemetry is **off** in local and desktop builds; it runs only in hosted mode.
- **Retention:** Compacted daily Parquet files live on the analysis machine; raw B2 chunks can be expired after $N$ months (e.g., 6 months) via a free B2 lifecycle rule.
- **Notice / opt-out:** Owner to confirm whether a minimal privacy notice or opt-out toggle is desired in the hosted UI.

## Design

### 1. Event envelope & ontology (Decision 3a)

Every event is emitted as one JSON object with a shared envelope:

- `schema_version` (integer, starting at `1`)
- `event_id` (UUID, for deduplicating client retries)
- `seq` (per-session monotonic integer, so gaps measure event loss)
- `ts` (server UTC timestamp) and `client_ts` (browser timestamp for client events, since client clocks are unreliable)
- `build` (git SHA, tying events to the deployed ontology version)
- `mode` (`hosted`), `session_id`, `world_id` (when applicable)
- `event` (verb + target or backend milestone)

**UI Ontology:**
- Interactive UI elements are registered in a single typed registry, each with a stable `element_id` (named by function, e.g., `start_run`, never by label or layout position), `kind` (`button`, `toggle`, `list_item`, `slider`, `modal`, `map`), and parent surface.
- IDs are **append-only**: retired, never renamed or reused. A CI check compares the registry against `origin/main` and fails on removal or rename.
- Fixed verb set: `shown`, `clicked`, `changed`, `dismissed`. Label text and layout are logged as attributes, not encoded in the ID.
- **Correctness enforcement:** `track()` accepts only the TypeScript union of registered IDs; interactive controls use a shared wrapper that requires an ID, paired with a lint rule that fails on raw click handlers outside the wrapper.
- **List & position bias:** Any list surface (saved worlds, view/projection menus, premade worlds) emits a `shown` event with the presented items and attaches `list_id`, `index`, `length`, `ordering`, and `in_view` on `clicked`.

**Backend & ingest:**
- **Browser ingest endpoint:** `POST /telemetry/events` accepts batched client events and uses `navigator.sendBeacon` on `pagehide` so end-of-session events are not lost. It enforces a per-request size limit, per-session rate limit, basic bot/crawler filtering, and schema/registry validation (rejecting and counting unknown IDs or malformed payloads).
- **Backend events:** Lifecycle (`session_started`, `world_created`, `world_loaded`, `world_deleted`), compute (`step_completed`/`step_failed` with `client_disconnect` status, `run_started`/`run_stopped`, `animation_started`/`animation_stopped`/`animation_recovered`, `worker_lifecycle` for cold/idle wakes), rendering (`map_rendered`, `cross_section_rendered`), and errors (`client_error`).
- High-frequency gestures (e.g., globe drag-rotation frames) are not logged per frame.

### 2. Process lifecycle, rotation, B2 export & disk safety (Decision 1a)

- **Single writer & orphan recovery:** The hosted front process owns the telemetry directory and receives forwarded events from the simulation worker (or sweeps any `active-<pid>` file whose PID is dead at startup and periodically). Files are **not** gzipped during the 2-second graceful shutdown window; the startup/periodic sweep handles any unsealed file left by a crash or 90-second worker idle stop.
- **Rotation:** Active files seal and gzip (off the event loop in a background thread) every 5 minutes (shrinking the crash loss window), at 10 MB, or on worker idle. Readers skip any torn last line from a mid-write crash and increment a torn-line counter.
- **Two-phase B2 upload:** A background task uploads sealed `.jsonl.gz` files to `telemetry/v1/dt=YYYY-MM-DD/<filename>` using B2's `X-Bz-Content-Sha1` header and verifies the returned `contentSha1` before calling `unlink()`. A crash between upload and `unlink()` re-uploads the same deterministic filename (creating a harmless B2 version); the dashboard deduplicates by `event_id`.
- **Least-privilege credentials:** Railway uses a B2 application key restricted to the `telemetry/` prefix with write-only capability, so a web-service leak cannot read or delete world saves.
- **Disk safety & silent drops:** Local telemetry is capped at 50 MB (<\~\$0.01/month at ~\$0.15/GB-month) **and** checks free volume space so telemetry never starves checkpoints. Under disk pressure, new events are dropped (preserving pre-outage context). All telemetry failures are swallowed so user requests never fail; drop, eviction, and reject counters are emitted and uploaded so data gaps are distinguishable from quiet days.

### 3. Off-box cron & dashboard generation (Decision 2a)

A standalone script runs via cron on the analysis machine:
1. Tracks a high-water mark (last closed day synced) and lists only day-prefixes from there on, staying well within B2's 2,500/day free listing/download allowance.
2. Compacts each closed day into a single Parquet file locally and deletes the raw local chunks.
3. Queries the Parquet files with DuckDB to generate a static HTML dashboard covering activation funnels, time-to-first-step, view/inspector adoption, position-bias tables, `p50`/`p95` step/render latencies, worker cold-start tax, `seq` gaps, drop/reject counters, and experiment readouts.

### 4. Typed flags & plain A/B experimentation (Decision 4a)

- **Manifest & `Flags`:** A versioned JSON manifest defines active flags and plain A/B splits (hashing `session_id` or `world_id` with a per-experiment salt). A typed `Flags` object provides safe production defaults for local runs/tests plus env/CLI overrides (`MB_FLAGS` / `--flag`) to reproduce any arm locally.
- **Exposures:** `experiment_exposed` is emitted when code reaches the flagged decision point. Any in-memory deduplication is purely a volume optimization (since workers idle out after 90 s); the dashboard deduplicates on `(unit_id, experiment_id)` taking the first exposure as the source of truth.
- **Save compatibility & replay:** Under Decision 4a, online experiments do not alter simulation physics. If a future flag does affect world state, it must be stored on `World` (with a default value so `npz_save.thaw` loads older v5 saves cleanly, and `__setstate__` only for legacy v1–v4 pickle saves), recorded by `replay_log`, and loudly flagged/refused in debug and RFC evidence runs.

### Runtime cost & knobs

- **Runtime cost:** Event emission (JSON encoding, IPC/buffering, and threaded gzip) is assumed to add negligible latency to request handlers; this is unmeasured and gated by the falsification criterion on `/step` and `/render` `p95` latency.
- **Knobs (design choices):**
  - `MB_TELEMETRY_CHUNK_BYTES = 10_485_760` (10 MB)
  - `MB_TELEMETRY_CHUNK_AGE_S = 300` (5 min)
  - `MB_TELEMETRY_MAX_LOCAL_BYTES = 52_428_800` (50 MB)
  - `MB_TELEMETRY_MIN_FREE_BYTES` (volume free-space floor protecting checkpoints)

## Detailed design (provisional touchpoints)

- **New:** `backend/app/telemetry.py` — event schema, buffered JSONL writer, PID-aware rotation, disk/free-space guard, drop counters, and B2 upload/delete worker (sharing or extracting B2 client logic from `bin/b2_saves.py`).
- **New:** `backend/app/experiments.py` — manifest loader, per-experiment salted A/B bucketing, `Flags` dataclass, and exposure helper.
- **New:** `frontend/src/uiElements.ts` — append-only element registry, typed `track()` helper, and tracked UI wrapper.
- **New:** `bin/telemetry_dashboard.py` — off-box B2 day-prefix sync, daily Parquet compaction, and DuckDB static HTML generator.
- **Existing:** `backend/app/hosted.py` (`_IDLE_SECONDS = 90.0`), `backend/app/hosted_worker.py`, `backend/app/main.py` — `POST /telemetry/events` ingest endpoint, lifecycle/compute/render event emission, and `Flags` forwarding.
- **Existing:** `frontend/src/api.ts`, `frontend/src/App.tsx`, `frontend/src/cookies.ts` — guest `session_id`, batched/beacon client reporting, and wrapped controls.
- **Existing (if physics flags are ever enabled):** `backend/app/world.py`, `backend/app/npz_save.py` (`thaw`), `backend/app/persistence.py`, `backend/app/replay_log.py`.
- **Docs:** `docs/architecture.md`, `docs/api-reference.md`, `docs/save-compatibility.md`.

## Assumptions and uncertainty

- Railway volume free-space headroom is sufficient for a 50 MB telemetry buffer alongside checkpoints.
- Hosted traffic remains below ~4,500 sessions/day; higher traffic requires B2 lifecycle expiration on raw chunks after daily Parquet compaction.
- Under low early traffic, A/B tests are best reserved for high-frequency UI/latency metrics (e.g., list ordering, render caching) rather than rare conversion events.

## Acceptance criteria

- [ ] **Ontology & UI coverage:** CI fails if a registered `element_id` in `uiElements.ts` is removed or renamed vs. `origin/main`; lint fails on raw click handlers outside the tracked wrapper; a frontend test clicks every registered element and verifies the emitted event; list events include `list_id`, `index`, `length`, `ordering`, and `in_view`.
- [ ] **Ingest & abuse protection:** `POST /telemetry/events` rejects and counts oversized, malformed, and unknown-ID payloads, enforces per-session rate limits, and filters bots/crawlers.
- [ ] **Crash & loss accounting:** An end-to-end test kills the worker mid-session and verifies the orphaned active file is sealed and uploaded; torn last lines are skipped and counted; the dashboard reports `seq` gaps and drop/eviction/reject counters.
- [ ] **B2 & disk safety:** Unit tests verify `X-Bz-Content-Sha1` / `contentSha1` verification before `unlink()`, and new-event dropping when `MAX_LOCAL_BYTES` or free-space floor is hit without failing requests.
- [ ] **Cost ceiling:** Documented monthly cost stays \$0 on B2 (and <\$0.01/mo worst-case on Railway volume) up to ~1,000 sessions/day, with ~4,500 sessions/day identified as the free-tier threshold.
- [ ] **Flags & save compatibility:** Unit tests verify `Flags` defaults, `MB_FLAGS`/`--flag` overrides, and loading a real pre-change `.mbworld` save (via `npz_save.thaw`).
- [ ] **Dashboard:** `bin/telemetry_dashboard.py` compacts sample `.jsonl.gz` fixtures to Parquet and renders the static HTML report.
- [ ] **Repo gates:** `./bin/affected_test.sh --base origin/main`, `ruff check`, `./bin/quality_check.sh`, and frontend lint/build pass for each implementation PR.

## Falsification criteria

- If measured event loss (`seq` gaps + dropped/torn events) exceeds 5% under normal hosted operation, or if telemetry adds >2 ms `p95` latency to `/step` or `/render`, move serialization/compression fully out-of-process or sample high-volume events.
- If hosted traffic cannot power even a two-arm A/B test on primary UI metrics within 8 weeks, defer online experimentation (Option 4c) and keep only typed flags + logging.

## References

- [`bin/b2_saves.py`](../bin/b2_saves.py) — existing B2 `X-Bz-Content-Sha1` upload and credential conventions.
- [`docs/architecture.md`](../docs/architecture.md) — backend/frontend request flow and `HostedWorker` lifecycle.
- [`docs/save-compatibility.md`](../docs/save-compatibility.md) — `.mbworld` backward-compatibility requirements.
