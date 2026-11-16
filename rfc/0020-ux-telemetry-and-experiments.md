# RFC 0020: Logging & Experimentation

- **Status:** Draft
- **Date:** 2026-11-16
- **Related issues / RFCs:** None yet

## Problem

We currently log only backend operational errors (checkpoint failures, worker startup timeouts, replay-log warnings, and animation stream failures). There is no structured record of how people use the product or how it performs for them:

1. **No product or UX visibility:** We can't see which entry path users take, how long until their first simulation step, which views and inspectors they open, what they click in lists, or whether hosted cold starts make them leave.
2. **No link between UI intent and simulation outcomes:** Even if UI clicks are recorded, we have no authoritative record of which simulation choice the backend accepted or what run/step outcome followed.
3. **No durable log pipeline or off-box analysis:** Hosted storage is sized for world checkpoints, not unbounded logs, and aggregation should run off-box rather than consuming serving CPU or memory.
4. **No online experimentation:** We have no mechanism to gate behavior behind typed flags, run controlled experiments, or keep simulation replays and saved worlds reproducible when experimental behavior is active.

### Constraints & cost objective

- **Minimize upfront engineering and operational burden**, while keeping ongoing usage costs bounded by an explicit budget ceiling (target **\$0–\$50/month** at early traffic, revisited if adoption/monetization grows).
- **No meaningful increase in hosting RAM or persistent disk cost:** Avoid continuously running heavy database containers on the app host, and respect Railway's volume constraint (volumes cannot be attached to replicated services, and Railway does not provide sticky sessions).

## Proposed solutions

*(Note: Because this is a product/infrastructure RFC rather than a physical-simulation model change, the "what happens in nature" and scientific-literature sections of `rfc-authoring` do not apply.)*

This RFC is structured around **five architectural decisions**. Accepting the RFC records the chosen option for each decision.

### Decision 1: Ingestion, storage, and stack

| Option | Cost (at 1,000 sessions/day ≈ 6M events/mo) | Operational / scaling trade-offs | Upfront effort |
|---|---|---|---|
| **1a. Append-only JSONL on volume → B2** *(Recommended for single-replica)* | **\$0/mo** on B2 (within 10 GB free storage & free uploads); **<\$0.01/mo** Railway volume (50 MB cap) | **Single-replica only:** Railway volumes cannot be used with >1 replica. Requires rotation, orphan sweeps, and B2 upload/delete logic. | Medium–High |
| **1b. Managed pay-as-you-go (e.g., PostHog)** | **~\$187/mo** at 6M anonymous events/mo (1M/mo free tier, then usage pricing); **\$0/mo** below 1M events/mo | Minimal upfront engineering and native horizontal scaling, but third-party scripts, data leaves our control, and usage bills scale linearly with event volume. | Low |
| **1c. Small self-hosted analytics DB (e.g., GoatCounter/Umami/Plausible on Postgres or SQLite)** | **\$0** license + Railway DB RAM/disk cost (or SQLite on single-replica volume) | Common pattern for small self-hosted analytics (Postgres for multi-replica, SQLite for single-node), but adds DB operations/migrations and continuous Railway RAM/disk billing for Postgres/ClickHouse. | Medium |
| **1d. Structured stdout → Railway logs** | **\$0/mo** (included in Railway plan) | Short retention window; no reliable bulk export for historical analysis. | Low |

**Why not a NoSQL store with cheap appends?** A self-hosted store (MongoDB, Redis, ClickHouse) bills RAM and disk continuously on Railway, while a hosted free tier (DynamoDB, Firestore, Atlas) puts a network call on the write path and imposes hard caps. Local file append is the cheapest append for a single replica.

**Scale boundary & multi-replica migration path:** Railway's documentation states that volumes cannot be used with multiple replicas and scaling does not support sticky sessions. Therefore, **1a (and any SQLite-on-volume fallback) is strictly a single-replica architecture**. Before scaling the hosted service beyond 1 replica, we must migrate ingestion off the local volume—either by buffering batches in memory per replica and flushing unique `<replica_id>` chunks directly to B2, or by switching to a shared central sink (Option 1b or 1c with Postgres).

**Volume & retention arithmetic (Option 1a):**
- At ~300 bytes/event and ~200 events/session (~60 KB/session raw, ~6 KB/session gzipped), **1,000 sessions/day** produces ~60 MB/day raw and **~6 MB/day gzipped** (~180 MB/month, **~2.2 GB/year compressed**).
- Under **1-year compressed retention**, B2's 10 GB free storage tier supports up to **~4,500 sessions/day** (~10 GB/year compressed).
- Under **6-month compressed retention** (via a B2 lifecycle rule that also expires old hidden object versions left by upload retries), B2's 10 GB free tier supports up to **~9,000 sessions/day** (~10 GB per 6 months compressed).

**Recommendation:** **1a (JSONL on volume → B2)** while single-replica and under 1M–6M events/month, with **1b (managed pay-as-you-go, capped by budget)** as the low-effort alternative if engineering time is prioritized over usage cost, and direct-to-B2 / central sink as the required path before enabling >1 replica.

### Decision 2: Where analysis runs

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **2a. Off-box cron + DuckDB → static HTML** *(Recommended)* | **\$0** (runs on existing analysis machine; B2 downloads within 2,500/day free allowance) | Periodic (hourly/daily), not real-time; DuckDB's single-writer model fits offline batch jobs, not replicated live ingest. | Medium |
| **2b. In-place queries on Railway** | **\$0** incremental, but risks CPU/RAM spikes | Competes with simulation workers for memory and CPU. | Medium |

**Recommendation:** **2a (Off-box cron + DuckDB)** with daily compaction to Parquet so queries stay fast as object count grows.

### Decision 3: UI ontology & backend domain events

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **3a. Element ontology + generic verbs + backend domain events** *(Recommended)* | **\$0** | Requires upfront registry, lint enforcement, and additive lifecycle discipline | Medium |
| **3b. Ad-hoc per-feature event names** | **\$0** | Time series breaks on UI redesigns; no link between UI intent and accepted backend simulation choices | Low |

**Recommendation:** **3a**. Stable element IDs keep UI metrics intact across redesigns, while separate backend domain events capture accepted simulation choices and their resulting outcomes.

### Decision 4: Experimentation scope

| Option | Monthly cost | Failure modes | Effort |
|---|---|---|---|
| **4a. Typed flags + plain A/B split** *(Recommended)* | **\$0** | Cannot run overlapping multi-layer or multi-arm sweeps without extension | Low–Medium |
| **4b. Orthogonal layers + Latin Hypercube (LHS)** | **\$0** | Severely underpowered at low traffic; requires multiple-comparison adjustment and per-experiment re-salting | High |
| **4c. None yet (logging only)** | **\$0** | Cannot test UI or caching changes online | Low |

**Power calculation:** For a baseline 20% conversion metric, detecting a 5 percentage-point lift (20% → 25%) at 80% power and $\alpha = 0.05$ requires **~1,100 sessions per arm** (~2,200 total for a 2-arm A/B test). With $K = 10$ LHS arms, each arm receives $1/K$ of traffic (~11,000 sessions total). At early traffic (e.g., 100–300 sessions/week), a 2-arm test takes ~2–5 months, and a 10-arm LHS would take years.

**Recommendation:** **4a (Typed flags + plain A/B split)**, keeping **simulation physics out of online experiments** so saved worlds are never contaminated. Orthogonal layers and LHS remain a future extension once traffic justifies them.

### Decision 5: Identity, privacy, and retention

Separate identifiers with explicit lifetimes ensure cross-tab, cross-day, and multi-replica correctness:
- **`event_id`:** UUID per event (for deduplicating client retries).
- **`session_id` + scoped `seq`:** Random ID per browser tab/session (`sessionStorage`), with `seq` scoped to `(session_id, source)` (frontend vs. backend) so concurrent tabs or replicas never collide.
- **`subject_id` (pseudonymous browser ID):** Random first-party ID persisted in a dedicated cookie/`localStorage` across sessions and replicas (with a reset/opt-out mechanism), enabling cross-day retention analysis without logging GitHub logins or IPs. Treat `subject_id` as pseudonymous data (event timelines can still link behavior); if longitudinal linkage is deferred, omit `subject_id` and defer cross-day user dashboards.
- **`world_id` / `run_id`:** Simulation domain identifiers.
- **Experiment assignment unit:** `subject_id` for cross-session UI/UX tests, `session_id` for single-session latency/caching tests, and `world_id` for world-scoped behavior.
- **Never collected:** GitHub logins, IP addresses, free text, rendered UI label strings, custom world names, or raw seeds.
- **Scope & retention:** Telemetry is **off** in local/desktop builds (hosted mode only). Raw B2 chunks expire after 6 months via a B2 lifecycle rule after off-box daily Parquet compaction. Owner to confirm whether a minimal privacy notice or opt-out toggle is desired in the hosted UI.

## Design

### 1. Event envelope, UI ontology & domain events (Decision 3a)

Every event shares a compact envelope:
- `schema_version` (integer, starting at `1`), `ontology_version`
- `event_id`, `seq`, `source` (`frontend` | `backend`)
- `ts` (server UTC timestamp) and `client_ts` (browser timestamp for client events)
- `build` (git SHA), `mode` (`hosted`)
- `subject_id` (if enabled), `session_id`, `world_id`, `run_id` (when applicable)
- `event` (verb + target, or backend domain event name)

**UI Ontology & Lifecycle:**
- Interactive UI elements are registered in a single typed registry with a stable `element_id` (named by function, e.g., `start_run`), `kind` (`button`, `toggle`, `list_item`, `slider`, `modal`, `map`), and parent surface.
- Fixed verb set: `shown`, `clicked`, `changed`, `dismissed`.
- **No free-text labels:** Send a code-owned `label_key` / `surface_revision` and a coarse `layout_slot`—never rendered strings (which could leak custom world names and create high-cardinality dimensions).
- **Ontology evolution:** Registry entries are additive in the UI PR; CI fails on removal or rename vs. `origin/main`. When an element's semantics change meaningfully, record an explicit deprecation/replacement mapping and bump `ontology_version` / `build` metadata so dashboards can preserve or deliberately split the historical series.
- **Correctness enforcement:** `track()` accepts only the TypeScript union of registered IDs; interactive controls use a shared wrapper that requires an ID, paired with a lint rule that fails on raw click handlers outside the wrapper.
- **List & position bias:** List surfaces emit a `shown` event with presented item IDs and attach `list_id`, `index`, `length`, `ordering`, and `in_view` on `clicked`.

**Authoritative backend domain events:**
- UI clicks record *intent*, not what the backend accepted. The backend emits separate domain events (`simulation_choice_accepted`, `simulation_outcome`) with controlled `choice_key`/`choice_value`, `world_id`, `run_id`/`command_id`, `build`/`config_revision`, and `experiment_assignment_id`, linking each accepted choice to the resulting step/run outcome (without raw seeds or world names).
- Backend also emits lifecycle (`session_started`, `world_created`, `world_loaded`, `world_deleted`), compute (`step_completed`/`step_failed` with `client_disconnect`, `run_started`/`run_stopped`, `animation_started`/`animation_stopped`/`animation_recovered`, `worker_lifecycle`), rendering (`map_rendered`, `cross_section_rendered`), and `client_error` events.

**Browser ingest endpoint:**
- `POST /telemetry/events` accepts batched client events and uses `navigator.sendBeacon` on `pagehide` so end-of-session events are not lost. It enforces a per-request size limit, per-session rate limit, basic bot/crawler filtering, and schema/registry validation (rejecting and counting unknown IDs or malformed payloads).

### 2. Process lifecycle, rotation, B2 export & disk safety (Decision 1a)

- **Single writer & orphan recovery:** Under single-replica hosting, the front process owns the telemetry directory and receives forwarded events from the simulation worker (or sweeps any `active-<pid>` file whose PID is dead at startup and periodically). Files are **not** gzipped during the 2-second graceful shutdown window; the startup/periodic sweep handles any unsealed file left by a crash or 90-second worker idle stop.
- **Rotation:** Active files seal and gzip (off the event loop in a background thread) every 5 minutes, at 10 MB, or on worker idle. Readers skip any torn last line from a mid-write crash and increment a torn-line counter.
- **Two-phase B2 upload:** A background task uploads sealed `.jsonl.gz` files to `telemetry/v1/dt=YYYY-MM-DD/<filename>` using B2's `X-Bz-Content-Sha1` header and verifies the returned `contentSha1` before calling `unlink()`. A crash between upload and `unlink()` re-uploads the same deterministic filename (creating a harmless B2 version, expired by the lifecycle rule); the dashboard deduplicates by `event_id`.
- **Least-privilege credentials:** Railway uses a B2 application key restricted to the `telemetry/` prefix with write-only capability, so a web-service leak cannot read or delete world saves.
- **Disk safety & silent drops:** Local telemetry is capped at 50 MB (<\~\$0.01/month at ~\$0.15/GB-month) **and** checks free volume space so telemetry never starves checkpoints. Under disk pressure, new events are dropped (preserving pre-outage context). All telemetry failures are swallowed so user requests never fail; drop, eviction, and reject counters are emitted and uploaded.

### 3. Off-box cron & dashboard generation (Decision 2a)

A standalone script runs via cron on the analysis machine:
1. Tracks a high-water mark (last closed day synced) and lists only day-prefixes from there on, staying within B2's 2,500/day free listing/download allowance.
2. Compacts each closed day into a single Parquet file locally and deletes the raw local chunks.
3. Queries the Parquet files with DuckDB to generate a static HTML dashboard covering activation funnels, cross-day retention (if `subject_id` is enabled), choice-to-outcome conversion, view/inspector adoption, position-bias tables, `p50`/`p95` step/render latencies, worker cold-start tax, `seq` gaps, drop/reject counters, and experiment readouts.

### 4. Typed flags & plain A/B experimentation (Decision 4a)

- **Manifest & `Flags`:** A versioned JSON manifest defines active flags and plain A/B splits (hashing the designated assignment unit with a per-experiment salt). A typed `Flags` object provides safe production defaults for local runs/tests plus env/CLI overrides (`MB_FLAGS` / `--flag`) to reproduce any arm locally.
- **Exposures:** `experiment_exposed` is emitted when code reaches the flagged decision point. In-memory deduplication is purely a volume optimization; the dashboard deduplicates on `(unit_id, experiment_id)` taking the first exposure as the source of truth.
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
- **New:** `frontend/src/uiElements.ts` — append-only element registry, `label_key` / `ontology_version` metadata, typed `track()` helper, and tracked UI wrapper.
- **New:** `bin/telemetry_dashboard.py` — off-box B2 day-prefix sync, daily Parquet compaction, and DuckDB static HTML generator.
- **Existing:** `backend/app/hosted.py` (`_IDLE_SECONDS = 90.0`), `backend/app/hosted_worker.py`, `backend/app/main.py` — `POST /telemetry/events` ingest endpoint, authoritative domain events (`simulation_choice_accepted`, `simulation_outcome`), lifecycle/compute/render event emission, and `Flags` forwarding.
- **Existing:** `frontend/src/api.ts`, `frontend/src/App.tsx`, `frontend/src/cookies.ts` — `session_id` / `subject_id`, batched/beacon client reporting, and wrapped controls.
- **Existing (if physics flags are ever enabled):** `backend/app/world.py`, `backend/app/npz_save.py` (`thaw`), `backend/app/persistence.py`, `backend/app/replay_log.py`.
- **Docs:** `docs/architecture.md`, `docs/api-reference.md`, `docs/save-compatibility.md`.

## Assumptions and uncertainty

- Hosted service remains single-replica while using local-volume JSONL buffering (Option 1a); enabling >1 replica requires migrating ingestion off the volume first.
- Railway volume free-space headroom is sufficient for a 50 MB telemetry buffer alongside checkpoints.
- Hosted traffic remains below ~9,000 sessions/day compressed (at 6-month B2 lifecycle retention) and within B2's 2,500/day free download/listing calls.

## Acceptance criteria

- [ ] **Ontology & UI coverage:** CI fails if a registered `element_id` in `uiElements.ts` is removed or renamed vs. `origin/main`; lint fails on raw click handlers outside the tracked wrapper; no rendered label strings are sent (`label_key` only); a frontend test clicks every registered element and verifies the emitted event; list events include `list_id`, `index`, `length`, `ordering`, and `in_view`.
- [ ] **Domain event linkage:** Accepted simulation choices emit `simulation_choice_accepted` linked to `simulation_outcome` by `world_id` and `run_id`/`command_id` without logging raw seeds or free-form world names.
- [ ] **Identity & sequencing:** `event_id`, `session_id`, `(session_id, source)`-scoped `seq`, and `subject_id` (if longitudinal linkage is enabled) are generated and preserved across tabs/restarts without depending on a specific replica.
- [ ] **Ingest & abuse protection:** `POST /telemetry/events` rejects and counts oversized, malformed, and unknown-ID payloads, enforces per-session rate limits, and filters bots/crawlers.
- [ ] **Crash, replica & loss accounting:** An end-to-end test kills the worker mid-session and verifies the orphaned active file is sealed and uploaded; torn last lines are skipped and counted; the dashboard reports `seq` gaps and drop/eviction/reject counters. Single-replica guard (or multi-replica central-sink migration test before enabling >1 replica) is verified.
- [ ] **B2 & disk safety:** Unit tests verify `X-Bz-Content-Sha1` / `contentSha1` verification before `unlink()`, and new-event dropping when `MAX_LOCAL_BYTES` or free-space floor is hit without failing requests.
- [ ] **Cost ceiling:** Documented monthly cost stays within the owner's ceiling (\$0/mo incremental on B2 and <\$0.01/mo Railway volume up to ~4,500 sessions/day at 1-year retention or ~9,000 sessions/day at 6-month retention).
- [ ] **Flags & save compatibility:** Unit tests verify `Flags` defaults, `MB_FLAGS`/`--flag` overrides, and loading a real pre-change `.mbworld` save (via `npz_save.thaw`).
- [ ] **Dashboard:** `bin/telemetry_dashboard.py` compacts sample `.jsonl.gz` fixtures to Parquet and renders the static HTML report.
- [ ] **Repo gates:** `./bin/affected_test.sh --base origin/main`, `ruff check`, `./bin/quality_check.sh`, and frontend lint/build pass for each implementation PR.

## Falsification criteria

- If measured event loss (`seq` gaps + dropped/torn events) exceeds 5% under normal hosted operation, or if telemetry adds >2 ms `p95` latency to `/step` or `/render`, move serialization/compression fully out-of-process or sample high-volume events.
- If the bespoke pipeline's maintenance burden exceeds the owner's threshold or Railway requires >1 replica before the central-sink migration is built, switch to managed pay-as-you-go (Option 1b) within the agreed monthly budget ceiling.
- If hosted traffic cannot power even a two-arm A/B test on primary UI metrics within 8 weeks, defer online experimentation (Option 4c) and keep only typed flags + logging.

## References

- [`bin/b2_saves.py`](../bin/b2_saves.py) — existing B2 `X-Bz-Content-Sha1` upload and credential conventions.
- [`docs/architecture.md`](../docs/architecture.md) — backend/frontend request flow and `HostedWorker` lifecycle.
- [`docs/save-compatibility.md`](../docs/save-compatibility.md) — `.mbworld` backward-compatibility requirements.
