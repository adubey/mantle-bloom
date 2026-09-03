# Publishing mantle-bloom as a public website

A plan for taking mantle-bloom from a single-user localhost dev tool to a public site,
assuming **low traffic (< 1,000 monthly active users)**. Covers hosting, a low-cost
marketing plan, banner-ad integration and ad-network analysis, and a feature backlog
drawn from what worldbuilding communities actually ask for.

> Status: planning document. Nothing here is built yet. Order of work is roughly the
> order of the sections: make it safe to host (§3.4) → host it (§3) → market it (§4) →
> monetize if it's worth it (§5) → grow it (§6).

---

## 1. Executive summary

- **The core problem is not the frontend, it's the backend.** The React bundle is ~300 KB
  of static files and can go on any CDN for free. The backend is a CPU-heavy Python
  simulation (`numpy`/`scipy`/`numba`) that holds **exactly one world in memory at a
  time**, has **no database**, **no auth**, a **CORS allowlist pinned to localhost**, and
  a **`/world/load` route that unpickles arbitrary bytes (remote code execution)**. It was
  built for one trusted user on `127.0.0.1`. Every hosting decision flows from fixing that.
- **Recommended architecture:** static frontend on **Cloudflare Pages** (free) + a single
  **Hetzner dedicated-vCPU VPS** (~€14–24/mo) running the backend in Docker behind
  **Caddy** (automatic HTTPS) behind **Cloudflare** (free CDN/DDoS/caching), with a
  **precomputed gallery of curated worlds** served as static assets so most visitors never
  touch the compute box.
- **Required before launch (§3.4):** remove/replace the pickle load route, add
  multi-session world storage with TTL eviction, cap resolution / step count / world size
  for anonymous users, add a concurrency gate + queue, add rate limiting, lock CORS to the
  real domain.
- **Marketing:** content-led and community-led, ~$20–120 one-time (domain + a few assets).
  The killer asset is already in the codebase: `/world/animate` exports H.264 timelapses of
  a planet evolving — that is native front-page material for [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration),
  [r/worldbuilding](https://old.reddit.com/r/worldbuilding), Hacker News, YouTube, and short-form video.
- **Banner ads:** at < 1,000 MAU, programmatic banner ads will earn **roughly $1–15/month**
  and will visually degrade an image-first tool. Recommendation: **do not run a
  programmatic network yet.** Instead do (a) one tasteful **EthicalAds** or **Carbon Ads**
  slot in the sidebar if you want passive income now, and (b) pursue **direct sponsorship**
  from adjacent worldbuilding tools. Keep AdSense as a documented fallback. Details and a
  full network comparison in §5.
- **Features:** the single most-requested thing from worldbuilders is *"tell me the
  climate at this exact spot, across the year"* plus **useful exports** (16-bit heightmap,
  GeoTIFF, Azgaar/Wonderdraft-friendly) and **shareable seed URLs**. Full backlog in §6.

---

## 2. What we are actually deploying

From `docs/architecture.md` and `docs/profiling.md`:

| Piece | Reality | Hosting implication |
| --- | --- | --- |
| Frontend | React + TS, Vite build, **~300 KB total** (`frontend/dist`), no server-side rendering, talks to backend via `VITE_API_BASE` baked at build time | Trivial. Any static host / CDN, free tier. |
| Backend | Python, FastAPI + uvicorn, `numpy`/`scipy`/`numba`/`Pillow`/`av` | Needs a real always-on machine with CPU headroom. Not serverless-friendly (see §3.3). |
| World state | **One `World` in a module-level dict in `main.py`.** `/world/generate` replaces it. ~131 K nodes, large arrays, hundreds of MB resident for one big world. | Must become per-session with eviction before it's public, or accept "one shared world for everyone" (§3.1 Option A′). |
| Persistence | None. "Save World" = the client downloads an opaque pickle; "Load World" = the client uploads one and the server **unpickles it** = arbitrary code execution. | `/world/load` must be removed or hard-sandboxed before public exposure. |
| Compute cost | Generate: a few seconds. Step: **~3–5 s** at default resolution. Render: **~1.8–2 s**. `/world/animate` (60 frames): **~4.6–5 s/frame ≈ 5 min wall**, and permanently advances the world. | A single visitor clicking "Play" pins a core for seconds at a time. Concurrency must be gated. Animation is a heavy endpoint that needs its own queue + caps. |
| Rendering | Entirely server-side PNG (`render_image.py`); client just draws the image. Renders are a pure function of `(world state, view, projection, size, rotation)`. | Renders are **very cacheable** — same seed + same step + same view = same PNG. Big lever for keeping the compute box idle. |
| Auth / multi-tenancy | None. CORS = `localhost`. | Add a session cookie, real CORS, rate limiting. |

**Traffic sizing.** < 1,000 MAU is roughly 30–60 sessions/day, maybe **3–10 concurrent
users at a peak**, a handful of whom are actively stepping a simulation at any instant. A
single 4–8 vCPU box with a concurrency gate handles this comfortably. The risk is not
average load, it's (a) a traffic spike from a front-page post and (b) one user kicking off
several `/world/animate` calls.

---

## 3. Hosting options

### 3.1 Option A — Precomputed gallery + static site (cheapest, always do this part)

Run mantle-bloom **offline** on your own machine (or a cheap spot VM) to generate a
catalogue of 20–50 curated worlds: for each seed, render every view (elevation, biome,
combined, temperature, precipitation, wind, currents, resources, soil) at multiple
projections, plus a timelapse MP4 and a downloadable hex-grid export. Commit the outputs
(or push to object storage). The public "site" is then a **static gallery + world
explorer**: pick a world, flip views, scrub the timelapse, download assets. No live
backend at all.

- **Cost:** ~$0–5/mo (static hosting + a little object storage/bandwidth). Effectively
  free under 1,000 MAU on Cloudflare Pages / R2.
- **Pros:** zero attack surface, infinite scalability, no security work, survives a
  Hacker News spike without noticing.
- **Cons:** no "generate *your* world" — which is the entire appeal for worldbuilders.
- **Verdict:** **Build this regardless.** It is the front page, the marketing engine, and
  the fallback that renders when the compute box is busy or down. Options B/C add live
  generation *on top* of it.

**Option A′ — one shared live world.** A middle ground: keep the current single-world
backend, expose it read-mostly, let anyone step/generate but accept that everyone shares
one world (last click wins). Cheap, but griefable and confusing for multiple simultaneous
visitors. Only worth it as a very short-term "is there any interest?" probe.

### 3.2 Option B — Single VPS backend + static frontend on CDN (recommended)

```
                    ┌──────────────────────────┐
   visitor ───────► │ Cloudflare (free)        │  CDN cache for static assets +
                    │  DNS, TLS, DDoS, caching  │  rendered PNGs; WAF rate-limit rules
                    └────────────┬─────────────┘
                                 │
              ┌──────────────────┴───────────────────┐
              ▼                                      ▼
   ┌────────────────────┐               ┌──────────────────────────────┐
   │ Cloudflare Pages   │               │ Hetzner VPS (Docker)         │
   │  frontend/dist     │               │  Caddy → uvicorn (FastAPI)   │
   │  + gallery assets  │  API calls →  │  + Redis (sessions, queue)   │
   │  (free)            │               │  + gallery fallback assets   │
   └────────────────────┘               │  CCX23 ~€14–24/mo            │
                                        └──────────────────────────────┘
```

- **Frontend:** `frontend/dist` + the Option A gallery, deployed to **Cloudflare Pages**
  (free, global CDN, unlimited bandwidth, preview deploys per PR). Alternatives: Netlify,
  Vercel, GitHub Pages, Cloudflare R2 + a worker.
- **Backend:** one **Hetzner Cloud** VPS. Recommended: **CCX23** (4 dedicated vCPU / 16 GB,
  ~€24/mo) for predictable performance under the CPU-bound step loop; **CAX31** (Ampere ARM,
  8 vCPU / 16 GB, ~€14/mo) is cheaper and `numpy`/`scipy` wheels are fine on ARM64 — verify
  `numba` + `av` build cleanly in the Docker image first. Hetzner is roughly 2–3× the
  price/performance of the US hyperscalers for this kind of always-on CPU box.
- **Reverse proxy / TLS:** **Caddy** in the same compose file — automatic Let's Encrypt,
  one-line config, gzip/zstd, easy request timeouts.
- **In front:** **Cloudflare** (free plan) for the apex domain: CDN, TLS, basic DDoS
  protection, and — the important part — **cache rules on `/world/render`** keyed by the
  full query string, so a popular seed's PNGs are served from the edge and never hit
  Python twice. Add a WAF rate-limiting rule (free plan allows one) on the API paths.
- **Session/queue store:** **Redis** (or just in-process structures to start). Holds the
  per-session world registry, the step/animate job queue, and rate-limit counters.
- **Process model:** `uvicorn` with **1 worker** (the simulation is CPU-bound and the
  world objects are big; multiple workers each with their own copy of a world multiplies
  memory and doesn't help a single-core-bound step). Concurrency is handled by the job
  queue, not by worker count. Put heavy work (`/world/step`, `/world/animate`) on a
  background task / worker thread pool of size = (vCPU − 1) with a hard queue length limit;
  return `503 + Retry-After` (and a gallery world) when full.
- **Cost:** **~€14–30/mo** all-in (VPS + domain amortized). Everything else free-tier.
- **Pros:** real "generate your own world," full control, cheap, simple mental model,
  one machine to reason about.
- **Cons:** you now run a server — patching, monitoring, backups, an on-call-ish
  responsibility. Single point of failure (mitigated by the gallery fallback).

### 3.3 Option C — Managed PaaS (Fly.io / Railway / Render)

Same shape as B but the platform runs the box.

- **Fly.io:** `fly launch` from a Dockerfile, scale-to-zero possible, `shared-cpu-2x` /
  `performance-2x` machines. Good DX, pay-as-you-go; a `performance-2x` (2 dedicated
  vCPU / 4 GB) is ~$40–50/mo if always on, less with autostop. Watch memory for big worlds.
- **Railway:** simplest deploy, usage-based (~$5 credit + ~$10–20/mo realistic), 8 GB/
  8 vCPU ceiling on the hobby plan. Fine for this scale.
- **Render:** predictable flat tiers; a Standard instance (2 GB) is ~$25/mo, but 2 GB is
  tight for multiple resident worlds — you'd need the 4 GB tier (~$85/mo). Least
  cost-effective of the three here.
- **Cost:** ~$10–50/mo depending on platform and autostop.
- **Pros:** no server administration, easy rollbacks, built-in metrics/logs, scale knobs.
- **Cons:** 1.5–3× the cost of Hetzner for equivalent CPU; cold starts if you use
  scale-to-zero (the `numba` JIT warm-up and first-render cost make the first request
  after a cold start slow); memory ceilings on cheap tiers are a real constraint given
  world size.
- **Verdict:** good if you'd rather pay ~$20/mo extra than ever `ssh` into a box. Start
  here if server ops is not your thing; migrate to B if the bill or the constraints bite.

### 3.4 Serverless (Lambda / Cloud Run / Workers) — why it doesn't fit

- **Stateful in-memory world.** Serverless instances are ephemeral and load-balanced;
  "the world" would vanish between requests or live on a random instance. You'd need to
  serialize the entire multi-hundred-MB world to external storage on every request and
  reload it — dominant cost, and slow.
- **Execution limits.** A 5-minute `/world/animate` blows past typical request/CPU limits;
  even a single step plus render can approach short timeouts under load.
- **Cold starts.** `scipy` + `numba` + `av` is a large, slow-initializing runtime.
- **Cost inversion.** Per-second CPU billing for a workload that wants a core pinned for
  seconds per interaction is *more* expensive than a flat VPS, not less.
- Cloud Run with min-instances=1 and a big CPU allocation is the least-bad serverless
  option, but at that point it's just Option C with worse ergonomics for this app.

### 3.5 Comparison

| Option | ~$/mo @ <1k MAU | Live "generate your world"? | Ops burden | Spike resilience | Notes |
| --- | --- | --- | --- | --- | --- |
| A — Static gallery only | $0–5 | No | ~None | Excellent | Always build this. |
| A′ — One shared live world | $5–15 | Shared/griefable | Low | Poor | Interest probe only. |
| **B — VPS + CDN (recommended)** | **€14–30** | **Yes** | Medium | Good (w/ gallery fallback) | Best $/perf, full control. |
| C — Managed PaaS | $10–50 | Yes | Low | Good | Pay more to skip server admin. |
| D — Serverless | $30–100+ | Awkward | Medium | Good | Architecture fights you. Skip. |

### 3.6 Required code changes before going public

These are **blockers**, roughly in priority order:

1. **Kill the pickle RCE.** `POST /world/load` unpickles request bytes
   (`persistence.py`). On a public server this is remote code execution. Options: remove
   the route entirely for the public deployment (feature-flag it), or replace the
   save/load format with a safe explicit schema (JSON/msgpack of the fields that matter) —
   see `docs/api-reference.md`'s own note that it's "not a stable interchange format."
   Shipping the hex-grid export (already JSON) as the shareable artifact is a good interim.
2. **Multi-session worlds.** Replace the single module-level world with a
   `dict[session_id -> World]` plus:
   - session id from a `Secure`/`HttpOnly` cookie,
   - **TTL eviction** (e.g. drop worlds idle > 20–30 min),
   - a **hard cap on concurrent resident worlds** (e.g. 20–40 depending on box RAM);
     evict LRU when full, or refuse new generations with a friendly "server at capacity,
     explore a gallery world" response.
3. **Resource caps for anonymous users.** Clamp:
   - render `width`/`height` (public max well below the current 4000 px),
   - `node_density` / `climate_density` / `fluid_density` (offer only the lighter presets
     publicly; reserve `4.0` for a future account tier),
   - `/world/step` `years` per call and steps per minute per session,
   - `/world/animate` `num_frames` (public max ~24–30) and **one animation per session at
     a time**, globally queued.
4. **Concurrency gate + queue.** A bounded worker pool for `step`/`animate`/`generate`;
   `503 + Retry-After` + a gallery fallback when the queue is full. Prevents a spike from
   OOM-ing or thrashing the box.
5. **Rate limiting.** Per-IP and per-session, at Caddy and/or in-app (`slowapi` or a small
   Redis token bucket). Separate, stricter budget for `/world/animate`.
6. **Lock down CORS** (`main.py`) to the real frontend origin(s); drop the `localhost`
   entries in the public config via env var.
7. **Request timeouts** at Caddy so a wedged computation can't hold a connection forever.
8. **Observability:** structured logs, a `/healthz`, basic metrics (Prometheus text
   endpoint or just Cloudflare + `netdata`), and an error tracker (GlitchTip/Sentry free
   tier). Privacy-light analytics: **Plausible** (self-host or ~$9/mo) or **Cloudflare Web
   Analytics** (free, no cookie → no consent banner needed for analytics).
9. **Abuse/cost ceiling:** a global daily compute budget; when exceeded, the site goes
   gallery-only until it resets. Protects you from a bad day.

### 3.7 Nice-to-have, not blockers

- **Persistent render cache** on disk/R2 keyed by `(seed, params, elapsed_years, view,
  projection, size, rotation)` so restarts and popular seeds stay cheap.
- **Deterministic seed → shareable URL** (`/w/<seed>?step=…&view=…&proj=…&rot=…`) that
  reconstructs a world by replaying generation + N steps. Powers both sharing (§4) and a
  cheap cache key.
- **Pre-warm** the `numba` JIT on boot with a tiny throwaway world so the first real
  request isn't the slow one.
- **CI/CD:** GitHub Actions → build frontend to Pages, build+push backend image, deploy
  to the VPS over SSH (or `flyctl deploy` for Option C).

### 3.8 When to scale past one box

Add a second backend instance behind Cloudflare load-balancing (with **sticky sessions** —
each world lives on one instance) when *any* of: sustained queue wait > ~10 s at peak,
p95 step latency creeping up from contention, or RAM regularly > 80%. That's unlikely
below a few thousand MAU. The gallery absorbs read traffic long before the compute box
does.

---

## 4. Low-cost marketing plan

### 4.1 Positioning

> **mantle-bloom — a physically grounded planet simulator for worldbuilders.**
> Not a doodle tool. Plates actually collide and raise mountains; wind actually blows
> down real pressure gradients; rivers actually follow the terrain erosion carved. Give it
> a seed, watch a world grow, and pull out the maps you need.

Differentiators vs. the field (Azgaar's Fantasy Map Generator, Wonderdraft, Inkarnate,
World Anvil, Worldengine, ProcGenesis, Gleba, planet-hex, World Climate Lab):

- **Time dimension.** Most tools generate a static map. mantle-bloom *runs history* — you
  can watch and export the evolution.
- **Coupled physics.** Tectonics ↔ climate ↔ hydrology ↔ erosion feed back on each other,
  documented honestly in `docs/simulation-model.md`.
- **Depth of inspection.** Plate / River / Lake inspectors, Köppen biomes, pelagic ocean
  provinces, resource and soil maps.
- **Free and open source (GPLv3).**

### 4.2 Audience segments

| Segment | Where they are | What they want from this |
| --- | --- | --- |
| Fantasy/SF writers & worldbuilders | [r/worldbuilding](https://old.reddit.com/r/worldbuilding), [r/worldbuilding](https://old.reddit.com/r/worldbuilding) Discord, World Anvil, Bluesky #worldbuilding, Tumblr | A believable base map + climate they can defend; "why is there a desert here?" |
| TTRPG GMs / homebrew setting makers | [r/DnD](https://old.reddit.com/r/DnD), [r/DMAcademy](https://old.reddit.com/r/DMAcademy), [r/rpg](https://old.reddit.com/r/rpg), [r/dndmaps](https://old.reddit.com/r/dndmaps), [r/UnearthedArcana](https://old.reddit.com/r/UnearthedArcana), Foundry/Roll20 communities | A campaign continent with rivers, biomes, resources; export to VTT / battlemap tools |
| Procedural-generation hobbyists & devs | [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration), Hacker News, lobste.rs, #procjam, Bluesky #procgen | The method and the writeup; the animations; the code |
| Map-making / cartography hobbyists | [r/mapmaking](https://old.reddit.com/r/mapmaking), [r/cartography](https://old.reddit.com/r/cartography), CartographersGuild | High-res exports, projections, heightmaps to paint over |
| Earth-science / geography educators & enthusiasts | [r/geology](https://old.reddit.com/r/geology), [r/geography](https://old.reddit.com/r/geography), science-teacher communities | "What if Earth's tectonics were different" as a teaching sandbox |
| Game developers | [r/gamedev](https://old.reddit.com/r/gamedev), itch.io, Godot/Unity communities | Heightmap / hex-grid / biome exports as a world seed for their game |

### 4.3 Channels (all free or near-free)

**Reddit** (respect each sub's self-promo rules — participate first, post artifacts not
links-only, use the flair):
- [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration) — the most receptive; timelapse GIF + short method note.
- [r/worldbuilding](https://old.reddit.com/r/worldbuilding) — "I built a tool that simulates your world's tectonics and climate,"
  lead with a gorgeous biome map + the "click anywhere for its climate" feature.
- [r/mapmaking](https://old.reddit.com/r/mapmaking), [r/cartography](https://old.reddit.com/r/cartography) — the export story (heightmaps, projections).
- [r/dndmaps](https://old.reddit.com/r/dndmaps), [r/DMAcademy](https://old.reddit.com/r/DMAcademy) — "generate a plausible campaign continent in 2 minutes."
- [r/simulated](https://old.reddit.com/r/simulated), [r/EarthPorn](https://old.reddit.com/r/EarthPorn)-adjacent, [r/geology](https://old.reddit.com/r/geology) — the physics angle.
- [r/InternetIsBeautiful](https://old.reddit.com/r/InternetIsBeautiful) — once the site is polished and fast.

**Hacker News** — a **Show HN** with a title like *"Show HN: A from-scratch planet
simulator — plate tectonics, climate, erosion"* linking a page that loads instantly
(gallery first!) and a "How it works" writeup. This audience rewards honest engineering
writeups; `docs/simulation-model.md` is already 90% of a great post. Also **lobste.rs**
(`programming`, `graphics`) and **tildes.net**.

**Short-form video** — the `/world/animate` MP4s are tailor-made for **YouTube Shorts,
TikTok, Instagram Reels, Bluesky video**. "Watch a planet form in 30 seconds." Post a
new seed weekly. Low effort (the tool makes the content), potentially high reach.

**YouTube (long form)** — one 8–15 min "I simulated a planet from scratch" build-log /
explainer. Evergreen, feeds search, and is the thing other creators cite.

**Communities & fediverse** — worldbuilding and procgen Discords, Bluesky (#worldbuilding,
#procgen, #proceduralgeneration, #madewithscience), Mastodon (fosstodon, mastodon.gamedev),
the ProcJam community.

**itch.io** — publish the tool (or a link page) with devlogs; it has a real
worldbuilding-tools audience and its own discovery.

**Directories & wikis** — AlternativeTo, "best fantasy map generator" listicles
(Reedsy, GM Hub, etc. — reach out), the [r/worldbuilding](https://old.reddit.com/r/worldbuilding) wiki resources list, awesome-lists
(awesome-proceduralgeneration, awesome-gamedev), World Anvil / Cartographers Guild
resource threads.

**SEO** — target long-tail intent: "fantasy world map generator with climate,"
"procedural planet generator," "plate tectonics simulator online," "worldbuilding climate
map maker," "realistic river placement worldbuilding," "Köppen climate map generator."
Each gallery world is an indexable page; the "How it works" doc targets the method
queries. Add `sitemap.xml`, Open Graph / Twitter card images per world, JSON-LD.

### 4.4 The viral loop (highest-leverage single feature)

**Shareable world URLs + auto-generated share images.** `/w/<seed>?step=…&view=…` that
anyone can open to the exact world/view, plus a "Share this world" button that produces a
1200×630 OG image (the current render + a small caption). Every share on Reddit/Discord/
Bluesky then renders a picture and links back. Pair with a **"Seed of the week"** featured
on the front page and posted to socials. This is the difference between a launch bump and
compounding growth.

### 4.5 Marketing materials checklist

Cheap to make; several are one command away given the existing export features.

- [ ] **Domain** — something short and memorable (~$10–15/yr; `.world`, `.earth`, `.cc`,
      `.app` all fit). Biggest single spend.
- [ ] **Logo / wordmark** — simple; a stylized planet/terminator line. DIY in Figma/
      Inkscape, or ~$0 with a careful type choice, or the `design` skill.
- [ ] **Landing page** — gallery-first, hero = an autoplaying muted timelapse loop, one
      clear "Generate your own" CTA, a strip of 6 curated worlds, a short "how it's
      different" section, link to "How it works."
- [ ] **Curated world set** — 20–50 named seeds ("Kes­trel," "The Long Coast," …), each
      with all views + timelapse + downloads. Doubles as SEO surface and cache warm-set.
- [ ] **Timelapse video library** — 10–20 clips at launch, then weekly. Multiple aspect
      ratios (16:9, 9:16, 1:1) from the same seed.
- [ ] **"How it works" page** — adapted from `docs/simulation-model.md`, with diagrams.
      This is the Hacker News / lobste.rs artifact.
- [ ] **Comparison figure** — same world across elevation / biome / climate / resources,
      one image. Very shareable, communicates depth instantly.
- [ ] **Press/creator kit** — logo pack, 8–10 hi-res screenshots, 3–4 GIFs, one-liner,
      100-word blurb, feature list, contact. A single page.
- [ ] **OG/Twitter card images** — default + per-world (auto-generated, §4.4).
- [ ] **Animated GIFs** — for Reddit/forum posts that don't autoplay video.
- [ ] **README / repo polish** — screenshots at the top, a "Try it live" link, GPLv3
      badge, a CONTRIBUTING note. The repo is itself a marketing channel.
- [ ] **A short FAQ** — "Is it accurate?" (honest answer), "Can I use maps commercially?"
      (yes — clarify the output license vs. the GPLv3 code), "Will you add X?"

### 4.6 Launch sequence

1. **Weeks −4 to −1 (private):** ship §3.4 blockers; build the gallery; make the curated
   set + first video batch; write "How it works"; set up analytics + error tracking; get
   a few worldbuilder friends to try it and break it.
2. **Soft launch:** post in 1–2 Discords and [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration) only. Fix what breaks.
   Confirm the box survives modest concurrent load.
3. **Launch week:** Show HN (Tue–Thu morning US) → same day, [r/worldbuilding](https://old.reddit.com/r/worldbuilding) + [r/mapmaking](https://old.reddit.com/r/mapmaking)
   with *tailored* posts and native media → lobste.rs → Bluesky/Mastodon thread with the
   best timelapse. Watch the box; keep the gallery fallback armed.
4. **Sustain (ongoing):** "Seed of the week" post every week (cross-post video everywhere),
   respond to every comment/issue, ship one visible feature every 2–4 weeks and announce
   it, submit to directories and listicles, reach out to 2–3 worldbuilding YouTubers with
   a creator kit and a personalized world.

### 4.7 Budget

| Item | One-time | Recurring |
| --- | --- | --- |
| Domain | — | ~$12/yr |
| Hosting (§3, Option B) | — | ~€14–30/mo |
| Logo / brand | $0–80 (DIY vs. cheap commission) | — |
| Analytics | $0 (Cloudflare) or ~$9/mo (Plausible) | optional |
| Error tracking | $0 (free tier) | — |
| Optional: a boosted [r/worldbuilding](https://old.reddit.com/r/worldbuilding) or Reddit ad experiment | $50–150 | — |
| **Total** | **~$0–230** | **~$15–40/mo** |

Everything load-bearing here is free-tier + the VPS + the domain. The rest is time.

---

## 5. Banner ads

### 5.1 Revenue reality check

At **< 1,000 MAU** you're looking at maybe **3,000–15,000 pageviews/month** (worldbuilders
flip views and scrub timelapses, so pageviews-per-session is decent). Realistic display
RPM for a niche hobby/tools audience is **$1–5** (session RPM lower). That's **~$5–50/month
gross**, and small-publisher networks with fees/minimums net far less — often **$1–15/mo**
in practice. This will **not** cover the VPS, and a programmatic ad unit sitting next to a
carefully-rendered map materially cheapens the product. **Weigh that before adding any
network.**

### 5.2 Ad network analysis

| Network | Min. traffic to join | Fit for this site | Est. RPM (this niche) | UX / performance cost | Notes |
| --- | --- | --- | --- | --- | --- |
| **Google AdSense** | None | OK | $1–3 | Medium (needs a CMP + consent banner for EU; layout shift if not careful) | Easiest to get, near-100% fill, needs "sufficient content" (the gallery pages help). Policy: don't place next to interactive UI that could drive accidental clicks — keep it away from the map. The default fallback. |
| **Ezoic** | None (dropped its minimum) | OK-ish | $3–8 claimed (AI-optimized) | **High** — heavy script, added latency, aggressive placements by default, wants to sit on your DNS or via a Cloudflare integration | Best *revenue* option that accepts small sites, worst *fit* for a fast image tool. Only if ad income becomes a real goal and you cap placements hard. |
| **Media.net** | Low but wants meaningful US traffic | OK | $1–3 | Medium | Yahoo/Bing contextual. Approval is pickier than AdSense; underwhelming for non-US traffic. |
| **EthicalAds** | Accepts small/independent + open-source sites | **Good** | $1–2 (flat, honest) | **Very low** — one static image/text ad, no tracking, no consent banner needed | Dev/technical + open-source audience, privacy-respecting, matches the project's GPLv3 ethos. Lower ceiling, but the right *kind* of ad. Strong first choice for a passive slot. |
| **Carbon Ads** (BuySellAds) | Curated; wants a design/dev/creative audience and steady traffic | **Good aesthetically** | $0.10–0.30/click, varies | **Very low** — a single tasteful unit | Beautiful single-ad format used by design/dev tools. Hard to get in at low traffic; apply anyway once the site is polished. |
| **BuySellAds direct / direct sponsorship** | You set terms | **Best** | Whatever you negotiate ($50–300/mo for one slot is plausible in this niche) | You control it entirely | Sell one slot to an adjacent tool (Inkarnate, World Anvil, Foundry VTT, Dungeondraft/Wonderdraft, dice/map shops, Kanka, LegendKeeper, Azgaar donations). One relevant sponsor > a wall of programmatic junk. |
| **Monumetric** | 10k pageviews/mo + **$99 setup fee** under 80k | Poor | $3–6 | Medium | Fee makes no sense at this scale. Revisit at ~50k+ pv/mo. |
| **Mediavine** | 50k sessions/mo | N/A yet | $10–25 | Medium-High | The real goal *if* traffic 10×s. Not now. |
| **Raptive (AdThrive)** | 100k pageviews/mo | N/A yet | $12–30 | Medium-High | Same — a "someday" target. |
| **PropellerAds / Adsterra / PopAds / pop-under & push networks** | None | **Bad** | varies, "high" but predatory | **Severe** — pop-unders, redirects, malvertising risk, destroys trust | Accept anyone, and it shows. Do not use on a site you want worldbuilders to bookmark and share. |

### 5.3 Recommendation

1. **Launch with no ads.** Get the product and the traffic right first. An ad on day one
   signals "content farm" to exactly the communities you're courting.
2. **First money (passive):** a **single EthicalAds slot** in the sidebar/below the map,
   or **Carbon Ads** if accepted. One unit, never over the canvas, never on the active
   simulation view.
3. **Best money (active):** pitch **one direct sponsor** from the worldbuilding-tools
   ecosystem. A 250×250 or leaderboard "Sponsored by …" in the sidebar. Revisit quarterly.
4. **Documented fallback:** **AdSense** with a strict single-placement config + a free CMP
   (Cloudflare Zaraz consent, CookieYes, or Osano free tier) + a privacy policy page. Turn
   it on only if 2–3 don't materialize and you specifically want the pennies.
5. **Revisit programmatic seriously only at ~25k–50k pageviews/month**, at which point
   Monumetric/Ezoic math changes and Mediavine comes into view.

### 5.4 Better-than-ads monetization (higher ROI at this scale)

- **Print-on-demand posters** of your own world — the render pipeline already outputs
  high-res PNGs; wire a "Order a print" button to Printful/Prodigi/Gelato. Worldbuilders
  and GMs genuinely want their world on the wall. Likely out-earns ads per visitor.
- **Ko-fi / GitHub Sponsors / "buy me a VPS month"** — a small, honest donation link.
  Open-source + a visible hosting cost converts surprisingly well.
- **Optional "supporter" tier** (Patreon or one-off): higher-resolution worlds, longer
  animations, no queue, cloud-saved worlds, batch seed generation, private galleries.
  This also maps cleanly onto the §3.4 resource-cap tiers — free users get the light
  presets, supporters unlock `4.0`.
- **Commissioned worlds** — a paid "give me a world matching this brief" service, even at
  low volume, is real money and great marketing (each one becomes a showcase).
- These stack with a single tasteful sponsor slot without making the site feel like adware.

### 5.5 Implementation sketch (whichever slot you pick)

- One `<AdSlot>` React component, **lazy-mounted** (IntersectionObserver), fixed
  reserved height to prevent layout shift, rendered only in the sidebar and only on
  non-interactive views (gallery, world detail, "how it works") — **never while a
  simulation step or animation is running**.
- Respect `navigator.connection.saveData` and `prefers-reduced-data`; skip the slot on
  small viewports.
- Consent: if AdSense/Ezoic → a CMP; if EthicalAds/Carbon/direct → no tracking, no banner
  needed (still ship a privacy policy).
- Keep all ad scripts out of the critical render path; they must not delay first paint of
  the map. Cloudflare Zaraz can load third-party tags off the main thread.
- Add a `robots`/`ads.txt` as required by the chosen network.

---

## 6. Feature ideas from worldbuilding communities

Synthesized from [r/worldbuilding](https://old.reddit.com/r/worldbuilding), [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration), [r/mapmaking](https://old.reddit.com/r/mapmaking), the Azgaar
community, and comparable tools (Azgaar's FMG, Wonderdraft, Worldengine, Gleba,
planet-hex, World Climate Lab, ProcGenesis). Grouped by theme; rough effort (S/M/L) and
impact (★–★★★) are guesses to help sequencing.

### 6.1 Make the existing data *usable* elsewhere — highest impact

The recurring complaint about physics-based generators is "beautiful, but I can't get my
world *out* of it." mantle-bloom already has a hex-grid export; extend the export story.

| Feature | Effort | Impact | Notes |
| --- | --- | --- | --- |
| **16-bit greyscale heightmap PNG** (equirectangular) | S–M | ★★★ | The single most-requested export across every map community. Paint-over in Wonderdraft/Photoshop/Wilbur/World Machine, import to Unreal/Unity/Godot/Blender. |
| **GeoTIFF / raster with real georeferencing** | M | ★★ | The GIS/QGIS crowd explicitly asks for this. Elevation, plus optional biome/temp/precip bands. |
| **Blender/SpaceEngine-ready pack** — equirectangular color + height + normal + biome id maps at 4k/8k | M | ★★★ | Turns a world into a renderable 3D globe. Very shareable output → marketing flywheel. |
| **Azgaar / Wonderdraft-friendly outputs** — coastline SVG, heightmap in their expected ranges, a `.json` of continents/rivers | M–L | ★★ | Meet the incumbent tools where they are; "start in mantle-bloom, finish in Azgaar." |
| **Vector coastline + rivers + lake polygons (GeoJSON/SVG)** | M | ★★ | Cartographers want editable linework, not just pixels. |
| **Per-region crop + high-res re-render** | M | ★★★ | "Zoom into this continent for my campaign map." Currently only whole-planet. |
| **Data download for a clicked cell / region** (elevation, climate, biome, soil, resources as CSV/JSON) | S | ★★ | GMs statting out a region; worldbuilders defending choices. |
| **Named-feature export** — if labeling ships (§6.4), include labels in every export | S | ★ | Azgaar users specifically ask for labels-in-export. |

### 6.2 Climate & environment inspection — the worldbuilder's core need

The #1 question in [r/worldbuilding](https://old.reddit.com/r/worldbuilding) is *"what's the climate here and is it plausible?"*

| Feature | Effort | Impact | Notes |
| --- | --- | --- | --- |
| **Click-anywhere climate readout** — monthly temperature & precipitation graph, Köppen code, hemisphere-aware seasons, "why" (latitude, elevation, rain shadow, continentality, current) | M–L | ★★★ | The feature that would make worldbuilders switch. Some data exists in `climate.py`; needs a seasonal cycle and a point-query UI. |
| **Prevailing wind & ocean-current overlay with labeled gyres** | S–M | ★★ | Already have wind/current arrow views; add gyre detection + labels + a cleaner cartographic style. Community teaches "draw your gyres" by hand — do it for them. |
| **Latitude / climate-band guide overlay** (Hadley/Ferrel/Polar, ITCZ, horse latitudes) | S | ★★ | Educational and reassuring; helps users sanity-check deserts and rainforests. |
| **Seasonal scrubber** — step through a year, watch the ITCZ, monsoons, sea ice, snow line move | L | ★★★ | Huge "wow," directly asked for ("show me summer vs winter"). Depends on a seasonal climate solve. |
| **Hemisphere/seasons toggle** on climate maps | S | ★★ | "It's summer in the north" — cheap clarity win. |
| **Biome hover/click identify with legend** everywhere (not just the `combined` view) | S | ★★ | Search results show this is explicitly wanted; you already do it in one view. |
| **"Is this river navigable / how big is this watershed"** in the River Inspector | S–M | ★ | GMs love navigable-river data for trade routes. Flow data already exists. |
| **Habitability / agriculture overlay** — growing season length, frost dates, arable land | M | ★★ | Feeds the "where would people settle / farm" question. Soil-quality view is a start. |

### 6.3 Generation controls & scenarios — "what if"

| Feature | Effort | Impact | Notes |
| --- | --- | --- | --- |
| **Scenario presets** — Pangaea / archipelago world / ice age / hothouse / waterworld / young violent tectonics / old quiet world | S–M | ★★★ | Low effort over existing sliders, high discoverability, great for social posts ("I generated an ice-age world"). |
| **Planet parameters** — radius / gravity / day length / stronger or weaker sun / higher axial tilt / low obliquity | M–L | ★★ | [r/worldbuilding](https://old.reddit.com/r/worldbuilding) and [r/SpaceEngine](https://old.reddit.com/r/SpaceEngine) love "super-earth" and "tilted 40°" experiments. Axial tilt + solar multiplier already exist. |
| **Tidally-locked / binary-star / eyeball planet** mode | L | ★★ | A whole content genre on its own. Big climate-model change. |
| **Guided континент placement / "draw your landmasses" seed** | L | ★★★ | The bridge between "roll a random planet" and "I have a specific map in mind." The most common reason people bounce off physics-based generators. |
| **Constrain generation** — "give me a world with a big equatorial continent and a polar ocean," reroll until it fits | M | ★★ | Rejection-sampling over existing params; big UX payoff. |
| **Reproducible seed + full parameter string in the URL** | S | ★★★ | Also the §4.4 viral loop and the §3.7 cache key. Do this early. |

### 6.4 Names, labels & lore hooks

| Feature | Effort | Impact | Notes |
| --- | --- | --- | --- |
| **Auto-name features** — oceans, seas, continents, major mountain ranges, big rivers, major lakes (procedural name generator, swappable language styles) | M | ★★★ | Azgaar's naming is a huge part of its appeal. Even rough names make a world feel real and screenshots shareable. |
| **Pin & annotate** — drop labeled markers, export coordinates | S–M | ★★ | Minimum viable "this is my world" layer without a full political-map editor. |
| **Editable label layer + styled map export** ("poster mode": pick projection, palette, labels, legend, title cartouche) | L | ★★★ | The output GMs actually frame/print. Ties into print-on-demand (§5.4). |
| **Suggested settlement sites** — coastal + river + arable + defensible scoring | M | ★★ | "Where would cities be?" is a top-5 worldbuilding question. |

### 6.5 Presentation, sharing & community

| Feature | Effort | Impact | Notes |
| --- | --- | --- | --- |
| **Shareable world URLs + auto OG images** | S–M | ★★★ | §4.4. Single highest-leverage growth feature. |
| **Public gallery with "Seed of the week" + user-submitted seeds** | M | ★★ | Community surface; also the §3.1 fallback content. |
| **3D globe / spinning-planet view** (orthographic is already a projection candidate) | M–L | ★★★ | Enormously shareable; "show me the globe" is constant. Even a slow client-side WebGL sphere textured with the equirectangular render. |
| **More projections** — Robinson, orthographic/globe, Mercator-ish, Winkel tripel, polar azimuthal, hex | S–M each | ★★ | Behrmann + Eckert IV today; cartographers want options, especially globe and polar. |
| **Embeddable viewer (iframe)** for World Anvil / campaign wikis / Notion | M | ★★ | Distribution into the tools worldbuilders already pay for. |
| **Timelapse export presets** (already have `/world/animate`) — aspect-ratio picker, "cinematic" easing, view-crossfade | S | ★★ | Turns the existing feature into a content machine (§4.3). |
| **Side-by-side seed compare** | S | ★ | "Which of these two worlds should I use?" |
| **Accounts + cloud-saved worlds** (replaces the pickle download) | M–L | ★★ | Also resolves the §3.4 security blocker and enables the supporter tier (§5.4). |
| **Mobile-friendly read-only viewer** | M | ★★ | A lot of Reddit traffic is on phones; at least the gallery + world detail should work well. |

### 6.6 Longer shots (interesting, lower priority for < 1k MAU)

- **Plate-history playback** — "this mountain range came from this collision 40 Myr ago."
- **Paleo-map export** — the world at N intermediate epochs, for deep-history settings.
- **Star field / moons / tides / eclipses** — asked for, but a big scope jump.
- **Simple civilization/biome-carrying-capacity pass** — population potential, trade
  routes; borders on being a different product (that's Azgaar's turf).
- **API access** for game devs (rate-limited, key-gated) — natural supporter-tier perk.
- **Localization** of the UI — later, once there's a reason.

### 6.7 Suggested near-term slice

If the goal is "launch a public site that worldbuilders actually adopt," the
highest–return-per-effort bundle is:

1. Reproducible **seed+params in URL** + **shareable OG images** (§6.5, §4.4).
2. **16-bit heightmap** + **clicked-region data** export (§6.1).
3. **Click-anywhere climate readout** with a monthly graph and Köppen code (§6.2).
4. **Scenario presets** (§6.3).
5. **Auto-names** for oceans/continents/major ranges and rivers (§6.4).
6. **3D globe view** (§6.5) — mostly for the marketing flywheel.

That set is squarely in the existing model's wheelhouse, needs no new physics except the
seasonal cycle for (3), and each item doubles as marketing material.

---

## 7. Rough monthly cost summary

| Scenario | Monthly |
| --- | --- |
| Gallery-only (Option A) | **~$0** (Cloudflare Pages/R2 free tiers) + ~$1/yr domain amortized |
| Recommended (Option B: Hetzner VPS + Cloudflare + Pages) | **~€15–30** + ~$1/mo domain |
| Managed PaaS (Option C) | **~$15–50** |
| Add Plausible analytics (optional) | +$0 (Cloudflare) or +$9 (Plausible cloud) |
| Ads at < 1k MAU | **−$1 to −$15 income** (does not cover hosting; see §5) |
| Print-on-demand / donations / supporter tier | variable; realistically the first thing that actually offsets the VPS |

**Bottom line:** budget **~$15–30/month** to run this properly for < 1,000 MAU, treat ad
revenue as negligible, and lean on the gallery + shareable worlds + timelapse videos as
both the scaling strategy and the marketing engine.
