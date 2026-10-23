# Hosting: current setup and next decisions

This page records Mantle Bloom's current hosted shape and the boundary of that design. The
full launch goals and open product decisions are in [launch readiness](launch-readiness.md).
Local development and localhost sign-in are covered by [local hosting](../docs/local-hosting.md).

## Current deployment

The public service runs on Railway from the repository's root [Dockerfile](../Dockerfile),
configured as code in [.railway/railway.ts](../.railway/railway.ts). The multi-stage image
builds the React frontend and installs the pinned Python runtime, then one Uvicorn process
serves both from the same origin. Railway builds from `main` and uses `GET /healthz` as its
health check. The Docker image is portable; Railway is the only deployment configured in
this repository.

The service is a small hosted preview, not a multi-user world service:

- One configured GitHub account has the full app and API.
- Logged-out visitors can use the public viewer and its limited step allowance.
- Everyone sees and steps the same active world. The world lives in one process's memory;
  there are no separate user-owned worlds or general user accounts.
- A Railway volume stores the latest world checkpoint and run status. On restart or wake, the
  process restores that state. If the service sleeps, simulation pauses until a request wakes
  it; elapsed time is the time simulated, not wall-clock time.
- One replica and one Uvicorn worker are required by the in-memory world and single-writer
  checkpoint. This setup does not provide horizontal scaling or independent world workers.

The route and access details are documented in the [API reference](../docs/api-reference.md).
The checkpoint and startup behavior are in [architecture](../docs/architecture.md#run-and-checkpoint).

## Configuration

Infrastructure settings belong in `.railway/railway.ts`; review them with `railway config plan`
and apply them with `railway config apply`. Secrets and environment-specific values are set in
Railway's variables, not committed to the repository:

- GitHub OAuth client ID and secret
- The numeric GitHub user ID allowed into the full app
- Session signing secret
- Public service URL

The Railway service mounts its persistent volume at `/data`. `RAILWAY_VOLUME_MOUNT_PATH` is
provided by Railway, and the app uses it for checkpoints. The service definition requires the
mount and its health check verifies that the checkpoint directory can be written. See
[checkpoint.py](../backend/app/checkpoint.py) and [architecture](../docs/architecture.md#run-and-checkpoint)
for the behavior.

## Hosting boundary

The Docker image keeps the application build portable, while the deployment contract is still
Railway-specific: one replica, one shared world, and one persistent checkpoint volume. Adding
accounts, paid access, per-user worlds, concurrent workers, a database or object store, or
sharding requires an explicit storage and execution design. The launch roadmap tracks those
decisions under [security and public service](launch-readiness.md#4-security-and-public-service)
and [performance and capacity](launch-readiness.md#5-performance-and-capacity).

The current service has no published capacity or availability target. Set those targets and
choose which hosting responsibilities to keep before designing scale-out; see the launch
roadmap's [decisions to settle](launch-readiness.md#decisions-to-settle-before-calling-this-a-launch-plan).

## Marketing and advertising notes

The sections below preserve the earlier audience, channel, launch, and monetization research.
They are planning hypotheses, not a launch commitment. Audience fit should be tested with
small outreach; ad-network eligibility, policies, pricing, and revenue estimates can change
and must be checked with providers before acting. Hosting cost figures are scenario estimates,
not a quote for the current Railway service. Recheck product claims against the current app
before using them in public copy.

## Marketing research and planning

### Positioning

> **mantle-bloom — a planet simulator for worldbuilders.**
> Explore how tectonics, climate, hydrology and erosion shape a seeded world over time, inspect
> its maps, and export what you need. The model is physically motivated and remains under
> development; describe its current limits alongside its strengths.

Differentiators vs. the field (Azgaar's Fantasy Map Generator, Wonderdraft, Inkarnate,
World Anvil, Worldengine, ProcGenesis, Gleba, planet-hex, World Climate Lab):

- **Time dimension.** Most tools generate a static map. mantle-bloom *runs history* — you
  can watch and export the evolution.
- **Coupled physics.** Tectonics ↔ climate ↔ hydrology ↔ erosion feed back on each other,
  documented honestly in `docs/simulation-model.md`.
- **Depth of inspection.** Plate / River / Lake inspectors, Köppen biomes, pelagic ocean
  provinces, resource and soil maps.
- **Free and open source (GPLv3).**

### Audience segments

| Segment | Where they are | What they want from this |
| --- | --- | --- |
| Fantasy/SF writers & worldbuilders | [r/worldbuilding](https://old.reddit.com/r/worldbuilding), r/worldbuilding Discord, World Anvil, Bluesky #worldbuilding, Tumblr | A believable base map + climate they can defend; "why is there a desert here?" |
| TTRPG GMs / homebrew setting makers | [r/DnD](https://old.reddit.com/r/DnD), [r/DMAcademy](https://old.reddit.com/r/DMAcademy), [r/rpg](https://old.reddit.com/r/rpg), [r/dndmaps](https://old.reddit.com/r/dndmaps), [r/UnearthedArcana](https://old.reddit.com/r/UnearthedArcana), Foundry/Roll20 communities | A campaign continent with rivers, biomes, resources; export to VTT / battlemap tools |
| Procedural-generation hobbyists & devs | [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration), Hacker News, lobste.rs, #procjam, Bluesky #procgen | The method and the writeup; the animations; the code |
| Map-making / cartography hobbyists | [r/mapmaking](https://old.reddit.com/r/mapmaking), [r/cartography](https://old.reddit.com/r/cartography), CartographersGuild | High-res exports, projections, heightmaps to paint over |
| Earth-science / geography educators & enthusiasts | [r/geology](https://old.reddit.com/r/geology), [r/geography](https://old.reddit.com/r/geography), science-teacher communities | "What if Earth's tectonics were different" as a teaching sandbox |
| Game developers | [r/gamedev](https://old.reddit.com/r/gamedev), itch.io, Godot/Unity communities | Heightmap / hex-grid / biome exports as a world seed for their game |

### Channels

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

### Shareable worlds

**Shareable world URLs + auto-generated share images.** `/w/<seed>?step=…&view=…` that
anyone can open to the exact world/view, plus a "Share this world" button that produces a
1200×630 OG image (the current render + a small caption). Every share on Reddit/Discord/
Bluesky then renders a picture and links back. Pair with a **"Seed of the week"** featured
on the front page and posted to socials. This is the difference between a launch bump and
compounding growth.

### Marketing materials checklist

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
- [ ] **OG/Twitter card images** — default + per-world (auto-generated, the shareable-world section).
- [ ] **Animated GIFs** — for Reddit/forum posts that don't autoplay video.
- [ ] **README / repo polish** — screenshots at the top, a "Try it live" link, GPLv3
      badge, a CONTRIBUTING note. The repo is itself a marketing channel.
- [ ] **A short FAQ** — "Is it accurate?" (honest answer), "Can I use maps commercially?"
      (yes — clarify the output license vs. the GPLv3 code), "Will you add X?"

### Launch sequence

1. **Weeks −4 to −1 (private):** ship launch-readiness blockers; build the gallery; make the curated
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

### Initial budget sketch

| Item | One-time | Recurring |
| --- | --- | --- |
| Domain | — | ~$12/yr |
| Hosting (hypothetical VPS option from the earlier comparison) | — | ~€14–30/mo |
| Logo / brand | $0–80 (DIY vs. cheap commission) | — |
| Analytics | $0 (Cloudflare) or ~$9/mo (Plausible) | optional |
| Error tracking | $0 (free tier) | — |
| Optional: a boosted [r/worldbuilding](https://old.reddit.com/r/worldbuilding) or Reddit ad experiment | $50–150 | — |
| **Total** | **~$0–230** | **~$15–40/mo** |

This estimate assumed a VPS plus a static gallery. It is not a current hosting quote.

---

## Advertising and monetization research

### Revenue estimates

At **< 1,000 MAU** you're looking at maybe **3,000–15,000 pageviews/month** (worldbuilders
flip views and scrub timelapses, so pageviews-per-session is decent). Realistic display
RPM for a niche hobby/tools audience is **$1–5** (session RPM lower). That's **~$5–50/month
gross**, and small-publisher networks with fees/minimums net far less — often **$1–15/mo**
in practice. This will **not** cover the VPS, and a programmatic ad unit sitting next to a
carefully-rendered map materially cheapens the product. **Weigh that before adding any
network.**

### Ad network comparison

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

### Recommendation

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

### Other monetization ideas

- **Print-on-demand posters** of your own world — the render pipeline already outputs
  high-res PNGs; wire a "Order a print" button to Printful/Prodigi/Gelato. Worldbuilders
  and GMs genuinely want their world on the wall. Likely out-earns ads per visitor.
- **Ko-fi / GitHub Sponsors / "buy me a VPS month"** — a small, honest donation link.
  Open-source + a visible hosting cost converts surprisingly well.
- **Optional "supporter" tier** (Patreon or one-off): higher-resolution worlds, longer
  animations, no queue, cloud-saved worlds, batch seed generation, private galleries.
  This also maps cleanly onto the future free and paid resource tiers — free users get the light
  presets, supporters unlock `4.0`.
- **Commissioned worlds** — a paid "give me a world matching this brief" service, even at
  low volume, is real money and great marketing (each one becomes a showcase).
- These stack with a single tasteful sponsor slot without making the site feel like adware.

### Ad implementation notes

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
