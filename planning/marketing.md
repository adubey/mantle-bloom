# Marketing plan

These are exploratory audience, channel, and launch ideas, not a committed schedule. Test
audience fit with small outreach and recheck claims about current product capabilities before
using them in public copy. Hosting and launch dependencies are summarized in
[launch readiness](launch-readiness.md) and [the hosting overview](../docs/hosting.md).

## Positioning

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

## Audience segments

| Segment | Where they are | What they want from this |
| --- | --- | --- |
| Fantasy/SF writers & worldbuilders | [r/worldbuilding](https://old.reddit.com/r/worldbuilding), r/worldbuilding Discord, World Anvil, Bluesky #worldbuilding, Tumblr | A believable base map + climate they can defend; "why is there a desert here?" |
| TTRPG GMs / homebrew setting makers | [r/DnD](https://old.reddit.com/r/DnD), [r/DMAcademy](https://old.reddit.com/r/DMAcademy), [r/rpg](https://old.reddit.com/r/rpg), [r/dndmaps](https://old.reddit.com/r/dndmaps), [r/UnearthedArcana](https://old.reddit.com/r/UnearthedArcana), Foundry/Roll20 communities | A campaign continent with rivers, biomes, resources; export to VTT / battlemap tools |
| Procedural-generation hobbyists & devs | [r/proceduralgeneration](https://old.reddit.com/r/proceduralgeneration), Hacker News, lobste.rs, #procjam, Bluesky #procgen | The method and the writeup; the animations; the code |
| Map-making / cartography hobbyists | [r/mapmaking](https://old.reddit.com/r/mapmaking), [r/cartography](https://old.reddit.com/r/cartography), CartographersGuild | High-res exports, projections, heightmaps to paint over |
| Earth-science / geography educators & enthusiasts | [r/geology](https://old.reddit.com/r/geology), [r/geography](https://old.reddit.com/r/geography), science-teacher communities | "What if Earth's tectonics were different" as a teaching sandbox |
| Game developers | [r/gamedev](https://old.reddit.com/r/gamedev), itch.io, Godot/Unity communities | Heightmap / hex-grid / biome exports as a world seed for their game |

## Channels

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

## Shareable worlds

**Shareable world URLs + auto-generated share images.** `/w/<seed>?step=…&view=…` that
anyone can open to the exact world/view, plus a "Share this world" button that produces a
1200×630 OG image (the current render + a small caption). Every share on Reddit/Discord/
Bluesky then renders a picture and links back. Pair with a **"Seed of the week"** featured
on the front page and posted to socials. This is the difference between a launch bump and
compounding growth.

## Marketing materials checklist

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

## Launch sequence

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

## Initial budget sketch

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
