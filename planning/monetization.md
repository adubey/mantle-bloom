# Advertising and monetization plan

This preserves the earlier ad-network comparison and monetization hypotheses. Traffic,
eligibility, policy, pricing, and revenue estimates can change; verify them with providers
before committing money or implementation time. The hosting and launch context is in
[the hosting overview](../docs/HOSTING.md) and [launch readiness](launch-readiness.md).

## Revenue estimates

At **< 1,000 MAU** you're looking at maybe **3,000–15,000 pageviews/month** (worldbuilders
flip views and scrub timelapses, so pageviews-per-session is decent). Realistic display
RPM for a niche hobby/tools audience is **$1–5** (session RPM lower). That's **~$5–50/month
gross**, and small-publisher networks with fees/minimums net far less — often **$1–15/mo**
in practice. This will **not** cover the VPS, and a programmatic ad unit sitting next to a
carefully-rendered map materially cheapens the product. **Weigh that before adding any
network.**

## Ad network comparison

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

## Recommendation

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

## Other monetization ideas

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

## Ad implementation notes

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