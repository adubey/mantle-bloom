---
name: documentation
description: Write or revise Mantle-Bloom project documentation in docs/ or rfc/. Use for documentation edits; RFC proposals and reviews should also follow their specialized skills.
---

# Mantle-Bloom documentation

Write for two readers: people learning or checking the system, and agents that need to find
an authoritative description while changing code. Documentation should describe current
behavior accurately and make the connection to implementation easy to verify.

## Choose the home

- `docs/` describes behavior that exists now. `rfc/` describes a proposed change and its
evidence; use `rfc-authoring` or `rfc-review` for RFC work.
- `docs/simulation-model.md` is a navigation page. Put model details in the relevant topic
page under `docs/simulation-model/`; its topic list maps modules to pages.
- `docs/architecture.md` maps modules and request flow; `docs/api-reference.md` documents
routes; `docs/debugging.md` holds reusable diagnostic procedures. Avoid maintaining the
same explanation in several places.
- Check the relevant page, implementation and linked sources before editing. Verify symbols,
defaults, units, limits and control flow in code rather than inferring them from old prose.

## Write for current behavior

Lead with what the system does now, then explain why, its limits and how a reader can inspect
or reproduce it. Keep historical rationale only when it explains a current choice; link to
the RFC or issue for the development history. Avoid dated section headings, “since removed”
and “an earlier version did” narratives, repeated claims of development-time confirmation,
stacked em-dash asides, and intensifiers such as “genuinely”, “exactly” or “actually”.

For agent readers, preserve precise names, units, defaults, formulas, state ownership and
links to the implementing symbols. When code renames or removes a symbol, changes a default
or changes a value, search the docs (for example, `rg 'OldSymbol|old_value' docs rfc`) and
update affected explanations and links. Search code comments too when a page moves. Do not
leave a plausible but stale description that an agent could mistake for current behavior.

## Diagrams and figures

Use a diagram when it clarifies a spatial relation, sequence, hierarchy or data transfer.
Prefer Mermaid for small relationships and flows because it is readable as source by both
people and agents. Use an SVG for a stable spatial schematic and rendered maps or plots for
actual simulation output.

Prose must still explain the claim without requiring the reader to open an image. A caption
and alt text should state what the figure shows and the point it supports; label conceptual
figures as conceptual. Include seed, time, projection and other context for simulation output.
For map images from a saved world, use `GET /world/render?view=...` or a save-aware renderer
such as `bin/debug/render_craton_map.py --load`; the route and examples are indexed in
`docs/debugging.md`. Record the save or seed, elapsed time and view. Keep
figures compact and legible at page size, and keep source data or generation steps for
quantitative claims. Store assets near their docs; keep large saves and replay outputs out of
the repository.

## Finish

Check changed links, anchors and image paths, then search for moved headings and renamed
symbols/values that may now be stale. Read the Markdown as a page: headings, lists, tables,
code and diagrams should remain easy to scan. Follow `AGENTS.md` for project workflow and
simulation-issue requirements; those requirements do not apply to descriptions of existing
behavior merely because they concern the simulation.
