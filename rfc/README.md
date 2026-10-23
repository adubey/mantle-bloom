# RFCs

An RFC proposes a change to the simulation model that is too large, too uncertain or too
tied to the science to start as a plain issue. It explains the problem, compares solutions
against the literature, works out the design to the level of units, knobs and source files,
and ends with criteria that say when the change is done. Once accepted, an RFC is
implemented through ordinary issues and PRs (see [AGENTS.md](../AGENTS.md)), each citing it.

## Index

| RFC | Title | Status |
|---|---|---|
| [0001](0001-merge-cap-craton-loss.md) | Plate merges must not destroy cratons | Accepted |
| [0002](0002-craton-trench-attrition.md) | Cratons resist sinking: margin-dependent attrition and collision sinks | Draft |
| [0003](0003-rifting-around-cratons.md) | Rifts go around cratons | Withdrawn (superseded by 0009, 0010) |
| [0004](0004-proto-cratons.md) | Proto-cratons: restarting craton formation | Withdrawn (superseded by 0009) |
| [0005](0005-continental-recycling-balance.md) | Balance continental recycling at trenches | Draft |
| [0006](0006-intra-oceanic-arcs.md) | Intra-oceanic arcs: new continent from the ocean | Draft |
| [0007](0007-mantle-plumes-and-plateaus.md) | Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks | Draft |
| [0008](0008-craton-collision-jamming.md) | Cratons jam collisions: braking and fusion | Accepted |
| [0009](0009-craton-mosaic.md) | Cratons first: continents built from cratonic cores, plates placed around them | Accepted (revised after G1) |
| [0010](0010-rifts-follow-weak-lithosphere.md) | Rifts follow weak lithosphere | Accepted |
| [0011](https://github.com/adubey/mantle-bloom/blob/claude/rfc-0011-contact-strength/rfc/0011-contact-strength.md) | Collisions are held by lithospheric strength | Draft |
| [0012](0012-physics-first-simulation.md) | Physics-first simulation | Accepted |
| [0014](0014-lake-water-volume-and-freezing.md) | Lakes conserve their water: in-step fill-spill-merge, seasonal freezing and water-surface temperature | Accepted |
| [0015](0015-frozen-water-in-physical-units.md) | Frozen water in physical units: glacier mass balance, ice flow and meltwater, freezing lakes and rivers, ocean ice and polar snow | Draft |

RFCs 0001–0007 come from one investigation: a 19.2 Gyr run (seed 456862899) that lost
almost all of its continental crust and all of its cratons. The run is used as a stress
test of the model's long-term behaviour, not as a calibration interval. Its generating
commit is unrecorded (first published at `d795a09`); it predates #358 and used a 200 Myr
split cooldown for most of its 10 Myr-step history. The current evidence is the #403
refresh on `2562f63`, with 5 Myr steps and the year-based 2 Myr cooldown; see
[`analysis/issue403/`](../analysis/issue403/README.md). Each RFC preserves the old figures
as historical evidence and identifies the new comparator.

The RFCs are numbered by how directly they act on that loss, not by priority. The
refreshed implementation order (subject to each RFC's acceptance and diagnosis gates):
1. 0001's merge-placement diagnosis and selected fix first. Frequent mergers make
   `delaminated_m3` a leading craton sink on the new baseline; #395 and #403 attribute
   that account to fusions. Run 0005's broadened sink-ranking measurement alongside it:
   continental collision-placement loss can rival deep subduction
   ([#408](https://github.com/adubey/mantle-bloom/issues/408)).
2. 0002 after resolving attrition contact ownership
   ([#405](https://github.com/adubey/mantle-bloom/issues/405)).
3. 0009's gates. G1 (done, #426) found a mosaic inside today's plates doesn't give rifts a
   way around, so 0009 now generates cratons first and places plates around them. G3 (#472)
   initially selected RFC 0010's static minimum cut with a resistance fallback, but #480's
   routing diagnostics missed the per-seed target. RFC 0010 now proposes local damage-driven
   cuts that may link into a rift; its load and accumulation diagnosis gates must pass before
   implementation issues are opened. RFCs 0009 and 0010 replace 0003, whose cut steering
   couldn't avoid cratons (PR #418), and 0004.
4. 0006 and 0007, after separating their proposed production from existing elevation-driven
   continental resets ([#406](https://github.com/adubey/mantle-bloom/issues/406),
   [#407](https://github.com/adubey/mantle-bloom/issues/407)). They share `terranes.py`
   (docking, buoyant-crust tracers and ledger),
   which whichever lands first builds. RFC 0010's plume weakening also depends on 0007.

RFC 0008 came later, from the same craton work. Its gate can run at any time. Its fusion
lands with or after 0001's merge placement, since it makes craton-bearing merges more
frequent.

RFC 0009 keeps RFC 0004's deferral of a collisional formation path to a future RFC on
mantle depletion.

RFCs 0009 and 0010 replace RFC 0003, whose cut steering couldn't avoid cratons (PR #418),
and RFC 0004. An offline probe ([`analysis/rfc0010/`](../analysis/rfc0010/README.md))
found that even a cut free to bend rarely misses a craton, because each continent's
interior is one craton. RFC 0009 reverses the generation order: cratonic cores first,
continents assembled around them with lasting weak belts, then plates placed around them.
G3 initially selected RFC 0010's min-cut with a resistance fallback, but the #480 1 Gyr
same-state routing diagnostics showed the min-cut did not reach its per-seed 50% target.
RFC 0010 is revised to reject a static whole-plate minimum cut and propose incremental,
damage-driven local cuts that may link into rift segments; see
[`analysis/issue480/`](../analysis/issue480/README.md).

## Writing and reviewing

Use [`start-rfc`](../.agents/skills/start-rfc/SKILL.md) to begin or advance an RFC
workstream, including acceptance, diagnosis gates, implementation issues, and their PRs. Use
[`continue-rfc`](../.agents/skills/continue-rfc/SKILL.md) to reconstruct the state, evidence
and worktree after a handoff, then return to the RFC workflow.
Use [`project-management`](../.agents/skills/project-management/SKILL.md) for project
board status transitions during RFC and issue work. [`rfc-authoring`](../.agents/skills/rfc-authoring/SKILL.md)
owns RFC content and numbering and defines its sections (problem, solutions, design,
provisional touchpoints, assumptions and uncertainty, acceptance and falsification criteria,
references), citation format, and source-checking procedure.
RFC proposals are reviewed with
[`rfc-review`](../.agents/skills/rfc-review/SKILL.md), and completed or withdrawn with
[`complete-rfc`](../.agents/skills/complete-rfc/SKILL.md).
`.agents/skills/rfc-authoring/check_citations.py` checks citations and DOIs.

## Lifecycle

1. **Draft:** opened as a PR that adds the RFC; revised in place as review comes in.
2. **Accepted:** the project owner accepts an RFC by merging its RFC PR. If it has already
   been merged while still marked Draft, only the project owner's explicit confirmation
   counts. Before gates or implementation, verify acceptance from GitHub or an explicit
   owner statement. If neither is clear, ask the project owner; a generic request to start
   work is not acceptance. Any agent may ask, but only the owner decides. Set the status to
   Accepted and create or update the tracking issue immediately, listing diagnosis gates
   and planned implementation steps. Run each gate next and record its result in the RFC.
3. **Implemented:** implementation issues may be opened at acceptance when their scope is
   independent of diagnosis-gate results. If a gate determines the implementation approach
   or issue scope, wait for its result before opening those sub-issues. In either case,
   implementation does not begin until required gates pass. Each implementation issue has
   its own worktree and branch from `origin/main`; the RFC proposal worktree remains
   separate (see [`start-issue`](../.agents/skills/start-issue/SKILL.md)). An RFC with
   sub-issues may have multiple implementation PRs, one per issue. The status is updated
   when the last one lands and the acceptance criteria are met, with
   [`complete-rfc`](../.agents/skills/complete-rfc/SKILL.md). From then on the code and
   `docs/` are the source of truth, and the RFC stays as the record of why.
4. **Withdrawn:** the project stops pursuing it before it is implemented, because its
   falsification criteria were met, a later RFC supersedes it, or the project owner
   declines it. A draft or an accepted RFC can be withdrawn. It is kept, with "(withdrawn)"
   in its title and an Outcome section giving the reason, what landed and stays, and a link
   to any superseding RFC. Use [`complete-rfc`](../.agents/skills/complete-rfc/SKILL.md).
