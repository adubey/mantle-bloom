# RFCs

An RFC proposes a change to the simulation model that is too large, too uncertain or too
tied to the science to start as a plain issue. It explains the problem, compares solutions
against the literature, works out the design to the level of units, knobs and source files,
and ends with criteria that say when the change is done. Once accepted, an RFC is
implemented through ordinary issues and PRs (see [AGENTS.md](../AGENTS.md)), each citing it.

## Index

| RFC | Title | Status |
|---|---|---|
| [0001](0001-merge-cap-craton-loss.md) | Plate merges must not destroy cratons | Draft |
| [0002](0002-craton-trench-attrition.md) | Cratons resist sinking: margin-dependent attrition and collision sinks | Draft |
| [0003](0003-rifting-around-cratons.md) | Rifts go around cratons | Accepted |
| [0004](0004-proto-cratons.md) | Proto-cratons: restarting craton formation | Draft |
| [0005](0005-continental-recycling-balance.md) | Balance continental recycling at trenches | Draft |
| [0006](0006-intra-oceanic-arcs.md) | Intra-oceanic arcs: new continent from the ocean | Draft |
| [0007](0007-mantle-plumes-and-plateaus.md) | Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks | Draft |
| [0008](0008-craton-collision-jamming.md) | Cratons jam collisions: braking and fusion | Accepted |
| [0009](0009-craton-mosaic.md) | Cratons as a mosaic: discrete cores and lasting weak belts | Draft |
| [0010](0010-rifts-follow-weak-lithosphere.md) | Rifts follow weak lithosphere | Draft |

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
   ([#405](https://github.com/adubey/mantle-bloom/issues/405)), and 0003's selected
   combination of cut steering and stronger resistance. RFC 0003's phase 1b is already
   on the year-based pace (#388, PR #392); #389 tracks its implementation.
3. 0004 only after the sink repairs, and only once the baseline passes its evaluation gate:
   a qualifying window after 500 Myr on at least two runs. `2562f63` has one. Most of its
   formation is a 200 Myr generation pulse, and its share rises through loss of
   non-cratonic crust, not formation
   ([#409](https://github.com/adubey/mantle-bloom/issues/409),
   [`analysis/issue409/`](../analysis/issue409/README.md)).
4. 0006 and 0007, after separating their proposed production from existing elevation-driven
   continental resets ([#406](https://github.com/adubey/mantle-bloom/issues/406),
   [#407](https://github.com/adubey/mantle-bloom/issues/407)). They share `terranes.py`
   (docking, buoyant-crust tracers and ledger),
   which whichever lands first builds. RFC 0003's plume weakening also depends on 0007.

RFC 0008 came later, from the same craton work. Its gate can run at any time. Its fusion
lands with or after 0001's merge placement, since it makes craton-bearing merges more
frequent.

RFC 0004 defers a collisional formation path to a future RFC on mantle depletion.

RFCs 0009 and 0010 replace RFC 0003, whose cut steering couldn't avoid cratons (PR #418),
and RFC 0004. An offline probe ([`analysis/rfc0010/`](../analysis/rfc0010/README.md))
found that even a cut free to bend rarely misses a craton, because each continent's
interior is one craton. 0009 makes cratons a mosaic of cores and lasting weak belts; its
gate G1 also decides whether 0010's strength-following cut is built.

## Writing and reviewing

Use [`start-rfc`](../.agents/skills/start-rfc/SKILL.md) to set up a new RFC workstream and
plan its transition into one or more implementation issues. RFCs are written with the
[`rfc-authoring`](../.agents/skills/rfc-authoring/SKILL.md) skill.
It sets the sections (problem, solutions, design, provisional touchpoints, assumptions and
uncertainty, acceptance and falsification criteria, references), the citation format, and
how sources are checked. They are reviewed with
[`rfc-review`](../.agents/skills/rfc-review/SKILL.md).
`.agents/skills/rfc-authoring/check_citations.py` checks citations and DOIs.

## Lifecycle

1. **Draft:** opened as a PR that adds the RFC; revised in place as review comes in.
2. **Accepted:** the project owner merges it, or, for an RFC already merged as a draft,
   says they accept it; the status is then set to Accepted. Any diagnosis gate the RFC sets (a measurement that
   decides between fixes) is run next, and the RFC is updated with the result.
3. **Implemented:** the accepted design is split into issues and PRs that cite the RFC,
   listed on the RFC's tracking issue, opened when it is accepted. The
   issues use the RFC worktree directory and keep one issue branch and PR per issue (see
   [`start-rfc`](../.agents/skills/start-rfc/SKILL.md) and
   [`start-issue`](../.agents/skills/start-issue/SKILL.md)). The status is updated when the
   last one lands. From then on the code and `docs/` are the source of truth, and the RFC
   stays as the record of why.
4. **Rejected or Superseded:** kept, with the reason, or a link to what replaced it. This
   includes an RFC whose falsification criteria were met.
