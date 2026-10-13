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
| [0002](0002-craton-trench-attrition.md) | Cratons are rarely subducted: margin-dependent attrition | Draft |
| [0003](0003-rifting-around-cratons.md) | Rifts go around cratons | Draft |
| [0004](0004-proto-cratons.md) | Proto-cratons: restarting craton formation | Draft |
| [0005](0005-continental-recycling-balance.md) | Balance continental recycling at trenches | Draft |
| [0006](0006-intra-oceanic-arcs.md) | Intra-oceanic arcs: new continent from the ocean | Draft |
| [0007](0007-mantle-plumes-and-plateaus.md) | Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks | Draft |

RFCs 0001–0007 come from one investigation: a 19.2 Gyr run (seed 456862899) that lost
almost all of its continental crust and all of its cratons. The run is used as a stress
test of the model's long-term behaviour, not as a calibration interval. The RFCs are
numbered by how directly they act on that loss, not by priority. The suggested
implementation order:
1. 0001;
2. 0002, together with 0005's measurement phase;
3. 0003 and 0004;
4. 0006 and 0007. They share `terranes.py` (docking, buoyant-crust tracers and ledger),
   which whichever lands first builds. RFC 0003's plume weakening also depends on 0007.

RFC 0004 defers a collisional formation path to a future RFC on mantle depletion.

## Writing and reviewing

RFCs are written with the [`rfc-authoring`](../.agents/skills/rfc-authoring/SKILL.md) skill.
It sets the sections (problem, solutions, design, provisional touchpoints, assumptions and
uncertainty, acceptance and falsification criteria, references), the citation format, and
how sources are checked. They are reviewed with
[`rfc-review`](../.agents/skills/rfc-review/SKILL.md).
`.agents/skills/rfc-authoring/check_citations.py` checks citations and DOIs.

## Lifecycle

1. **Draft:** opened as a PR that adds the RFC; revised in place as review comes in.
2. **Accepted:** the project owner merges it. Any diagnosis gate the RFC sets (a measurement that
   decides between fixes) is run next, and the RFC is updated with the result.
3. **Implemented:** the accepted design is split into issues and PRs that cite the RFC.
   The status is updated when the last one lands. From then on the code and `docs/` are
   the source of truth, and the RFC stays as the record of why.
4. **Rejected or Superseded:** kept, with the reason, or a link to what replaced it. This
   includes an RFC whose falsification criteria were met.
