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
| [0004](0004-proto-cratons.md) | Proto-cratons stabilised by collision | Draft |
| [0005](0005-continental-recycling-balance.md) | Balance continental recycling at trenches | Draft |
| [0006](0006-intra-oceanic-arcs.md) | Intra-oceanic arcs: new continent from the ocean | Draft |
| [0007](0007-mantle-plumes-and-plateaus.md) | Mantle plumes: oceanic plateaus, flood basalts and hotspot tracks | Draft |

RFCs 0001-0007 come from one investigation: a 19.2 Gyr run (seed 456862899) that lost
almost all of its continental crust and all of its cratons. They are numbered by how
directly they act on that loss, not by priority. The suggested implementation order:
1. 0001;
2. 0002, together with 0005's measurement phase;
3. 0003, 0006 and 0004;
4. 0007, the largest. RFC 0003's plume weakening depends on it.

## Format

One Markdown file per RFC, `rfc/NNNN-short-slug.md`, numbered in order of creation. Each
has these sections, in this order:

1. **Header**: title, status (Draft, Accepted, Implemented, Rejected or Superseded by
   NNNN), date, and related issues and RFCs.
2. **Problem**: what goes wrong, in general terms, with the evidence. An RFC may cite local
   files such as saved worlds (`.mbworld`) or replay outputs, but must say they are not
   available for download, and give enough numbers that the argument stands without them.
   Prefer evidence anyone can regenerate: a seed and a script in `bin/debug/`.
3. **Proposed solutions**: one or more, each grounded in the geology, computational
   geodynamics, climatology, ecology or biology literature as the problem calls for.
   Cite as `[[Author et al., 2020]](#ref-author2020)`: the citation links to its entry in a
   **References** section at the end, and that entry links to the paper itself (a DOI
   link where one exists). One author is `[Author, 2020]`, two are
   `[Author & Other, 2020]`, three or more are `[Author et al., 2020]`. Say which solution
   is recommended and why.
4. **Design**: how the recommended solution works in the model: the quantities and their
   physical units (metres, m^3, Myr, km^2/Myr, ... rather than per-step or per-node
   factors), the runtime cost per step, save compatibility, and every knob that needs
   tuning, with its starting value, the source of that value, and how it will be tuned
   (which seeds, which script, which metric).
5. **Detailed design**: the actual source files and functions that change, the tests that
   go with them, and the docs that need updating.
6. **Acceptance criteria**: a checklist of measurable outcomes (replay metrics with
   thresholds, tests, ledger closure, runtime budget) that says when the RFC is done.

A **References** section closes the file. Keep the RFC current while it is a Draft. Once
it is Implemented, the code and `docs/` are the source of truth; the RFC stays as the record
of why.
