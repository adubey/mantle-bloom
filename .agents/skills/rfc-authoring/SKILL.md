---
name: rfc-authoring
description: Write or revise the scientific content of a Mantle-Bloom RFC. Use when asked to author a proposal, turn investigation findings into an RFC, or address technical review comments; this skill also chooses the NNNN for a new RFC. To split off a revision into a new session, use edit-rfc.
---

# Write an RFC

An RFC proposes; it doesn't commit to code. It records the evidence, alternatives,
assumptions and uncertainty, rough costs, and *provisional* touchpoints in the source. This
skill owns the document: for a new RFC, choose its next `NNNN` by checking the index and
existing filenames, then write and validate the content. Use `start-rfc` to set up or
advance the RFC workstream, including acceptance, gate sequencing, and implementation
planning. Open the RFC proposal PR from this skill; implementation work is proposed through
one or more implementation PRs coordinated by `start-rfc`, generally one per sub-issue.
Follow the [RFC lifecycle](../../../rfc/README.md#lifecycle) for issue timing and gates.

## 1. Gather the evidence

- Prefer evidence anyone can regenerate: a seed plus a script in `bin/debug/` (with a usage
  docstring). `bin/debug/ledger_history.py` turns any save's `stats_history` into craton
  and continental-material budgets.
- A local `.mbworld` save or replay output may be cited by file name, but say it is **not
  available for download**, and quote enough numbers that the argument stands without it.
- **A long run is a stress test of model behaviour, not a calibration interval.** Runs of
  hundreds of Myr are typical; a multi-Gyr run is rare. Long runs expose ratchets and
  missing sources and sinks. Don't fit Earth rates to them, and don't read their late,
  degenerate states (no continents, no cratons) as geology.
- Say what the evidence does *not* show. A ledger account that lumps several code paths
  can't attribute a loss to one of them.
- **Read the code before hypothesising about it.** Trace the control flow of every
  function a hypothesis or fix depends on: what it already searches, clips or books, and
  which writers run between it and where the effect is measured. A fix that repeats
  something the code already does can't help.
- **Reproduction needs every input.** A replay or captured call needs everything the real
  call used (third-party geometry, not just positions), plus the original step schedule
  and controls. A seed alone doesn't reproduce a particular event, because topology and
  per-step random streams depend on step length.

## 2. Find and check the literature

- Ground each proposed mechanism in whatever literature the problem calls for: geology,
  geodynamics, climatology, hydrology, ecology, biology, numerical methods, or another
  field.
- **Verify every reference:** run
  `python3 .agents/skills/rfc-authoring/check_citations.py --crossref`. It checks that
  every DOI resolves and prints title, authors and year to compare against the citation.
- **Verify every number you quote in the paper's own text** (the abstract, or the full
  text when the abstract doesn't have it), not from a search summary or memory. If you
  can't confirm a figure, don't use it, or label it unverified. Check the paper's
  *mechanism* the same way. A paper on cratons formed in extension doesn't support
  cratons formed in compression.
- Respect each source's scope.
  - Published fluxes are usually process-specific, global, and focused on the Cenozoic or
    Phanerozoic. Use them as **sanity checks** with that scope stated, not as targets for
    a model ledger.
  - Statistical tendencies ("erosion is favoured where…") are not per-cell rules. Keep
    them probabilistic, or averaged in time and space, unless you justify a hard
    threshold.
  - Keep **gross production** (what a process makes), **net preservation** (what
    survives) and **model-ledger fluxes** (what an account books) distinct. Never compare
    one kind against another.
  - Match definitions before comparing numbers: a total thickness isn't an excess over the
    model's reference, and a porous sediment column isn't a rock-equivalent volume.
    Convert explicitly.

## 3. Structure

`rfc/NNNN-short-slug.md`, numbered in order of creation, in Markdown, with these sections
in this order:

1. **Header:** title; status (Draft, Accepted, Implemented, or Withdrawn with any superseding RFC; see `complete-rfc`);
   date; related issues and RFCs. Revising an accepted RFC doesn't change its Accepted status
   or add an "under revision" state. Add each revision's date and trigger to the status line,
   newest first, as in `Accepted (revised 2026-10-30 after G2; revised 2026-10-17 after G1)`.
2. **Problem:** what goes wrong, in general terms, with the evidence (section 1 above).
3. **Proposed solutions:** what happens in nature, with citations; then one or more
   alternatives, each with its trade-offs; then the recommendation and why.
4. **Design:** how the recommended solution works.
   - Quantities in **physical units** (m, m³, Myr, km²/Myr, m³/yr), not per-step or
     per-node factors. Work the units through, and do a scale check against the
     literature where one exists.
   - **Bounds:** every factor applied to a physical field stays in its valid range (a
     fraction in [0, 1]; a thickness ≥ 0), and composed factors are clamped.
   - **Timestep robustness:** results don't depend on the step length.
     - An event shorter than a step integrates over the overlap of its interval with the
       step, and never delivers more than its total.
     - A rate that is nonlinear in a quantity changing during the step is integrated
       along that quantity's path, not evaluated at the end and multiplied by the step.
     - Random events are keyed by event index or continuous time, not by step, so the same
       seed gives the same events under any step partition.
     - Test with steps longer than every duration in the design, comparing accumulated
       results, not just final states.
   - **Conservation and ownership:**
     - Sum only disjoint quantities: check which fields are shares of which (for example,
       mobile cover inside Hc) before adding them.
     - Give each physical volume one owner, so two code paths don't both add it.
     - A new tracer that is a share of another field needs its whole lifecycle: every
       writer of the parent field (erosion, thinning, resets, caps, transfers, merges)
       must update it or be caught by an audit. Registering it as EXTENSIVE covers only
       remaps.
     - A ledger's reservoirs must be disjoint; count repeatable events outside the
       balance.
   - **Diagnosis gates:** if the cause isn't confirmed, a measurement phase comes first. Say
     what will be measured and which outcome leads to which fix. Fixes that depend on it
     are chosen only after the gate, and the RFC is updated with the result.
   - **Rough runtime cost** per step, and **save compatibility** (new pickled state, and
     how old saves load).
   - **Knobs:** each with its starting value, the source of that value (a citation or
     "design choice"), and how it will be tuned: which seeds, which script, which metric.
     Separate literature-constrained values from calibration assumptions, and plan
     sensitivity runs for the latter. Keep tuning and sensitivity seeds separate from
     held-out acceptance seeds; do not tune against held-out results. Only the project owner
     may change an agreed diagnosis gate or acceptance criterion. If evidence suggests a
     criterion is invalid, propose a revision in the RFC with the old criterion, proposed
     change and reason, and get the owner's decision before using it for a rerun. Under
     [`edit-rfc`](../edit-rfc/SKILL.md), the owner's merge of a PR that sets this out is that
     decision.
5. **Detailed design (provisional touchpoints):** the source files and functions likely
   to change, plus their tests and docs. Check that every file, function and constant you
   name exists, and mark new ones as new.
6. **Assumptions and uncertainty:** what the design takes as given; which numbers are
   weakly constrained; which proxies stand in for physics the model lacks (a thickness
   standing in for depletion, for example), and how each will be validated.
7. **Acceptance criteria:** a checklist of measurable outcomes: replay metrics with
   thresholds on named seeds, unit tests, ledger closure, runtime and save-size budgets,
   old saves loading. Put the diagnosis gates first.
8. **Falsification criteria:** observations that would show the RFC's premise or
   mechanism is wrong, and what follows from each (reject, revise, or switch to an
   alternative). They differ from acceptance criteria: acceptance says the change works,
   falsification says the idea was mistaken. Make them checkable with the same scripts.
9. **References.**

Two more sections go between the header and **Problem** when they apply, outside the numbered
list. `## Outcome` comes first and records how an RFC ended; see `complete-rfc`.
`## Revision note` comes next, for an accepted RFC revised through
[`edit-rfc`](../edit-rfc/SKILL.md). Add one dated entry per revision, newest first. Each entry
gives the revision PR, what the earlier version proposed, what the evidence showed with links,
and what changed, as in RFC 0009. For each changed gate or acceptance criterion, it gives the
old criterion, new criterion and reason, and it links the impact table. A draft RFC doesn't
need a note. Revise the body in place.

## 4. Citations

- In text: `[[Author et al., 2020]](#ref-author2020)`. One author is `[Author, 2020]`, two
  are `[Author & Other, 2020]`, three or more `[Author et al., 2020]`. Add a/b suffixes
  for same-author, same-year pairs.
- Each reference entry starts with `<a id="ref-author2020"></a>` and ends with a DOI link
  (`https://doi.org/...`), or the publisher's URL when there is no DOI.
- `check_citations.py` (no flags) checks that every in-text anchor has an entry and every
  entry is cited.

## 5. Finish

- Add the RFC to the index in `rfc/README.md`.
- Run `check_citations.py --crossref` and `backend/.venv/bin/ruff check` (for any new
  `bin/debug/` script).
- Open the **RFC proposal PR** that adds the RFC, with a `## Summary` of the problem and
  recommendation. This is separate from implementation PRs, which `start-rfc` coordinates
  according to the [RFC lifecycle](../../../rfc/README.md#lifecycle).
  When addressing review comments, revise the RFC in place and add a `## Review fixes`
  section to the PR body listing each comment and what changed.
- When RFC proposal work starts or resumes, follow
  [Starting or resuming an RFC](../project-management/SKILL.md#starting-or-resuming-an-rfc)
  to update its project status.
