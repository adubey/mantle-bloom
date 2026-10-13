---
name: rfc-review
description: Review a mantle-bloom RFC (rfc/NNNN-slug.md) for scientific grounding, design soundness and completeness, and post the review on its PR. Use when asked to review an RFC or an RFC pull request.
---

# Review an RFC

The aim is an RFC that can be accepted and then implemented without surprises. Check the
science against its sources, the design against the code, and the plan against what could
go wrong. Read [rfc-authoring](../rfc-authoring/SKILL.md) first: it defines the sections
and rules this review checks against.

## 1. Structure

Run `python3 .agents/skills/rfc-authoring/check_citations.py <rfc> --crossref`. Then check
that every section the authoring skill requires is present: problem, solutions, design,
provisional touchpoints, assumptions and uncertainty, acceptance criteria, falsification
criteria, references.

## 2. Evidence

- Can the problem's numbers be regenerated (seed plus a `bin/debug/` script)? Local files
  should be marked as not downloadable.
- Is a long run treated as a stress test, not as a calibration interval?
- Does the RFC claim more than the evidence shows? For example, attributing a lumped
  ledger account to one code path, or naming a cause before reproducing the event.
  Unconfirmed causes need a **diagnosis gate** that comes before fix selection, not
  alongside it.

## 3. Science

For each load-bearing citation, open the paper (abstract or full text) and check:
- the number quoted is in it;
- the RFC uses it within the paper's scope: a Cenozoic global budget, a classification by
  net behaviour over more than 10 Myr, a regional case study;
- a statistical tendency hasn't become a deterministic per-cell rule;
- gross production, net preservation and model-ledger fluxes are kept distinct;
- a proxy is not taken for the property itself (thickness for depletion, age alone for
  subduction polarity, mafic arc crust for continental crust). The RFC must state the
  assumption and how it will be validated.

## 4. Design

- **Units and scale:** quantities physical, and units worked through.
- **Bounds:** every factor on a physical field stays in range, including where factors
  compose (resistance × weakening × taper). Look for unbounded proxies feeding a `1 − k·s`
  form.
- **Timestep:** try a step longer than every duration in the design. Nothing should
  deliver more than its total, nonlinear rates should be integrated over the step, and
  random events should be keyed independently of the step partition.
- **Conservation vs. validity:** a design can conserve every ledger and still put material
  somewhere implausible. Ask for physical checks (distance limits, contiguity,
  thickness distributions) alongside ledger closure.
- **Accounting:** quantities summed are disjoint; each volume has one owner; every new
  tracer has a lifecycle through all writers of its parent field; ledger reservoirs are
  disjoint.
- **Scale checks:** redo the key arithmetic (volume over area, isostatic relief, rate ×
  band width × units). Do the stated inputs actually produce the stated outcome?
- **Side effects:** a change aimed at one sink can move others (plate counts, overlaps,
  split frequency). The acceptance criteria should watch for that.
- **Code:** trace the control flow the design relies on, and check that the touchpoints
  named exist (or are marked new). Check the step order between the writers the design
  touches; that save compatibility is handled; and that the runtime estimate covers the
  whole path, not just its cheapest part.
- **Knobs:** literature-constrained values are kept apart from calibration assumptions,
  and the latter have sensitivity runs.

## 5. Criteria

- Acceptance criteria are measurable, name their seeds and scripts, and put gates first.
- Falsification criteria name observations that would show the premise or mechanism
  wrong, and say what follows from each.

## 6. Post the review

Review the PR with `gh`: one inline comment per concern, anchored at the line it concerns,
and a summary comment for cross-cutting issues. Each comment states the problem, why it
matters, and a concrete change, with a link to the source when the science is at issue.
Mark the review "comment" or "request changes"; approval and merging are the project
owner's. On a later pass, say so if a previous comment was wrong, and withdraw it.
