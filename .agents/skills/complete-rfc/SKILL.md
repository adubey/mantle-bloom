---
name: complete-rfc
description: Complete or withdraw a mantle-bloom RFC (rfc/NNNN-slug.md). Close it as Implemented when its acceptance criteria are met, or as Withdrawn when it is falsified, superseded or dropped, and settle its landed code, open issues, dependent RFCs and tracking issue. Use when asked to finish, close out, withdraw, abandon, supersede or reject an RFC.
---

# Complete or withdraw an RFC

Every RFC ends in one of two states. Either way it stays in `rfc/` as the record of why.

- **Implemented:** its accepted design has landed and its acceptance criteria are met. From
  then on the code and `docs/` are the source of truth.
- **Withdrawn:** the project stops pursuing it before it is implemented. The reason is one
  of these:
  - its falsification criteria were met;
  - a later RFC supersedes it;
  - review or the project owner declined it.

  This covers drafts and accepted RFCs alike, including ones with work already landed.

Withdrawing is the project owner's decision. An agent that finds a reason to withdraw (a
falsification criterion met, say) reports the evidence and recommends withdrawal, but
doesn't withdraw the RFC until the owner says so.

## 1. Gather the record

- Read the RFC, its tracking issue with comments (`gh issue view <N> --comments`), and
  every issue and PR the tracking issue or the RFC's **Related** line names. Note which
  merged, which are open, and which are still being worked on (open PRs, recent commits on
  their branches, other worktrees: `git worktree list`).
- Note every gate result and replay measurement, with its commit and `analysis/`
  directory.
- Find the RFC's dependents: `grep -rn "NNNN" rfc/ docs/ AGENTS.md .agents/`. Other RFCs
  may depend on its design, for example a hook it promised to provide.

## 2. Check the outcome

- **Implemented:** every acceptance box can be ticked against a named measurement or test,
  on the seeds and baseline the RFC names. If one can't be ticked, the RFC isn't
  implemented yet. Either finish the work, or ask the owner to waive that criterion and
  record the waiver with its reason.
- **Withdrawn:** state the reason in one or two sentences, with its evidence: the
  falsification criterion met and the measurement that met it, the superseding RFC, or the
  owner's decision. If a replacement RFC is planned, it should exist at least as a draft
  PR, so the withdrawal can link to it.

## 3. Update the RFC

- **Title:** a withdrawn RFC's title gets the outcome, so it shows wherever the title is
  quoted: `# RFC 0003 (withdrawn): Rifts go around cratons`. An implemented RFC's title
  doesn't change.
- **Status line:** `Implemented (YYYY-MM-DD)`, or `Withdrawn (YYYY-MM-DD)`, followed by
  `; superseded by [RFC NNNN](NNNN-slug.md)` when one exists. Add the completion date to
  the **Date** line too.
- **Outcome section:** add `## Outcome` directly after the header, before **Problem**, so a
  reader learns the outcome before reading the proposal. It covers:
  - what happened, and why, with evidence links (`analysis/` directories, PRs, commits);
  - what landed and stays, and where its behaviour is now documented;
  - for a withdrawal, what carries over to the superseding RFC (evidence, gates still
    valid, design pieces) and what doesn't;
  - questions left open, and where they went.
- **Leave the body as written.** It is the record of what was proposed and why. Don't
  rewrite the problem or design to match the outcome. Tick an acceptance box only with
  evidence; on a withdrawn RFC, leave unmet boxes unticked.

## 4. Settle landed code

For each merged implementation PR, decide whether its change stays.

- **Implemented:** it stays. Check that `docs/simulation-model.md` (or the relevant doc)
  describes the behaviour without depending on the RFC.
- **Withdrawn:** keep a landed change if it stands on its own (a consistency fix, a
  diagnostic account), and say so in the Outcome. If a superseding RFC builds on it, say
  that RFC now owns it. If it only makes sense as part of the withdrawn design, open an
  issue to remove it rather than reverting it in the completion PR. Removal is a model
  change, with its own tests and replays.
- Code comments citing the RFC stay: the RFC remains the record of why the code is there.

## 5. Settle the tracking issue and sub-issues

- Inspect the tracking issue's GitHub sub-issues as well as its checklist. Every checklist
  step and every child issue must have a disposition before the RFC is complete.
- **Implemented:** verify each child issue's work is merged or otherwise confirmed complete.
  If a child remains open because its work is still needed, the RFC is not complete. Close
  children whose work has landed with reason `completed`; close any made obsolete by the
  outcome with reason `not planned` and a comment linking the RFC's Outcome.
- **Withdrawn:** close the tracking issue with reason `not planned` and a comment linking
  the RFC's Outcome. Keep open any child issue that is no longer pursued and set its project
  status to `Depreciated`; do not close it as not planned. Move relevant children to the superseding RFC's
  tracking issue and plan. If still relevant but no RFC covers them, remove them as
  sub-issues of the withdrawn tracking issue and keep them as standalone issues. For an
  existing issue, use `gh issue edit <N> --parent <new-tracking-issue>` to move it or
  `gh issue edit <N> --remove-parent` to make it standalone, then verify with
  `gh issue view <N> --json parent`. Do not close work that remains needed; keep its
  project status aligned with its disposition.
- For either outcome, settle the child issues when the completion PR merges, or earlier if
  the project owner has already confirmed the outcome. Do not close or push to another
  agent's or the owner's PR without the owner's go-ahead.
- **Open PRs:** comment with the outcome and what you recommend (close, or salvage its
  analysis). Don't close or push to another agent's or the owner's PR without the owner's
  go-ahead.
- **Tracking issue:** for an implemented RFC, the completion PR closes it (`Closes #N`).
  For a withdrawn RFC, close it with `gh issue close <N> --reason "not planned"` after the
  completion PR merges; don't use `Closes #N`, which would mark it completed. Tick or
  strike checklist items first, with the reason for each strike.
- After the completion PR merges or the owner confirms the outcome, follow
  [Completing an RFC](../project-management/SKILL.md#completing-an-rfc) to set both the RFC
  PR and tracking issue's project status, and reconcile child issue statuses.

## 6. Update dependents and the index

- `rfc/README.md`:
  - set the index row's status: `Implemented`, or `Withdrawn`, plus `(superseded by
    NNNN)` when one exists;
  - update the implementation-order prose that names the RFC.
- Other RFCs that depend on it: repoint their **Related** lines and the text that relies on
  its design to the superseding RFC. Leave any section that is historical evidence as it
  is. Don't change another RFC's design beyond the repointing; open a follow-up issue if
  it needs more.
- `docs/` that cite the RFC for behaviour that stays need no change.

## 7. Open the PR

- One PR per completion, branched off `origin/main` (`<agent>/rfc-NNNN-complete` or
  `<agent>/rfc-NNNN-withdraw`). Several RFCs withdrawn for the same reason can share one.
- `## Summary`: the outcome and its reason, the evidence, what stays, what carries over,
  and what was closed or moved. End with `Closes #<tracking issue>` for an implemented RFC;
  for a withdrawn RFC, link the tracking issue without an auto-close directive so it can be
  closed with reason `not planned` after merge.
- Run `python3 .agents/skills/rfc-authoring/check_citations.py rfc/NNNN-slug.md` on every
  RFC you edited.
