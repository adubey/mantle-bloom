---
name: create-issue
description: Create a new Mantle-Bloom GitHub issue. Use for requests to create an issue, an RFC sub-issue, a standalone task, or an exploration tracker.
---

# Create a Mantle-Bloom issue

Use this when asked to create an issue or when an investigation identifies issue-sized
follow-up work, including an exploration tracking issue. This skill creates and links the
GitHub issue; use
[`project-management`](../project-management/SKILL.md#creating-an-issue) for project-board
membership and its initial status.

## 1. Decide where the work belongs

- Read `AGENTS.md`, the relevant issue comments, and the RFC index or RFC when applicable.
- Check existing issues and RFC tracking checklists for duplicate or already-planned work.
- An exploration tracker, when requested, is a new issue linked to its source item; do not
  reuse or reopen the issue or RFC that prompted the exploration. Under
  [`edit-rfc`](../edit-rfc/SKILL.md), the project owner's invocation approves the tracker.
- An exploration tracker requested by `explore-problem` is a bounded investigation issue,
  even when the underlying scientific question is large or uncertain. It does not authorize
  model implementation or replace an RFC that may be needed for a resulting model change.
- For an RFC implementation sub-issue, follow the acceptance, gate and issue-timing rules
  in the [canonical RFC lifecycle](../../../rfc/README.md#lifecycle). If the work changes
  the RFC's design or plan, revise the RFC first with
  [`rfc-authoring`](../rfc-authoring/SKILL.md). Create an eligible issue as a child of that
  RFC's tracking issue and update its checklist in dependency order.
- If no RFC covers proposed implementation work, create a standalone issue only for small,
  bounded changes. For a large, uncertain or science-led model change that meets
  `AGENTS.md`'s RFC criteria, start an RFC instead of treating an issue as authorization to
  implement it.

## 2. Write and create the issue

Use the repository's usual issue format. State the problem or finding, relevant evidence,
bounded deliverable, acceptance checks, and dependencies. For an exploration tracking
issue, make the investigation and its findings summary the deliverable; it is not an
implementation issue and does not require an RFC merely because the scientific question is
large or uncertain. Follow [`explore-problem`](../explore-problem/SKILL.md) for its progress
comments and summary. For an RFC sub-issue, cite the RFC and tracking issue, and keep the
deliverable consistent with the accepted design.
Avoid including personal filesystem paths or machine-specific artifact locations; link to
durable repository evidence or use neutral run labels.

Create the GitHub issue, supplying `--parent` only for an RFC sub-issue:

```bash
gh issue create --title "<title>" --body-file <issue-body-file> \
  --parent <tracking-issue-number>
```

After creation, add its URL to the RFC tracking checklist when applicable. Then follow
[Creating an issue](../project-management/SKILL.md#creating-an-issue) to add it to project
2 and set its initial status. If it is being started immediately, follow
[Starting an issue](../project-management/SKILL.md#starting-an-issue) for its status.
