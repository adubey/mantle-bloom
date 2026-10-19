---
name: create-issue
description: Create a Mantle Bloom GitHub issue as a standalone task or an additional step under an RFC tracking issue.
---

# Create a Mantle Bloom issue

Use this when asked to create an issue or when an investigation identifies issue-sized
follow-up work. This skill creates and links the GitHub issue; use
[`project-management`](../project-management/SKILL.md#creating-an-issue) for project-board
membership and its initial status.

## 1. Decide where the work belongs

- Read `AGENTS.md`, the relevant issue comments, and the RFC index or RFC when applicable.
- Check existing issues and RFC tracking checklists for duplicate or already-planned work.
- If this is an additional step for an RFC, confirm the RFC is accepted and its required
  diagnosis gates have passed. If the work changes the RFC's design or plan, revise the
  RFC first with [`rfc-authoring`](../rfc-authoring/SKILL.md). Create the issue as a child
  of that RFC's tracking issue and update its checklist in dependency order.
- If no RFC covers the work, create a standalone issue. For a large, uncertain or
  science-led model change that meets `AGENTS.md`'s RFC criteria, start an RFC instead of
  treating the issue as authorization to implement it.
- If the RFC is still a draft or its gates are pending, do not create an implementation
  sub-issue. Record the needed proposal or gate work in the RFC workflow first.

## 2. Write and create the issue

Use the repository's usual issue format. State the problem or finding, relevant evidence,
bounded deliverable, acceptance checks, and dependencies. For an RFC sub-issue, cite the
RFC and tracking issue, and keep the deliverable consistent with the accepted design.
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
