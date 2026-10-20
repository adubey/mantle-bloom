---
name: project-management
description: Update or inspect RFC and issue statuses in Mantle Bloom GitHub Projects project 2. Workflow skills link here for status how-tos instead of duplicating board procedures.
---

# Manage the Mantle Bloom project board

The board is the user's GitHub Projects project 2. Keep issue and RFC work visible using
the board's existing single-select status field. This skill changes project item fields;
it does not close GitHub issues, merge PRs or change project fields/options.

## 1. Verify access and the board schema

- At the start of a session that will read or change project data, check `gh auth status`.
  The token must be valid
  and have the `project` scope. If it is invalid or lacks scope, stop project mutations and
  tell the user to reauthenticate or run `gh auth refresh -s project`; don't attempt a
  different account or silently skip the update.
- Inspect the project's fields and current items once per session, and reuse the discovered
  field and option IDs for later changes unless an error or ambiguity requires a refresh:

  ```bash
  gh project field-list 2 --owner @me --format json
  gh project item-list 2 --owner @me --limit 500 --format json
  ```

- Find the existing single-select field used for the board columns and read its option
  names. Do not assume its name or capitalization. Match the requested state to an
  existing option, preserving the exact option spelling. The intended states are
  `In Progress`, `Done`, `Ready`, `Deprecated`, `On Hold`, and `Backlog`. Use the exact
  existing spelling `Deprecated`; do not create a missing option. If the
  field or option is missing or ambiguous, report that and ask before changing the board
  schema.
- Check whether each issue or PR is already in the project. Add it only when absent; don't
  create duplicate project items.

## 2. Set a status

For an item already in project 2, set the existing status field using its exact name and
the exact option from `field-list`. The `--url`, `--field` and `--value` form is documented
in [`gh project item-edit`](https://cli.github.com/manual/gh_project_item-edit). Add an
absent issue or PR with
[`gh project item-add`](https://cli.github.com/manual/gh_project_item-add).

```bash
gh project item-edit 2 --owner @me --url <issue-or-pr-url> \
  --field "<status-field-name>" --value "<exact-option-name>"
```

If it is not in the project, add it first, then set its status:

```bash
gh project item-add 2 --owner @me --url <issue-or-pr-url>
gh project item-edit 2 --owner @me --url <issue-or-pr-url> \
  --field "<status-field-name>" --value "<exact-option-name>"
```

Verify the result with `gh project item-list 2 --owner @me --limit 500 --format json` and
confirm the item URL and status. Don't infer success only from a successful command exit.

## 3. Lifecycle how-to

### Starting an issue

When issue work begins, set its project item to `In Progress`. If a prerequisite prevents
work from starting, use `On Hold` and record the dependency in the issue.

### Starting or resuming an RFC

Set the RFC project item to `In Progress`. Use the RFC tracking issue once it exists;
before then, use the RFC PR. If both are project items, keep both `In Progress` while
RFC-level work continues. An active implementation sub-issue has its own item in
`In Progress`.

### Completing an issue

Set an issue to `Done` only after its acceptance criteria are met and its implementation is
merged or the project owner has otherwise confirmed completion. An open PR or local work
alone is not completion. If the issue belongs to an RFC, inspect its tracking checklist in
dependency order and tell the user the next step. If the next step is an eligible issue,
set it to `In Progress` only when work on it starts immediately; otherwise leave it `Ready`.
If blocked, set it to `On Hold`. If it is a gate, tell the user what gate to run; if no issue
exists for the next implementation step, say one needs to be created.

### Investigating additional work

Check whether the work belongs under this RFC or another one. If so, add it to the tracking
checklist and link the issue there. If no RFC covers it, create an issue for ordinary
issue-sized work; for a large, uncertain or science-led model change that meets
`AGENTS.md`'s RFC criteria, start an RFC first. An issue that blocks another issue puts the
blocked item `On Hold` until the dependency is resolved.

### Creating an issue

Use [`create-issue`](../create-issue/SKILL.md) to decide whether the issue is standalone
or an RFC sub-issue, prepare its scope and create it. Add every new issue to project 2 and
set it to `Backlog` by default. Use `Ready` only if the work is urgent and can start now.
If an eligible issue is being started immediately, follow [Starting an issue](#starting-an-issue)
to move it to `In Progress`; otherwise leave it `Ready` until work begins.

### Choosing not to pursue an issue

When an issue will not be worked on any further without a fix, set its project status to
`Deprecated` and close it with reason `not planned`. Preserve the issue and its discussion
as the record of the decision.

### Completing an RFC

After the completion PR merges or the project owner confirms the outcome, set the RFC PR
and tracking issue to `Done` when implemented or `Deprecated` when withdrawn, if both are
project items. Follow the [RFC completion workflow](../complete-rfc/SKILL.md#5-settle-the-tracking-issue-and-sub-issues)
to settle child issues, then update their project statuses to match their disposition.

## 4. Connect the other workflow skills

- `start-issue` follows [Starting an issue](#starting-an-issue) and
  [Completing an issue](#completing-an-issue).
- `start-rfc` follows [Starting or resuming an RFC](#starting-or-resuming-an-rfc) and uses
  [`create-issue`](../create-issue/SKILL.md) when follow-up work needs an issue.
- `continue-rfc` reconstructs handoff context, then routes lifecycle work to `start-rfc`.
- `complete-rfc` follows [Completing an RFC](#completing-an-rfc).
