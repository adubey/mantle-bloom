---
name: project-management
description: Keep Mantle Bloom RFCs and issues in GitHub Projects project 2 aligned with their work state. Use when starting, completing, blocking, deprioritizing or creating RFC and issue work.
---

# Manage the Mantle Bloom project board

The board is the user's GitHub Projects project 2. Keep issue and RFC work visible using
the board's existing single-select status field. This skill changes project item fields;
it does not close GitHub issues, merge PRs or change project fields/options.

## 1. Verify access and the board schema

- Check `gh auth status` before reading or changing project data. The token must be valid
  and have the `project` scope. If it is invalid or lacks scope, stop project mutations and
  tell the user to reauthenticate or run `gh auth refresh -s project`; don't attempt a
  different account or silently skip the update.
- Inspect the project's fields and current items:

  ```bash
  gh project field-list 2 --owner @me --format json
  gh project item-list 2 --owner @me --limit 500 --format json
  ```

- Find the existing single-select field used for the board columns and read its option
  names. Do not assume its name or capitalization. Match the requested state to an
  existing option, preserving the exact option spelling. The intended states are
  `In Progress`, `Done`, `Ready`, `Depreciated`, `On Hold`, and `Backlog`. In particular,
  do not silently change `Depreciated` to `Deprecated` or create a missing option. If the
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

## 3. Apply the lifecycle rules

- **Starting an issue:** move that issue to `In Progress` when work begins. If a required
  prerequisite prevents work from starting, set it to `On Hold` instead and record the
  dependency in the issue.
- **Starting or resuming an RFC:** move its project item to `In Progress`. Use the RFC
  tracking issue as the canonical item once it exists; before that, use the RFC PR. If
  both the RFC PR and tracking issue are already project items, keep both `In Progress`
  while RFC-level work continues. An active implementation issue also has its own item in
  `In Progress`.
- **Completing an issue:** mark it `Done` only after its acceptance criteria are met and
  its implementation is merged or the project owner has otherwise confirmed completion.
  An open PR or an agent finishing its local changes is not completion. If it belongs to
  an RFC, inspect the tracking checklist in dependency order and tell the user the next
  RFC step. If the next step has an issue and its prerequisites are met, move that issue
  to `In Progress`; if it is blocked, move it to `On Hold`. If the next step is a gate
  without an issue, tell the user what gate to run. If no issue exists for the next
  implementation step, say one needs to be created.
- **Investigation finds additional work:** check whether the work belongs under the current
  RFC or another RFC. If so, add it to that RFC's tracking checklist and link the issue
  there. If no RFC covers it and the change is ordinary issue-sized work, create an issue;
  if it meets `AGENTS.md`'s RFC criteria for a large, uncertain or science-led model
  change, start an RFC first. Add newly created issues to project 2. Set ordinary new work
  to `Backlog`; set it to `Ready` only when it is urgent and can start now. If it blocks
  another issue, put the blocked issue `On Hold` until the dependency is resolved.
- **Choosing not to pursue an issue:** set the project item to `Depreciated` when the
  decision is to stop without fixing the issue. Follow the relevant RFC or issue workflow
  for any GitHub issue closure; changing a project status does not close the issue.
- **A new issue:** add it to project 2 and set it to `Backlog` by default. The more specific
  rules above take precedence when the issue is urgent (`Ready`) or is the next eligible
  RFC step (`In Progress`).

## 4. Create an investigation issue

Before creating an issue, check the RFC tracking checklist and existing issues so work is
not duplicated. Use the repository's usual issue format and include the finding, evidence,
scope, acceptance checks and any dependency. When it is an RFC step, identify the RFC and
tracking issue; add it as a child of the tracking issue when appropriate. Create the GitHub
issue, using `--parent <tracking-issue-number>` when it is a child step:

```bash
gh issue create --title "<title>" --body-file <issue-body-file> \
  --parent <tracking-issue-number>
```

Omit `--parent` when the issue is not a child of an RFC tracking issue. Add the created
issue URL to project 2 with `item-add`, set its status according to section 3, and update
the RFC checklist when applicable. For an investigation that needs an issue, issue
creation is authorized by the project workflow; still follow the RFC requirement for
changes that need a new RFC.

## 5. Connect the other workflow skills

- `start-issue` calls this skill when issue work starts and when a completed issue advances
  the RFC checklist.
- `start-rfc` and `continue-rfc` call this skill when RFC work starts or resumes, and when
  they create a follow-up issue from an investigation.
- `complete-rfc` calls this skill to set the RFC's project item to `Done` when implemented,
  or `Depreciated` when withdrawn, and to settle any remaining active issue statuses.
