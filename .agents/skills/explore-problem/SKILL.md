---
name: explore-problem
description: Investigate a Mantle Bloom problem interactively, keep findings in an exploration issue, and recommend scientific hypotheses and next artifacts.
---

# Explore a Mantle Bloom problem

Use this before authoring an RFC or implementation issue when the problem, cause, scope,
or best next step is still uncertain. Also use it when an RFC or issue was abandoned
because a gate failed, acceptance criteria were not met, or the proposed work appeared
infeasible, and the user wants to understand what happened or whether another path exists.
The exploration is collaborative: involve the user as evidence suggests useful directions,
and ask which artifacts to create before making them.

## 1. Frame the problem with the user

- Clarify the observed behavior, who or what is affected, and what decision the exploration
  should inform. Read the relevant issue, RFC, PR and discussion when provided. Start by
  summarizing the problem and main unknowns; gather initial evidence before proposing
  avenues so early guesses do not anchor the investigation.
- For simulation or model questions, read `AGENTS.md`, `docs/architecture.md`,
  `docs/simulation-model.md`, `docs/design-principles.md`, and the relevant topic page.
  Read `docs/debugging.md` for saved-world diagnostics.
- Assess the relevant design or modeling principles as criteria: identify whether the
  behavior appears to violate one, is consistent with one, or has insufficient evidence.
  Separately formulate causal hypotheses that could explain the observations.
- Create a new exploration tracking issue early, after the initial framing, using
  [`create-issue`](../create-issue/SKILL.md); never reuse or reopen the issue or RFC that
  prompted the exploration. Link related RFCs or issues; make the new issue a child only
  when the investigation is an authorized step under an active accepted RFC. In its
  description, state the question and include a `## Exploration summary` section with this
  note: “Intermediate comments record selected findings and may be disregarded. Read the
  Exploration summary for the consolidated results.” Start the summary as pending. Let
  [`create-issue`](../create-issue/SKILL.md) add it to project 2; because exploration starts
  immediately, follow
  [`project-management`](../project-management/SKILL.md#starting-an-issue) to set it to
  `In Progress`. Follow the same guide for other transitions; ask the user only when the
  requested status or appropriate transition is unclear.
- After gathering initial evidence, present a few evidence-based avenues, alternatives and
  a recommendation. Ask a focused question when the user's knowledge or preference could
  materially change which avenue is useful.

## 2. Investigate interactively

- Trace current behavior through code and durable project records before proposing a fix.
  Check related issues, RFCs, analysis reports, and existing runs to avoid repeating work.
- When exploring an abandoned or failed RFC or issue, verify its current GitHub state and
  inspect the stated reason, relevant acceptance criteria, gate definition, PR discussion,
  and available evidence. Treat failure as evidence about the tested design and conditions;
  distinguish that from evidence that the whole problem is infeasible. Check whether the
  result could instead reflect an invalid assumption, input or setup mismatch, measurement
  problem, incomplete evidence, or criteria that need reconsideration. Do not rerun an
  expensive failed gate until existing evidence and a meaningful change to the run are
  understood. The exploration goal is to find plausible, interesting hypotheses about the
  failure and possible mechanisms, not to test a candidate implementation.
- Share concise findings as they emerge. When evidence opens materially different avenues,
  explain the evidence, alternatives and a recommendation, then ask the user which direction
  to explore. Continue independent investigation while waiting for input; don't interrupt
  for routine choices.
- The exploration phase finds and sharpens hypotheses. Inspect code, literature, existing
  runs, and observations, and use non-invasive diagnostics of the current model when they
  help distinguish explanations. For each promising hypothesis, state its proposed
  mechanism, which observations it explains, a discriminating prediction, and important
  uncertainties or alternative explanations. Do not prototype or test a candidate model
  implementation here; that belongs in an RFC investigation with a selected implementation
  and explicit diagnosis gates.
- Distinguish observations and exploratory evidence from implementation verification and
  model validation. For scientific evidence, account for controls, initial conditions,
  units, seeds, resolution and relevant uncertainty; use representative runs when stochastic
  or parameter-sensitive behavior matters. Follow `replay-ab` for its controlled before/after
  comparisons and avoid repeating matching prior runs.
- Use `git worktree list --porcelain` to locate relevant registered worktrees and resolve
  repository-relative files beneath each listed `worktree` root. Treat other worktrees as
  read-only evidence unless the user explicitly directs otherwise. The listing does not
  identify an owner or expose worktrees from another clone or host.
- If relevant evidence or uncommitted work is in an inaccessible worktree, clone, or host,
  record it as a missing dependency, continue independent investigation, and ask the user
  for access or a durable handoff only when that evidence is needed to decide the next step.
- Do not implement a model change or create outcome artifacts during exploration. The
  exploration tracking issue is the exception: add comments containing concise, interesting
  findings, evidence links, and changed hypotheses—not a transcript of the conversation.
  Keep its `## Exploration summary` pending until the evidence supports a useful handoff,
  then summarize findings, hypotheses, uncertainties, and recommended next steps. Update
  that summary when conclusions materially change or before ending the exploration. Keep
  the exploration issue as the investigation record; when the user selects a downstream
  RFC or issue, link this summary from it and carry over only findings it needs rather than
  duplicating the full exploration.
- Put reusable evidence in `analysis/issue<N>/` in line with `AGENTS.md`; keep only compact
  JSON and reports that belong in the repository. For temporary evidence used mainly to
  develop or test an RFC, record how to regenerate it in the exploration issue and, once
  selected, in the RFC's diagnosis-gate record. Include the command or script, commit,
  neutral input/save label and checksum, seed, schedule, controls, and measured values as
  applicable. Raw saves and large run outputs may remain ephemeral and can be removed when
  the RFC is completed or withdrawn, once no active follow-up depends on them. Never remove
  another task's files.
- Never publish personal usernames, home-directory paths, absolute machine paths, private
  mount paths, saves, or raw run output. If same-host recovery needs a local path, use an
  ignored local-only handoff note and verify it is ignored.

## 3. Present outcomes and ask what to create

When the evidence is sufficient or the user wants to decide, provide a concise summary of
the problem, evidence, applicable principles, remaining uncertainty, and a prioritized
list of next steps. Outcomes may be combined; do not force the investigation into one
category:

1. One or more new RFCs for large, uncertain, or science-led model changes.
2. One or more issues for small, bounded changes, either standalone or as sub-issues of an
   accepted RFC after its diagnosis gates pass.
3. A revision or reassessment of one or more existing RFCs.
4. No viable path identified within the explored scope, with the reason and evidence. If an
   issue is abandoned for this reason, keep it open and set its project status to
   `Depreciated`. If an RFC or issue was already abandoned, preserve its history; propose
   separate follow-up work if the investigation identifies a viable new scope.

For each proposed next step, name the artifact, its purpose, dependencies and the status
prescribed by `project-management`. Ask which outcome artifacts the user wants created or
revised. Apply established status transitions from that skill; ask for input only when a
status is not covered or the user's preference is unclear. Do not infer approval to create
or revise outcome artifacts from approval to investigate. Do not close an issue or withdraw
an RFC without the user's direction.

## 4. Carry out the user's choices

After the user selects artifacts:

- For a new RFC, follow [`rfc-authoring`](../rfc-authoring/SKILL.md) and
  [`start-rfc`](../start-rfc/SKILL.md).
- For a standalone issue or an RFC sub-issue, follow
  [`create-issue`](../create-issue/SKILL.md). Respect RFC acceptance and diagnosis gates.
- For an existing RFC revision, follow [`rfc-authoring`](../rfc-authoring/SKILL.md); use
  [`rfc-review`](../rfc-review/SKILL.md) when the selected work is a review.
- For an RFC the user chooses to withdraw, follow [`complete-rfc`](../complete-rfc/SKILL.md).
  For an issue the user chooses to deprioritize or abandon, follow
  [Choosing not to pursue an issue](../project-management/SKILL.md#choosing-not-to-pursue-an-issue);
  keep the issue open with status `Depreciated`; a project status change does not close it.
- When the user confirms the exploration summary is complete, follow
  [Completing an issue](../project-management/SKILL.md#completing-an-issue) for the
  exploration tracking issue. Its deliverable is the findings summary; implementation does
  not need to be merged for the user to confirm the investigation is complete.
- Follow the relevant how-to in [`project-management`](../project-management/SKILL.md) and
  verify each requested board update.

Do not start implementation merely because an issue or RFC was created. Use
[`start-issue`](../start-issue/SKILL.md) only when the user also wants implementation to
begin.
