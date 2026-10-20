---
name: explore-problem
description: Explore a Mantle Bloom problem interactively to find scientific hypotheses. Use before an RFC or small issue, or after an RFC or issue fails a gate or appears infeasible.
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
- After gathering initial evidence, present a few evidence-based avenues, alternatives and
  a recommendation. Ask which direction to investigate and whether to create a new
  exploration issue; never reuse or reopen the issue or RFC that prompted the exploration.
  A tracker is especially useful when the investigation will span sessions. Do not create
  it or update the project board before the user chooses. If approved, use
  [`create-issue`](../create-issue/SKILL.md), link related RFCs or issues, and make it a
  child only when the investigation is an authorized step under an active accepted RFC.
  Include a `## Exploration summary` section with this note: “Intermediate comments record
  selected findings and may be disregarded. Read the Exploration summary for the
  consolidated results.” Start the summary as pending, then use
  [`project-management`](../project-management/SKILL.md#starting-an-issue) to set its
  initial status.

## 2. Investigate interactively

- Trace current behavior through code and durable project records before proposing a fix.
  Check related issues, RFCs, analysis reports, and existing runs to avoid repeating work.
- When exploring an abandoned or failed RFC or issue, verify its current GitHub state and
  inspect the stated reason, relevant acceptance criteria, gate definition, PR discussion,
  and available evidence. Treat failure as evidence about the tested design and conditions;
  distinguish that from evidence that the whole problem is infeasible. Check whether the
  result could instead reflect an invalid assumption, input or setup mismatch, measurement
  problem or incomplete evidence. If evidence suggests an agreed gate or acceptance criterion
  needs to change, present the evidence and ask the project owner; do not change it yourself.
  Do not rerun an
  expensive failed gate until existing evidence and a meaningful change to the run are
  understood. The exploration goal is to find plausible, interesting hypotheses about the
  failure and possible mechanisms, not to test a candidate implementation.
- At each new material decision point, share concise findings, explain the alternatives and
  recommendation, and ask which direction the user wants to explore. Pause until they choose.
- The exploration phase finds and sharpens hypotheses. Inspect code, literature, existing
  runs, and observations. Clearly labelled, temporary scratch perturbations of the current
  model (for example, a knockout or parameter change) may test a causal hypothesis if they
  help distinguish explanations. Follow steps 4–6 of [`replay-ab`](../replay-ab/SKILL.md)
  for pre-registration, run ordering and interpretation; record every probe outcome, including
  null results, in the exploration issue or another user-approved investigation artifact.
  If no record was approved, ask before running a probe. Follow the one-at-a-time limit in
  `AGENTS.md`. Probes must not tune a candidate fix, become a PR, or count as acceptance
  evidence. For
  each promising hypothesis, state its proposed
  mechanism, which observations it explains, a discriminating prediction, and important
  uncertainties or alternative explanations. Do not prototype or test a candidate model
  implementation here; that belongs in an RFC investigation with a selected implementation
  and explicit diagnosis gates.
- Distinguish observations and exploratory evidence from implementation verification and
  model validation. For scientific evidence, account for controls, initial conditions,
  units, seeds, resolution and relevant uncertainty; use representative runs when stochastic
  or parameter-sensitive behavior matters. Follow `replay-ab` for its controlled before/after
  comparisons and avoid repeating matching prior runs.
- Use the worktree discovery and path-resolution procedure in
  [`continue-rfc`](../continue-rfc/SKILL.md#find-rfc-and-evidence-worktrees) to locate
  relevant registered trees. Treat other worktrees as read-only evidence unless the user
  explicitly directs otherwise.
- If relevant evidence or uncommitted work is in an inaccessible worktree, clone, or host,
  record it as a missing dependency, continue independent investigation, and ask the user
  for access or a durable handoff only when that evidence is needed to decide the next step.
- Do not implement a model change or create outcome artifacts during exploration. If the user
  approved an exploration tracker, use it as the investigation record: add comments with
  concise findings, evidence links, and changed hypotheses—not a transcript of the
  conversation. Keep its `## Exploration summary` pending until the evidence supports a useful handoff,
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
  applicable. Keep raw saves and large outputs outside Git. Do not delete saved `.mbworld`
  files without the project owner's direction; never remove another task's files. Follow
  the local-only evidence handoff procedure in
  [`AGENTS.md`](../../../AGENTS.md#workflow).

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
   issue is abandoned for this reason, set its project status to `Deprecated` and close it
   with reason `not planned`. If an RFC or issue was already abandoned, preserve its history; propose
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
  [Choosing not to pursue an issue](../project-management/SKILL.md#choosing-not-to-pursue-an-issue)
  to set its project status to `Deprecated` and close it with reason `not planned`.
- If a tracker was created, when the user confirms its exploration summary is complete, follow
  [Completing an issue](../project-management/SKILL.md#completing-an-issue) for the
  exploration tracking issue. Its deliverable is the findings summary; implementation does
  not need to be merged for the user to confirm the investigation is complete.
- Follow the relevant how-to in [`project-management`](../project-management/SKILL.md) and
  verify each requested board update.

Do not start implementation merely because an issue or RFC was created. Use
[`start-issue`](../start-issue/SKILL.md) only when the user also wants implementation to
begin.
