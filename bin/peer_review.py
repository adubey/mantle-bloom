#!/usr/bin/env python3
"""Coordinate independent multi-agent PR reviews across open windows or CLIs (#829).

An author agent requests two (or more) independent reviews on a PR, specifying each
reviewer's agent (`antigravity`, `claude`, `codex`, `muse`), model, and reasoning level.
Reviewer sessions in already-open windows either wait in standby (`listen`) and wake as
soon as a matching slot appears, or claim a slot on demand (`claim`) when prompted. With
`--spawn`, CLI-capable reviewers can also be launched headlessly in detached worktrees.

Each reviewer inspects the PR independently, posts its findings on GitHub with a standard
`## Peer Review` header, and marks its slot complete (`complete`). Meanwhile, the author
agent runs `wait` in the background so it remains available and wakes automatically once
all requested reviews have landed (detected from local state or GitHub PR comments/reviews).
See .agents/skills/peer-review/SKILL.md for the full procedure.

    bin/peer_review.py request --pr N --reviewer AGENT[:MODEL[:REASONING]] [...]
    bin/peer_review.py listen --agent AGENT [--model M] [--reasoning R]
    bin/peer_review.py claim --pr N --slot SLOT_OR_AGENT [--model M] [--reasoning R]
    bin/peer_review.py complete --pr N --slot SLOT_OR_AGENT [--url URL] [--summary TEXT]
    bin/peer_review.py wait --pr N [--timeout SEC]
    bin/peer_review.py status [--pr N] [--json]

Standard library only, so any agent can run it without activating a worktree venv.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterator

AGENT_ALIASES = {
    "agy": "antigravity",
    "antigravity": "antigravity",
    "claude": "claude",
    "claude-code": "claude",
    "codex": "codex",
    "muse": "muse",
    "muse-code": "muse",
}

HEADER_RE = re.compile(
    r"^##\s+(?:Peer\s+Review|Code\s+review)\s*(?:—|--|-|:|\()\s*([a-zA-Z0-9_-]+)",
    re.IGNORECASE | re.MULTILINE,
)
MARKER_RE = re.compile(
    r"<!--\s*mantle-bloom-peer-review:\s*pr=(\d+)\s+slot=([^\s>]+)\s+agent=([^\s>]+)\s*-->"
)


@dataclass
class ReviewSlot:
    slot: str
    agent: str
    model: str
    reasoning: str
    mode: str
    status: str = "pending"
    claimed_at: str = ""
    completed_at: str = ""
    review_url: str = ""
    summary: str = ""

    def to_dict(self) -> dict[str, str]:
        return {
            "slot": self.slot,
            "agent": self.agent,
            "model": self.model,
            "reasoning": self.reasoning,
            "mode": self.mode,
            "status": self.status,
            "claimed_at": self.claimed_at,
            "completed_at": self.completed_at,
            "review_url": self.review_url,
            "summary": self.summary,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ReviewSlot:
        return cls(
            slot=str(data["slot"]),
            agent=str(data["agent"]),
            model=str(data.get("model", "default")),
            reasoning=str(data.get("reasoning", "default")),
            mode=str(data.get("mode", "window")),
            status=str(data.get("status", "pending")),
            claimed_at=str(data.get("claimed_at", "")),
            completed_at=str(data.get("completed_at", "")),
            review_url=str(data.get("review_url", "")),
            summary=str(data.get("summary", "")),
        )


@dataclass
class ReviewRequest:
    pr: int
    requested_at: str
    author: str
    branch: str
    head_sha: str
    slots: list[ReviewSlot]

    def to_dict(self) -> dict[str, Any]:
        return {
            "pr": self.pr,
            "requested_at": self.requested_at,
            "author": self.author,
            "branch": self.branch,
            "head_sha": self.head_sha,
            "slots": [s.to_dict() for s in self.slots],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ReviewRequest:
        return cls(
            pr=int(data["pr"]),
            requested_at=str(data["requested_at"]),
            author=str(data.get("author", "")),
            branch=str(data.get("branch", "")),
            head_sha=str(data.get("head_sha", "")),
            slots=[ReviewSlot.from_dict(item) for item in data.get("slots", [])],
        )

    def all_completed(self) -> bool:
        return bool(self.slots) and all(s.status == "completed" for s in self.slots)


def now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_iso(ts: str) -> datetime | None:
    if not ts:
        return None
    cleaned = ts.strip()
    if cleaned.endswith("Z"):
        cleaned = cleaned[:-1] + "+00:00"
    try:
        dt = datetime.fromisoformat(cleaned)
        if dt.tzinfo is None:
            dt = dt.replace(tzinfo=timezone.utc)
        return dt
    except ValueError:
        return None


def canonical_agent(name: str) -> str:
    cleaned = name.strip().lower()
    if not cleaned:
        raise ValueError("Agent name cannot be empty")
    return AGENT_ALIASES.get(cleaned, cleaned)


def parse_reviewer_spec(spec: str, mode: str = "window") -> tuple[str, str, str, str]:
    """Parse 'agent[:model[:reasoning]]' into (agent, model, reasoning, mode)."""
    parts = [p.strip() for p in spec.split(":")]
    if not parts or not parts[0]:
        raise ValueError(f"Invalid reviewer specification: {spec!r}")
    if len(parts) > 3:
        raise ValueError(
            f"Reviewer specification {spec!r} has too many colons; expected AGENT[:MODEL[:REASONING]]"
        )
    agent = canonical_agent(parts[0])
    model = parts[1] if len(parts) >= 2 and parts[1] else "default"
    reasoning = parts[2] if len(parts) == 3 and parts[2] else "default"
    return agent, model, reasoning, mode


def assign_slots(specs: list[tuple[str, str, str, str]]) -> list[ReviewSlot]:
    counts: dict[str, int] = {}
    for agent, _, _, _ in specs:
        counts[agent] = counts.get(agent, 0) + 1
    seen: dict[str, int] = {}
    slots: list[ReviewSlot] = []
    for agent, model, reasoning, mode in specs:
        seen[agent] = seen.get(agent, 0) + 1
        slot_name = agent if counts[agent] == 1 else f"{agent}-{seen[agent]}"
        slots.append(
            ReviewSlot(
                slot=slot_name,
                agent=agent,
                model=model,
                reasoning=reasoning,
                mode=mode,
            )
        )
    return slots


def default_state_dir(cwd: Path | None = None) -> Path:
    env_dir = os.environ.get("MANTLE_BLOOM_REVIEWS_DIR")
    if env_dir:
        return Path(env_dir).expanduser().resolve()
    proc = subprocess.run(
        ["git", "rev-parse", "--git-common-dir"],
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Could not resolve git common dir; run inside a git checkout or pass --state-dir"
        )
    common = Path(proc.stdout.strip())
    if not common.is_absolute():
        base = cwd if cwd is not None else Path.cwd()
        common = (base / common).resolve()
    return common / "mantle-bloom-handoff" / "reviews"


def pr_file(state_dir: Path, pr: int) -> Path:
    return state_dir / f"pr-{pr}.json"


def lock_file(state_dir: Path, pr: int) -> Path:
    return state_dir / f"pr-{pr}.lock"


@contextmanager
def locked_pr(state_dir: Path, pr: int) -> Iterator[None]:
    state_dir.mkdir(parents=True, exist_ok=True)
    lpath = lock_file(state_dir, pr)
    with lpath.open("a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def load_request(state_dir: Path, pr: int) -> ReviewRequest:
    path = pr_file(state_dir, pr)
    if not path.exists():
        raise FileNotFoundError(f"No review request found for PR #{pr} at {path.name}")
    data = json.loads(path.read_text(encoding="utf-8"))
    return ReviewRequest.from_dict(data)


def save_request(state_dir: Path, req: ReviewRequest) -> Path:
    state_dir.mkdir(parents=True, exist_ok=True)
    target = pr_file(state_dir, req.pr)
    fd, tmp_name = tempfile.mkstemp(prefix=f".pr-{req.pr}-", suffix=".tmp", dir=str(state_dir))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(req.to_dict(), handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(tmp_name, target)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise
    return target


def find_matching_slot(
    req: ReviewRequest,
    target: str,
    require_pending: bool = False,
) -> ReviewSlot | None:
    norm = canonical_agent(target)
    # 1. Exact match on slot name
    for slot in req.slots:
        if slot.slot.lower() == target.lower():
            if not require_pending or slot.status == "pending":
                return slot

    # 2. Match by agent name
    matching_slots = [s for s in req.slots if s.agent == norm]
    if len(matching_slots) == 1:
        s = matching_slots[0]
        if not require_pending or s.status == "pending":
            return s
    elif len(matching_slots) > 1:
        # If multiple slots belong to the same agent, pick the first pending slot
        for s in matching_slots:
            if s.status == "pending":
                return s
        if not require_pending:
            # If not requiring pending, exact slot name was needed above
            return None
    return None


def review_header(slot: ReviewSlot) -> str:
    return (
        f"## Peer Review — {slot.agent} (model: {slot.model}, reasoning: {slot.reasoning})\n"
        f"<!-- mantle-bloom-peer-review: pr={{pr}} slot={slot.slot} agent={slot.agent} -->"
    )


def format_review_brief(req: ReviewRequest, slot: ReviewSlot) -> str:
    header = review_header(slot).format(pr=req.pr)
    lines = [
        f"=== Peer Review Ticket: PR #{req.pr} (slot: {slot.slot}) ===",
        f"Agent:     {slot.agent}",
        f"Model:     {slot.model}",
        f"Reasoning: {slot.reasoning}",
        f"Branch:    {req.branch or '(see PR)'}",
        f"Head SHA:  {req.head_sha or '(see PR)'}",
        "",
        "Instructions:",
        "1. Review the PR independently in a detached read-only worktree or via `gh pr diff`.",
        "   Do not read or rely on other peer review comments posted after this request.",
        "2. When the author reports tests passing, do not rerun the full suite; focus on",
        "   model invariants, conservation, units, save compatibility, and targeted checks.",
        f"3. Post your review on GitHub (`gh pr review {req.pr} --comment` or `gh pr comment {req.pr}`)",
        "   beginning with these exact two lines:",
        f"   {header.splitlines()[0]}",
        f"   {header.splitlines()[1]}",
        "4. Notify the author agent by marking your slot complete:",
        f"   python3 bin/peer_review.py complete --pr {req.pr} --slot {slot.slot} "
        '--url "<comment-url>" --summary "<one-line verdict>"',
    ]
    return "\n".join(lines)


def format_open_window_prompt(req: ReviewRequest, slot: ReviewSlot) -> str:
    return (
        f"Use the peer-review skill to review PR #{req.pr} as slot '{slot.slot}' "
        f"(agent: {slot.agent}, model: {slot.model}, reasoning: {slot.reasoning})."
    )


def create_review_request(
    state_dir: Path,
    pr: int,
    reviewer_specs: list[str],
    author: str = "",
    branch: str = "",
    head_sha: str = "",
    spawn: bool = False,
    force: bool = False,
    requested_at: str | None = None,
) -> ReviewRequest:
    if not reviewer_specs:
        raise ValueError("At least one --reviewer is required")
    mode = "cli" if spawn else "window"
    parsed = [parse_reviewer_spec(spec, mode=mode) for spec in reviewer_specs]
    slots = assign_slots(parsed)
    req = ReviewRequest(
        pr=pr,
        requested_at=requested_at or now_iso(),
        author=canonical_agent(author) if author else "",
        branch=branch,
        head_sha=head_sha,
        slots=slots,
    )
    with locked_pr(state_dir, pr):
        path = pr_file(state_dir, pr)
        if path.exists() and not force:
            existing = load_request(state_dir, pr)
            has_active = any(s.status in ("claimed", "completed") for s in existing.slots)
            if has_active:
                raise RuntimeError(
                    f"Review request for PR #{pr} already exists with active/completed reviews. "
                    "Use --force to overwrite."
                )
        save_request(state_dir, req)
    return req


def claim_slot(
    state_dir: Path,
    pr: int,
    target: str,
    model: str = "",
    reasoning: str = "",
    require_pending: bool = False,
) -> tuple[ReviewRequest, ReviewSlot] | None:
    with locked_pr(state_dir, pr):
        req = load_request(state_dir, pr)
        slot = find_matching_slot(req, target, require_pending=require_pending)
        if slot is None:
            return None
        if slot.status == "completed":
            # Reopen or update model/reasoning if requested
            if model:
                slot.model = model
            if reasoning:
                slot.reasoning = reasoning
            save_request(state_dir, req)
            return req, slot
        slot.status = "claimed"
        slot.claimed_at = now_iso()
        if model:
            slot.model = model
        if reasoning:
            slot.reasoning = reasoning
        save_request(state_dir, req)
        return req, slot


def claim_next_pending_for_agent(
    state_dir: Path,
    agent: str,
    model: str = "",
    reasoning: str = "",
) -> tuple[ReviewRequest, ReviewSlot] | None:
    if not state_dir.exists():
        return None
    candidates = sorted(state_dir.glob("pr-*.json"))
    for path in candidates:
        match = re.match(r"^pr-(\d+)\.json$", path.name)
        if not match:
            continue
        pr_num = int(match.group(1))
        claimed = claim_slot(
            state_dir,
            pr_num,
            agent,
            model=model,
            reasoning=reasoning,
            require_pending=True,
        )
        if claimed is not None:
            return claimed
    return None


def listen_for_slot(
    state_dir: Path,
    agent: str,
    model: str = "",
    reasoning: str = "",
    timeout: float = 3600.0,
    poll_interval: float = 1.0,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> tuple[ReviewRequest, ReviewSlot]:
    deadline = time.monotonic() + timeout
    while True:
        claimed = claim_next_pending_for_agent(
            state_dir, agent, model=model, reasoning=reasoning
        )
        if claimed is not None:
            return claimed
        if time.monotonic() >= deadline:
            raise TimeoutError(f"Timed out after {timeout:.0f}s waiting for a {agent} review slot")
        sleep_fn(poll_interval)


def complete_slot(
    state_dir: Path,
    pr: int,
    target: str,
    url: str = "",
    summary: str = "",
) -> tuple[ReviewRequest, ReviewSlot]:
    with locked_pr(state_dir, pr):
        req = load_request(state_dir, pr)
        slot = find_matching_slot(req, target, require_pending=False)
        if slot is None:
            raise ValueError(f"No review slot matching {target!r} on PR #{pr}")
        slot.status = "completed"
        slot.completed_at = now_iso()
        if url:
            slot.review_url = url
        if summary:
            slot.summary = summary
        save_request(state_dir, req)
        return req, slot


def match_gh_entry_to_slot(
    req: ReviewRequest, body: str, created_at: str, allow_completed: bool = False
) -> ReviewSlot | None:
    dt_created = parse_iso(created_at)
    dt_requested = parse_iso(req.requested_at)
    if dt_created and dt_requested and dt_created < dt_requested:
        return None
    marker = MARKER_RE.search(body)
    if marker:
        marker_pr = int(marker.group(1))
        if marker_pr != req.pr:
            return None
        slot_name = marker.group(2)
        agent_name = marker.group(3)
        slot = find_matching_slot(req, slot_name, require_pending=False)
        if slot is not None and (allow_completed or slot.status != "completed"):
            return slot
        slot = find_matching_slot(req, agent_name, require_pending=False)
        if slot is not None and (allow_completed or slot.status != "completed"):
            return slot
    header = HEADER_RE.search(body)
    if header:
        agent_token = header.group(1)
        for slot in req.slots:
            if (allow_completed or slot.status != "completed") and (
                slot.slot.lower() == agent_token.lower()
                or slot.agent == canonical_agent(agent_token)
            ):
                return slot
    return None


def extract_first_summary_line(body: str) -> str:
    return next(
        (
            line.strip()
            for line in body.splitlines()
            if line.strip()
            and not line.strip().startswith("## Peer Review")
            and not line.strip().startswith("## Code review")
            and not line.strip().startswith("<!--")
        ),
        "Posted peer review on GitHub",
    )[:240]


def sync_from_github_payload(
    state_dir: Path,
    pr: int,
    payload: dict[str, Any],
) -> ReviewRequest:
    """Reconcile PR reviews and comments from `gh pr view --json reviews,comments`."""
    entries: list[tuple[str, str, str]] = []
    for item in payload.get("reviews", []) or []:
        body = str(item.get("body") or "").strip()
        if body:
            entries.append(
                (body, str(item.get("submittedAt") or ""), str(item.get("url") or ""))
            )
    for item in payload.get("comments", []) or []:
        body = str(item.get("body") or "").strip()
        if body:
            entries.append(
                (body, str(item.get("createdAt") or ""), str(item.get("url") or ""))
            )

    with locked_pr(state_dir, pr):
        req = load_request(state_dir, pr)
        changed = False
        for body, timestamp, url in entries:
            slot = match_gh_entry_to_slot(req, body, timestamp, allow_completed=False)
            if slot is None:
                # Check if it matches an already completed slot needing backfill
                slot = match_gh_entry_to_slot(req, body, timestamp, allow_completed=True)
                if slot is not None:
                    if url and not slot.review_url:
                        slot.review_url = url
                        changed = True
                    if not slot.summary:
                        slot.summary = extract_first_summary_line(body)
                        changed = True
                continue

            slot.status = "completed"
            slot.completed_at = timestamp or now_iso()
            if url and not slot.review_url:
                slot.review_url = url
            if not slot.summary:
                slot.summary = extract_first_summary_line(body)
            changed = True
        if changed:
            save_request(state_dir, req)
        return req


def fetch_github_pr_payload(pr: int) -> dict[str, Any] | None:
    proc = subprocess.run(
        ["gh", "pr", "view", str(pr), "--json", "reviews,comments"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or not proc.stdout.strip():
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None


def wait_for_reviews(
    state_dir: Path,
    pr: int,
    timeout: float = 3600.0,
    poll_interval: float = 5.0,
    gh_fetcher: Callable[[int], dict[str, Any] | None] | None = fetch_github_pr_payload,
    sleep_fn: Callable[[float], None] = time.sleep,
) -> ReviewRequest:
    deadline = time.monotonic() + timeout
    while True:
        req = load_request(state_dir, pr)
        if req.all_completed():
            return req
        if gh_fetcher is not None:
            payload = gh_fetcher(pr)
            if payload is not None:
                req = sync_from_github_payload(state_dir, pr, payload)
                if req.all_completed():
                    return req
        if time.monotonic() >= deadline:
            pending = ", ".join(s.slot for s in req.slots if s.status != "completed")
            raise TimeoutError(
                f"Timed out after {timeout:.0f}s waiting for PR #{pr} reviews (pending: {pending})"
            )
        sleep_fn(poll_interval)


def build_cli_command(req: ReviewRequest, slot: ReviewSlot, worktree: Path) -> list[str]:
    prompt = format_review_brief(req, slot)
    if slot.agent == "codex":
        cmd = ["codex", "exec", "--cd", str(worktree)]
        if slot.model != "default":
            cmd.extend(["-m", slot.model])
        if slot.reasoning != "default":
            cmd.extend(["-c", f'model_reasoning_effort="{slot.reasoning}"'])
        cmd.append(prompt)
        return cmd
    if slot.agent == "claude":
        cmd = ["claude", "-p"]
        if slot.model != "default":
            cmd.extend(["--model", slot.model])
        if slot.reasoning != "default":
            cmd.extend(["--effort", slot.reasoning])
        cmd.append(prompt)
        return cmd
    if slot.agent == "antigravity":
        cmd = ["agy"]
        if slot.model != "default":
            cmd.extend(["--model", slot.model])
        cmd.append(prompt)
        return cmd
    raise ValueError(
        f"Agent {slot.agent!r} does not have a headless CLI builder; use open-window mode instead"
    )


def detect_git_metadata() -> tuple[str, str]:
    branch_proc = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    sha_proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    branch = branch_proc.stdout.strip() if branch_proc.returncode == 0 else ""
    sha = sha_proc.stdout.strip() if sha_proc.returncode == 0 else ""
    return branch, sha


def resolve_primary_checkout(state_dir: Path | None = None) -> Path:
    if state_dir is not None and ".git" in state_dir.parts:
        idx = state_dir.parts.index(".git")
        return Path(*state_dir.parts[:idx])
    proc = subprocess.run(
        ["git", "rev-parse", "--git-common-dir"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode == 0 and proc.stdout.strip():
        common = Path(proc.stdout.strip()).resolve()
        if common.name == ".git":
            return common.parent
        if ".git" in common.parts:
            idx = common.parts.index(".git")
            return Path(*common.parts[:idx])
    if state_dir is not None and len(state_dir.parents) >= 3:
        return state_dir.parents[2]
    return Path.cwd().resolve()


def review_worktree_path(
    state_dir: Path, pr: int, slot: str, worktree_root: Path | None = None
) -> Path:
    """Where a spawned reviewer's detached tree goes, per the worktree skill: a sibling of
    the primary checkout, or `<worktree_root>/review-<pr>-<slot>` when the caller's sandbox
    is scoped to the workspace root. The state dir sits in .git, never a place for a tree."""
    name = f"review-{pr}-{slot}"
    if worktree_root is not None:
        return worktree_root / name
    primary = resolve_primary_checkout(state_dir)
    return primary.parent / f"{primary.name}-{name}"


def spawn_cli_reviewer(
    state_dir: Path,
    req: ReviewRequest,
    slot: ReviewSlot,
    worktree_root: Path | None = None,
) -> subprocess.Popen[bytes]:
    wt_dir = review_worktree_path(state_dir, req.pr, slot.slot, worktree_root)
    # 1. Build command upfront to catch unsupported CLI builders (e.g. muse) before any side-effects
    cmd = build_cli_command(req, slot, wt_dir)

    created_wt = False
    log_handle = None
    try:
        if not wt_dir.exists():
            subprocess.run(
                ["git", "fetch", "origin", f"pull/{req.pr}/head"],
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                ["git", "worktree", "add", "--detach", str(wt_dir), "FETCH_HEAD"],
                check=True,
                capture_output=True,
                text=True,
            )
            created_wt = True

        claim_slot(state_dir, req.pr, slot.slot)
        log_path = state_dir / f"pr-{req.pr}-{slot.slot}.log"
        log_handle = log_path.open("ab")
        proc = subprocess.Popen(
            cmd,
            cwd=str(wt_dir),
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        return proc
    except Exception:
        if log_handle is not None:
            try:
                log_handle.close()
            except Exception:
                pass
        if created_wt and wt_dir.exists():
            subprocess.run(["git", "worktree", "remove", "--force", str(wt_dir)], check=False)
        # Release slot back to pending
        with locked_pr(state_dir, req.pr):
            cur_req = load_request(state_dir, req.pr)
            s = find_matching_slot(cur_req, slot.slot)
            if s and s.status == "claimed":
                s.status = "pending"
                s.claimed_at = ""
                save_request(state_dir, cur_req)
        raise


def format_status_report(req: ReviewRequest) -> str:
    lines = [
        f"PR #{req.pr} (requested {req.requested_at}, head {req.head_sha[:8] if req.head_sha else 'unknown'}):"
    ]
    for slot in req.slots:
        detail = f"  - {slot.slot}: agent={slot.agent} model={slot.model} reasoning={slot.reasoning} status={slot.status}"
        if slot.review_url:
            detail += f" url={slot.review_url}"
        lines.append(detail)
        if slot.summary:
            lines.append(f"    summary: {slot.summary}")
    return "\n".join(lines)


def cmd_request(args: argparse.Namespace, state_dir: Path) -> int:
    branch, sha = detect_git_metadata()
    req = create_review_request(
        state_dir=state_dir,
        pr=args.pr,
        reviewer_specs=args.reviewer,
        author=args.author,
        branch=args.branch or branch,
        head_sha=args.head_sha or sha,
        spawn=args.spawn,
        force=args.force,
    )
    if args.spawn:
        for slot in req.slots:
            spawn_cli_reviewer(state_dir, req, slot, args.worktree_root)
            print(
                f"Spawned headless CLI reviewer for slot '{slot.slot}' "
                f"({slot.agent}, model={slot.model}, reasoning={slot.reasoning})."
            )
        return 0

    if args.claim_wait > 0:
        deadline = time.monotonic() + args.claim_wait
        while time.monotonic() < deadline:
            req = load_request(state_dir, args.pr)
            if all(s.status != "pending" for s in req.slots):
                break
            time.sleep(0.25)
        req = load_request(state_dir, args.pr)

    print(f"Created review request for PR #{req.pr} with {len(req.slots)} slot(s):")
    for slot in req.slots:
        if slot.status != "pending":
            print(
                f"- [{slot.slot}] Claimed automatically by standby {slot.agent} window "
                f"(model={slot.model}, reasoning={slot.reasoning})."
            )
        else:
            print(
                f"- [{slot.slot}] In your open '{slot.agent}' window, select model "
                f"'{slot.model}' and reasoning '{slot.reasoning}', then paste:"
            )
            print(f"    {format_open_window_prompt(req, slot)}")
    return 0


def cmd_listen(args: argparse.Namespace, state_dir: Path) -> int:
    req, slot = listen_for_slot(
        state_dir=state_dir,
        agent=args.agent,
        model=args.model,
        reasoning=args.reasoning,
        timeout=args.timeout,
        poll_interval=args.poll_interval,
    )
    print(format_review_brief(req, slot))
    return 0


def cmd_claim(args: argparse.Namespace, state_dir: Path) -> int:
    claimed = claim_slot(
        state_dir=state_dir,
        pr=args.pr,
        target=args.slot,
        model=args.model,
        reasoning=args.reasoning,
        require_pending=False,
    )
    if claimed is None:
        print(f"No matching slot {args.slot!r} on PR #{args.pr}", file=sys.stderr)
        return 1
    req, slot = claimed
    print(format_review_brief(req, slot))
    return 0


def cmd_complete(args: argparse.Namespace, state_dir: Path) -> int:
    req, slot = complete_slot(
        state_dir=state_dir,
        pr=args.pr,
        target=args.slot,
        url=args.url,
        summary=args.summary,
    )
    remaining = sum(1 for s in req.slots if s.status != "completed")
    print(
        f"Marked slot '{slot.slot}' completed on PR #{req.pr} "
        f"({remaining} slot(s) remaining)."
    )
    return 0


def cmd_wait(args: argparse.Namespace, state_dir: Path) -> int:
    fetcher = None if args.no_gh else fetch_github_pr_payload
    req = wait_for_reviews(
        state_dir=state_dir,
        pr=args.pr,
        timeout=args.timeout,
        poll_interval=args.poll_interval,
        gh_fetcher=fetcher,
    )
    print(f"All {len(req.slots)} peer review(s) completed for PR #{req.pr}:")
    print(format_status_report(req))
    return 0


def cmd_status(args: argparse.Namespace, state_dir: Path) -> int:
    if args.pr is not None:
        req = load_request(state_dir, args.pr)
        if args.json:
            print(json.dumps(req.to_dict(), indent=2, sort_keys=True))
        else:
            print(format_status_report(req))
        return 0

    requests: list[ReviewRequest] = []
    if state_dir.exists():
        for path in sorted(state_dir.glob("pr-*.json")):
            if re.match(r"^pr-\d+\.json$", path.name):
                requests.append(ReviewRequest.from_dict(json.loads(path.read_text("utf-8"))))
    if args.json:
        print(json.dumps([r.to_dict() for r in requests], indent=2, sort_keys=True))
    elif not requests:
        print("No active peer review requests.")
    else:
        for req in requests:
            print(format_status_report(req))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=None,
        help="Override the shared handoff reviews directory",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_req = sub.add_parser("request", help="Create a peer review request for a PR")
    p_req.add_argument("--pr", type=int, required=True, help="Pull request number")
    p_req.add_argument(
        "--reviewer",
        action="append",
        required=True,
        help="Reviewer spec AGENT[:MODEL[:REASONING]] (repeatable)",
    )
    p_req.add_argument("--author", default="", help="Author agent identifier")
    p_req.add_argument("--branch", default="", help="PR branch name")
    p_req.add_argument("--head-sha", default="", help="PR head commit SHA")
    p_req.add_argument(
        "--spawn",
        action="store_true",
        help="Spawn headless CLI reviewers in detached worktrees instead of open windows",
    )
    p_req.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing review request even if slots are claimed or completed",
    )
    p_req.add_argument(
        "--worktree-root",
        type=Path,
        default=None,
        help="With --spawn, put review worktrees here instead of beside the primary checkout",
    )
    p_req.add_argument(
        "--claim-wait",
        type=float,
        default=2.0,
        help="Seconds to wait for standby open windows to auto-claim slots (default: 2.0)",
    )

    p_listen = sub.add_parser("listen", help="Wait in standby for a review slot for an agent")
    p_listen.add_argument("--agent", required=True, help="Agent name (antigravity, claude, codex, muse)")
    p_listen.add_argument("--model", default="", help="Override model active in this window")
    p_listen.add_argument("--reasoning", default="", help="Override reasoning level active in this window")
    p_listen.add_argument("--timeout", type=float, default=3600.0, help="Seconds to wait")
    p_listen.add_argument("--poll-interval", type=float, default=1.0, help="Poll interval in seconds")

    p_claim = sub.add_parser("claim", help="Claim a specific PR review slot in an open window")
    p_claim.add_argument("--pr", type=int, required=True, help="Pull request number")
    p_claim.add_argument("--slot", required=True, help="Slot name or agent name")
    p_claim.add_argument("--model", default="", help="Record model used by this window")
    p_claim.add_argument("--reasoning", default="", help="Record reasoning level used by this window")

    p_done = sub.add_parser("complete", help="Mark a review slot completed after posting to GitHub")
    p_done.add_argument("--pr", type=int, required=True, help="Pull request number")
    p_done.add_argument("--slot", required=True, help="Slot name or agent name")
    p_done.add_argument("--url", default="", help="GitHub comment or review URL")
    p_done.add_argument("--summary", default="", help="One-line review summary or verdict")

    p_wait = sub.add_parser("wait", help="Wait until all requested reviews for a PR complete")
    p_wait.add_argument("--pr", type=int, required=True, help="Pull request number")
    p_wait.add_argument("--timeout", type=float, default=3600.0, help="Seconds to wait")
    p_wait.add_argument("--poll-interval", type=float, default=5.0, help="Poll interval in seconds")
    p_wait.add_argument(
        "--no-gh",
        action="store_true",
        help="Wait on local handoff state only without polling GitHub PR comments",
    )

    p_status = sub.add_parser("status", help="Show review slot statuses")
    p_status.add_argument("--pr", type=int, default=None, help="Pull request number")
    p_status.add_argument("--json", action="store_true", help="Emit JSON output")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    state_dir = args.state_dir.resolve() if args.state_dir else default_state_dir()
    dispatch = {
        "request": cmd_request,
        "listen": cmd_listen,
        "claim": cmd_claim,
        "complete": cmd_complete,
        "wait": cmd_wait,
        "status": cmd_status,
    }
    return dispatch[args.command](args, state_dir)


if __name__ == "__main__":
    sys.exit(main())
