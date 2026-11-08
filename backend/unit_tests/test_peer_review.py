"""Tests for bin/peer_review.py (#829): multi-agent review slot assignment, atomic
standby claiming, GitHub comment reconciliation, and headless CLI command builders."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "peer_review.py"
spec = importlib.util.spec_from_file_location("peer_review", SCRIPT)
assert spec is not None and spec.loader is not None
pr_mod = importlib.util.module_from_spec(spec)
sys.modules["peer_review"] = pr_mod
spec.loader.exec_module(pr_mod)


def test_parse_reviewer_spec_and_aliases():
    assert pr_mod.parse_reviewer_spec("muse:opus:high") == (
        "muse",
        "opus",
        "high",
        "window",
    )
    assert pr_mod.parse_reviewer_spec("agy:pro") == (
        "antigravity",
        "pro",
        "default",
        "window",
    )
    assert pr_mod.parse_reviewer_spec("claude-code", mode="cli") == (
        "claude",
        "default",
        "default",
        "cli",
    )
    with pytest.raises(ValueError, match="Invalid reviewer"):
        pr_mod.parse_reviewer_spec("")
    with pytest.raises(ValueError, match="too many colons"):
        pr_mod.parse_reviewer_spec("codex:m:r:extra")


def test_assign_slots_numbers_duplicate_agents_only():
    distinct = pr_mod.assign_slots(
        [
            pr_mod.parse_reviewer_spec("muse:opus:high"),
            pr_mod.parse_reviewer_spec("codex:gpt-5.4:xhigh"),
        ]
    )
    assert [s.slot for s in distinct] == ["muse", "codex"]

    repeated = pr_mod.assign_slots(
        [
            pr_mod.parse_reviewer_spec("codex:gpt-5.4:high"),
            pr_mod.parse_reviewer_spec("codex:o3:medium"),
        ]
    )
    assert [s.slot for s in repeated] == ["codex-1", "codex-2"]


def test_standby_listeners_claim_distinct_slots(tmp_path):
    pr_mod.create_review_request(
        tmp_path,
        pr=812,
        reviewer_specs=["codex:gpt-5.4:high", "codex:o3:medium"],
        author="claude",
        branch="claude/issue-812-example",
        head_sha="deadbeef1234",
    )

    first = pr_mod.listen_for_slot(tmp_path, agent="codex", timeout=0.1, poll_interval=0.01)
    second = pr_mod.listen_for_slot(tmp_path, agent="codex", timeout=0.1, poll_interval=0.01)
    assert first[1].slot == "codex-1"
    assert second[1].slot == "codex-2"

    with pytest.raises(TimeoutError):
        pr_mod.listen_for_slot(tmp_path, agent="codex", timeout=0.01, poll_interval=0.01)


def test_local_complete_unblocks_wait(tmp_path):
    pr_mod.create_review_request(
        tmp_path,
        pr=815,
        reviewer_specs=["muse:opus:high", "codex:gpt-5.4:high"],
        author="antigravity",
    )
    pr_mod.claim_slot(tmp_path, pr=815, target="muse", model="opus-4.7", reasoning="high")
    pr_mod.complete_slot(
        tmp_path,
        pr=815,
        target="muse",
        url="https://github.com/adubey/mantle-bloom/pull/815#issuecomment-1",
        summary="Checked conservation and save compatibility; looks sound.",
    )
    with pytest.raises(TimeoutError):
        pr_mod.wait_for_reviews(
            tmp_path,
            pr=815,
            timeout=0.01,
            poll_interval=0.01,
            gh_fetcher=None,
        )

    pr_mod.complete_slot(
        tmp_path,
        pr=815,
        target="codex",
        url="https://github.com/adubey/mantle-bloom/pull/815#issuecomment-2",
        summary="Found one edge case in boundary flux.",
    )
    done = pr_mod.wait_for_reviews(
        tmp_path,
        pr=815,
        timeout=0.1,
        poll_interval=0.01,
        gh_fetcher=None,
    )
    assert done.all_completed()
    assert [s.review_url.endswith(("#issuecomment-1", "#issuecomment-2")) for s in done.slots] == [
        True,
        True,
    ]


def test_github_payload_reconciles_reviews_and_ignores_stale_comments(tmp_path):
    pr_mod.create_review_request(
        tmp_path,
        pr=820,
        reviewer_specs=["muse:opus:high", "codex:gpt-5.4:high"],
        requested_at="2026-11-07T19:00:00Z",
    )

    stale_and_partial_payload = {
        "reviews": [
            {
                "body": "## Peer Review — muse (model: opus, reasoning: high)\nOld pass before push.",
                "submittedAt": "2026-11-07T18:30:00Z",
                "url": "https://github.com/adubey/mantle-bloom/pull/820#pullrequestreview-10",
            },
            {
                "body": (
                    "## Peer Review — muse (model: opus, reasoning: high)\n"
                    "<!-- mantle-bloom-peer-review: pr=820 slot=muse agent=muse -->\n"
                    "No conservation issues found on the new commit."
                ),
                "submittedAt": "2026-11-07T19:05:00Z",
                "url": "https://github.com/adubey/mantle-bloom/pull/820#pullrequestreview-11",
            },
        ],
        "comments": [
            {
                "body": (
                    "## Peer Review — codex (model: gpt-5.4, reasoning: high)\n"
                    "Checked timestep invariance and bounds; one minor nit."
                ),
                "createdAt": "2026-11-07T19:06:00Z",
                "url": "https://github.com/adubey/mantle-bloom/pull/820#issuecomment-99",
            }
        ],
    }

    done = pr_mod.wait_for_reviews(
        tmp_path,
        pr=820,
        timeout=0.1,
        poll_interval=0.01,
        gh_fetcher=lambda _pr: stale_and_partial_payload,
    )
    assert done.all_completed()
    by_slot = {s.slot: s for s in done.slots}
    assert by_slot["muse"].review_url.endswith("#pullrequestreview-11")
    assert by_slot["muse"].summary == "No conservation issues found on the new commit."
    assert by_slot["codex"].review_url.endswith("#issuecomment-99")
    assert by_slot["codex"].summary == "Checked timestep invariance and bounds; one minor nit."


def test_build_cli_command_passes_model_and_reasoning(tmp_path):
    req = pr_mod.ReviewRequest(
        pr=825,
        requested_at="2026-11-07T19:00:00Z",
        author="antigravity",
        branch="antigravity/issue-825-test",
        head_sha="abcdef123456",
        slots=[],
    )
    codex_slot = pr_mod.ReviewSlot("codex", "codex", "gpt-5.4", "xhigh", "cli")
    codex_cmd = pr_mod.build_cli_command(req, codex_slot, tmp_path)
    assert codex_cmd[:7] == [
        "codex",
        "exec",
        "--cd",
        str(tmp_path),
        "-m",
        "gpt-5.4",
        "-c",
    ]
    assert codex_cmd[7] == 'model_reasoning_effort="xhigh"'
    assert "PR #825" in codex_cmd[8]

    claude_slot = pr_mod.ReviewSlot("claude", "claude", "opus", "high", "cli")
    claude_cmd = pr_mod.build_cli_command(req, claude_slot, tmp_path)
    assert claude_cmd[:6] == ["claude", "-p", "--model", "opus", "--effort", "high"]

    agy_slot = pr_mod.ReviewSlot("antigravity", "antigravity", "pro", "default", "cli")
    agy_cmd = pr_mod.build_cli_command(req, agy_slot, tmp_path)
    assert agy_cmd[:3] == ["agy", "--model", "pro"]

    muse_slot = pr_mod.ReviewSlot("muse", "muse", "opus", "high", "window")
    with pytest.raises(ValueError, match="open-window mode"):
        pr_mod.build_cli_command(req, muse_slot, tmp_path)


def test_cli_end_to_end_workflow(tmp_path, capsys):
    base = ["--state-dir", str(tmp_path)]
    assert (
        pr_mod.main(
            [
                *base,
                "request",
                "--pr",
                "829",
                "--author",
                "antigravity",
                "--reviewer",
                "muse:opus:high",
                "--reviewer",
                "codex:gpt-5.4:high",
                "--claim-wait",
                "0",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out
    assert "Use the peer-review skill to review PR #829" in out

    assert pr_mod.main([*base, "claim", "--pr", "829", "--slot", "muse"]) == 0
    assert "## Peer Review — muse" in capsys.readouterr().out

    assert (
        pr_mod.main(
            [
                *base,
                "complete",
                "--pr",
                "829",
                "--slot",
                "muse",
                "--summary",
                "LGTM from muse",
            ]
        )
        == 0
    )
    assert (
        pr_mod.main(
            [
                *base,
                "complete",
                "--pr",
                "829",
                "--slot",
                "codex",
                "--summary",
                "LGTM from codex",
            ]
        )
        == 0
    )
    capsys.readouterr()

    assert pr_mod.main([*base, "wait", "--pr", "829", "--no-gh", "--timeout", "1"]) == 0
    wait_out = capsys.readouterr().out
    assert "All 2 peer review(s) completed for PR #829" in wait_out

    assert pr_mod.main([*base, "status", "--pr", "829", "--json"]) == 0
    status_json = json.loads(capsys.readouterr().out)
    assert status_json["pr"] == 829
    assert [s["status"] for s in status_json["slots"]] == ["completed", "completed"]
