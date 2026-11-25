"""Unit tests for backend/app/animation_job.py (#1033)."""

from __future__ import annotations

import json

from app.animation_job import AnimationJob


def test_subscribe_receives_progress_heartbeat_and_done():
    job = AnimationJob(total=3, unbounded=False, projection="behrmann", view="elevation", elapsed_years=100.0)
    sub = job.subscribe(replay_latest=False)

    p1 = json.dumps({"type": "progress", "frame": 1, "total": 3, "elapsed_years": 200.0}) + "\n"
    hb = json.dumps({"type": "heartbeat", "frame": 2, "total": 3, "elapsed_years": 250.0}) + "\n"
    done = json.dumps({"type": "done", "mime": "video/mp4", "video_base64": "AAAA", "stopped_early": False}) + "\n"

    job.publish_progress(p1, frame=1, elapsed_years=200.0)
    assert job.summary()["frame"] == 1
    assert job.summary()["active"] is True

    job.publish_heartbeat(hb, elapsed_years=250.0)
    job.mark_stopping()
    assert job.summary()["stopping"] is True

    job.finish_done(done, stopped_early=True, error=None, elapsed_years=300.0)
    assert job.is_active() is False
    summary = job.summary()
    assert summary["has_result"] is True
    assert summary["stopped_early"] is True
    assert summary["stopping"] is False
    assert summary["elapsed_years"] == 300.0

    lines = list(job.iter_stream(sub))
    assert lines == [p1, hb, done]


def test_reconnect_mid_run_replays_latest_progress_and_follows_completion():
    job = AnimationJob(total=4, unbounded=True, projection="eckert4", view="biome")
    p1 = json.dumps({"type": "progress", "frame": 1, "total": 4, "elapsed_years": 10.0}) + "\n"
    p2 = json.dumps({"type": "progress", "frame": 2, "total": 4, "elapsed_years": 20.0}) + "\n"
    done = json.dumps({"type": "done", "mime": "video/mp4", "video_base64": "BBBB", "stopped_early": False}) + "\n"

    orig = job.subscribe(replay_latest=False)
    job.publish_progress(p1, frame=1, elapsed_years=10.0)
    # Original client disconnects mid-run; job keeps going.
    job.unsubscribe(orig)

    job.publish_progress(p2, frame=2, elapsed_years=20.0)
    reconnected = job.subscribe(replay_latest=True)
    job.finish_done(done, stopped_early=False, error=None, elapsed_years=30.0)

    # Reconnected client gets latest snapshot (p2, not stale p1) followed by done.
    assert list(job.iter_stream(reconnected)) == [p2, done]


def test_reconnect_after_completion_immediately_delivers_terminal_line():
    job = AnimationJob(total=2, unbounded=False, projection="behrmann", view="elevation")
    done = json.dumps({"type": "done", "mime": "video/mp4", "video_base64": "CCCC", "stopped_early": False}) + "\n"
    job.finish_done(done, stopped_early=False, error=None, elapsed_years=50.0)

    late = job.subscribe(replay_latest=True)
    assert list(job.iter_stream(late)) == [done]


def test_finish_error_and_finish_empty():
    err_job = AnimationJob(total=2, unbounded=False, projection="behrmann", view="elevation")
    err_line = json.dumps({"type": "error", "detail": "boom"}) + "\n"
    err_job.finish_error(err_line, detail="boom")
    assert err_job.summary()["error"] == "boom"
    assert err_job.summary()["has_result"] is False
    assert list(err_job.iter_stream(err_job.subscribe(replay_latest=True))) == [err_line]

    empty_job = AnimationJob(total=2, unbounded=False, projection="behrmann", view="elevation")
    sub = empty_job.subscribe(replay_latest=False)
    empty_job.finish_empty()
    assert empty_job.is_active() is False
    assert list(empty_job.iter_stream(sub)) == []
