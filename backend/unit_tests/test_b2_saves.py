"""Tests for bin/b2_saves.py (#708): the rules that decide what may be deleted, on a fake
bucket, and the B2 request shapes, against a local HTTP stub."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
import threading
import zipfile
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "b2_saves.py"
spec = importlib.util.spec_from_file_location("b2_saves", SCRIPT)
assert spec is not None and spec.loader is not None
b2 = importlib.util.module_from_spec(spec)
sys.modules["b2_saves"] = b2  # dataclasses looks the module up by name
spec.loader.exec_module(b2)

NOW = 2_000_000_000.0
DAY = 86400


class FakeBucket:
    def __init__(self):
        self.files = {}  # name -> (Remote, bytes)
        self.corrupt_downloads = False
        self.fail_listing = False

    def list_files(self):
        if self.fail_listing:
            raise b2.B2Error(503, "service_unavailable", "down")
        return [remote for remote, _ in self.files.values()]

    def upload(self, path, name, sha1, info):
        remote = b2.Remote(name, path.stat().st_size, sha1, 0, info)
        self.files[name] = (remote, path.read_bytes())
        return remote

    def download(self, name, dest):
        data = self.files[name][1] + (b"x" if self.corrupt_downloads else b"")
        dest.write_bytes(data)
        return hashlib.sha1(data).hexdigest()


def save(directory, name, content=None, age_days=0.0):
    path = directory / name
    path.write_bytes((content or name).encode())
    when = NOW - age_days * DAY
    os.utime(path, (when, when))
    return path


def config(saves_dir, keep_days=14.0, keep_newest=0):
    return b2.Config("id", "key", "bucket", saves_dir, keep_days, keep_newest)


def test_purge_deletes_old_saves_held_by_the_bucket(tmp_path):
    old = save(tmp_path, "mantle-bloom-seed1-1y.bloom", age_days=30)
    recent = save(tmp_path, "mantle-bloom-seed2-1y.bloom", age_days=1)
    bucket = FakeBucket()
    b2.sync(bucket, tmp_path, NOW)

    purged = b2.purge(bucket, config(tmp_path), NOW, dry_run=False, pins=set())

    assert purged == [old]
    assert not old.exists() and recent.exists()
    assert "mantle-bloom-seed1-1y.bloom" in bucket.files


def test_purge_keeps_an_old_save_missing_from_the_bucket(tmp_path):
    old = save(tmp_path, "mantle-bloom-seed1-1y.bloom", age_days=30)
    assert b2.purge(FakeBucket(), config(tmp_path), NOW, dry_run=False, pins=set()) == []
    assert old.exists()


def test_purge_compares_contents_not_names(tmp_path):
    """A same-named, same-sized bucket file isn't a backup of different bytes."""
    old = save(tmp_path, "mantle-bloom-seed1-1y.bloom", content="local!", age_days=30)
    bucket = FakeBucket()
    bucket.files[old.name] = (b2.Remote(old.name, 6, "0" * 40, 0), b"cloud!")
    assert b2.purge(bucket, config(tmp_path), NOW, dry_run=False, pins=set()) == []
    assert old.exists()


def test_purge_deletes_nothing_when_the_listing_fails(tmp_path):
    old = save(tmp_path, "mantle-bloom-seed1-1y.bloom", age_days=30)
    bucket = FakeBucket()
    b2.sync(bucket, tmp_path, NOW)
    bucket.fail_listing = True
    with pytest.raises(b2.B2Error):
        b2.purge(bucket, config(tmp_path), NOW, dry_run=False, pins=set())
    assert old.exists()


def test_purge_honours_dry_run_pins_and_keep_newest(tmp_path):
    names = [f"mantle-bloom-seed{i}-1y.bloom" for i in range(4)]
    paths = [save(tmp_path, name, age_days=30 + i) for i, name in enumerate(names)]
    bucket = FakeBucket()
    b2.sync(bucket, tmp_path, NOW)

    cfg = config(tmp_path, keep_newest=1)
    dry = b2.purge(bucket, cfg, NOW, dry_run=True, pins={names[1]})
    assert dry == paths[2:]
    assert all(path.exists() for path in paths)

    b2.purge(bucket, cfg, NOW, dry_run=False, pins={names[1]})
    assert [path.exists() for path in paths] == [True, True, False, False]


def test_purge_ignores_files_that_are_not_saves(tmp_path):
    other = save(tmp_path, "tax-return.pdf", age_days=400)
    assert b2.purge(FakeBucket(), config(tmp_path), NOW, dry_run=False, pins=set()) == []
    assert other.exists()


def test_an_unreadable_saves_directory_is_an_error_not_an_empty_list(tmp_path, monkeypatch):
    """macOS blocks cron from ~/Downloads; that must not look like having no saves."""

    def denied(_path):
        raise PermissionError("Operation not permitted")

    monkeypatch.setattr(b2.os, "listdir", denied)
    with pytest.raises(PermissionError):
        b2.local_saves(tmp_path)


def test_a_missing_saves_directory_has_no_saves(tmp_path):
    assert b2.local_saves(tmp_path / "absent") == []


def test_sync_skips_a_save_still_being_written(tmp_path):
    save(tmp_path, "mantle-bloom-seed1-1y.bloom", age_days=0)
    bucket = FakeBucket()
    assert b2.sync(bucket, tmp_path, NOW) == 0
    assert bucket.files == {}


def test_push_keeps_different_saves_with_one_name_apart(tmp_path):
    first = save(tmp_path, "mantle-bloom-seed1-1y.bloom", content="first")
    bucket = FakeBucket()
    remotes = []
    b2.push_file(bucket, first, remotes)
    second = save(tmp_path, "mantle-bloom-seed1-1y.bloom", content="other")
    b2.push_file(bucket, second, remotes)
    assert len(bucket.files) == 2


def test_push_skips_an_exact_duplicate_under_another_name(tmp_path):
    first = save(tmp_path, "mantle-bloom-seed1-1y.bloom", content="same")
    copy = save(tmp_path, "mantle-bloom-seed1-1y(1).bloom", content="same")
    bucket = FakeBucket()
    remotes = []
    b2.push_file(bucket, first, remotes)
    b2.push_file(bucket, copy, remotes)
    assert list(bucket.files) == [first.name]


def test_pull_rejects_a_corrupt_download_and_leaves_nothing(tmp_path):
    source = save(tmp_path, "mantle-bloom-seed1-1y.bloom")
    bucket = FakeBucket()
    b2.push_file(bucket, source, [])
    dest_dir = tmp_path / "elsewhere"
    bucket.corrupt_downloads = True
    with pytest.raises(ValueError, match="SHA-1"):
        b2.pull(bucket, source.name, dest_dir)
    assert list(dest_dir.iterdir()) == []


def test_pull_never_overwrites_different_contents(tmp_path):
    source = save(tmp_path, "mantle-bloom-seed1-1y.bloom", content="cloud")
    bucket = FakeBucket()
    b2.push_file(bucket, source, [])
    save(tmp_path, source.name, content="local")
    with pytest.raises(ValueError, match="different contents"):
        b2.pull(bucket, source.name, tmp_path)
    assert source.read_text() == "local"


def test_pull_keeps_a_bucket_name_inside_the_saves_directory(tmp_path):
    source = save(tmp_path, "mantle-bloom-seed1-1y.bloom")
    bucket = FakeBucket()
    bucket.upload(source, "../../escape.bloom", b2.sha1_of(source), {})
    dest_dir = tmp_path / "saves"
    assert b2.pull(bucket, "../../escape.bloom", dest_dir) == dest_dir / "escape.bloom"


def test_metadata_reads_name_and_manifest(tmp_path):
    path = tmp_path / "mantle-bloom-seed42-1500000y.bloom"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("manifest.json", json.dumps({"version": 5, "name": "Gliese b"}))
    info = b2.save_metadata(path, note="n" * 1000, machine="studio")
    assert info["seed"] == "42"
    assert info["sim_years"] == "1500000"
    assert info["format_version"] == "5"
    assert info["world_name"] == "Gliese b"
    assert info["machine"] == "studio"
    assert len(info["note"]) == b2.MAX_NOTE_CHARS


def test_metadata_tolerates_a_pickled_save(tmp_path):
    path = save(tmp_path, "mantle-bloom-seed42-1500000y.mbworld")
    info = b2.save_metadata(path, machine="studio")
    assert info["kind"] == "mbworld"
    assert "format_version" not in info


def test_catalog_marks_where_each_file_lives():
    remotes = [
        b2.Remote("both.bloom", 3, "a" * 40, 1000, {"machine": "studio"}),
        b2.Remote("cloud.bloom", 3, "b" * 40, 2000),
    ]

    class Stat:
        st_size = 3

    class Local:
        def __init__(self, name):
            self.name = name

        def stat(self):
            return Stat

    rows = b2.catalog(remotes, [Local("both.bloom"), Local("local.bloom")], {"both.bloom"})
    where = {row["name"]: row["where"] for row in rows}
    assert where == {"both.bloom": "both", "cloud.bloom": "cloud", "local.bloom": "local"}
    assert next(row for row in rows if row["name"] == "both.bloom")["pinned"]


def test_config_reads_the_env_file_and_the_environment_wins(tmp_path):
    env = tmp_path / "b2.env"
    env.write_text("# comment\nB2_APPLICATION_KEY_ID=fromfile\nB2_BUCKET=bucket\n")
    cfg = b2.load_config({"B2_BUCKET": "fromenv", "MANTLE_BLOOM_SAVES_DIR": "/x"}, env)
    assert (cfg.key_id, cfg.bucket, cfg.saves_dir) == ("fromfile", "fromenv", Path("/x"))
    assert b2.missing_settings(cfg) == ["B2_APPLICATION_KEY"]


def test_cron_line_replacement_keeps_other_jobs():
    existing = "0 1 * * * backup\n" + b2.cron_line("/py", Path("/s.py"), Path("/l"))
    assert b2.without_cron_line(existing) == "0 1 * * * backup\n"


class Stub(BaseHTTPRequestHandler):
    """Just enough of B2's native API, recording every request."""

    requests: list = []
    stored: dict = {}
    fail_next_upload = False

    def log_message(self, *_args):
        pass

    def reply(self, body, status=200):
        data = body if isinstance(body, bytes) else json.dumps(body).encode()
        self.send_response(status)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def base(self):
        return f"http://127.0.0.1:{self.server.server_address[1]}"

    def do_GET(self):
        Stub.requests.append((self.path, dict(self.headers), b""))
        if self.path.endswith("b2_authorize_account"):
            self.reply(
                {
                    "accountId": "acct",
                    "authorizationToken": "token",
                    "apiInfo": {"storageApi": {"apiUrl": self.base(), "downloadUrl": self.base()}},
                }
            )
        else:
            self.reply(Stub.stored.get(self.path, b""))

    def do_POST(self):
        body = self.rfile.read(int(self.headers["Content-Length"]))
        Stub.requests.append((self.path, dict(self.headers), body))
        if self.path.endswith("b2_list_buckets"):
            self.reply({"buckets": [{"bucketId": "bid", "bucketType": "allPrivate"}]})
        elif self.path.endswith("b2_get_upload_url"):
            self.reply({"uploadUrl": f"{self.base()}/upload", "authorizationToken": "up"})
        elif self.path.endswith("b2_list_file_names"):
            self.reply({"files": [], "nextFileName": None})
        elif self.path == "/upload":
            if Stub.fail_next_upload:
                Stub.fail_next_upload = False
                self.reply({"code": "service_unavailable", "message": "busy"}, 503)
                return
            self.reply(
                {
                    "fileName": self.headers["X-Bz-File-Name"].replace("%28", "("),
                    "contentLength": len(body),
                    "contentSha1": self.headers["X-Bz-Content-Sha1"],
                    "uploadTimestamp": 1000,
                    "fileInfo": {"machine": self.headers["X-Bz-Info-machine"]},
                }
            )


@pytest.fixture
def client(monkeypatch):
    server = ThreadingHTTPServer(("127.0.0.1", 0), Stub)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    Stub.requests, Stub.stored, Stub.fail_next_upload = [], {}, False
    monkeypatch.setattr(b2.time, "sleep", lambda _seconds: None)
    yield b2.B2("id", "secret", "bucket", api_url=f"http://127.0.0.1:{server.server_address[1]}")
    server.shutdown()
    server.server_close()


def test_wire_upload_sends_the_headers_b2_requires(client, tmp_path):
    path = save(tmp_path, "mantle-bloom-seed1-1y(1).bloom", content="bytes")
    sha1 = b2.sha1_of(path)

    remote = client.upload(path, path.name, sha1, {"machine": "my studio"})

    assert remote.sha1 == sha1 and remote.size == 5
    url, headers, body = Stub.requests[-1]
    assert (url, body) == ("/upload", b"bytes")
    assert headers["Authorization"] == "up"
    assert headers["X-Bz-File-Name"] == "saves/mantle-bloom-seed1-1y%281%29.bloom"
    assert headers["X-Bz-Content-Sha1"] == sha1
    assert headers["Content-Length"] == "5"
    assert headers["X-Bz-Info-Machine"] == "my%20studio"
    assert "Transfer-Encoding" not in headers


def test_wire_upload_gets_a_new_url_after_a_503(client, tmp_path):
    path = save(tmp_path, "mantle-bloom-seed1-1y.bloom")
    Stub.fail_next_upload = True
    client.upload(path, path.name, b2.sha1_of(path), {"machine": "studio"})
    urls = [request[0] for request in Stub.requests]
    assert sum(url.endswith("b2_get_upload_url") for url in urls) == 2


def test_wire_authorizes_with_basic_auth_on_api_v4_and_finds_the_bucket(client):
    assert client.type_of_bucket() == "allPrivate"
    url, headers, _ = Stub.requests[0]
    assert url == "/b2api/v4/b2_authorize_account"
    assert headers["Authorization"] == "Basic aWQ6c2VjcmV0"  # id:secret
    assert json.loads(Stub.requests[1][2]) == {"accountId": "acct", "bucketName": "bucket"}


def test_wire_download_authorizes_and_hashes_what_arrived(client, tmp_path):
    Stub.stored["/file/bucket/saves/a%20b.bloom"] = b"payload"
    dest = tmp_path / "out"
    assert client.download("a b.bloom", dest) == hashlib.sha1(b"payload").hexdigest()
    assert dest.read_bytes() == b"payload"
    assert Stub.requests[-1][1]["Authorization"] == "token"
