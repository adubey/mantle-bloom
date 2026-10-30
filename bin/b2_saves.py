#!/usr/bin/env python3
"""Keep Mantle Bloom saves in a private Backblaze B2 bucket across machines (#708).

Recent saves live in the local saves directory and in the bucket. Older local copies are
purged once the bucket is confirmed to hold them; bucket copies stay until deleted by hand.
See .agents/skills/b2-saves/SKILL.md for the procedure.

    bin/b2_saves.py doctor                  # check credentials, bucket and cron
    bin/b2_saves.py list [--json]           # every file, local and in the bucket
    bin/b2_saves.py push FILE [--note TEXT]
    bin/b2_saves.py pull NAME
    bin/b2_saves.py sync                    # push every local save the bucket lacks
    bin/b2_saves.py purge [--dry-run]       # delete old local copies that are in the bucket
    bin/b2_saves.py pin NAME | unpin NAME   # exempt a local save from purge, on this machine
    bin/b2_saves.py delete-remote NAME --yes
    bin/b2_saves.py cron                    # sync then purge; what the cron job runs
    bin/b2_saves.py install-cron | uninstall-cron

Standard library only, so it runs on any machine without a worktree or venv. Credentials
come from ~/.config/mantle-bloom/b2.env (KEY=value lines) or the environment, never from the
repo. The local directory is MANTLE_BLOOM_SAVES_DIR (default ~/Downloads), the same place
the replay tools look, so saves are not moved anywhere.

Metadata is B2 file info, stored with each file, so one listing shows everything on any
machine and there's no shared catalog for two machines to race on.
"""

from __future__ import annotations

import argparse
import base64
import fcntl
import fnmatch
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

API_URL = "https://api.backblazeb2.com"
API_VERSION = "b2api/v4"
REMOTE_PREFIX = "saves/"
SAVE_PATTERNS = ("mantle-bloom-*.bloom", "mantle-bloom-*.mbworld")
# b2_upload_file takes one file up to 5 GB; saves are well below that, so no multipart path.
MAX_UPLOAD_BYTES = 5_000_000_000
# A browser may still be writing a file this recently changed.
SETTLE_SECONDS = 60
MAX_NOTE_CHARS = 400
RETRY_STATUSES = {401, 408, 429, 500, 503}
ATTEMPTS = 4
HASH_CHUNK = 1 << 20

CONFIG_DIR = Path.home() / ".config" / "mantle-bloom"
ENV_FILE = CONFIG_DIR / "b2.env"
PINS_FILE = CONFIG_DIR / "b2_pins.txt"
STATE_DIR = Path.home() / ".local" / "state" / "mantle-bloom"
LOCK_FILE = STATE_DIR / "b2_saves.lock"
LOG_FILE = STATE_DIR / "b2_saves.log"
# Cron must not depend on a worktree that can vanish, so it runs a private copy.
INSTALLED_SCRIPT = Path.home() / ".local" / "share" / "mantle-bloom" / "b2_saves.py"
CRON_MARKER = "# mantle-bloom-b2-saves"
# Hourly, not daily: a sleeping laptop skips cron, and a run with nothing to do is cheap.
CRON_SCHEDULE = "17 * * * *"

SEED_YEARS = re.compile(r"mantle-bloom-seed(\d+)-(\d+)y")


class B2Error(RuntimeError):
    def __init__(self, status: int, code: str, message: str):
        super().__init__(f"B2 error {status} {code}: {message}")
        self.status = status


@dataclass(frozen=True)
class Config:
    key_id: str
    key: str
    bucket: str
    saves_dir: Path
    keep_days: float
    keep_newest: int


@dataclass(frozen=True)
class Remote:
    name: str  # without REMOTE_PREFIX
    size: int
    sha1: str
    uploaded_ms: int
    info: dict[str, str] = field(default_factory=dict)


def read_env_file(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    values = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def load_config(environ=os.environ, env_file: Path = ENV_FILE) -> Config:
    values = {**read_env_file(env_file), **environ}
    return Config(
        key_id=values.get("B2_APPLICATION_KEY_ID", ""),
        key=values.get("B2_APPLICATION_KEY", ""),
        bucket=values.get("B2_BUCKET", ""),
        saves_dir=Path(values.get("MANTLE_BLOOM_SAVES_DIR") or Path.home() / "Downloads"),
        keep_days=float(values.get("B2_SAVES_KEEP_DAYS") or 14),
        keep_newest=int(values.get("B2_SAVES_KEEP_NEWEST") or 5),
    )


def missing_settings(cfg: Config) -> list[str]:
    settings = {
        "B2_APPLICATION_KEY_ID": cfg.key_id,
        "B2_APPLICATION_KEY": cfg.key,
        "B2_BUCKET": cfg.bucket,
    }
    return [name for name, value in settings.items() if not value]


def sha1_of(path: Path) -> str:
    digest = hashlib.sha1()
    with path.open("rb") as handle:
        while chunk := handle.read(HASH_CHUNK):
            digest.update(chunk)
    return digest.hexdigest()


def remote_key(name: str) -> str:
    return REMOTE_PREFIX + name


def to_remote(entry: dict) -> Remote:
    info = {k: urllib.parse.unquote(v) for k, v in (entry.get("fileInfo") or {}).items()}
    return Remote(
        name=entry["fileName"].removeprefix(REMOTE_PREFIX),
        size=entry["contentLength"],
        sha1=entry["contentSha1"],
        uploaded_ms=entry["uploadTimestamp"],
        info=info,
    )


class B2:
    """The few B2 native API calls this tool needs."""

    def __init__(self, key_id: str, key: str, bucket: str, api_url: str = API_URL):
        self.key_id, self.key, self.bucket = key_id, key, bucket
        self.api_url = api_url
        self.token = self.account_id = self.storage_url = self.download_url = None
        self.bucket_id = self.bucket_type = None

    @staticmethod
    def _open(request: urllib.request.Request) -> bytes:
        try:
            with urllib.request.urlopen(request, timeout=600) as response:
                return response.read()
        except urllib.error.HTTPError as err:
            try:
                body = json.loads(err.read())
            except ValueError:
                body = {}
            raise B2Error(err.code, body.get("code", "unknown"), body.get("message", "")) from None

    def _authorize(self) -> None:
        basic = base64.b64encode(f"{self.key_id}:{self.key}".encode()).decode()
        request = urllib.request.Request(
            f"{self.api_url}/{API_VERSION}/b2_authorize_account",
            headers={"Authorization": f"Basic {basic}"},
        )
        auth = json.loads(self._open(request))
        storage = auth["apiInfo"]["storageApi"]
        self.token = auth["authorizationToken"]
        self.account_id = auth["accountId"]
        self.storage_url = storage["apiUrl"]
        self.download_url = storage["downloadUrl"]

    def _ensure_bucket(self) -> None:
        if self.bucket_id:
            return
        # bucketName is required for a bucket-restricted key and harmless for any other.
        found = self._call(
            "b2_list_buckets",
            {"accountId": self._account(), "bucketName": self.bucket},
        )["buckets"]
        if not found:
            raise B2Error(404, "no_such_bucket", f"bucket {self.bucket!r} not found")
        self.bucket_id = found[0]["bucketId"]
        self.bucket_type = found[0]["bucketType"]

    def _account(self) -> str:
        if not self.token:
            self._authorize()
        return self.account_id

    def _call(self, endpoint: str, payload: dict) -> dict:
        self._account()
        for attempt in range(ATTEMPTS):
            request = urllib.request.Request(
                f"{self.storage_url}/{API_VERSION}/{endpoint}",
                data=json.dumps(payload).encode(),
                headers={"Authorization": self.token},
            )
            try:
                return json.loads(self._open(request))
            except B2Error as err:
                if err.status not in RETRY_STATUSES or attempt == ATTEMPTS - 1:
                    raise
                if err.status == 401:
                    self._authorize()  # tokens last a day
                else:
                    time.sleep(2**attempt)
        raise AssertionError("unreachable")

    def type_of_bucket(self) -> str:
        self._ensure_bucket()
        return self.bucket_type

    def list_files(self) -> list[Remote]:
        self._ensure_bucket()
        files, start = [], None
        while True:
            query = {"bucketId": self.bucket_id, "prefix": REMOTE_PREFIX, "maxFileCount": 1000}
            if start:
                query["startFileName"] = start
            page = self._call("b2_list_file_names", query)
            files += [
                to_remote(entry)
                for entry in page["files"]
                if entry.get("action", "upload") == "upload"
            ]
            start = page.get("nextFileName")
            if not start:
                return files

    def upload(self, path: Path, name: str, sha1: str, info: dict[str, str]) -> Remote:
        self._ensure_bucket()
        size = path.stat().st_size
        for attempt in range(ATTEMPTS):
            target = self._call("b2_get_upload_url", {"bucketId": self.bucket_id})
            headers = {
                "Authorization": target["authorizationToken"],
                "X-Bz-File-Name": urllib.parse.quote(remote_key(name), safe="/"),
                "Content-Type": "application/octet-stream",
                "Content-Length": str(size),  # B2 refuses chunked uploads
                "X-Bz-Content-Sha1": sha1,
            }
            for key, value in info.items():
                headers[f"X-Bz-Info-{key}"] = urllib.parse.quote(value, safe="")
            try:
                with path.open("rb") as body:
                    request = urllib.request.Request(target["uploadUrl"], body, headers)
                    return to_remote(json.loads(self._open(request)))
            except (B2Error, OSError) as err:
                retriable = not isinstance(err, B2Error) or err.status in RETRY_STATUSES
                if not retriable or attempt == ATTEMPTS - 1:
                    raise
                time.sleep(2**attempt)
        raise AssertionError("unreachable")

    def download(self, name: str, dest: Path) -> str:
        """Stream a file to `dest` and return the SHA-1 of what arrived."""
        self._account()
        url = "/".join(
            [
                self.download_url,
                "file",
                urllib.parse.quote(self.bucket, safe=""),
                urllib.parse.quote(remote_key(name), safe="/"),
            ]
        )
        request = urllib.request.Request(url, headers={"Authorization": self.token})
        digest = hashlib.sha1()
        try:
            with urllib.request.urlopen(request, timeout=600) as response, dest.open("wb") as out:
                while chunk := response.read(HASH_CHUNK):
                    digest.update(chunk)
                    out.write(chunk)
        except urllib.error.HTTPError as err:
            raise B2Error(err.code, "download_failed", name) from None
        return digest.hexdigest()

    def delete(self, name: str) -> int:
        """Delete every version of a file and return how many went."""
        self._ensure_bucket()
        key = remote_key(name)
        versions = self._call(
            "b2_list_file_versions",
            {"bucketId": self.bucket_id, "startFileName": key, "maxFileCount": 100},
        )["files"]
        doomed = [version for version in versions if version["fileName"] == key]
        for version in doomed:
            self._call(
                "b2_delete_file_version",
                {"fileName": key, "fileId": version["fileId"]},
            )
        return len(doomed)


def save_metadata(path: Path, note: str = "", machine: str | None = None) -> dict[str, str]:
    """What the bucket records about a save, read from its name and zip manifest."""
    stat = path.stat()
    info = {
        "src_last_modified_millis": str(int(stat.st_mtime * 1000)),
        "machine": machine or socket.gethostname(),
        "kind": path.suffix.lstrip(".") or "file",
    }
    if match := SEED_YEARS.search(path.name):
        info["seed"], info["sim_years"] = match.groups()
    try:
        with zipfile.ZipFile(path) as archive:
            manifest = json.loads(archive.read("manifest.json"))
        info["format_version"] = str(manifest.get("version", ""))
        info["world_name"] = str(manifest.get("name", ""))
    except (zipfile.BadZipFile, KeyError, ValueError):
        pass  # pickled saves before format 5, and files that aren't saves
    if note:
        info["note"] = note[:MAX_NOTE_CHARS]
    return {key: value for key, value in info.items() if value}


def local_saves(saves_dir: Path) -> list[Path]:
    # listdir, not glob: glob hides a PermissionError, so cron blocked from ~/Downloads by
    # macOS privacy controls would see no saves and report success forever.
    try:
        names = os.listdir(saves_dir)
    except FileNotFoundError:
        return []
    return sorted(
        saves_dir / name
        for name in names
        if any(fnmatch.fnmatchcase(name, pattern) for pattern in SAVE_PATTERNS)
    )


def read_pins(path: Path = PINS_FILE) -> set[str]:
    return set(path.read_text().split("\n")) - {""} if path.exists() else set()


def write_pins(pins: set[str], path: Path = PINS_FILE) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{name}\n" for name in sorted(pins)))


def unique_remote_name(name: str, sha1: str, taken: set[str]) -> str:
    """Keep two different saves that share a name from becoming versions of one file."""
    if name not in taken:
        return name
    path = Path(name)
    return f"{path.stem}-{sha1[:8]}{path.suffix}"


def push_file(client, path: Path, remotes: list[Remote], note: str = "") -> None:
    size = path.stat().st_size
    if size > MAX_UPLOAD_BYTES:
        raise ValueError(f"{path.name} is {size} bytes, over the {MAX_UPLOAD_BYTES} upload limit")
    sha1 = sha1_of(path)
    if existing := next((r for r in remotes if r.sha1 == sha1), None):
        print(f"already in bucket: {path.name} as {existing.name}")
        return
    name = unique_remote_name(path.name, sha1, {r.name for r in remotes})
    remote = client.upload(path, name, sha1, save_metadata(path, note))
    if remote.sha1 != sha1:
        raise ValueError(f"{path.name}: bucket reports SHA-1 {remote.sha1}, expected {sha1}")
    remotes.append(remote)
    print(f"pushed {path.name} as {name} ({size / 1e6:.0f} MB)")


def sync(client, saves_dir: Path, now: float) -> int:
    """Push every settled local save the bucket lacks; returns how many failed."""
    remotes = client.list_files()
    # Name and size spare hashing every save each hour. A false match only leaves that save
    # unpushed, and purge then keeps it, because purge compares SHA-1s.
    known = {(r.name, r.size) for r in remotes}
    failures = 0
    for path in local_saves(saves_dir):
        stat = path.stat()
        if (path.name, stat.st_size) in known or now - stat.st_mtime < SETTLE_SECONDS:
            continue
        try:
            push_file(client, path, remotes)
        except (B2Error, OSError, ValueError) as err:
            failures += 1
            print(f"FAILED {path.name}: {err}", file=sys.stderr)
    return failures


def purge_candidates(
    paths: list[Path], pins: set[str], now: float, keep_days: float, keep_newest: int
) -> list[Path]:
    """Saves old enough to purge, before checking the bucket holds them."""
    newest_first = sorted(paths, key=lambda path: path.stat().st_mtime, reverse=True)
    cutoff = now - keep_days * 86400
    return [
        path
        for path in newest_first[keep_newest:]
        if path.stat().st_mtime < cutoff and path.name not in pins
    ]


def purge(client, cfg: Config, now: float, dry_run: bool, pins: set[str]) -> list[Path]:
    """Delete old local saves, each only if a fresh listing shows the same bytes in the bucket."""
    candidates = purge_candidates(
        local_saves(cfg.saves_dir), pins, now, cfg.keep_days, cfg.keep_newest
    )
    if not candidates:
        return []
    in_bucket = {remote.sha1 for remote in client.list_files()}
    purged = []
    for path in candidates:
        if sha1_of(path) not in in_bucket:
            print(f"keeping {path.name}: not in bucket")
            continue
        print(f"{'would purge' if dry_run else 'purged'} {path.name}")
        if not dry_run:
            path.unlink()
        purged.append(path)
    return purged


def pull(client, name: str, saves_dir: Path) -> Path:
    remotes = {remote.name: remote for remote in client.list_files()}
    if name not in remotes:
        raise ValueError(f"{name} is not in the bucket")
    remote = remotes[name]
    dest = saves_dir / Path(name).name  # a bucket name must not escape the saves directory
    if dest.exists():
        if sha1_of(dest) == remote.sha1:
            print(f"already here: {dest}")
            return dest
        raise ValueError(f"{dest} exists with different contents; move it first")
    saves_dir.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".part")
    try:
        sha1 = client.download(name, partial)
        if sha1 != remote.sha1:
            raise ValueError(f"{name}: downloaded SHA-1 {sha1}, expected {remote.sha1}")
        os.replace(partial, dest)
    finally:
        partial.unlink(missing_ok=True)
    print(f"pulled {dest}")
    return dest


def catalog(remotes: list[Remote], paths: list[Path], pins: set[str]) -> list[dict]:
    """One row per file, joining the bucket's metadata with what this machine holds."""
    local = {path.name: path.stat().st_size for path in paths}
    rows = []
    for remote in remotes:
        here = local.get(remote.name) == remote.size
        rows.append(
            {
                "name": remote.name,
                "where": "both" if here else "cloud",
                "size": remote.size,
                "sha1": remote.sha1,
                "uploaded": datetime.fromtimestamp(
                    remote.uploaded_ms / 1000, timezone.utc
                ).isoformat(timespec="seconds"),
                "pinned": remote.name in pins,
                **remote.info,
            }
        )
    in_rows = {row["name"] for row in rows if row["where"] == "both"}
    for name, size in sorted(local.items()):
        if name not in in_rows:
            rows.append({"name": name, "where": "local", "size": size, "pinned": name in pins})
    return sorted(rows, key=lambda row: (row.get("uploaded", ""), row["name"]))


def print_catalog(rows: list[dict]) -> None:
    for row in rows:
        myr = f"{int(row['sim_years']) / 1e6:.1f} Myr" if "sim_years" in row else "-"
        parts = [
            f"{row['where']:5}{'*' if row['pinned'] else ' '}",
            f"{row['size'] / 1e6:6.0f} MB",
            f"{myr:>12}",
            f"{row.get('uploaded', '-')[:10]:10}",
            f"{row.get('machine', '-'):20}",
            row["name"],
        ]
        extra = " | ".join(row[key] for key in ("world_name", "note") if key in row)
        print("  ".join(parts) + (f"  [{extra}]" if extra else ""))
    print(f"{len(rows)} files; * is pinned on this machine")


def client_for(cfg: Config) -> B2:
    if missing := missing_settings(cfg):
        sys.exit(f"Set {', '.join(missing)} in {ENV_FILE}; see the b2-saves skill.")
    return B2(cfg.key_id, cfg.key, cfg.bucket)


def read_crontab() -> str:
    result = subprocess.run(["crontab", "-l"], capture_output=True, text=True)
    return result.stdout if result.returncode == 0 else ""  # "no crontab" exits non-zero


def write_crontab(text: str) -> None:
    subprocess.run(["crontab", "-"], input=text, text=True, check=True)


def without_cron_line(crontab: str) -> str:
    return "".join(
        line + "\n" for line in crontab.splitlines() if not line.endswith(CRON_MARKER)
    )


def cron_line(python: str, script: Path, log: Path) -> str:
    return f'{CRON_SCHEDULE} "{python}" "{script}" cron >> "{log}" 2>&1 {CRON_MARKER}\n'


def cmd_doctor(cfg: Config, _args: argparse.Namespace) -> int:
    print(f"saves directory: {cfg.saves_dir}")
    print(f"policy: keep {cfg.keep_newest} newest and anything under {cfg.keep_days:g} days")
    print(f"cron: {'installed' if CRON_MARKER in read_crontab() else 'not installed'}")
    if missing := missing_settings(cfg):
        print(f"missing in {ENV_FILE}: {', '.join(missing)}")
        return 1
    client = client_for(cfg)
    bucket_type = client.type_of_bucket()
    print(f"bucket {cfg.bucket}: {bucket_type}, {len(client.list_files())} files")
    if bucket_type != "allPrivate":
        print("WARNING: this bucket is not private, so anyone could download these saves")
        return 1
    return 0


def cmd_list(cfg: Config, args: argparse.Namespace) -> int:
    paths = local_saves(cfg.saves_dir)
    rows = catalog(client_for(cfg).list_files(), paths, read_pins())
    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        print_catalog(rows)
    return 0


def cmd_push(cfg: Config, args: argparse.Namespace) -> int:
    client = client_for(cfg)
    remotes = client.list_files()
    for file in args.files:
        push_file(client, Path(file), remotes, args.note)
    return 0


def cmd_pull(cfg: Config, args: argparse.Namespace) -> int:
    client = client_for(cfg)
    for name in args.names:
        pull(client, name, cfg.saves_dir)
    return 0


def cmd_sync(cfg: Config, _args: argparse.Namespace) -> int:
    return 1 if sync(client_for(cfg), cfg.saves_dir, time.time()) else 0


def cmd_purge(cfg: Config, args: argparse.Namespace) -> int:
    purged = purge(client_for(cfg), cfg, time.time(), args.dry_run, read_pins())
    print(f"{len(purged)} saves {'would be ' if args.dry_run else ''}purged")
    return 0


def cmd_pin(_cfg: Config, args: argparse.Namespace) -> int:
    write_pins(read_pins() | {args.name})
    return 0


def cmd_unpin(_cfg: Config, args: argparse.Namespace) -> int:
    write_pins(read_pins() - {args.name})
    return 0


def cmd_delete_remote(cfg: Config, args: argparse.Namespace) -> int:
    if not args.yes:
        sys.exit("Deleting from the bucket can't be undone; pass --yes to confirm.")
    print(f"deleted {client_for(cfg).delete(args.name)} versions of {args.name}")
    return 0


def cmd_cron(cfg: Config, _args: argparse.Namespace) -> int:
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_FILE.open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print("previous run still going")
            return 0
        print(datetime.now().isoformat(timespec="seconds"))
        client = client_for(cfg)
        now = time.time()
        # A failed push can't lose anything: purge checks the bucket for each file itself.
        failures = sync(client, cfg.saves_dir, now)
        purge(client, cfg, now, dry_run=False, pins=read_pins())
    return 1 if failures else 0


def cmd_install_cron(_cfg: Config, _args: argparse.Namespace) -> int:
    INSTALLED_SCRIPT.parent.mkdir(parents=True, exist_ok=True)
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(Path(__file__).resolve(), INSTALLED_SCRIPT)
    line = cron_line(sys.executable, INSTALLED_SCRIPT, LOG_FILE)
    write_crontab(without_cron_line(read_crontab()) + line)
    print(f"installed: {line.strip()}")
    return 0


def cmd_uninstall_cron(_cfg: Config, _args: argparse.Namespace) -> int:
    write_crontab(without_cron_line(read_crontab()))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)

    def add(name: str, handler) -> argparse.ArgumentParser:
        sub = commands.add_parser(name)
        sub.set_defaults(handler=handler)
        return sub

    add("doctor", cmd_doctor)
    add("list", cmd_list).add_argument("--json", action="store_true")
    push_parser = add("push", cmd_push)
    push_parser.add_argument("files", nargs="+")
    push_parser.add_argument("--note", default="")
    add("pull", cmd_pull).add_argument("names", nargs="+")
    add("sync", cmd_sync)
    add("purge", cmd_purge).add_argument("--dry-run", action="store_true")
    add("pin", cmd_pin).add_argument("name")
    add("unpin", cmd_unpin).add_argument("name")
    delete_parser = add("delete-remote", cmd_delete_remote)
    delete_parser.add_argument("name")
    delete_parser.add_argument("--yes", action="store_true")
    add("cron", cmd_cron)
    add("install-cron", cmd_install_cron)
    add("uninstall-cron", cmd_uninstall_cron)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(load_config(), args)
    except PermissionError as err:
        # macOS privacy controls block cron from ~/Downloads until it has Full Disk Access.
        print(f"error: {err}; on macOS, give /usr/sbin/cron Full Disk Access", file=sys.stderr)
    except (B2Error, ValueError, OSError) as err:
        print(f"error: {err}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
