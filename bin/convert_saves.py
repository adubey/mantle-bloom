#!/usr/bin/env python3
"""Convert legacy `.mbworld` saves to `.bloom`, keeping every old hash resolvable.

Most files named `.mbworld` are already v5 zips (the extension is all that is legacy). Those
are renamed: same bytes, same SHA-1, so nothing needs recording. The rest are pickles
(format 1-4); each is loaded, written as a v5 `.bloom`, reloaded, and accepted only if the
reloaded world has the same `surface_parity.state_hash`. A converted file has new bytes, so
the old name and hashes are recorded where lookups can find them:

* locally, in `LEGACY_HASHES.json` beside the file (legacy SHA-1 -> new name and hashes);
* in the bucket, as `legacy_name` and `legacy_sha1` file info (plus `legacy_sha256` for a
  converted pickle; a server-side copy never reads the bytes, and the bucket stores only
  SHA-1), which `bin/b2_saves.py` resolves in `pull`, `sync`, `purge` and `list`.

Originals are never deleted. Run with the backend venv (it imports `backend/app`):

    backend/.venv/bin/python -I bin/convert_saves.py local DIR [DIR ...] [--dry-run]
    backend/.venv/bin/python -I bin/convert_saves.py remote [--name NAME] [--dry-run]
    backend/.venv/bin/python -I bin/convert_saves.py recompress DIR [DIR ...] [--dry-run]
    backend/.venv/bin/python -I bin/convert_saves.py upload DIR [DIR ...] [--dry-run]
    backend/.venv/bin/python -I bin/convert_saves.py zips DIR [DIR ...] [--dry-run] [--verify-step]

`upload` sends one local `.bloom` per run, the next whose original `.mbworld` is still the only
copy in the bucket, with `legacy_name`/`legacy_sha1`/`legacy_sha256` file info. It finds the
`.bloom` through the `LEGACY_HASHES.json` files under DIR, so nothing is downloaded and the
converted bytes are exactly what was verified locally. The bucket's `.mbworld` stays.

`zips` converts each `*.mbworld.zip` (a zip wrapping one `.mbworld`) to a deflated `.bloom`
named after the inner file. The inner file's SHA-1 and SHA-256, which are what the bucket holds,
and the zip's own are recorded in `LEGACY_HASHES.json`, with the world's true seed and years
(a file name can disagree with its world). The zip is kept.

`recompress` rewrites stored `.bloom` zips with deflate (what `save_world_bytes` now writes,
#860) in place. It copies entries without loading the world, replaces the file only after the
new zip's entry names, sizes and CRCs match the old one, and then records the old SHA-1 in
`LEGACY_HASHES.json` and corrects the file's `SHA256SUMS` line.

`remote` converts exactly one bucket file per run (the next unconverted one, or `--name`),
so a shell loop controls the download/upload budget. A bucket file that is already a v5 zip is
renamed with a server-side copy, which moves no bytes through this machine.
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "backend"))

LEGACY_FILE = "LEGACY_HASHES.json"
SUMS_FILE = "SHA256SUMS"


def _load_b2():
    spec = importlib.util.spec_from_file_location("b2_saves", REPO / "bin" / "b2_saves.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["b2_saves"] = module
    spec.loader.exec_module(module)
    return module


def digests(data: bytes) -> tuple[str, str]:
    return hashlib.sha1(data).hexdigest(), hashlib.sha256(data).hexdigest()


def is_zip(data: bytes) -> bool:
    return data[:2] == b"PK"


def file_is_zip(path: Path) -> bool:
    with path.open("rb") as handle:
        return is_zip(handle.read(2))


def new_name(name: str) -> str:
    return str(Path(name).with_suffix(".bloom")) if name.endswith(".mbworld") else name


def backfill_dataclass_defaults(world) -> list[str]:
    """Give cached dataclass objects the fields added since the pickle was written.

    A pickle restores `__dict__` only, so a field added later (`HydrologyFields.lake_ice_m`)
    is simply absent and the v5 writer fails on it. Fill it from the field's own default,
    which is what a freshly built object holds. Returns the filled field names."""
    filled: list[str] = []
    seen: set[int] = set()

    def visit(value) -> None:
        if id(value) in seen:
            return
        if isinstance(value, dict):
            seen.add(id(value))
            for child in list(value.values()):
                visit(child)
        elif isinstance(value, (list, tuple, set, frozenset)):
            seen.add(id(value))
            for child in list(value):
                visit(child)
        elif hasattr(value, "__dict__") and not isinstance(value, type):
            # Any object, not just a dataclass: a plate (plain class) can hold a cache
            # dataclass that gained a field.
            seen.add(id(value))
            if dataclasses.is_dataclass(value):
                fill(value)
            for child in list(vars(value).values()):
                visit(child)

    def fill(value) -> None:
        for field in dataclasses.fields(value):
            if field.name in value.__dict__ or hasattr(type(value), field.name):
                continue
            if field.default is not dataclasses.MISSING:
                setattr(value, field.name, field.default)
            elif field.default_factory is not dataclasses.MISSING:
                setattr(value, field.name, field.default_factory())
            else:
                raise ValueError(f"{type(value).__name__}.{field.name} has no default")
            filled.append(f"{type(value).__name__}.{field.name}")

    visit(world)
    return filled


def convert_pickle_bytes(data: bytes, *, verify_step: bool = False) -> tuple[bytes, dict]:
    """The v5 bytes for a pickled save, and what was checked. Raises if the round trip moves state."""
    from app import persistence, surface_parity, world as world_mod

    world = persistence.load_world_bytes(data)
    name = persistence.load_world_name(data)
    filled = backfill_dataclass_defaults(world)
    bloom = persistence.save_world_bytes(world, name=name)
    again = persistence.load_world_bytes(bloom)
    before, after = surface_parity.state_hash(world), surface_parity.state_hash(again)
    if before != after:
        raise ValueError("reloaded .bloom differs from the pickle-loaded world")
    report = {"state_hash": before, "backfilled": sorted(set(filled))}
    if verify_step:
        # One step from each side: the conversion must not change how the world evolves.
        world_mod.step_world(world, 100_000.0)
        world_mod.step_world(again, 100_000.0)
        if surface_parity.state_hash(world) != surface_parity.state_hash(again):
            raise ValueError("a step from the .bloom lands on different state than from the pickle")
        report["step_checked"] = True
    return bloom, report


def read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def write_json(path: Path, value: dict) -> None:
    partial = path.with_name(path.name + ".part")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(partial, path)


def record_legacy(directory: Path, legacy: dict, current: dict) -> None:
    table = read_json(directory / LEGACY_FILE)
    table[legacy["sha1"]] = {
        "legacy_name": legacy["name"],
        "legacy_sha256": legacy["sha256"],
        "name": current["name"],
        "sha1": current["sha1"],
        "sha256": current["sha256"],
        "how": current["how"],
    }
    write_json(directory / LEGACY_FILE, table)


def update_sums(directory: Path, old: str, new: str, new_sha256: str | None, keep_old: bool) -> None:
    """Keep a directory's SHA256SUMS true: rename a line (same bytes), or add one (new bytes)."""
    path = directory / SUMS_FILE
    if not path.exists():
        return
    lines, out = path.read_text().splitlines(), []
    for line in lines:
        digest, _, name = line.partition("  ")
        if name.lstrip("*") == old:
            if keep_old:
                out.append(line)
                out.append(f"{new_sha256}  {new}")
            else:
                out.append(f"{digest}  {new}")
        else:
            out.append(line)
    path.write_text("\n".join(out) + "\n")


def convert_local(path: Path, *, dry_run: bool, verify_step: bool) -> str:
    data = path.read_bytes()
    sha1, sha256 = digests(data)
    target = path.with_name(new_name(path.name))
    legacy = {"name": path.name, "sha1": sha1, "sha256": sha256}
    if target.exists():
        return f"SKIP {path}: {target.name} already exists"
    if is_zip(data):
        if not dry_run:
            os.rename(path, target)
            update_sums(path.parent, path.name, target.name, None, keep_old=False)
            record_legacy(path.parent, legacy, {"name": target.name, "sha1": sha1, "sha256": sha256, "how": "renamed"})
        return f"renamed {path.name} -> {target.name}"
    if dry_run:
        return f"would convert pickle {path.name} ({len(data) / 1e6:.0f} MB)"
    bloom, report = convert_pickle_bytes(data, verify_step=verify_step)
    new_sha1, new_sha256 = digests(bloom)
    partial = target.with_name(target.name + ".part")
    partial.write_bytes(bloom)
    os.replace(partial, target)
    record_legacy(path.parent, legacy, {"name": target.name, "sha1": new_sha1, "sha256": new_sha256, "how": "converted"})
    update_sums(path.parent, path.name, target.name, new_sha256, keep_old=True)
    return f"converted {path.name} -> {target.name} ({len(bloom) / 1e6:.0f} MB; {report}); original kept"


def set_sum(directory: Path, name: str, sha256: str) -> None:
    """Point a directory's SHA256SUMS line for `name` at new bytes."""
    path = directory / SUMS_FILE
    if not path.exists():
        return
    out = []
    for line in path.read_text().splitlines():
        _digest, _, listed = line.partition("  ")
        out.append(f"{sha256}  {listed}" if listed.lstrip("*") == name else line)
    path.write_text("\n".join(out) + "\n")


def file_digests(path: Path) -> tuple[str, str]:
    sha1, sha256 = hashlib.sha1(), hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            sha1.update(chunk)
            sha256.update(chunk)
    return sha1.hexdigest(), sha256.hexdigest()


def is_stored(path: Path) -> bool:
    with zipfile.ZipFile(path) as archive:
        return any(info.compress_type == zipfile.ZIP_STORED for info in archive.infolist())


def recompress_zip(src: Path, dst: Path) -> None:
    """Copy `src` to `dst` with every entry deflated, one entry in memory at a time, and
    raise unless `dst` holds the same names, sizes and CRCs."""
    from app import npz_save

    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst, "w", zipfile.ZIP_DEFLATED, compresslevel=npz_save.COMPRESSLEVEL) as zout:
        for info in zin.infolist():
            zout.writestr(info.filename, zin.read(info.filename))
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(dst) as zout:
        before = [(i.filename, i.file_size, i.CRC) for i in zin.infolist()]
        after = [(i.filename, i.file_size, i.CRC) for i in zout.infolist()]
        if before != after or zout.testzip() is not None:
            raise ValueError("recompressed zip differs from the original")


def record_recompressed(directory: Path, name: str, old: tuple[str, str], new: tuple[str, str]) -> None:
    """Keep the old SHA-1 resolvable: refresh every entry that pointed at the old bytes, and
    add one for the old SHA-1 itself when none had it as a key."""
    table = read_json(directory / LEGACY_FILE)
    for entry in table.values():
        if entry.get("sha1") == old[0]:
            entry.update(sha1=new[0], sha256=new[1], how=entry["how"] + "+recompressed")
    table.setdefault(old[0], {"legacy_name": name, "legacy_sha256": old[1], "name": name, "sha1": new[0], "sha256": new[1], "how": "recompressed"})
    write_json(directory / LEGACY_FILE, table)


def recompress_local(path: Path, *, dry_run: bool) -> str:
    if not is_stored(path):
        return f"SKIP {path.name}: already deflated"
    if dry_run:
        return f"would recompress {path.name} ({path.stat().st_size / 1e6:.0f} MB)"
    old = file_digests(path)
    before = path.stat().st_size
    partial = path.with_name(path.name + ".part")
    try:
        recompress_zip(path, partial)
        new = file_digests(partial)
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)
    record_recompressed(path.parent, path.name, old, new)
    set_sum(path.parent, path.name, new[1])
    return f"recompressed {path.name}: {before / 1e6:.0f} -> {path.stat().st_size / 1e6:.0f} MB"


def convert_zip(path: Path, *, dry_run: bool, verify_step: bool) -> str:
    with zipfile.ZipFile(path) as outer:
        (info,) = outer.infolist()
        inner_name = Path(info.filename).name
        data = outer.read(info.filename)
    target = path.parent / new_name(inner_name)
    if target.exists():
        return f"SKIP {path.name}: {target.name} already exists"
    if dry_run:
        return f"would convert {path.name} ({'v5' if is_zip(data) else 'pickle'}, {len(data) / 1e6:.0f} MB) -> {target.name}"
    inner_sha1, inner_sha256 = digests(data)
    partial = target.with_name(target.name + ".part")
    try:
        if is_zip(data):
            staged = target.with_name(target.name + ".stored")
            staged.write_bytes(data)
            try:
                recompress_zip(staged, partial)
            finally:
                staged.unlink(missing_ok=True)
            how, report = "unzipped+recompressed", {}
        else:
            bloom, report = convert_pickle_bytes(data, verify_step=verify_step)
            partial.write_bytes(bloom)
            how = "unzipped+converted"
        os.replace(partial, target)
    finally:
        partial.unlink(missing_ok=True)
    new_sha1, new_sha256 = file_digests(target)
    zip_sha1, zip_sha256 = file_digests(path)
    table = read_json(path.parent / LEGACY_FILE)
    table[inner_sha1] = {
        "legacy_name": inner_name, "legacy_sha256": inner_sha256, "legacy_zip": path.name,
        "legacy_zip_sha1": zip_sha1, "legacy_zip_sha256": zip_sha256,
        "name": target.name, "sha1": new_sha1, "sha256": new_sha256, "how": how,
        "state_hash": report.get("state_hash"),
    }
    write_json(path.parent / LEGACY_FILE, table)
    return f"{path.name} -> {target.name} ({len(data) / 1e6:.0f} -> {target.stat().st_size / 1e6:.0f} MB; {how})"


def cmd_zips(args: argparse.Namespace) -> int:
    paths = sorted(p for d in args.dirs for p in Path(d).expanduser().glob("*.mbworld.zip"))
    failed = 0
    for path in paths:
        try:
            print(convert_zip(path, dry_run=args.dry_run, verify_step=args.verify_step), flush=True)
        except (ValueError, OSError, TypeError, AttributeError, zipfile.BadZipFile) as err:
            failed += 1
            print(f"FAILED {path}: {err}", file=sys.stderr)
    print(f"{len(paths)} files, {failed} failed")
    return 1 if failed else 0


def cmd_recompress(args: argparse.Namespace) -> int:
    paths = sorted(p for d in args.dirs for p in Path(d).expanduser().rglob("*.bloom"))
    failed = 0
    for path in paths:
        try:
            print(recompress_local(path, dry_run=args.dry_run), flush=True)
        except (ValueError, OSError, zipfile.BadZipFile) as err:
            failed += 1
            print(f"FAILED {path}: {err}", file=sys.stderr)
    print(f"{len(paths)} files, {failed} failed")
    return 1 if failed else 0


def cmd_local(args: argparse.Namespace) -> int:
    paths = sorted(p for d in args.dirs for p in Path(d).expanduser().rglob("*.mbworld"))
    # Zips first: they are free, and pickles are the memory-heavy ones, done one at a time.
    paths.sort(key=lambda p: not file_is_zip(p))
    failed = 0
    for path in paths:
        try:
            print(convert_local(path, dry_run=args.dry_run, verify_step=args.verify_step), flush=True)
        except (ValueError, OSError, TypeError, AttributeError) as err:
            failed += 1
            print(f"FAILED {path}: {err}", file=sys.stderr)
    print(f"{len(paths)} files, {failed} failed")
    return 1 if failed else 0


def next_remote(remotes: list, name: str | None):
    by_name = {r.name: r for r in remotes}
    converted = {r.info.get("legacy_sha1") for r in remotes} | {r.sha1 for r in remotes if r.name.endswith(".bloom")}
    candidates = [r for r in remotes if r.name.endswith(".mbworld") and r.sha1 not in converted and new_name(r.name) not in by_name]
    if name:
        candidates = [r for r in candidates if r.name == name]
    return candidates[0] if candidates else None


def local_blooms_for_bucket(dirs: list[str], remotes: list) -> list[tuple]:
    """(remote .mbworld, local .bloom, legacy entry) for every bucket `.mbworld` whose
    SHA-1 some `LEGACY_HASHES.json` under `dirs` maps to an existing local file and that
    has no converted copy in the bucket yet."""
    converted = {r.info.get("legacy_sha1") for r in remotes}
    by_sha1 = {r.sha1: r for r in remotes if r.name.endswith(".mbworld") and r.sha1 not in converted}
    found: dict[str, tuple] = {}
    for directory in dirs:
        for table_path in sorted(Path(directory).expanduser().rglob(LEGACY_FILE)):
            for key, entry in read_json(table_path).items():
                local = table_path.parent / entry["name"]
                if key in by_sha1 and key not in found and local.exists():
                    found[key] = (by_sha1[key], local, entry)
    return sorted(found.values(), key=lambda item: item[0].name)


def cmd_upload(args: argparse.Namespace) -> int:
    b2 = _load_b2()
    client = b2.client_for(b2.load_config())
    remotes = client.list_files()
    pending = local_blooms_for_bucket(args.dirs, remotes)
    print(f"{len(pending)} bucket .mbworld files have a local .bloom and no converted copy")
    if not pending:
        return 0
    remote, local, entry = pending[0]
    if local.stat().st_size > b2.MAX_UPLOAD_BYTES:
        raise ValueError(f"{local.name} is over the upload limit")
    name = b2.unique_remote_name(new_name(remote.name), remote.sha1, {r.name for r in remotes})
    if args.dry_run:
        print(f"would upload {local} as {name} ({local.stat().st_size / 1e6:.0f} MB)")
        return 0
    sha1, sha256 = file_digests(local)
    if entry["sha1"] != sha1:
        raise ValueError(f"{local} no longer matches its LEGACY_HASHES.json entry; rerun the conversion")
    info = {
        **b2.save_metadata(local),
        **{k: v for k, v in remote.info.items() if k not in {"kind", "format_version", "world_name"}},
        "kind": "bloom",
        "format_version": "5",
        "legacy_name": remote.name,
        "legacy_sha1": remote.sha1,
    }
    if entry.get("legacy_sha256"):
        info["legacy_sha256"] = entry["legacy_sha256"]
    if entry.get("state_hash"):
        info["converted_state_hash"] = entry["state_hash"]
    uploaded = client.upload(local, name, sha1, info)
    if uploaded.sha1 != sha1:
        raise ValueError(f"{name}: bucket reports SHA-1 {uploaded.sha1}, expected {sha1}")
    print(f"uploaded {name} ({uploaded.size / 1e6:.0f} MB) for {remote.name}; {len(pending) - 1} left")
    return 0


def cmd_remote(args: argparse.Namespace) -> int:
    b2 = _load_b2()
    cfg = b2.load_config()
    client = b2.client_for(cfg)
    remotes = client.list_files()
    remote = next_remote(remotes, args.name)
    if remote is None:
        print("nothing left to convert")
        return 0
    zipped = remote.info.get("format_version") == "5"
    taken = {r.name for r in remotes}
    name = b2.unique_remote_name(new_name(remote.name), remote.sha1, taken)
    if zipped:
        # Already v5: a server-side copy gives it the new name with identical bytes.
        print(f"{'would copy' if args.dry_run else 'copying'} {remote.name} -> {name} server-side ({remote.size / 1e6:.0f} MB, no transfer)")
        if not args.dry_run:
            copied = client.copy(remote, name, {**remote.info, "kind": "bloom", "legacy_name": remote.name, "legacy_sha1": remote.sha1})
            assert copied.sha1 == remote.sha1
        return 0
    print(f"{'would convert' if args.dry_run else 'converting'} {remote.name} ({remote.size / 1e6:.0f} MB down, about as much up)")
    if args.dry_run:
        return 0
    with tempfile.TemporaryDirectory(dir=args.work_dir) as work:
        source = Path(work) / remote.name
        got = client.download(remote.name, source)
        if got != remote.sha1:
            raise ValueError(f"{remote.name}: downloaded SHA-1 {got}, expected {remote.sha1}")
        data = source.read_bytes()
        if is_zip(data):
            raise ValueError(f"{remote.name} is a zip but lacks format_version; rerun after fixing its info")
        bloom, report = convert_pickle_bytes(data, verify_step=args.verify_step)
        source.unlink()
        target = Path(work) / name
        target.write_bytes(bloom)
        sha1, sha256 = digests(bloom)
        info = {
            **b2.save_metadata(target),
            **{k: v for k, v in remote.info.items() if k not in {"kind", "format_version", "world_name"}},
            "kind": "bloom",
            "format_version": "5",
            "legacy_name": remote.name,
            "legacy_sha1": remote.sha1,
            "legacy_sha256": hashlib.sha256(data).hexdigest(),
            "converted_state_hash": report["state_hash"],
        }
        uploaded = client.upload(target, name, sha1, info)
        if uploaded.sha1 != sha1:
            raise ValueError(f"{name}: bucket reports SHA-1 {uploaded.sha1}, expected {sha1}")
    print(f"uploaded {name} ({uploaded.size / 1e6:.0f} MB); {report}; the original {remote.name} stays in the bucket")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    local = commands.add_parser("local")
    local.add_argument("dirs", nargs="+")
    local.add_argument("--dry-run", action="store_true")
    local.add_argument("--verify-step", action="store_true", help="also step both worlds once and compare")
    local.set_defaults(handler=cmd_local)
    recompress = commands.add_parser("recompress")
    recompress.add_argument("dirs", nargs="+")
    recompress.add_argument("--dry-run", action="store_true")
    recompress.set_defaults(handler=cmd_recompress)
    zips = commands.add_parser("zips")
    zips.add_argument("dirs", nargs="+")
    zips.add_argument("--dry-run", action="store_true")
    zips.add_argument("--verify-step", action="store_true")
    zips.set_defaults(handler=cmd_zips)
    upload = commands.add_parser("upload")
    upload.add_argument("dirs", nargs="+")
    upload.add_argument("--dry-run", action="store_true")
    upload.set_defaults(handler=cmd_upload)
    remote = commands.add_parser("remote")
    remote.add_argument("--name")
    remote.add_argument("--dry-run", action="store_true")
    remote.add_argument("--verify-step", action="store_true")
    remote.add_argument("--work-dir", default=None, help="scratch directory for the one download")
    remote.set_defaults(handler=cmd_remote)
    args = parser.parse_args(argv)
    return args.handler(args)


if __name__ == "__main__":
    sys.exit(main())
