"""Tests for bin/convert_saves.py (the .mbworld to .bloom move): renames keep bytes and hashes
recoverable, and the pickle path's default backfill. The pickle conversion itself loads a real
world and is exercised by hand on saved worlds."""

from __future__ import annotations

import dataclasses
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "convert_saves.py"
spec = importlib.util.spec_from_file_location("convert_saves", SCRIPT)
assert spec is not None and spec.loader is not None
cs = importlib.util.module_from_spec(spec)
sys.modules["convert_saves"] = cs
spec.loader.exec_module(cs)

ZIP = b"PK\x03\x04 stand-in for a v5 zip"


def test_a_zip_is_renamed_with_the_same_bytes_and_its_sums_line_follows(tmp_path):
    (tmp_path / "w.mbworld").write_bytes(ZIP)
    digest = hashlib.sha256(ZIP).hexdigest()
    (tmp_path / "SHA256SUMS").write_text(f"{digest}  w.mbworld\n{'0' * 64}  other.txt\n")
    message = cs.convert_local(tmp_path / "w.mbworld", dry_run=False, verify_step=False)
    assert message == "renamed w.mbworld -> w.bloom"
    assert (tmp_path / "w.bloom").read_bytes() == ZIP and not (tmp_path / "w.mbworld").exists()
    assert (tmp_path / "SHA256SUMS").read_text() == f"{digest}  w.bloom\n{'0' * 64}  other.txt\n"
    table = json.loads((tmp_path / cs.LEGACY_FILE).read_text())
    assert table[hashlib.sha1(ZIP).hexdigest()]["name"] == "w.bloom"


def test_dry_run_changes_nothing_and_an_existing_target_is_not_overwritten(tmp_path):
    (tmp_path / "w.mbworld").write_bytes(ZIP)
    assert cs.convert_local(tmp_path / "w.mbworld", dry_run=True, verify_step=False).startswith("renamed")
    assert (tmp_path / "w.mbworld").exists() and not (tmp_path / "w.bloom").exists()
    (tmp_path / "w.bloom").write_bytes(b"different")
    assert cs.convert_local(tmp_path / "w.mbworld", dry_run=False, verify_step=False).startswith("SKIP")
    assert (tmp_path / "w.bloom").read_bytes() == b"different"


def test_backfill_gives_missing_dataclass_fields_their_defaults():
    @dataclasses.dataclass
    class Cache:
        kept: int = 1
        added: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros(0))

    cache = Cache()
    del cache.__dict__["added"]  # what unpickling an older save leaves
    holder = {"cache": cache}
    assert cs.backfill_dataclass_defaults(holder) == ["Cache.added"]
    assert cache.added.shape == (0,) and cache.kept == 1


def test_next_remote_skips_files_already_converted_by_hash_or_name():
    b2 = cs._load_b2()
    remotes = [
        b2.Remote("a.mbworld", 1, "s1", 0, {}),
        b2.Remote("a.bloom", 1, "n1", 0, {"legacy_sha1": "s1"}),
        b2.Remote("b.mbworld", 1, "s2", 0, {}),
        b2.Remote("c.mbworld", 1, "s3", 0, {}),
        b2.Remote("c.bloom", 1, "other", 0, {}),
    ]
    assert cs.next_remote(remotes, None).name == "b.mbworld"
    assert cs.next_remote(remotes, "c.mbworld") is None


def test_backfill_reaches_a_dataclass_held_by_a_plain_object():
    @dataclasses.dataclass
    class Cache:
        added: list = dataclasses.field(default_factory=lambda: [7])

    class Plate:  # a plain class, like the real plates
        pass

    cache = Cache()
    del cache.__dict__["added"]
    plate = Plate()
    plate.cache = cache
    assert cs.backfill_dataclass_defaults([plate, {"again": plate}]) == ["Cache.added"]
    assert cache.added == [7]


def test_one_unsaveable_file_is_reported_and_does_not_stop_the_batch(tmp_path, monkeypatch):
    (tmp_path / "a.mbworld").write_bytes(b"not a zip")
    (tmp_path / "b.mbworld").write_bytes(ZIP)

    def boom(*_args, **_kwargs):
        raise AttributeError("a field the writer cannot find")

    monkeypatch.setattr(cs, "convert_pickle_bytes", boom)
    args = cs.argparse.Namespace(dirs=[str(tmp_path)], dry_run=False, verify_step=False)
    assert cs.cmd_local(args) == 1
    assert (tmp_path / "b.bloom").exists()  # the zip still moved
