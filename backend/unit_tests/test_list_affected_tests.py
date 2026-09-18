"""Tests for bin/list_affected_tests.py's import-graph and git-diff logic -- see its
docstring for what "affected" means. Runs the real script as a subprocess against a
throwaway git repo shaped like backend/app + backend/unit_tests, since the script's
REPO_ROOT/APP_DIR/etc. are derived from its own file location and its diffing needs a real
git repo -- there's no meaningful way to unit-test that in-process."""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "bin" / "list_affected_tests.py"


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=repo, check=True, capture_output=True, text=True
    ).stdout


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


@pytest.fixture
def fake_repo(tmp_path) -> Path:
    """backend/app/{a,b,c}.py chained a <- b <- c (c imports b, b imports a), each with its
    own backend/unit_tests/test_*.py importing only its own module, committed on `main`."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    _git(repo, "symbolic-ref", "HEAD", "refs/heads/main")  # git init -b needs git 2.28+
    _git(repo, "config", "user.email", "test@example.com")
    _git(repo, "config", "user.name", "Test")

    _write(repo / "backend/app/__init__.py", "")
    _write(repo / "backend/app/a.py", "VALUE = 1\n")
    _write(repo / "backend/app/b.py", "from . import a\nVALUE = a.VALUE + 1\n")
    _write(repo / "backend/app/c.py", "from .b import VALUE as b_value\nVALUE = b_value + 1\n")
    _write(repo / "backend/unit_tests/__init__.py", "")
    _write(repo / "backend/unit_tests/test_a.py", "from app import a\n\n\ndef test_a():\n    assert a.VALUE == 1\n")
    _write(repo / "backend/unit_tests/test_b.py", "from app.b import VALUE\n\n\ndef test_b():\n    assert VALUE == 2\n")
    _write(repo / "backend/unit_tests/test_c.py", "from app import c\n\n\ndef test_c():\n    assert c.VALUE == 3\n")

    script_dest = repo / "bin" / "list_affected_tests.py"
    script_dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(SCRIPT, script_dest)

    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "initial")
    return repo


def _run(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(repo / "bin" / "list_affected_tests.py"), *args],
        cwd=repo,
        capture_output=True,
        text=True,
    )


def _affected(result: subprocess.CompletedProcess) -> set[str]:
    assert result.returncode == 0, result.stderr
    return {line for line in result.stdout.splitlines() if line.strip()}


def test_no_changes_reports_nothing(fake_repo):
    result = _run(fake_repo)
    assert _affected(result) == set()
    assert "no changes" in result.stderr


def test_leaf_module_change_affects_only_its_own_test(fake_repo):
    (fake_repo / "backend/app/c.py").write_text("from .b import VALUE as b_value\nVALUE = b_value + 100\n")
    assert _affected(_run(fake_repo)) == {"unit_tests/test_c.py"}


def test_change_propagates_to_transitive_importers(fake_repo):
    # a is imported by b, b by c -- a change to a should reach every test in the chain
    (fake_repo / "backend/app/a.py").write_text("VALUE = 999\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_changed_test_file_itself_is_affected(fake_repo):
    (fake_repo / "backend/unit_tests/test_b.py").write_text(
        "from app.b import VALUE\n\n\ndef test_b():\n    assert VALUE == 2\n\n\ndef test_b_extra():\n    pass\n"
    )
    assert _affected(_run(fake_repo)) == {"unit_tests/test_b.py"}


def test_uncommitted_changes_are_picked_up(fake_repo):
    (fake_repo / "backend/app/c.py").write_text("from .b import VALUE as b_value\nVALUE = b_value + 5\n")
    assert _affected(_run(fake_repo)) == {"unit_tests/test_c.py"}


def test_new_leaf_module_with_no_importers_affects_nothing(fake_repo):
    # nothing imports it yet, so nothing currently exercises it -- correctly, not a fallback
    _write(fake_repo / "backend/app/new_module.py", "VALUE = 1\n")
    assert _affected(_run(fake_repo)) == set()


def test_new_module_wired_into_an_existing_one_affects_its_importers(fake_repo):
    _write(fake_repo / "backend/app/new_module.py", "VALUE = 1\n")
    (fake_repo / "backend/app/a.py").write_text("from . import new_module\nVALUE = new_module.VALUE\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_changed_app_init_falls_back_to_full_suite(fake_repo):
    (fake_repo / "backend/app/__init__.py").write_text("# no longer empty\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_syntax_error_in_app_module_falls_back_to_full_suite(fake_repo):
    (fake_repo / "backend/app/a.py").write_text("def broken(:\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_deleted_module_still_propagates_to_a_stale_importer(fake_repo):
    # b.py still says `from . import a` even though a.py is gone -- a real bug (b will
    # ImportError), but that's exactly when running b's test matters most
    (fake_repo / "backend/app/a.py").unlink()
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_deleted_test_file_is_not_listed(fake_repo):
    (fake_repo / "backend/unit_tests/test_c.py").unlink()
    assert _affected(_run(fake_repo)) == set()


def test_absolute_self_import_inside_app_is_tracked(fake_repo):
    # a couple of real app/ files import absolutely rather than relatively, e.g. desktop.py's
    # `from app.main import app` -- that must count as a dependency edge just like `from . import`
    _write(fake_repo / "backend/app/d.py", "from app.a import VALUE\n")
    _write(
        fake_repo / "backend/unit_tests/test_d.py",
        "from app import d\n\n\ndef test_d():\n    assert d.VALUE == 1\n",
    )
    _git(fake_repo, "add", "-A")
    _git(fake_repo, "commit", "-q", "-m", "add d")

    (fake_repo / "backend/app/a.py").write_text("VALUE = 2\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
        "unit_tests/test_d.py",
    }


def test_change_under_app_subdirectory_falls_back_to_full_suite(fake_repo):
    # e.g. the real backend/app/data/major_plates.json that real_plates.py reads -- not a
    # top-level app/*.py module, so it can't be mapped and must not be silently dropped
    _write(fake_repo / "backend/app/data/plates.json", "{}\n")
    assert _affected(_run(fake_repo)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }


def test_base_override(fake_repo):
    (fake_repo / "backend/app/a.py").write_text("VALUE = 2\n")
    _git(fake_repo, "add", "-A")
    _git(fake_repo, "commit", "-q", "-m", "bump a")
    previous = _git(fake_repo, "rev-parse", "HEAD~1").strip()

    # default base (main, == HEAD) sees no changes: the bump is already committed on main
    assert _affected(_run(fake_repo)) == set()
    # explicitly diffing against the commit before it surfaces the whole downstream
    assert _affected(_run(fake_repo, "--base", previous)) == {
        "unit_tests/test_a.py",
        "unit_tests/test_b.py",
        "unit_tests/test_c.py",
    }
