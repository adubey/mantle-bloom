#!/usr/bin/env python3
"""Prints the backend/unit_tests/ files affected by the changes on this branch, one per line
-- used by bin/affected_test.sh to run only the tests worth running instead of the whole
suite. "Affected" means: the test file's own module changed, or it imports (directly or
transitively, through app/'s own internal imports) a module that changed.

Change set = committed diff against the merge-base with origin/main (or main), unioned with
whatever's currently uncommitted (staged, unstaged, and untracked) -- so this is useful both
mid-edit and right after committing on a branch, without a --base flag to remember.

The dependency graph is built by statically parsing app/*.py and unit_tests/*.py with ast
(not by importing them -- app/ pulls in numba/scipy at import time, too slow to pay per
invocation) and only understands this codebase's two import styles: relative imports within
app/ (`from . import x, y` / `from .x import y`) and absolute imports from tests
(`from app import x` / `from app.x import y` / `import app.x`). A new leaf module with no
importers yet is handled directly by the graph (it affects nothing, correctly -- nothing
exercises it). What this can't map to a single module -- a changed app/__init__.py or
unit_tests/__init__.py, a parse error -- falls back to "affected == everything", since a
wrong "nothing affected" is a silent gap in coverage where a wrong "run everything" just
costs time.
"""

from __future__ import annotations

import argparse
import ast
import subprocess
import sys
from fnmatch import fnmatch
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
BACKEND = REPO_ROOT / "backend"
APP_DIR = BACKEND / "app"
TESTS_DIR = BACKEND / "unit_tests"


def run_git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=REPO_ROOT, check=True, capture_output=True, text=True
    ).stdout


def merge_base_with_main() -> str:
    for candidate in ("origin/main", "main"):
        result = subprocess.run(
            ["git", "merge-base", "HEAD", candidate],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    raise SystemExit("could not find a merge-base with origin/main or main")


def changed_backend_files(base: str) -> set[Path]:
    """Every backend/app or backend/unit_tests file touched since `base`, plus whatever's
    currently uncommitted -- as absolute paths, additions/deletions/renames included (a
    rename is a deletion of the old path and an addition of the new one)."""
    changed = set()
    committed = run_git("diff", "--name-only", base, "HEAD", "--", "backend/app", "backend/unit_tests")
    working_tree = run_git("status", "--porcelain", "--untracked-files=all", "--", "backend/app", "backend/unit_tests")
    for line in committed.splitlines():
        if line.strip():
            changed.add(REPO_ROOT / line.strip())
    for line in working_tree.splitlines():
        # e.g. " M backend/app/world.py", "?? backend/app/new_file.py", "R  old -> new"
        path_part = line[3:].strip()
        if "->" in path_part:
            path_part = path_part.split("->", 1)[1].strip()
        if path_part:
            changed.add(REPO_ROOT / path_part)
    return changed


def module_name(path: Path) -> str | None:
    if path.suffix != ".py" or path.name == "__init__.py":
        return None
    return path.stem


def imported_app_modules(path: Path) -> set[str] | None:
    """App-package module names this file imports, however deep in its body (app/world.py
    lazy-imports app/stats.py inside a method) -- or None if the file doesn't parse, which
    the caller treats as "assume it could depend on anything". Every name this resolves is
    trusted as a submodule name even if that module file no longer exists on disk (app's own
    __init__.py is empty and always-empty is enforced by the __init__.py fallback below, so
    `from . import x` / `from app import x` can only ever mean the submodule `x`) -- that
    matters for deletions: a module that still imports a just-deleted one must keep counting
    as depending on it, or removing a module would silently stop propagating to its callers'
    tests."""
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except SyntaxError:
        return None

    is_app_file = path.parent == APP_DIR
    deps: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if is_app_file and node.level == 1:
                if node.module is None:
                    deps.update(alias.name for alias in node.names)
                else:
                    deps.add(node.module.split(".")[0])
            elif not is_app_file and node.level == 0 and node.module:
                parts = node.module.split(".")
                if parts[0] == "app":
                    if len(parts) > 1:
                        deps.add(parts[1])
                    else:
                        deps.update(alias.name for alias in node.names)
        elif isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == "app" and len(parts) > 1:
                    deps.add(parts[1])
    return deps


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base", help="compare against this ref instead of the merge-base with main/origin-main"
    )
    args = parser.parse_args()

    app_files = sorted(APP_DIR.glob("*.py"))
    app_module_names = {m for f in app_files if (m := module_name(f))}
    test_files = sorted(TESTS_DIR.glob("test_*.py"))

    base = args.base or merge_base_with_main()
    changed = changed_backend_files(base)
    if not changed:
        print("# no changes under backend/app or backend/unit_tests", file=sys.stderr)
        return 0

    # app-internal dependency graph: module name -> module names it imports
    depends_on: dict[str, set[str]] = {}
    for f in app_files:
        mod = module_name(f)
        if mod is None:
            continue
        deps = imported_app_modules(f)
        if deps is None:
            print(f"# could not parse {f.relative_to(REPO_ROOT)}, falling back to full suite", file=sys.stderr)
            for t in test_files:
                print(t.relative_to(BACKEND))
            return 0
        depends_on[mod] = deps

    changed_modules = {module_name(f) for f in changed if f.parent == APP_DIR}
    changed_modules.discard(None)
    unresolved_app_changes = {f for f in changed if f.parent == APP_DIR and module_name(f) is None}
    # matched by name, not by current existence -- a *deleted* test_*.py has nothing left to
    # run and isn't unresolved, it just drops out of test_files on its own
    unresolved_test_changes = {
        f for f in changed if f.parent == TESTS_DIR and not fnmatch(f.name, "test_*.py")
    }
    unresolved = unresolved_app_changes | unresolved_test_changes
    if unresolved:
        names = ", ".join(str(f.relative_to(REPO_ROOT)) for f in sorted(unresolved))
        print(f"# {names} isn't a plain app/test module (e.g. __init__.py), falling back to full suite", file=sys.stderr)
        for t in test_files:
            print(t.relative_to(BACKEND))
        return 0

    # reverse-reachability: every module that (transitively) depends on a changed module is
    # itself considered changed, since its behavior can shift too
    reverse_deps: dict[str, set[str]] = {m: set() for m in app_module_names}
    for mod, deps in depends_on.items():
        for dep in deps:
            reverse_deps.setdefault(dep, set()).add(mod)

    affected_modules: set[str] = set()
    frontier = list(changed_modules)
    while frontier:
        mod = frontier.pop()
        if mod in affected_modules:
            continue
        affected_modules.add(mod)
        frontier.extend(reverse_deps.get(mod, ()))

    affected_tests: set[Path] = set()
    for t in test_files:
        if t in changed:
            affected_tests.add(t)
            continue
        deps = imported_app_modules(t)
        if deps is None:
            affected_tests.add(t)  # unparseable test file -- always include it, it's cheap
            continue
        if deps & affected_modules:
            affected_tests.add(t)

    for t in sorted(affected_tests):
        print(t.relative_to(BACKEND))
    print(
        f"# {len(affected_tests)}/{len(test_files)} test files affected by "
        f"{len(changed_modules)} changed app module(s)",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
