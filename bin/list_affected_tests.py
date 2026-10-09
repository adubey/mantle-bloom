#!/usr/bin/env python3
"""Prints the backend/unit_tests/ files affected by the changes on this branch, one per line
-- used by bin/affected_test.sh to run only the tests worth running instead of the whole
suite. "Affected" means: the test file's own module changed, or it imports (directly or
transitively, through app/'s own internal imports) a module that changed.

Change set = committed diff against the merge-base with origin/main (or main), unioned with
whatever's currently uncommitted (staged, unstaged, and untracked) -- so this is useful both
mid-edit and right after committing on a branch, without a --base flag to remember.

The dependency graph is built by statically parsing every app/**/*.py module and
unit_tests/test_*.py file with ast (not by importing them -- app/ pulls in numba/scipy at
import time, too slow to pay per invocation). Modules are keyed by their dotted name under
app/ ("world", "hydroclimate.climate"), and the parser understands this codebase's import
styles: relative imports at any depth (`from . import x`, `from .x import y`,
`from .. import geometry` inside app/hydroclimate/), and absolute imports of app's own
modules (`from app import x` / `from app.hydroclimate.climate import y` / `import app.x`)
wherever they appear -- tests use these exclusively, and a couple of app/ files do too (e.g.
desktop.py's `from app.main import app`). A new leaf module with no importers yet is handled
directly by the graph (it affects nothing, correctly -- nothing exercises it).

What this can't map to an app module or unit_tests/test_*.py file -- any changed
__init__.py, a non-Python file like app/data/*.json, a unit_tests/ helper such as
quad_fixtures.py, a parse error -- falls back to "affected == everything", as does a change to
the backend's dependency pins (FULL_SUITE_FILES), since a wrong "nothing affected" is a
silent gap in coverage where a wrong "run everything" just costs time.
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

# Changing these can change any test's outcome without touching an import edge.
FULL_SUITE_FILES = {
    BACKEND / "constraints.txt",
    BACKEND / "requirements.txt",
    BACKEND / "requirements-runtime.txt",
}
# Files outside backend/ that a test exercises without importing them.
EXTRA_TEST_DEPENDENCIES = {
    REPO_ROOT / "bin" / "list_affected_tests.py": TESTS_DIR / "test_list_affected_tests.py",
}


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
    """Every backend/app or backend/unit_tests file (plus FULL_SUITE_FILES and
    EXTRA_TEST_DEPENDENCIES) touched since `base`, plus whatever's currently uncommitted -- as
    absolute paths, additions/deletions/renames included (a rename is a deletion of the old
    path and an addition of the new one)."""
    changed = set()
    paths = ["backend/app", "backend/unit_tests"] + [
        str(f.relative_to(REPO_ROOT)) for f in sorted(FULL_SUITE_FILES | set(EXTRA_TEST_DEPENDENCIES))
    ]
    committed = run_git("diff", "--name-only", base, "HEAD", "--", *paths)
    working_tree = run_git("status", "--porcelain", "--untracked-files=all", "--", *paths)
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
    """Dotted name under app/ ("world", "hydroclimate.climate"), or None for an __init__.py or
    a non-Python file."""
    if path.suffix != ".py" or path.name == "__init__.py":
        return None
    return ".".join(path.relative_to(APP_DIR).with_suffix("").parts)


def imported_app_modules(path: Path, known_modules: set[str]) -> set[str] | None:
    """App modules this file imports, however deep in its body (app/world.py lazy-imports
    app/stats.py inside a method) -- or None if the file doesn't parse, which the caller
    treats as "assume it could depend on anything".

    Each imported name resolves to the longest prefix that's a module in `known_modules`:
    `from app.hydroclimate import climate` and `from app.hydroclimate.climate import X` both
    mean hydroclimate.climate. The caller puts deleted modules in `known_modules` too, so a
    module that still imports a just-deleted one keeps counting as depending on it -- otherwise
    removing a module would silently stop propagating to its callers' tests. A bare package
    import resolves to nothing; package __init__.py changes take the full-suite fallback."""
    try:
        tree = ast.parse(path.read_text(), filename=str(path))
    except SyntaxError:
        return None

    # the file's own package under app/ ([] for app/world.py, ["hydroclimate"] for
    # app/hydroclimate/climate.py), or None outside app/, where relative imports aren't ours
    package = list(path.relative_to(APP_DIR).parent.parts) if APP_DIR in path.parents else None
    deps: set[str] = set()

    def add(parts: list[str]) -> None:
        for n in range(len(parts), 0, -1):
            name = ".".join(parts[:n])
            if name in known_modules:
                deps.add(name)
                return

    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            if node.level:
                if package is None or node.level - 1 > len(package):
                    continue
                base = package[: len(package) - (node.level - 1)]
                if node.module:
                    base = base + node.module.split(".")
            elif node.module and node.module.split(".")[0] == "app":
                base = node.module.split(".")[1:]
            else:
                continue
            # each imported name is either a submodule (`from . import x`) or an attribute of
            # `base` itself, which the longest-prefix walk resolves to `base`
            for alias in node.names:
                add(base + [alias.name])
        elif isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                if parts[0] == "app":
                    add(parts[1:])
    return deps


def main() -> int:  # noqa: C901
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base", help="compare against this ref instead of the merge-base with main/origin-main"
    )
    args = parser.parse_args()

    app_files = sorted(f for f in APP_DIR.rglob("*.py") if "__pycache__" not in f.parts)
    app_module_names = {m for f in app_files if (m := module_name(f))}
    test_files = sorted(TESTS_DIR.glob("test_*.py"))

    base = args.base or merge_base_with_main()
    changed = changed_backend_files(base)
    if not changed:
        print("# no changes under backend/app or backend/unit_tests", file=sys.stderr)
        return 0

    full_suite_triggers = changed & FULL_SUITE_FILES
    if full_suite_triggers:
        names = ", ".join(str(f.relative_to(REPO_ROOT)) for f in sorted(full_suite_triggers))
        print(f"# {names} changed, falling back to full suite", file=sys.stderr)
        for t in test_files:
            print(t.relative_to(BACKEND))
        return 0

    # deleted modules stay resolvable, so their stale importers still count as affected
    known_modules = app_module_names | {
        m for f in changed if APP_DIR in f.parents and (m := module_name(f))
    }

    # app-internal dependency graph: module name -> module names it imports
    depends_on: dict[str, set[str]] = {}
    for f in app_files:
        mod = module_name(f)
        if mod is None:
            continue
        deps = imported_app_modules(f, known_modules)
        if deps is None:
            print(f"# could not parse {f.relative_to(REPO_ROOT)}, falling back to full suite", file=sys.stderr)
            for t in test_files:
                print(t.relative_to(BACKEND))
            return 0
        depends_on[mod] = deps

    # Walk every changed path against both directories recursively -- app/data/major_plates.json
    # is a real file real_plates.py reads, and a change to it (or to anything else this can't
    # resolve to an app module or unit_tests/test_*.py file) must fall into `unresolved` below,
    # not vanish silently.
    changed_modules: set[str] = set()
    unresolved: set[Path] = set()
    extra_tests: set[Path] = set()
    for f in changed:
        if f in EXTRA_TEST_DEPENDENCIES:
            extra_tests.add(EXTRA_TEST_DEPENDENCIES[f])
        elif APP_DIR in f.parents:
            mod = module_name(f)
            if mod is None:
                unresolved.add(f)
            else:
                changed_modules.add(mod)
        elif TESTS_DIR in f.parents:
            # matched by name, not by current existence -- a *deleted* test_*.py has nothing
            # left to run and isn't unresolved, it just drops out of test_files on its own
            if not (f.parent == TESTS_DIR and fnmatch(f.name, "test_*.py")):
                unresolved.add(f)

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

    affected_tests: set[Path] = {t for t in extra_tests if t.exists()}
    for t in test_files:
        if t in changed:
            affected_tests.add(t)
            continue
        deps = imported_app_modules(t, known_modules)
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
