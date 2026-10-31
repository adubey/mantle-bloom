"""Guard the Gate 3b probes' cross-script references, which no other test imports (#699)."""
import ast
import importlib
import sys
from pathlib import Path

import pytest

DEBUG_DIR = Path(__file__).resolve().parents[2] / "bin" / "debug"
PROBES = sorted(DEBUG_DIR.glob("probe_gate3b_*.py"))


@pytest.mark.parametrize("path", PROBES, ids=lambda path: path.stem)
def test_probe_references_exist_in_the_sibling_scripts_it_imports(path, monkeypatch):
    monkeypatch.setattr(sys, "path", [str(DEBUG_DIR), *sys.path])
    tree = ast.parse(path.read_text())
    siblings = {
        alias.asname or alias.name: alias.name
        for node in ast.walk(tree) if isinstance(node, ast.Import)
        for alias in node.names
        if (DEBUG_DIR / f"{alias.name}.py").exists()
    }
    modules = {name: importlib.import_module(module) for name, module in siblings.items()}
    missing = [
        f"{node.value.id}.{node.attr}"
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
        and node.value.id in modules and not hasattr(modules[node.value.id], node.attr)
    ]
    assert not missing
