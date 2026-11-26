#!/usr/bin/env python3
"""Build and install the `mantle_bloom_rust` PyO3 extension from `backend/rust/`.

Uses plain `cargo build --release` without `maturin`, keeping `backend/constraints.txt`
untouched. Installs the compiled native extension into the target Python's `platlib`
directory (or `--out-dir` for staged Docker builds).
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
RUST_DIR = REPO_ROOT / "backend" / "rust"
MODULE_NAME = "mantle_bloom_rust"


def _find_cargo() -> str:
    found = shutil.which("cargo")
    if found:
        return found
    for candidate in (
        Path.home() / ".cargo" / "bin" / "cargo",
        Path("/opt/homebrew/bin/cargo"),
        Path("/usr/local/bin/cargo"),
    ):
        if candidate.is_file():
            return str(candidate)
    raise SystemExit(
        "cargo not found. Install Rust locally (e.g. `brew install rust` or https://rustup.rs) "
        "or set MANTLE_BLOOM_SKIP_RUST=1 to skip building the Rust acceleration extension."
    )


def _target_filename() -> str:
    if sys.platform == "win32":
        return f"{MODULE_NAME}.pyd"
    return f"{MODULE_NAME}.so"


def _built_lib_path(target_dir: Path) -> Path:
    release_dir = target_dir / "release"
    if sys.platform == "darwin":
        candidates = [release_dir / f"lib{MODULE_NAME}.dylib", release_dir / f"lib{MODULE_NAME}.so"]
    elif sys.platform == "win32":
        candidates = [release_dir / f"{MODULE_NAME}.dll"]
    else:
        candidates = [release_dir / f"lib{MODULE_NAME}.so"]
    for path in candidates:
        if path.is_file():
            return path
    raise SystemExit(f"Compiled Rust library not found in {release_dir}")


def _resolve_platlib(python_exe: str) -> Path:
    if Path(python_exe).resolve() == Path(sys.executable).resolve():
        return Path(sysconfig.get_path("platlib"))
    out = subprocess.check_output(
        [python_exe, "-c", "import sysconfig; print(sysconfig.get_path('platlib'))"],
        text=True,
    ).strip()
    return Path(out)


def _newest_rust_mtime() -> float:
    watch_paths = [
        RUST_DIR / "Cargo.toml",
        RUST_DIR / "Cargo.lock",
        RUST_DIR / ".cargo" / "config.toml",
    ]
    src_dir = RUST_DIR / "src"
    if src_dir.is_dir():
        watch_paths.extend(src_dir.rglob("*.rs"))
    mtimes = [p.stat().st_mtime for p in watch_paths if p.is_file()]
    return max(mtimes, default=0.0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--python", default=sys.executable, help="Target Python interpreter")
    ap.add_argument("--out-dir", type=Path, help="Install directory (defaults to target Python's platlib)")
    ap.add_argument("--locked", action="store_true", help="Pass --locked to cargo build")
    ap.add_argument("--skip-if-unchanged", action="store_true", help="Skip cargo build if installed module is newer than Rust sources")
    args = ap.parse_args()

    out_dir = args.out_dir if args.out_dir is not None else _resolve_platlib(args.python)
    dest = out_dir / _target_filename()
    if args.skip_if_unchanged and dest.is_file() and dest.stat().st_mtime >= _newest_rust_mtime():
        return 0

    cargo = _find_cargo()
    target_dir = Path(os.environ.get("CARGO_TARGET_DIR", str(RUST_DIR / "target")))
    cmd = [cargo, "build", "--release"]
    if args.locked and (RUST_DIR / "Cargo.lock").is_file():
        cmd.append("--locked")

    env = os.environ.copy()
    env["PYO3_PYTHON"] = str(Path(args.python).resolve())
    subprocess.run(cmd, cwd=RUST_DIR, env=env, check=True)

    built = _built_lib_path(target_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, dest)
    print(f"Installed {MODULE_NAME} -> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
