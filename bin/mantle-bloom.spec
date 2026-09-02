# PyInstaller spec for the single-process desktop build.
#
#   cd backend && source .venv/bin/activate && pip install -r requirements-build.txt
#   cd frontend && npm ci && npm run build && cd ..
#   pyinstaller bin/mantle-bloom.spec --noconfirm
#
# Output: dist/mantle-bloom/ (a folder build -- onedir, not onefile: numba/llvmlite and the
# ffmpeg libs PyAV carries make onefile slow to start and awkward to sign). On macOS you also
# get dist/mantle-bloom.app. Build on the OS you're targeting -- PyInstaller does not
# cross-compile.
from pathlib import Path

from PyInstaller.utils.hooks import collect_dynamic_libs, collect_submodules

REPO = Path(SPECPATH).parent
BACKEND = REPO / "backend"
FRONTEND_DIST = REPO / "frontend" / "dist"

if not (FRONTEND_DIST / "index.html").is_file():
    raise SystemExit("build the frontend first: cd frontend && npm run build")

datas = [(str(FRONTEND_DIST), "frontend_dist")]

# numba compiles at runtime and needs llvmlite's shared lib; PyAV carries its own ffmpeg.
binaries = collect_dynamic_libs("llvmlite") + collect_dynamic_libs("av")

hiddenimports = (
    collect_submodules("uvicorn")
    + collect_submodules("app")
    + ["numba", "llvmlite", "scipy._lib.array_api_compat", "astropy_healpix"]
)

a = Analysis(
    [str(BACKEND / "app" / "desktop.py")],
    pathex=[str(BACKEND)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    excludes=["tkinter", "pytest", "httpx", "matplotlib", "IPython"],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="mantle-bloom",
    console=True,
    disable_windowed_traceback=False,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    name="mantle-bloom",
)
app = BUNDLE(
    coll,
    name="mantle-bloom.app",
    icon=None,
    bundle_identifier="com.mantlebloom.app",
    info_plist={"LSBackgroundOnly": False, "NSHighResolutionCapable": True},
)
