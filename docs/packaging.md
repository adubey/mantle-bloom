# Packaging mantle-bloom as a desktop binary

Turns the two-process dev setup (uvicorn on `:8000` + Vite on `:5173`) into **one
self-contained executable** that a non-developer can double-click. The backend serves the
built frontend from its own origin, picks a free port, and opens a browser at it.

This is a **local single-user app**, not a deployment. For putting mantle-bloom on the
public internet see [HOSTING.md](HOSTING.md) instead — the security blockers there
(pickle RCE in `/world/load`, single in-memory world, localhost CORS) do **not** matter
when the whole thing runs on one person's machine, which is the only thing this build
targets.

## How it works

| Piece | Dev (`bin/restart.sh --dev`) | Single-process (`bin/restart.sh` / `bin/package.sh`) |
| --- | --- | --- |
| Frontend | Vite dev server on `:5173` (HMR) | `frontend/dist`, served by FastAPI at `/` (bundled into the binary when frozen) |
| API base | `VITE_API_BASE=http://localhost:8000` baked in | `""` — same origin as the page (`api.ts`, `import.meta.env.PROD`) |
| Port | backend `8000` + frontend `5173` | one port: `8000` from `restart.sh`, a free ephemeral one from the binary (`MANTLE_BLOOM_PORT` overrides) |
| CORS | localhost allowlist | irrelevant — same origin |
| Entry point | `uvicorn app.main:app` | `app/desktop.py` → `uvicorn.run(app)` (+ open browser, unless `MANTLE_BLOOM_NO_BROWSER`) |
| Process count | 2 | 1 |

The two moving parts added for this:

- **`backend/app/main.py` → `mount_frontend()`** — mounts a `vite build` output dir at `/`,
  *after* every API route so `/world/*` and `/docs` keep priority, with an SPA fallback to
  `index.html`. Only active when `MANTLE_BLOOM_FRONTEND_DIST` is set, so the dev workflow is
  untouched.
- **`backend/app/desktop.py`** — the frozen-app entry point. Locates `frontend/dist`
  (bundled next to the exe when frozen, `../frontend/dist` from source), sets the env var,
  starts the server on a free port, opens a browser when it's up.

Running the single-process app from source without packaging anything is just the default
`./bin/restart.sh` (it builds the frontend and runs `python -m app.desktop` on `:8000`), or
by hand:

```bash
cd frontend && npm run build && cd ..
cd backend && source .venv/bin/activate && python -m app.desktop
```

## Building the binary

```bash
./bin/package.sh
```

Does: `npm ci && npm run build` in `frontend/`, `pip install -r requirements-build.txt`
in the backend venv, then `pyinstaller bin/mantle-bloom.spec`. Output:

- **`dist/mantle-bloom/`** — a onedir bundle (the folder *is* the app; ship the whole
  folder, zipped).
- **`dist/mantle-bloom.app`** — macOS app bundle wrapping the same thing.

Roughly **~290 MB** unpacked (~100 MB of that is LLVM, shipped for numba's JIT; the rest
is numpy/scipy/Pillow and the ffmpeg libraries PyAV carries).

### Why onedir, not onefile

onefile unpacks the whole ~290 MB to a temp dir on every launch (slow first paint) and the
bundled `.dylib`/`.so`/`.dll` files inside a single executable make code-signing and
notarization harder. onedir starts fast and signs cleanly.

## Per-platform

**PyInstaller does not cross-compile.** Build the macOS app on macOS, the Windows `.exe`
on Windows. Use a CI matrix (`.github/workflows/package.yml`) so both come out of one tag.

### macOS

- The `.app` runs locally as-is. To hand it to someone else without Gatekeeper blocking it,
  **codesign + notarize** (needs a paid Apple Developer account):
  ```bash
  codesign --deep --force --options runtime --sign "Developer ID Application: NAME (TEAMID)" dist/mantle-bloom.app
  ditto -c -k --keepParent dist/mantle-bloom.app mantle-bloom.zip
  xcrun notarytool submit mantle-bloom.zip --apple-id … --team-id … --password … --wait
  xcrun stapler staple dist/mantle-bloom.app
  ```
- The build is single-arch (whatever the build machine is). For a universal2 app you'd
  need universal2 wheels for numpy/scipy/numba/av — not worth it; ship separate
  arm64 / x86_64 builds or just arm64.
- A harmless build warning — `Library not found: @rpath/libomp.dylib` — is fine: numba
  falls back to its `workqueue` threading layer, which is the one `main.py` already
  documents relying on. (Bundle `libomp` only if you switch numba to the OpenMP layer.)

### Windows

- No signing needed to run, but SmartScreen will warn on first launch of an unsigned exe.
  An OV/EV code-signing cert removes that.
- Build inside a clean venv on a `windows-latest` runner; confirm `av` and `llvmlite`
  wheels install (they do, on CPython 3.12/3.13).

### Python version

The repo runs on **CPython 3.14**. PyInstaller 6.x, numba ≥ 0.67, and llvmlite ≥ 0.49
support it, but if a build machine has trouble, pin the build venv to 3.12/3.13 — the app
code is compatible and the wheels are unambiguously available there.

## Gotchas / verified

Smoke-tested from the frozen macOS binary:

- ✅ `/world/generate`, `/world/step` (numba JIT compiles on first step — adds ~1.5 s once,
  then cached), `/world/render` (Pillow), `/world/animate` (PyAV encodes a valid MP4).
- ✅ `/world/save` + `/world/load` round-trip (pickle).
- ✅ `scipy.spatial.cKDTree`, `astropy_healpix`.
- ✅ GET inspector routes (`/world/rivers` etc.) and `/docs` are not shadowed by the SPA
  catch-all; a real API 404 still returns its JSON `detail`, not `index.html`.
- The `.spec` explicitly `collect_submodules`es `uvicorn` (dynamic protocol imports) and
  `app`, and `collect_dynamic_libs`es `llvmlite` and `av`.

Not yet exercised: the Windows build (no Windows machine — relying on CI).

## Alternatives considered

- **python-build-standalone + a launcher script** — ship an embedded interpreter + the
  venv + a `.command`/`.bat`. Not a single file, but no PyInstaller hook-chasing and no
  code-signing fight. Fall back to this if the frozen build fights you.
- **Nuitka** — compiles to C; smaller, faster startup, but longer builds and its own set
  of numba/scipy edge cases. Worth trying if binary size matters.
- **Tauri / Electron shell** — only if you want a real app window and menu bar instead of
  a browser tab. Still needs the Python backend as a sidecar (i.e. everything above plus a
  shell).
