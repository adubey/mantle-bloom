# Issue 920: security second pass

Status on 2026-11-15: the dependency scan is done for Python and npm, the frontend XSS
audit is done, and the base-image scan (`python:3.14-slim` and `node:26-slim`) was completed
under [#926](https://github.com/adubey/mantle-bloom/issues/926). Most of the penetration
test is **not** done. One finding that a logged-out visitor could reach was fixed in the same
branch. The launch-readiness §4 security box stays open.

| Part of #920 | State |
|---|---|
| Python dependency scan | Done: 45 advisories in 3 packages, triaged below (cleared in #926) |
| npm dependency scan | Done: 3 advisories, all in dev tooling, cleared by a lockfile update |
| Base-image scan | Done (#926): Trivy 0.76.2 on `python:3.14-slim` and `node:26-slim` (0 critical) |
| Frontend XSS audit | Done: no HTML injection sink found |
| Penetration test | Partly done, on a local hosted copy only; see [what was not run](#what-was-not-run) |

## Findings

Severity is for the service as deployed today: one owner account, logged-out visitors on
`/public/*`. "Owner only" means the request needs the owner's session, so today the owner
can only harm their own server. Each of those becomes reachable by other people once
accounts exist ([#560](https://github.com/adubey/mantle-bloom/issues/560)).

| # | Finding | Who can reach it | Severity | State |
|---|---|---|---|---|
| S1 | Request bodies had no size limit; the worker buffered a visitor's step body whole | Anyone | Medium | Fixed in this branch |
| S2 | A v5 save's compressed entries expand without a limit on load | Owner only | Low now, High with accounts | [#925](https://github.com/adubey/mantle-bloom/issues/925) |
| S3 | Pillow 12.3.0 has 43 advisories; the sketch upload decodes any image format | Owner only | Low now, High with accounts | [#926](https://github.com/adubey/mantle-bloom/issues/926) |
| S4 | certifi 2026.7.22 bundles root CAs with compromised intermediates | Not reachable | Accepted until the next pin bump | [#926](https://github.com/adubey/mantle-bloom/issues/926) |
| S5 | The v5 loader imports whatever module a manifest names | Owner only | Low | [#925](https://github.com/adubey/mantle-bloom/issues/925) |
| S6 | A v5 manifest sets the HEALPix grid size with no upper bound | Owner only | Low | [#925](https://github.com/adubey/mantle-bloom/issues/925) |
| S7 | `POST /worlds/import` doesn't apply the 100-character name limit | Owner only | Low | [#925](https://github.com/adubey/mantle-bloom/issues/925) |
| S8 | Save loading runs on the worker's event loop | Owner only | Low | [#925](https://github.com/adubey/mantle-bloom/issues/925) |
| S9 | `node:20-slim` is past Node 20's end of life | Build stage only | Low | [#926](https://github.com/adubey/mantle-bloom/issues/926) |
| S10 | Logging out doesn't invalidate a copied session cookie | Needs a stolen cookie | Low | Accepted, recorded here |
| S11 | Generate takes unbounded `num_plates` and `num_mantle_centers` | Owner only | Low | [#925](https://github.com/adubey/mantle-bloom/issues/925) |

S1 and S2 were measured. S3 to S11 come from reading the code and the scan output and were
not exercised.

### S1: unbounded request bodies (fixed)

FastAPI reads a JSON body into memory before validating it. Nothing bounded that read, and
the hosted public process streams a step body over 1 KiB straight to the worker.

Measured on a local hosted copy (`uvicorn app.hosted:app`, loopback, throwaway OAuth values)
with the worker awake: a `POST /public/worlds/<any id>/step` with no session and a 200 MiB
body of spaces was read to the end (the `422` named byte 209,715,200) and the worker's
resident memory rose from 210 MiB to 463 MiB. The world id did not have to exist, because
the body is parsed before the route looks the world up.

Fix: `backend/app/body_limit.py`, installed by both entry points. After the fix, on the same
copy: a declared 8 KiB body gets `413`, a 200 MiB chunked body is stopped at the cap, a normal step body is answered as before, and the worker's memory did not move
(213 MiB before and after). Whether Railway's edge also caps bodies was not checked.

Verify: `backend/unit_tests/test_body_limit.py`.

### S2 and S6: a small save can expand without limit

`npz_save._load_world_v5_from_zip` reads `manifest.json` and every array with no limit on
the uncompressed size. Measured in-process: a 299 KiB zip whose manifest is 300 MiB of
whitespace raised the loader's peak memory from 653 MiB to 1,050 MiB before it was rejected
as corrupt. The upload cap is 512 MB of compressed bytes, so the same ratio allows far more
than a hosted container holds. Separately, the `__healpix_grid__` tag passes a
manifest-supplied `nside` to `healpix_grid.build`, which allocates per pixel (12 × nside²).

Suggested fix: on the upload routes, refuse a zip whose declared uncompressed sizes pass a
ceiling before reading any entry, and bound `nside` to the values the simulator produces.
Verify with a test that such a file raises `CorruptSaveError` without the memory growth.

### S3 and S4: Python advisories

Scanned all 46 pins in `backend/constraints.txt` against OSV on 2026-11-15. Results are in
[python-osv.json](python-osv.json).

| Package | Advisories | Fixed in | Reachable here |
|---|---|---|---|
| pillow 12.3.0 | 43 (5 critical, 24 high, 14 moderate) | 13.0.0 | Partly, owner only |
| certifi 2026.7.22 | 1 critical | 2026.10.18 | No |
| setuptools 84.0.0 | 1 low | 84.1.0 | No |

**Pillow.** The app writes PNG and MP4 frames from its own arrays. The one place it decodes
bytes it did not make is `worldsketch.parse_sketch_image`, which calls `Image.open` on the
sketch from `POST /worlds/{world_id}/generate` with no format restriction. The file picker
sends whatever image the owner chose (`accept="image/*"`), so restricting to PNG would break
JPEG sketches and was not done here.

| Group of advisories | Reachable |
|---|---|
| Decoders that `Image.open` selects by content (TIFF, GIF, BMP, JPEG, JPEG 2000, APNG, PNG eXIf, MSP, SPIDER, MPO, DCX, XBM) | Yes, by the owner through the sketch upload |
| Decoders needing something the image lacks: FlashPix and MIC (`olefile`), EPS (Ghostscript), WMF and EMF (Windows) | No |
| Encoders other than PNG (TIFF, PDF, GIF and JPEG save paths) | No, the app saves PNG only |
| Operations the app never calls (rank filters, morphology, `merge`, `transform`, `quantize`, Arrow export, stroked text, `TarIO`) | No |

Suggested fix: bump Pillow with `./bin/pin_deps.sh` under its own issue, with a before and
after replay as the pin rules require, and pass `formats=` to `Image.open` with the formats
the file picker should accept.

**certifi.** `httpx` imports it, but the only `httpx` client talks to the worker over a Unix
socket with no TLS. The GitHub OAuth calls use `urllib`, which reads the system trust store.
No TLS connection in the app verifies against certifi's bundle. Bump it with Pillow.

**setuptools.** The advisory is about installing legacy wheels on Windows. setuptools is a
build-time pin and is not in the runtime image.

### npm advisories

`npm audit` on 2026-11-15, results in [npm-audit.json](npm-audit.json).

| Lockfile | Before | After |
|---|---|---|
| `frontend/package-lock.json` | 3 high (vite 8.2.x, source-map-js, zod) | 0 |
| `frontend/package-lock.json`, production dependencies only | 0 | 0 |
| `.railway/package-lock.json` | 0 | 0 |

All three were in dev tooling: the Vite dev server, PostCSS's source-map reader, and knip's
validator. None is in the built bundle, which has four production dependencies.
`npm audit fix` moved vite to 8.3.8, source-map-js to 1.2.2 and zod to 4.6.7 within the
existing version ranges; `npm run lint`, `npm run build` and `npm test` pass afterwards.

### S5: the loader imports modules named by the manifest

`npz_save._resolve` calls `importlib.import_module` on a manifest-supplied dotted name and
only then checks that the result is a dataclass. Importing runs the module's top-level
code, and any importable dataclass can be instantiated with manifest-chosen arguments.
Suggested fix: accept only names under `app.` before importing.

### S7 to S11

- **S7.** `CreateWorldRequest`, `PatchWorldRequest` and `GenerateRequest` bound a name to
  1-100 characters. `POST /worlds/import` takes its `name` query parameter, or the name in
  the save's manifest, with no bound, and the name is rewritten to `worlds.json` on every
  record save.
- **S8.** `load_world` and `import_world` are `async` and call `persistence.load_world_file`
  directly, so a slow load stalls every request the worker is serving, public long-polls
  included.
- **S9.** The frontend build stage uses `node:20-slim`. Node 20 left support in April 2026,
  so the image stops receiving fixes. It builds the bundle and is not in the final image.
  CI pins the same version, so both move together.
- **S10.** Sessions are stateless signed cookies valid for 14 days. `POST /logout` clears
  the browser's copy; a copy taken earlier keeps working until it expires or the session
  secret is rotated. The cookie is `Secure`, `HttpOnly` and `SameSite=Lax`, so taking a copy
  needs access to the owner's browser. Accepted for a single-owner service; revisit with
  accounts.
- **S11.** `GenerateRequest.num_plates` and `num_mantle_centers` are unbounded integers,
  unlike `voronoi_points` and the density choices.

## Frontend XSS audit

Reviewed every `.ts` and `.tsx` file under `frontend/src` (13,732 lines) and the backend's
HTML responses.

**No HTML sink is used.** There is no `dangerouslySetInnerHTML`, `innerHTML`, `outerHTML`,
`insertAdjacentHTML`, `document.write`, `eval`, `new Function`, `srcdoc`, `DOMParser`,
`window.open` or `postMessage` anywhere in the frontend. All server-supplied and
user-supplied text (world names, manifest names, event-log messages, error `detail` strings,
phase labels) is rendered as React text children, which React escapes.

Sinks reviewed:

| Sink | Where | Data | Result |
|---|---|---|---|
| `<a href>` | `Viewer.tsx`, `NotFound.tsx`, `RunOverlay.tsx` | Constants, and the static link table in `phaseInfo.ts` | Safe |
| `<a href={debugUrl()}>` | `App.tsx` | The page's own URL, parsed by `URL` | Safe: the scheme can't change |
| `<form action>` | `AdminSession.tsx` | `${API_BASE}/logout`, build-time constant | Safe |
| `img.src`, `<img src>`, SVG `<image href>` | `MapCanvas.tsx`, `PlatesAndFaults.tsx`, `globeTexture.ts`, `VoronoiPreview.tsx`, `SketchEditor.tsx`, `MantleCrossSection.tsx`, `App.tsx` | `data:image/png;base64,` plus server base64, a blob URL, or the user's own sketch | Safe: an image context doesn't run script |
| `a.href` and `a.download` on a temporary link | `App.tsx`, `FileModal.tsx` | A blob URL; file name built from the numeric seed and age | Safe |
| `fetch("data:<mime>;base64,...")` | `App.tsx` (save animation) | Server-chosen MIME and base64 | Safe: read as a blob and downloaded, never rendered |
| `window.location.replace` | `AdminSession.tsx` | Constant `/view` | Safe |
| `history.replaceState` | `App.tsx`, `Viewer.tsx` | The page URL with `world` set through `URLSearchParams` | Safe |
| Request paths | `worldRoutes.ts` | `?world=` from the page URL, passed through `encodeURIComponent` | Safe: can't add path segments |
| `document.cookie` | `cookies.ts` | View-state names and values, both `encodeURIComponent`-encoded | Safe |
| Inline `style` values | `Legend.tsx`, `TimeSeriesChart.tsx`, `MantleCrossSection.tsx` | Colours from static palettes | Safe |

One note: `apiResponse.expandColumnarStatsHistory` assigns server-supplied field names onto
plain objects. A field named `__proto__` would change that one row's prototype, not
`Object.prototype`. The data comes from an owner-only route. Not a finding.

Backend HTML:

- `/login` is a fixed template. `?error=` selects one of five fixed messages, which is
  passed through `html.escape`; an unknown code shows nothing.
- `/docs` and `/redoc` are owner only and are the only pages with a relaxed policy.
- Error `detail` strings that echo input (an unknown `node_density`, for example) are sent
  as `application/json` with `nosniff`, so a browser won't render them as HTML.
- `GET /worlds/{world_id}/save` builds its `Content-Disposition` file name from the integer
  seed and age only.
- The SPA fallback resolves the requested path and serves it only when it is inside the
  build directory; anything else gets `index.html`.

What the content security policy from #918 does:

- Blocks inline and third-party scripts (`script-src 'self'`, no `unsafe-inline`, no
  `unsafe-eval`), plugins, a changed `<base>`, forms posting elsewhere, framing by other
  sites, and fetches to other origins.
- Does not block inline styles (`style-src 'unsafe-inline'`), so injected markup could
  still restyle the page. No injection point for markup was found.
- Does not restrict `data:` and `blob:` images.
- Does not apply its resource rules to `/docs` and `/redoc`.
- Has no `report-uri`, so a blocked load is visible only in the browser console.

## Penetration test

### Scope

Target: a local copy of the hosted entry point on loopback, with throwaway OAuth values and
a throwaway session secret. The owner withdrew the deployed service as a target during the
work, so nothing was sent to it. Planned areas, from #920: auth and session handling,
authorization on `/worlds/*`, `/debug/*` and `/public/*`, the worker proxy, the guest step
allowance, resource exhaustion, and the upload path.

### What was run

| Check | Result |
|---|---|
| Logged-out `GET /auth/session` | `admin: false`, `local: false` |
| Security headers on a logged-out response | All present |
| Anonymous oversized step body, before the fix | Buffered whole by the worker (S1) |
| The same, after the fix, declared and chunked | Refused at the cap; worker memory unchanged |
| A session cookie signed with the copy's own secret | Accepted, and wakes the worker |
| Small v5 file with a highly compressed manifest, in-process | Expands without limit (S2) |

### What was not run

- Everything against the deployed service: its TLS and edge behaviour, the cookie flags as
  served, and whether its variables match the repository.
- OAuth state and PKCE handling, session forgery and replay, and the cross-site write
  checks. The code for each was read (see `auth.py`) and nothing wrong was seen, but none
  was exercised.
- Authorization bypass attempts on `/worlds/*` and `/debug/*`, and how the public process
  and the worker each interpret an unusual request path. The proxy forwards the raw path to
  the worker; whether the two always agree on which route a path names is untested.
- Render and animate bounds, the concurrent-world cap, and the 512 MB upload cap under load.
- Legacy pickles and malformed v5 saves through the upload routes.

## Base-image scan (#926)

Scanned `python:3.14-slim` (`python@sha256:9aaa474097d094040649f2162a627c0049b408f52138d980916788d80f04802f`)
and `node:26-slim` (`node@sha256:e6bedef217677b276c591795e169840edd9096dbf13f29df8ef8d1856b8820bc`,
replacing `node:20-slim`) with Trivy 0.76.2 on 2026-11-15. Both images are based on Debian 13.7
and have **0 critical** advisories.

| Image | Target | Critical | High | Medium | Low | Unknown | Fixable High/Critical |
|---|---|---|---|---|---|---|---|
| `python:3.14-slim` | Debian 13.7 OS | 0 | 49 | 93 | 65 | 5 | 2 (`zlib1g` CVE-2026-76844, CVE-2026-85091) |
| `python:3.14-slim` | Bundled Python (`pip`/`setuptools`/`urllib3`/`msgpack`) | 0 | 5 | 2 | 1 | 0 | 5 (`msgpack` 1.1.2, `setuptools` 70.3.0, `urllib3` 2.7.0) |
| `node:26-slim` (build stage only) | Debian 13.7 OS | 0 | 48 | 94 | 65 | 5 | 2 (`zlib1g` CVE-2026-76844, CVE-2026-85091) |
| `node:26-slim` (build stage only) | Bundled Node.js (`npm`'s `brace-expansion`, `undici`) | 0 | 6 | 12 | 1 | 0 | 6 (`brace-expansion` 5.0.9, `undici` 6.28.0) |

- `node:26-slim` is used only in the multi-stage `frontend` build stage (`Dockerfile`) to run
  `npm ci` and `npm run build`; neither Node.js nor its Debian packages enter the final image.
- In `python:3.14-slim`, the 47 non-`zlib1g` Debian HIGH advisories (`util-linux`/`bsdutils`/
  `libblkid1`/`libmount1`, `libacl1`, `libncursesw6`, etc.) have no upstream Debian 13 fix yet;
  the 5 Python-package HIGH advisories are in `pip`/base-image tooling (`msgpack`, `setuptools`,
  `urllib3`), not in `backend/constraints.txt` (which uses `setuptools==84.2.1` at build time and
  does not include `msgpack` or `urllib3`).

## Reproducing the scans

- Python: post each `name==version` from `backend/constraints.txt` to
  `https://api.osv.dev/v1/querybatch` with ecosystem `PyPI`. Only package names and versions
  leave the machine. `pip-audit -r backend/constraints.txt --no-deps --disable-pip` reads
  the same data.
- npm: `npm audit` and `npm audit --omit=dev` in `frontend/` and `.railway/`.
- Base images: `trivy image python:3.14-slim` and `trivy image node:26-slim`.
