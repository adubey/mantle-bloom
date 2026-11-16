# Frontend tests

Tests run under Node's built-in `node:test` runner, typechecked by `tsc -p tsconfig.test.json` and executed directly with `tsx`:

```bash
cd frontend && npm test
```

## Component testing setup

React component tests use `@testing-library/react` with `happy-dom` under `node:test` (rather than adding Vitest or Playwright Component Testing, which would introduce a second runner or browser binary).

- **`tests/domSetup.ts`** registers a `happy-dom` `Window` on `globalThis` and exports `resetDom(url)` for `beforeEach` isolation. Because `@testing-library/react` reads `globalThis.document` at module evaluation time, **import `./domSetup` before `@testing-library/react` or any component**.
- **`window.setInterval` / `window.clearInterval`** on the `happy-dom` window are forwarded to `globalThis` so `t.mock.timers.enable({ apis: ["setInterval"] })` drives polling loops deterministically.
- **Network requests** (`fetch`) are stubbed per test; long-poll requests (`?after=...`) should stay pending until their `AbortSignal` fires.

## Limits

- `happy-dom` does not perform CSS layout or 2D canvas rasterization (`canvas.getContext("2d")` returns `null`), so pixel rendering, canvas sizing, and pointer geometry cannot be asserted here.
- Backend API contracts and cookie/auth flows are mocked at `fetch`; end-to-end verification still requires a running server and browser.
