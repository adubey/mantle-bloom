"""Serve the packaged SPA from either application entry point."""

import os
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles


def mount_frontend(app: FastAPI, dist_dir: str | os.PathLike[str]) -> None:
    """Serve a `vite build` output directory (its `index.html` + `assets/`) at `/`.

    Mounted last, after every API route, so `/worlds/{world_id}/*` and `/docs` keep priority; anything
    else falls through to the bundle, with unknown non-API paths served `index.html` so a
    hard refresh on a client-side view doesn't 404."""
    dist = Path(dist_dir)
    index = dist / "index.html"
    if not index.is_file():
        raise RuntimeError(f"no index.html under {dist!r} -- run `npm run build` in frontend/ first")

    app.mount("/assets", StaticFiles(directory=dist / "assets"), name="assets")

    @app.get("/", include_in_schema=False)
    def _index() -> FileResponse:
        return FileResponse(index)

    @app.get("/{path:path}", include_in_schema=False)
    def _spa_fallback(path: str) -> FileResponse:
        candidate = (dist / path).resolve()
        if dist.resolve() in candidate.parents and candidate.is_file():
            return FileResponse(candidate)
        return FileResponse(index)
