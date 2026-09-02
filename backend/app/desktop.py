"""Entry point for the packaged single-process desktop build.

`python -m app.desktop` (or the PyInstaller binary that bundles it) starts the FastAPI app on
a free local port, has it serve the built frontend from that same origin, and opens a browser
at it. No separate Vite server, no fixed :8000/:5173, no CORS.

Running from source needs the frontend built once (`cd frontend && npm run build`); the
frozen binary carries `frontend/dist` inside it (see bin/mantle-bloom.spec).
"""

from __future__ import annotations

import os
import socket
import sys
import threading
import time
import webbrowser
from pathlib import Path


def _frontend_dist() -> Path:
    """Locate the `vite build` output, whether we're frozen or running from a source tree."""
    if getattr(sys, "frozen", False):
        # PyInstaller unpacks bundled data under sys._MEIPASS (onefile) or next to the exe.
        base = Path(getattr(sys, "_MEIPASS", Path(sys.executable).parent))
        return base / "frontend_dist"
    # backend/app/desktop.py -> repo root -> frontend/dist
    return Path(__file__).resolve().parents[2] / "frontend" / "dist"


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _open_when_up(url: str, host: str, port: int) -> None:
    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                webbrowser.open(url)
                return
        except OSError:
            time.sleep(0.2)


def main() -> None:
    dist = _frontend_dist()
    if not (dist / "index.html").is_file():
        sys.exit(
            f"frontend build not found at {dist}\n"
            "Run `cd frontend && npm run build` first (or use the packaged binary)."
        )
    os.environ["MANTLE_BLOOM_FRONTEND_DIST"] = str(dist)

    host = "127.0.0.1"
    port = int(os.environ.get("MANTLE_BLOOM_PORT") or _free_port())
    url = f"http://{host}:{port}/"

    import uvicorn

    from app.main import app

    print(f"mantle-bloom running at {url}  (Ctrl+C to quit)", flush=True)
    # MANTLE_BLOOM_NO_BROWSER=1 keeps this from spawning a tab -- bin/restart.sh sets it so a
    # dev restart doesn't open a new browser window every time.
    if not os.environ.get("MANTLE_BLOOM_NO_BROWSER"):
        threading.Thread(target=_open_when_up, args=(url, host, port), daemon=True).start()
    uvicorn.run(app, host=host, port=port, log_level="warning")


if __name__ == "__main__":
    main()
