# The hosted app (#376): one container serving the built frontend and the FastAPI backend from
# one origin, as the desktop build does. Railway builds it as .railway/railway.ts configures;
# see docs/hosting.md for the variables and volume it needs.

# Node 26 matches CI's frontend job (.github/workflows/test.yml).
FROM node:26-slim AS frontend
WORKDIR /src/frontend
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY frontend/ ./
# Empty: the API is on the same origin as the page (api.ts).
ENV VITE_API_BASE=""
RUN npm run build

# Build the Rust glacier-march PyO3 extension targeting Python 3.14 on Linux.
FROM python:3.14-slim AS rust-builder
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates build-essential \
    && rm -rf /var/lib/apt/lists/* \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable
ENV PATH="/root/.cargo/bin:${PATH}"
WORKDIR /src
COPY backend/rust ./backend/rust
COPY bin/build_rust.py ./bin/build_rust.py
RUN python3 bin/build_rust.py --locked --out-dir /rust-ext

# Python 3.14 matches CI, which tests against the same backend/constraints.txt pins.
FROM python:3.14-slim
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1
WORKDIR /app/backend
# Runtime deps only, at the pinned versions: no pytest/ruff/vulture in the image.
COPY backend/requirements-runtime.txt backend/constraints.txt ./
RUN pip install -r requirements-runtime.txt -c constraints.txt
COPY --from=rust-builder /rust-ext/mantle_bloom_rust.so /tmp/mantle_bloom_rust.so
RUN python3 -c "import pathlib, shutil, sysconfig; shutil.move('/tmp/mantle_bloom_rust.so', str(pathlib.Path(sysconfig.get_path('platlib')) / 'mantle_bloom_rust.so'))" \
    && python3 -c "import mantle_bloom_rust"
COPY backend/app ./app
COPY --from=frontend /src/frontend/dist /app/frontend/dist

# Every request reaches the container through Railway's proxy, so nothing is local: turn the
# loopback no-login rule off rather than rely on the proxy's forwarding headers (auth.py).
ENV MANTLE_BLOOM_FRONTEND_DIST=/app/frontend/dist \
    MANTLE_BLOOM_LOCAL_NO_LOGIN=0

# Runs as root: Railway mounts the checkpoint volume root-owned, and a non-root user can't
# write to it (docs/hosting.md).
#
# The public process serves the frontend and persisted snapshots. It starts one
# private simulation worker when needed. Railway injects PORT.
CMD ["sh", "-c", "exec uvicorn app.hosted:app --host 0.0.0.0 --port \"${PORT:-8080}\" --timeout-graceful-shutdown 2"]
