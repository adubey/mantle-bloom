# The hosted app (#376): one container serving the built frontend and the FastAPI backend from
# one origin, as the desktop build does. Railway builds it from railway.toml; see
# docs/HOSTING.md for the variables and volume it needs.

# Node 20 matches CI's frontend job (.github/workflows/test.yml).
FROM node:20-slim AS frontend
WORKDIR /src/frontend
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY frontend/ ./
# Empty: the API is on the same origin as the page (api.ts).
ENV VITE_API_BASE=""
RUN npm run build

# Python 3.14 matches CI, which tests against the same backend/constraints.txt pins.
FROM python:3.14-slim
ENV PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1
WORKDIR /app/backend
# Runtime deps only, at the pinned versions: no pytest/ruff/vulture in the image.
COPY backend/requirements-runtime.txt backend/constraints.txt ./
RUN pip install -r requirements-runtime.txt -c constraints.txt
COPY backend/app ./app
COPY --from=frontend /src/frontend/dist /app/frontend/dist

# Every request reaches the container through Railway's proxy, so nothing is local: turn the
# loopback no-login rule off rather than rely on the proxy's forwarding headers (auth.py).
ENV MANTLE_BLOOM_FRONTEND_DIST=/app/frontend/dist \
    MANTLE_BLOOM_LOCAL_NO_LOGIN=0

# Runs as root: Railway mounts the checkpoint volume root-owned, and a non-root user can't
# write to it (docs/HOSTING.md).
#
# One process: the world lives in its memory, so more workers would each hold a different
# world. Railway injects PORT. The graceful-shutdown cap ends viewers' 25 s long polls on a
# redeploy instead of waiting them out; the checkpoint is already on disk.
CMD ["sh", "-c", "exec uvicorn app.main:app --host 0.0.0.0 --port \"${PORT:-8080}\" --timeout-graceful-shutdown 5"]
