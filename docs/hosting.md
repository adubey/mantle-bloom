# Hosting: current setup and next decisions

This page records Mantle Bloom's current hosted shape and the boundary of that design. The
full launch goals and open product decisions are in [launch readiness](../planning/launch-readiness.md).
Local development and localhost sign-in are covered by [local hosting](local-hosting.md).

## Current deployment

The public service runs on Railway from the repository's root [Dockerfile](../Dockerfile),
configured as code in [.railway/railway.ts](../.railway/railway.ts). The multi-stage image
builds the React frontend and installs the pinned Python runtime, then one Uvicorn process
serves both from the same origin. Railway builds from `main` and uses `GET /healthz` as its
health check. The Docker image is portable; Railway is the only deployment configured in
this repository.

The service is a small hosted preview, not a multi-user world service:

- One configured GitHub account has the full app and API.
- Logged-out visitors can use the public viewer and its limited step allowance.
- Everyone sees and steps the same active world. The world lives in one process's memory;
  there are no separate user-owned worlds or general user accounts.
- A Railway volume stores the latest world checkpoint and run status. On restart or wake, the
  process restores that state. If the service sleeps, simulation pauses until a request wakes
  it; elapsed time is the time simulated, not wall-clock time.
- One replica and one Uvicorn worker are required by the in-memory world and single-writer
  checkpoint. This setup does not provide horizontal scaling or independent world workers.

The route and access details are documented in the [API reference](api-reference.md).
The checkpoint and startup behavior are in [architecture](architecture.md#run-and-checkpoint).

## Configuration

Infrastructure settings belong in `.railway/railway.ts`; review them with `railway config plan`
and apply them with `railway config apply`. Secrets and environment-specific values are set in
Railway's variables, not committed to the repository:

- GitHub OAuth client ID and secret
- The numeric GitHub user ID allowed into the full app
- Session signing secret
- Public service URL

The Railway service mounts its persistent volume at `/data`. `RAILWAY_VOLUME_MOUNT_PATH` is
provided by Railway, and the app uses it for checkpoints. The service definition requires the
mount and its health check verifies that the checkpoint directory can be written. See
[checkpoint.py](../backend/app/checkpoint.py) and [architecture](architecture.md#run-and-checkpoint)
for the behavior.

## Hosting boundary

The Docker image keeps the application build portable, while the deployment contract is still
Railway-specific: one replica, one shared world, and one persistent checkpoint volume. Adding
accounts, paid access, per-user worlds, concurrent workers, a database or object store, or
sharding requires an explicit storage and execution design. The launch roadmap tracks those
decisions under [security and public service](../planning/launch-readiness.md#4-security-and-public-service)
and [performance and capacity](../planning/launch-readiness.md#5-performance-and-capacity).

Marketing and monetization notes are in [planning/marketing.md](../planning/marketing.md)
and [planning/monetization.md](../planning/monetization.md).

The current service has no published capacity or availability target. Set those targets and
choose which hosting responsibilities to keep before designing scale-out; see the launch
roadmap's [decisions to settle](../planning/launch-readiness.md#decisions-to-settle-before-calling-this-a-launch-plan).
