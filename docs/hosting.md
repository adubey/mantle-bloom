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
- The owner keeps several named worlds in one process, each addressed by world ID; the
  public viewer opens whichever ones the owner marks public. There are no separate
  user-owned worlds or general user accounts.
- A Railway volume stores the registry index plus every world's checkpoint and run status
  (`worlds.json` and `worlds/<id>/` per world; a default world pickles to about 60 MB).
  On restart or wake, the process restores all of that state. If the service sleeps,
  simulation pauses until a request wakes it; elapsed time is the time simulated, not
  wall-clock time.
- One replica and one Uvicorn worker are required by the in-memory worlds and single-writer
  checkpoints. This setup does not provide horizontal scaling or parallel physics workers:
  independently running worlds take turns at one compute gate. `MANTLE_BLOOM_MAX_WORLDS`
  caps the registry (default 4 records); size the volume and memory for the worlds kept.

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
- User database file (`MANTLE_BLOOM_USER_DB_PATH`, a SQLite file on the volume;
  unset means the in-memory store, the desktop default -- see #620)

### Reaching the admin dashboard

The full app at `/` needs the owner's session. Behind Railway's proxy the
local-request bypass never applies (forwarding headers mark the traffic as relayed),
so GitHub OAuth is the only remote admin path, and without it `/view` is all anyone
gets. Checklist when the dashboard is unreachable (#626):

1. Create a GitHub OAuth app whose authorization callback URL is
   `https://<service-domain>/auth/github/callback`.
2. Set all five Railway variables: `MANTLE_BLOOM_GITHUB_CLIENT_ID`,
   `MANTLE_BLOOM_GITHUB_CLIENT_SECRET`, `MANTLE_BLOOM_GITHUB_OWNER_ID` (the numeric
   user ID, e.g. `gh api user -q .id`), `MANTLE_BLOOM_SESSION_SECRET` (32+ random
   bytes), and `MANTLE_BLOOM_PUBLIC_URL` (the public `https://` URL -- the OAuth
   callback and the cross-site checks are derived from it).
3. Redeploy, open `/view`, and use its Sign in entry, which lands on `/login` (with
   setup help while login is still unconfigured) and then GitHub.

Rotating the session secret signs everyone out; changing the owner ID invalidates
existing sessions.

The Railway service mounts its persistent volume at `/data`. `RAILWAY_VOLUME_MOUNT_PATH` is
provided by Railway, and the app uses it for checkpoints. The service definition requires the
mount and its health check verifies that the checkpoint directory can be written. See
[checkpoint.py](../backend/app/checkpoint.py) and [architecture](architecture.md#run-and-checkpoint)
for the behavior.

## Hosting boundary

The Docker image keeps the application build portable, while the deployment contract is still
Railway-specific: one replica, one multi-world registry, and one persistent checkpoint volume.
Adding accounts, paid access, per-user worlds, parallel physics workers, a database or object
store, or sharding requires an explicit storage and execution design. The launch roadmap
tracks those decisions under [security and public
service](../planning/launch-readiness.md#4-security-and-public-service) and [performance and
capacity](../planning/launch-readiness.md#5-performance-and-capacity).

Marketing and monetization notes are in [planning/marketing.md](../planning/marketing.md)
and [planning/monetization.md](../planning/monetization.md).

The current service has no published capacity or availability target. Set those targets and
choose which hosting responsibilities to keep before designing scale-out; see the launch
roadmap's [decisions to settle](../planning/launch-readiness.md#decisions-to-settle-before-calling-this-a-launch-plan).
