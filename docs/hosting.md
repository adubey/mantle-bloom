# Hosting: current setup and next decisions

This page records Mantle-Bloom's current hosted shape and the boundary of that design. The
full launch goals and open product decisions are in [launch readiness](../planning/launch-readiness.md).
Local development and localhost sign-in are covered by [local hosting](local-hosting.md).

## Current deployment

The public service runs on Railway from the repository's root [Dockerfile](../Dockerfile),
configured as code in [.railway/railway.ts](../.railway/railway.ts). The multi-stage image
builds the React frontend and installs the pinned Python runtime. One public Uvicorn process
serves the frontend and API from the same origin and starts a private simulation worker on demand. Railway builds from `main` and uses `GET /healthz` as its
health check. The Docker image is portable; Railway is the only deployment configured in
this repository.

The service is a small hosted preview, not a multi-user world service:

- One configured GitHub account has the full app and API.
- Logged-out visitors can use the public viewer and its limited step allowance.
- The owner keeps several named worlds in one process, each addressed by world ID; the
  public viewer opens whichever ones the owner marks public. There are no separate
  user-owned worlds or general user accounts. On Railway the registry holds at most two
  worlds (`MANTLE_BLOOM_MAX_WORLDS: "2"` in `.railway/railway.ts`); the app default is 4.
- A Railway volume stores the registry index, world checkpoints, run status, and public
  previews. One simulation worker owns checkpoint writes. It serializes compute across
  worlds, so independently running worlds take turns at the same compute gate.
- `MANTLE_BLOOM_MAX_WORLDS` caps the registry (default 4 records). Size memory for active
  simulation and rendering at the detail you use.

## Idle worker and memory

On `/view`, the public process reads the last committed map from the volume. After painting
the map, the browser asks the worker to restore the selected world in the background.
Admin requests for a world and authorized step requests start or join the same worker.
The first open of a world saved before snapshot support may need to restore and render it.
Public warm hints share a cooldown of twice the idle timeout, starting when a hint
finishes. Explicit admin and step requests can wake the worker immediately.

`MANTLE_BLOOM_IDLE_SECONDS` sets how long the worker stays up after its last active request
(default 90 seconds). Routine public polling can continue through sleep. Play, requests
with active simulation work, and unsaved checkpoint changes keep the worker alive.
If the worker exits unexpectedly, the public process restarts it for any checkpointed
Play session. A failed hosted step returns the world to its last committed checkpoint.
The public list omits worlds beyond the registry cap. The admin list shows their load error.
Every completed change to a public world can trigger a canonical 2200 × 1222 render,
including worlds other than the default. Allow memory for that render when raising the cap.
Railway Serverless is disabled so the public process can always serve the saved map and
Play can continue without traffic. The [memory profile](../analysis/exploration-railway-sleep/README.md)
records the original idle footprint and the local measurements that informed this design.

The hosted entry point requires a signed owner session for admin routes, including when
run on loopback. Use the desktop or development entry point for local admin access without OAuth.
Opening the API schema starts the worker on the first request so the hosted docs include
the admin API.

To avoid generating on the host at all, generate locally and upload the result: File > Save
World on the local app, then File > Load World on the hosted app (admin only, like every
`/worlds/*` route). The upload streams to a temp file and unpickles from disk rather than
holding the whole body in RAM, and bodies past 512 MB (`MAX_LOAD_WORLD_BYTES` in
`backend/app/main.py`) are refused with `413` while still streaming.

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

### Hosted hydrology stop-gap

The hosted service doesn't run the production lake and river model. `.railway/railway.ts` sets
`MANTLE_BLOOM_RIVER_ARM=outlet-bound` and `MANTLE_BLOOM_DEPTH_CAP_M=20000`, two experiment
knobs from #867 (see [hydrology](simulation-model/hydrology.md)). The first lets a basin whose
water reached the ocean on the previous step cut its outlet open. The second raises the
per-node notch and channel-depth ceilings from 2,000 m to 20 km. Neither is physical: both
hold back the oversized inland water the production model makes, and that excess comes from
the water balance, not incision. Lake outcomes on the hosted world therefore differ from local
runs, which leave both variables unset. [#882](https://github.com/adubey/mantle-bloom/issues/882)
tracks removing them once lakes are fixed. Both are read once at startup, so a change takes
effect on the next deploy.

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
