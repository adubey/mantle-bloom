// The hosted mantle-bloom (#376): one Railway service built from the root Dockerfile, serving
// the frontend and the API from one origin, with a volume for the world checkpoint.
// `railway config plan` / `railway config apply` from the repo root; see docs/hosting.md.
import { defineRailway, github, preserve, project, service, volume } from "railway/iac";

// checkpoint.py finds the volume through RAILWAY_VOLUME_MOUNT_PATH, which Railway sets to this.
const CHECKPOINT_MOUNT = "/data";

export default defineRailway(() => {
  const checkpoint = volume("mantle-bloom-checkpoint");

  const app = service("mantle-bloom", {
    source: github("adubey/mantle-bloom", { branch: "main" }),
    build: { builder: "DOCKERFILE", dockerfilePath: "Dockerfile" },
    deploy: {
      healthcheckPath: "/healthz",
      // Covers loading a large world checkpoint at startup, before /healthz can answer.
      healthcheckTimeout: 300,
      // The world lives in one process's memory and its checkpoint on one volume, which
      // Railway can't attach to more than one replica anyway.
      numReplicas: 1,
      // Fails a deploy with no volume here, rather than running one that forgets the world
      // whenever it sleeps.
      requiredMountPath: CHECKPOINT_MOUNT,
      // Sleep when idle. The checkpoint restores the world and its run status on the next
      // request; a running world is paused while asleep, not caught up.
      sleepApplication: true,
      restartPolicyType: "ON_FAILURE",
      restartPolicyMaxRetries: 5,
      // SIGTERM to SIGKILL. uvicorn closes viewers' long polls within 5 s (Dockerfile), then
      // waits 1 s for the run loop (run_control.WorldRunner.close). A step still under way
      // after that is cut off: the world restarts from the step before, whose checkpoint is
      // intact, and the cut-off write's temporary file is removed at startup.
      drainingSeconds: 10,
    },
    volumeMounts: { [CHECKPOINT_MOUNT]: checkpoint },
    // Secrets are set with `railway variables --set` (docs/hosting.md), never in this file.
    env: {
      // SQLite user database on the persistent volume (see #620). Not a secret, so it
      // lives here rather than in dashboard variables, derived from the mount path above.
      MANTLE_BLOOM_USER_DB_PATH: `${CHECKPOINT_MOUNT}/database.db`,
      MANTLE_BLOOM_GITHUB_CLIENT_ID: preserve(),
      MANTLE_BLOOM_GITHUB_CLIENT_SECRET: preserve(),
      MANTLE_BLOOM_GITHUB_OWNER_ID: preserve(),
      MANTLE_BLOOM_SESSION_SECRET: preserve(),
      MANTLE_BLOOM_PUBLIC_URL: preserve(),
      // Two worlds fit the container's RAM: about 135 MB steady state per default-detail
      // world on top of a ~170 MB baseline, with a fresh generate peaking near 480 MB (see
      // docs/hosting.md). Raise this only alongside the service's memory.
      MANTLE_BLOOM_MAX_WORLDS: "2",
    },
  });

  return project("mantle-bloom", { resources: [app, checkpoint] });
});
