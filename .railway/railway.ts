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
      // The public process checks the checkpoint volume without loading worlds.
      healthcheckTimeout: 30,
      // The world lives in one process's memory and its checkpoint on one volume, which
      // Railway can't attach to more than one replica anyway.
      numReplicas: 1,
      // Fails a deploy with no volume here, rather than running one that forgets the world
      // whenever it sleeps.
      requiredMountPath: CHECKPOINT_MOUNT,
      // The app retires its simulation worker on idle and keeps the public snapshot
      // available. Platform sleep would also stop a local-only Play session.
      sleepApplication: false,
      restartPolicyType: "ON_FAILURE",
      restartPolicyMaxRetries: 5,
      // The public process drains requests for 2 s, then allows its worker up to 7 s
      // to exit. The last committed checkpoint survives an interrupted write.
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
      // Stop-gap for the hosted world's oversized lakes, not a physical model: let a pit
      // that drained to the ocean last step keep its outlet open (#867's outlet-preservation
      // bound), with the notch and channel-depth ceilings raised from 2,000 m to 20 km. In
      // #867 this cut lake volume by about a third, but lakes stayed ~15x Earth's; the cause
      // is the water balance. Remove both once that is fixed (#882).
      MANTLE_BLOOM_RIVER_ARM: "outlet-bound",
      MANTLE_BLOOM_DEPTH_CAP_M: "20000",
    },
  });

  return project("mantle-bloom", { resources: [app, checkpoint] });
});
