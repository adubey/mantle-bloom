// The hosted mantle-bloom (#376): one Railway service built from the root Dockerfile, serving
// the frontend and the API from one origin, with a volume for the world checkpoint.
// `railway config plan` / `railway config apply` from the repo root; see docs/HOSTING.md.
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
      // Time for uvicorn's 5 s graceful shutdown (Dockerfile) before the old deploy is killed.
      drainingSeconds: 10,
    },
    volumeMounts: { [CHECKPOINT_MOUNT]: checkpoint },
    // Secrets are set with `railway variables --set` (docs/HOSTING.md), never in this file.
    env: {
      MANTLE_BLOOM_GITHUB_CLIENT_ID: preserve(),
      MANTLE_BLOOM_GITHUB_CLIENT_SECRET: preserve(),
      MANTLE_BLOOM_GITHUB_OWNER_ID: preserve(),
      MANTLE_BLOOM_SESSION_SECRET: preserve(),
      MANTLE_BLOOM_PUBLIC_URL: preserve(),
    },
  });

  return project("mantle-bloom", { resources: [app, checkpoint] });
});
