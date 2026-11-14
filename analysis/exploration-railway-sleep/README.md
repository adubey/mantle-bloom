# Railway idle memory and map startup

Measured November 13, 2026. Production stayed near **3.02 GB (2.81 GiB) of reported memory** with almost no CPU use. A local reproduction on the deployed commit went from **424 MiB after checkpoint restore** to **2,482 MiB after the first public render**. After an admin open, releasing the world and clearing the render grids left **2,526 MiB** resident. A fresh HTTP process served the saved image in about **3 ms** at **57–64 MiB RSS**. Exiting the process that ran the render is the strongest memory-reclamation option these measurements support.

## Production observations

- Railway production is running commit `56c7159156beb5f512aa12db1ed55a78362aa214`. The service instance reports `sleepApplication: false`. The repository has `sleepApplication: true` in [.railway/railway.ts](../../.railway/railway.ts). Railway applies a Serverless setting change when a new container is deployed. Sleeping containers incur no compute usage. [Railway Serverless](https://docs.railway.com/deployments/serverless), [Railway deployment behavior](https://docs.railway.com/platform/advanced-concepts#deployments)
- Over seven days, 169 hourly buckets reported 3.016–3.019 GB in the raw `MEMORY_USAGE_GB` series and 0.00084–0.00226 vCPU. A 24-hour edge-log window had 180 requests, no step, run, or animation POSTs, and a maximum request gap of about 100 minutes. Over the observed 23.7 hours, 67–82% of time lay beyond 10- and 5-minute gaps respectively. This is an HTTP-only estimate of time eligible for sleep. Railway considers all outbound traffic and samples inactivity. An open viewer's world-list polling and state long polling can keep the deployment active. [Railway Serverless](https://docs.railway.com/deployments/serverless)
- One already-awake public request sequence took 259 ms for the HTML, 147 ms for the list, 136 ms for state, and 407 ms for the 1.60 MB render response. The sum is about 0.95 seconds for those serial requests on this client. This sample excludes frontend assets, browser decoding and painting, and Railway cold boot.
- The public listing showed one stopped default world at 6.3 Myr. Railway SSH has no registered key for this account, so the production process RSS, cgroup file-cache split, and private checkpoint count could not be inspected.

## Local profile

Used the deployed backend source in an isolated scratch directory on a 16 GB Apple Silicon Mac, Python 3.14.6 and pinned NumPy, SciPy, Numba, FastAPI, and Uvicorn. FastAPI's in-process HTTP client exercised the real routes. The fixture was an existing v5 save of the same public world at 6.2 Myr, one 100 kyr step earlier, with point, climate, and fluid densities 4, 4, and 2 and diagnostic wind. Its SHA-256 is `87359964b55a8e82bab47131ab2bd9ae26d87b9e9ab7963dd6541ae2187f3992`. The route run shared the host with one probe and two test workers. These are single-run macOS RSS and timing observations. Railway reports container memory, and Linux allocator and filesystem-cache behavior can differ.

| Phase | Local wall time | Process RSS after phase |
|---|---:|---:|
| FastAPI and HTTP test client | 0.17–0.31 s | 49 MiB |
| Import deployed backend after HTTP imports, first run | 15.1 s | 225 MiB |
| Import deployed backend in two later fresh Python processes with files already read | 1.38–1.49 s | 215 MiB in the route run |
| Restore one world through FastAPI lifespan | 0.29–0.32 s | 415–424 MiB |
| First public Combined render, 2200 × 1222, cache miss | 3.16 s | 2,482 MiB |
| Same public render, RAM hit / disk hit | 2.1 / 3.5 ms | 2,482 MiB |
| Admin list, summary, and stats history | 1.6–3.4 ms each | 2,482 MiB |
| Admin opens controls via empty POST | 2.29 s | 2,746 MiB |
| Admin Combined render | 1.57 s | 2,613 MiB |
| Admin stats | 0.19 s | 2,601 MiB |
| Release world, clear public render cache, and collect garbage | | 2,556 MiB |
| Also clear projected-grid and biome-grid caches | | 2,526 MiB |
| Fresh HTTP process serves saved JSON or raw PNG | 2.8–3.9 ms | 57–64 MiB |
| First 100 kyr guest step, physics in a fresh worker | 7.24 s | 912 MiB |
| Write the new checkpoint after the guest step | 0.19 s | 1,031 MiB |

The local process high-water mark was 2,559 MiB on the first render and 3,169 MiB during admin open. A weak reference confirmed the world object was collected. The remaining RSS includes allocations outside the world and the cleared caches. This run did not attribute retained bytes among native allocators, compiled code, and other module state. The spread between local imports warrants direct measurement on Railway before setting a wake-up expectation.

The step run used a new process and the original 6.2 Myr checkpoint copy. It omitted the automatic preview render to measure physics and checkpoint writing independently. The real guest route returned `computed: true` at 6.3 Myr and wrote a 127 MB checkpoint. Its measured 7.47-second request excludes import, restore, and preview generation. The deployed path also regenerates the default preview after a step. The step run retained 1,026 MiB after releasing the world.

The save contained roughly 48 MiB of hydrology entries, 35 MiB of plate entries, and 30 MiB of climate entries before compression. The Combined renderer resamples onto an 801 × 1601 internal grid and builds multiple fields, projected geometry, and color arrays before encoding. The biome blender also creates a full grid of float64 weights per class and copies the weight array to select the runner-up. Those allocations are candidates for a separate render-peak profile. The rendered PNG was 1.22 MB, compared with 1.63 MB in its JSON/base64 wrapper. Serving the original PNG bytes preserves the alpha channel that the Combined and Biome maps use for interactive class highlights. The test client returned identical bytes on disk and RAM render hits.

## Sleep and wake behavior

Use a lightweight HTTP process to serve the frontend and a revisioned public snapshot, with the simulation in a child process that can exit after its idle timeout. In the local profile the post-admin process and snapshot server differed by about **2.48 GiB RSS**. A **100 MiB** snapshot-only service would save about **2.7 GiB** against the current Railway baseline if production reproduces that target. A Railway Serverless sleep can remove the entire running compute charge for that interval. Use an edge-served viewer shell and snapshot if the first map must also stay fast while Railway itself has stopped the HTTP process.

- **New `/view` navigation.** Return the last committed map and its metadata without waiting on simulation imports. The current viewer discovers the world, fetches state, and then requests the image. A snapshot URL in initial state or the page response lets the browser request the image earlier. After serving the browser's initial snapshot, start one shared warm-up for the default or explicitly selected world. Load the full checkpoint into RAM and reuse the saved PNG. Route asset fetches, list polling, and `after=` long polls through the lightweight process so routine polling does not rehydrate the worker. Persist a canonical preview for each public world that needs this path. The existing disk cache covers the default world. A snapshot miss or request for an uncached render can wait on the same warm-up.
- **Authenticated admin.** Start or join warm-up when a world is selected, then restore controls, statistics, and the requested admin map. The current `POST /worlds/{id}/controls` with `{}` recomputes climate, publishes a public revision, renders the default preview and writes the checkpoint. Give admin open a read-only controls path and reuse the existing map only when world revision, view, rotation, and size match.
- **Step, Play, and animation.** Require the restored simulation and use the existing compute gate. Keep the worker up while a run or write is active. If Play must continue with no viewers, the hosting policy also needs to keep the container running during Play. Railway's traffic-based sleep can pause a local-only run. Preserve the default checkpoint's climate, hydrology, CFD state and collision polarity data, which the next step can read before refreshing them. Rebuild derived spatial indexes and load or compile Numba kernels on demand. Report waking and stepping as distinct phases, with guest allowance charged only for an actual computation.
- **Return to sleep.** Publish snapshot metadata and the PNG under the same checkpoint generation, using an atomic pointer or equivalent commit record. Gate access on current visibility, including unpublish and delete. The current preview key uses seed, age, view, size and rotation, and the image can be written before the checkpoint succeeds. A durable revision link is needed for snapshot serving while the worker is down. Keep the worker resident if checkpointing failed or if a step, run, animation or write is active. Deduplicate overlapping wake requests and keep a single checkpoint writer.

Production validation should record cold first paint and peak/idle container memory over repeated sleep cycles, with one and two worlds. Include simultaneous viewer/admin arrivals, failed checkpoint writes, unpublish, and a warm-up interrupted by a restart. Verify that a resumed step from the same committed checkpoint matches uninterrupted execution under the existing save-continuation contract.
