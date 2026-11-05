# #801: where the simulation spends its time

Measured with [`bin/debug/profile_phases.py`](../../bin/debug/profile_phases.py) at
`17baaca5` (`origin/main`), 6-core Apple Silicon, 8 GB, load 2-3 while running. Wall-clock
seconds, exclusive times that add up (see the script's docstring). Inputs are the UI defaults:
node density 4, climate density 4, fluid density 2, 100 kyr steps. Three runs, compact trees
in the JSON files beside this one.

| File | World | Scenarios |
|---|---|---|
| `default-seed0.json` | seed 0, generated, 123 K nodes, 18 plates | generate x2, render x2, step x3 |
| `default-seed0-cfd-step.json` | same, `wind_model="cfd"` | step x3 |
| `synthetic-seed0.json` | seed 0 with `synthetic_relief=True` (pseudo-history) | generate x2, render x2, step x3 |
| `seed443034896-67.4Myr.json` | save at 67.4 Myr, 128 K nodes, 8 plates | render x2, step x3 (peak RSS 1.9 GB) |

## Stepping

| Phase | Young world, diagnostic wind | Young world, CFD wind | 67 Myr world, diagnostic wind |
|---|---|---|---|
| **Step total** | **6.3 s** (6.2-6.6) | **8.8 s** | **14.8 s** (12.2-18.5) |
| Eroding and depositing sediment (all of `apply_erosion`) | 2.75 (43%) | 2.46 (28%) | 6.15 (41%) |
| &emsp;Resolving rivers and lakes | 1.28 (20%) | 1.07 (12%) | 4.60 (31%) |
| &emsp;&emsp;Resolving lakes | 0.19 (3%) | 0.19 (2%) | 2.31 (16%) |
| &emsp;&emsp;Marching glaciers | 0.42 (7%) | 0.36 (4%) | 1.51 (10%) |
| &emsp;Computing climate | 0.51 (8%) | 0.43 (5%) | 0.54 (4%) |
| &emsp;Routing flow downstream (erosion's own calls) | 0.40 (6%) | 0.40 (5%) | 0.46 (3%) |
| Advancing the wind and ocean solver | skipped | 2.87 (33%) | skipped |
| Updating intraplate faults | 1.03 (16%) | 0.99 (11%) | 1.87 (13%) |
| Solving plate motion (mantle flow) | 0.54 (9%) | 0.54 (6%) | 1.07 (7%) |
| Filling seafloor gaps | 0.53 (8%) | 0.49 (6%) | 1.48 (10%) |
| Deforming plates | 0.45 (7%) | 0.42 (5%) | 1.42 (10%) |
| Reading collision polarity | 0.33 (5%) | 0.33 (4%) | 1.25 (8%) |
| Moving magma | 0.07 (1%) | < 1% | 0.58 (4%) |
| Merging and splitting plates (finding rifts) | < 1% | < 1% | 0.27 (2%) |
| Raising mantle plumes (mantle blooms) | < 1% | < 1% | < 1% |
| Balancing sea level | 0.19 (3%) | 0.18 (2%) | 0.17 (1%) |
| Untraced | 2% | 1.5% | 1% |

Rows under 1% of a step are hidden by the script, so "< 1%" means "fired, but too small to list".
Rivers, lakes and glaciers are nested inside erosion, so they are part of its row, not in
addition to it.

## Building world (young world, 11.2 s cold and 11.6 s warm)

| Phase | Seconds | Share |
|---|---|---|
| Placing sutures and cratons | 7.5 | 67% |
| &emsp;Laying out continents and orogenic belts | 4.4 | 39% |
| &emsp;Partitioning the sphere into plates | 3.1 | 27% |
| Seeding rivers (initial hydrology, climate, glaciers) | 1.15 | 10% |
| Raising initial orogens | 0.71 | 6% |
| Shaping initial bathymetry | 0.54 | 5% |
| Computing climate (baseline) | 0.46 | 4% |
| Dating the initial rocks | 0.33 | 3% |

Cold and warm runs are within noise of each other: numba's on-disk cache was already warm.
The very first run in a fresh worktree (under cProfile, so not comparable) spent about 3 s in
numba compilation, so a first-ever start pays that once.

### With synthetic history (`--synthetic-relief`, `synthetic-seed0.json`)

Generation takes 16.0-16.3 s against 11.2-11.6 s. "Building pseudo-history" is 5.1-5.3 s
(32%) and replaces "Raising initial orogens" (0.71 s) and the lowland relief pass, so the net
cost is about 4.6 s. Everything else is unchanged: layout 4.5 s, partitioning 3.1 s, rivers
1.35 s. Stepping (6.1 s) and rendering (3.0 s cold, 2.0 s warm) match the default world, as
expected since synthetic history only changes the starting terrain.

## First map render after a step

2.9-3.0 s cold against 1.8-2.1 s warm, so the cache fill costs about 1 s: looking up the
nearest node for every pixel (0.52-0.56 s cold, 0.23-0.25 s warm) and resampling climate onto
the map grid (0.45-0.49 s, absent when warm). Biome classification costs 1.2-1.5 s on every
render, cached or not.

## What this says about the draft list in the issue

- **Big enough to show:** resolving rivers and lakes, marching glaciers, computing climate,
  updating faults, plate motion, filling gaps, deforming plates, reading collision polarity.
  Erosion's own work (0.5 s) is the only thing not covered by a nested label.
- **"Mantle Flow, Mantle Blooms"** are separate and both small: plate motion is 7-9% and
  plumes are under 1%. Splitting them is fine, but the bloom bar will flash by.
- **"Finding rifts"** is under 1% on a young world and 2% on the 67 Myr one. It will not be
  visible as its own trace.
- **"Finding stationary distribution for wind & rain"** is the CFD solver, 33% of a step with
  `wind_model="cfd"` and zero at the UI default. The humidity and rainfall sweeps in the
  diagnostic path are under 1%. The trace should appear only when CFD is on.
- **Step cost grows with world age** (6.3 s to 14.8 s here), and the growth is concentrated in
  lakes (12x), glaciers (3.6x) and collision polarity (3.8x). On the mature world, lakes and
  glaciers alone are over a quarter of the step.
- **"Building cache"** is not a state I could find in the UI or backend (the overlay headline is
  only "Generating world", "Advancing world" or "Recording animation"). The closest measured
  thing is the render-side fill above, about 1 s, which seems too short to deserve a bar.

## Caveats

- Three steps per run; the 67 Myr steps ranged 12.2-18.5 s, so the mature-world shares are
  indicative, not precise. Percentages of one run use that run's own total.
- One young and one mature seed. Step cost depends on plate count, ice cover and lake count,
  so other worlds will shift the shares.
- "Building pseudo-history" is off by default (`SYNTHETIC_RELIEF_DEFAULT = False`), so it is
  measured only in the `--synthetic-relief` run above.
- Generation was only measured on the young world; the 67 Myr save skips it.

## Caches after loading a save

Derived caches are not saved and are rebuilt lazily on first use (`persistence._drop_derived_caches`
clears the plate k-d trees, `node_kdtree_cache` and `node_position_tree_cache` on load;
`npz_save` omits the rest). Measured on a 6.2 Myr save: load 1.6 s; the first render right
after load 3.1 s against 1.9 s for the next, so the render-side fill costs about 1.2 s; step 1
after load 6.1 s against 5.5 s for steps 2-3, about 0.5 s extra. A step also drops the node
k-d trees itself (`step_world_progress`), so those are rebuilt every step regardless of loading.
That save had no climate, hydrology or CFD cache, so its first step recomputed them too.

## Refined hooks (`hooks-v2-*.json`)

The first set of files above used a coarse hook table. The `hooks-v2` files use the expanded
table in `profile_phases.py` (about 70 hooks, each labelled with its user-visible group), which
leaves 0.7-1.3% of a step untraced. Measured group shares of a step:

| Group | Young | Young, CFD wind | 67 Myr |
|---|---|---|---|
| Moving plates | 50.5% | 34.7% | 59.0% |
| Hydrometeorology (climate, water & ice) | 34.2% | 54.2% | 34.4% |
| Reshaping the surface | 12.7% | 9.3% | 5.3% |

Details and per-function numbers are in [traces.md](traces.md). The young and mature runs are
the same worlds as above; the young run uses synthetic history, which only changes generation.

Follow-up (`hooks-v3-seed443034896-67.4Myr.json`): hooking the neighbour search attributes the
mature world's unexplained plate costs to one function. Plate outlines
(`PlateWithSparseQuadPatch.outline_world`) take 2.92 s of a 14.9 s step (19.6%), shared across
gap fill, collision polarity, deformation and the plate-motion solve. After this, no hook sits on
more than 12% of a step on that world (lakes, 11.4%, is next). See [traces.md](traces.md).

## Reproducing

Run from the repository root with this worktree's `backend/.venv`, one heavy run at a time
(compute-budget skill).

| Files | Command |
|---|---|
| `hooks-v2-seed0-synthetic.json` | `bin/debug/profile_phases.py --synthetic-relief --generate-runs 1 --steps 3 --min-pct 0.7 --json <out>` (young world generated from seed 0 at the UI defaults; the group shares are regenerable to within timing noise) |
| `hooks-v2-seed0-cfd-step.json` | `bin/debug/profile_phases.py --wind-model cfd --scenarios step --steps 3 --min-pct 1.5 --json <out>` |
| `hooks-v2-seed443034896-67.4Myr.json`, `hooks-v3-...json` | `bin/debug/profile_phases.py --load <save> --scenarios step --steps 3 --warmup 1 --min-pct 1.5 --json <out>` |
| `default-seed0*.json`, `synthetic-seed0.json`, `seed443034896-67.4Myr.json` | the same script with the coarse hook table it had at the time (it has since gained hooks), so use the `hooks-v*` files for current numbers |

The mature runs load `mantle-bloom-seed443034896-67400000y.mbworld` (79 MB, SHA-256
`3d3000d872e68fee55645361b2686028c3faccd6edf03ce9d27e26ce7f286942`). It is a local save of a
user-run world and has not been uploaded to B2, so the mature numbers cannot be regenerated from
the repository alone; treat them as indicative, not as a baseline (3 steps per run, 12.2-19.0 s
between steps). Any save past ~60 Myr will show the same shape: lakes and glacier shaping grow
the most with age.

Step times quoted elsewhere come from different runs and so differ slightly: 6.06 s is the young
world with synthetic history (`hooks-v2-seed0-synthetic.json`), 6.3 s the earlier run with the
coarse table (`default-seed0.json`), 8.6 s with CFD wind, and 14.9 s the 67 Myr save.
