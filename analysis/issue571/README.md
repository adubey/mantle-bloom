# Issue 571: sea level falls without the ice changing

On `origin/main`, seed 997271774 replayed from 87.4 Myr loses 117 m of sea level in 10 Myr
while its ice barely changes (#554's `main` arm). Spread over the ocean, that is about 40 M km³
of water.

**Result: the seafloor sinks.** Water isn't moving onto land, and #550 isn't involved. Over
2 Myr:
- the ocean basins gain 32 m of sea level's worth of capacity below the waterline, at constant
  area;
- oceanic mantle lithosphere thickens by 1.6 km per Myr while oceanic crust stays the same, so
  the seafloor subsides isostatically by about 32 m per Myr;
- the commit before #550 gives the same result to within 0.5 m.

## Method

[`measure_sea_level_budget.py`](../../bin/debug/measure_sea_level_budget.py) replays a save and
splits the change in sea level into two parts:
- **water held on land:** ice as water, plus lakes at least 1 m deep (eustasy's own debit);
- **basin capacity:** how much water the connected ocean holds below the *opening* sea level
  (`eustasy.ocean_water_volume_m3` at that fixed level).

Each checkpoint also records:
- area and mean depth below the opening level;
- for oceanic and continental columns (by node crust density): area, mean elevation, crust
  thickness (Hc) and mantle-lithosphere thickness (Hm).

Seed 997271774 at 87.4 Myr (`08ad9461…ed762b`), 2 Myr in 0.25 Myr steps, run one at a time:
- `origin/main` at `8ed179f`;
- the commit before #550, `009d63a`.

Every checkpoint is in [`sea-level-budget.json`](sea-level-budget.json).

## Results

| After 2 Myr | `origin/main` | Before #550 |
|---|---:|---:|
| Sea-level change | −32.0 m | −31.9 m |
| …from water held on land | +0.1 m | +0.7 m |
| …from basin capacity | −32.0 m | −32.5 m |
| Trapped lake water, range | 0.09–0.32 M km³ | 0.03–0.32 M km³ |
| Water budget | constant | constant |

Lake water moves up and down but doesn't accumulate. Its whole range is under 1 m of sea
level.

| `origin/main`, 87.4 → 89.4 Myr | Opening | After 2 Myr |
|---|---:|---:|
| Area below the opening sea level | 362.6 M km² | 363.0 M km² |
| Mean depth below it | 3,981 m | 4,008 m |
| Oceanic columns: area | 242.8 M km² | 242.4 M km² |
| …mean elevation | −4,803 m | −4,867 m |
| …crust thickness (Hc) | 10.23 km | 10.22 km |
| …mantle lithosphere (Hm) | 21.1 km | 24.3 km |
| Continental columns: mean elevation | −919 m | −910 m |
| …Hc, Hm | 32.0, 107.8 km | 32.0, 107.8 km |

The basins deepen at constant area. Oceanic Hm grows steadily, about 0.4 km per 0.25 Myr
step, and the oceanic surface sinks with it. Continental columns don't change.

## Interpretation

Seafloor sinking as its mantle lithosphere cools and thickens is real physics, but here it is
fast and one-sided:
- **Rate.** About 16 m of sea level per Myr, 5–8× the fastest rates Earth gets from changes in
  ridge volume (a few metres per Myr).
- **No balance.** On Earth, new high-standing seafloor at ridges roughly balances old seafloor
  sinking, so mean ocean depth changes slowly.
- **A young ocean.** A mean oceanic Hm of 21 km corresponds to a few Myr of half-space cooling.
  Earth's mean seafloor age is about 60 Myr, with mantle lithosphere near 100 km.

Either too little new seafloor is being made to balance the cooling, or oceanic Hm grows faster
than its age implies. This is lithosphere and tectonics behaviour, not RFC 0014's lake water or
RFC 0015's ice.

**Principle:** #6, account for material through every process. The water is accounted for;
the basin geometry is what moves.

**Next:** measure oceanic Hm growth against the half-space cooling law, and the rate of new
seafloor production against the area that subsides, on more than one save. Seed 649658380
falls only 14 m in 10 Myr. The fix belongs in the lithosphere model.

## Reproduce

```bash
cd backend
.venv/bin/python ../bin/debug/measure_sea_level_budget.py <save> --myr 2 --step-years 250000 \
    --checkpoint-myr 0.25 --out <scratch>/main.json
.venv/bin/python ../bin/debug/measure_sea_level_budget.py <save> --myr 2 --step-years 250000 \
    --checkpoint-myr 0.25 --backend <pre-550 tree>/backend --out <scratch>/before.json
```
