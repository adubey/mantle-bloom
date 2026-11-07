"""Long-march glacial carving (review of #806): the realistic 5 Myr default-density run
that the bounded 750 yr unit case in unit_tests/test_erosion.py stands in for. Too
slow for the PR gate (tens of minutes); nightly stress only. Logs substeps, wall time
and march phase timings with the result."""

import time

import numpy as np

from app import trace
from app.hydroclimate import erosion, ice_sheet
from app.world import generate_world


def test_glacial_erosion_carves_channels_over_five_myr():
    world = generate_world(seed=3, num_plates=8)
    world.rain_erosion_multiplier = 0.0
    world.river_erosion_multiplier = 0.0
    for plate in world.plates:
        plate.set_fields_on_plate(glacier_depth=np.where(plate.collect("elevation") > world.sea_level_m, 500.0, 0.0))
    _, _, before, _, _, _ = erosion._gather_nodes(world)
    ice_sheet.MARCH_LOG.clear()
    with trace.recording() as rec:
        t0 = time.time()
        erosion.apply_erosion(world, years=5_000_000)
        wall_s = time.time() - t0
    _, _, after, _, _, _ = erosion._gather_nodes(world)
    carved = float(np.max(after - before))
    entry = ice_sheet.MARCH_LOG[-1]
    print(
        f"5 Myr march: {entry['substeps']} substeps, settled={entry['settled']}, "
        f"wall_s={wall_s:.1f}, carved_m={carved:.3f}"
    )
    for path, (incl, _excl, calls) in sorted(rec.snapshot().items(), key=lambda kv: -kv[1][0]):
        label = " > ".join(path)
        if "march" in label or "Marching" in label:
            print(f"  {incl:.1f}s x{calls}  {label}")
    # No river erosion at all, so all new channel is a glacial trough.
    assert carved > 0.1
