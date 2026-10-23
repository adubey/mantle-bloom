"""A hand-built river whose lake spills over a rim several nodes away (issue #520), shared by
the hydrology, render and River Inspector tests.

Eight nodes on the equator, 0.01 rad apart, each a k-NN neighbour of the next. River
0 -> 1 -> sink 2, whose basin spans 2-4 behind rim 5; filled, sink 2 spills straight to node 6
beyond the rim (a non-neighbour hop, as `lakes.compute_spill_routing` routes it), then
6 -> ocean 7. Nodes 3 and 4 hold `lake_m` of water; the sink itself is shallower than the
visible depth unless `sink_m` says otherwise.
"""

import numpy as np

from app.hydroclimate import hydrology

N = 8
# Every edge of the dry spill path from sink 2 over rim 5, and on to the ocean.
DRY_SPILL_PATH = {(2, 3), (3, 4), (4, 5), (5, 6)}


def spill_hop_fields(lake_m: float, store_should_spill: bool = True, sink_m: float = 0.5) -> hydrology.HydrologyFields:
    points = np.array([[np.cos(0.01 * i), np.sin(0.01 * i), 0.0] for i in range(N)])
    neighbor_idx = np.array([[max(i - 1, 0), min(i + 1, N - 1)] for i in range(N)])
    lake_depth = np.zeros(N)
    lake_depth[2] = sink_m
    lake_depth[[3, 4]] = lake_m
    spill_target = np.full(N, -1, dtype=np.int64)
    spill_target[2] = 6
    should_spill = np.zeros(N, dtype=bool)
    should_spill[2] = True
    zeros = np.zeros(N)
    return hydrology.HydrologyFields(
        points=points,
        elevation=np.array([100.0, 80.0, 10.0, 20.0, 25.0, 40.0, 30.0, -10.0]),
        is_ocean=np.arange(N) == 7,
        neighbor_idx=neighbor_idx,
        flow_target=np.array([1, 2, 6, 2, 3, 6, 7, -1]),
        flow_accum=np.array([5.0, 6.0, 8.0, 0.0, 0.0, 0.0, 9.0, 0.0]),
        water_deposited=zeros.copy(),
        filled_elevation=zeros.copy(),
        spill_target=spill_target,
        is_river=np.isin(np.arange(N), [0, 1, 2, 6]),
        lake_depth=lake_depth,
        glacier_depth=zeros.copy(),
        plates_in_order=[],
        should_spill=should_spill if store_should_spill else None,
        channel_width=np.full(N, 5000.0),
        discharge_m3_s=np.array([500.0, 600.0, 800.0, 0.0, 0.0, 0.0, 900.0, 0.0]),
    )


def segment_pairs(fields: hydrology.HydrologyFields, sources) -> set[tuple[int, int]]:
    from_idx, to_idx, _, _ = hydrology.flow_segments(fields, np.asarray(sources))
    return set(zip(from_idx.tolist(), to_idx.tolist()))
