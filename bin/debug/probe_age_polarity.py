#!/usr/bin/env python3
"""Why ocean-ocean age verdicts stay undecided and why they reverse (RFC 0010, #737). Read-only.

Usage, one process per seed:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_age_polarity.py <seed> <out.json>
"""
import json, sys
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0, 'bin/debug')
import probe_rift_settings as p

seed, out = int(sys.argv[1]), sys.argv[2]
undecided, rev = Counter(), []
frozen = {}
original_fronts = p.polarity_fronts

def fronts(world, plate, state, by):
    rows = original_fronts(world, plate, state, by)
    now = world.elapsed_years / 1e6
    for r in rows:
        other = by.get(r['neighbour'])
        own = plate.collect('node_created_years')[r['nodes']] / 1e6
        oth = other.collect('node_created_years')[state['inputs'].neighbor_node_index[r['nodes']]] / 1e6 if other else own
        r['stats'] = dict(own_med=now - max(np.median(own), 0), oth_med=now - max(np.median(oth), 0),
                          own_start=float((own <= 0).mean()), oth_start=float((oth <= 0).mean()))
    return rows
p.polarity_fronts = fronts

original_observe = p.Census.observe_polarity
def observe(self, rows, now, by):
    real = self.polarity_tracker.update
    captured = {}
    def update(*a):
        result = real(*a)
        captured.update(result[1])
        return result
    self.polarity_tracker.update = update
    original_observe(self, rows, now, by)
    self.polarity_tracker.update = real
    for i, track in captured.items():
        plate, state, row = rows[i]
        if not row['counted'] or row['kind'] != 'ocean_ocean':
            continue
        s = row['stats']
        gap = s['own_med'] - s['oth_med']
        local = p.geometry.to_local(plate.frame, state['points'][row['nodes']])
        if track.verdict is None:
            both = s['own_med'] >= now - 1e-6 and s['oth_med'] >= now - 1e-6
            undecided['both_medians_start' if both else 'gap_under_10_with_younger_cells'] += 1
            undecided['myr_0_50' if now <= 50 else 'myr_50_100' if now <= 100 else 'myr_100_200'] += 1
        elif track.track_id not in frozen:
            frozen[track.track_id] = dict(at=now, gap=gap, local=local, own=s['own_med'], oth=s['oth_med'], done=False)
        f = frozen.get(track.track_id)
        if f and not f['done'] and row['age'] not in ('tied', track.verdict):
            f['done'] = True
            d = cKDTree(f['local']).query(local)[0]
            rev.append(dict(frozen_myr=f['at'], reversed_myr=now, gap_then=round(f['gap'], 1), gap_now=round(gap, 1),
                            own_shift=round(s['own_med'] - f['own'] - (now - f['at']), 1),
                            other_shift=round(s['oth_med'] - f['oth'] - (now - f['at']), 1),
                            same_nodes_share=round(float((d <= self.spacing).mean()), 2),
                            front_nodes=len(row['nodes']), frozen_nodes=len(f['local'])))
p.Census.observe_polarity = observe
report = p.run_seed(seed, 200.0, 1.0)
json.dump(dict(seed=seed, fingerprint=report['fingerprint'], undecided=undecided, reversals=rev, frozen=len(frozen)), open(out, 'w'))
