"""Print G1 scores side by side (#975).

Usage: python3 bin/debug/tabulate_g1_scores.py analysis/issue975/scores b123-brake-off b123 [<setting> ...]
Each <setting>.json is the output of `calibrate_plate_forces.py score`.
"""

import json, sys
S = sys.argv[1]; names = sys.argv[2:]
def g(d, *ks):
    for k in ks:
        d = d.get(k) if isinstance(d, dict) else None
    return d
rows = [("T1 NNR cm/yr", "T1_speed_nnr_cm_per_yr", "{:.2f}"), ("T2 ratio (+10 Myr)", "T2_slab_attached_ratio", "{:.2f}"),
        ("T2 same time", "T2_same_time_ratio", "{:.2f}"), ("T3 pull share", "T3_slab_pull_driving_share", "{:.1%}"),
        ("T4 MPa", "T4_slabless_basal_traction_mpa", "{:.2f}"), ("T5 N/m", "T5_old_ridge_push_n_per_m", "{:.2e}"),
        ("T7 fast/slow world", "T7_fastest_over_slowest", "{:.2f}")]
data = {n: json.load(open(f"{S}/{n}.json")) for n in names}
print("measure".ljust(26), *(n[:14].rjust(15) for n in names))
for label, key, fmt in rows:
    print(label.ljust(26), *(fmt.format(data[n][key]).rjust(15) if data[n].get(key) is not None else "-".rjust(15) for n in names))
extra = [("net slab work", ("breakdown","slab","power_over_pull_toward_trench","net"), "{:+.1%}"),
         ("drag normal power", ("breakdown","slab","power_over_pull_toward_trench","drag_trench_normal"), "{:+.1%}"),
         ("drag parallel power", ("breakdown","slab","power_over_pull_toward_trench","drag_trench_parallel"), "{:+.1%}"),
         ("pull N/m", ("breakdown","slab","pull_n_per_m"), "{:.2e}"), ("pull column km", ("breakdown","slab","hm_km"), "{:.1f}"),
         ("closing pulled", ("breakdown","closing_down_going_length_pulled"), "{:.1%}"),
         ("ridge push drive share", ("breakdown","driving_share","ridge_push"), "{:.1%}"),
         ("slab drag resist share", ("breakdown","resisting_share","slab_drag"), "{:.1%}"),
         ("splits /world/100Myr", ("T9_rifting","splits_per_world_100myr"), "{:.1f}"),
         ("T8 jammed drops", ("T8_contacts","collision_time_drops_jammed"), "{}"),
         ("T8 largest plate share", ("T8_contacts","largest_plate_share_mean"), "{:.2f}"),
         ("T8 plate count", ("T8_contacts","plate_count_mean"), "{:.1f}"),
         ("solve residual max", ("closure","residual_over_driving_max"), "{:.1e}")]
for label, ks, fmt in extra:
    print(label.ljust(26), *(fmt.format(g(data[n], *ks)).rjust(15) if g(data[n], *ks) is not None else "-".rjust(15) for n in names))
