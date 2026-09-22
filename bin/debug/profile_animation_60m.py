"""Profile a 60-frame, 60-Myr animation and emit per-frame high-level timings."""
from __future__ import annotations

import cProfile
import csv
import io
import json
import pstats
from collections import defaultdict
from pathlib import Path
from time import perf_counter

import av
import numpy as np
from PIL import Image

from app import erosion, eustasy, faults, geology, merge_split, render_image, volcanism, world as world_mod
from app.lithosphere_plate import LithospherePlate

FRAMES = 60
STEPS_PER_FRAME = 10
STEP_YEARS = 100_000
OUT = Path("analysis/issue147-profile-20260922")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frame_bucket: dict[str, float] = defaultdict(float)
    current_frame = [0]

    def wrap(obj, name: str, label: str):
        original = getattr(obj, name)
        def timed(*args, **kwargs):
            started = perf_counter()
            try:
                return original(*args, **kwargs)
            finally:
                frame_bucket[label] += perf_counter() - started
        setattr(obj, name, timed)

    # These wrappers match step_world's major sequential phases. shift/deform are aggregated
    # over all plates; the remaining step overhead is reported explicitly as unclassified.
    wrap(LithospherePlate, "shift", "plate_shift")
    wrap(LithospherePlate, "deform", "plate_deform")
    wrap(faults, "update_faults", "faults")
    wrap(merge_split, "apply_topology_changes", "topology_changes")
    wrap(erosion, "apply_erosion", "climate_erosion_hydrology")
    wrap(volcanism, "apply_volcanic_activity", "volcanism")
    wrap(geology, "apply_resource_formation", "resource_formation")
    wrap(eustasy, "update_sea_level", "sea_level")

    metadata = {
        "frames": FRAMES, "steps_per_frame": STEPS_PER_FRAME,
        "step_years": STEP_YEARS, "total_years": FRAMES * STEPS_PER_FRAME * STEP_YEARS,
        "seed": 0, "node_density": 4.0, "climate_density": 4.0,
        "fluid_density": 2.0, "projection": "eckert4", "view": "combined",
        "width": 2200, "height": 1222,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    t0 = perf_counter()
    sim_world = world_mod.generate_world(seed=0, node_density=4.0, climate_density=4.0, fluid_density=2.0)
    generation_s = perf_counter() - t0
    metadata.update(generation_s=generation_s, initial_nodes=sum(p.node_count() for p in sim_world.plates), initial_plates=len(sim_world.plates))
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"generated in {generation_s:.3f}s: {metadata['initial_nodes']} nodes, {metadata['initial_plates']} plates", flush=True)

    mp4_buf = io.BytesIO()
    container = av.open(mp4_buf, mode="w", format="mp4")
    stream = container.add_stream("libx264", rate=render_image.ANIMATION_FPS)
    stream.width, stream.height, stream.pix_fmt = 2200, 1222, "yuv420p"
    stream.options = {"crf": render_image.ANIMATION_CRF, "movflags": "+faststart"}

    fields = ["frame", "start_myr", "end_myr", "nodes", "plates", "step_total", "plate_shift", "plate_deform", "faults", "topology_changes", "climate_erosion_hydrology", "volcanism", "resource_formation", "sea_level", "step_other", "render", "decode_encode", "frame_total"]
    rows = []
    profiler = cProfile.Profile()
    profiler.enable()
    run_started = perf_counter()
    for frame in range(1, FRAMES + 1):
        current_frame[0] = frame
        frame_bucket.clear()
        frame_started = perf_counter()
        start_years = sim_world.elapsed_years
        step_started = perf_counter()
        for _ in range(STEPS_PER_FRAME):
            world_mod.step_world(sim_world, STEP_YEARS)
        step_total = perf_counter() - step_started
        render_started = perf_counter()
        png = render_image.render_png(sim_world, "eckert4", "combined", 2200, 1222, None)
        render_s = perf_counter() - render_started
        encode_started = perf_counter()
        image = Image.open(io.BytesIO(png)).convert("RGB")
        video_frame = av.VideoFrame.from_ndarray(np.asarray(image), format="rgb24")
        for packet in stream.encode(video_frame):
            container.mux(packet)
        encode_s = perf_counter() - encode_started
        measured = sum(frame_bucket.values())
        row = {
            "frame": frame, "start_myr": start_years / 1e6, "end_myr": sim_world.elapsed_years / 1e6,
            "nodes": sum(p.node_count() for p in sim_world.plates), "plates": len(sim_world.plates),
            "step_total": step_total, **{k: frame_bucket.get(k, 0.0) for k in fields[6:14]},
            "step_other": max(0.0, step_total - measured), "render": render_s,
            "decode_encode": encode_s, "frame_total": perf_counter() - frame_started,
        }
        rows.append(row)
        with (OUT / "frames.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
        print("frame", frame, json.dumps({k: round(row[k], 3) for k in fields[4:] if isinstance(row[k], float)}), flush=True)

    for packet in stream.encode():
        container.mux(packet)
    container.close()
    profiler.disable()
    metadata.update(profiled_run_s=perf_counter() - run_started, final_nodes=rows[-1]["nodes"], final_plates=rows[-1]["plates"], mp4_bytes=len(mp4_buf.getvalue()))
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (OUT / "animation.mp4").write_bytes(mp4_buf.getvalue())
    profiler.dump_stats(OUT / "profile.pstats")
    with (OUT / "profile-cumulative.txt").open("w") as f:
        pstats.Stats(profiler, stream=f).strip_dirs().sort_stats("cumulative").print_stats(120)
    with (OUT / "profile-internal.txt").open("w") as f:
        pstats.Stats(profiler, stream=f).strip_dirs().sort_stats("tottime").print_stats(120)
    print("done", json.dumps(metadata), flush=True)


if __name__ == "__main__":
    main()
