import base64
import io
import math
import threading
import pytest
from fastapi.testclient import TestClient
from PIL import Image
from app.lithosphere_plate import generate_plates
from app.main import app
from app.plates import MAX_AUTO_PLATES, MIN_AUTO_PLATES


@pytest.fixture
def client():
    return TestClient(app)


def _decode_image(body: dict) -> Image.Image:
    return Image.open(io.BytesIO(base64.b64decode(body["image_base64"])))


def test_render_before_generate_returns_404(client):
    resp = client.get("/world/render")
    assert resp.status_code == 404


def test_step_before_generate_returns_404(client):
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 404


def test_plates_before_generate_returns_404(client):
    assert client.get("/world/plates").status_code == 404


def test_plate_at_before_generate_returns_404(client):
    assert client.get("/world/plate_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_rivers_before_generate_returns_404(client):
    assert client.get("/world/rivers").status_code == 404


def test_river_at_before_generate_returns_404(client):
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_lakes_before_generate_returns_404(client):
    assert client.get("/world/lakes").status_code == 404


def test_lake_at_before_generate_returns_404(client):
    assert client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": 0}).status_code == 404


def test_save_before_generate_returns_404(client):
    assert client.get("/world/save").status_code == 404


def test_animate_before_generate_returns_404(client):
    resp = client.post("/world/animate", json={"years_per_frame": 1_000_000, "num_frames": 2})
    assert resp.status_code == 404


def test_export_hexgrid_before_generate_returns_404(client):
    assert client.post("/world/export_hexgrid", json={}).status_code == 404


def test_stats_before_generate_returns_404(client):
    assert client.get("/world/stats").status_code == 404


def test_overlapping_step_returns_503(client, monkeypatch):
    client.post("/world/generate", json={"seed": 1, "num_plates": 6})

    # Makes the first /world/step's own critical section deterministically overlap the
    # second's: the blocking replacement only runs (and lets app.main._step_lock be released)
    # once entered_step_world confirms the first request is already holding the lock, so the
    # second request's 503 isn't a race on request scheduling order.
    entered_step_world = threading.Event()
    release_step_world = threading.Event()

    def blocking_step_world(world, years):
        entered_step_world.set()
        release_step_world.wait(timeout=5)

    import app.main as main_module

    monkeypatch.setattr(main_module, "step_world", blocking_step_world)

    results: list[int] = []
    t1 = threading.Thread(target=lambda: results.append(client.post("/world/step", json={"years": 1_000_000}).status_code))
    t1.start()
    assert entered_step_world.wait(timeout=5)

    t2_response = client.post("/world/step", json={"years": 1_000_000})
    release_step_world.set()
    t1.join()

    assert t2_response.status_code == 503
    assert results == [200]


def test_render_waits_for_an_in_progress_step_instead_of_racing_it(client, monkeypatch):
    # Regression test: GET /world/render used to run with no lock at all, so it could read
    # plates.collect_all_points/collect_all_lake_depth/etc. (each independently walking
    # world.plates) while an in-progress step's own deform() was mid-mutation, desyncing
    # those separately-gathered arrays and crashing with an IndexError -- confirmed directly
    # by deliberately racing a step against repeated renders before this fix existed. Same
    # blocking_step_world technique as test_overlapping_step_returns_503 above, but this
    # time proving render *waits out* the lock (and then succeeds) rather than either racing
    # past it or being rejected with a 503 the way the write endpoints are -- a render is a
    # read, called far more often, so it should wait a moment, not fail.
    client.post("/world/generate", json={"seed": 1, "num_plates": 6})

    entered_step_world = threading.Event()
    release_step_world = threading.Event()

    def blocking_step_world(world, years):
        entered_step_world.set()
        release_step_world.wait(timeout=5)

    import app.main as main_module

    monkeypatch.setattr(main_module, "step_world", blocking_step_world)

    step_results: list[int] = []
    t1 = threading.Thread(target=lambda: step_results.append(client.post("/world/step", json={"years": 1_000_000}).status_code))
    t1.start()
    assert entered_step_world.wait(timeout=5)

    render_results: list[int] = []
    t2 = threading.Thread(target=lambda: render_results.append(client.get("/world/render").status_code))
    t2.start()
    # The render thread must actually be blocked on _step_lock right now, not racing through
    # -- confirmed by it still not having finished shortly after starting, while the step
    # still holds the lock.
    t2.join(timeout=0.5)
    assert t2.is_alive()
    assert render_results == []

    release_step_world.set()
    t1.join(timeout=5)
    t2.join(timeout=5)

    assert step_results == [200]
    assert render_results == [200]


def test_generate_returns_summary(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_plates"] == 6
    assert body["elapsed_years"] == 0.0
    assert body["seed"] == 1


def test_generate_returns_a_generation_event(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "continental_fraction": 0.5})
    body = resp.json()
    assert len(body["events"]) == 1
    assert "6 plates" in body["events"][0]["message"]
    assert "3 continental" in body["events"][0]["message"]
    assert body["events"][0]["elapsed_years"] == 0.0


def test_generate_with_continental_fraction_gives_exact_count(client):
    # crust_type per plate isn't part of the render response (see render_image.py) -- the
    # generation event log is the documented way to confirm the exact count (also covered by
    # test_generate_returns_a_generation_event above; this locks in a second seed/count pair).
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 10, "continental_fraction": 0.4})
    assert "4 continental" in resp.json()["events"][0]["message"]


def test_generate_with_land_fraction(client):
    resp = client.post(
        "/world/generate",
        json={"seed": 1, "num_plates": 10, "continental_fraction": 0.7, "land_fraction": 0.29},
    )
    assert resp.status_code == 200


def test_generate_with_node_density(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "node_density": 4.0})
    assert resp.status_code == 200


def test_generate_with_the_coarsest_node_density(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "node_density": 0.5})
    assert resp.status_code == 200


def test_generate_with_unknown_node_density_returns_400(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "node_density": 2.5})
    assert resp.status_code == 400


def test_generate_with_climate_density(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "climate_density": 2.0})
    assert resp.status_code == 200


def test_generate_with_the_finest_climate_density(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "climate_density": 4.0})
    assert resp.status_code == 200


def test_generate_with_unknown_climate_density_returns_400(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "climate_density": 3.0})
    assert resp.status_code == 400


def test_render_defaults_to_elevation_view_at_1100x611(client):
    client.post("/world/generate", json={"seed": 3, "num_plates": 6})
    resp = client.get("/world/render")
    assert resp.status_code == 200
    assert _decode_image(resp.json()).size == (1100, 611)


def test_render_different_views_produce_different_images(client):
    """Not a pixel-exact check (that would just re-derive render_image.py's own math) --
    just confirms the `view` param actually changes what's drawn, at the HTTP layer."""
    client.post("/world/generate", json={"seed": 3, "num_plates": 8, "continental_fraction": 0.5})
    bodies = {
        view: client.get("/world/render", params={"view": view, "width": 300, "height": 200}).json()["image_base64"]
        for view in ("elevation", "plates", "platesDetail")
    }
    assert len(set(bodies.values())) == 3


def test_render_different_resolutions_produce_different_size_images(client):
    client.post("/world/generate", json={"seed": 3, "num_plates": 6})
    small = client.get("/world/render", params={"width": 200, "height": 100}).json()
    large = client.get("/world/render", params={"width": 800, "height": 400}).json()
    assert _decode_image(small).size == (200, 100)
    assert _decode_image(large).size == (800, 400)


def test_render_unknown_projection_returns_400(client):
    client.post("/world/generate", json={"seed": 4, "num_plates": 6})
    resp = client.get("/world/render", params={"projection": "mercator"})
    assert resp.status_code == 400


def test_render_unknown_view_returns_400(client):
    client.post("/world/generate", json={"seed": 4, "num_plates": 6})
    resp = client.get("/world/render", params={"view": "topographic"})
    assert resp.status_code == 400


def test_render_out_of_range_dimensions_return_400(client):
    client.post("/world/generate", json={"seed": 4, "num_plates": 6})
    assert client.get("/world/render", params={"width": 0, "height": 100}).status_code == 400
    assert client.get("/world/render", params={"width": 100, "height": 0}).status_code == 400
    assert client.get("/world/render", params={"width": 100_000, "height": 100}).status_code == 400


def test_generate_replaces_previous_world(client):
    client.post("/world/generate", json={"seed": 5, "num_plates": 6})
    resp = client.post("/world/generate", json={"seed": 6, "num_plates": 8})
    assert resp.json()["num_plates"] == 8


def test_generate_without_num_plates_picks_a_plausible_count(client):
    resp = client.post("/world/generate", json={"seed": 7})
    assert resp.status_code == 200
    assert MIN_AUTO_PLATES <= resp.json()["num_plates"] <= MAX_AUTO_PLATES


def test_render_omitted_rotation_matches_explicit_identity(client):
    client.post("/world/generate", json={"seed": 8, "num_plates": 8})
    identity = "1,0,0,0,1,0,0,0,1"
    default_body = client.get("/world/render", params={"width": 300, "height": 200}).json()
    explicit_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": identity}).json()
    assert default_body["image_base64"] == explicit_body["image_base64"]


def test_render_nontrivial_rotation_changes_the_image(client):
    client.post("/world/generate", json={"seed": 8, "num_plates": 8})
    identity = "1,0,0,0,1,0,0,0,1"
    rotated = "0,0,1,0,1,0,-1,0,0"  # a valid 90-degree rotation matrix
    identity_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": identity}).json()
    rotated_body = client.get("/world/render", params={"width": 300, "height": 200, "rotation": rotated}).json()
    assert identity_body["image_base64"] != rotated_body["image_base64"]


def test_render_malformed_rotation_returns_400(client):
    client.post("/world/generate", json={"seed": 8, "num_plates": 8})
    assert client.get("/world/render", params={"rotation": "1,0,0,0,1,0,0,0"}).status_code == 400  # only 8 values
    assert client.get("/world/render", params={"rotation": "not,a,valid,rotation,matrix,at,all,here,either"}).status_code == 400
    assert client.get("/world/render", params={"rotation": "1,0,0,0,1,0,0,0,nan"}).status_code == 400


def test_plates_endpoint_matches_directly_generated_plate_state(client):
    resp = client.post("/world/generate", json={"seed": 11, "num_plates": 9, "continental_fraction": 0.5})
    assert resp.status_code == 200

    ground_truth = {p.plate_id: p for p in generate_plates(seed=11, num_plates=9, continental_fraction=0.5)}
    body = client.get("/world/plates").json()
    assert len(body["plates"]) == len(ground_truth)

    for entry in body["plates"]:
        truth = ground_truth[entry["plate_id"]]
        assert entry["crust_type"] == truth.crust_type
        assert entry["num_points"] == truth.node_count()
        assert entry["num_rows"] == sum(1 for line in truth.lines if len(line.theta) > 0)
        assert len(entry["outline"]) == len(truth.outline_world())
        assert len(entry["points"]) == truth.node_count()
        if truth.node_count() > 0:
            assert entry["bounding_ellipse"] is not None
            assert entry["bounding_ellipse"]["diameter_a_km"] >= entry["bounding_ellipse"]["diameter_b_km"] >= 0.0
        else:
            assert entry["bounding_ellipse"] is None


def test_plate_at_returns_the_owning_plate_id(client):
    client.post("/world/generate", json={"seed": 12, "num_plates": 8})
    plates = client.get("/world/plates").json()["plates"]
    target = next(p for p in plates if p["num_points"] > 0)
    x, y, z = target["outline"][0]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/plate_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    assert resp.json()["plate_id"] == target["plate_id"]


def test_plate_at_rejects_non_finite_query(client):
    client.post("/world/generate", json={"seed": 12, "num_plates": 8})
    assert client.get("/world/plate_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/plate_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_rivers_and_river_at_are_empty_before_the_first_step(client):
    # hydrology_cache is None until erosion.py runs once (see World.hydrology_cache) --
    # /world/rivers should degrade to an empty list rather than erroring, same spirit as
    # /world/plates always having *something* right after generate.
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/rivers")
    assert resp.status_code == 200
    assert resp.json()["rivers"] == []
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": 0}).json()["river_id"] is None


def test_river_at_rejects_non_finite_query(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/river_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_lakes_and_lake_at_are_empty_before_the_first_step(client):
    # hydrology_cache is None until erosion.py runs once, same degrade-to-empty contract
    # /world/rivers already has -- see test_rivers_and_river_at_are_empty_before_the_first_step.
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    resp = client.get("/world/lakes")
    assert resp.status_code == 200
    assert resp.json()["lakes"] == []
    at = client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": 0}).json()
    assert at == {"kind": "no_basin", "basin": None}


def test_lake_at_rejects_non_finite_query(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/lake_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/lake_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


def test_stats_returns_expected_shape(client):
    client.post("/world/generate", json={"seed": 13, "num_plates": 8, "continental_fraction": 0.5})
    resp = client.get("/world/stats")
    assert resp.status_code == 200
    body = resp.json()

    assert body["elapsed_years"] == 0.0
    assert 0.0 <= body["land_fraction"] <= 1.0
    assert 0.0 <= body["ocean_fraction"] <= 1.0
    assert math.isclose(body["land_fraction"] + body["ocean_fraction"], 1.0, abs_tol=1e-9)
    assert body["elevation_min_m"] <= body["elevation_mean_m"] <= body["elevation_max_m"]
    # A seed/continental_fraction generating both crust types should have both land and
    # ocean grid cells, so none of these should fall back to their None (empty-mask) case.
    for key in (
        "land_temperature_min_c", "land_temperature_max_c", "land_temperature_mean_c",
        "air_temperature_min_c", "air_temperature_max_c", "air_temperature_mean_c",
        "ocean_temperature_min_c", "ocean_temperature_max_c", "ocean_temperature_mean_c",
    ):
        assert body[key] is not None
    assert body["precipitation_min_mm"] <= body["precipitation_mean_mm"] <= body["precipitation_max_mm"]


def test_save_then_load_round_trips_the_exact_world_state(client):
    client.post("/world/generate", json={"seed": 5, "num_plates": 8, "continental_fraction": 0.5})
    client.post("/world/step", json={"years": 2_000_000})
    generated_summary = client.get("/world/plates").json()

    saved = client.get("/world/save")
    assert saved.status_code == 200
    assert saved.headers["content-type"] == "application/octet-stream"
    assert "attachment" in saved.headers["content-disposition"]

    # Advance further so the live world visibly differs from the snapshot just saved.
    client.post("/world/step", json={"years": 2_000_000})
    assert client.get("/world/render").json()["elapsed_years"] == 4_000_000.0

    loaded = client.post("/world/load", content=saved.content)
    assert loaded.status_code == 200
    body = loaded.json()
    assert body["seed"] == 5
    assert body["elapsed_years"] == 2_000_000.0
    assert client.get("/world/plates").json() == generated_summary


def test_load_with_malformed_bytes_returns_400(client):
    resp = client.post("/world/load", content=b"definitely not a saved world")
    assert resp.status_code == 400


def test_animate_advances_the_world_and_returns_a_multi_frame_gif(client):
    client.post("/world/generate", json={"seed": 9, "num_plates": 6})
    resp = client.post(
        "/world/animate",
        json={"projection": "eckert4", "view": "elevation", "width": 200, "height": 110, "years_per_frame": 1_000_000, "num_frames": 3},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["elapsed_years"] == 2_000_000.0  # (num_frames - 1) * years_per_frame

    image = Image.open(io.BytesIO(base64.b64decode(body["image_base64"])))
    assert image.format == "GIF"
    frame_count = 0
    try:
        while True:
            image.seek(frame_count)
            frame_count += 1
    except EOFError:
        pass
    assert frame_count == 3


def test_animate_rejects_out_of_range_frame_counts(client):
    client.post("/world/generate", json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate", json={"years_per_frame": 1_000_000, "num_frames": 0})
    assert resp.status_code == 400
    resp = client.post("/world/animate", json={"years_per_frame": 1_000_000, "num_frames": 10_000})
    assert resp.status_code == 400


def test_animate_rejects_unknown_view(client):
    client.post("/world/generate", json={"seed": 9, "num_plates": 6})
    resp = client.post("/world/animate", json={"view": "not-a-real-view", "years_per_frame": 1_000_000, "num_frames": 2})
    assert resp.status_code == 400


def test_export_hexgrid_returns_the_requested_tile_count(client):
    client.post("/world/generate", json={"seed": 11, "num_plates": 8})
    resp = client.post("/world/export_hexgrid", json={"frequency": 8})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_tiles"] == 642 == len(body["tiles"])
    assert sum(1 for t in body["tiles"] if t["is_pentagon"]) == 12


def test_export_hexgrid_rejects_unknown_frequency(client):
    client.post("/world/generate", json={"seed": 11, "num_plates": 8})
    resp = client.post("/world/export_hexgrid", json={"frequency": 5})
    assert resp.status_code == 400


def test_generate_with_unknown_fluid_density_returns_400(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "fluid_density": 3.0})
    assert resp.status_code == 400


def test_generate_with_the_finest_allowed_fluid_density(client):
    # fluid_density is capped at 2.0 ("High") -- unlike climate_density, which allows 4.0 --
    # since Ocean/Atmospheric Fluid Dynamics now runs every step, not just opt-in.
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "fluid_density": 2.0})
    assert resp.status_code == 200


def test_retired_cfd_sediment_views_are_not_available(client):
    # The ocean CFD (and its sediment concentration/deposition views) was retired -- see
    # render_image.py's VIEWS, which no longer lists them at all.
    client.post("/world/generate", json={"seed": 12, "num_plates": 6, "fluid_density": 0.5})
    for view in ["oceanCfdSediment", "oceanCfdDeposition"]:
        assert client.get("/world/render", params={"view": view}).status_code == 400


def test_climate_views_render_natively_off_the_healpix_grid(client):
    # Every CLIMATE_VIEWS member renders directly off the world's own HEALPix grid (see
    # render_image._render_climate_view) rather than resampling down to an equirectangular
    # grid first.
    from app import render_image

    client.post("/world/generate", json={"seed": 21, "num_plates": 6, "node_density": 0.5, "climate_density": 0.5, "fluid_density": 0.5})
    client.post("/world/step", json={"years": 2_000_000})
    for view in render_image.CLIMATE_VIEWS:
        resp = client.get("/world/render", params={"view": view, "width": 100, "height": 50})
        assert resp.status_code == 200
        assert len(resp.json()["image_base64"]) > 0


def test_wind_and_ocean_currents_views_are_also_always_renderable(client):
    # `wind` draws the CFD-sourced wind; `oceanCurrents` draws climate.py's diagnostic
    # currents (see climate.py's own module docstring) -- both always renderable regardless of
    # what else the world is doing.
    client.post("/world/generate", json={"seed": 12, "num_plates": 6, "fluid_density": 0.5})
    assert client.get("/world/render", params={"view": "wind"}).status_code == 200
    assert client.get("/world/render", params={"view": "oceanCurrents"}).status_code == 200


def test_step_advances_atmosphere_cfd_by_its_own_fixed_seconds(client):
    # Regardless of the tectonic `years` requested, the atmospheric wind solve always advances
    # by its own fixed real-time increment per step (atmosphere_cfd.SECONDS_PER_TECTONIC_STEP)
    # -- reaches into main.py's internal `_state` since this is not surfaced through any HTTP
    # response.
    from app import atmosphere_cfd, main

    client.post("/world/generate", json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    world = main._state["world"]
    assert world.atmosphere_cfd_state.elapsed_seconds == 0.0

    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 200
    assert world.atmosphere_cfd_state.elapsed_seconds == pytest.approx(atmosphere_cfd.SECONDS_PER_TECTONIC_STEP)


def test_step_does_not_advance_fluid_dynamics_when_climate_biomes_paused(client):
    from app import main

    client.post("/world/generate", json={"seed": 12, "num_plates": 6, "climate_density": 0.5, "fluid_density": 0.5})
    client.post("/world/controls", json={"simulate_climate_biomes": False})
    world = main._state["world"]

    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 200
    assert world.atmosphere_cfd_state.elapsed_seconds == 0.0
