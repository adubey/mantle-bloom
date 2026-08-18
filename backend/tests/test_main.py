import base64
import io
import math

import pytest
from fastapi.testclient import TestClient
from PIL import Image

from app.main import app
from app.plates import MAX_AUTO_PLATES, MIN_AUTO_PLATES, generate_plates


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


def test_stats_before_generate_returns_404(client):
    assert client.get("/world/stats").status_code == 404


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


def test_step_response_includes_growing_event_log(client):
    client.post("/world/generate", json={"seed": 2, "num_plates": 6})
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert len(resp.json()["events"]) >= 1  # at least the generation event


def test_step_advances_elapsed_years(client):
    client.post("/world/generate", json={"seed": 2, "num_plates": 6})
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 200
    assert resp.json()["elapsed_years"] == 1_000_000


def test_render_returns_a_decodable_png_at_the_requested_size(client):
    client.post("/world/generate", json={"seed": 3, "num_plates": 6})
    for projection in ("behrmann", "eckert4"):
        for view in ("elevation", "plates", "platesDetail"):
            resp = client.get(
                "/world/render",
                params={"projection": projection, "view": view, "width": 400, "height": 300},
            )
            assert resp.status_code == 200
            body = resp.json()
            assert body["projection"] == projection
            assert "elapsed_years" in body

            image = _decode_image(body)
            assert image.format == "PNG"
            assert image.size == (400, 300)


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


def test_rivers_endpoint_returns_networks_with_expected_shape_after_stepping(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(3):
        client.post("/world/step", json={"years": 2_000_000})

    body = client.get("/world/rivers").json()
    rivers = body["rivers"]
    assert len(rivers) > 0  # seed 20 reliably produces rivers after 3 steps

    seen_ids = set()
    for river in rivers:
        assert river["river_id"] not in seen_ids
        seen_ids.add(river["river_id"])
        assert river["num_nodes"] >= 1
        assert len(river["mouth_xyz"]) == 3
        assert river["mouth_type"] in ("ocean", "lake", "other")
        assert river["flow_rate"] >= 0.0
        assert river["speed"] >= 0.0
        assert river["num_tributaries"] >= 0
        for a, b in river["segments"]:
            assert len(a) == 3 and len(b) == 3


def test_river_at_finds_the_river_nearest_its_own_mouth(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(3):
        client.post("/world/step", json={"years": 2_000_000})

    rivers = client.get("/world/rivers").json()["rivers"]
    target = rivers[0]
    x, y, z = target["mouth_xyz"]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/river_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    assert resp.json()["river_id"] == target["river_id"]


def test_river_at_rejects_non_finite_query(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    assert client.get("/world/river_at", params={"lat_deg": "nan", "lon_deg": 0}).status_code == 400
    assert client.get("/world/river_at", params={"lat_deg": 0, "lon_deg": "inf"}).status_code == 400


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


def test_stats_reflects_world_stepping(client):
    client.post("/world/generate", json={"seed": 14, "num_plates": 8})
    before = client.get("/world/stats").json()
    client.post("/world/step", json={"years": 5_000_000})
    after = client.get("/world/stats").json()
    assert after["elapsed_years"] > before["elapsed_years"]
