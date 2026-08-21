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


def test_lakes_endpoint_returns_lakes_with_expected_shape_after_stepping(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(6):
        client.post("/world/step", json={"years": 2_000_000})

    body = client.get("/world/lakes").json()
    lakes = body["lakes"]
    assert len(lakes) > 0  # seed 20 reliably produces lakes after 6 steps

    seen_ids = set()
    for lake in lakes:
        assert lake["lake_id"] not in seen_ids
        seen_ids.add(lake["lake_id"])
        assert lake["is_lake"] is True
        assert lake["member_count"] >= 1
        assert len(lake["member_xyz"]) == lake["member_count"]
        assert len(lake["floor_xyz"]) == 3
        assert lake["water_elevation_m"] is not None
        assert lake["water_elevation_m"] >= lake["floor_elevation_m"]
        if lake["outlet_elevation_m"] is not None:
            assert len(lake["outlet_xyz"]) == 3
            assert lake["outlet_elevation_m"] >= lake["floor_elevation_m"]
        else:
            assert lake["outlet_xyz"] is None
        for river in lake["inflow_rivers"]:
            assert len(river["mouth_xyz"]) == 3
            assert river["flow_rate"] >= 0.0


def test_lake_at_finds_the_lake_at_its_own_floor(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(6):
        client.post("/world/step", json={"years": 2_000_000})

    lakes = client.get("/world/lakes").json()["lakes"]
    target = lakes[0]
    x, y, z = target["floor_xyz"]
    lat_deg = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    lon_deg = math.degrees(math.atan2(y, x))

    resp = client.get("/world/lake_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg})
    assert resp.status_code == 200
    body = resp.json()
    assert body["kind"] == "lake"
    assert body["basin"]["lake_id"] == target["lake_id"]


def test_lake_at_reports_a_dry_basin_nonetheless(client):
    client.post("/world/generate", json={"seed": 20, "num_plates": 10, "continental_fraction": 0.5})
    for _ in range(6):
        client.post("/world/step", json={"years": 2_000_000})

    found_basin = False
    for lat_deg in range(-80, 81, 10):
        for lon_deg in range(-170, 171, 10):
            body = client.get("/world/lake_at", params={"lat_deg": lat_deg, "lon_deg": lon_deg}).json()
            assert body["kind"] in ("lake", "basin", "no_basin", "ocean")
            if body["kind"] == "basin":
                found_basin = True
                basin = body["basin"]
                assert basin["is_lake"] is False
                assert basin["lake_id"] is None
                assert len(basin["floor_xyz"]) == 3
    assert found_basin  # seed 20 reliably has some never-flooded catchment after 6 steps


def test_stats_reflects_world_stepping(client):
    client.post("/world/generate", json={"seed": 14, "num_plates": 8})
    before = client.get("/world/stats").json()
    client.post("/world/step", json={"years": 5_000_000})
    after = client.get("/world/stats").json()
    assert after["elapsed_years"] > before["elapsed_years"]
