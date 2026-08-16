import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.plates import MAX_AUTO_PLATES, MIN_AUTO_PLATES


@pytest.fixture
def client():
    return TestClient(app)


def test_render_before_generate_returns_404(client):
    resp = client.get("/world/render")
    assert resp.status_code == 404


def test_step_before_generate_returns_404(client):
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 404


def test_generate_returns_summary(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6})
    assert resp.status_code == 200
    body = resp.json()
    assert body["num_plates"] == 6
    assert body["elapsed_years"] == 0.0
    assert body["seed"] == 1


def test_generate_returns_a_generation_event(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 6, "num_continents": 3})
    body = resp.json()
    assert len(body["events"]) == 1
    assert "6 plates" in body["events"][0]["message"]
    assert "3 continental" in body["events"][0]["message"]
    assert body["events"][0]["elapsed_years"] == 0.0


def test_generate_with_num_continents_gives_exact_count(client):
    resp = client.post("/world/generate", json={"seed": 1, "num_plates": 10, "num_continents": 4})
    render_resp = client.get("/world/render", params={"projection": "behrmann"})
    continental = [p for p in render_resp.json()["plates"] if p["crust_type"] == "continental"]
    assert len(continental) == 4


def test_step_response_includes_growing_event_log(client):
    client.post("/world/generate", json={"seed": 2, "num_plates": 6})
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert len(resp.json()["events"]) >= 1  # at least the generation event


def test_step_advances_elapsed_years(client):
    client.post("/world/generate", json={"seed": 2, "num_plates": 6})
    resp = client.post("/world/step", json={"years": 1_000_000})
    assert resp.status_code == 200
    assert resp.json()["elapsed_years"] == 1_000_000


def test_render_behrmann_and_eckert4(client):
    client.post("/world/generate", json={"seed": 3, "num_plates": 6})
    for projection in ("behrmann", "eckert4"):
        resp = client.get("/world/render", params={"projection": projection})
        assert resp.status_code == 200
        body = resp.json()
        assert body["projection"] == projection
        assert len(body["plates"]) == 6
        first_plate = body["plates"][0]
        assert "crust_type" in first_plate
        assert len(first_plate["lines"]) > 0
        first_line = first_plate["lines"][0]
        assert len(first_line["points"]) == len(first_line["elevation"])
        assert len(first_line["points"][0]) == 2

        for plate in body["plates"]:
            assert "boundary" in plate and isinstance(plate["boundary"], list)
            assert "rotation_rate_deg_per_myr" in plate
            if plate["pole"] is not None:
                assert len(plate["pole"]) == 2
            if plate["velocity_arrow"] is not None:
                assert set(plate["velocity_arrow"].keys()) == {"start", "end"}


def test_render_unknown_projection_returns_400(client):
    client.post("/world/generate", json={"seed": 4, "num_plates": 6})
    resp = client.get("/world/render", params={"projection": "mercator"})
    assert resp.status_code == 400


def test_generate_replaces_previous_world(client):
    client.post("/world/generate", json={"seed": 5, "num_plates": 6})
    resp = client.post("/world/generate", json={"seed": 6, "num_plates": 8})
    assert resp.json()["num_plates"] == 8
    render_resp = client.get("/world/render", params={"projection": "behrmann"})
    assert len(render_resp.json()["plates"]) == 8


def test_generate_without_num_plates_picks_a_plausible_count(client):
    resp = client.post("/world/generate", json={"seed": 7})
    assert resp.status_code == 200
    assert MIN_AUTO_PLATES <= resp.json()["num_plates"] <= MAX_AUTO_PLATES
