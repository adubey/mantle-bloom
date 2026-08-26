from fastapi.testclient import TestClient

from app.main import app

client = TestClient(app)

_GENERATE_BODY = {"seed": 21, "node_density": 0.5, "climate_density": 0.5, "fluid_density": 0.5, "num_plates": 6}


def test_v2_mounted_alongside_v1():
    paths = {getattr(r, "path", None) for r in app.routes}
    assert "/v2" in paths


def test_v2_generate_and_step():
    r = client.post("/v2/world/generate", json=_GENERATE_BODY)
    assert r.status_code == 200
    body = r.json()
    assert body["num_plates"] == 6
    assert body["elapsed_years"] == 0.0

    r = client.post("/v2/world/step", json={"years": 2_000_000})
    assert r.status_code == 200
    assert r.json()["elapsed_years"] == 2_000_000.0


def test_v2_render_elevation_view():
    client.post("/v2/world/generate", json=_GENERATE_BODY)
    r = client.get("/v2/world/render", params={"view": "elevation", "width": 100, "height": 50})
    assert r.status_code == 200
    assert len(r.json()["image_base64"]) > 0


def test_v2_render_rejects_sediment_views():
    # Sediment has no HEALPix port yet -- the only two FLUID_VIEWS entries left after the
    # rest (oceanCfdVelocity/oceanCfdTemperature/atmosphereCfdVelocity/
    # atmosphereCfdTemperature/atmosphereCfdHumidity) were removed as duplicates of
    # CLIMATE_VIEWS (see test_v2_climate_views_render_natively_below).
    client.post("/v2/world/generate", json=_GENERATE_BODY)
    r = client.get("/v2/world/render", params={"view": "oceanCfdSediment"})
    assert r.status_code == 400


def test_v2_climate_views_render_natively():
    # Once climate.py started sourcing temperature/humidity/precipitation straight from the
    # CFD state, every CLIMATE_VIEWS member renders directly off the world's own HEALPix grid
    # for a v2 world (see render_image._render_climate_view_healpix) instead of resampling
    # down to an equirectangular grid first -- still succeeds via the same endpoint/view names
    # v1 uses, just a different render path server-side.
    from app import render_image

    client.post("/v2/world/generate", json=_GENERATE_BODY)
    client.post("/v2/world/step", json={"years": 2_000_000})
    for view in render_image.CLIMATE_VIEWS:
        r = client.get("/v2/world/render", params={"view": view, "width": 100, "height": 50})
        assert r.status_code == 200
        assert len(r.json()["image_base64"]) > 0


def test_v2_and_v1_worlds_are_independent():
    r1 = client.post("/world/generate", json={"seed": 1})
    assert r1.status_code == 200
    v1_plates = r1.json()["num_plates"]

    r2 = client.post("/v2/world/generate", json=_GENERATE_BODY)
    assert r2.status_code == 200
    assert r2.json()["num_plates"] == 6

    # v1's own world is untouched by the v2 generate call above.
    r1_again = client.get("/world/summary")
    assert r1_again.status_code == 200
    assert r1_again.json()["num_plates"] == v1_plates


def test_v2_plates_stats_rivers_endpoints():
    client.post("/v2/world/generate", json=_GENERATE_BODY)
    assert client.get("/v2/world/plates").status_code == 200
    assert client.get("/v2/world/stats").status_code == 200
    assert client.get("/v2/world/rivers").status_code == 200
    assert client.get("/v2/world/lakes").status_code == 200


def test_v2_returns_404_before_generate():
    # main_v2.py's `_state` is a module-level singleton shared by every TestClient pointing
    # at the same `v2_app` object (this file's own earlier tests have already generated a
    # world into it) -- reset it directly to exercise the "nothing generated yet" 404 path,
    # then restore a world so later tests in this module aren't affected by ordering.
    from app.v2 import main_v2

    previous = main_v2._state["world"]
    main_v2._state["world"] = None
    try:
        r = client.get("/v2/world/summary")
        assert r.status_code == 404
    finally:
        main_v2._state["world"] = previous
