"""Sea-ice extent statistics shared by RFC 0015 gate G5's model driver and its ERA5 reference
(#548), so both sides measure the same way. Pure numpy; imported by `gate_sea_ice.py` and
`era5_sea_ice_climatology.py`.

On a regular latitude-longitude grid of cell areas `area` (H, W), with boolean `ocean`,
`perennial` and `seasonal` masks:
- areas by hemisphere and by `BAND_DEG` latitude band poleward of `EDGE_MIN_LAT_DEG`;
- the *equivalent-area edge* (RFC 0015 G5): the latitude poleward of which the ocean area
  equals the ice area. It is reported only where the hemisphere has ocean poleward of
  `EDGE_MIN_LAT_DEG`, and computed for perennial ice and for all ice (perennial + seasonal);
- the per-longitude edge: in each longitude column with ocean poleward of
  `EDGE_MIN_LAT_DEG`, the most equatorward latitude of ice in that column.
"""

from __future__ import annotations

import numpy as np

EARTH_RADIUS_KM = 6371.0
EDGE_MIN_LAT_DEG = 60.0
BAND_DEG = 10
BAND_MIN_LAT_DEG = 40


def cell_areas_km2(lat_deg: np.ndarray, n_lon: int) -> np.ndarray:
    """(H, n_lon) exact areas of cells centred on `lat_deg` (ascending or descending, evenly
    spaced), bounded half-way between centres and clipped at the poles."""
    lat = np.asarray(lat_deg, dtype=float)
    half = 0.5 * abs(lat[1] - lat[0])
    north = np.radians(np.clip(np.maximum(lat + half, lat - half), -90.0, 90.0))
    south = np.radians(np.clip(np.minimum(lat + half, lat - half), -90.0, 90.0))
    band = 2.0 * np.pi * EARTH_RADIUS_KM**2 * (np.sin(north) - np.sin(south))
    return np.repeat((band / n_lon)[:, None], n_lon, axis=1)


def equivalent_area_edge(lat_deg: np.ndarray, area: np.ndarray, ocean: np.ndarray, ice: np.ndarray, north: bool) -> float | None:
    """|Latitude| enclosing, poleward of it, as much ocean as there is ice in the hemisphere;
    linear within a row. None where the hemisphere has no ocean poleward of EDGE_MIN_LAT_DEG."""
    lat = np.asarray(lat_deg, dtype=float)
    side = lat > 0 if north else lat < 0
    polar = side & (np.abs(lat) >= EDGE_MIN_LAT_DEG)
    if not (ocean[polar].any()):
        return None
    ice_area = float(area[ice & ocean & side[:, None]].sum())
    rows = np.flatnonzero(side)
    rows = rows[np.argsort(-np.abs(lat[rows]))]  # pole first
    half = 0.5 * abs(lat[1] - lat[0])
    total = 0.0
    for row in rows:
        row_ocean = float(area[row][ocean[row]].sum())
        if total + row_ocean >= ice_area and row_ocean > 0.0:
            frac = (ice_area - total) / row_ocean
            edge = min(abs(lat[row]) + half, 90.0) - frac * 2.0 * half
            return float(edge)
        total += row_ocean
    return 0.0


def per_longitude_edge(lat_deg: np.ndarray, ocean: np.ndarray, ice: np.ndarray, north: bool) -> np.ndarray:
    """(W,) most equatorward |latitude| of ice in each column: NaN for a column with no ocean
    poleward of EDGE_MIN_LAT_DEG, 90 for one that has such ocean but no ice."""
    lat = np.asarray(lat_deg, dtype=float)
    side = (lat > 0) if north else (lat < 0)
    polar = side & (np.abs(lat) >= EDGE_MIN_LAT_DEG)
    has_ocean = ocean[polar].any(axis=0)
    abs_lat = np.where((ice & ocean & side[:, None]), np.abs(lat)[:, None], np.inf)
    edge = abs_lat.min(axis=0)
    edge = np.where(np.isfinite(edge), edge, 90.0)
    return np.where(has_ocean, edge, np.nan)


def _edge_spread(edges: np.ndarray) -> dict | None:
    valid = edges[np.isfinite(edges)]
    if valid.size == 0:
        return None
    return {
        "columns": int(valid.size),
        "min": round(float(valid.min()), 2),
        "p05": round(float(np.percentile(valid, 5)), 2),
        "median": round(float(np.median(valid)), 2),
        "p95": round(float(np.percentile(valid, 95)), 2),
        "max": round(float(valid.max()), 2),
    }


def hemisphere_summary(lat_deg: np.ndarray, area: np.ndarray, ocean: np.ndarray, perennial: np.ndarray, seasonal: np.ndarray) -> dict:
    """Areas (million km2), bands and edges for each hemisphere. `seasonal` excludes perennial."""
    lat = np.asarray(lat_deg, dtype=float)
    any_ice = perennial | seasonal
    out = {}
    for name, north in (("north", True), ("south", False)):
        side = (lat > 0) if north else (lat < 0)
        rows = side[:, None] & ocean

        def mkm2(mask: np.ndarray, rows: np.ndarray = rows) -> float:
            return round(float(area[mask & rows].sum()) / 1e6, 3)

        bands = []
        for lo in range(BAND_MIN_LAT_DEG, 90, BAND_DEG):
            band = side & (np.abs(lat) >= lo) & (np.abs(lat) < lo + BAND_DEG)
            band_ocean = band[:, None] & ocean
            ocean_mkm2 = float(area[band_ocean].sum()) / 1e6
            bands.append({
                "band": f"{lo}-{lo + BAND_DEG}",
                "ocean_mkm2": round(ocean_mkm2, 3),
                "perennial_share": round(float(area[perennial & band_ocean].sum()) / 1e6 / ocean_mkm2, 3) if ocean_mkm2 else None,
                "seasonal_share": round(float(area[seasonal & band_ocean].sum()) / 1e6 / ocean_mkm2, 3) if ocean_mkm2 else None,
            })
        perennial_edge = equivalent_area_edge(lat, area, ocean, perennial, north)
        all_edge = equivalent_area_edge(lat, area, ocean, any_ice, north)
        out[name] = {
            "ocean_mkm2": mkm2(np.ones_like(ocean)),
            "perennial_mkm2": mkm2(perennial),
            "seasonal_mkm2": mkm2(seasonal),
            "ocean_poleward_50_perennial_share": round(
                float(area[perennial & rows & (np.abs(lat) >= 50.0)[:, None]].sum())
                / max(float(area[rows & (np.abs(lat) >= 50.0)[:, None]].sum()), 1e-9), 3),
            "perennial_edge_deg": None if perennial_edge is None else round(perennial_edge, 2),
            "all_ice_edge_deg": None if all_edge is None else round(all_edge, 2),
            "perennial_edge_by_longitude": _edge_spread(per_longitude_edge(lat, ocean, perennial, north)),
            "bands": bands,
        }
    return out
