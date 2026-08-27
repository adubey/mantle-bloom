import numpy as np
from app import geometry
from app.lithosphere_plate import generate_plates
from app.plates import (
    ELLIPSE_OUTLINE_POINTS,
    MAX_AUTO_PLATES,
    MIN_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    NODE_DENSITY_CHOICES,
    TARGET_LINE_SPACING_RAD,
    collect_all_points,
    iter_local_lattice,
    line_spacing_rad,
    nearest_plate_id,
    plate_bounding_ellipse,
)
from app.world import generate_world


def _measured_land_fraction(plates_list) -> float:
    total = sum(p.node_count() for p in plates_list)
    land = sum(int(np.sum(line.elevation > 0)) for p in plates_list for line in p.lines)
    return land / total if total else 0.0


def _measured_continental_area_fraction(plates_list) -> float:
    total = sum(p.node_count() for p in plates_list)
    continental = sum(p.node_count() for p in plates_list if p.crust_type == "continental")
    return continental / total if total else 0.0


def test_generate_plates_without_num_plates_picks_a_plausible_count():
    for seed in range(20):
        plates = generate_plates(seed=seed)
        assert MIN_AUTO_PLATES <= len(plates) <= MAX_AUTO_PLATES
