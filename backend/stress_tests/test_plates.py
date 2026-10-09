from app.lithosphere_plate import generate_plates
from app.elevation_lines import (
    NODE_DENSITY_CHOICES,
    TARGET_LINE_SPACING_RAD,
    iter_local_lattice,
    line_spacing_rad,
)
from app.plates import (
    MAX_AUTO_PLATES,
    MIN_AUTO_PLATES,
)


def _measured_land_fraction(plates_list) -> float:
    total = sum(float(p.node_areas_m2().sum()) for p in plates_list)
    land = sum(float(p.node_areas_m2()[p.collect("elevation") > 0].sum()) for p in plates_list)
    return land / total if total else 0.0


def _measured_continental_area_fraction(plates_list) -> float:
    total = sum(float(p.node_areas_m2().sum()) for p in plates_list)
    continental = sum(float(p.node_areas_m2().sum()) for p in plates_list if p.crust_type == "continental")
    return continental / total if total else 0.0


def test_generate_plates_without_num_plates_picks_a_plausible_count():
    for seed in range(20):
        plates = generate_plates(seed=seed)
        assert MIN_AUTO_PLATES <= len(plates) <= MAX_AUTO_PLATES
