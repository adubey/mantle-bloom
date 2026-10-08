import numpy as np
from app.elevation_lines import (
    CRUST_TYPE_CONTINENTAL,
    CRUST_TYPE_INHERIT,
    CRUST_TYPE_OCEANIC,
    TARGET_LINE_SPACING_RAD,
    effective_is_continental_from_codes,
    iter_local_lattice,
    line_spacing_rad,
)


def test_iter_local_lattice_default_spacing_matches_target():
    rows = list(iter_local_lattice(np.eye(3)))
    assert len(rows) > 0
    phis = [phi for phi, _, _ in rows]
    assert np.allclose(np.diff(sorted(phis)), TARGET_LINE_SPACING_RAD, atol=1e-9)


def test_iter_local_lattice_custom_spacing_changes_row_count():
    coarse = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD * 4))
    fine = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD / 2))
    assert len(fine) > len(coarse)
    total_fine = sum(len(theta) for _, theta, _ in fine)
    total_coarse = sum(len(theta) for _, theta, _ in coarse)
    assert total_fine > total_coarse


def test_line_spacing_rad_matches_default_at_density_one():
    assert line_spacing_rad(1.0) == TARGET_LINE_SPACING_RAD


def test_line_spacing_rad_halves_at_4x_density():
    # Node count scales with the square of resolution, so 4x the nodes needs the spacing
    # *halved*, not quartered.
    assert np.isclose(line_spacing_rad(4.0), TARGET_LINE_SPACING_RAD / 2)


def test_effective_is_continental_resolves_inherit_against_the_plate():
    codes = np.array([CRUST_TYPE_INHERIT, CRUST_TYPE_OCEANIC, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_INHERIT], dtype=np.int8)
    # On a continental plate: INHERIT reads continental, explicit codes read at face value.
    assert list(effective_is_continental_from_codes(codes, True)) == [True, False, True, True]
    # On an oceanic plate: only INHERIT flips; the explicit codes are unchanged.
    assert list(effective_is_continental_from_codes(codes, False)) == [False, False, True, False]
