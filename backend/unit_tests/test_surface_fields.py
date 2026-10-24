"""The surface-field registry's shared remap rules (surface_fields.py)."""

import numpy as np

from app.surface_fields import BELT_UNMARKED_YEARS, SURFACE_FIELDS, RemapClass, majority_history

NEVER = BELT_UNMARKED_YEARS


def test_belt_marks_remap_by_majority_and_default_to_never():
    spec = SURFACE_FIELDS["belt_formed_years"]
    assert spec.remap_class is RemapClass.MAJORITY_HISTORY
    assert spec.default == spec.sentinel == NEVER


def test_a_cell_is_marked_only_if_most_of_its_area_was_and_keeps_the_latest_date():
    # Cell 0: marked area 3 of 4. Cell 1: exactly half is not most. Cell 2: nothing maps to it.
    cell = np.array([0, 0, 0, 1, 1])
    values = np.array([1e8, 3e8, NEVER, 2e8, NEVER])
    weights = np.array([1.0, 2.0, 1.0, 1.0, 1.0])

    out = majority_history(cell, values, weights, 3, NEVER)

    np.testing.assert_array_equal(out, [3e8, NEVER, NEVER])


def test_an_existing_cells_own_mark_votes_with_its_area():
    cell = np.array([0, 1])
    values = np.array([NEVER, 5e8])
    weights = np.array([1.0, 1.0])
    own = np.array([1e8, 1e8])

    # Cell 0: its own mark outweighs the unmarked crust it takes in; cell 1 keeps the later
    # date of the two.
    out = majority_history(cell, values, weights, 2, NEVER, own, np.array([2.0, 1.0]))

    np.testing.assert_array_equal(out, [1e8, 5e8])
