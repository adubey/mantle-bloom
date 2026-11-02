"""date_saved_world.py -- dating the starting rock of an existing save (#747)."""

import copy
import functools
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "bin" / "debug"))

import date_saved_world as script  # noqa: E402

from app import persistence, surface_parity  # noqa: E402
from app.world import generate_world  # noqa: E402

BORN_LATER_YEARS = 3.0e6


def _created(world):
    return np.concatenate([plate.collect("node_created_years") for plate in world.plates])


@functools.cache
def _generated_world():
    return generate_world(747, num_plates=6, node_density=0.5)


def _old_world():
    """A world as a pre-#744 save holds it: undated, but for a few cells created since t=0."""
    world = copy.deepcopy(_generated_world())
    world.elapsed_years = 5.0e6
    world.initial_rock_ages = {}
    for plate in world.plates:
        created = np.full(plate.node_count(), -1.0)
        created[:3] = BORN_LATER_YEARS
        plate.set_fields_on_plate(node_created_years=created)
    return world


def _save(world, path):
    path.write_bytes(persistence.save_world_bytes(world, name="Old world"))
    return path


def test_only_undated_cells_change_and_ages_vary_within_a_plate():
    world = _old_world()
    before = _created(world)
    summary = script.date_world(world)
    after = _created(world)

    assert not np.any(after == -1.0)
    np.testing.assert_array_equal(after[before != -1.0], before[before != -1.0])
    assert np.all(after[before == -1.0] < 0.0)
    assert summary["cells_dated"] == np.count_nonzero(before == -1.0)
    for plate in world.plates:
        record = world.initial_rock_ages[plate.plate_id]
        ages = -plate.collect("node_created_years")[3:]
        assert len(np.unique(ages)) > 1
        assert record["youngest_age_years"] == ages.min()
        assert record["dated_at_years"] == world.elapsed_years
        assert record["primary_crust_type"] == plate.crust_type


def test_dating_changes_no_other_field():
    world, dated = _old_world(), _old_world()
    script.date_world(dated)
    for plate, source in zip(dated.plates, world.plates):
        plate.set_fields_on_plate(node_created_years=source.collect("node_created_years"))
    assert surface_parity.state_hash(dated) == surface_parity.state_hash(world)


def test_a_plate_with_nothing_undated_is_left_alone():
    world = _old_world()
    dated_plate = world.plates[0]
    original = np.full(dated_plate.node_count(), BORN_LATER_YEARS)
    dated_plate.set_fields_on_plate(node_created_years=original)
    script.date_world(world)
    np.testing.assert_array_equal(dated_plate.collect("node_created_years"), original)
    assert dated_plate.plate_id not in world.initial_rock_ages


def test_the_same_save_dates_the_same_way():
    first, second = _old_world(), _old_world()
    script.date_world(first)
    script.date_world(second)
    np.testing.assert_array_equal(_created(first), _created(second))


def test_the_cli_writes_a_loadable_copy_and_keeps_the_source_and_name(tmp_path, capsys):
    source = _save(_old_world(), tmp_path / "old.mbworld")
    original = source.read_bytes()
    assert script.main([str(source)]) == 0

    out = tmp_path / "old-dated.mbworld"
    assert json.loads(capsys.readouterr().out)["written"] == str(out)
    assert source.read_bytes() == original
    data = out.read_bytes()
    assert persistence.load_world_name(data) == "Old world"
    world = persistence.load_world_bytes(data)
    assert not np.any(_created(world) == -1.0)
    assert set(world.initial_rock_ages) == {plate.plate_id for plate in world.plates}


def test_a_dated_save_writes_nothing(tmp_path, capsys):
    dated = tmp_path / "dated.mbworld"
    world = _old_world()
    script.date_world(world)
    _save(world, dated)
    out = tmp_path / "again.mbworld"
    assert script.main([str(dated), "--out", str(out)]) == 0
    assert json.loads(capsys.readouterr().out)["written"] is None
    assert not out.exists()


def test_the_source_can_never_be_the_output(tmp_path):
    source = _save(_old_world(), tmp_path / "old.mbworld")
    with pytest.raises(SystemExit):
        script.main([str(source), "--out", str(source), "--force"])


def test_an_existing_output_needs_force(tmp_path, capsys):
    source = _save(_old_world(), tmp_path / "old.mbworld")
    out = tmp_path / "taken.mbworld"
    out.write_bytes(b"keep me")
    with pytest.raises(SystemExit):
        script.main([str(source), "--out", str(out)])
    assert out.read_bytes() == b"keep me"
    assert script.main([str(source), "--out", str(out), "--force"]) == 0
    assert out.read_bytes() != b"keep me"
