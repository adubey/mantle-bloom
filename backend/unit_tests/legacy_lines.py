"""Test-only stand-ins for line-backed plates, the surface retired in #251.

A line save holds `LithospherePlate`s (`app.lithosphere_plate`) whose `_lines` are
`ElevationLine`s (`app.elevation_lines`). This build no longer has either class: it reads a
line save structurally (`legacy_conversion.read_line_plate`) and converts it to quads. These
stand-ins carry exactly the state those classes pickled, so a test can build a line world in
memory and, through `retired_class_names`, write it as a genuine line save.
"""

from contextlib import contextmanager

import numpy as np

from app import geometry, lithosphere
from app.elevation_lines import iter_local_lattice
from app.surface_fields import SURFACE_FIELDS


class LegacyLine:
    """One row of nodes at plate-local latitude `phi`: `_phi`, `_theta` and one `_<field>`
    array per surface field, as `ElevationLine` pickled them."""

    def __init__(self, phi: float, theta, **fields):
        theta = np.asarray(theta, dtype=float)
        self._phi = float(phi)
        self._theta = theta
        for name, spec in SURFACE_FIELDS.items():
            value = fields.get(name)
            if value is None:
                value = np.full(len(theta), spec.default, dtype=spec.dtype)
            setattr(self, f"_{name}", np.asarray(value, dtype=spec.dtype))

    @property
    def phi(self) -> float:
        return self._phi

    @property
    def theta(self) -> np.ndarray:
        return self._theta

    def __len__(self) -> int:
        return len(self._theta)


class LegacyLinePlate:
    """A line-backed plate's pickled state (`_plate_id`, `_frame`, `_crust_type`, `_omega`,
    `_age_steps`, `_internal_stress`, `_lines`), with just enough behaviour for a test to
    read and edit it."""

    def __init__(self, plate_id, frame, crust_type, lines, omega=None, age_steps=0, internal_stress=0.0):
        self._plate_id = plate_id
        self._frame = np.asarray(frame, dtype=float)
        self._crust_type = crust_type
        self._omega = np.zeros(3) if omega is None else np.asarray(omega, dtype=float)
        self._age_steps = age_steps
        self._internal_stress = internal_stress
        self._lines = list(lines)

    @property
    def plate_id(self) -> int:
        return self._plate_id

    @property
    def frame(self) -> np.ndarray:
        return self._frame

    @property
    def crust_type(self) -> str:
        return self._crust_type

    @property
    def omega(self) -> np.ndarray:
        return self._omega

    @property
    def age_steps(self) -> int:
        return self._age_steps

    @property
    def lines(self) -> list[LegacyLine]:
        return self._lines

    def accounting_areas_m2(self, spacing_rad: float) -> np.ndarray:
        """The nominal footprint per node, as line plates accounted for area."""
        return np.full(self.node_count(), lithosphere.node_area_m2(spacing_rad))

    def node_count(self) -> int:
        return sum(len(line) for line in self._lines)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        phi = np.concatenate([np.full(len(line), line.phi) for line in self._lines])
        theta = np.concatenate([line.theta for line in self._lines])
        return geometry.to_world(self._frame, geometry.local_xyz(phi, theta)), self.collect("elevation")

    def collect(self, name: str) -> np.ndarray:
        return np.concatenate([getattr(line, f"_{name}") for line in self._lines])

    def set_fields_on_plate(self, **fields) -> None:
        offset = 0
        for line in self._lines:
            n = len(line)
            for name, values in fields.items():
                setattr(line, f"_{name}", np.asarray(values[offset : offset + n], dtype=SURFACE_FIELDS[name].dtype))
            offset += n

    def rotate(self, increment: np.ndarray) -> None:
        self._frame = increment @ self._frame


def legacy_line_plate(plate_id, frame, crust_type, spacing_rad, is_owned, omega=None) -> LegacyLinePlate:
    """A line plate built the way the line engine generated one: every node of the plate-local
    lattice `is_owned(world_pts)` keeps, at its crust type's reference column with elevation
    from isostasy."""
    hc, hm = lithosphere.reference_thickness(crust_type)
    lines = []
    for phi, theta, world_pts in iter_local_lattice(frame, spacing_rad):
        owned = np.asarray(is_owned(world_pts), dtype=bool)
        if not np.any(owned):
            continue
        n = int(owned.sum())
        hc_row, hm_row = np.full(n, hc), np.full(n, hm)
        elevation = lithosphere.isostatic_elevation(hc_row, hm_row, lithosphere.crust_density(crust_type))
        lines.append(LegacyLine(phi, theta[owned], elevation=elevation, crustal_thickness_m=hc_row, mantle_lithosphere_thickness_m=hm_row))
    return LegacyLinePlate(plate_id, frame, crust_type, lines, omega=omega)


def line_world_like(world, spacing_rad: float):
    """A line world tiling the sphere as `world`'s quad plates do: each plate's lattice nodes
    inside its own cells, carrying the nearest cell's fields. So the line plates abut, carry
    real relief, and convert back to roughly the world they came from."""
    from scipy.spatial import cKDTree

    from app.world import World

    plates = []
    for plate in world.plates:
        line_plate = legacy_line_plate(plate.plate_id, plate.frame, plate.crust_type, spacing_rad, plate.contains_batch, omega=plate.omega)
        nearest = cKDTree(plate.all_points_and_elevation()[0]).query(line_plate.all_points_and_elevation()[0])[1]
        line_plate.set_fields_on_plate(**{name: plate.collect(name)[nearest] for name in SURFACE_FIELDS})
        plates.append(line_plate)
    lines = World(seed=world.seed, plates=plates, next_plate_id=world.next_plate_id, node_density=world.node_density, mantle_centers=world.mantle_centers)
    lines.sea_level_m = world.sea_level_m
    return lines


class LegacyRowLookup:
    """A line plate's cached row index (`_row_lookup_cache`), which a line save pickled
    whenever the plate's containment fast path had run."""

    def __init__(self, phi):
        self.phi = np.asarray(phi, dtype=float)


# The names a real line save pickled its classes under.
_RETIRED = (
    (LegacyLinePlate, "app.lithosphere_plate", "LithospherePlate"),
    (LegacyLine, "app.elevation_lines", "ElevationLine"),
    (LegacyRowLookup, "app.plates", "_RowLookup"),
)


@contextmanager
def retired_class_names():
    """Within this block, pickling a stand-in writes it under its retired class's name, as a
    line save from before #251 did."""
    import importlib

    saved = []
    for cls, module_name, name in _RETIRED:
        module = importlib.import_module(module_name)
        saved.append((cls, cls.__module__, cls.__qualname__, module, name))
        cls.__module__, cls.__qualname__ = module_name, name
        setattr(module, name, cls)
    try:
        yield
    finally:
        for cls, old_module, old_qualname, module, name in saved:
            cls.__module__, cls.__qualname__ = old_module, old_qualname
            delattr(module, name)
