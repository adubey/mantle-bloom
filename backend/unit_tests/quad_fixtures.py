"""Small hand-built quad plates for tests: a block of cells around the centre of cube face 0.

With an identity frame, face 0 is centred on world (lat 0, lon 0); column offsets run east
and row offsets north, one cell (about one node spacing) apart.
"""

import numpy as np

from app import lithosphere
from app.elevation_lines import line_spacing_rad
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys


def cells_per_edge(density: float = 1.0) -> int:
    return cells_per_face_edge(line_spacing_rad(density))


def block_keys(columns=range(4), rows=range(1), density: float = 1.0, n: int | None = None) -> np.ndarray:
    """Cell keys for every (row, column) offset from the face centre, in the plate's own node
    order: row by row, then column by column. `columns` and `rows` must be ascending. `n`
    (cells per face edge) overrides the lattice `density` gives."""
    centre = (cells_per_edge(density) if n is None else n) // 2
    jj, ii = np.meshgrid(centre + np.asarray(rows), centre + np.asarray(columns), indexing="ij")
    return pack_cell_keys(np.zeros(ii.size, dtype=np.int64), ii.ravel(), jj.ravel())


def quad_plate(
    plate_id: int,
    crust_type: str = "oceanic",
    columns=range(4),
    rows=range(1),
    density: float = 1.0,
    frame: np.ndarray | None = None,
    keys: np.ndarray | None = None,
    n: int | None = None,
    **fields,
) -> PlateWithSparseQuadPatch:
    """A plate of `len(rows) * len(columns)` cells (see `block_keys`), or of `keys` when given,
    at its crust type's reference column. `fields` override any field per cell (an array in
    node order) or for every cell (a scalar); elevation is synced from Hc/Hm unless given."""
    keys = block_keys(columns, rows, density, n) if keys is None else np.sort(keys)
    count = len(keys)
    hc, hm = lithosphere.reference_thickness(crust_type)
    values = {"crustal_thickness_m": np.full(count, hc), "mantle_lithosphere_thickness_m": np.full(count, hm)}
    for name, value in fields.items():
        values[name] = np.full(count, value) if np.ndim(value) == 0 else np.asarray(value)
    plate = PlateWithSparseQuadPatch(
        plate_id, np.eye(3) if frame is None else frame, crust_type, cells_per_edge(density) if n is None else n, keys,
        fields=values,
    )
    if "elevation" not in fields:
        lithosphere.sync_plate_elevation(plate)
    return plate


def row_plate(plate_id: int, crust_type: str, theta, frame: np.ndarray | None = None, rows=range(1), **fields) -> PlateWithSparseQuadPatch:
    """One cell per entry of `theta` (uniformly spaced plate-local longitudes), on a lattice
    whose cells are exactly that spacing wide -- so a row of cells stands in for a row of
    nodes at the same positions. `rows` stacks more rows, the same spacing apart."""
    theta = np.asarray(theta, dtype=float)
    step = float(theta[1] - theta[0]) if len(theta) > 1 else 0.01
    n = int(round((np.pi / 2) / step))
    columns = int(np.round(theta[0] / step)) + np.arange(len(theta))
    return quad_plate(plate_id, crust_type, columns=columns, rows=rows, n=n, frame=frame, **fields)


def lobed_plate(lobes, plate_id: int = 0, rows: int = 12, per_row: int = 8, crust_type: str = "oceanic", density: float = 1.0, frame=None, **plate_kwargs):
    """A plate whose cells form one connected blob per entry in `lobes`, each a starting
    column offset or an `(offset, cells_per_row)` tuple. Offsets must sit far enough apart to
    read as separate connected components (2.5 spacings) and stay on the face (within about
    +-35 columns at density 1). Each lobe spans `rows` rows, so it contributes
    `rows * cells_per_row` cells."""
    n = cells_per_edge(density)
    i0, j0 = n // 2 - 20, n // 2 - rows // 2
    columns = np.concatenate(
        [i0 + offset + np.arange(count) for offset, count in ((lobe if isinstance(lobe, tuple) else (lobe, per_row)) for lobe in lobes)]
    )
    jj, ii = np.meshgrid(np.arange(j0, j0 + rows), columns, indexing="ij")
    keys = np.sort(pack_cell_keys(np.zeros(ii.size, dtype=np.int64), ii.ravel(), jj.ravel()))
    return PlateWithSparseQuadPatch(plate_id, np.eye(3) if frame is None else frame, crust_type, n, keys, **plate_kwargs)


def globe_plate(crust_type: str, density: float = 0.5, **fields) -> PlateWithSparseQuadPatch:
    """One plate covering the whole sphere. Each `fields` value is a scalar for every cell or a
    function of the cells' world positions."""
    cells = PlateWithSparseQuadPatch.from_lattice(
        0, np.eye(3), crust_type, line_spacing_rad(density), lambda pts: np.ones(len(pts), dtype=bool)
    )
    points = node_points(cells)
    values = {name: value(points) if callable(value) else value for name, value in fields.items()}
    return quad_plate(0, crust_type, density=density, keys=cells.cell_keys, **values)


def node_points(plate) -> np.ndarray:
    """Every node's world position, in node order."""
    return plate.all_points_and_elevation()[0]
