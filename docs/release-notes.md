# Release notes

## Elevation lines retired

The legacy elevation-line terrain model is gone (#251). Generate World's Advanced Settings no
longer offers **Legacy elevation lines**, and every world uses sparse quad surfaces.

Line-backed `.mbworld` files still load: they are converted to quads on load, one way (see
[save-compatibility.md](save-compatibility.md)), and save as quad worlds afterwards. The
`convert_lines` load option is gone, since conversion is now automatic.

The Points debug view's arrow-key stepping along a line, the plate summary's `rows` count and
the corner-notch log panel went with the line model. Points still reports the clicked node's
plate, plate-local coordinates and elevation, and marks it on the map. The Controls window's
debug checkbox now reads **Record debug diagnostics**; it gates the phase budget.

## Sparse-quad world cutover

New worlds now use sparse adaptive quad surfaces by default. Generate World's Advanced
Settings keeps a **Legacy elevation lines (for comparison)** option for short-term rollback
and parity work; it is not a second supported production format. The choice resets to
**Standard** after generation, while a legacy-lines badge identifies the current world.

Existing line-backed `.mbworld` files remain loadable under the current compatibility policy.
New saves declare their authoritative surface in the versioned envelope, and a save whose
declaration conflicts with its actual plates is rejected. Line worlds remain line-backed when
loaded normally; the existing explicit `convert_lines` load option can migrate one to quads.
