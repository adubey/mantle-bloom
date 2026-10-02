# Release notes

## Sparse-quad world cutover

New worlds now use sparse adaptive quad surfaces by default. Generate World's Advanced
Settings keeps a **Legacy elevation lines (for comparison)** option for short-term rollback
and parity work; it is not a second supported production format. The choice resets to
**Standard** after generation, while a legacy-lines badge identifies the current world.

Existing line-backed `.mbworld` files remain loadable under the current compatibility policy.
New saves declare their authoritative surface in the versioned envelope, and a save whose
declaration conflicts with its actual plates is rejected. Line worlds remain line-backed when
loaded normally; the existing explicit `convert_lines` load option can migrate one to quads.
