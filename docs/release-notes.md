# Release notes

## Sparse-quad world cutover

New worlds now use sparse adaptive quad surfaces by default. The Generate World dialog keeps
an **Elevation lines (legacy / diagnostic)** option for short-term rollback and parity work;
it is not a second supported production format.

Existing line-backed `.mbworld` files remain loadable under the current compatibility policy.
New saves carry explicit `surface.kind` and `surface.version` metadata in their versioned
envelope (`mixed` is reserved for synthetic diagnostic worlds). A save whose metadata
conflicts with its actual plate representation is rejected rather than guessed or converted.
There is no automatic conversion between line and quad worlds: the selected representation
remains authoritative for that world's lifetime.
