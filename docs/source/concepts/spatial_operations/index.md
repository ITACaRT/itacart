# Spatial Operations

Operations that take you from geography to cells, and move you around the grid.

Everything in this section reduces to two questions: which cell is a position
in, and which cells are near or inside something. The first has a short answer
and a long one about edges. The second splits into moving on the lattice — by
neighbour, ring and disk — and covering a geometry, which is where the choice
of what "inside" means stops being obvious.

```{toctree}
:maxdepth: 1

point_to_cell_indexing
neighbors_and_grid_traversal
grid_distance
cell_filling_polyfill
containment_modes
```
