---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Traverse the Grid

Move around by neighbour, ring and disk — and see what happens when a step
reaches a seam.

```{code-cell} python
import itacart
from itacart import topology

cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
cell
```

## One step

```{code-cell} python
north = topology.get_neighbor(cell, "N")
north, topology.are_neighbor_cells(cell, north)
```

Directions are lattice compass codes — `N`, `S`, `E`, `W` and the four
diagonals. Every direction at once:

```{code-cell} python
{d: topology.get_neighbor(cell, d) is not None
 for d in ("N", "S", "E", "W", "NE", "NW", "SE", "SW")}
```

None of them is `None` here, because this cell is deep in the interior.

## Rings and disks

```{code-cell} python
[(k, len(topology.grid_disk(cell, k)), len(topology.grid_ring(cell, k)))
 for k in (1, 2, 3)]
```

A disk is everything within `k` steps including the origin; a ring is the
hollow shell at exactly `k`. The counts follow from the metric:

```{code-cell} python
(len(topology.grid_disk(cell, 1)),
 len(topology.grid_disk(cell, 1, metric="manhattan")))
```

`chebyshev`, the default, counts a diagonal as one step and gives nine.
`manhattan` allows only axis steps and gives five.

Do not compute a disk's size from `k` — near a seam some of the shell names no
cell, and the count comes back smaller.

## Distance

```{code-cell} python
far = itacart.geo_to_cell(-46.6300, -23.5500, resolution=9)
(topology.grid_distance(cell, far),
 topology.grid_distance(cell, far, metric="manhattan"))
```

Steps, not metres. Multiplying by a cell edge gives a figure that looks like
metres and is not one, because the path is a staircase on the lattice rather
than an arc on the Earth.
[Grid Distance](../../concepts/spatial_operations/grid_distance.md).

## Edges as objects

```{code-cell} python
edge = topology.cells_to_directed_edge(cell, north)
edge
```

```{code-cell} python
topology.directed_edge_to_cells(edge) == (cell, north), len(topology.cell_to_edges(cell))
```

## At the prime meridian

Now take a resolution-1 cell in the meridian column, which is a triangle:

```{code-cell} python
from itacart import boundary
tri = "NE(0000/0500)"
boundary.cell_shape(tri)
```

```{code-cell} python
{d: topology.get_neighbor(tri, d) for d in ("W", "E", "N", "S")}
```

Two things to notice. Stepping **west** leaves the quadrant entirely and lands
in `NW` — the quadrants mirror one another, so a step that would go negative
deflects instead, which is why ITACaRT has no negative indices. And stepping
**north** returns `None`: there is no cell there.

```{code-cell} python
topology.deflect(tri, "W")
```

`deflect` is the explicit form of the same resolution, and it is defined on the
resolution-1 lattice, where the step is arithmetic on the `XXXX/YYYY` pair.

## At the equator

```{code-cell} python
{d: topology.get_neighbor("NE(0500/0000)", d) for d in ("N", "S")}
```

Stepping south from row 0 crosses into the southern quadrant. Stepping north
moves one row up **and one column west** — that is the shear of the lattice
showing through, not an error.
[Meridian and Equator](../../concepts/boundaries/meridian_and_equator.md).

## A walk that crosses a seam

```{code-cell} python
here, path = "NE(0003/0500)", ["NE(0003/0500)"]
for _ in range(4):
    here = topology.get_neighbor(here, "W")
    path.append(here)
path
```

Four steps west from column 3, and the walk crosses the prime meridian into the
northwestern quadrant without any special handling on the caller's side.
