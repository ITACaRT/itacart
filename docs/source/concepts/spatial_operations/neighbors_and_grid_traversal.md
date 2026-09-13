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

# Neighbors & Grid Traversal

Adjacency is computed from the index, not from coordinates. That is a
requirement the standard makes of a discrete global grid, and it is what keeps
traversal cheap.

```{code-cell} python
import itacart
from itacart import topology
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
north = topology.get_neighbor(cell, "N")
north
```

Directions are lattice compass codes — `N`, `S`, `E`, `W` and the four
diagonals. Two cells are neighbours when they share an edge:

```{code-cell} python
topology.are_neighbor_cells(cell, north)
```

## How the arithmetic works

The rule depends on which kind of level you are on, and all three cases are
arithmetic on the string.

At **resolutions 0 and 1** the step is integer arithmetic on the
`XXXX/YYYY` pair: the Y index decreases for the cell above and increases for
the cell below, the X index decreases to the west and increases to the east.

At **even resolutions** the four children sit in a two-by-two block. Vertical
adjacency adds or subtracts two; horizontal adjacency adds or subtracts one.

At **odd resolutions** the twenty-five children sit in a five-by-five block.
The numeric part of the code moves horizontally — `C2` to `C3` — wrapping
between columns 5 and 1; the alphabetic part moves vertically — `C2` to `B2` —
wrapping between rows E and A.

```{figure} ../../../_static/f6/f6_02_refinement_grids.png
:alt: The two refinement blocks neighbours are computed in
:width: 100%

`docs/_static/f6/f6_02_refinement_grids.png`
```

## Crossing into another parent

When the neighbour is not a sibling, the wrap-around has carried you out of the
parent block. The procedure is to ascend, step at the parent's level, and
descend into the corresponding child:

```{figure} ../../../_static/f6/f6_03_ascend_and_descend.png
:alt: Ascending to the parent, stepping, and descending again
:width: 100%

`docs/_static/f6/f6_03_ascend_and_descend.png`
```

That is why the wrap is safe rather than wrong: a code that wraps is the signal
to go up a level, not an answer in itself.

## When the step falls off the lattice

At a quadrant edge or the domain border, the lexical target of a step can name
no cell. `deflect` resolves those:

```{code-cell} python
from itacart.exceptions import ITACaRTError
try:
    print(topology.deflect("NE(0000/0500)", "W"))
except ITACaRTError as exc:
    print(type(exc).__name__, "--", exc)
```

```{figure} ../../../_static/f6/f6_04_boundary_deflections.png
:alt: Deflections at the quadrant boundaries
:width: 100%

`docs/_static/f6/f6_04_boundary_deflections.png`
```

Because the quadrants mirror one another, ITACaRT never uses negative indices —
a step that would go negative deflects into the neighbouring quadrant instead.
[Meridian and Equator](../boundaries/meridian_and_equator.md) works through the
geometry behind those deflections.

## Rings and disks

A disk is every cell within `k` steps, the origin included; a ring is the
hollow shell at exactly `k`.

```{code-cell} python
[(k, len(topology.grid_disk(cell, k)), len(topology.grid_ring(cell, k)))
 for k in (1, 2, 3)]
```

```{figure} ../../../_static/f6/f6_06_disk_cardinality.png
:alt: How many cells a disk holds at each radius
:width: 100%

`docs/_static/f6/f6_06_disk_cardinality.png`
```

The counts follow from the metric. Under the default `chebyshev` diagonal steps
count as one, so the unit disk is a filled rhombus of nine cells. Under
`manhattan` only axis steps are allowed and the unit disk is a diamond of five:

```{code-cell} python
(len(topology.grid_disk(cell, 1)),
 len(topology.grid_disk(cell, 1, metric="manhattan")))
```

Near a boundary the counts fall below these figures, because some of the shell
names no cell. Do not compute a disk's size from `k` — ask for the disk.

## Edges as first-class objects

A step can be named rather than taken:

```{code-cell} python
edge = topology.cells_to_directed_edge(cell, north)
topology.directed_edge_to_cells(edge) == (cell, north)
```

The identifier is the two cells and a direction, and it round-trips. Every edge
leaving a cell is available at once:

```{code-cell} python
len(topology.cell_to_edges(cell))
```
