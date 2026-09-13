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

# Quick Start

A first session. Every cell below runs when this page is built, so what you
read is what the installed package actually returned.

```{code-cell} python
import itacart
```

## Address a position

Give a longitude, a latitude and a resolution. Resolution 13 is one centimetre:

```{code-cell} python
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
cell
```

The nesting is the hierarchy: each parenthesis is one step down the resolution
ladder, and the whole string is the cell's address at every level at once.

## Go back to geography

```{code-cell} python
lon, lat = itacart.cell_to_centroid(cell)
lon, lat
```

The boundary comes back as a ring of coordinates:

```{code-cell} python
ring = itacart.cell_to_boundary(cell, close=True)
len(ring), ring[0]
```

## Move up and down the hierarchy

The parent is one resolution coarser:

```{code-cell} python
parent = itacart.get_parent(cell)
parent
```

Descending gives every child. Pass `flatten=True` to get the children
themselves rather than a grouping:

```{code-cell} python
children = list(itacart.get_children(parent, flatten=True))
len(children), children[0]
```

Twenty-five, because resolution 13 is an odd level. Even levels refine into
four. The count is obtained by enumeration rather than by multiplication, which
matters near the domain border, where a cell can hold fewer children than the
ratio predicts —
[Parent, Children & Descendants](../concepts/index_and_hierarchy/parent_children_and_descendants.md)
covers that case.

## Move sideways

```{code-cell} python
neighbours = itacart.grid_disk(cell, k_distance=1)
len(neighbours)
```

## Measure an area

Area is a count multiplied by a constant, because the cells are equal-area:

```{code-cell} python
itacart.nominal_cell_area(13)
```

For a real parcel, fill it and count what fell inside:

```{code-cell} python
from shapely.geometry import Polygon

parcel = Polygon([
    (-46.6329, -23.5509), (-46.6327, -23.5509),
    (-46.6327, -23.5507), (-46.6329, -23.5507),
])
covered = itacart.count_internal_cells(parcel, 9)
covered, covered * itacart.nominal_cell_area(9)
```

The second number is square metres. No conversion factor appears anywhere,
which is the practical consequence of the areas being decimal.

## Where to go next

[Highlights](highlights.md) states what the grid guarantees.
[Cell Filling / Polyfill](../concepts/spatial_operations/cell_filling_polyfill.md)
covers covering a polygon properly, including the choice of containment mode,
and [Example Notebooks](examples/index.md) works each of these through at
length.
