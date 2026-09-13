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

# Cell Geometry

Everything geometric about a cell is derived from its index. Nothing is stored,
and nothing is looked up.

```{code-cell} python
import itacart
cell = "SE(1400/0374)"
itacart.cell_to_boundary(cell, close=True)
```

Five positions: four distinct vertices and a repeat of the first, which closes
the ring. Longitude first, then latitude, in EPSG:4326.

## From index to vertices

Descending the index arrives at the anchor — the lattice corner the address
inverts to exactly. The other three vertices follow from the construction: base
length $l$ along the parallel, height $l$ along the meridian, and the lean that
carries the upper edge one base-length to the west.
[The Sinusoidal Construction](../grid_fundamentals/the_sinusoidal_construction.md)
is where those coordinates come from.

```{code-cell} python
itacart.cell_to_anchor(cell), itacart.cell_to_centroid(cell)
```

The two points are different and both are meaningful — the anchor is the
lattice corner, the centroid is the cell's representative position.
[Cells, Vertices, and Spatial Identity](../grid_fundamentals/cells_vertices_and_spatial_identity.md)
is why the distinction is load-bearing.

```{figure} ../../../_static/f4/f4_10_anchor_and_centroid.png
:alt: Anchor and centroid of a cell
:width: 100%

`docs/_static/f4/f4_10_anchor_and_centroid.png`
```

## As a Shapely polygon

When you want to intersect, buffer or export, use `cell_to_polygon` polygon
rather than the boundary ring:

```{code-cell} python
poly = itacart.cell_to_polygon(cell)
poly.is_valid, poly.geom_type
```

`cell_to_polygon` returns a valid Shapely polygon for the ordinary,
border-absorbing, and polar cell families.

There is one important convention near the antimeridian. Cells that
belong to an extension zone may use longitudes greater than 180°
instead of wrapping immediately to the opposite side of the globe.
Keeping those coordinates continuous prevents the polygon from
being split or folded across the ±180° seam.

For example, a cell extending eastward from 179° may continue to
181° rather than being represented as 179° to −179°. In that
continuous representation, Shapely sees the intended local polygon
instead of a geometry spanning almost the entire longitude range.

See [Extension Zones](../boundaries/extension_zones.md) covers where those
cells are.

## The projected side

The same cell exists on the sinusoidal plane, and that is where the lattice is
regular:

```{code-cell} python
itacart.cell_to_sinusoidal(cell)
```

Going the other way, `sinusoidal_to_cell` quantizes a projected position. The
pair is the seam between the regular grid and the curved surface, and most of
the package's arithmetic happens on the plane side of it.

## Vertex count is not constant

```{code-cell} python
from itacart import boundary
[(c, boundary.cell_shape(c), len(itacart.cell_to_boundary(c, close=True)))
 for c in ("SE(1400/0374)", "NE(0000/0500)", "SE(1930/0196)")]
```

A parallelogram and a trapezoid close with five positions; the meridian
triangle closes with four. Code that assumes four corners will be wrong at the
boundary, which is why `cell_shape` exists and why
[Polar & Border Cell Families](polar_and_border_cell_families.md) is a page
rather than a footnote.
