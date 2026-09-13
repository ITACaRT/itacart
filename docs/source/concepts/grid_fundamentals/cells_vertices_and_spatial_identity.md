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

# Cells, Vertices, and Spatial Identity

A cell is an area. A vertex is a point. An index names a cell, but the number
the index literally encodes is a vertex. Keeping those apart is the difference
between an address that works and one that quietly means two things.

```{code-cell} python
import itacart
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
cell
```

## The anchor: the point the index encodes

Descending the index arrives at a lattice corner. That corner is the cell's
anchor, and it is what the address inverts to exactly:

```{code-cell} python
itacart.cell_to_anchor(cell)
```

The anchor is the lower-left vertex (at NE quadrant), and choosing it is a usability decision.
Addressing by a corner is how plane grids have always worked, so a surveyor
reading an ITACaRT index is reading something already familiar — row, column,
then refinement. The paper names this vertex for the standard's requirement
that every cell carry a designated representative position.

## The centroid: the point the cell answers with

The requirement asks for a position *within* the cell. A vertex lies on the
boundary, not inside, so this package does not use the anchor for that role. It
designates the centroid:

```{code-cell} python
itacart.cell_to_centroid(cell)
```

This is the point every single-point route answers with, and it is what the
engine declares when it describes a cell to a standards-facing client. The
paper anticipated the move — it offers the vertex for Cartesian familiarity and
notes that an implementation may represent the centroid instead — and this
implementation takes that option.

So both points exist and both are meaningful. The anchor is the lattice corner
the descent inverts and the point at which corner-defined metrics are measured;
the centroid is the cell's representative position. Neither is a fallback for
the other.

## Why the distinction is load-bearing

Three things go wrong when the two are collapsed.

**A boundary point would belong to no cell, or to four.** The anchor of one
cell is a corner shared with its neighbours. If the address pointed at the
anchor *as* the cell's location, then asking what is at that location has four
equally good answers. The cell is the area; the anchor merely labels it.

**Area computations would be defined at the wrong place.** A corner is where
some shape metrics are evaluated, because that is where the lattice is exact.
An area is a property of the enclosed figure. Mixing them produces numbers that
look right and are measured at the wrong point.

**Interoperability would carry the wrong geometry.** An exported feature has to
say which of the two it is giving you.
[Geographic vs Grid Geometry](../data_and_interoperability/geographic_vs_grid_geometry.md)
is where that is stated, and it is stated because the two are not the same.

## The whole cell

The boundary comes back as a ring, and the ring is the cell:

```{code-cell} python
ring = itacart.cell_to_boundary(cell, close=True)
len(ring)
```

Four distinct vertices plus a repeat of the first, which is what closes the
ring. [Cell Geometry](../geometry/cell_geometry.md) derives each of them from
the index.
