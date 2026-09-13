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

# Why a Cadastral DGGS?

Cadastral mapping defines where one person's property stops and another's
begins. It is the foundation a Land Administration System rests on — legal
ownership, taxation, urban planning — and it is measurement with consequences.

So what does a discrete global grid buy that a coordinate pair does not?

## A coordinate is a position; a cell is a piece of the world

A parcel stored as a ring of coordinates is a description of a boundary. To
answer "do these two parcels overlap", "what is this parcel's area", "is this
point inside", you compute — every time, against floating-point coordinates,
with whatever tolerance your library happens to use.

A parcel stored as a set of cells is a set of named pieces. Overlap is set
intersection. Containment is membership. Area is a count.

```{code-cell} python
import itacart
from itacart import geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
cover = geometry.polyfill(parcel, resolution=7)
itacart.count_cells(cover), itacart.count_cells(cover) * itacart.nominal_cell_area(7)
```

Four hundred and forty-eight cells, and forty-four thousand eight hundred
square metres — a count multiplied by a constant, with no conversion factor
and no tolerance parameter.

## Why the areas being decimal matters

They are round metric quantities on purpose: one square kilometre at resolution
3, one square metre at 9, one square centimetre at 13. That makes one cell one
unit of area, which is what lets a cell count stand in for an area in a
register, in an assessment, or in a token.

The paper's blockchain case rests on exactly this — a token can correspond to a
standard metric area, so a parcel becomes a quantity of tokens rather than a
geometry someone has to re-measure. The simpler and far more common case is a
tax assessment that whoever checks it can reproduce by counting.

[Equal Area](../grid_fundamentals/equal_area.md) is the guarantee this rests on.

## Why the ellipsoid rather than a sphere

Most discrete global grids project a polyhedron onto a sphere. It is cheaper
and it is right for most applications. It is wrong here twice: the sphere is
not the figure survey control is referred to, and the projection distorts area
in a way that varies across the globe.

A parcel's area carries tax and legal weight, so area is the quantity that has
to be preserved rather than approximated.
[Design Decisions](../standards_and_conformance/design_decisions.md) is that
trade written out.

## Why the addressing looks Cartesian

Because surveyors already think in rows and columns. The resolution-1 component
of an index is literally a pair of integers, and the refinements are codes on a
two-by-two or five-by-five block. Nothing has to be relearned:

```{code-cell} python
itacart.geo_to_cell(-46.6329, -23.5509, resolution=1)
```

Usability was a design criterion, not a nicety — a representation that
practitioners will not adopt does not administer any land.

## What it does not buy

A DGGS does not make an imprecise survey precise. The cells are exact; what you
fill them from is whatever your source mapping was, and the cover inherits its
error. [Choosing Resolution from Mapping Scale](choosing_resolution_from_mapping_scale.md)
is about not pretending otherwise.
