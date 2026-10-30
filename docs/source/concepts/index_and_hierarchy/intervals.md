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

# Intervals

A region as sorted ranges of integers, with set operations. The Multi-Order
Coverage maps of the IVOA store a region of the sky as half-open ranges of
cell numbers at one depth, and combine regions by merging those lists. The
same works on the [integer index](integer_index.md): the overlap of two
parcels, or of an obstacle and a flight corridor, is one intersection of two
lists, in time linear in their lengths.

```{code-cell} python
import itacart

parcel = "NE(0625/0451(1,2))"
corridor = "NE(0625/0451(2,3))"
a = itacart.cells_to_intervals(parcel)
b = itacart.cells_to_intervals(corridor)
a, itacart.intervals_to_cells(itacart.intervals_intersection(a, b))
```

## The form

An interval list is a list of `(start, stop)` pairs, sorted, disjoint and with
no two touching, every integer at one resolution: the list's *depth*. The
integers of a cell's descendants are consecutive, so away from the border a
cell coarser than the depth is a single range:

```{code-cell} python
itacart.cells_to_intervals("NE(0625/0451)", resolution=4)
```

Here a hundred cells of resolution 4 make one range. The functions work on that
form and return it, with no new type: an interval list is another
representation of a region, not another spelling of an index.

```{code-cell} python
itacart.intervals_count(itacart.cells_to_intervals("NE(0625/0451(3))", resolution=4))
```

## The depth is chosen

The depth defaults to the finest resolution in the index, and it can be set.
It is not fixed at resolution 13, because the region's semantics is physical,
and on the border that has a price. A cell whose outer side was carried onto
the border has children spelled under the next column's prefix, outside the
range its own integers span. Such a cell is expanded through its physical
children, down to the depth. Under one resolution-1 border cell, the
absorbing descendants number about one per row at every resolution: a
thousand at resolution 7, and a million at 13. Pinned to resolution 13, one
border cell would become a million ranges.

```{code-cell} python
border = "NE(2003/0000)"
len(itacart.cells_to_intervals(border, resolution=2)), list(
    itacart.get_children(border, flatten=True)
)
```

The trapezoid takes two ranges, because one of its five children is spelled
under `NE(2004/0000)`. Two lists of different depths are brought to the finer
one before they are combined, which is exact, because a coarser list names
whole cells.

## Back to an index

From intervals back to an index, the answer is what
{func}`itacart.compact_cells` writes, so its TreeBlob is the one
{func}`itacart.serialize_to_blob` writes for the region:

```{code-cell} python
from itacart import serialization

region = "NE(0625/0451(1(A1,A2,A3)),0626/0451)"
back = itacart.intervals_to_cells(itacart.cells_to_intervals(region))
back == itacart.compact_cells(region), (
    serialization.serialize_to_blob(back) == serialization.serialize_to_blob(region)
)
```

The operations are union, intersection, difference, containment, overlap and
count. Each takes two lists and returns a list or an answer, never a string,
and an interval that crosses resolutions, or a list of mixed depths, is
refused rather than guessed.
