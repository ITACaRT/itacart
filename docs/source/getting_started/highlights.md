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

# Highlights

Four properties distinguish this grid. Every number on this page is computed
when the documentation is built, so a claim that stops being true stops the
build rather than surviving in print.

```{code-cell} python
import itacart
```

## Cells of equal area, on the ellipsoid

A cell's area is the area the specification assigns to its resolution, not an
average that holds somewhere and fails elsewhere. The value is a round metric
quantity, which is what makes one cell a usable unit of account:

```{code-cell} python
itacart.nominal_cell_area(3), itacart.nominal_cell_area(9), itacart.nominal_cell_area(13)
```

One square kilometre, one square metre, one square centimetre. Counting cells
is therefore measuring area, with no conversion factor in between.

The guarantee is not universal, and the exceptions are specified rather than
incidental. Cells along the antemeridian extension edge are trapezoidal and do
not hold the area, so an area-sensitive computation asks first:

```{code-cell} python
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
itacart.is_equal_area_cell(cell)
```

[Equal Area](../concepts/grid_fundamentals/equal_area.md) establishes the
property and states where it is qualified;
[Effective vs Nominal Geometry](../concepts/geometry/effective_vs_nominal_geometry.md)
covers what to do when it does not hold.

## A decimal hierarchy

The resolution ladder is decimal, so a resolution corresponds to a length a
surveyor already thinks in rather than to an arbitrary subdivision count:

```{code-cell} python
table = itacart.resolution_table()
len(table), table[1]["cell_size_m"], table[-1]["cell_size_m"]
```

Fourteen entries. Resolution 0 addresses the quadrant and carries no metric
size; the sized levels run from ten kilometres at resolution 1 down to one
centimetre at resolution 13, the two values printed above. Levels alternate
between refining one cell into four and into twenty-five, which is what keeps
the areas decimal:

```{code-cell} python
itacart.constants.REFINEMENT_RATIO[2:8]
```

[Resolution & Scale](../concepts/grid_fundamentals/resolution_and_scale.md)
relates a resolution to a mapping scale.

## Cartesian-like addressing

Within a quadrant the base grid is addressed by row and column, so the
intuition a surveyor already has about a plane grid keeps working. The
resolution-1 component of an index is literally a pair of numbers:

```{code-cell} python
itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=1)
```

## One string for a whole feature

An index addresses a set of cells as easily as a single one, and functions
accept either. A composite index keeps its cells in order, and results come
back aligned with that order:

```{code-cell} python
region = "NE(0001/0002(1(A1,A2,A3)))"
itacart.count_cells(region), itacart.get_parent(region)
```

That alignment is the whole point of the notation: a parcel stays one object
through a pipeline instead of decomposing into a list that has to be carried
alongside its own metadata.
[Composition & Decomposition](../concepts/index_and_hierarchy/composition_and_decomposition.md)
covers the round trip.
