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

# Equal Area

Every cell at a given resolution encloses the same area — except at the
antemeridian, where it does not, and the exception is written into the
specification rather than discovered later.

```{code-cell} python
import itacart
itacart.nominal_cell_area(3), itacart.nominal_cell_area(9), itacart.nominal_cell_area(13)
```

One square kilometre, one square metre, one square centimetre.

## What the guarantee rests on

Not on measurement. The cells have constant base and height, so equal
area is a property of the construction itself —
[The Sinusoidal Construction](the_sinusoidal_construction.md) provides
the geometric argument. Measurement does not establish the equal-area
property; it checks that the implementation delivers the construction
it claims.

That distinction matters for how you read a number. `nominal_cell_area` returns
the value the specification assigns to a resolution; it takes a resolution, not
a cell, because under the guarantee the cell cannot matter:

```{code-cell} python
import inspect
inspect.signature(itacart.nominal_cell_area)
```

## Where it stops

Three cell shapes exist in this grid, and only two of them hold the area:

```{code-cell} python
itacart.constants.CELL_SHAPES
```

Interior cells are parallelograms. Cells straddling the prime meridian are
isosceles triangles, chosen because a triangle has the same base-and-height
behaviour as a parallelogram — the area survives. Cells along the antemeridian
extension edge are trapezoids, and they do not.

The paper is explicit about this in its conformance table: the equal-size
constraint is met only partially, because the trapezoidal cells at the
antemeridian do not hold their dimensions, and it argues that for land cadastre
this is acceptable. It is acceptable because the affected cells lie along one
meridian across two inhabited landmasses, and because the grid says so rather
than hiding it.

[Extension Zones](../boundaries/extension_zones.md) covers where those cells
are; [Polar & Border Cell Families](../geometry/polar_and_border_cell_families.md)
covers what they look like.

## Asking before you count

Because the guarantee is conditional, an area-sensitive computation should test
rather than assume:

```{code-cell} python
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
itacart.is_equal_area_cell(cell)
```

When the answer is `False`, the nominal value is the wrong number and the
measured one is available instead:

```{code-cell} python
itacart.effective_cell_area(cell), itacart.nominal_cell_area(13)
```

For an interior cell the two agree, which is the guarantee doing its job. Where
they diverge, [Effective vs Nominal Geometry](../geometry/effective_vs_nominal_geometry.md)
explains which to use and why the package refuses to quietly substitute one for
the other.

## The consequence worth having

Counting cells is measuring area, with no conversion factor:

```{code-cell} python
from shapely.geometry import Polygon
parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
n = itacart.count_internal_cells(parcel, 9)
n, n * itacart.nominal_cell_area(9)
```

Square metres, directly. That is what makes a cell usable as a unit of account
— the paper's tokenization case rests on it, and so does the simpler and more
common case of a tax assessment that has to be reproducible by whoever checks
it. [Area Estimation & Error Bounds](../cadastral_workflows/area_estimation_and_error_bounds.md)
bounds how close that count is to the parcel's true area.
