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

# Area Estimation & Error Bounds

Area from a cell set is a count times a constant. The question this page
answers is how close that count is to the parcel's true area, and what bounds
the difference.

```{code-cell} python
import itacart
from itacart import geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
n = geometry.count_internal_cells(parcel, 9)
n, n * itacart.nominal_cell_area(9)
```

## The estimate converges

```{code-cell} python
[(r, geometry.count_internal_cells(parcel, r) * itacart.nominal_cell_area(r))
 for r in (5, 6, 7, 8, 9, 10, 11)]
```

Sixty thousand square metres at resolution 5, forty thousand at 6, then 44 800,
44 900, 45 127, 45 231 and 45 231.05. The coarse levels are wild — the parcel
is only a few cells across — and by resolution 9 the estimate has settled to
within a quarter of a per cent of where it lands at 11.

```{figure} ../../../_static/f7/f7_03_area_convergence.png
:alt: Area estimate converging with resolution
:width: 100%

`docs/_static/f7/f7_03_area_convergence.png`
```

Note the non-monotonicity between 5 and 6. A cover is not a shrinking
approximation from above or below; cells are won and lost along the boundary as
the grid shifts under it, and at coarse resolutions a single cell is a large
fraction of the answer.

## What bounds the error

The error lives entirely on the boundary. Interior cells are exactly right, so
the difference between the estimate and the truth is at most the area of the
cells the outline passes through — the count of boundary cells times the cell
area.

That gives the practical rule: **error scales with the perimeter, not the
area.** Doubling a parcel's linear size quadruples its area and only doubles
its boundary, so a large parcel is estimated proportionally better than a small
one at the same resolution.

And refining one odd level divides the cell area by twenty-five while
multiplying the boundary count by about five, so the absolute error falls by
roughly a factor of five per odd level.

```{figure} ../../../_static/f7/f7_04_area_error_by_latitude.png
:alt: Area error across latitude
:width: 100%

`docs/_static/f7/f7_04_area_error_by_latitude.png`
```

## Why `center` is the mode for area

The default containment mode is the unbiased one: a cell whose centre is inside
is counted whole, a cell whose centre is outside is dropped whole, and along a
boundary the two errors trade against each other.

```{code-cell} python
{mode: itacart.count_cells(geometry.polyfill(parcel, 7, containment=mode))
 for mode in ("contains", "center", "intersects")}
```

`contains` understates systematically and `intersects` overstates
systematically. Both are useful — as bounds. The pair brackets the true area,
and the bracket is a defensible thing to report when a single number is not:

```{code-cell} python
lo = itacart.count_cells(geometry.polyfill(parcel, 9, containment="contains"))
hi = itacart.count_cells(geometry.polyfill(parcel, 9, containment="intersects"))
lo * itacart.nominal_cell_area(9), hi * itacart.nominal_cell_area(9)
```

## Where the arithmetic stops being a multiplication

Everything above assumes every cell carries the nominal area. Near the
antemeridian and at the poles they do not, and the count-times-constant rule
silently gives a wrong number.

```{code-cell} python
def parcel_area(cells, resolution):
    """Nominal where the guarantee holds, measured where it does not."""
    return sum(itacart.nominal_cell_area(resolution) if itacart.is_equal_area_cell(c)
               else itacart.effective_cell_area(c)
               for c in cells)

parcel_area(itacart.decompose(geometry.polyfill(parcel, 7)), 7)
```

For an interior parcel this is the same number the multiplication gives, at
more cost. Near a seam it is the only right one.
[Known Boundary Exceptions](../standards_and_conformance/known_boundary_exceptions.md).

## What none of this bounds

The error of the source. A cover inherits whatever the surveyed boundary got
wrong, and no resolution repairs it —
[Choosing Resolution from Mapping Scale](choosing_resolution_from_mapping_scale.md).
The bounds on this page are the grid's contribution to the error, not the
total.
