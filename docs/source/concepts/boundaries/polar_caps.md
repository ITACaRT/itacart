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

# Polar Caps

The grid does not stop at high latitude — it thins. Rows continue towards the
pole, and each holds fewer cells than the one below it, until almost nothing is
left.

```{code-cell} python
from itacart import boundary
[(row, boundary.last_lattice_column("SE", row, 10000.0))
 for row in (0, 300, 600, 900, 950, 999)]
```

## Why it thins

The construction holds base and height constant on the ellipsoid. Towards the
pole a parallel is shorter, so fewer cells of that base fit along it. The
lattice narrowing is the sinusoidal projection doing exactly what it is
supposed to — [The Sinusoidal Construction](../grid_fundamentals/the_sinusoidal_construction.md)
is where that comes from.

```{figure} ../../../_static/f9/f9_23_lattice_narrowing_to_the_pole.png
:alt: The last addressable column per row, from the equator to the pole
:width: 100%

Columns per row, with the last fifteen rows inset — the tail is the whole
point and is invisible at full scale.
`docs/_static/f9/f9_23_lattice_narrowing_to_the_pole.png`
```

```{figure} ../../../_static/f4/f4_16_quadrant_limits.png
:alt: The limits of a quadrant's lattice
:width: 100%

`docs/_static/f4/f4_16_quadrant_limits.png`
```

## The last row

At resolution 1 the final row holds **six cells in the entire world**:

```{code-cell} python
import itacart
census = []
for quadrant in ("NE", "NW", "SE", "SW"):
    last = boundary.last_lattice_column(quadrant, 999, 10000.0)
    for column in (0, last):
        candidate = f"{quadrant}({column:04d}/0999)"
        try:
            census.append((candidate, boundary.cell_shape(candidate),
                           boundary.absorbs_border(candidate)))
        except itacart.exceptions.ITACaRTError:
            pass
census
```

Four trapezoidal absorbers, one per quadrant in the last existing column, and
two prime-meridian triangles — eastern quadrants only, because the western
quadrants have no meridian column at all.

```{figure} ../../../_static/f4/f4_13_polar_row.png
:alt: The last row of the lattice
:width: 100%

`docs/_static/f4/f4_13_polar_row.png`
```

Six is the enumerated total of that row, not a sample of it.

## Row 1000: triangles with a different area rule

One row remains above it, and it holds two cells in the whole world:

```{code-cell} python
[(c, boundary.cell_shape(c), boundary.absorbs_border(c), itacart.is_equal_area_cell(c))
 for c in ("NE(0000/1000)", "SE(0000/1000)")]
```

Both are **triangles**, both **absorb the border**, and neither carries the
nominal area. Everywhere else a triangle is a prime-meridian cell that keeps
its area exactly; here the same shape name attaches to a cell that does not.

For code that needs to distinguish these cases, use `is_equal_area_cell`,
rather then inferring equal area from `shape != "trapezoid"`.

## The poles themselves are addressable

Thinning is not stopping. The poles themselves can be indexed:

```{code-cell} python
north = itacart.geo_to_cell(0.0, 90.0, resolution=9)
south = itacart.geo_to_cell(0.0, -90.0, resolution=9)
```

These positions belong to the special polar family rather than to the
ordinary parallelogram lattice. Their geometry and refinement should
therefore be handled through the same public cell APIs used elsewhere,
without assuming the ordinary-cell shape or nominal area.

## Why document the polar family?

The polar caps are outside the ordinary terrestrial cadastral use for
which ITACaRT is designed. They still matter because ITACaRT is a global
grid: every geographic position must have a well-defined representation,
and global operations, validation sweeps, and conformance tests eventually
reach the poles.

The purpose of documenting this family is therefore completeness of the
DGGS and correctness of the API, not a special cadastral workflow for the
polar regions.
