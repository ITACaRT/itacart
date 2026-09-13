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

# Absorbing Boundary Cells

Where the lattice meets the edge of its domain, the outermost cell does not fit.
The rule is that it takes up the remainder: its outer side is constrained to
the boundary line, and the parallelogram becomes a trapezoid with a longer base.

```{code-cell} python
import itacart
from itacart import boundary
cell = "SE(1930/0196)"
(boundary.cell_shape(cell), boundary.absorbs_border(cell),
 itacart.is_equal_area_cell(cell))
```

```{figure} ../../../_static/f4/f4_08_trapezoid_absorption.png
:alt: A parallelogram absorbing the border into a trapezoid
:width: 100%

`docs/_static/f4/f4_08_trapezoid_absorption.png`
```

```{figure} ../../../_static/f9/f9_22_absorbing_cell_against_interior.png
:alt: An absorbing trapezoid beside an interior parallelogram of the same row
:width: 100%

The same row, the same nominal area, drawn on the sinusoidal plane where the
shapes are regular. `docs/_static/f9/f9_22_absorbing_cell_against_interior.png`
```

```{figure} ../../../_static/fpaper/fpaper_06_trapezoidal_cells.png
:alt: Examples of trapezoidal cells
:width: 80%

Figure 6 of the paper: the trapezoidal cells the rule produces.
`docs/_static/fpaper/fpaper_06_trapezoidal_cells.png`
```

The rule is stated narrowly in the specification: it applies when a vertex of
the parallelogram on the side opposite the prime meridian exceeds the boundary,
and it applies **to that cell only**, not to all subsequent resolutions.

## They do not carry the nominal area

This is the one family that breaks the equal-area guarantee, and it is the
reason the guarantee is stated conditionally everywhere else:

```{code-cell} python
round(itacart.effective_cell_area(cell)), round(itacart.nominal_cell_area(1))
```

About half. Not a rounding error — half a cell, which in cadastral terms is
fifty square kilometres at resolution 1.

```{figure} ../../../_static/f4/f4_17_trapezoid_areas.png
:alt: The spread of trapezoidal cell areas
:width: 100%

`docs/_static/f4/f4_17_trapezoid_areas.png`
```

```{figure} ../../../_static/f4/f4_09_absorption_across_latitude.png
:alt: How absorption varies with latitude
:width: 100%

`docs/_static/f4/f4_09_absorption_across_latitude.png`
```

The guard is never the shape name. Ask the construction:
[Effective vs Nominal Geometry](../geometry/effective_vs_nominal_geometry.md).

## Where the edge is

The lattice narrows towards the poles, so the last existing column is a
function of the row:

```{code-cell} python
[(row, boundary.last_lattice_column("SE", row, 10000.0))
 for row in (0, 196, 500, 900, 999)]
```

Two thousand and three columns at the equator, one at the last row. That
narrowing is what [Polar Caps](polar_caps.md) is about.

## The trapezoid is not inherited

A meridian triangle refines into triangles. An absorbing trapezoid does not
refine into trapezoids, and its children do not number the refinement ratio:

```{code-cell} python
len(list(itacart.get_children(cell, flatten=True)))
```

Three where the rule says four, because the border truncates the parent's own
footprint. A cell that absorbs has been measured to hold between two and six
children, which is why the package enumerates children rather than multiplying.

```{figure} ../../../_static/f5/f5_07_trapezoid_children.png
:alt: The children of a trapezoidal cell
:width: 100%

`docs/_static/f5/f5_07_trapezoid_children.png`
```

```{figure} ../../../_static/f5/f5_06_border_children_counts.png
:alt: Child counts along the border
:width: 100%

`docs/_static/f5/f5_06_border_children_counts.png`
```

## Children reach past their parent

The border is a smooth convex curve, and an absorbing cell replaces its outer
side by the chord of that curve across its own height. A parent spans one cell
height and takes one chord; each child spans half that height and takes two —
and a chord over a shorter span lies further out.

So a child of an absorbing cell lies partly outside its parent's effective
ring, by construction. This is a declared limitation of this version rather
than a defect awaiting repair, and it is stated in full, with its measurements,
at [Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md).

```{figure} ../../../_static/f5/f5_03_parent_prefix_divergence.png
:alt: Where a child diverges from its parent's ring
:width: 100%

`docs/_static/f5/f5_03_parent_prefix_divergence.png`
```

## Working near the border

Test `is_equal_area_cell` before any area arithmetic. Enumerate children rather
than assuming a count. Expect grid distances to stop behaving like a square
lattice — [Grid Distance](../spatial_operations/grid_distance.md) says why. And
do not reingest an exported border outline as fresh geographic geometry:
[Geographic vs Grid Geometry](../data_and_interoperability/geographic_vs_grid_geometry.md)
measures what that costs.
