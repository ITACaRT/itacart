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

# Known Boundary Exceptions

Every qualification on the equal-area claim, in one place, each with its scope.

```{code-cell} python
import itacart
itacart.describe()["constraints"]
```

`cellEqualSized` is asserted, and the three exception families are named in the
same breath. That is the pattern this whole section follows: the claim and its
scope travel together.

## The three families

**Trapezoids.** Cells whose outer side has been carried onto the domain border
— the antemeridian, or the meridian edge of an extension zone. Their longer
base is clipped, so they do not carry the nominal area. There are 3 941 of them
at resolution 1. [Absorbing Boundary Cells](../boundaries/absorbing_boundary_cells.md).

**Polar caps.** Two of them, one per pole. [Polar Caps](../boundaries/polar_caps.md).

**The polar row.** The row above the last full one, where the lattice has
thinned to almost nothing.

Together: 3 943 cells, **0.077 per cent of the grid** at resolution 1.

## What the exception costs, per family

The trapezoid family runs from 0.032 to 2.067 of the theoretical average cell
area. The caps sit at 0.121. Everything else — parallelograms and meridian
triangles alike — sits at 0.999915, which is 0.0085 per cent off the
theoretical average.

So the spread within the exception families is large, and the spread outside
them is negligible. An average taken across the whole grid would hide both
facts.

```{figure} ../../../_static/f4/f4_17_trapezoid_areas.png
:alt: The spread of trapezoidal cell areas
:width: 100%

`docs/_static/f4/f4_17_trapezoid_areas.png`
```

## The meridian triangle is not an exception

Worth stating because the shape suggests otherwise. A prime-meridian cell is a
triangle rather than a parallelogram, and it carries the nominal area exactly:

```{code-cell} python
from itacart import boundary
[(c, boundary.cell_shape(c), itacart.is_equal_area_cell(c))
 for c in ("NE(0000/0500)", "SE(1930/0196)", "NE(0000/1000)")]
```

A triangle at row 500 keeps its area; a trapezoid at the border does not; and a
triangle at row 1000 — a polar cap — does not either. The shape name is not the
exception; the construction is.

## Where the exceptions surface

They are not confined to area arithmetic. The same cells are where:

- child counts stop matching the refinement ratio —
  [Parent, Children & Descendants](../index_and_hierarchy/parent_children_and_descendants.md);
- a child reaches past its parent's ring —
  [Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md);
- grid distance stops behaving like a square lattice —
  [Grid Distance](../spatial_operations/grid_distance.md);
- a reingested outline covers a third of the cell —
  [Geographic vs Grid Geometry](../data_and_interoperability/geographic_vs_grid_geometry.md).

One family of cells, four consequences. If your work touches the antemeridian
or the poles, read all four before trusting a result.

## The guard

```{code-cell} python
itacart.is_equal_area_cell("SE(1930/0196)")
```

One call, and it is the only reliable test. Not the shape name, not the row
index, not a latitude threshold.
