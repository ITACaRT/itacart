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

# Polar & Border Cell Families

Three cell shapes exist in this grid, and knowing which one you hold changes
the vertex count, the area computation and whether the equal-area guarantee
applies.

```{code-cell} python
import itacart
itacart.constants.CELL_SHAPES
```

## Parallelogram — the interior

The ordinary case, and overwhelmingly the common one. Four vertices, nominal
area, the lean measured by
[Base Angle & Shape](base_angle_and_shape.md).

## Triangle — the prime meridian

A cell whose index falls on the prime meridian is an isosceles triangle,
mirrored across the meridian, with a base twice its height. The shape is chosen
because a triangle keeps the base-and-height behaviour of a parallelogram, so
the area survives:

```{code-cell} python
from itacart import boundary
tri = "NE(0000/0500)"
boundary.cell_shape(tri), itacart.is_equal_area_cell(tri), len(itacart.cell_to_boundary(tri, close=True))
```

Four positions rather than five, because a triangle closes with three distinct
vertices.

```{figure} ../../../_static/f4/f4_01_grid_triangle.png
:alt: The grid triangle at the prime meridian
:width: 80%

`docs/_static/f4/f4_01_grid_triangle.png`
```

**The triangle is inherited.** A meridian triangle refines into triangles at
every resolution, because the indexing does not create separate western cells
along that boundary — they are part of the hierarchical subdivision of the
adjacent eastern cells. That is also why a western cell with an X index of zero
does not exist at all:
[Cell Validity & Existence](../index_and_hierarchy/cell_validity_and_existence.md).

## Trapezoid — the absorbing edge

Where the lattice meets the domain border, the outer cell absorbs the
remainder. Its outer side becomes a chord of the border curve, and it is the
one family that does not carry the nominal area:

```{code-cell} python
abs_cell = "SE(1930/0196)"
(boundary.cell_shape(abs_cell), boundary.absorbs_border(abs_cell),
 itacart.is_equal_area_cell(abs_cell))
```

**The trapezoid is not inherited.** Unlike the triangle, a trapezoid's children
are not all trapezoids, and their count is not the refinement ratio — this cell
holds three children where the rule says four. It is also why a child of an
absorbing cell reaches past its parent's ring, which is stated in full at
[Lexical Ancestry vs Physical Refinement](../index_and_hierarchy/lexical_ancestry_vs_physical_refinement.md).

```{figure} ../../../_static/f4/f4_08_trapezoid_absorption.png
:alt: How a trapezoidal cell absorbs the border
:width: 80%

`docs/_static/f4/f4_08_trapezoid_absorption.png`
```

## The poles

At high latitude the lattice runs out. The last row of the grid is very nearly
empty — at resolution 1 it holds six cells in the whole world: four absorbers,
one per quadrant in the last existing column, and two prime-meridian triangles,
which exist in the eastern quadrants only because the western quadrants have no
meridian column.

```{figure} ../../../_static/f4/f4_13_polar_row.png
:alt: The last row of the lattice
:width: 80%

`docs/_static/f4/f4_13_polar_row.png`
```

A polar cell is still *named* a triangle by `cell_shape`, and it still absorbs
the border — which is exactly why the guard everywhere in this documentation is
`is_equal_area_cell` and never `shape != "trapezoid"`.

## The rule to carry away

Ask `cell_shape` when you care about vertices. Ask `is_equal_area_cell` when
you care about area. They answer different questions and the answers do not
line up.
