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

# Meridian and Equator

The prime meridian and the equator divide the globe into the four quadrants,
and each quadrant's lattice is a mirror of the others. That symmetry is why
ITACaRT never uses a negative index — but it leaves a seam down the middle.

## The problem the mirror creates

Applying the same parallelogram geometry to each quadrant produces
discontinuities along the quadrant's own meridian. The cells on either side do
not meet: reflected across the meridian, the lean of one is the mirror of the
other, and the two leave a gap that neither fills.

## The triangle

The solution keeps the area and changes the shape. A cell straddling the prime
meridian is an **isosceles triangle** whose base is twice its height, with the
cell's index sitting at the midpoint of the base. The figure is therefore
mirrored across the meridian, and it tiles the gap exactly.

```{figure} ../../../_static/f4/f4_01_grid_triangle.png
:alt: The grid triangle at the prime meridian
:width: 80%

`docs/_static/f4/f4_01_grid_triangle.png`
```

```{figure} ../../../_static/fpaper/fpaper_04_grid_triangle_and_prime_meridian.png
:alt: The grid triangle and prime meridian behaviour under both refinements
:width: 100%

Figure 4 of the paper: the triangle, and how it refines on each side of
the meridian.
`docs/_static/fpaper/fpaper_04_grid_triangle_and_prime_meridian.png`
```

A triangle keeps the base-and-height behaviour of a parallelogram, so the area
survives the change of shape:

```{code-cell} python
import itacart
from itacart import boundary
tri = "NE(0000/0500)"
boundary.cell_shape(tri), itacart.is_equal_area_cell(tri)
```

Equal-area, like any interior cell. The vertex count is what differs:

```{code-cell} python
len(itacart.cell_to_boundary(tri, close=True))
```

Four positions — three distinct vertices plus the closing repeat.

```{figure} ../../../_static/f4/f4_02_triangle_tiles_its_row.png
:alt: The triangle tiling its row
:width: 100%

`docs/_static/f4/f4_02_triangle_tiles_its_row.png`
```

## Indexed once, in the east

The triangle spans both sides of the meridian, so indexing it in both quadrants
would name the same ground twice. It is indexed in the **eastern** quadrants
only. At resolution 1 and coarser, the western quadrants simply have no
meridian column:

```{code-cell} python
from itacart import index as ix
from itacart.exceptions import NonExistentCellError
ix.is_valid_index("NW(0000/0500)")
```

Valid, and nonexistent:

```{code-cell} python
try:
    itacart.cell_to_centroid("NW(0000/0500)")
except NonExistentCellError as exc:
    print(exc)
```

That is the case [Cell Validity & Existence](../index_and_hierarchy/cell_validity_and_existence.md)
is built around, and this is where it comes from.

## The triangle is inherited

At finer resolutions the indexing does not create separate western cells along
this boundary either — they are part of the hierarchical subdivision of the
adjacent eastern cells. A meridian triangle therefore refines into triangles,
all the way down:

```{figure} ../../../_static/f4/f4_03_triangle_refinement.png
:alt: A meridian triangle refining into triangles
:width: 100%

`docs/_static/f4/f4_03_triangle_refinement.png`
```

This is the opposite of what happens at the outer edge, where the trapezoid is
*not* inherited — [Absorbing Boundary Cells](absorbing_boundary_cells.md) is
the contrast.

## The equator

The equator divides north from south and needs no special cell. Rows are
counted outward from it in each quadrant, so row 0 sits against the equator and
the row index grows towards the pole. The mirror across the equator is exact,
which is why a northern cell and its southern twin share a base angle:

```{code-cell} python
from itacart import metrics
metrics.cell_base_angle("NE(1400/0374)"), metrics.cell_base_angle("SE(1400/0374)")
```

## Stepping across a seam

A step whose lexical target falls outside its quadrant's range does not fail —
it deflects into the neighbouring quadrant:

```{code-cell} python
from itacart import topology
topology.deflect("NE(0000/0500)", "W")
```

[Neighbors & Grid Traversal](../spatial_operations/neighbors_and_grid_traversal.md)
covers when a traversal needs this and when the ordinary arithmetic suffices.
