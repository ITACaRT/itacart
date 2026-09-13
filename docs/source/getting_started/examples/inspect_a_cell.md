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

# Inspect a Cell

Take a cell apart — its address, its place in the hierarchy, its shape family
and its measured geometry.

```{code-cell} python
import itacart
from itacart import index as ix, boundary, metrics

cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
cell
```

## The address, component by component

```{code-cell} python
ix.split_components(cell)
```

Quadrant, base cell, then one refinement code per level. The two shortcuts save
parsing it yourself:

```{code-cell} python
ix.quadrant_of(cell), ix.base_cell_of(cell)
```

Note the alternating alphabets: single digits on even levels, letter-digit
pairs on odd ones. Which alphabet is legal at a position is fixed by the
position, not chosen:

```{code-cell} python
ix.is_valid_index("NE(0001/0002(1))"), ix.is_valid_index("NE(0001/0002(A1))")
```

## Where it sits in the hierarchy

```{code-cell} python
parent = itacart.get_parent(cell)
parent, itacart.child_position(cell)
```

The parent is the address with its last component removed — ascent is
truncation, not computation. The ancestors, coarse to fine:

```{code-cell} python
itacart.get_ancestors(cell)
```

Descending gives the siblings this cell belongs to. Pass `flatten=True`, or you
get one list per input cell rather than the cells themselves:

```{code-cell} python
siblings = list(itacart.get_children(parent, flatten=True))
len(siblings), siblings[0]
```

Twenty-five, because resolution 9 is an odd level.

## Which family it belongs to

```{code-cell} python
boundary.cell_shape(cell), boundary.is_boundary_cell(cell), boundary.absorbs_border(cell)
```

An ordinary interior parallelogram. Compare against one of each other family:

```{code-cell} python
[(c, boundary.cell_shape(c), itacart.is_equal_area_cell(c))
 for c in (cell, "NE(0000/0500)", "SE(1930/0196)", "NE(0000/1000)")]
```

A meridian triangle keeps the nominal area; a border trapezoid does not; and a
polar cap is *called* a triangle and does not either. This is why the guard is
always `is_equal_area_cell` and never the shape name.
[Polar & Border Cell Families](../../concepts/geometry/polar_and_border_cell_families.md).

## Its geometry

```{code-cell} python
itacart.cell_to_anchor(cell), itacart.cell_to_centroid(cell)
```

```{code-cell} python
itacart.cell_to_boundary(cell, close=True)
```

And on the sinusoidal plane, where the lattice is regular:

```{code-cell} python
itacart.cell_to_sinusoidal(cell)
```

## Its measurements

```{code-cell} python
(itacart.nominal_cell_area(9),
 itacart.effective_cell_area(cell),
 itacart.is_equal_area_cell(cell))
```

Nominal and effective agree, which is the guarantee doing its job.

```{code-cell} python
(round(metrics.cell_base_angle(cell), 6),
 round(metrics.compactness(cell), 6),
 round(metrics.normalized_cell_area(cell), 6))
```

The base angle is what the sinusoidal plane's forty-five degrees becomes on the
ellipsoid — it varies across the lattice and is measured, not predicted. The
other two are transcriptions of published conventions, for comparing this grid
against others; each states where its convention stops being trustworthy.
[Area, Perimeter & Compactness](../../concepts/geometry/area_perimeter_and_compactness.md).

## Is it where you think it is?

```{code-cell} python
from itacart import interop
interop.cell_to_wkt(cell)
```
