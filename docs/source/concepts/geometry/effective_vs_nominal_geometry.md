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

# Effective vs Nominal Geometry

Two numbers describe a cell's area. The **nominal** area is what the
specification assigns to the resolution. The **effective** area is what the
cell's own boundary encloses. They agree almost everywhere, and the places they
do not are the reason both exist.

```{code-cell} python
import itacart
cell = "SE(1400/0374)"
itacart.nominal_cell_area(1), itacart.effective_cell_area(cell)
```

Identical, for an interior cell. Note the signatures differ in a way that says
what each one is: nominal takes a *resolution*, because under the guarantee the
cell cannot matter; effective takes a *cell*, because it is a measurement of
that cell.

## Where they part

```{code-cell} python
absorber = "SE(1930/0196)"
itacart.effective_cell_area(absorber), itacart.nominal_cell_area(1)
```

Roughly half. This cell absorbs the domain border: its longer base is clipped
at the boundary line, so the area is computed from the actual vertices rather
than assumed from the table.

```{figure} ../../../_static/f4/f4_06_effective_vs_nominal_area.png
:alt: Effective area against nominal area across the border
:width: 100%

`docs/_static/f4/f4_06_effective_vs_nominal_area.png`
```

For clipped and other special boundary cells, the effective area is the
area of the geometry that actually exists, not the nominal area of a
regular cell at the same resolution.

This distinction matters for geometric correctness and for any computation
that relies on the actual cell footprint. Returning the nominal area would
overstate the area of these cells, in some cases substantially.

In ITACaRT, however, the affected cells occur in polar or oceanic boundary
regions rather than in the cadastral areas for which the system is intended.
The issue is therefore one of geometric and API correctness, not a practical parcel-area or taxation error.

## Asking the construction, not the shape name

The natural guard looks like `shape != "trapezoid"`. It is wrong.

```{code-cell} python
from itacart import boundary
[(c, boundary.cell_shape(c), boundary.absorbs_border(c), itacart.is_equal_area_cell(c))
 for c in ("SE(1400/0374)", "NE(0000/0500)", "SE(1930/0196)")]
```

A prime-meridian triangle is not a parallelogram and still carries the nominal
area — a triangle has the same base-and-height behaviour, so the area survives
the change of shape. Meanwhile a polar cell is still *named* a triangle and
still absorbs the border, so it does not carry the nominal area.

Shape and area-honesty are two different questions. `is_equal_area_cell` asks
the second one directly, and it is the guard to use:

```{code-cell} python
def parcel_area(cells, resolution):
    """Nominal where the guarantee holds, measured where it does not."""
    return sum(itacart.nominal_cell_area(resolution) if itacart.is_equal_area_cell(c)
               else itacart.effective_cell_area(c)
               for c in cells)

parcel_area(["SE(1400/0374)", "SE(1930/0196)"], 1)
```

## Why the package refuses to choose for you

It would be easy for `nominal_cell_area` to quietly return the effective
area when the two diverge. ITACaRT deliberately keeps the quantities separate:
the nominal area describes the regular cell implied by the resolution, while
the effective area describes the geometry that actually exists.

Silently substituting one for the other would hide an important distinction
in the grid model and could make downstream geometric calculations appear
consistent while using the wrong quantity.

For that reason the two functions remain separate, and `is_equal_area_cell`
lets the caller determine whether the nominal and effective interpretations
coincide for a particular cell. In ordinary cadastral regions they normally
do; the distinction becomes relevant mainly in the special polar and boundary
families where the effective geometry departs from the regular cell.
