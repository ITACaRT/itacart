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

# Choosing Resolution from Mapping Scale

The resolution is the single most consequential parameter in a cadastral
pipeline, and it should be chosen from the scale of the source mapping rather
than from how precise you would like the answer to be.

```{code-cell} python
import itacart
table = itacart.resolution_table()
[(r["resolution"], r["cell_size_m"], r["visualization_scale"], r["analysis_scale"])
 for r in table if r["resolution"] in (3, 5, 7, 9, 11, 13)]
```

## Two scales, not one

Each resolution carries two suggested scales, and they are different numbers
for a reason.

The **visualisation scale** comes from the minimum visible line on a
cartographic representation: below it, the cells are too small to see. The
**analysis scale** comes from sampling theory: below it, the cells are finer
than the data can support, and the extra resolution is invented detail.

The analysis scale is consistently finer than the visualisation scale. A cell
you can see is coarser than a cell you should compute with, which is the
opposite of the intuition most people bring.

## The range cadastral work spans

A cadastral hierarchy has to cover roughly 1:500 for urban plots out to
1:10 000 for large rural regions, and at the fine end it has to reach the
centimetre accuracy a modern GNSS field survey delivers. That is the span the
fourteen levels were designed against.

Reading the table the other way — from a scale to a resolution — a 1:1 000
urban survey sits around resolution 11, and a 1:10 000 rural sheet around
resolution 9.

## Do not over-refine

The temptation is to take resolution 13 always, on the grounds that finer
cannot hurt. It can.

```{code-cell} python
from itacart import geometry
from shapely.geometry import Polygon
parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
[(r, geometry.count_internal_cells(parcel, r)) for r in (7, 9, 11)]
```

Two levels of refinement multiply the cell count by a hundred. Four levels
multiply it by ten thousand. Storage, hashing and transmission all scale with
it, and none of that extra detail is knowable if the parcel boundary came from
a 1:2 000 sheet.

The honest resolution is the one the source supports. Refining past it produces
a number with more digits and no more truth.

## The ceiling and the floor

```{code-cell} python
from itacart.exceptions import ITACaRTError
try:
    itacart.geo_to_cell(-46.6329, -23.5509, resolution=14)
except ITACaRTError as exc:
    print(exc)
```

Thirteen is the finest level, and one centimetre is below what any
GNSS-referenced cadastral survey delivers — so the ceiling is where the grid
stops promising more than the measurement can support rather than a limitation
in practice.

At the coarse end, a parcel smaller than a cell covers nothing under the
default containment mode, and the package raises rather than returning an empty
cover. [Containment Modes](../spatial_operations/containment_modes.md) is what
to do about that.
