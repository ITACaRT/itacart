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

# Fill a Parcel

Cover a polygon with cells, under each definition of "inside", and see what the
choice costs.

```{code-cell} python
import itacart
from itacart import geometry, index as ix
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6330, -23.5520),
                  (-46.6330, -23.5500), (-46.6340, -23.5500)])
parcel.bounds
```

About a hundred metres by two hundred.

## Densify first

A straight edge on the plane is not a geodesic, so a long edge that has not
gained intermediate vertices fills the wrong cells in the middle — and does it
without raising anything.

```{code-cell} python
dense = geometry.densify_orthodromic(parcel, max_segment_m=5.0)
len(parcel.exterior.coords), len(dense.exterior.coords)
```

`polyfill` does this itself at a threshold it derives from the resolution.
Doing it explicitly makes the threshold a decision instead of a default, and
the step is idempotent:

```{code-cell} python
len(geometry.densify_orthodromic(dense, max_segment_m=5.0).exterior.coords)
```

## Three covers, three answers

```{code-cell} python
covers = {mode: geometry.polyfill(dense, resolution=7, containment=mode)
          for mode in ("contains", "center", "intersects")}
{mode: ix.count_cells(c) for mode, c in covers.items()}
```

They nest — `contains` inside `center` inside `intersects` — and the nesting is
a property of how acceptance is decided, not a coincidence of this parcel:

```{code-cell} python
sets = {m: set(ix.decompose(c)) for m, c in covers.items()}
sets["contains"] <= sets["center"] <= sets["intersects"]
```

Between the tightest and the loosest is more than a quarter of the cells. A
stored cover that does not record which mode produced it cannot be compared
with another one.
[Containment Modes](../../concepts/spatial_operations/containment_modes.md).

## Counting without naming

If the answer you want is an area, naming the cells is wasted work:

```{code-cell} python
n = geometry.count_internal_cells(dense, 7)
n, n * itacart.nominal_cell_area(7)
```

Square metres directly — no conversion factor, because the cells are
equal-area and the areas are decimal.

## Refining, and knowing when to stop

```{code-cell} python
[(res,
  geometry.count_internal_cells(dense, res),
  geometry.count_internal_cells(dense, res) * itacart.nominal_cell_area(res))
 for res in (5, 7, 9, 11)]
```

The estimate settles. The cell count does not — it is multiplied by a hundred
every two levels — so past the resolution your source mapping supports you are
buying storage rather than accuracy.
[Choosing Resolution from Mapping Scale](../../concepts/cadastral_workflows/choosing_resolution_from_mapping_scale.md).

## Compact before storing

```{code-cell} python
cover = covers["center"]
compact = itacart.compact_cells(cover)
ix.count_cells(cover), ix.count_cells(compact)
```

Identical ground, a fraction of the cells: the solid interior folds into coarse
cells and the fine ones stay where the information is, along the boundary.

Compaction is lexical — it manipulates the index tree and never consults
geometry — so it cannot introduce a geometric error.

## When a parcel is smaller than a cell

```{code-cell} python
from itacart.exceptions import ITACaRTError
tiny = Polygon([(-46.63400, -23.55200), (-46.63399, -23.55200),
                (-46.63399, -23.55199), (-46.63400, -23.55199)])
try:
    geometry.polyfill(tiny, resolution=5)
except ITACaRTError as exc:
    print(type(exc).__name__, "--", exc)
```

The package raises rather than returning an empty cover. An empty cover read as
"no overlap" is the kind of error that reaches a land registry.
