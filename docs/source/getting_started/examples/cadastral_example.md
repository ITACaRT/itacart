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

# Cadastral Example

Two adjoining parcels, registered properly: covered, measured, identified,
stored and exported. This is the session the whole library is for.

```{code-cell} python
import itacart
from itacart import geometry, serialization, interop, index as ix
from shapely.geometry import Polygon

lot_a = Polygon([(-46.6340, -23.5520), (-46.6330, -23.5520),
                 (-46.6330, -23.5500), (-46.6340, -23.5500)])
lot_b = Polygon([(-46.6330, -23.5520), (-46.6320, -23.5520),
                 (-46.6320, -23.5500), (-46.6330, -23.5500)])
RESOLUTION, CONTAINMENT = 7, "center"
```

The two lots share a boundary along −46.6330. That shared edge is where a
representation usually goes wrong, in one of two ways: the boundary is counted
twice, or it falls through the gap.

## Cover both

```{code-cell} python
dense_a = geometry.densify_orthodromic(lot_a, max_segment_m=5.0)
dense_b = geometry.densify_orthodromic(lot_b, max_segment_m=5.0)
cover_a = geometry.polyfill(dense_a, RESOLUTION, containment=CONTAINMENT)
cover_b = geometry.polyfill(dense_b, RESOLUTION, containment=CONTAINMENT)
ix.count_cells(cover_a), ix.count_cells(cover_b)
```

## The shared boundary is neither double-counted nor lost

```{code-cell} python
cells_a, cells_b = set(ix.decompose(cover_a)), set(ix.decompose(cover_b))
len(cells_a & cells_b), len(cells_a | cells_b)
```

No cell belongs to both. And filling the two lots together gives exactly the
union:

```{code-cell} python
together = geometry.polyfill(
    geometry.densify_orthodromic(lot_a.union(lot_b), max_segment_m=5.0),
    RESOLUTION, containment=CONTAINMENT)
ix.count_cells(together), ix.count_cells(together) == len(cells_a | cells_b)
```

The parts sum to the whole, with no overlap and no gap. Under `center` a cell
on the shared edge goes to whichever lot holds its centre — exactly one of
them — which is the property that makes a cadastral partition work.

## Measure

```{code-cell} python
def parcel_area(cover, resolution):
    """Nominal where the equal-area guarantee holds, measured where it does not."""
    return sum(itacart.nominal_cell_area(resolution) if itacart.is_equal_area_cell(c)
               else itacart.effective_cell_area(c)
               for c in ix.decompose(cover))

area_a = parcel_area(cover_a, RESOLUTION)
area_b = parcel_area(cover_b, RESOLUTION)
area_a, area_b, area_a + area_b
```

The guard matters even though these lots are nowhere near a seam: a pipeline
that assumes the nominal figure everywhere is a pipeline that will be wrong the
first time it is handed a parcel in Fiji.

Report the figure with a bracket, not alone:

```{code-cell} python
lo = ix.count_cells(geometry.polyfill(dense_a, RESOLUTION, containment="contains"))
hi = ix.count_cells(geometry.polyfill(dense_a, RESOLUTION, containment="intersects"))
(lo * itacart.nominal_cell_area(RESOLUTION),
 area_a,
 hi * itacart.nominal_cell_area(RESOLUTION))
```

## Capture identity

The cover says *where* the lot is. The vertex sequence says *which* lot it is —
and only the second survives as an identity, because a cover is lossy about the
boundary.

```{code-cell} python
ring_a = geometry.vertex_to_cell(dense_a, resolution=13)
blob_a = serialization.encode_geometry([ring_a], geometry_type="POLYGON")
digest_a = serialization.geometry_hash(blob_a)
len(ring_a), len(blob_a), digest_a[:12]
```

The hash is taken over the canonical form, so the same lot encoded again —
starting the ring at a different vertex — produces the same digest. That is
what lets two authorities compare their copies without sending geometry.
[Canonical Geometry](../../concepts/cadastral_workflows/canonical_geometry.md).

## Store

```{code-cell} python
compact_a = itacart.compact_cells(cover_a)
tree_a = serialization.serialize_to_blob(compact_a)
ix.count_cells(cover_a), ix.count_cells(compact_a), len(tree_a)
```

## Verify before writing anything down

```{code-cell} python
(serialization.validate_geometry(blob_a) is None,
 serialization.validate_tree(tree_a) is None,
 serialization.deserialize_from_blob(tree_a) == ix.normalize(compact_a))
```

Both blobs are structurally sound and the tree round-trips to the index it
stores. A blob that survives validation decodes to cells that exist.

## Export for someone who has never heard of this grid

```{code-cell} python
collection = interop.cells_to_geojson(compact_a, properties={"lot": "A"})
(collection["type"],
 len(collection["features"]),
 collection["features"][0]["properties"])
```

Each Feature carries more than your own properties: the index, the resolution,
the shape family, and both the nominal and the effective area of that cell. A
consumer that knows nothing about ITACaRT can still see which cells are not
equal-area, without having to ask.

```{code-cell} python
interop.recover_from_geojson(collection) == ix.decompose(compact_a)
```

Recovery is exact because the index travelled in each Feature's `id`. What you
must not do is rebuild the cells from the exported coordinates — that is
reingestion of a rendering, and nothing promises it returns the cell.
[Geographic vs Grid Geometry](../../concepts/data_and_interoperability/geographic_vs_grid_geometry.md).

## What the record holds

Six fields, and none of them is optional:

```{code-cell} python
record = {
    "lot": "A",
    "identity": digest_a.hex()[:16],
    "resolution": RESOLUTION,
    "containment": CONTAINMENT,
    "cells": ix.count_cells(compact_a),
    "area_m2": round(area_a, 3),
}
record
```

Without the resolution and the containment mode, the cover is a set of cells
whose meaning cannot be reconstructed and whose area cannot be audited by
anyone who was not in the room when it was computed.
