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

# End-to-End Cadastral Workflow

The whole thing, from a surveyed boundary to a stored, verifiable record. Every
step below runs when this page is built.

```{code-cell} python
import itacart
from itacart import geometry, serialization, interop, index as ix
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
```

## 1. Choose the resolution from the source

The survey sits around 1:2 000, so the analysis column points at resolution 9.
Taking 13 would add four levels of invented detail —
[Choosing Resolution from Mapping Scale](choosing_resolution_from_mapping_scale.md).

```{code-cell} python
RESOLUTION = 9
```

## 2. Densify before covering

A straight edge on the plane is not a geodesic, so a long edge must gain
intermediate vertices or the cover is wrong in the middle, silently.

```{code-cell} python
dense = geometry.densify_orthodromic(parcel, max_segment_m=5.0)
len(parcel.exterior.coords), len(dense.exterior.coords)
```

`polyfill` does this itself at a threshold derived from the resolution; doing
it explicitly makes the parameter a decision rather than a default.

## 3. Cover the parcel, and record the mode

```{code-cell} python
CONTAINMENT = "center"
cover = geometry.polyfill(dense, resolution=RESOLUTION, containment=CONTAINMENT)
itacart.count_cells(cover)
```

The mode is not recoverable from the cover, so it is stored beside it —
[Containment Modes](../spatial_operations/containment_modes.md).

## 4. Report the area, with a bracket

```{code-cell} python
def parcel_area(cells, resolution):
    return sum(itacart.nominal_cell_area(resolution) if itacart.is_equal_area_cell(c)
               else itacart.effective_cell_area(c)
               for c in cells)

area = parcel_area(ix.decompose(cover), RESOLUTION)
lo = ix.count_cells(geometry.polyfill(dense, RESOLUTION, containment="contains"))
hi = ix.count_cells(geometry.polyfill(dense, RESOLUTION, containment="intersects"))
area, (lo * itacart.nominal_cell_area(RESOLUTION),
       hi * itacart.nominal_cell_area(RESOLUTION))
```

The single figure, and the bracket the two other modes put around it.

## 5. Capture the boundary as identity

The cover says where; the vertex sequence says which. Both are needed, and the
second is the one a registry compares.

```{code-cell} python
ring = geometry.vertex_to_cell(dense, resolution=13)
blob = serialization.encode_geometry([ring], geometry_type="POLYGON")
digest = serialization.geometry_hash(blob)
len(ring), len(blob), digest[:12]
```

The boundary is captured at resolution 13 even though the cover is at 9 — the
vertices are positions, not an extent, and there is no reason to coarsen them.

## 6. Compact the cover for storage

```{code-cell} python
compact = itacart.compact_cells(cover)
tree = serialization.serialize_to_blob(compact)
itacart.count_cells(cover), itacart.count_cells(compact), len(tree)
```

## 7. Verify before storing

```{code-cell} python
(serialization.validate_geometry(blob) is None,
 serialization.validate_tree(tree) is None,
 serialization.deserialize_from_blob(tree) == ix.normalize(compact))
```

Both blobs validate, and the tree round-trips to the index it stores.

## 8. Export for a consumer that does not know the grid

```{code-cell} python
collection = interop.cells_to_geojson(compact)
(collection["type"], len(collection["features"]),
 interop.recover_from_geojson(collection) == ix.decompose(compact))
```

The index travels in each Feature's `id`, which is what makes recovery exact.
Do not plan to rebuild the cells from the exported coordinates —
[Geographic vs Grid Geometry](../data_and_interoperability/geographic_vs_grid_geometry.md).

## What the record holds

| Field | Value from above | Why |
|---|---|---|
| Geometry hash | `geometry_hash(blob)` | the parcel's identity |
| GeometryBlob | `blob` | the boundary, canonical |
| TreeBlob | `tree` | the cover, compacted |
| Resolution | 9 | what the cover means |
| Containment | `center` | how the cover was decided |
| Area | computed above | with its bracket |

The resolution and the containment mode are important. Without them the
cover is a set of cells whose meaning cannot be reconstructed, and the area
cannot be audited by anyone who was not in the room.
