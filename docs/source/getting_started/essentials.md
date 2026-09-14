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

# Essentials

Eight things you meet in the first hour, one screen each. Every output below
is produced when this page is built.

This is a map, not an argument. Each section links to the page where the
reasoning lives.

```{code-cell} python
import itacart
from itacart import geometry, interop, topology, index as ix
from shapely.geometry import Polygon
```

## 1. A point becomes a cell

```{code-cell} python
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=9)
cell
```

Longitude, latitude, resolution. Resolution 9 is one metre; the ladder runs 1
to 13, and 0 addresses the quadrant.
→ [Point-to-Cell Indexing](../concepts/spatial_operations/point_to_cell_indexing.md)

## 2. A cell becomes geometry

```{code-cell} python
(itacart.cell_to_centroid(cell),
 len(itacart.cell_to_boundary(cell, close=True)),
 itacart.cell_to_polygon(cell).geom_type)
```

The **centroid** is the cell's representative position and the point every
single-point route answers with. `cell_to_anchor` gives the lattice corner the
index encodes — a different point, on the boundary, not a substitute.
→ [Cell Geometry](../concepts/geometry/cell_geometry.md)

## 3. A polygon becomes a cover

```{code-cell} python
parcel = Polygon([(-46.6340, -23.5520), (-46.6330, -23.5520),
                  (-46.6330, -23.5500), (-46.6340, -23.5500)])
{mode: ix.count_cells(geometry.polyfill(parcel, 7, containment=mode))
 for mode in ("contains", "center", "intersects")}
```

Three definitions of *inside*, three covers, nested in that order. `center` is
the default and the unbiased one for area. **Store the mode beside the cover**
— it is not recoverable from the cells.
→ [Containment Modes](../concepts/spatial_operations/containment_modes.md)

## 4. A cover becomes GeoJSON

```{code-cell} python
cover = geometry.polyfill(parcel, 7)
collection = interop.cells_to_geojson(cover)
(len(collection["features"]),
 interop.recover_from_geojson(collection) == ix.decompose(cover))
```

The index travels in each Feature's `id`, which is what makes recovery exact.
Rebuilding cells from the exported **coordinates** is a different operation and
is not promised to return them.
→ [Geographic vs Grid Geometry](../concepts/data_and_interoperability/geographic_vs_grid_geometry.md)

## 5. Moving up and down

```{code-cell} python
parent = itacart.get_parent(cell)
(parent,
 len(list(itacart.get_children(parent, flatten=True))),
 ix.count_cells(itacart.compact_cells(cover)))
```

Ascent is truncation of the string. Descent is enumeration, not
multiplication — at the domain border a cell can hold fewer children than the
ratio predicts. Compaction folds a solid interior into coarse cells without
touching geometry.
→ [Parent, Children & Descendants](../concepts/index_and_hierarchy/parent_children_and_descendants.md)

## 6. What comes back: one cell or many

An index names one cell or a whole set, and the return shape follows. This is
the single largest source of surprise in the library.

The table below is built when this page is compiled, by calling each function
twice — once with a single cell, once with a composite index of three — and
reporting what came back.

```{code-cell} python
:tags: [remove-input]

import types

one = cell
three = ix.compose([cell, topology.get_neighbor(cell, "N"),
                    topology.get_neighbor(cell, "E")])

def shape(value):
    if isinstance(value, bool): return "bool"
    if isinstance(value, str): return "str"
    if isinstance(value, int): return "int"
    if isinstance(value, tuple): return f"tuple[{len(value)}]"
    if isinstance(value, list):
        if not value: return "list[]"
        kinds = {type(item).__name__ for item in value}
        if kinds == {"list"}: return f"list[{len(value)}] of list"
        return f"list[{len(value)}] of {kinds.pop()}"
    return type(value).__name__

def call(fn, arg, **kw):
    out = fn(arg, **kw)
    return shape(list(out) if isinstance(out, types.GeneratorType) else out)

calls = [("get_parent", itacart.get_parent, {}),
         ("get_ancestors", itacart.get_ancestors, {}),
         ("get_children", itacart.get_children, {}),
         ("get_children(flatten=True)", itacart.get_children, {"flatten": True}),
         ("grid_disk(k_distance=1)", topology.grid_disk, {"k_distance": 1}),
         ("cell_to_centroid", itacart.cell_to_centroid, {}),
         ("cell_to_boundary", itacart.cell_to_boundary, {}),
         ("is_equal_area_cell", itacart.is_equal_area_cell, {}),
         ("count_cells", ix.count_cells, {}),
         ("compact_cells", itacart.compact_cells, {})]

print(f"{'':<28}{'one cell':<22}one composite index of three")
for name, fn, kw in calls:
    print(f"{name:<28}{call(fn, one, **kw):<22}{call(fn, three, **kw)}")
```

Two rules and one exception.

A function answering **one value per cell** gives the bare value for a single
cell and a list aligned with `decompose()` order for a composite index. A
function answering **many values per cell** gives the flat list for a single
cell and a list of lists for a composite one.

`get_children` is the exception: it groups **always**, so
`list(get_children(cell))` has length one and the children are inside it. Pass
`flatten=True` to read them directly.
→ [Composition & Decomposition](../concepts/index_and_hierarchy/composition_and_decomposition.md)

## 7. Coordinate order

```{code-cell} python
(itacart.cell_to_centroid(cell),        # (longitude, latitude)
 itacart.cell_to_sinusoidal(cell))      # (x, y) in metres on the plane
```

**Longitude first, then latitude**, everywhere geodetic — the GeoJSON order,
and the opposite of how most people say it aloud. Shapely geometries you hand
in follow the same convention. `cell_to_sinusoidal` is the other space
entirely: metres on the projection plane, where the lattice is regular.
→ [The Sinusoidal Construction](../concepts/grid_fundamentals/the_sinusoidal_construction.md)

## 8. Limits, and where the area is safe

```{code-cell} python
from itacart.exceptions import ITACaRTError
for lon, lat, res in [(181.0, 0.0, 9), (-46.6, -23.5, 14)]:
    try:
        itacart.geo_to_cell(lon, lat, res)
    except ITACaRTError as exc:
        print(type(exc).__name__, "--", exc)
```

Resolutions run 1 to 13. Longitudes past 180 are refused unless they are the
eastern spelling of an extension zone. The package raises rather than returning
a plausible value, because a plausible wrong answer is the one that reaches
stored data.

Area is a count times a constant — **except** at the antemeridian and the
poles, where cells do not carry the nominal area. The guard is never the shape
name:

```{code-cell} python
[(c, itacart.is_equal_area_cell(c))
 for c in (cell, "NE(0000/0500)", "SE(1930/0196)", "NE(0000/1000)")]
```

An interior cell and a meridian triangle carry it; a border trapezoid and a
polar cap do not — and the polar cap is *called* a triangle. Ask
`is_equal_area_cell`, and use `effective_cell_area` when it answers `False`.
→ [Known Boundary Exceptions](../concepts/standards_and_conformance/known_boundary_exceptions.md)

## Where to go next

[Example Notebooks](examples/index.md) work each of these through at length.
[Concepts & Guides](../concepts/index.md) argues them.
