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

# Index a Point

Turn a geodetic coordinate into a cell index, then get the geography back.

```{code-cell} python
import itacart
itacart.__version__
```

## One coordinate, one cell

A position in São Paulo, addressed at one-centimetre resolution:

```{code-cell} python
LON, LAT = -46.6328862, -23.5508962
cell = itacart.geo_to_cell(LON, LAT, resolution=13)
cell
```

Longitude first, then latitude — the order GeoJSON uses, and the opposite of
the order most people say out loud.

The string is not opaque. `SW` is the quadrant, `0476/0260` is the
ten-kilometre cell in that quadrant's Cartesian grid, and each parenthesis
after it is one step finer.
[Index Structure](../../concepts/index_and_hierarchy/index_structure.md) reads
it properly.

## The same point at every resolution

```{code-cell} python
for res in (1, 3, 5, 7, 9, 11, 13):
    print(res, "\t", itacart.geo_to_cell(LON, LAT, resolution=res))
```

Each line is a prefix of the next. That is the hierarchy — a coarser address is
literally the start of the finer one, so no lookup is needed to go up a level.

## Back to geography

```{code-cell} python
lon, lat = itacart.cell_to_centroid(cell)
lon, lat
```

The centroid is the cell's representative position — the one point every
single-point route answers with. The round trip closes:

```{code-cell} python
itacart.geo_to_cell(lon, lat, resolution=13) == cell
```

## The anchor is a different point

```{code-cell} python
itacart.cell_to_anchor(cell)
```

Close to the centroid at this resolution, and not the same thing. The anchor is
the lattice corner the index literally encodes; it sits on the cell's boundary,
which is why it is not used as the representative position.
[Cells, Vertices, and Spatial Identity](../../concepts/grid_fundamentals/cells_vertices_and_spatial_identity.md).

## The whole cell

```{code-cell} python
ring = itacart.cell_to_boundary(cell, close=True)
ring
```

Four vertices and a repeat of the first. As a Shapely polygon:

```{code-cell} python
poly = itacart.cell_to_polygon(cell)
poly.geom_type, poly.is_valid
```

## How big is it, really

```{code-cell} python
itacart.nominal_cell_area(13), itacart.is_equal_area_cell(cell)
```

One square centimetre, and the equal-area guarantee holds for this cell. That
second call is not ceremony — near the antemeridian and the poles it returns
`False`, and then the nominal figure is the wrong number.
[Equal Area](../../concepts/grid_fundamentals/equal_area.md).

## What the grid refuses

```{code-cell} python
from itacart.exceptions import ITACaRTError
for lon_, lat_, res in [(181.0, 0.0, 9), (LON, LAT, 14)]:
    try:
        itacart.geo_to_cell(lon_, lat_, res)
    except ITACaRTError as exc:
        print(type(exc).__name__, "--", exc)
```

Both refusals are deliberate. A position the grid cannot address is almost
always a fault upstream, and answering it with the nearest plausible cell
buries that fault in stored data.
