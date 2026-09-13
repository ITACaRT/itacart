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

# GeoJSON

Exporting a cell set to GeoJSON, and reading it back.

```{code-cell} python
import itacart
from itacart import interop, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
cover = geometry.polyfill(parcel, resolution=7)
collection = interop.cells_to_geojson(cover)
collection["type"], len(collection["features"])
```

One Feature per cell. The cell's vertices become the Feature's geometry, and
the index goes into the Feature's `id` member:

```{code-cell} python
collection["features"][0]["id"]
```

## The index travels with the Feature

That `id` is what makes recovery possible:

```{code-cell} python
interop.recover_from_geojson(collection) == itacart.decompose(cover)
```

Exactly the cells that went in, in `decompose` order. The recovery works
because the index travelled in the Feature, not because it was reconstructed
from the coordinates — and that distinction is the whole subject of
[Geographic vs Grid Geometry](geographic_vs_grid_geometry.md). Read it before
relying on any round trip.

## Going the other way

`from_geojson` fills the geometries of a GeoJSON object into a cover. It is
covering, not recovery: it takes the shapes and rasterises them, and the
containment mode is a parameter with the same three choices as
[polyfill](../spatial_operations/cell_filling_polyfill.md).

The two directions are not symmetric and are not meant to be. `cells_to_geojson`
and `recover_from_geojson` are the inverse pair; `from_geojson` is a filling
operation that happens to accept the same file format.

## WKT and GeoDataFrame

```{code-cell} python
interop.cell_to_wkt(itacart.decompose(cover)[0])[:60] + " ..."
```

`cells_to_wkt` does a whole index at once, with an optional `dissolve` that
merges the cells into one geometry instead of a collection. With the `geo`
extra installed, `to_geodataframe` and `from_geodataframe` are the geopandas
side of the same pair.

## Longitudes past 180

A cell inside an extension zone carries longitudes past 180 degrees rather than
wrapping, because a wrapped polygon would be self-intersecting and its area
meaningless — [Antimeridian](../boundaries/antimeridian.md) explains why.

Strict GeoJSON wants longitudes within −180 to 180, so a consumer that
validates against the specification will object. This is a genuine
incompatibility between a format designed for a wrapped world and a grid that
refuses to wrap, and it affects only the two extension zones. If you are
exporting from Fiji or Chukotka to a strict consumer, you will need to decide
what to do about it — and measuring area after normalising the longitudes is
not one of the options.
