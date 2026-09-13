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

# Parcel Representation

A parcel can be represented three ways in this grid, and they answer different
questions. Choosing badly is not a performance problem — it is a correctness
problem, because two of the three are lossy in different directions.

```{code-cell} python
import itacart
from itacart import geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6340, -23.5520), (-46.6320, -23.5520),
                  (-46.6320, -23.5500), (-46.6340, -23.5500)])
```

## As a cover — the cells the parcel occupies

```{code-cell} python
cover = geometry.polyfill(parcel, resolution=7)
itacart.count_cells(cover)
```

This is the representation for area, overlap and containment. It says *where*
the parcel is, to within one cell, and it is what makes area a count.

It is lossy about the boundary. Every cell in the cover is entirely in or
entirely out, so the true outline is approximated to the cell size — and which
approximation you get depends on the containment mode, which must be recorded
alongside the cover.

```{figure} ../../../_static/fpaper/fpaper_07_cells_over_an_osm_base.png
:alt: ITACaRT cells at resolutions 6 and 7 over an OpenStreetMap base
:width: 80%

Figure 7 of the paper: a real parcel covered at two resolutions — coarse
cells inside, fine cells along the boundary, which is what compaction
produces.
`docs/_static/fpaper/fpaper_07_cells_over_an_osm_base.png`
```

## As a compacted cover — the same ground, fewer cells

```{code-cell} python
compact = itacart.compact_cells(cover)
itacart.count_cells(cover), itacart.count_cells(compact)
```

Identical ground, far fewer cells. The interior folds into coarse cells and the
fine ones stay along the boundary, which is where the information actually is.

Compaction is lexical, so it cannot introduce a geometric error —
[Compaction & Normalization](../index_and_hierarchy/compaction_and_normalization.md).
Use it before storing, hashing or transmitting. Note that an area computed from
a compacted cover must weight each cell by its own resolution, because the
cells are no longer all the same size.

## As a vertex sequence — the boundary itself

```{code-cell} python
ring = geometry.vertex_to_cell(parcel, resolution=13)
len(ring), ring[0]
```

Each vertex of the parcel mapped to the cell containing it, **in sequence**.
This is the representation for identity: it preserves the outline, the ring
topology and the order, so two records of the same parcel can be compared as
the same object.

It says nothing about area. A vertex list is a boundary, not a region.

```{code-cell} python
rebuilt = geometry.cells_to_geometry(ring, geometry_type="Polygon")
rebuilt.geom_type, rebuilt.is_valid
```

`cells_to_geometry` is the inverse, which is what makes the sequence a
representation rather than a summary.

## Which to store

Most cadastral records need **both**, and they are not substitutes.

Store the **vertex sequence** — encoded as a GeometryBlob — as the parcel's
identity. It is what a hash is taken over, what two registries compare, what
survives an export and comes back the same.
[GeometryBlob](../data_and_interoperability/geometry_blob.md).

Store or derive the **compacted cover** for anything areal or spatial. It is
what a tax base is computed from and what a spatial index is built on.

The bridge between them is one-way: a GeometryBlob yields the TreeBlob of its
vertex set, and order, ring topology and type vanish in the derivation.
Coverage survives; identity does not. So the cover can always be rebuilt from
the boundary, and the boundary can never be rebuilt from the cover.

## The mistake to avoid

Do not store only the cell cover and treat it as the parcel geometry.

A cover is a discrete representation of the parcel at a chosen resolution.
It preserves which cells were selected and can be used for indexing, filling,
comparison, and area estimation, but it does not preserve the original
boundary inside the cells crossed by that boundary.

Once the vector outline is discarded, that information cannot be
reconstructed from the cover. A later conversion back to GeoJSON can draw
the boundary of the selected cells, but that is the boundary of the cover,
not the cadastral boundary from which the cover was derived.

Resolution controls how finely the parcel is represented, but increasing
resolution does not turn a cell cover into the original survey geometry.
Keep the authoritative parcel geometry and use the ITACaRT cover as its
discrete spatial representation.
