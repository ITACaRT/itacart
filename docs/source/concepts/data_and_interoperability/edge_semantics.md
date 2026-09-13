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

# Edge Semantics

Two positions do not determine a line on a curved surface. What runs between
them is a choice, and an encoded geometry declares which choice it made.

```{code-cell} python
from itacart import serialization, geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
ring = geometry.vertex_to_cell(parcel, resolution=13)
blob = serialization.encode_geometry([ring], geometry_type="POLYGON")
profile = serialization.decode_geometry(blob)
{k: profile[k] for k in ("edge_model", "densification_model", "max_segment_m")}
```

## The two edge models

**`WGS84_GEODESIC`** — an edge is the geodesic between its endpoints on the
ellipsoid. This is what a surveyor means by a straight line between two marks,
and it is the default.

**`PLANAR_SINUSOIDAL_STRAIGHT`** — an edge is a straight segment on the
sinusoidal plane. This is what the grid's own cells are: a cell's boundary is
straight in the projection, which is why the exported vertices rebuild the cell
to within a unit in the last place.

They are different lines between the same endpoints, and over a long edge they
are noticeably different:

```{figure} ../../../_static/f7/f7_05_chord_versus_geodesic.png
:alt: A chord against the geodesic between the same endpoints
:width: 100%

`docs/_static/f7/f7_05_chord_versus_geodesic.png`
```

## Densification is a separate field

An edge model says what the line *is*. Densification says how many intermediate
vertices were inserted so that a consumer drawing straight segments gets close
to it: `NONE`, `ORTHODROMIC_VINCENTY` with a maximum segment length, or
`ASSUMED_BY_PRODUCER` when the producer asserts the geometry is already dense
enough.

## Two combinations are contradictions

Declaring straight planar edges *and* geodesic densification says two
incompatible things about the same geometry. The format refuses rather than
storing the contradiction:

```{code-cell} python
from itacart.exceptions import IncompatibleProfileError
try:
    serialization.encode_geometry(
        [ring], geometry_type="POLYGON",
        edge_model="PLANAR_SINUSOIDAL_STRAIGHT",
        densification_model="ORTHODROMIC_VINCENTY")
except IncompatibleProfileError as exc:
    print(exc)
```

The full compatibility matrix is in the
[GeometryBlob binary encoding](geometry_encoding_spec.md) specification.

## Why it is declared rather than inferred

A consumer cannot tell from coordinates alone which model produced them. Dense
vertices along a geodesic and dense vertices along a planar chord look the
same — a list of positions — and the difference only appears when someone
measures an area or tests a containment.

Carrying the model in the header means the question has an answer that does not
depend on guessing. It is the same reasoning as carrying the containment mode
with a cover: [Containment Modes](../spatial_operations/containment_modes.md).

## Which to use

Use `WGS84_GEODESIC` for geometry that came from the world — surveyed
boundaries, GNSS traces, imported cadastral records. Use
`PLANAR_SINUSOIDAL_STRAIGHT` for geometry that came from the grid, because that
is what grid cells actually are. Do not convert between them by relabelling the
header: the vertices mean different shapes under the two models.
