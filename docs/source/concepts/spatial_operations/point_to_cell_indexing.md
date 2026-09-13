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

# Point-to-Cell Indexing

One function, three arguments, one string back.

```{code-cell} python
import itacart
itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
```

Longitude first, then latitude, then the resolution. Both angles in decimal
degrees, longitude in the closed range from −180 to 180.

## What happens inside

The position is projected onto the ellipsoidal parallels plane, mirrored into
its quadrant's positive octant, sheared onto the square lattice, and then
descended: the resolution-1 Cartesian pair first, then alternating one-to-four
and one-to-twenty-five refinements down to the level you asked for.

The descent stays in plane metres the whole way. No trigonometry runs below
resolution 1 — the projection is evaluated once, at the top, and everything
after it is integer and metric arithmetic on a regular lattice.

```{figure} ../../../_static/f6/f6_01_lattice_step.png
:alt: The step from geography onto the square lattice
:width: 100%

`docs/_static/f6/f6_01_lattice_step.png`
```

That structure is why indexing is cheap and why it is exact. The expensive,
approximate part happens once; the rest is arithmetic that cannot drift.

## Consistency with the inverse

Indexing and inverse geometry are a pair, and the round trip lands inside the
cell you started in:

```{code-cell} python
cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
lon, lat = itacart.cell_to_centroid(cell)
itacart.geo_to_cell(lon, lat, resolution=13) == cell
```

The centroid re-indexes to its own cell. The anchor would not be a fair test —
it lies on the boundary, which is the subject of the next section.

## Edges, corners and the extension zones

A position exactly on a shared edge belongs to one cell, not to two, and a
position on a shared corner belongs to one, not four. The grid resolves this
rather than leaving it to floating point: the descent assigns a boundary point
deterministically, so the same coordinates always produce the same index.

Inside the latitude band of an extension zone a position has two spellings, and
the eastern one is the one the grammar admits.
[Extension Zones](../boundaries/extension_zones.md) covers why those zones
exist and where they are.

## What it refuses

Coordinates outside the addressable domain raise rather than clamping:

```{code-cell} python
from itacart.exceptions import ITACaRTError
for lon, lat, res in [(181.0, 0.0, 9), (-46.6, -23.5, 14)]:
    try:
        itacart.geo_to_cell(lon, lat, res)
    except ITACaRTError as exc:
        print(type(exc).__name__, "--", exc)
```

The first message is worth reading twice: a longitude past 180 is refused
*unless* it is the eastern spelling of an extension zone, which is the one case
where coordinates legitimately run past the antemeridian.

The refusal is the same principle as everywhere else in the package: a position
the grid cannot address is almost always a fault upstream, and answering it
with the nearest plausible cell buries that fault in stored data. Note that the
poles themselves are addressable — the lattice thins towards them rather than
stopping, which is [Polar Caps](../boundaries/polar_caps.md).

## Vertices of a whole geometry

For a line or a polygon, `vertex_to_cell` maps the vertices in sequence rather
than covering the figure:

```{code-cell} python
from shapely.geometry import LineString
from itacart import geometry
line = LineString([(-46.6329, -23.5509), (-46.6327, -23.5507)])
geometry.vertex_to_cell(line, resolution=9)
```

Sequence is preserved, which is what lets `cells_to_geometry` rebuild the
figure afterwards. Covering the figure rather than its vertices is
[Cell Filling / Polyfill](cell_filling_polyfill.md), and the two are not
interchangeable.
