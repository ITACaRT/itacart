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

# Antimeridian

The grid does not wrap around the world. Each quadrant's lattice runs outward
from the prime meridian and stops, and where it stops is the antemeridian at
180 degrees — except in two named zones, which are
[Extension Zones](extension_zones.md).

## Why it stops

Maintaining uniform cell areas all the way around is impractical: the lattice
that fits one hemisphere does not close cleanly against itself on the other
side. Rather than distort every cell slightly to force a closure, ITACaRT
accepts a discontinuity at one meridian and keeps every cell elsewhere exact.

That trade is defensible because of what lies on that meridian. Almost all of
it is ocean. The land it crosses is Fiji, a piece of the Russian mainland in
the Chukotka Autonomous Okrug with Wrangel Island nearby, and a stretch of
Antarctica — and Antarctica is excluded on the grounds that it has little
cadastral use. So the cells that pay for the discontinuity are in the ocean and
in Antarctica, and the inhabited land is handled by extending two zones.

## What it means for geometry

A geometry that spans the 180th meridian outside a defined extension zone
cannot be covered, and the package says so rather than producing a cover that
wraps:

```{code-cell} python
from itacart import boundary
from shapely.geometry import LineString
spanning = LineString([(179.0, 10.0), (-179.0, 10.0)])
boundary.crosses_antemeridian(spanning)
```

`crosses_antemeridian` is the test to run before filling anything whose extent
you do not control. Inside an extension zone the situation is different, and
the zone's own page covers it.

## Longitudes past 180

Inside an extension zone, cells carry longitudes that run **past** 180 degrees
rather than wrapping back to negative values:

```{code-cell} python
import itacart
cell = itacart.geo_to_cell(-179.0, -18.0, resolution=1)
poly = itacart.cell_to_polygon(cell)
poly.is_valid, round(poly.bounds[0], 3), round(poly.bounds[2], 3)
```

This is deliberate and it matters. A polygon whose coordinates wrapped would be
self-intersecting, its area would be meaningless, and every areal operation on
it would silently return nonsense. Carrying the longitude past the line keeps
the polygon simple and its area measurable.

The cost is that such a polygon is not directly valid GeoJSON, which requires
longitudes in the range from −180 to 180.
[GeoJSON](../data_and_interoperability/geojson.md) covers how export handles
it.

## Indexing a position near the line

A position inside the latitude band of an extension zone has two spellings —
the eastern one past 180, and the ordinary negative one. The grammar admits
only the eastern spelling, and the package's refusal message says so:

```{code-cell} python
from itacart.exceptions import ITACaRTError
try:
    itacart.geo_to_cell(181.0, 0.0, resolution=9)
except ITACaRTError as exc:
    print(exc)
```

At latitude 0 there is no extension zone, so the coordinate is simply out of
domain. Inside a zone's band the same longitude would be accepted.

## The rule to carry away

Test with `crosses_antemeridian` before covering geometry of unknown extent.
Expect longitudes past 180 from cells inside an extension zone, and do not
normalise them before measuring.
