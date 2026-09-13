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

# Extension Zones

Two regions where the eastern quadrants are carried past the antemeridian, so
that inhabited land crossed by the 180th meridian falls inside the grid rather
than on its seam.

```{code-cell} python
from itacart import constants
list(constants.EXTENSION_ZONES)
```

## Where they are

```{code-cell} python
from itacart import boundary
{z: boundary.extension_bounds(z) for z in ("FIJI", "CHUKOTKA")}
```

Bounds are west, south, east, north in degrees.

**Fiji** extends the southeastern quadrant to 178 degrees west, between 15.5
and 21.5 degrees south.

**Chukotka** extends the northeastern quadrant to 169.5 degrees west, between
64 and 72 degrees north, covering the Russian mainland, Wrangel Island and the
islands nearby.

```{figure} ../../../_static/f4/f4_05_extension_zones.png
:alt: The limits of the two extension zones
:width: 100%

`docs/_static/f4/f4_05_extension_zones.png`
```

```{figure} ../../../_static/fpaper/fpaper_05_extension_area_limits.png
:alt: The limits of the eastern quadrants' extension areas
:width: 100%

Figure 5 of the paper: the two extension areas and their limits.
`docs/_static/fpaper/fpaper_05_extension_area_limits.png`
```

## Why those limits

The longitude limits are given to half a degree, and the precision is the
point: half a degree is the coarsest, easiest-to-use value that does not cut
across any land other than Antarctica. A tighter figure would buy nothing and
would be harder to state; a looser one would clip a coastline.

Antarctica is excluded from the whole arrangement. The paper's reasoning is
that its cadastral applications are limited, and the consequence is explicit —
some cells in the oceans and in Antarctica have unequal areas.

```{figure} ../../../_static/f4/f4_11_extension_zone_latitude_edges.png
:alt: The latitude edges of an extension zone
:width: 100%

`docs/_static/f4/f4_11_extension_zone_latitude_edges.png`
```

```{figure} ../../../_static/f4/f4_12_extension_zone_meridian_edge.png
:alt: The meridian edge of an extension zone
:width: 100%

`docs/_static/f4/f4_12_extension_zone_meridian_edge.png`
```

## Asking whether you are in one

By position, before you have a cell:

```{code-cell} python
[(p, boundary.extension_zone_for_point(*p))
 for p in ((-179.0, -18.0), (-179.0, 0.0), (-175.0, 68.0))]
```

By cell, once you have one:

```{code-cell} python
import itacart
cell = itacart.geo_to_cell(-179.0, -18.0, resolution=1)
cell, boundary.extension_zone(cell)
```

## Inside a zone, cells are ordinary

A cell well inside an extension zone is a parallelogram and carries the nominal
area — the zone moves where the lattice stops, it does not change the cells
within it:

```{code-cell} python
(boundary.cell_shape(cell), boundary.absorbs_border(cell),
 itacart.is_equal_area_cell(cell))
```

The cells that pay are the ones at the zone's own meridian edge. There, the
same absorbing rule applies as at any other domain border:
[Absorbing Boundary Cells](absorbing_boundary_cells.md).

## What the zones are not

They are **not** a general antemeridian crossing. Outside their latitude bands
the 180th meridian is still the edge of the grid, and
[Antimeridian](antimeridian.md) still applies.

They are also **not** adjustable. The bounds are constants of the specification,
transcribed from the paper — a grid whose extent depended on configuration
would not give the same index for the same position in two installations, which
would defeat the point of a global reference system.
