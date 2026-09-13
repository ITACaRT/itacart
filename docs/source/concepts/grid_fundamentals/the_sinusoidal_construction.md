# The Sinusoidal Construction

The ITACaRT grid is constructed directly on the ellipsoid rather than
by first projecting the Earth onto a plane. Its geometry is defined
from physical distances on the ellipsoidal surface: cell height
follows the true north–south distance associated with latitude, while
the east–west coordinate follows the true longitudinal distance from the
Greenwich meridian at that latitude.

These two quantities have exactly the geometric behavior expressed by an
equal-area sinusoidal projection. The sinusoidal plane is therefore a
convenient representation of the grid's intrinsic ellipsoidal construction,
not the surface on which the grid is originally defined.

This correspondence is what makes the equal-area property especially
useful: the planar representation preserves the area of the cells
already defined on the ellipsoid, rather than creating equal-area cells
by projection.

## Base and height, taken from the ellipsoid

Start with a cell's lower-left vertex in the northeast quadrant and call its
projected coordinates $\{x, y\}$. Let $l$ be the length of the horizontal base.
In this construction the vertical height is also $l$, and the remaining three
vertices follow:

| Vertex | Coordinates |
|---|---|
| lower-left (origin) | $\{x,\; y\}$ |
| lower-right | $\{x + l,\; y\}$ |
| upper-right | $\{x,\; y + l\}$ |
| upper-left | $\{x - l,\; y + l\}$ |

```{figure} ../../../_static/fpaper/fpaper_01_quadrant_and_parallelogram_arrangement.png
:alt: Quadrant division and the parallelogram's coordinate arrangement
:width: 100%

Figure 1 of the paper: the quadrant division and the coordinates of a
cell's four vertices.
`docs/_static/fpaper/fpaper_01_quadrant_and_parallelogram_arrangement.png`
```

The upper-right vertex sits directly above the origin, so the figure leans one
base-length to the left over one base-length of rise — the forty-five degree
lean. Cells in the other quadrants are obtained by reflection: across the
x-axis for the southern quadrants, across the y-axis for the western ones.

The base length and the height are not chosen in the plane and then projected
down. They are measured on the ellipsoid — the base along a parallel, the
height along a meridian — and the projection is what results from holding both
fixed across the grid. That ordering is the whole trick.

```{figure} ../../../_static/f1/f1_01_graticule.png
:alt: The graticule of the parallels-plane projection
:width: 100%

The graticule the construction produces.
`docs/_static/f1/f1_01_graticule.png`
```

## Why equal area follows

Holding base and height fixed while taking both from the ellipsoid is the
defining property of the sinusoidal projection — the parallel-planes
projection, in the paper's terms. A sinusoidal projection is equal-area by
construction: it preserves the area of every region it maps, paying for that
with shape distortion that grows away from the central meridian.

For a grid, that trade is exactly the right way round. Shape distortion changes
what a cell looks like; area distortion would change what a cell is worth.
Since the cells here are defined as figures of constant base and height in an
equal-area projection, every cell encloses the same area on the ellipsoid, and
the distortion is absorbed into the lean of the parallelogram rather than into
its size.

[Equal Area](equal_area.md) states the guarantee precisely, including where it
stops holding.

## The projection is computed, not delegated

The projection and the geodesic solutions are evaluated from the paper's
equations inside `itacart.geodesy`. PROJ is never called at run time, and there
is no compiled geospatial stack behind the package — `shapely` is the only hard
dependency.

That is a verifiability decision as much as a packaging one. The equations that
define the grid are in the source, at the resolution the paper states them, and
the round trip between geography and the projection is exercised directly:

```{figure} ../../../_static/f1/f1_02_roundtrip.png
:alt: Round-trip error between geography and the projection
:width: 100%

The round trip from geodetic coordinates into the projection and back.
`docs/_static/f1/f1_02_roundtrip.png`
```

## What this page does not settle

That a cell has a constant area says nothing yet about which point of the cell
its address names, or about how a cell divides into smaller cells. Those are
[Cells, Vertices, and Spatial Identity](cells_vertices_and_spatial_identity.md)
and [Resolution & Scale](resolution_and_scale.md) respectively.
