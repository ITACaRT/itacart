# The ITACaRT Grid

The globe is divided into four quadrants — northeast, northwest, southeast,
southwest — by the prime meridian and the equator. Each quadrant is tessellated
into parallelogram cells, and every cell has an address that locates it within
its quadrant at every resolution at once.

## Direct surface tessellation

A discrete global grid has to get from a curved planet to a finite set of named
pieces, and almost every system does it in two steps: approximate the planet
with a polyhedron, then subdivide the polyhedron's faces. The approach is
computationally cheap and it has a large literature behind it.

ITACaRT does it in one step. The cells are laid out on the WGS84 ellipsoid
itself, from base and height measurements taken on that ellipsoid. There is no
polyhedron and no sphere in between.

```{figure} ../../../_static/f1/f1_07_ellipsoid_vs_sphere.png
:alt: Ellipsoid compared against a spherical approximation
:width: 100%

The divergence between the ellipsoid and a spherical approximation of it,
measured across latitude. `docs/_static/f1/f1_07_ellipsoid_vs_sphere.png`
```

The cost is arithmetic: distances and areas come from geodesic solutions
evaluated numerically rather than from closed forms on a sphere. The gain is
that survey data, which is already referred to WGS84, meets a grid built on the
same figure — no reprojection sits between the measurement and the cell.

This choice is why ITACaRT does not meet the parts of the OGC standard that
presume a polyhedral interface. That is a deliberate divergence rather than an
omission, and [Conformance Status](../standards_and_conformance/conformance_status.md)
records it as one.

## Why parallelograms

The natural cell shape for a grid of this kind is a square, and square cells in
this projection were tried in earlier work. They distort badly: the projection
does not preserve angles, so a square drawn in the projection plane is not a
square on the ellipsoid, and the further from the quadrant origin the worse it
gets.

ITACaRT takes the shape the projection actually produces. Cells originating at
the centre of a quadrant approximate a rhombus, standing at roughly forty-five
degrees, and the angle changes across the surface as the cells follow the
curvature. Accepting that shape rather than fighting it is what keeps the areas
exact. [Base Angle & Shape](../geometry/base_angle_and_shape.md) measures how
the angle varies; [The Sinusoidal Construction](the_sinusoidal_construction.md)
is where the shape comes from.

Parallelograms are also simple polygons, convex and four-sided, which is what
the standard asks of a cell geometry and what a surveyor's tooling expects.

## What the grid is for

The design criteria were chosen for cadastral work specifically, and they pull
against each other: geodetic precision against computational convenience, equal
area against uniform shape, standards conformance against usability. Where they
conflict, the paper resolves in favour of the cadastral requirement — area has
legal and fiscal weight, so it is the quantity that is preserved exactly and
the others that bend.

The consequences of that ordering are visible throughout the grid: in the
decimal areas, which make one cell one unit of account; in the Cartesian-like
addressing, which keeps a surveyor's intuitions working; and in the boundary
cells, which are specified departures rather than unhandled edge cases.
