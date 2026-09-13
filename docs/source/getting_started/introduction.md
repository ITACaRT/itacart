# Introduction

ITACaRT divides the surface of the Earth into cells of equal area, gives every
cell a name, and makes that name sufficient to recover where the cell is and
what shape it has. A coordinate pair says where something is; a cell index says
which piece of the world it occupies, and pieces can be compared, nested,
counted and stored.

## The problem it was built for

Cadastral work is measurement with consequences. A parcel boundary decides
where one person's property stops and another's begins, and a parcel's area
carries tax and legal weight. That makes area a quantity the representation has
to preserve rather than approximate.

Most discrete global grid systems project a polyhedron — an icosahedron, a
cube — onto a sphere, then subdivide its faces. The approach is computationally
convenient and it is the right trade for most applications. It is the wrong
trade here twice over: the sphere is not the figure that survey control is
referred to, and the projection distorts area in a way that varies across the
globe.

ITACaRT tessellates the WGS84 ellipsoid directly. There is no intermediate
polyhedron and no spherical approximation. The cost is computational: the grid
is built from the ellipsoidal sinusoidal projection and from geodesic
solutions, evaluated numerically. The gain is that a cell's area is the area it
is specified to have, and that the grid sits on the same figure the survey data
already refers to.

## What a cell name looks like

An index is a string, and it nests:

```
SE(1400/0374(3(C2(3))))
```

Parentheses descend a resolution; the leading pair of letters names the
quadrant. The same notation addresses more than one cell at a time, so a whole
vector feature — a parcel, a lot, a block — is one string rather than an
unstructured list of identifiers. [Index Structure](../concepts/index_and_hierarchy/index_structure.md)
takes the notation apart.

## What this documentation assumes

That you know what a coordinate reference system is and roughly what a
projection does to area, and that you do not already know this grid.
[Highlights](highlights.md) states what the grid guarantees and points at the
page that establishes each guarantee; [Quick Start](quick_start.md) gets you to
a first cell in a dozen lines. The conceptual material in
[Concepts & Guides](../concepts/index.md) is where each claim is argued rather
than asserted.

The grid is described in full in the paper the package implements; see
[Citing](../index.rst) on the front page for the reference.
