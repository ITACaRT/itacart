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

# Cell Filling / Polyfill

Covering a geometry with cells. This is the operation a cadastral pipeline runs
most, and the one with the most decisions hidden inside it.

```{code-cell} python
import itacart
from itacart import geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
cover = geometry.polyfill(parcel, resolution=9)
itacart.count_cells(cover)
```

The result is a compositional index — one string naming the whole cover, not a
list of identifiers.

## The cost tracks the boundary, not the area

`polyfill` descends the hierarchy over the sinusoidal plane instead of testing
cells one at a time. A coarse cell wholly inside the geometry is accepted whole
and never subdivided; a coarse cell wholly outside is discarded whole. Only
cells straddling the outline are refined further.

```{figure} ../../../_static/f7/f7_02_descent_cost.png
:alt: How the descent cost scales
:width: 100%

`docs/_static/f7/f7_02_descent_cost.png`
```

So the work is proportional to the length of the boundary rather than to the
enclosed area. Doubling a parcel's size does not double the cost of covering
it, which is what makes fine resolutions usable on large parcels at all.

## Densification is applied, not assumed

A straight edge on the plane is not a geodesic. An undensified long edge
therefore fills the wrong cells in the middle — and the failure is silent,
producing a plausible cover with no error raised.

```{figure} ../../../_static/f7/f7_05_chord_versus_geodesic.png
:alt: A chord against the geodesic between the same endpoints
:width: 100%

`docs/_static/f7/f7_05_chord_versus_geodesic.png`
```

`polyfill` densifies areal input first, at a threshold derived from the target
resolution. If you have already densified the geometry yourself, the step is
idempotent and costs you nothing:

```{code-cell} python
dense = geometry.densify_orthodromic(parcel, max_segment_m=5.0)
twice = geometry.densify_orthodromic(dense, max_segment_m=5.0)
len(parcel.exterior.coords), len(dense.exterior.coords), len(twice.exterior.coords)
```

Five positions become twenty-one at a five-metre threshold, and densifying the
result again changes nothing.

```{figure} ../../../_static/f7/f7_06_densification_idempotence.png
:alt: Densifying an already densified geometry
:width: 100%

`docs/_static/f7/f7_06_densification_idempotence.png`
```

## Counting without naming

When the answer you want is an area, naming the cells is wasted work:

```{code-cell} python
n = geometry.count_internal_cells(parcel, 9)
n, n * itacart.nominal_cell_area(9)
```

`count_internal_cells` runs the same descent and returns the count. It matches
the default containment mode exactly, which is why the number above and the
cover computed at the top of this page agree.

## A centre on the outline

Under `center` a cell is kept when its centre lies inside the geometry, and a
centre lying exactly on the outline is not inside. That is the common case for
a polygon drawn through cell centroids: every vertex is a centre, and so is
every lattice point along an edge. Such a centre used to be kept by no region
at all — not by the polygon, and not by its neighbour across the edge — and
whether a centre landed exactly on the line, or a unit in the last place to
either side, was decided by the C library that computed the projection.

So the point tested is the centre moved **one micrometre east and ten
nanometres north** on the sinusoidal plane. A centre on an edge then lies on
one side of it, and of two regions sharing the edge exactly **one** keeps the
cell: **west and south edges keep** the centres on them, **east and north
edges do not**. A square of 100 by 100 resolution-13 cells, drawn through the
centroids of its corner cells, holds 10 000 of them, and its two halves add up
to it:

```{code-cell} python
def walk(cell, direction, steps):
    for _ in range(steps):
        cell = itacart.get_neighbor(cell, direction)
    return cell

def through(*cells):
    return Polygon([itacart.cell_to_centroid(cell) for cell in cells])

sw = itacart.geo_to_cell(-47.0, -22.0, 13)
se, nw = walk(sw, "E", 100), walk(sw, "N", 100)
ne = walk(se, "N", 100)
west, east = walk(sw, "N", 50), walk(se, "N", 50)

whole = geometry.count_internal_cells(through(sw, se, ne, nw), 13)
south = geometry.count_internal_cells(through(sw, se, east, west), 13)
north = geometry.count_internal_cells(through(west, east, ne, nw), 13)
whole, south, north
```

The centres on the shared row went to the northern half. A centre farther from
the outline than the shift, about a micrometre, never changes side.

### The vertices drawn on the lattice

The shift settles a centre on the outline only if the outline is where it was
drawn. The vertices reach the plane through the projection, whose last bits
differ between C libraries by a few nanometres, and the shift itself puts some
centres that close to an edge: a centre the shift's own size from an edge lands
on it. So, under `center`, before anything is tested, a projected vertex within
**a micrometre**, in each coordinate, of a point of the **five-millimetre plane
grid** is moved onto that point, and one within a micrometre of the centroid of
a prime-meridian triangle — a third of a side above the triangle's base — onto
that centroid.

The grid holds **every corner and every centre of an ordinary cell, at every
resolution**, so a polygon drawn through the centroids or the corners of any
resolution's cells is read as drawn, with the same numbers on every machine: the
move uses only division, rounding to a whole number, multiplication and
addition. A corner or centroid this package returns in degrees projects back at
most 7.5 nanometres from its point, as measured, against the micrometre
allowed. The count and the fill under `center` take the same step; `intersects`
and `contains` read the vertices as projected.

Not moved:

- the centroid of a cell that absorbs the domain border, and of a polar cap,
  which the geometry engine computes from the cell's ring;
- any vertex not drawn on the lattice, which includes a lattice vertex written
  with **too few decimals**: rounding the latitude moves the abscissa too, in
  proportion to the longitude, and eleven decimals of a degree left a vertex
  1.6 micrometres off its point near the antemeridian. Twelve decimals, or a
  float written out in full, keep every lattice vertex within reach.

A polygon the move makes invalid — two vertices a micrometre apart about one
grid point become one — is refused, as any invalid polygon is, and not repaired.

### What the rule decides otherwise than the exact count

The exact count is that of the lattice points a polygon keeps when each is moved
an infinitesimal distance in the shift's direction. With the vertices on the
lattice and no edge densified, the count departs from it in two cases only, and
**in both it is the same on every machine**. Measured on one Linux machine:

- **An edge along the shift.** The shift points 1:100 north of east, and a
  centre on an edge running that way stays on it after the shift; its side is
  settled by the rounding of the shifted point, not by the rule. At 22 degrees
  of latitude, in three quadrants, on edges carrying nine lattice points between
  their ends, the count departed from the exact one by up to **ten cells** for
  edges **within 0.053 degrees** of that direction, and by none from 0.06
  degrees on.
- **A lattice point nearer an edge than the shift.** A centre that is not on
  the outline but nearer to it than the shift's component across the edge — at
  most a micrometre — changes side; one exactly that far lands on the edge, a
  near tie the arithmetic decides. On a polygon drawn through cell centroids
  that takes an edge whose direction on the lattice has no short step: edges of
  130 to 500 m with no common factor in their steps moved the count by up to
  **four cells** each, and edges up to 80 m by none. On a polygon whose
  vertices are not cell centres, this is all the shift changes: against the
  rule without it, the count of one of 160 random polygons at resolutions 6 to
  10 moved by one cell, and of 34 of 39 at resolution 13, where a cell is a
  centimetre wide, by **at most six cells** in counts of up to two billion.

### What still depends on the machine

Both are points the move does not reach:

- **An edge the fill or the count densifies.** A polygon arrives in longitude
  and latitude, and an edge at least as long as the densification threshold —
  **505.1 m at resolution 13, one kilometre at the others** — is read as a
  geodesic and gains points that are not on the lattice. Near them the decision
  is the platform's again. A polygon whose edges are straight on the plane is
  counted as drawn only while every edge stays below the threshold, and with a
  margin, because the length is computed with the C library too. On a square of
  20 km at resolution 13, with the C library's results moved by a unit in the
  last place to stand for four other platforms, a vertex every 400 m kept the
  count exact on all four; the corners alone let it move by one cell on three.
- **A vertex not on the lattice.** A vertex the move does not reach is where the
  projection puts it, give or take the last bits, and a centre that lands within
  a few nanometres of an edge through it is decided by them. Nothing makes that
  systematic, as a tie does, but at resolution 13 nothing makes it rare either:
  on average a centre lies within a nanometre of an edge once in fifty
  kilometres of outline. Of 39 random polygons of 100 to 400 m at resolution 13,
  with their coordinates moved by a unit in the last place to stand for six
  other platforms, **two moved by one cell**, on one or two of the six. The
  centroids of the border-absorbing cells and of the polar caps are vertices of
  this kind.

## The prime-meridian column

Cells on the prime meridian are triangles straddling the line, so the ordinary
descent cannot test them and a separate walk handles them. They carry the
eastern spelling, which is the only one the grammar admits, and their
containment is decided against the triangle rather than against a parallelogram
that is not there. Nothing about this is visible from the calling side —
it is here because a cover that quietly skipped that column would be wrong
along one meridian.

## Compact the result

A cover at a fine resolution is mostly solid interior, and the interior folds:

```{code-cell} python
compact = itacart.compact_cells(cover)
itacart.count_cells(cover), itacart.count_cells(compact)
```

[Compaction & Normalization](../index_and_hierarchy/compaction_and_normalization.md)
is the general treatment;
[Parcel Representation](../cadastral_workflows/parcel_representation.md) is the
cadastral one.

## Workers

`polyfill` runs in the calling thread whatever `n_jobs` says, and no extra
changes that. The parameter is accepted and validated, so calls that pass it
keep working. Until release 1.2.0 the fill spread its descent over threads, and
two threads made it slower than one: every Shapely call in the descent handed
the interpreter lock to the other thread.

Counting is where workers pay. `count_internal_cells` asks a whole level of the
descent at a time, in arrays, and with `n_jobs` above 1 it deals a large
descent out to processes. They are started by `spawn` on every system, so a
script whose count reaches them needs an `if __name__ == "__main__":` guard
around its entry point. A notebook does not, and neither does this page, which
counts a one-degree box on two workers:

```{code-cell} python
from shapely.geometry import box

region = box(-48.0, -16.0, -47.0, -15.0)
geometry.count_internal_cells(region, 9, n_jobs=2) == geometry.count_internal_cells(region, 9)
```

A count too small to repay starting the processes never leaves the calling
one, so `n_jobs` changes the time a count takes and never the count.
