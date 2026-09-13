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

## Parallel filling

With the `parallel` extra installed, the descent can be spread across workers
by passing `n_jobs`. It is worth it for large parcels at fine resolutions and
not worth it otherwise — the descent is already cheap relative to the area it
covers.
