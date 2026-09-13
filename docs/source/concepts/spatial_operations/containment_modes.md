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

# Containment Modes

"Inside" is three different questions, and the answer to each is a different
cover. Choosing between them is the most consequential decision in a cadastral
pipeline that uses this grid.

```{code-cell} python
import itacart
from itacart import geometry
from shapely.geometry import Polygon

parcel = Polygon([(-46.6329, -23.5509), (-46.6327, -23.5509),
                  (-46.6327, -23.5507), (-46.6329, -23.5507)])
{mode: itacart.count_cells(geometry.polyfill(parcel, 9, containment=mode))
 for mode in ("contains", "center", "intersects")}
```

**`contains`** keeps only cells wholly inside the geometry. **`center`**, the
default, keeps cells whose centre falls inside. **`intersects`** keeps every
cell the geometry touches at all.

## They nest, by construction

`contains` ⊆ `center` ⊆ `intersects`, and the nesting is a property of how
acceptance is decided rather than something that happens to hold on the
examples above:

```{code-cell} python
covers = {m: set(itacart.decompose(geometry.polyfill(parcel, 9, containment=m)))
          for m in ("contains", "center", "intersects")}
(covers["contains"] <= covers["center"], covers["center"] <= covers["intersects"])
```

```{figure} ../../../_static/f7/f7_01_containment_chain.png
:alt: The three containment modes nesting
:width: 100%

`docs/_static/f7/f7_01_containment_chain.png`
```

## The centre, not the anchor

`center` tests the cell's centre, and the choice is forced rather than
stylistic. A cell's anchor is a vertex and lies on its own border, so an anchor
test would award a cell on the strength of a point it shares with three
neighbours — and `contains` would stop being a subset of `center` for every
cell whose anchor happened to sit on the outline. The centre is interior to the
cell and has neither problem.

That is the same distinction
[Cells, Vertices, and Spatial Identity](../grid_fundamentals/cells_vertices_and_spatial_identity.md)
draws, arriving here as a concrete consequence.

## Which one to use

**`center` for area.** It is the default because it is the unbiased one: cells
lost along the boundary are traded against cells gained, and the estimate
converges on the true area as resolution increases.
[Area Estimation & Error Bounds](../cadastral_workflows/area_estimation_and_error_bounds.md)
bounds that convergence.

**`contains` when a false positive is expensive.** Everything in the cover is
certainly inside the parcel. Use it when the cover will be treated as a claim —
a registered footprint, a taxable extent — and understate rather than overstate.

**`intersects` when a false negative is expensive.** Nothing touching the
parcel is missed. Use it for search and indexing, where the cover is a filter
that a precise test will run behind.

## An empty cover is not always an answer

A geometry small enough to hold no cell centre produces an empty `center`
cover. The package raises rather than returning it, because an empty cover
silently treated as "no overlap" is the kind of error that reaches a land
registry. Under `intersects` the same geometry returns the cells it touches.

That case is not hypothetical: it is exactly what happens to a reingested
border outline in the last row of the lattice, which is documented in full
under [Geographic vs Grid Geometry](../data_and_interoperability/geographic_vs_grid_geometry.md).

## Carry the mode with the data

A cell set is not self-describing: the same parcel at the same resolution has
three legitimate covers, and they differ here by more than a quarter. A stored
cover that does not record which mode produced it cannot be compared with
another, and cannot be audited at all.
